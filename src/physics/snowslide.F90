!>----------------------------------------------------------
!! Standalone SNOWSLIDE gravitational snow redistribution module
!!
!! Implementation of the SnowSlide model (Bernhardt and Schulz, 2010)
!! Adapted from FSM2-specific SNOWSLIDE_interface.F90 to work with
!! any snow model (FSM, SNOWPACK, etc.) using HICAR's native grid indexing.
!!
!! The algorithm routes excess snow from steep cells to lower neighbors
!! using an elevation-weighted D8 scheme, iterating multiple times with
!! MPI halo exchanges to allow cascading avalanches across MPI boundaries.
!!
!! @author Dylan Reynolds (adapted from Louis Quéno's FSM implementation)
!!----------------------------------------------------------
module module_snowslide
    use icar_constants
    use domain_interface,    only : domain_t
    use options_interface,   only : options_t
    use meta_data_interface, only : meta_data_t
    use variable_interface,  only : variable_t

    implicit none
    private

    public :: snowslide_var_request, snowslide_init, snowslide_step

    ! Module state
    integer :: Nx_local, Ny_local         ! Local grid dimensions (interior + 2 halo cells)
    integer :: its, ite, jts, jte         ! Tile/interior bounds
    integer :: ims, ime, jms, jme         ! Memory bounds (tile + halo); used to
                                          ! give assumed-shape helper dummies the
                                          ! correct lower bound so `i_dom` / `j_dom`
                                          ! remain valid indices inside the helpers.

    integer, allocatable  :: index_sorted(:,:)  ! (Nx*Ny, 2) — grid points sorted by decreasing elevation
    real, allocatable     :: dem_local(:,:)      ! Local DEM (m), in local (transposed) indexing
    real, allocatable     :: slope_local(:,:)    ! Local slope (degrees)
    real, allocatable     :: shd_local(:,:)      ! Local snow holding depth (m)

    ! Mirror the four module-scope allocatables on device. Populated once
    ! in snowslide_init; read-only thereafter from GPU regions. The
    ! `!$acc declare create` attaches the descriptor; `!$acc update device`
    ! in snowslide_init pushes the data over after allocation + fill.
    !$acc declare create(dem_local, slope_local, shd_local, index_sorted)

    ! Working arrays for routing (persist between calls for timing logic)
    type(variable_t) :: SD_0_var, Sice_0_var    ! For halo exchange

    ! Number of routing iterations per snowslide call. Cadence (how often
    ! snowslide_step is invoked) is managed by the caller — see
    ! sm_driver.F90's last_snowslide_time / next_snowslide_time gate.
    integer :: n_slide_iters = 10

    ! Layer info (detected at init)
    logical :: use_snowpack = .false.              ! Using SNOWPACK-style volume fractions

    ! Parameters (from PARAM_SNOWSLIDE)
    real, parameter :: dyn_ratio    = 0.7         ! Dynamic snow holding depth ratio
    real, parameter :: rho_deposit  = 300.0       ! Avalanche deposit density (kg/m³)
    real, parameter :: slope_min    = 25.0        ! Minimum slope for initial avalanche (degrees)
    real, parameter :: Shd_min      = 0.05        ! Minimum snow holding depth (m)
    real, parameter :: rho_ice_c    = 917.0       ! Ice density (kg/m³)
    real, parameter :: rho_water_c  = 1000.0      ! Water density (kg/m³)

    real :: Ds_min_fsm = 0.02   ! Min layer thickness (m)
    real :: Ds_surflay = 0.5    ! Max surface fine-layer zone (m)

    ! Union of all 3D per-layer snow kVARS across FSM + SNOWPACK. Any not
    ! allocated in the current build is skipped via vidx<=0 checks, so this
    ! list is safe for both snow-model builds. Kept parallel with
    ! snow_3d_fresh_defaults so the SNOWPACK reconcile can initialize a
    ! fresh surface element with a single loop.
    ! NOTE: kVARS%X is not a constant expression (kVARS is a module variable,
    ! not a parameter), so snow_3d_vars can't be declared `parameter`. It is
    ! populated at runtime in snowslide_init.
    integer, parameter :: n_snow_3d = 20
    integer :: snow_3d_vars(n_snow_3d)

    ! Fresh-snow defaults for a SNOWPACK avalanche deposit element. Values
    ! chosen to match SNOWPACK's new-snowfall init path at
    ! build/external/SNOWPACK/fortran/snowpack_driver.F90:666-697. Entries
    ! marked PATCH are overwritten per-cell in snowslide_reconcile_snowpack
    ! (Ds, Vol_Frac_I/A from rho_deposit, depositionDate from sim_time,
    ! snow_temperature/_i from skin_temperature, N3 from density piecewise).
    ! Rg/Rb are in millimeters, matching snowpack_driver.F90's new_rg =
    ! new_snow_grain_size/2 convention (default 0.2/2 = 0.1 mm) and the
    ! bond-radius rule new_rb = new_rg/4.
    real, parameter :: snow_3d_fresh_defaults(n_snow_3d) = &
        [0.0,          &  ! Ds                 PATCH
         0.0,          &  ! Sice               (FSM, skipped in SNOWPACK build)
         0.0,          &  ! Sliq               (FSM, skipped in SNOWPACK build)
         273.15,       &  ! snow_temperature   PATCH (falls back to 273.15)
         0.0,          &  ! Vol_Frac_I         PATCH
         0.0,          &  ! Vol_Frac_W
         0.0,          &  ! Vol_Frac_A         PATCH
         0.0,          &  ! Vol_Frac_S
         0.0,          &  ! Vol_Frac_WP
         273.15,       &  ! snow_temperature_i PATCH (falls back to 273.15)
         0.1,          &  ! Rg (mm) = new_snow_grain_size / 2
         0.025,        &  ! Rb (mm) = Rg / 4
         1.0,          &  ! Dd (dendricity)
         0.5,          &  ! Sp (sphericity)
         0.0,          &  ! mk (SN_NEW_SNOW_MARKER = 0)
         0.0,          &  ! mass_hoar
         0.0,          &  ! CDot
         0.0,          &  ! snow_stress
         1.75,         &  ! N3                 PATCH (from rho_deposit piecewise)
         0.0]             ! depositionDate     PATCH

contains

    !>----------------------------------------------------------
    !! Request allocation of variables needed by SNOWSLIDE
    !!----------------------------------------------------------
    subroutine snowslide_var_request(options)
        implicit none
        type(options_t), intent(inout) :: options

        call options%alloc_vars( &
            [kVARS%dSWE_slide, kVARS%shd, kVARS%terrain, kVARS%slope_angle, &
             kVARS%snow_height, kVARS%snow_water_equivalent])

        call options%restart_vars([kVARS%dSWE_slide])
    end subroutine snowslide_var_request


    !>----------------------------------------------------------
    !! Initialize SNOWSLIDE module
    !!----------------------------------------------------------
    subroutine snowslide_init(domain, options)
        implicit none
        type(domain_t), intent(inout) :: domain
        type(options_t), intent(in)   :: options

        integer :: i, j
        real, parameter :: piconst = 3.14159265358979323846

        ! Store grid bounds
        its = domain%grid%its; ite = domain%grid%ite
        jts = domain%grid%jts; jte = domain%grid%jte
        ims = domain%grid%ims; ime = domain%grid%ime
        jms = domain%grid%jms; jme = domain%grid%jme

        ! Local grid: interior + 1 halo cell on each side (HICAR native indexing)
        Nx_local = (ite - its + 1) + 2
        Ny_local = (jte - jts + 1) + 2

        ! Allocate local arrays
        if (allocated(dem_local))   deallocate(dem_local)
        if (allocated(slope_local)) deallocate(slope_local)
        if (allocated(shd_local))   deallocate(shd_local)
        if (allocated(index_sorted)) deallocate(index_sorted)

        allocate(dem_local(Nx_local, Ny_local))
        allocate(slope_local(Nx_local, Ny_local))
        allocate(shd_local(Nx_local, Ny_local))
        allocate(index_sorted(Nx_local * Ny_local, 2))

        ! Copy terrain data from domain to local transposed arrays
        ! Fill local arrays including halo cells from domain (halo populated by MPI exchange).
        ! This ensures buffer cells see correct neighbor elevations at MPI boundaries.
        dem_local = 0.0
        slope_local = 0.0
        shd_local = 10000.0  ! Default: very large holding depth (no slides)

        dem_local(1:Nx_local, 1:Ny_local) = &
            domain%vars_2d(domain%var_indx(kVARS%terrain)%v)%data_2d(its-1:ite+1, jts-1:jte+1)

        if (domain%var_indx(kVARS%slope_angle)%v > 0) then
            slope_local(1:Nx_local, 1:Ny_local) = &
                domain%vars_2d(domain%var_indx(kVARS%slope_angle)%v)%data_2d(its-1:ite+1, jts-1:jte+1) &
                * 180.0 / piconst   ! Convert radians to degrees
        endif

        shd_local = min(shd_local, domain%vars_2d(domain%var_indx(kVARS%shd)%v)%data_2d(its-1:ite+1, jts-1:jte+1))

        call sort_dem_local()

        ! Push the populated module arrays to device. They are read-only
        ! from every subsequent snowslide GPU region.
        !$acc update device(dem_local, slope_local, shd_local, index_sorted)

        ! Initialize variable_t objects for halo exchange
        call SD_0_var%initialize(kVARS%density, domain%grid2d)
        call Sice_0_var%initialize(kVARS%density, domain%grid2d)

        ! Detect layer representation and max layers
        use_snowpack = (options%physics%snowmodel == kSM_SNOWPACK)

        Ds_min_fsm = options%sm%fsm_ds_min
        Ds_surflay = options%sm%fsm_ds_surflay

        ! Populate the union var list (can't be a parameter — kVARS%X is not
        ! a constant expression). Order must stay parallel with snow_3d_fresh_defaults.
        snow_3d_vars = [kVARS%Ds,                                                           &
                        kVARS%Sice, kVARS%Sliq, kVARS%snow_temperature,                     &
                        kVARS%Vol_Frac_I, kVARS%Vol_Frac_W, kVARS%Vol_Frac_A,                &
                        kVARS%Vol_Frac_S, kVARS%Vol_Frac_WP, kVARS%snow_temperature_i,      &
                        kVARS%Rg, kVARS%Rb, kVARS%Dd, kVARS%Sp,                             &
                        kVARS%mk, kVARS%mass_hoar, kVARS%CDot, kVARS%snow_stress,           &
                        kVARS%N3, kVARS%depositionDate]

    end subroutine snowslide_init


    !>----------------------------------------------------------
    !! Run one SNOWSLIDE step. Call cadence is controlled by the caller
    !! (sm_driver.F90) — this routine always executes when invoked.
    !!----------------------------------------------------------
    subroutine snowslide_step(domain, options)
        implicit none
        type(domain_t), intent(inout)  :: domain
        type(options_t), intent(in)    :: options

        real, allocatable :: SD_0(:,:), Sice_0(:,:)
        real, allocatable :: dSWE_slide_local(:,:)
        logical, allocatable :: snow_depo(:,:)
        integer :: iter
        integer :: i_s, i_e, j_s, j_e

        i_s = 2;  i_e = Nx_local - 1
        j_s = 2;  j_e = Ny_local - 1

        allocate(SD_0(Nx_local, Ny_local));      SD_0 = 0.0
        allocate(Sice_0(Nx_local, Ny_local));    Sice_0 = 0.0
        allocate(dSWE_slide_local(Nx_local, Ny_local)); dSWE_slide_local = 0.0
        allocate(snow_depo(Nx_local, Ny_local)); snow_depo = .False.

        call sync_snow_to_host(domain)

        ! ================================================================
        ! Iteration loop. Each iteration:
        !   Step 1: Route interior only — erode + deposit, buffer accumulates
        !   Step 2: Exchange domain Ds halos — neighbor sees our state
        !   Step 3: Exchange SD_0 halos — neighbor sees buffer deposits
        !   Step 4: Route buffer+halo — buffer erodes, halo deposits cascade in
        ! After all iterations: SNOW_LAYERING applies remaining SD_0
        ! ================================================================
        do iter = 1, n_slide_iters

            ! ---- Step 1: Route interior only ----
            call snowslide_route(domain, dem_local, slope_local, shd_local, &
                                 index_sorted, Nx_local, Ny_local, &
                                 SD_0, Sice_0, dSWE_slide_local, &
                                 snow_depo, .True.)

            ! ---- Step 2: Exchange domain Ds halos ----
            call sync_snow_halo_to_device(domain)
            call domain%halo%exch_var(domain%vars_2d(domain%var_indx(kVARS%snow_nlayers)%v),corners=.True.)
            call domain%halo%exch_var(domain%vars_3d(domain%var_indx(kVARS%Ds)%v),corners=.True.)
            if (use_snowpack) then
                call domain%halo%exch_var(domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_I)%v),corners=.True.)
                call domain%halo%exch_var(domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_W)%v),corners=.True.)
            else
                call domain%halo%exch_var(domain%vars_3d(domain%var_indx(kVARS%Sice)%v),corners=.True.)
                call domain%halo%exch_var(domain%vars_3d(domain%var_indx(kVARS%Sliq)%v),corners=.True.)
            end if
            call sync_snow_halo_to_host(domain)

            ! ---- Step 3: Exchange SD_0/Sice_0 halos ----
            call exch_slide_halo(domain, SD_0, Sice_0)

            ! Flag halo cells that received deposits from neighbor
            where(SD_0 > 0.0) snow_depo = .True.
            snow_depo(i_s:i_e, j_s:j_e) = .False.

            call snowslide_route(domain, dem_local, slope_local, shd_local, &
                                 index_sorted, Nx_local, Ny_local, &
                                 SD_0, Sice_0, dSWE_slide_local, &
                                 snow_depo, .False.)

        end do

        domain%vars_2d(domain%var_indx(kVARS%dSWE_slide)%v)%data_2d(its:ite, jts:jte) = &
            domain%vars_2d(domain%var_indx(kVARS%dSWE_slide)%v)%data_2d(its:ite, jts:jte) + &
            dSWE_slide_local(i_s:i_e, j_s:j_e)
        !$acc update device(domain%vars_2d(domain%var_indx(kVARS%dSWE_slide)%v)%data_2d(its:ite, jts:jte))

        ! ================================================================
        ! After all iterations: apply remaining deposits via SNOW_LAYERING
        ! ================================================================
        if (options%physics%snowmodel == kSM_FSM) then
#ifdef FSM
            ! Push FSM-relevant arrays to device, run the reconcile as a
            ! gang-parallel kernel, leave state on device. SD_0/Sice_0 are
            ! scoped onto device for the call only.
            !$acc enter data copyin(SD_0, Sice_0)
            !$acc update device( &
            !$acc   domain%vars_2d(domain%var_indx(kVARS%snow_nlayers)%v)%data_2di, &
            !$acc   domain%vars_2d(domain%var_indx(kVARS%fsnow)%v)%data_2d, &
            !$acc   domain%vars_2d(domain%var_indx(kVARS%snow_height)%v)%data_2d, &
            !$acc   domain%vars_2d(domain%var_indx(kVARS%snow_water_equivalent)%v)%data_2d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Ds)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Sice)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Sliq)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%snow_temperature)%v)%data_3d)
            call snowslide_reconcile_fsm(domain, SD_0, Sice_0)
            !$acc exit data delete(SD_0, Sice_0)
#endif
        else if (options%physics%snowmodel == kSM_SNOWPACK) then
#ifdef SNOWPACK
            ! snowslide_reconcile_snowpack runs on GPU as a `!$acc parallel loop`.
            ! Push all per-layer + per-cell arrays it reads/writes to device first;
            ! the kernel leaves the updated state on device. SD_0/Sice_0 are local
            ! host-side arrays so they get a scoped copyin for the call only.
            !$acc enter data copyin(SD_0, Sice_0)
            !$acc update device( &
            !$acc   domain%vars_2d(domain%var_indx(kVARS%snow_nlayers)%v)%data_2di, &
            !$acc   domain%vars_2d(domain%var_indx(kVARS%skin_temperature)%v)%data_2d, &
            !$acc   domain%vars_2d(domain%var_indx(kVARS%snow_height)%v)%data_2d, &
            !$acc   domain%vars_2d(domain%var_indx(kVARS%snow_water_equivalent)%v)%data_2d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Ds)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%snow_temperature)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%snow_temperature_i)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_I)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_W)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_A)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_S)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_WP)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Rg)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Rb)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Dd)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Sp)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%mk)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%mass_hoar)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%CDot)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%snow_stress)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%N3)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%depositionDate)%v)%data_3d)
            call snowslide_reconcile_snowpack(domain, SD_0, Sice_0)
            !$acc exit data delete(SD_0, Sice_0)
#endif
        end if

        deallocate(SD_0, Sice_0, dSWE_slide_local, snow_depo)

    end subroutine snowslide_step


    !>----------------------------------------------------------
    !! Sync domain snow arrays from device to host (before routing)
    !!----------------------------------------------------------
    subroutine sync_snow_to_host(domain)
        implicit none
        type(domain_t), intent(inout) :: domain

        ! Called ONCE at the top of snowslide_step. Pulls the full set of
        ! per-layer snow vars that the route phase will read/write on host
        ! via ablate_snow_layers (which shifts all 18 SNOWPACK vars or the 4
        ! FSM vars). Once on host, they persist through the iter loop —
        ! only halo-exchanged vars need re-pulling, via sync_snow_halo_to_host.
        !$acc update host( &
        !$acc   domain%vars_2d(domain%var_indx(kVARS%snow_nlayers)%v)%data_2di, &
        !$acc   domain%vars_2d(domain%var_indx(kVARS%skin_temperature)%v)%data_2d, &
        !$acc   domain%vars_3d(domain%var_indx(kVARS%Ds)%v)%data_3d)
        if (use_snowpack) then
            !$acc update host( &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%snow_temperature)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%snow_temperature_i)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_I)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_W)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_A)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_S)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_WP)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Rg)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Rb)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Dd)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Sp)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%mk)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%mass_hoar)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%CDot)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%snow_stress)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%N3)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%depositionDate)%v)%data_3d)
        else
            !$acc update host( &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Sice)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Sliq)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%snow_temperature)%v)%data_3d)
        end if
    end subroutine sync_snow_to_host


    !>----------------------------------------------------------
    !! Lighter variant for use *inside* the route iteration loop, after a
    !! halo exchange. Only re-pulls the arrays that the halo exchange just
    !! wrote (snow_nlayers, Ds, Sice/Sliq for FSM or Vol_Frac_I/W for
    !! SNOWPACK). Using the full sync_snow_to_host here would clobber
    !! route-side modifications to the non-exchanged vars (Tsn, Tsni, VFA,
    !! VFS, Rg, Rb, Dd, Sp, mk, ...) that ablate_snow_layers has shifted on
    !! host but not yet pushed to device.
    !!----------------------------------------------------------
    subroutine sync_snow_halo_to_host(domain)
        implicit none
        type(domain_t), intent(inout) :: domain

        !$acc update host( &
        !$acc   domain%vars_2d(domain%var_indx(kVARS%snow_nlayers)%v)%data_2di, &
        !$acc   domain%vars_3d(domain%var_indx(kVARS%Ds)%v)%data_3d)
        if (use_snowpack) then
            !$acc update host( &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_I)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_W)%v)%data_3d)
        else
            !$acc update host( &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Sice)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Sliq)%v)%data_3d)
        end if
    end subroutine sync_snow_halo_to_host

    subroutine sync_snow_halo_to_device(domain)
        implicit none
        type(domain_t), intent(inout) :: domain

        !$acc update device( &
        !$acc   domain%vars_2d(domain%var_indx(kVARS%snow_nlayers)%v)%data_2di, &
        !$acc   domain%vars_3d(domain%var_indx(kVARS%Ds)%v)%data_3d)
        if (use_snowpack) then
            !$acc update device( &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_I)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_W)%v)%data_3d)
        else
            !$acc update device( &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Sice)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Sliq)%v)%data_3d)
        end if
    end subroutine sync_snow_halo_to_device

    !>----------------------------------------------------------
    !! Sync modified domain snow arrays from host back to device
    !!----------------------------------------------------------
    subroutine sync_snow_to_device(domain)
        implicit none
        type(domain_t), intent(inout) :: domain

        !$acc update device( &
        !$acc   domain%vars_2d(domain%var_indx(kVARS%snow_height)%v)%data_2d, &
        !$acc   domain%vars_2d(domain%var_indx(kVARS%snow_water_equivalent)%v)%data_2d, &
        !$acc   domain%vars_2d(domain%var_indx(kVARS%dSWE_slide)%v)%data_2d, &
        !$acc   domain%vars_2d(domain%var_indx(kVARS%snow_nlayers)%v)%data_2di, &
        !$acc   domain%vars_3d(domain%var_indx(kVARS%Ds)%v)%data_3d)
        if (use_snowpack) then
            !$acc update device( &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_I)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_W)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_A)%v)%data_3d)
        else
            !$acc update device( &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Sice)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Sliq)%v)%data_3d)
        end if
    end subroutine sync_snow_to_device


    !>----------------------------------------------------------
    !! Core SNOWSLIDE routing algorithm
    !! Adapted from SNOWSLIDE_interface.F90 with layer-aware erosion
    !!----------------------------------------------------------
    subroutine snowslide_route(domain, dem, slope, shd, &
                               idx_sorted, Nx, Ny, &
                               snowdepth0, Sice0, dSWE_slide, &
                               snow_depo, interior_only)
        implicit none

        type(domain_t), intent(inout) :: domain
        integer, intent(in) :: Nx, Ny
        real, intent(in)    :: dem(Nx, Ny)              ! Terrain elevation (m)
        real, intent(in)    :: slope(Nx, Ny)            ! Slope (degrees)
        real, intent(in)    :: shd(Nx, Ny)              ! Snow holding depth (m)
        integer, intent(in) :: idx_sorted(Nx*Ny, 2)    ! Sorted grid point indices

        real, intent(inout)    :: snowdepth0(Nx, Ny)    ! Fresh avalanche deposit depth (m)
        real, intent(inout)    :: Sice0(Nx, Ny)         ! Fresh avalanche deposit mass (kg/m²)
        real, intent(inout)    :: dSWE_slide(Nx, Ny)    ! SWE change from slides (kg/m²)
        logical, intent(inout) :: snow_depo(Nx, Ny)     ! Cells receiving deposits
        logical, intent(in)    :: interior_only         ! If true, skip buffer and halo cells

        ! D8 neighbor offsets (ordering matches original S,N,W,E,SW,SE,NW,NE layout)
        integer, parameter :: di(8) = [-1,  1,  0,  0, -1, -1,  1,  1]
        integer, parameter :: dj(8) = [ 0,  0, -1,  1, -1,  1, -1,  1]

        integer :: i, j, n, d, ni, nj
        real    :: elev, neighbor_elev, delev_tot, dswe
        real    :: snowdepth_available, swe_available, snowdepth_updated
        real    :: w(8)
        real    :: Shd_corr(Nx, Ny)

        Shd_corr = max(shd, Shd_min)

        ! Process grid points from highest to lowest elevation
        do n = 1, Nx * Ny
            i = idx_sorted(n, 1)
            j = idx_sorted(n, 2)

            ! Skip buffer and halo cells when routing interior only
            if (interior_only .and. (i <= 2 .or. i >= Nx-1 .or. j <= 2 .or. j >= Ny-1)) cycle

            if (slope(i,j) < slope_min) cycle

            snowdepth_updated = get_snow_depth(domain, i, j) + snowdepth0(i,j)
            elev = dem(i,j) + snowdepth_updated

            ! Dynamic reduction of holding depth for cells receiving avalanche
            if (snow_depo(i,j)) then
                Shd_corr(i,j) = max(Shd_corr(i,j) * dyn_ratio, Shd_min)
            end if

            snowdepth_available = max(0.0, snowdepth_updated - Shd_corr(i,j))
            if (snowdepth_available <= epsilon(snowdepth_available)) cycle

            ! Compute elevation-weighted routing to 8 neighbors
            w = 0.0
            do d = 1, 8
                ni = i + di(d)
                nj = j + dj(d)
                if (ni < 1 .or. ni > Nx .or. nj < 1 .or. nj > Ny) cycle
                neighbor_elev = dem(ni,nj) + get_snow_depth(domain, ni, nj) + snowdepth0(ni,nj)
                w(d) = max(0.0, elev - neighbor_elev)
            end do

            delev_tot = sum(w)
            if (delev_tot <= epsilon(delev_tot)) cycle

            w = w / delev_tot

            ! Erode snow: first from fresh deposits, then from pack layers
            call erode_snow(domain, i, j, snowdepth_available, snowdepth0(i,j), Sice0(i,j), swe_available)
            dSWE_slide(i,j) = dSWE_slide(i,j) - swe_available

            ! Transport to 8 neighbors weighted by elevation difference
            do d = 1, 8
                if (w(d) <= epsilon(w(d))) cycle
                ni = i + di(d)
                nj = j + dj(d)
                dswe = w(d) * swe_available
                Sice0(ni,nj)      = Sice0(ni,nj) + dswe
                snowdepth0(ni,nj) = snowdepth0(ni,nj) + dswe / rho_deposit
                snow_depo(ni,nj)  = .TRUE.
                dSWE_slide(ni,nj) = dSWE_slide(ni,nj) + dswe
            end do
        end do

    end subroutine snowslide_route


    !>----------------------------------------------------------
    !! Compute snow depth from domain Ds at local grid position (i,j)
    !! Maps local→domain and sums layer thicknesses, matching FSM's
    !! sum(Ds(:,i,j)) * fsnow(i,j) (we assume fsnow=1 for snowslide)
    !!----------------------------------------------------------
    function get_snow_depth(domain, i_loc, j_loc) result(depth)
        implicit none
        type(domain_t), intent(in) :: domain
        integer, intent(in) :: i_loc, j_loc
        real :: depth
        integer :: i_dom, j_dom, nsnow

        i_dom = its + (i_loc - 2)
        j_dom = jts + (j_loc - 2)
        nsnow = domain%vars_2d(domain%var_indx(kVARS%snow_nlayers)%v)%data_2di(i_dom, j_dom)
        depth = sum(domain%vars_3d(domain%var_indx(kVARS%Ds)%v)%data_3d(i_dom, 1:nsnow, j_dom))
        ! FSM computes snow depth as sum(Ds) * fsnow
        if (domain%var_indx(kVARS%fsnow)%v > 0) then
            depth = depth * domain%vars_2d(domain%var_indx(kVARS%fsnow)%v)%data_2d(i_dom, j_dom)
        end if
    end function get_snow_depth


    !>----------------------------------------------------------
    !! Erode snow from a cell: first from fresh deposits (snowdepth0/Sice0),
    !! then from domain pack layers via ablate_snow_layers.
    !! Returns total SWE eroded for depositing on neighbors.
    !!----------------------------------------------------------
    subroutine erode_snow(domain, i_loc, j_loc, &
                          snowdepth_available, snowdepth0, Sice0, swe_eroded)
        implicit none
        type(domain_t), intent(inout) :: domain
        integer, intent(in)    :: i_loc, j_loc
        real, intent(in)       :: snowdepth_available  ! Total depth to erode (m)
        real, intent(inout)    :: snowdepth0            ! Fresh deposit depth (consumed first)
        real, intent(inout)    :: Sice0                 ! Fresh deposit mass
        real, intent(out)      :: swe_eroded            ! Total SWE eroded (kg/m²)

        integer :: i_dom, j_dom
        real :: depth_from_pack, swe_from_pack

        swe_eroded = 0.0
        i_dom = its + (i_loc - 2)
        j_dom = jts + (j_loc - 2)

        if (snowdepth0 >= snowdepth_available) then
            ! All excess comes from fresh deposits
            swe_eroded = snowdepth_available / max(snowdepth0, 1.0e-10) * Sice0
            snowdepth0 = snowdepth0 - snowdepth_available
            Sice0 = Sice0 - swe_eroded
        else
            ! Consume all deposits, then erode from pack layers
            swe_eroded = Sice0
            depth_from_pack = snowdepth_available - snowdepth0
            snowdepth0 = 0.0
            Sice0 = 0.0

            call ablate_snow_layers(domain, i_dom, j_dom, depth_from_pack, swe_from_pack)
            swe_eroded = swe_eroded + swe_from_pack
        end if

    end subroutine erode_snow


    !>----------------------------------------------------------
    !! Top-down layer removal operating directly on domain arrays.
    !! Erodes dhs depth from surface at domain cell (i_dom, j_dom).
    !! Shifts ALL allocated 3D snow variables, zeros above new nlay.
    !! Returns the SWE eroded.
    !!----------------------------------------------------------
    subroutine ablate_snow_layers(domain, i_dom, j_dom, dhs, swe_eroded)
        implicit none
        type(domain_t), intent(inout) :: domain
        integer, intent(in) :: i_dom, j_dom
        real, intent(in)    :: dhs              ! Depth to erode (m)
        real, intent(out)   :: swe_eroded       ! SWE removed (kg/m²)

        integer :: k, k_eaten, nlay, nlay_new
        integer :: vidx_Ds, vidx_nlay
        integer :: vidx_Sice, vidx_Sliq, vidx_Tsn
        integer :: vidx_Tsni, vidx_VFI, vidx_VFW, vidx_VFA, vidx_VFS, vidx_VFWP
        integer :: vidx_Rg, vidx_Rb, vidx_Dd, vidx_Sp
        integer :: vidx_mk, vidx_mh, vidx_CDot, vidx_sns, vidx_N3, vidx_dep
        real :: depth_accum, Ds_k, rho_k, ds_old, ds_new, scale

        swe_eroded = 0.0
        vidx_nlay = domain%var_indx(kVARS%snow_nlayers)%v
        nlay = domain%vars_2d(vidx_nlay)%data_2di(i_dom, j_dom)
        if (nlay <= 0 .or. dhs <= 0.0) return

        vidx_Ds = domain%var_indx(kVARS%Ds)%v
        if (use_snowpack) then
            vidx_Sice = -1; vidx_Sliq = -1; vidx_Tsn = -1
        else
            vidx_Sice = domain%var_indx(kVARS%Sice)%v
            vidx_Sliq = domain%var_indx(kVARS%Sliq)%v
            vidx_Tsn  = domain%var_indx(kVARS%snow_temperature)%v
        end if
        if (use_snowpack) then
            vidx_Tsn   = domain%var_indx(kVARS%snow_temperature)%v
            vidx_Tsni  = domain%var_indx(kVARS%snow_temperature_i)%v
            vidx_VFI   = domain%var_indx(kVARS%Vol_Frac_I)%v
            vidx_VFW   = domain%var_indx(kVARS%Vol_Frac_W)%v
            vidx_VFA   = domain%var_indx(kVARS%Vol_Frac_A)%v
            vidx_VFS   = domain%var_indx(kVARS%Vol_Frac_S)%v
            vidx_VFWP  = domain%var_indx(kVARS%Vol_Frac_WP)%v
            vidx_Rg    = domain%var_indx(kVARS%Rg)%v
            vidx_Rb    = domain%var_indx(kVARS%Rb)%v
            vidx_Dd    = domain%var_indx(kVARS%Dd)%v
            vidx_Sp    = domain%var_indx(kVARS%Sp)%v
            vidx_mk    = domain%var_indx(kVARS%mk)%v
            vidx_mh    = domain%var_indx(kVARS%mass_hoar)%v
            vidx_CDot  = domain%var_indx(kVARS%CDot)%v
            vidx_sns   = domain%var_indx(kVARS%snow_stress)%v
            vidx_N3    = domain%var_indx(kVARS%N3)%v
            vidx_dep   = domain%var_indx(kVARS%depositionDate)%v
        end if

        ! Walk from the surface consuming whole layers, accumulating the real
        ! mass removed from each one. This is the *only* correct way to get
        ! swe_eroded — the previous column-average (rho_bulk) form was wrong
        ! whenever top-layer density differs from the column mean, creating
        ! or destroying mass in snowslide_route's source→neighbor transfer.
        k = 0
        depth_accum = 0.0
        do while (k < nlay .and. depth_accum < dhs)
            k = k + 1
            Ds_k = domain%vars_3d(vidx_Ds)%data_3d(i_dom, k, j_dom)
            if (use_snowpack) then
                ! SNOWPACK: volume fractions are already mass per unit volume
                ! when multiplied by their phase density, no Ds involved.
                rho_k = domain%vars_3d(vidx_VFI)%data_3d(i_dom, k, j_dom) * rho_ice_c &
                      + domain%vars_3d(vidx_VFW)%data_3d(i_dom, k, j_dom) * rho_water_c
            else
                ! FSM: Sice/Sliq are mass per unit area (kg/m²), so
                ! density = (Sice + Sliq) / Ds.
                rho_k = (domain%vars_3d(vidx_Sice)%data_3d(i_dom, k, j_dom) &
                       + domain%vars_3d(vidx_Sliq)%data_3d(i_dom, k, j_dom)) &
                      / max(Ds_k, 1.0e-10)

            end if
            if (depth_accum + Ds_k <= dhs) then
                ! Whole layer consumed
                swe_eroded  = swe_eroded  + Ds_k * rho_k
            else
                ! Partial layer — take only (dhs - depth_accum) meters off the top.
                ! Advance depth_accum by the full Ds_k so the loop exits *and* the
                ! residual thickness `depth_accum - dhs = Ds_k - (dhs - prev_accum)`
                ! falls out of the formula below.
                swe_eroded  = swe_eroded  + (dhs - depth_accum) * rho_k
            end if
            depth_accum = depth_accum + Ds_k
        end do
        ! If the partial branch fired, layer k has a residue (depth_accum > dhs).
        ! Otherwise every touched layer was fully consumed.
        if (depth_accum > dhs + 1.0e-10) then
            k_eaten = k - 1
        else
            k_eaten = k
        end if
        nlay_new = max(nlay - k_eaten, 0)

        ! Partial erosion of the remaining top layer + shift up.
        !
        ! SNOWPACK-correctness note: Vol_Frac_* are intensive (volume fractions
        ! preserved when a partial slice is removed from the surface), so they
        ! fall through to the plain-copy branch. Microstructure vars (Rg, Rb,
        ! Dd, Sp, mk, ...) also travel with the layer — same copy branch. Only
        ! FSM's Sice/Sliq need extensive-mass scaling.
        if (nlay_new > 0) then
            ds_old = domain%vars_3d(vidx_Ds)%data_3d(i_dom, 1+k_eaten, j_dom)
            ds_new = max(0.0, depth_accum - dhs)
            scale  = ds_new / max(ds_old, 1.0e-10)

            domain%vars_3d(vidx_Ds)%data_3d(i_dom, 1, j_dom) = ds_new
            do k = 2, nlay_new
                domain%vars_3d(vidx_Ds)%data_3d(i_dom, k, j_dom) = &
                    domain%vars_3d(vidx_Ds)%data_3d(i_dom, k+k_eaten, j_dom)
            end do
        end if

        ! Tail: always zero every slot above nlay_new (covers k=1..kSNOW_GRID_Z
        ! when the full column was consumed).
        do k = nlay_new + 1, kSNOW_GRID_Z
            domain%vars_3d(vidx_Ds)%data_3d(i_dom, k, j_dom) = 0.0
        end do


        if (use_snowpack) then
            ! SNOWPACK: Vol_Frac_* and microstructure are intensive and
            ! ride along with the layer — plain copy via the 17-var
            ! helpers. No rescale needed for a partial top layer because
            ! density is independent of element thickness. Ds is already
            ! handled by the branch-agnostic block above.
            associate( &
                Tsn_3d  => domain%vars_3d(vidx_Tsn  )%data_3d, &
                Tsni_3d => domain%vars_3d(vidx_Tsni )%data_3d, &
                VFI_3d  => domain%vars_3d(vidx_VFI  )%data_3d, &
                VFW_3d  => domain%vars_3d(vidx_VFW  )%data_3d, &
                VFA_3d  => domain%vars_3d(vidx_VFA  )%data_3d, &
                VFS_3d  => domain%vars_3d(vidx_VFS  )%data_3d, &
                VFWP_3d => domain%vars_3d(vidx_VFWP )%data_3d, &
                Rg_3d   => domain%vars_3d(vidx_Rg   )%data_3d, &
                Rb_3d   => domain%vars_3d(vidx_Rb   )%data_3d, &
                Dd_3d   => domain%vars_3d(vidx_Dd   )%data_3d, &
                Sp_3d   => domain%vars_3d(vidx_Sp   )%data_3d, &
                mk_3d   => domain%vars_3d(vidx_mk   )%data_3d, &
                mh_3d   => domain%vars_3d(vidx_mh   )%data_3d, &
                CDot_3d => domain%vars_3d(vidx_CDot )%data_3d, &
                sns_3d  => domain%vars_3d(vidx_sns  )%data_3d, &
                N3_3d   => domain%vars_3d(vidx_N3   )%data_3d, &
                dep_3d  => domain%vars_3d(vidx_dep  )%data_3d, &
                ims => domain%ims, ime => domain%ime, jms => domain%jms, jme => domain%jme)
                if (nlay_new > 0) then
                    call shift_snowpack_layer_slot(                        &
                        Tsn_3d, Tsni_3d,                                   &
                        VFI_3d, VFW_3d, VFA_3d, VFS_3d, VFWP_3d,           &
                        Rg_3d, Rb_3d, Dd_3d, Sp_3d,                        &
                        mk_3d, mh_3d, CDot_3d, sns_3d, N3_3d, dep_3d,      &
                        i_dom, j_dom, 1, 1+k_eaten, ims, ime, jms, jme, kSNOW_GRID_Z)
                    do k = 2, nlay_new
                        call shift_snowpack_layer_slot(                    &
                            Tsn_3d, Tsni_3d,                               &
                            VFI_3d, VFW_3d, VFA_3d, VFS_3d, VFWP_3d,       &
                            Rg_3d, Rb_3d, Dd_3d, Sp_3d,                    &
                            mk_3d, mh_3d, CDot_3d, sns_3d, N3_3d, dep_3d,  &
                            i_dom, j_dom, k, k+k_eaten, ims, ime, jms, jme, kSNOW_GRID_Z)
                    end do
                endif
                do k = nlay_new + 1, kSNOW_GRID_Z
                    call zero_snowpack_layer_slot(                     &
                        Tsn_3d, Tsni_3d,                               &
                        VFI_3d, VFW_3d, VFA_3d, VFS_3d, VFWP_3d,       &
                        Rg_3d, Rb_3d, Dd_3d, Sp_3d,                    &
                        mk_3d, mh_3d, CDot_3d, sns_3d, N3_3d, dep_3d,  &
                        i_dom, j_dom, k, ims, ime, jms, jme, kSNOW_GRID_Z)
                end do
            end associate
        else
            ! FSM: Sice/Sliq are extensive (mass per unit area) and must
            ! be rescaled for the partial top layer; snow_temperature
            ! is intensive (plain copy).
            associate(Sice_3d => domain%vars_3d(vidx_Sice)%data_3d, &
                        Sliq_3d => domain%vars_3d(vidx_Sliq)%data_3d, &
                        Tsn_3d  => domain%vars_3d(vidx_Tsn )%data_3d)
                if (nlay_new > 0) then
                    Sice_3d(i_dom, 1, j_dom) = Sice_3d(i_dom, 1+k_eaten, j_dom) * scale
                    Sliq_3d(i_dom, 1, j_dom) = Sliq_3d(i_dom, 1+k_eaten, j_dom) * scale
                    Tsn_3d (i_dom, 1, j_dom) = Tsn_3d (i_dom, 1+k_eaten, j_dom)
                    do k = 2, nlay_new
                        Sice_3d(i_dom, k, j_dom) = Sice_3d(i_dom, k+k_eaten, j_dom)
                        Sliq_3d(i_dom, k, j_dom) = Sliq_3d(i_dom, k+k_eaten, j_dom)
                        Tsn_3d (i_dom, k, j_dom) = Tsn_3d (i_dom, k+k_eaten, j_dom)
                    end do
                endif
                do k = nlay_new + 1, kSNOW_GRID_Z
                    Sice_3d(i_dom, k, j_dom) = 0.0
                    Sliq_3d(i_dom, k, j_dom) = 0.0
                    Tsn_3d (i_dom, k, j_dom) = 0.0
                end do
            end associate
        end if

        domain%vars_2d(vidx_nlay)%data_2di(i_dom, j_dom) = nlay_new

    end subroutine ablate_snow_layers



#ifdef FSM
    !>----------------------------------------------------------
    !! Shift column arrays up by one starting at k_from (destructive: caller
    !! must zero/decrement Nsnow afterward). Used by the FSM relayering steps.
    !!----------------------------------------------------------
    subroutine shift_column_up(Ds, Sice, Sliq, U, k_from, Nsnow)
        !$acc routine seq
        implicit none
        real, intent(inout) :: Ds(:), Sice(:), Sliq(:), U(:)
        integer, intent(in) :: k_from, Nsnow
        integer :: k
        do k = k_from, Nsnow - 1
            Ds(k) = Ds(k+1); Sice(k) = Sice(k+1)
            Sliq(k) = Sliq(k+1); U(k) = U(k+1)
        end do
    end subroutine shift_column_up


    !>----------------------------------------------------------
    !! FSM deposit reconciliation: SNOLAY=1 density-dependent layering
    !! Reads current (eroded) domain layer state, adds deposits from
    !! SD_0/Sice_0 as new surface layers, then relayers.
    !! Matches SNOW_LAYERING.F90 (Quéno et al. 2024) steps 1-7.
    !! Erosion is already synced to domain before this is called.
    !!----------------------------------------------------------
    subroutine snowslide_reconcile_fsm(domain, SD_0, Sice_0)
        implicit none
        type(domain_t), intent(inout) :: domain
        real, intent(in)    :: SD_0(Nx_local, Ny_local)
        real, intent(in)    :: Sice_0(Nx_local, Ny_local)

        ! FSM constants
        real, parameter :: hcap_ice = 2100.0, hcap_wat = 4180.0, Tm = 273.15
        real, parameter :: hfsn = 0.1   ! Snowcover fraction depth scale (m), FSM default

        integer :: i, j, k, i_dom, j_dom, Nsnow_loc
        integer :: kmin_idx, kmerge, kup, kdown, kmax_idx, k_surflay
        integer :: vidx_nlay, vidx_fsnow, vidx_Ds, vidx_Sice, vidx_Sliq, vidx_Tsn
        integer :: vidx_height, vidx_swe
        real :: deposit_depth, deposit_swe, csnow_k, fsnow_loc, fsnow_old, fsnow_new
        real :: snowdepth_loc, Dtemp_surflay, Ds_excess, frac
        real :: Ds_old, Sice_old, Sliq_old, U_old
        real :: rho_kmin, rho_kup, rho_kdown
        real :: Ds_loc(kSNOW_GRID_Z+1), Sice_loc(kSNOW_GRID_Z+1), Sliq_loc(kSNOW_GRID_Z+1), U_loc(kSNOW_GRID_Z+1)
        real :: Tsnow_loc(kSNOW_GRID_Z+1), rho_k(kSNOW_GRID_Z+1), diff_rho(kSNOW_GRID_Z)

        vidx_nlay   = domain%var_indx(kVARS%snow_nlayers)%v
        vidx_fsnow  = domain%var_indx(kVARS%fsnow)%v
        vidx_Ds     = domain%var_indx(kVARS%Ds)%v
        vidx_Sice   = domain%var_indx(kVARS%Sice)%v
        vidx_Sliq   = domain%var_indx(kVARS%Sliq)%v
        vidx_Tsn    = domain%var_indx(kVARS%snow_temperature)%v
        vidx_height = domain%var_indx(kVARS%snow_height)%v
        vidx_swe    = domain%var_indx(kVARS%snow_water_equivalent)%v

        associate( &
            Ds_3d     => domain%vars_3d(vidx_Ds   )%data_3d, &
            Sice_3d   => domain%vars_3d(vidx_Sice )%data_3d, &
            Sliq_3d   => domain%vars_3d(vidx_Sliq )%data_3d, &
            Tsn_3d    => domain%vars_3d(vidx_Tsn  )%data_3d, &
            nlay_2di  => domain%vars_2d(vidx_nlay )%data_2di, &
            fsnow_2d  => domain%vars_2d(vidx_fsnow)%data_2d,  &
            height_2d => domain%vars_2d(vidx_height)%data_2d, &
            swe_2d    => domain%vars_2d(vidx_swe  )%data_2d)

        !$acc parallel loop gang collapse(2) default(present) &
        !$acc   private(i_dom, j_dom, k, Nsnow_loc, kmin_idx, kmerge, kup, kdown, kmax_idx, k_surflay) &
        !$acc   private(deposit_depth, deposit_swe, csnow_k, fsnow_loc, fsnow_old, fsnow_new) &
        !$acc   private(snowdepth_loc, Dtemp_surflay, Ds_excess, frac) &
        !$acc   private(Ds_old, Sice_old, Sliq_old, U_old, rho_kmin, rho_kup, rho_kdown) &
        !$acc   private(Ds_loc, Sice_loc, Sliq_loc, U_loc, Tsnow_loc, rho_k, diff_rho) firstprivate(kSNOW_GRID_Z)
        do j = 2, Ny_local - 1
            do i = 2, Nx_local - 1
                i_dom = its + (i - 2)
                j_dom = jts + (j - 2)

                deposit_depth = SD_0(i,j)
                deposit_swe   = Sice_0(i,j)

                ! Load current (post-erosion) domain state into column arrays
                Nsnow_loc = nlay_2di(i_dom, j_dom)
                fsnow_loc = fsnow_2d(i_dom, j_dom)
                !$acc loop seq
                do k = 1, kSNOW_GRID_Z+1
                    Ds_loc(k) = 0.0; Sice_loc(k) = 0.0; Sliq_loc(k) = 0.0; U_loc(k) = 0.0
                    Tsnow_loc(k) = Tm
                end do

                !$acc loop seq
                do k = 1, Nsnow_loc
                    Ds_loc(k)    = Ds_3d(i_dom, k, j_dom)
                    Sice_loc(k)  = Sice_3d(i_dom, k, j_dom)
                    Sliq_loc(k)  = Sliq_3d(i_dom, k, j_dom)
                    Tsnow_loc(k) = Tsn_3d(i_dom, k, j_dom)
                    csnow_k      = (Sice_loc(k)*hcap_ice + Sliq_loc(k)*hcap_wat) / max(fsnow_loc, 1.0e-10)
                    U_loc(k)     = csnow_k * (Tsnow_loc(k) - Tm)
                end do

                ! --- Step 1: Add deposit as new surface layer (SNOLAY=1) ---
                if (deposit_depth > 1.0e-6) then
                    Nsnow_loc = Nsnow_loc + 1
                    !$acc loop seq
                    do k = Nsnow_loc, 2, -1
                        Ds_loc(k)   = Ds_loc(k-1)
                        Sice_loc(k) = Sice_loc(k-1)
                        Sliq_loc(k) = Sliq_loc(k-1)
                        U_loc(k)    = U_loc(k-1)
                    end do
                    Ds_loc(1)    = deposit_depth
                    Sice_loc(1)  = deposit_swe
                    Sliq_loc(1)  = 0.0
                    Tsnow_loc(1) = Tm    ! Fresh deposit at freezing point ⇒ U = 0
                    U_loc(1)     = 0.0
                end if

                ! --- Step 2: Update fsnow and rescale Ds (SNFRAC=4: tanh model) ---
                snowdepth_loc = sum(Ds_loc(1:Nsnow_loc))

                fsnow_old = fsnow_2d(i_dom, j_dom)

                if (snowdepth_loc > epsilon(snowdepth_loc)) then
                    fsnow_new = min(1.0, tanh(snowdepth_loc / hfsn))
                else
                    fsnow_new = 0.0
                end if

                if (fsnow_new > epsilon(fsnow_new) .and. fsnow_old > epsilon(fsnow_old)) then
                    !$acc loop seq
                    do k = 1, Nsnow_loc
                        Ds_loc(k) = Ds_loc(k) * fsnow_old / fsnow_new
                    end do
                end if

                fsnow_2d(i_dom, j_dom) = fsnow_new

                ! --- Step 3: Restrict surface fine layering ---
                if (Nsnow_loc > 1) then
                    Dtemp_surflay = 0.0
                    k_surflay = 0
                    !$acc loop seq
                    do k = 1, Nsnow_loc - 1
                        Dtemp_surflay = Dtemp_surflay + Ds_loc(k)
                        if (Dtemp_surflay > Ds_surflay) then
                            k_surflay = k
                            exit
                        end if
                    end do
                    if (k_surflay > 0) then
                        Ds_excess = Dtemp_surflay - Ds_surflay
                        Ds_old   = Ds_loc(k_surflay)
                        Sice_old = Sice_loc(k_surflay)
                        Sliq_old = Sliq_loc(k_surflay)
                        U_old    = U_loc(k_surflay)
                        Ds_loc(k_surflay)   = Ds_old - Ds_excess
                        frac                = Ds_loc(k_surflay) / max(Ds_old, 1.0e-10)
                        Sice_loc(k_surflay) = Sice_old * frac
                        Sliq_loc(k_surflay) = Sliq_old * frac
                        U_loc(k_surflay)    = U_old    * frac
                        ! Merge all layers below k_surflay into one — explicit
                        ! reduction loops over temp accumulators instead of
                        ! whole-array `sum(...)`. Accumulators are loaded with
                        ! the remainder from the excess/partial split first,
                        ! then the actual layer contents are added in seq.
                        snowdepth_loc = Ds_excess           ! reuse scalar as tmp
                        Dtemp_surflay = Sice_old * (1.0 - frac)
                        Ds_old   = Sliq_old * (1.0 - frac)  ! reuse Ds_old as Sliq tmp
                        Sice_old = U_old    * (1.0 - frac)  ! reuse Sice_old as U tmp
                        !$acc loop seq
                        do k = k_surflay + 1, Nsnow_loc
                            snowdepth_loc = snowdepth_loc + Ds_loc(k)
                            Dtemp_surflay = Dtemp_surflay + Sice_loc(k)
                            Ds_old        = Ds_old        + Sliq_loc(k)
                            Sice_old      = Sice_old      + U_loc(k)
                        end do
                        Ds_loc(k_surflay+1)   = snowdepth_loc
                        Sice_loc(k_surflay+1) = Dtemp_surflay
                        Sliq_loc(k_surflay+1) = Ds_old
                        U_loc(k_surflay+1)    = Sice_old
                        !$acc loop seq
                        do k = k_surflay + 2, Nsnow_loc
                            Ds_loc(k) = 0.0; Sice_loc(k) = 0.0; Sliq_loc(k) = 0.0; U_loc(k) = 0.0
                        end do
                        Nsnow_loc = k_surflay + 1
                    end if
                end if

                ! --- Step 4: Merge thin layers (SNOLAY=1) ---
                do while (Nsnow_loc > 1)
                    kmin_idx = 1
                    !$acc loop seq
                    do k = 2, Nsnow_loc
                        if (Ds_loc(k) < Ds_loc(kmin_idx)) kmin_idx = k
                    end do
                    if (Ds_loc(kmin_idx) >= Ds_min_fsm) exit

                    if (kmin_idx == 1) then
                        Ds_loc(1)   = Ds_loc(1)   + Ds_loc(2)
                        Sice_loc(1) = Sice_loc(1) + Sice_loc(2)
                        Sliq_loc(1) = Sliq_loc(1) + Sliq_loc(2)
                        U_loc(1)    = U_loc(1)    + U_loc(2)
                        call shift_column_up(Ds_loc, Sice_loc, Sliq_loc, U_loc, 2, Nsnow_loc)
                    else if (kmin_idx == Nsnow_loc) then
                        Ds_loc(Nsnow_loc-1)   = Ds_loc(Nsnow_loc-1)   + Ds_loc(Nsnow_loc)
                        Sice_loc(Nsnow_loc-1) = Sice_loc(Nsnow_loc-1) + Sice_loc(Nsnow_loc)
                        Sliq_loc(Nsnow_loc-1) = Sliq_loc(Nsnow_loc-1) + Sliq_loc(Nsnow_loc)
                        U_loc(Nsnow_loc-1)    = U_loc(Nsnow_loc-1)    + U_loc(Nsnow_loc)
                    else
                        ! Merge with closest-density neighbor
                        kup = kmin_idx - 1; kdown = kmin_idx + 1
                        rho_kup   = (Sice_loc(kup)+Sliq_loc(kup))     / max(Ds_loc(kup), 1.0e-10)
                        rho_kdown = (Sice_loc(kdown)+Sliq_loc(kdown)) / max(Ds_loc(kdown), 1.0e-10)
                        rho_kmin  = (Sice_loc(kmin_idx)+Sliq_loc(kmin_idx)) / max(Ds_loc(kmin_idx), 1.0e-10)
                        if (abs(rho_kmin - rho_kup) < abs(rho_kmin - rho_kdown)) then
                            Ds_loc(kup)   = Ds_loc(kup)   + Ds_loc(kmin_idx)
                            Sice_loc(kup) = Sice_loc(kup) + Sice_loc(kmin_idx)
                            Sliq_loc(kup) = Sliq_loc(kup) + Sliq_loc(kmin_idx)
                            U_loc(kup)    = U_loc(kup)    + U_loc(kmin_idx)
                            call shift_column_up(Ds_loc, Sice_loc, Sliq_loc, U_loc, kmin_idx, Nsnow_loc)
                        else
                            Ds_loc(kmin_idx)   = Ds_loc(kmin_idx)   + Ds_loc(kdown)
                            Sice_loc(kmin_idx) = Sice_loc(kmin_idx) + Sice_loc(kdown)
                            Sliq_loc(kmin_idx) = Sliq_loc(kmin_idx) + Sliq_loc(kdown)
                            U_loc(kmin_idx)    = U_loc(kmin_idx)    + U_loc(kdown)
                            call shift_column_up(Ds_loc, Sice_loc, Sliq_loc, U_loc, kdown, Nsnow_loc)
                        end if
                    end if
                    Ds_loc(Nsnow_loc) = 0.0; Sice_loc(Nsnow_loc) = 0.0
                    Sliq_loc(Nsnow_loc) = 0.0; U_loc(Nsnow_loc) = 0.0
                    Nsnow_loc = Nsnow_loc - 1
                end do

                ! --- Step 5: Merge excess layers (SNOLAY=1) ---
                do while (Nsnow_loc > kSNOW_GRID_Z)
                    !$acc loop seq
                    do k = 1, Nsnow_loc
                        rho_k(k) = (Sice_loc(k)+Sliq_loc(k)) / max(Ds_loc(k), 1.0e-10)
                    end do
                    !$acc loop seq
                    do k = 1, Nsnow_loc - 1
                        diff_rho(k) = abs(rho_k(k) - rho_k(k+1))
                    end do
                    kmerge = 1
                    !$acc loop seq
                    do k = 2, Nsnow_loc - 1
                        if (diff_rho(k) < diff_rho(kmerge)) kmerge = k
                    end do
                    Ds_loc(kmerge)   = Ds_loc(kmerge)   + Ds_loc(kmerge+1)
                    Sice_loc(kmerge) = Sice_loc(kmerge) + Sice_loc(kmerge+1)
                    Sliq_loc(kmerge) = Sliq_loc(kmerge) + Sliq_loc(kmerge+1)
                    U_loc(kmerge)    = U_loc(kmerge)    + U_loc(kmerge+1)
                    call shift_column_up(Ds_loc, Sice_loc, Sliq_loc, U_loc, kmerge+1, Nsnow_loc)
                    Ds_loc(Nsnow_loc) = 0.0; Sice_loc(Nsnow_loc) = 0.0
                    Sliq_loc(Nsnow_loc) = 0.0; U_loc(Nsnow_loc) = 0.0
                    Nsnow_loc = Nsnow_loc - 1
                end do

                ! --- Step 6: Split thick layers (SNOLAY=1) ---
                do while (Nsnow_loc > 0 .and. Nsnow_loc < kSNOW_GRID_Z)
                    kmax_idx = 1
                    !$acc loop seq
                    do k = 2, Nsnow_loc
                        if (Ds_loc(k) > Ds_loc(kmax_idx)) kmax_idx = k
                    end do
                    if (Ds_loc(kmax_idx) < 2.0 * Ds_min_fsm) exit
                    Nsnow_loc = Nsnow_loc + 1
                    !$acc loop seq
                    do k = Nsnow_loc, kmax_idx + 2, -1
                        Ds_loc(k) = Ds_loc(k-1); Sice_loc(k) = Sice_loc(k-1)
                        Sliq_loc(k) = Sliq_loc(k-1); U_loc(k) = U_loc(k-1)
                    end do
                    Ds_loc(kmax_idx+1)   = Ds_loc(kmax_idx)   / 2.0
                    Ds_loc(kmax_idx)     = Ds_loc(kmax_idx)   / 2.0
                    Sice_loc(kmax_idx+1) = Sice_loc(kmax_idx) / 2.0
                    Sice_loc(kmax_idx)   = Sice_loc(kmax_idx) / 2.0
                    Sliq_loc(kmax_idx+1) = Sliq_loc(kmax_idx) / 2.0
                    Sliq_loc(kmax_idx)   = Sliq_loc(kmax_idx) / 2.0
                    U_loc(kmax_idx+1)    = U_loc(kmax_idx)    / 2.0
                    U_loc(kmax_idx)      = U_loc(kmax_idx)    / 2.0
                end do

                ! --- Step 7: Diagnose temperatures from internal energy ---
                !$acc loop seq
                do k = 1, Nsnow_loc
                    csnow_k = (Sice_loc(k)*hcap_ice + Sliq_loc(k)*hcap_wat)
                    if (csnow_k > 1.0e-10) then
                        Tsnow_loc(k) = Tm + U_loc(k) / csnow_k
                    else
                        Tsnow_loc(k) = Tm
                    end if
                end do

                ! --- Write layers back to domain + zero above nlay ---
                nlay_2di(i_dom, j_dom) = Nsnow_loc
                !$acc loop seq
                do k = 1, Nsnow_loc
                    Ds_3d  (i_dom, k, j_dom) = Ds_loc(k)
                    Sice_3d(i_dom, k, j_dom) = Sice_loc(k)
                    Sliq_3d(i_dom, k, j_dom) = Sliq_loc(k)
                    Tsn_3d (i_dom, k, j_dom) = Tsnow_loc(k)
                end do
                !$acc loop seq
                do k = Nsnow_loc + 1, kSNOW_GRID_Z
                    Ds_3d  (i_dom, k, j_dom) = 0.0
                    Sice_3d(i_dom, k, j_dom) = 0.0
                    Sliq_3d(i_dom, k, j_dom) = 0.0
                    Tsn_3d (i_dom, k, j_dom) = 0.0
                end do

                ! Recompute snow_height = sum(Ds) * fsnow, and SWE from layers.
                ! Explicit seq accumulation (no whole-array sum()) for gang safety.
                snowdepth_loc = 0.0
                Dtemp_surflay = 0.0                    ! reused as SWE accumulator
                !$acc loop seq
                do k = 1, Nsnow_loc
                    snowdepth_loc = snowdepth_loc + Ds_3d  (i_dom, k, j_dom)
                    Dtemp_surflay = Dtemp_surflay + Sice_3d(i_dom, k, j_dom) &
                                                  + Sliq_3d(i_dom, k, j_dom)
                end do
                height_2d(i_dom, j_dom) = snowdepth_loc * fsnow_2d(i_dom, j_dom)
                swe_2d   (i_dom, j_dom) = Dtemp_surflay
            end do
        end do
        !$acc end parallel loop

        end associate

    end subroutine snowslide_reconcile_fsm
#endif

    !>----------------------------------------------------------
    !! Copy one SNOWPACK layer slot (k_src → k_dst) at cell (i_dom, j_dom),
    !! across the 17 SNOWPACK per-layer variables OTHER than Ds. Replaces
    !! the dynamic `do iv = 1, n_snow_3d` copy pattern with an explicit
    !! unrolled write so the body is reachable from inside `!$acc parallel`
    !! regions. Callers are responsible for Ds separately — Ds handling
    !! often needs a per-call correction (e.g. partial erosion writes
    !! ds_new instead of shifting), so keeping it out of this helper
    !! avoids a double-write.
    !!----------------------------------------------------------
    subroutine shift_snowpack_layer_slot(Tsn_3d, Tsni_3d,                            &
                                         VFI_3d, VFW_3d, VFA_3d, VFS_3d, VFWP_3d,    &
                                         Rg_3d, Rb_3d, Dd_3d, Sp_3d,                 &
                                         mk_3d, mh_3d, CDot_3d, sns_3d, N3_3d, dep_3d, &
                                         i_dom, j_dom, k_dst, k_src, ims, ime, jms, jme, k_max)
        !$acc routine seq
        implicit none
        ! Fully explicit bounds — matches HICAR's domain memory allocation
        ! (ims:ime, 1:k_max, jms:jme). Caller passes i_dom/j_dom in
        ! original domain coordinates; the dummy preserves that indexing
        ! directly. Fortran does not allow mixing explicit bounds with
        ! assumed-shape `:` in the same array declaration.
        real, intent(inout) :: Tsn_3d(ims:ime, 1:k_max, jms:jme), Tsni_3d(ims:ime, 1:k_max+1, jms:jme)
        real, intent(inout) :: VFI_3d(ims:ime, 1:k_max, jms:jme), VFW_3d(ims:ime, 1:k_max, jms:jme)
        real, intent(inout) :: VFA_3d(ims:ime, 1:k_max, jms:jme), VFS_3d(ims:ime, 1:k_max, jms:jme)
        real, intent(inout) :: VFWP_3d(ims:ime, 1:k_max, jms:jme)
        real, intent(inout) :: Rg_3d(ims:ime, 1:k_max, jms:jme), Rb_3d(ims:ime, 1:k_max, jms:jme)
        real, intent(inout) :: Dd_3d(ims:ime, 1:k_max, jms:jme), Sp_3d(ims:ime, 1:k_max, jms:jme)
        real, intent(inout) :: mk_3d(ims:ime, 1:k_max, jms:jme), mh_3d(ims:ime, 1:k_max, jms:jme)
        real, intent(inout) :: CDot_3d(ims:ime, 1:k_max, jms:jme)
        real, intent(inout) :: sns_3d(ims:ime, 1:k_max, jms:jme), N3_3d(ims:ime, 1:k_max, jms:jme)
        real, intent(inout) :: dep_3d(ims:ime, 1:k_max, jms:jme)
        integer, value, intent(in) :: i_dom, j_dom, k_dst, k_src, ims, ime, jms, jme, k_max

        Tsn_3d(i_dom, k_dst, j_dom)  = Tsn_3d(i_dom, k_src, j_dom)
        Tsni_3d(i_dom, k_dst, j_dom) = Tsni_3d(i_dom, k_src, j_dom)
        VFI_3d(i_dom, k_dst, j_dom)  = VFI_3d(i_dom, k_src, j_dom)
        VFW_3d(i_dom, k_dst, j_dom)  = VFW_3d(i_dom, k_src, j_dom)
        VFA_3d(i_dom, k_dst, j_dom)  = VFA_3d(i_dom, k_src, j_dom)
        VFS_3d(i_dom, k_dst, j_dom)  = VFS_3d(i_dom, k_src, j_dom)
        VFWP_3d(i_dom, k_dst, j_dom) = VFWP_3d(i_dom, k_src, j_dom)
        Rg_3d(i_dom, k_dst, j_dom)   = Rg_3d(i_dom, k_src, j_dom)
        Rb_3d(i_dom, k_dst, j_dom)   = Rb_3d(i_dom, k_src, j_dom)
        Dd_3d(i_dom, k_dst, j_dom)   = Dd_3d(i_dom, k_src, j_dom)
        Sp_3d(i_dom, k_dst, j_dom)   = Sp_3d(i_dom, k_src, j_dom)
        mk_3d(i_dom, k_dst, j_dom)   = mk_3d(i_dom, k_src, j_dom)
        mh_3d(i_dom, k_dst, j_dom)   = mh_3d(i_dom, k_src, j_dom)
        CDot_3d(i_dom, k_dst, j_dom) = CDot_3d(i_dom, k_src, j_dom)
        sns_3d(i_dom, k_dst, j_dom)  = sns_3d(i_dom, k_src, j_dom)
        N3_3d(i_dom, k_dst, j_dom)   = N3_3d(i_dom, k_src, j_dom)
        dep_3d(i_dom, k_dst, j_dom)  = dep_3d(i_dom, k_src, j_dom)
    end subroutine shift_snowpack_layer_slot


    !>----------------------------------------------------------
    !! Zero one SNOWPACK layer slot at cell (i_dom, j_dom) across the 17
    !! SNOWPACK per-layer variables OTHER than Ds. Used for "zero slots
    !! above nlay" cleanup. Callers handle Ds separately.
    !!----------------------------------------------------------
    subroutine zero_snowpack_layer_slot(Tsn_3d, Tsni_3d,                            &
                                        VFI_3d, VFW_3d, VFA_3d, VFS_3d, VFWP_3d,    &
                                        Rg_3d, Rb_3d, Dd_3d, Sp_3d,                 &
                                        mk_3d, mh_3d, CDot_3d, sns_3d, N3_3d, dep_3d, &
                                        i_dom, j_dom, k, ims, ime, jms, jme, k_max)
        !$acc routine seq
        implicit none
        real, intent(inout) :: Tsn_3d(ims:ime, 1:k_max, jms:jme), Tsni_3d(ims:ime, 1:k_max+1, jms:jme)
        real, intent(inout) :: VFI_3d(ims:ime, 1:k_max, jms:jme), VFW_3d(ims:ime, 1:k_max, jms:jme)
        real, intent(inout) :: VFA_3d(ims:ime, 1:k_max, jms:jme), VFS_3d(ims:ime, 1:k_max, jms:jme)
        real, intent(inout) :: VFWP_3d(ims:ime, 1:k_max, jms:jme)
        real, intent(inout) :: Rg_3d(ims:ime, 1:k_max, jms:jme), Rb_3d(ims:ime, 1:k_max, jms:jme)
        real, intent(inout) :: Dd_3d(ims:ime, 1:k_max, jms:jme), Sp_3d(ims:ime, 1:k_max, jms:jme)
        real, intent(inout) :: mk_3d(ims:ime, 1:k_max, jms:jme), mh_3d(ims:ime, 1:k_max, jms:jme)
        real, intent(inout) :: CDot_3d(ims:ime, 1:k_max, jms:jme)
        real, intent(inout) :: sns_3d(ims:ime, 1:k_max, jms:jme), N3_3d(ims:ime, 1:k_max, jms:jme)
        real, intent(inout) :: dep_3d(ims:ime, 1:k_max, jms:jme)
        integer, value, intent(in) :: i_dom, j_dom, k, ims, ime, jms, jme, k_max

        Tsn_3d(i_dom, k, j_dom)  = 0.0
        Tsni_3d(i_dom, k, j_dom) = 0.0
        VFI_3d(i_dom, k, j_dom)  = 0.0
        VFW_3d(i_dom, k, j_dom)  = 0.0
        VFA_3d(i_dom, k, j_dom)  = 0.0
        VFS_3d(i_dom, k, j_dom)  = 0.0
        VFWP_3d(i_dom, k, j_dom) = 0.0
        Rg_3d(i_dom, k, j_dom)   = 0.0
        Rb_3d(i_dom, k, j_dom)   = 0.0
        Dd_3d(i_dom, k, j_dom)   = 0.0
        Sp_3d(i_dom, k, j_dom)   = 0.0
        mk_3d(i_dom, k, j_dom)   = 0.0
        mh_3d(i_dom, k, j_dom)   = 0.0
        CDot_3d(i_dom, k, j_dom) = 0.0
        sns_3d(i_dom, k, j_dom)  = 0.0
        N3_3d(i_dom, k, j_dom)   = 0.0
        dep_3d(i_dom, k, j_dom)  = 0.0
    end subroutine zero_snowpack_layer_slot




    !>----------------------------------------------------------
    !! SNOWPACK deposit reconciliation: insert a fresh surface element for
    !! every cell that received an avalanche deposit. Erosion is already
    !! applied by ablate_snow_layers during routing, so this routine only
    !! handles the deposition side.
    !!
    !! Max-element edge case: if a column is already at kmax elements, the
    !! two deepest elements are merged first (SNOWPACK's combineElements
    !! philosophy: protect the top, coarsen the bottom) before the shift +
    !! insert runs. This guarantees the fresh deposit always lands at k=1.
    !!----------------------------------------------------------
#ifdef SNOWPACK
    subroutine snowslide_reconcile_snowpack(domain, SD_0, Sice_0)
        implicit none
        type(domain_t), intent(inout) :: domain
        real, intent(in)    :: SD_0(Nx_local, Ny_local)
        real, intent(in)    :: Sice_0(Nx_local, Ny_local)

        integer :: i, j, k, i_dom, j_dom, nlay_new
        integer :: vidx_nlay, vidx_Ds, vidx_Tsn, vidx_Tsni
        integer :: vidx_VFI, vidx_VFW, vidx_VFA, vidx_VFS, vidx_VFWP
        integer :: vidx_Rg, vidx_Rb, vidx_Dd, vidx_Sp
        integer :: vidx_mk, vidx_mh, vidx_CDot, vidx_sns, vidx_N3, vidx_dep
        integer :: vidx_skin, vidx_height, vidx_swe
        real :: deposit_depth, theta_i_new, t_surf, n3_fresh, julian_now
        ! Scratch variables for inlined merge_snowpack_bottom_elements
        integer :: k_up, k_lo
        real    :: L_u, L_l, L_new, theta_I_u, theta_I_l, wi_u, wi_l, air_frac

        ! SNOWPACK fresh-element coordination number for rho_deposit = 300 kg/m³,
        ! evaluated from the piecewise polynomial in
        ! build/external/SNOWPACK/fortran/snowpack_driver.F90:1751-1760. Computed
        ! once outside the i,j loop because rho_deposit is a module parameter.
        n3_fresh = 1.4153 - 7.5580e-5 * rho_deposit &
                         + 5.1495e-5 * rho_deposit**2 &
                         - 1.7345e-7 * rho_deposit**3 &
                         + 1.8082e-10 * rho_deposit**4

        julian_now = real(domain%sim_time%mjd())

        vidx_nlay   = domain%var_indx(kVARS%snow_nlayers)%v
        vidx_Ds     = domain%var_indx(kVARS%Ds)%v
        vidx_Tsn    = domain%var_indx(kVARS%snow_temperature)%v
        vidx_Tsni   = domain%var_indx(kVARS%snow_temperature_i)%v
        vidx_VFI    = domain%var_indx(kVARS%Vol_Frac_I)%v
        vidx_VFW    = domain%var_indx(kVARS%Vol_Frac_W)%v
        vidx_VFA    = domain%var_indx(kVARS%Vol_Frac_A)%v
        vidx_VFS    = domain%var_indx(kVARS%Vol_Frac_S)%v
        vidx_VFWP   = domain%var_indx(kVARS%Vol_Frac_WP)%v
        vidx_Rg     = domain%var_indx(kVARS%Rg)%v
        vidx_Rb     = domain%var_indx(kVARS%Rb)%v
        vidx_Dd     = domain%var_indx(kVARS%Dd)%v
        vidx_Sp     = domain%var_indx(kVARS%Sp)%v
        vidx_mk     = domain%var_indx(kVARS%mk)%v
        vidx_mh     = domain%var_indx(kVARS%mass_hoar)%v
        vidx_CDot   = domain%var_indx(kVARS%CDot)%v
        vidx_sns    = domain%var_indx(kVARS%snow_stress)%v
        vidx_N3     = domain%var_indx(kVARS%N3)%v
        vidx_dep    = domain%var_indx(kVARS%depositionDate)%v
        vidx_skin   = domain%var_indx(kVARS%skin_temperature)%v
        vidx_height = domain%var_indx(kVARS%snow_height)%v
        vidx_swe    = domain%var_indx(kVARS%snow_water_equivalent)%v

        associate( &
            Ds_3d       => domain%vars_3d(vidx_Ds  )%data_3d, &
            Tsn_3d      => domain%vars_3d(vidx_Tsn )%data_3d, &
            Tsni_3d     => domain%vars_3d(vidx_Tsni)%data_3d, &
            VFI_3d      => domain%vars_3d(vidx_VFI )%data_3d, &
            VFW_3d      => domain%vars_3d(vidx_VFW )%data_3d, &
            VFA_3d      => domain%vars_3d(vidx_VFA )%data_3d, &
            VFS_3d      => domain%vars_3d(vidx_VFS )%data_3d, &
            VFWP_3d     => domain%vars_3d(vidx_VFWP)%data_3d, &
            Rg_3d       => domain%vars_3d(vidx_Rg  )%data_3d, &
            Rb_3d       => domain%vars_3d(vidx_Rb  )%data_3d, &
            Dd_3d       => domain%vars_3d(vidx_Dd  )%data_3d, &
            Sp_3d       => domain%vars_3d(vidx_Sp  )%data_3d, &
            mk_3d       => domain%vars_3d(vidx_mk  )%data_3d, &
            mh_3d       => domain%vars_3d(vidx_mh  )%data_3d, &
            CDot_3d     => domain%vars_3d(vidx_CDot)%data_3d, &
            sns_3d      => domain%vars_3d(vidx_sns )%data_3d, &
            N3_3d       => domain%vars_3d(vidx_N3  )%data_3d, &
            dep_3d      => domain%vars_3d(vidx_dep )%data_3d, &
            nlay_2di    => domain%vars_2d(vidx_nlay  )%data_2di, &
            skin_2d     => domain%vars_2d(vidx_skin  )%data_2d,  &
            height_2d   => domain%vars_2d(vidx_height)%data_2d,  &
            swe_2d      => domain%vars_2d(vidx_swe   )%data_2d,  &
            ims         => domain%ims, &
            ime         => domain%ime, &
            jms         => domain%jms, &
            jme         => domain%jme)

        !$acc parallel loop gang collapse(2) default(present) &
        !$acc   firstprivate(kSNOW_GRID_Z, ims, ime, jms, jme, n3_fresh, julian_now) &
        !$acc   private(i_dom, j_dom, k, nlay_new, deposit_depth, theta_i_new, t_surf, &
        !$acc           k_up, k_lo, L_u, L_l, L_new, theta_I_u, theta_I_l, wi_u, wi_l, air_frac)
        do j = 2, Ny_local - 1
            do i = 2, Nx_local - 1
                i_dom = its + (i - 2)
                j_dom = jts + (j - 2)
                deposit_depth = SD_0(i,j)
                nlay_new = nlay_2di(i_dom, j_dom)

                if (deposit_depth > 1.0e-6) then

                    ! Edge case: column is full — merge the two deepest elements
                    ! to free a slot, then fall through to the regular shift+insert.
                    ! Inlined from merge_snowpack_bottom_elements: passing arrays
                    ! through `!$acc routine seq` boundaries corrupts device
                    ! pointers (compute-sanitizer shows OOB at the call site).
                    if (nlay_new >= kSNOW_GRID_Z) then
                        k_lo = kSNOW_GRID_Z
                        k_up = k_lo - 1
                        L_u  = Ds_3d(i_dom, k_up, j_dom)
                        L_l  = Ds_3d(i_dom, k_lo, j_dom)
                        L_new = L_u + L_l
                        if (L_new > 0.0) then
                            ! Ice content captured before blending (microstructure weights)
                            theta_I_u = VFI_3d(i_dom, k_up, j_dom)
                            theta_I_l = VFI_3d(i_dom, k_lo, j_dom)
                            wi_u = theta_I_u * L_u
                            wi_l = theta_I_l * L_l

                            ! L-weighted volume fractions + temperatures + diagnostics
                            VFI_3d (i_dom, k_up, j_dom) = (L_u * VFI_3d (i_dom, k_up, j_dom) + L_l * VFI_3d (i_dom, k_lo, j_dom)) / L_new
                            VFW_3d (i_dom, k_up, j_dom) = (L_u * VFW_3d (i_dom, k_up, j_dom) + L_l * VFW_3d (i_dom, k_lo, j_dom)) / L_new
                            VFWP_3d(i_dom, k_up, j_dom) = (L_u * VFWP_3d(i_dom, k_up, j_dom) + L_l * VFWP_3d(i_dom, k_lo, j_dom)) / L_new
                            Tsn_3d (i_dom, k_up, j_dom) = (L_u * Tsn_3d (i_dom, k_up, j_dom) + L_l * Tsn_3d (i_dom, k_lo, j_dom)) / L_new
                            Tsni_3d(i_dom, k_up, j_dom) = (L_u * Tsni_3d(i_dom, k_up, j_dom) + L_l * Tsni_3d(i_dom, k_lo, j_dom)) / L_new
                            sns_3d (i_dom, k_up, j_dom) = (L_u * sns_3d (i_dom, k_up, j_dom) + L_l * sns_3d (i_dom, k_lo, j_dom)) / L_new
                            N3_3d  (i_dom, k_up, j_dom) = (L_u * N3_3d  (i_dom, k_up, j_dom) + L_l * N3_3d  (i_dom, k_lo, j_dom)) / L_new

                            ! Ice-mass-weighted microstructure (no-op if both layers ice-free)
                            if (wi_u + wi_l > 0.0) then
                                Dd_3d  (i_dom, k_up, j_dom) = (wi_u * Dd_3d  (i_dom, k_up, j_dom) + wi_l * Dd_3d  (i_dom, k_lo, j_dom)) / (wi_u + wi_l)
                                Sp_3d  (i_dom, k_up, j_dom) = (wi_u * Sp_3d  (i_dom, k_up, j_dom) + wi_l * Sp_3d  (i_dom, k_lo, j_dom)) / (wi_u + wi_l)
                                Rg_3d  (i_dom, k_up, j_dom) = (wi_u * Rg_3d  (i_dom, k_up, j_dom) + wi_l * Rg_3d  (i_dom, k_lo, j_dom)) / (wi_u + wi_l)
                                Rb_3d  (i_dom, k_up, j_dom) = (wi_u * Rb_3d  (i_dom, k_up, j_dom) + wi_l * Rb_3d  (i_dom, k_lo, j_dom)) / (wi_u + wi_l)
                                CDot_3d(i_dom, k_up, j_dom) = (wi_u * CDot_3d(i_dom, k_up, j_dom) + wi_l * CDot_3d(i_dom, k_lo, j_dom)) / (wi_u + wi_l)
                            end if

                            ! Derived air fraction
                            air_frac = 1.0 - VFI_3d(i_dom, k_up, j_dom) &
                                           - VFW_3d(i_dom, k_up, j_dom) &
                                           - VFWP_3d(i_dom, k_up, j_dom) &
                                           - VFS_3d(i_dom, k_up, j_dom)
                            VFA_3d(i_dom, k_up, j_dom) = max(0.0, air_frac)

                            ! Additive: mass_hoar
                            mh_3d(i_dom, k_up, j_dom) = mh_3d(i_dom, k_up, j_dom) + mh_3d(i_dom, k_lo, j_dom)

                            ! Ds: combined length (assigned LAST)
                            Ds_3d(i_dom, k_up, j_dom) = L_new
                            ! mk and dep keep k_up value (no action needed)
                        end if
                        nlay_new = kSNOW_GRID_Z - 1
                    end if

                    ! Shift existing elements down by 1 to make room at the top.
                    ! Inlined from shift_snowpack_layer_slot for the same reason.
                    if (nlay_new > 0) then
                        !$acc loop seq
                        do k = nlay_new + 1, 2, -1
                            Ds_3d  (i_dom, k, j_dom) = Ds_3d  (i_dom, k-1, j_dom)
                            Tsn_3d (i_dom, k, j_dom) = Tsn_3d (i_dom, k-1, j_dom)
                            Tsni_3d(i_dom, k, j_dom) = Tsni_3d(i_dom, k-1, j_dom)
                            VFI_3d (i_dom, k, j_dom) = VFI_3d (i_dom, k-1, j_dom)
                            VFW_3d (i_dom, k, j_dom) = VFW_3d (i_dom, k-1, j_dom)
                            VFA_3d (i_dom, k, j_dom) = VFA_3d (i_dom, k-1, j_dom)
                            VFS_3d (i_dom, k, j_dom) = VFS_3d (i_dom, k-1, j_dom)
                            VFWP_3d(i_dom, k, j_dom) = VFWP_3d(i_dom, k-1, j_dom)
                            Rg_3d  (i_dom, k, j_dom) = Rg_3d  (i_dom, k-1, j_dom)
                            Rb_3d  (i_dom, k, j_dom) = Rb_3d  (i_dom, k-1, j_dom)
                            Dd_3d  (i_dom, k, j_dom) = Dd_3d  (i_dom, k-1, j_dom)
                            Sp_3d  (i_dom, k, j_dom) = Sp_3d  (i_dom, k-1, j_dom)
                            mk_3d  (i_dom, k, j_dom) = mk_3d  (i_dom, k-1, j_dom)
                            mh_3d  (i_dom, k, j_dom) = mh_3d  (i_dom, k-1, j_dom)
                            CDot_3d(i_dom, k, j_dom) = CDot_3d(i_dom, k-1, j_dom)
                            sns_3d (i_dom, k, j_dom) = sns_3d (i_dom, k-1, j_dom)
                            N3_3d  (i_dom, k, j_dom) = N3_3d  (i_dom, k-1, j_dom)
                            dep_3d (i_dom, k, j_dom) = dep_3d (i_dom, k-1, j_dom)
                        end do
                    end if
                    nlay_new = nlay_new + 1

                    ! Set fresh surface element — defaults from table (inlined
                    ! from set_fresh_snowpack_element_defaults), then per-cell
                    ! patches. Ds is patched below (always per-cell).
                    Tsn_3d (i_dom, 1, j_dom) = snow_3d_fresh_defaults(4)   ! snow_temperature (PATCH)
                    VFI_3d (i_dom, 1, j_dom) = snow_3d_fresh_defaults(5)   ! Vol_Frac_I (PATCH)
                    VFW_3d (i_dom, 1, j_dom) = snow_3d_fresh_defaults(6)   ! Vol_Frac_W
                    VFA_3d (i_dom, 1, j_dom) = snow_3d_fresh_defaults(7)   ! Vol_Frac_A (PATCH)
                    VFS_3d (i_dom, 1, j_dom) = snow_3d_fresh_defaults(8)   ! Vol_Frac_S
                    VFWP_3d(i_dom, 1, j_dom) = snow_3d_fresh_defaults(9)   ! Vol_Frac_WP
                    Tsni_3d(i_dom, 1, j_dom) = snow_3d_fresh_defaults(10)  ! snow_temperature_i (PATCH)
                    Rg_3d  (i_dom, 1, j_dom) = snow_3d_fresh_defaults(11)  ! Rg
                    Rb_3d  (i_dom, 1, j_dom) = snow_3d_fresh_defaults(12)  ! Rb
                    Dd_3d  (i_dom, 1, j_dom) = snow_3d_fresh_defaults(13)  ! Dd
                    Sp_3d  (i_dom, 1, j_dom) = snow_3d_fresh_defaults(14)  ! Sp
                    mk_3d  (i_dom, 1, j_dom) = snow_3d_fresh_defaults(15)  ! mk
                    mh_3d  (i_dom, 1, j_dom) = snow_3d_fresh_defaults(16)  ! mass_hoar
                    CDot_3d(i_dom, 1, j_dom) = snow_3d_fresh_defaults(17)  ! CDot
                    sns_3d (i_dom, 1, j_dom) = snow_3d_fresh_defaults(18)  ! snow_stress
                    N3_3d  (i_dom, 1, j_dom) = snow_3d_fresh_defaults(19)  ! N3 (PATCH, overwritten below)
                    dep_3d (i_dom, 1, j_dom) = snow_3d_fresh_defaults(20)  ! depositionDate (PATCH, overwritten below)

                    theta_i_new           = rho_deposit / rho_ice_c
                    Ds_3d (i_dom, 1, j_dom) = deposit_depth
                    VFI_3d(i_dom, 1, j_dom) = theta_i_new
                    VFA_3d(i_dom, 1, j_dom) = max(0.0, 1.0 - theta_i_new)
                    dep_3d(i_dom, 1, j_dom) = julian_now

                    ! Temperature: clamp to surface skin temperature the way
                    ! snowpack_driver.F90:670 does for new snowfall. A fresh
                    ! element set to Tm collapses the visc_temp term in
                    ! comp_snow_viscosity_default and produces an oversized
                    ! first-call creep.
                    t_surf = min(273.15, skin_2d(i_dom, j_dom))
                    Tsn_3d (i_dom, 1, j_dom) = t_surf
                    Tsni_3d(i_dom, 1, j_dom) = t_surf

                    ! Coordination number: overwrite the table default (1.75)
                    ! with the density-based piecewise from snowpack_driver.F90.
                    ! The metamorphism step will re-derive this each call, but
                    ! the creep block runs BEFORE metamorphism, so the first
                    ! call must already have the right value.
                    N3_3d(i_dom, 1, j_dom) = n3_fresh

                    nlay_2di(i_dom, j_dom) = nlay_new
                end if

                ! Zero all 3D vars above nlay_new (defensive; ablate should
                ! already have). Inlined from zero_snowpack_layer_slot.
                !$acc loop seq
                do k = nlay_new + 1, kSNOW_GRID_Z
                    Ds_3d  (i_dom, k, j_dom) = 0.0
                    Tsn_3d (i_dom, k, j_dom) = 0.0
                    Tsni_3d(i_dom, k, j_dom) = 0.0
                    VFI_3d (i_dom, k, j_dom) = 0.0
                    VFW_3d (i_dom, k, j_dom) = 0.0
                    VFA_3d (i_dom, k, j_dom) = 0.0
                    VFS_3d (i_dom, k, j_dom) = 0.0
                    VFWP_3d(i_dom, k, j_dom) = 0.0
                    Rg_3d  (i_dom, k, j_dom) = 0.0
                    Rb_3d  (i_dom, k, j_dom) = 0.0
                    Dd_3d  (i_dom, k, j_dom) = 0.0
                    Sp_3d  (i_dom, k, j_dom) = 0.0
                    mk_3d  (i_dom, k, j_dom) = 0.0
                    mh_3d  (i_dom, k, j_dom) = 0.0
                    CDot_3d(i_dom, k, j_dom) = 0.0
                    sns_3d (i_dom, k, j_dom) = 0.0
                    N3_3d  (i_dom, k, j_dom) = 0.0
                    dep_3d (i_dom, k, j_dom) = 0.0
                end do

                ! Recompute snow_height and SWE from layers. Inlined as a
                ! seq loop so the gang-parallel region captures no reductions.
                height_2d(i_dom, j_dom) = 0.0
                swe_2d(i_dom, j_dom) = 0.0
                !$acc loop seq
                do k = 1, nlay_new
                    height_2d(i_dom, j_dom) = height_2d(i_dom, j_dom) + Ds_3d(i_dom, k, j_dom)
                    swe_2d(i_dom, j_dom)    = swe_2d(i_dom, j_dom)    &
                        + (VFI_3d(i_dom, k, j_dom) * rho_ice_c        &
                         + VFW_3d(i_dom, k, j_dom) * rho_water_c)     &
                        * Ds_3d(i_dom, k, j_dom)
                end do
            end do
        end do
        !$acc end parallel loop

        end associate

    end subroutine snowslide_reconcile_snowpack
#endif


    !>----------------------------------------------------------
    !! MPI halo exchange for snowslide deposit arrays
    !!----------------------------------------------------------
    subroutine exch_slide_halo(domain, SD_0, Sice_0)
        implicit none
        type(domain_t), intent(inout) :: domain
        real, dimension(Nx_local, Ny_local), intent(inout) :: SD_0, Sice_0

        ! Copy interior to variable_t for halo exchange
        SD_0_var%data_2d(its:ite, jts:jte) = SD_0(2:Nx_local-1, 2:Ny_local-1)
        Sice_0_var%data_2d(its:ite, jts:jte) = Sice_0(2:Nx_local-1, 2:Ny_local-1)

        !$acc update device(SD_0_var%data_2d, Sice_0_var%data_2d)

        call domain%halo%exch_var(SD_0_var, corners=.True.)
        call domain%halo%exch_var(Sice_0_var, corners=.True.)

        !$acc update host(SD_0_var%data_2d, Sice_0_var%data_2d)

        ! Copy back full domain INCLUDING halo to local arrays
        ! its-1:ite+1 includes the halo cells populated by MPI exchange
        SD_0(1:Nx_local, 1:Ny_local) = SD_0_var%data_2d(its-1:ite+1, jts-1:jte+1)
        Sice_0(1:Nx_local, 1:Ny_local) = Sice_0_var%data_2d(its-1:ite+1, jts-1:jte+1)

    end subroutine exch_slide_halo


    !>----------------------------------------------------------
    !! Sort grid points by decreasing DEM elevation (quicksort)
    !!----------------------------------------------------------
    subroutine sort_dem_local()
        implicit none

        integer :: i, j, k, ntotal
        real, allocatable :: dem_vals(:)
        integer, allocatable :: idx_i(:), idx_j(:)

        ntotal = Nx_local * Ny_local

        allocate(dem_vals(ntotal))
        allocate(idx_i(ntotal))
        allocate(idx_j(ntotal))

        ! Flatten 2D DEM into 1D array
        k = 0
        do j = 1, Ny_local
            do i = 1, Nx_local
                k = k + 1
                dem_vals(k) = dem_local(i,j)
                idx_i(k) = i
                idx_j(k) = j
            end do
        end do

        ! Sort by descending elevation
        call quicksort_desc(dem_vals, idx_i, idx_j, 1, ntotal)

        ! Copy results to module index array
        do k = 1, ntotal
            index_sorted(k, 1) = idx_i(k)
            index_sorted(k, 2) = idx_j(k)
        end do

        deallocate(dem_vals, idx_i, idx_j)

    end subroutine sort_dem_local


    !>----------------------------------------------------------
    !! Iterative quicksort in descending order
    !!----------------------------------------------------------
    recursive subroutine quicksort_desc(vals, ii, jj, lo_in, hi_in)
        implicit none
        real, intent(inout)    :: vals(:)
        integer, intent(inout) :: ii(:), jj(:)
        integer, intent(in)    :: lo_in, hi_in

        integer :: iq, jq, lo, hi, ti, sp
        real    :: pivot, tv
        integer :: stack_lo(64), stack_hi(64)

        if (lo_in >= hi_in) return

        sp = 1
        stack_lo(1) = lo_in
        stack_hi(1) = hi_in

        do while (sp > 0)
            lo = stack_lo(sp)
            hi = stack_hi(sp)
            sp = sp - 1

            pivot = vals((lo + hi) / 2)
            iq = lo
            jq = hi

            do while (iq <= jq)
                do while (vals(iq) > pivot)
                    iq = iq + 1
                end do
                do while (vals(jq) < pivot)
                    jq = jq - 1
                end do
                if (iq <= jq) then
                    tv = vals(iq); vals(iq) = vals(jq); vals(jq) = tv
                    ti = ii(iq); ii(iq) = ii(jq); ii(jq) = ti
                    ti = jj(iq); jj(iq) = jj(jq); jj(jq) = ti
                    iq = iq + 1
                    jq = jq - 1
                end if
            end do

            if (jq - lo > hi - iq) then
                if (lo < jq) then
                    sp = sp + 1; stack_lo(sp) = lo; stack_hi(sp) = jq
                end if
                if (iq < hi) then
                    sp = sp + 1; stack_lo(sp) = iq; stack_hi(sp) = hi
                end if
            else
                if (iq < hi) then
                    sp = sp + 1; stack_lo(sp) = iq; stack_hi(sp) = hi
                end if
                if (lo < jq) then
                    sp = sp + 1; stack_lo(sp) = lo; stack_hi(sp) = jq
                end if
            end if
        end do

    end subroutine quicksort_desc

end module module_snowslide
