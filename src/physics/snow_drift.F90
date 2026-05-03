!>----------------------------------------------------------
!! CRYOWRF-style blowing snow redistribution module
!!
!! Implements wind-driven snow transport including:
!!   - Saltation (Liston et al. 1998)
!!   - Near-surface suspension on a fine vertical mesh
!!   - Sublimation of blowing snow (Thorpe & Mason 1966)
!!   - Coupling to 3D atmospheric tracers
!!   - Atmospheric feedback (moisture + temperature)
!!
!! Based on: Sharma et al. (2023) GMD, CRYOWRF
!!
!! @author
!! Dylan Reynolds (dylan.reynolds@epfl.ch)
!!
!!----------------------------------------------------------
module snow_drift
    use domain_interface,   only : domain_t
    use options_interface,  only : options_t
    use variable_interface, only : variable_t
    use mod_wrf_constants,  only : KARMAN, gravity, cp, R_d, XLV, XLS
    use wind,               only : calc_divergence
    use icar_constants
    use data_structures, only : index_type
    use mod_atm_utilities, only : relative_humidity
    implicit none
    private
    public :: snow_drift_var_request, snow_drift_init, snow_drift_step, snow_drift_apply_feedback

    ! -----------------------------------------------
    ! Physical parameters
    ! -----------------------------------------------
    real, parameter :: PI_CONST     = 3.14159265358979323846
    real, parameter :: LOWEST_LAYER       = 0.5        ! Lowest layer height (m)
    real, parameter :: BS_QS_MIN    = 0.0!1.0E-12     ! Floor for mass mixing ratio (kg/kg)
    real, parameter :: BS_NS_MIN    = 0.0!1.0E-6      ! Floor for number concentration (#/kg)
    real, parameter :: RHO_ICE      = 917.0       ! Density of ice (kg/m^3)
    real, parameter :: D_MEAN_SALT  = 200.0E-6    ! Mean saltation particle diameter (m)
    real, parameter :: SIGMA_SALT   = 1.5         ! Log-normal std dev for saltation size dist
    real, parameter :: KA_AIR       = 0.024       ! Thermal conductivity of air (W/m/K)
    real, parameter :: DV_AIR       = 2.2E-5      ! Diffusivity of water vapor in air (m^2/s)
    real, parameter :: NU_AIR       = 1.5E-5      ! Kinematic viscosity of air (m^2/s)
    real, parameter :: RV           = 461.5       ! Gas constant for water vapor (J/kg/K)
    real, parameter :: SETTLE_MIN   = 0.001       ! Minimum settling velocity (m/s)
    real, parameter :: AT_COEFF     = 0.02        ! Geometrical coeff for weight term (Lehning et al. 2000)
    real, parameter :: BT_COEFF     = 0.0015      ! Geometrical coeff for binding term (Lehning et al. 2000)
    real, parameter :: SIGMA_REF    = 300.0       ! Reference shear strength (Pa) (Lehning et al. 2000)
    real, parameter :: VQ_RATIO     = 2.303       ! Mass-weighted Stokes ratio (gamma alpha=3)
    real, parameter :: VN_RATIO     = 0.822       ! Number-weighted Stokes ratio (gamma alpha=3)
    real, parameter :: alpha        = 3.0         ! Exponent for gamma distribution of particle sizes (Sharma et al. 2023)

    ! Mitchell (1996) Re-X regime boundaries and coefficients
    real, parameter :: MITCHELL_X1 = 10.0       ! Stokes / intermediate boundary
    real, parameter :: MITCHELL_X2 = 585.0      ! Intermediate / turbulent boundary
    real, parameter :: MITCHELL_AM_S  = 0.04167 ! a_m Stokes (1/24)
    real, parameter :: MITCHELL_BM_S  = 1.0     ! b_m Stokes
    real, parameter :: MITCHELL_AM_I  = 0.3285  ! a_m intermediate
    real, parameter :: MITCHELL_BM_I  = 0.5842  ! b_m intermediate
    real, parameter :: MITCHELL_AM_T  = 0.8388  ! a_m turbulent
    real, parameter :: MITCHELL_BM_T  = 0.5044  ! b_m turbulent

    real, parameter :: LAMBDA_COEFF = 3.91487   ! 60^(1/3) for lambda = LAMBDA_COEFF / D_m

    ! Precomputed gamma ratios for Eq. 13 sublimation (Sharma et al. 2023)
    ! Gamma(1+1.5*b_m+alpha)/Gamma(alpha) for each Mitchell regime (computed in init)
    real, save :: GAMMA_RATIO_S, GAMMA_RATIO_I, GAMMA_RATIO_T
    real, save :: SC_THIRD   ! Sc^(1/3) = (NU_AIR/DV_AIR)^(1/3)

    ! -----------------------------------------------
    ! Module-level fine mesh arrays (allocated in snow_drift_init)
    ! -----------------------------------------------
    integer, save :: snc_N_loc                                  ! Local copy of fine mesh levels
    real, allocatable, save :: snc_Z(:,:,:)                     ! (ims:ime, snc_N, jms:jme) heights AGL
    real, allocatable, save :: snc_dz(:,:,:)                        ! (snc_N) layer thickness
    real, allocatable, save :: u_fm(:,:,:)                      ! (ims:ime, snc_N, jms:jme) U-wind at fine mesh
    real, allocatable, save :: v_fm(:,:,:)                      ! (ims:ime, snc_N, jms:jme) V-wind at fine mesh
    real, allocatable, save :: w_fm(:,:,:)                      ! (ims:ime, snc_N, jms:jme) W-wind at fine mesh
    real, allocatable, save :: wind_fm(:,:,:)                   ! (ims:ime, snc_N, jms:jme) wind speed at fine mesh
    real, allocatable, save :: t_fm(:,:,:)                      ! (ims:ime, snc_N, jms:jme) temperature at fine mesh
    real, allocatable, save :: qv_fm(:,:,:)                     ! (ims:ime, snc_N, jms:jme) water vapor mixing ratio at fine mesh
    real, allocatable, save :: Kh_fm(:,:,:)                     ! (ims:ime, snc_N, jms:jme) eddy diffusivity
    real, allocatable, save :: subl_mass_2d(:,:)                ! (ims:ime, jms:jme) accumulated sublimation mass
    real, allocatable, save :: qs_flux_cache(:,:)            ! (ims:ime, jms:jme) cached interface flux tendency to qs
    real, allocatable, save :: ns_flux_cache(:,:)            ! (ims:ime, jms:jme) cached interface flux tendency to ns

    ! Halo exchange variable indices for batch exchange
    type(index_type), save :: exch_vars_fm(2)

    ! Domain index copies
    integer, save :: ids, ide, jds, jde, kds, kde
    integer, save :: ims, ime, jms, jme, kms, kme
    integer, save :: its, ite, jts, jte, kts, kte

    logical, save :: module_initialized = .false.

contains

    !>----------------------------------------------------------
    !! Request allocation of variables needed by snow drift
    !!----------------------------------------------------------
    subroutine snow_drift_var_request(options)
        implicit none
        type(options_t), intent(inout) :: options

        ! Request 2D diagnostic variables
        call options%alloc_vars( &
            [kVARS%bs_threshold_ustar,  &
             kVARS%bs_saltation_flux,   &
             kVARS%bs_saltation_height, &
             kVARS%bs_saltation_concentration, &
             kVARS%bs_swe_exchange,    &
             kVARS%bs_sublimation_flux, &
             kVARS%bs_suspension_flux, &
             kVARS%bs_drift_swe_salt,   &
             kVARS%bs_drift_swe_susp,   &
             kVARS%bs_drift_swe_subl,   &
             kVARS%qs_fm, kVARS%ns_fm, &
             kVARS%ustar, kVARS%density, kVARS%dz_interface, kVARS%exner, &
             kVARS%water_vapor, kVARS%potential_temperature, kVARS%temperature, &
             kVARS%u_10m, kVARS%v_10m, kVARS%snow_water_equivalent, &
             kVARS%snow_height, kVARS%skin_temperature, kVARS%land_mask])

        ! Request restart variables for 3D tracers
        call options%restart_vars( &
            [kVARS%qs_fm,    &
             kVARS%ns_fm,  &
             kVARS%bs_drift_swe_salt,   &
             kVARS%bs_drift_swe_susp,   &
             kVARS%bs_drift_swe_subl])

    end subroutine snow_drift_var_request


    !>----------------------------------------------------------
    !! Initialize fine mesh arrays and module state
    !!----------------------------------------------------------
    subroutine snow_drift_init(domain, options)
        implicit none
        type(domain_t),  intent(inout) :: domain
        type(options_t), intent(in)    :: options

        integer :: i, j, k
        real    :: z_top, dlog_z, z_bot, fm_top

        if (module_initialized) then
            !$acc exit data delete(snc_dz,snc_Z,subl_mass_2d,kH_fm,wind_fm,v_fm,u_fm, qv_fm,t_fm,w_fm,qs_flux_cache,ns_flux_cache)
            deallocate(snc_Z)
            deallocate(snc_dz)
            deallocate(u_fm)
            deallocate(v_fm)
            deallocate(w_fm)
            deallocate(wind_fm)
            deallocate(t_fm)
            deallocate(qv_fm)
            deallocate(Kh_fm)
            deallocate(subl_mass_2d)
            deallocate(qs_flux_cache)
            deallocate(ns_flux_cache)
            module_initialized = .False.
        endif

        ! Store domain indices
        ids = domain%ids ; ide = domain%ide ; jds = domain%jds ; jde = domain%jde
        kds = domain%kds ; kde = domain%kde
        ims = domain%ims ; ime = domain%ime ; jms = domain%jms ; jme = domain%jme
        kms = domain%kms ; kme = domain%kme
        its = domain%its ; ite = domain%ite ; jts = domain%jts ; jte = domain%jte
        kts = domain%kts ; kte = domain%kte

        snc_N_loc = options%sm%suspension_fine_mesh_levels
        if (snc_N_loc > domain%kme .and. STD_OUT_PE) then
            write(*,*) "WARNING: number of fine mesh levels specified for suspension scheme greater than number of atmospheric levels: ",kte
            write(*,*) "WARNING: clipping fine mesh levels to: ",kte
        endif
        snc_N_loc = min(snc_N_loc,kte)

        ! Allocate fine mesh arrays
        allocate(snc_Z(ims:ime, snc_N_loc, jms:jme))
        allocate(snc_dz(ims:ime, snc_N_loc, jms:jme))
        allocate(u_fm(ims:ime, snc_N_loc, jms:jme))
        allocate(v_fm(ims:ime, snc_N_loc, jms:jme))
        allocate(w_fm(ims:ime, snc_N_loc, jms:jme))
        allocate(wind_fm(ims:ime, snc_N_loc, jms:jme))
        allocate(t_fm(ims:ime, snc_N_loc, jms:jme))
        allocate(qv_fm(ims:ime, snc_N_loc, jms:jme))
        allocate(Kh_fm(ims:ime, snc_N_loc, jms:jme))
        allocate(subl_mass_2d(ims:ime, jms:jme))
        allocate(qs_flux_cache(ims:ime, jms:jme))
        allocate(ns_flux_cache(ims:ime, jms:jme))

        ! Initialize to zero
        u_fm        = 0.0
        v_fm        = 0.0
        wind_fm     = 0.0
        Kh_fm       = 0.0
        subl_mass_2d = 0.0

        ! Build logarithmic fine mesh from LOWEST_LAYER to approximate first model level height
        ! Heights are AGL; the actual coupling height varies per cell and is handled at runtime
        do j = jms, jme
            do i = ims, ime
                fm_top  = domain%vars_3d(domain%var_indx(kVARS%dz_interface)%v)%data_3d(i,1,j) * 0.5  ! Approximate first model level; adjusted per cell at runtime

                dlog_z = log(fm_top / LOWEST_LAYER) / real(snc_N_loc)
                do k = 1, snc_N_loc
                    snc_Z(i, k, j) = LOWEST_LAYER * exp((real(k) - 0.5) * dlog_z)
                end do
                snc_dz(i,1,j) = LOWEST_LAYER * (exp(dlog_z) - 1.0)
                do k = 2, snc_N_loc
                    snc_dz(i,k,j) = LOWEST_LAYER * (exp(real(k) * dlog_z) - exp(real(k-1) * dlog_z))
                enddo
            enddo
        end do

        !$acc enter data copyin(snc_dz,snc_Z,subl_mass_2d,kH_fm,wind_fm,v_fm,u_fm) create(qv_fm,t_fm,w_fm,qs_flux_cache,ns_flux_cache)
        ! Initialize fine mesh exchange variable indices for batch halo exchange
        exch_vars_fm(1)%v = domain%var_indx(kVARS%qs_fm)%v
        exch_vars_fm(1)%id = domain%vars_3d(exch_vars_fm(1)%v)%id
        exch_vars_fm(2)%v = domain%var_indx(kVARS%ns_fm)%v
        exch_vars_fm(2)%id = domain%vars_3d(exch_vars_fm(2)%v)%id

        ! Precompute gamma ratios for Eq. 13 sublimation (Sharma et al. 2023)
        GAMMA_RATIO_S = gamma(1.0 + 1.5*MITCHELL_BM_S + alpha) / gamma(alpha)
        GAMMA_RATIO_I = gamma(1.0 + 1.5*MITCHELL_BM_I + alpha) / gamma(alpha)
        GAMMA_RATIO_T = gamma(1.0 + 1.5*MITCHELL_BM_T + alpha) / gamma(alpha)
        SC_THIRD = (NU_AIR / DV_AIR)**(1.0/3.0)

        module_initialized = .true.

        !if (STD_OUT_PE) write(*,*) "  Snow drift (CRYOWRF-style) initialized with", snc_N_loc, "fine mesh levels"

    end subroutine snow_drift_init


    !>----------------------------------------------------------
    !! Main blowing snow physics step
    !!----------------------------------------------------------
    subroutine snow_drift_step(domain, options, dt)
        implicit none
        type(domain_t),  intent(inout) :: domain
        type(options_t), intent(in)    :: options
        real,            intent(in)    :: dt

        integer :: i, j
        real    :: dx

        if (options%physics%snowmodel == 0 .or. options%sm%suspension_layer /= 1) return

        dx = domain%dx

        ! 1. Build fine mesh wind profiles
        call build_fine_mesh_winds(domain)

        ! 2. Compute threshold friction velocity from the currently-exposed
        ! snowpack layer. Works for both FSM (generic Li & Pomeroy) and
        ! SNOWPACK (Schmidt/Lehning).
        ! self-limits once erosion exposes a harder layer.
        if (options%physics%snowmodel /= kSM_SNOWPACK) then
            call compute_threshold_ustar_generic(domain)
        else
            call compute_threshold_ustar_snowpack(domain, options%sm%saltation_model)
        endif

        ! saltation mass flux and concentration are computed by snowpack fortran driver
#ifndef SNOWPACK_FORTRAN
        ! 3. Saltation concentration calculation
        call saltation_step(domain, dt, dx)
#endif

        ! 4. Fine mesh suspension (horizontal advection + vertical diffusion/settling/sublimation)
        call snow_drift_integrate(domain, dt, dx, options)

        ! 5. Couple fine mesh top to 3D atmospheric grid
        call couple_to_3d_grid(domain, dt)

    end subroutine snow_drift_step


    !>----------------------------------------------------------
    !! Compute threshold friction velocity from SNOWPACK grain properties.
    !!
    !! Lehning et al. (2000) Eq. 15 / SNOWPACK SnowDrift.cc
    !! (Sharma et al. 2023, GMD — CRYOWRF).
    !!----------------------------------------------------------
    subroutine compute_threshold_ustar_snowpack(domain, salt_model)
        implicit none
        type(domain_t), intent(inout) :: domain
        integer, intent(in) :: salt_model

        integer :: i, j, n_snow, k
        real :: Sp_top, Rg_top, Rb_top, N3_top
        real :: rho_air, rg_m
        real :: weight, binding, tau_thresh, ustar_t, ustar_t_out

        associate( &
            bs_ustar_t  => domain%vars_2d(domain%var_indx(kVARS%bs_threshold_ustar)%v)%data_2d, &
            ustar       => domain%vars_2d(domain%var_indx(kVARS%ustar)%v)%data_2d, &
            swe         => domain%vars_2d(domain%var_indx(kVARS%snow_water_equivalent)%v)%data_2d, &
            land_mask   => domain%vars_2d(domain%var_indx(kVARS%land_mask)%v)%data_2di, &
            density     => domain%vars_3d(domain%var_indx(kVARS%density)%v)%data_3d, &
            Ds_3d       => domain%vars_3d(domain%var_indx(kVARS%Ds)%v)%data_3d, &
            VFI_3d      => domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_I)%v)%data_3d, &
            Sp_3d       => domain%vars_3d(domain%var_indx(kVARS%Sp)%v)%data_3d, &
            Rg_3d       => domain%vars_3d(domain%var_indx(kVARS%Rg)%v)%data_3d, &
            Rb_3d       => domain%vars_3d(domain%var_indx(kVARS%Rb)%v)%data_3d, &
            N3_3d       => domain%vars_3d(domain%var_indx(kVARS%N3)%v)%data_3d, &
            bs_salt_conc => domain%vars_2d(domain%var_indx(kVARS%bs_saltation_concentration)%v)%data_2d, &
            bs_salt_flux => domain%vars_2d(domain%var_indx(kVARS%bs_saltation_flux)%v)%data_2d, &
            bs_salt_height => domain%vars_2d(domain%var_indx(kVARS%bs_saltation_height)%v)%data_2d, &
            bs_swe_exch => domain%vars_2d(domain%var_indx(kVARS%bs_swe_exchange)%v)%data_2d, &
            n_layers    => domain%vars_2d(domain%var_indx(kVARS%snow_nlayers)%v)%data_2di &
        )

        !$acc parallel loop gang vector collapse(2)
        do j = jts, jte
            do i = its, ite
                n_snow = n_layers(i,j)

                if (swe(i,j) < 1.0 .or. n_snow < 1 .or. land_mask(i,j) == 2) then
                    bs_ustar_t(i,j) = 9999.0  ! No snow or water; suppress drift
                    cycle
                endif

                ! Top snow layer properties (layer index 1 = top in SNOWPACK)
                Sp_top  = max(Sp_3d(i,1,j), 0.0)
                Rg_top  = max(Rg_3d(i,1,j), 0.01)   ! min 0.05 mm (stored in mm)
                Rb_top  = max(Rb_3d(i,1,j), 0.0)
                N3_top  = max(N3_3d(i,1,j), 0.0)     ! coordination number from SNOWPACK

                ! Air density at surface
                rho_air = density(i,kts,j)

                ! Grain radius in meters for weight term (Rg stored in mm)
                rg_m = Rg_top * 1.0E-3

                ! Threshold shear stress (Lehning et al. 2000, Eq. 15)
                ! tau = At * rho_ice * (Sp+1) * g * rg  +  Bt * sigma_ref * N3 * (rb/rg)^2
                weight  = AT_COEFF * RHO_ICE * (Sp_top + 1.0) * gravity * rg_m
                binding = BT_COEFF * SIGMA_REF * N3_top * (Rb_top / max(Rg_top, 0.01))**2

                tau_thresh = weight + binding

                ! Threshold friction velocity: u*_t = sqrt(tau / rho_air)
                ustar_t = sqrt(tau_thresh / rho_air)

                ! Clamp to physical range
                bs_ustar_t(i,j) = max(0.05, min(ustar_t, 5.0))

            enddo
        enddo

        if (salt_model==kSALTATION_SORENSEN) then
            !$acc parallel loop gang vector collapse(2)
            do j = jts, jte
                do i = its, ite
                    if (n_layers(i,j) < 1) cycle

                    ! Code copied from snowpack_laws.F90
                    bs_salt_height(i,j) = 0.15 ! hard coded based on Armin's analysis (See Viaro et al., 2026)

                    if (ustar(i,j) > bs_ustar_t(i,j)) then
                        bs_salt_flux(i,j) = 0.0014 * density(i,kts,j) * ustar(i,j) &
                            * (ustar(i,j) - bs_ustar_t(i,j)) &
                            * (ustar(i,j) + 7.6 * bs_ustar_t(i,j) + 205.0)
                        ! Concentration: arbitrary scaling chosen by C++ to match Doorschot
                        ! magnitudes (Saltation.cc:411). Units are nominally kg/m^3 once the
                        ! suspension scheme consumes c_salt / rho_air.
                        if (ustar(i,j) > 0.0) bs_salt_conc(i,j) = bs_salt_flux(i,j) / ustar(i,j) * 0.001
                    else
                        bs_salt_flux(i,j) = 0.0
                        bs_salt_conc(i,j) = 0.0
                    endif
                enddo
            enddo
        endif
        end associate

    end subroutine compute_threshold_ustar_snowpack


    !>----------------------------------------------------------
    !! Compute threshold friction velocity from snow temperature
    !! Li & Pomeroy (1997)
    !!----------------------------------------------------------
    subroutine compute_threshold_ustar_generic(domain)
        implicit none
        type(domain_t), intent(inout) :: domain

        integer :: i, j
        real :: T_snow

        associate( &
            bs_ustar_t  => domain%vars_2d(domain%var_indx(kVARS%bs_threshold_ustar)%v)%data_2d, &
            swe         => domain%vars_2d(domain%var_indx(kVARS%snow_water_equivalent)%v)%data_2d, &
            skin_temp   => domain%vars_2d(domain%var_indx(kVARS%skin_temperature)%v)%data_2d, &
            land_mask   => domain%vars_2d(domain%var_indx(kVARS%land_mask)%v)%data_2di &
        )

        !$acc parallel loop gang vector collapse(2)
        do j = jts, jte
            do i = its, ite
                if (swe(i,j) < 1.0 .or. land_mask(i,j) == 2) then
                    bs_ustar_t(i,j) = 9999.0
                    cycle
                endif

                T_snow = skin_temp(i,j)
                ! Li & Pomeroy (1997): u*_t = 0.15 + 0.00167 * max(0, T - 258.15)^2.5
                bs_ustar_t(i,j) = 0.15 + 0.00167 * max(0.0, T_snow - 258.15)**2.5

                ! Clamp to physical range
                bs_ustar_t(i,j) = max(0.05, min(bs_ustar_t(i,j), 5.0))
            enddo
        enddo

        end associate

    end subroutine compute_threshold_ustar_generic


    !>----------------------------------------------------------
    !! Build Monin-Obukhov wind profiles on the fine mesh
    !!----------------------------------------------------------
    subroutine build_fine_mesh_winds(domain)
        implicit none
        type(domain_t), intent(inout) :: domain

        integer :: i, j, k
        real    :: z0_loc, z_ref, log_ratio

        associate( &
            w   => domain%vars_3d(domain%var_indx(kVARS%w)%v)%data_3d, &
            u        => domain%vars_2d(domain%var_indx(kVARS%u_10m)%v)%data_2d, &
            v        => domain%vars_2d(domain%var_indx(kVARS%v_10m)%v)%data_2d, &
            temperature => domain%vars_3d(domain%var_indx(kVARS%temperature)%v)%data_3d, &
            qv       => domain%vars_3d(domain%var_indx(kVARS%water_vapor)%v)%data_3d, &
            dz       => domain%vars_3d(domain%var_indx(kVARS%dz_interface)%v)%data_3d, &
            z0       => domain%vars_2d(domain%var_indx(kVARS%roughness_z0)%v)%data_2d, &
            Kh       => domain%vars_3d(domain%var_indx(kVARS%coeff_heat_exchange_3d)%v)%data_3d &
        )

        z_ref = 10.0  ! Reference height for 10m wind

        ! Compute over full memory domain (including halo) for advection stencil access
        !$acc parallel loop gang vector collapse(2)
        do j = jms, jme
            do i = ims, ime
                z0_loc    = max(z0(i,j), 1.0E-5)
                !$acc loop
                do k = 1, snc_N_loc
                    ! vertical wind
                    w_fm(i,k,j) = w(i,1,j) * (snc_Z(i,k,j) / (dz(i,1,j) * 0.5))
                    ! Diffusivity coefficient
                    Kh_fm(i,k,j) = Kh(i,2,j) * (snc_Z(i,k,j) / (dz(i,1,j) * 0.5))

                    log_ratio = log(max(snc_Z(i,k,j), z0_loc * 1.1) / z0_loc) &
                                / log(z_ref / z0_loc)

                    ! Monin-Obukhov logarithmic profiles for horizontal wind, temperature, and moisture
                    u_fm(i,k,j) = u(i,j) * log_ratio
                    v_fm(i,k,j) = v(i,j) * log_ratio
                    t_fm(i,k,j) = temperature(i,1,j) * log_ratio
                    qv_fm(i,k,j) = qv(i,1,j) * log_ratio

                    wind_fm(i,k,j) = sqrt(u_fm(i,k,j)**2 + v_fm(i,k,j)**2)

                enddo
            enddo
        enddo

        end associate

    end subroutine build_fine_mesh_winds


    !>----------------------------------------------------------
    !! Saltation step — now a no-op.
    !!
    !! The saltation mass flux, saltation height, and saltation-top concentration
    !! are computed by the SNOWPACK Fortran driver (Doorschot & Lehning 2002) and
    !! written directly into bs_saltation_flux / bs_saltation_height /
    !! bs_saltation_concentration BEFORE snow_drift_step is called (lsm runs at
    !! time_step.F90:392, snow_drift_step at line 426). This subroutine is kept
    !! as a hook in case future logic (halo exchanges, diagnostics) needs to run
    !! between the SNOWPACK write and the suspension integrate.

    !! Saltation step with iterative halo exchange
    !! Sorensen (2004) / Sharma et al. (2023) Eq. 17 + upwind divergence
    !!----------------------------------------------------------
    subroutine saltation_step(domain, dt, dx)
        implicit none
        type(domain_t), intent(inout) :: domain
        real,           intent(in)    :: dt, dx

        integer :: i, j, iter
        real    :: ustar_loc, ustar_t, rho_air, ratio_t
        real    :: flux_u, flux_v, div_salt, swe_change

        associate( &
            bs_ustar_t     => domain%vars_2d(domain%var_indx(kVARS%bs_threshold_ustar)%v)%data_2d, &
            ustar_2d       => domain%vars_2d(domain%var_indx(kVARS%ustar)%v)%data_2d,              &
            density        => domain%vars_3d(domain%var_indx(kVARS%density)%v)%data_3d,             &
            swe            => domain%vars_2d(domain%var_indx(kVARS%snow_water_equivalent)%v)%data_2d, &
            bs_salt_mass   => domain%vars_2d(domain%var_indx(kVARS%bs_saltation_flux)%v)%data_2d,   &
            bs_swe_salt    => domain%vars_2d(domain%var_indx(kVARS%bs_drift_swe_salt)%v)%data_2d    &
        )

        ! Compute saltation mass flux (Qsalt in kg/m/s)
        !$acc parallel loop gang vector collapse(2)
        do j = jts, jte
            do i = its, ite
                ustar_loc = max(ustar_2d(i,j), 0.01)
                ustar_t   = bs_ustar_t(i,j)
                rho_air   = density(i,kts,j)

                if (ustar_loc > ustar_t .and. swe(i,j) > 1.0) then
                    ! This comes from the SNOWPACK code, which is based on Doorschot & Lehning (2002) and Sorensen (1991):
                    bs_salt_mass(i,j) = 0.0014 * rho_air * ustar_loc &
                        * (ustar_loc - ustar_t) &
                        * (ustar_loc + 7.6 * ustar_t + 205.0)
                else
                    bs_salt_mass(i,j) = 0.0
                endif
            enddo
        enddo

        end associate

    end subroutine saltation_step


    !>----------------------------------------------------------
    !! Operator-split blowing snow transport (Sharma et al. 2023, Eqs. 20-22)
    !!
    !! Eq. 20: dq_bs/dt = dq_bs/dt|_(m+s) + dq_bs/dt|_(p) + dq_bs/dt|_(advec)
    !!
    !! Operator A: Horizontal advection (Eq. 21c) via adv_std_advect_horiz
    !! Operator B: Vertical mixing + sedimentation (Eq. 21a) semi-implicit Thomas
    !! Operator C: Phase change / sublimation (Eq. 21b) explicit
    !! Surface mass balance (Eq. 22)
    !!----------------------------------------------------------
    subroutine snow_drift_integrate(domain, dt, dx, options)
        use adv_std,      only: adv_std_compute_wind_2d_fm, flux_2d_fm, sum_kernel_2d_fm, adv_std_clean_wind_arrays_fm, &
                                flux_x_fm, flux_y_fm, flux_z_fm
        use adv_fluxcorr, only: init_fluxcorr_fm, set_sign_arrays_fm, WRF_flux_corr_fm
        implicit none
        type(domain_t), intent(inout) :: domain
        real,           intent(in)    :: dt, dx
        type(options_t), intent(in)   :: options

        integer :: i, j, k, iter
        real    :: rho_air, T_air, qv_air, e_sat, rh, sigma_i
        real    :: D_m, V_coef
        real    :: S_q, S_N, h_salt_loc
        real    :: Kh_half, dep_mass_salt, dep_mass_susp, dep_num_susp
        real    :: dz_below, dz_above
        real    :: q_salt_val, n_salt_val, ghost_coeff
        real    :: dq_entrain, dz_salt_interface
        real    :: q_ghost_limited, n_ghost_limited
        real    :: ghost_coeff_top, dz_top_interface
        real    :: lambda, A_thermo, eta, phi_best, X_best
        real    :: am_loc, bm_loc, gr_loc
        ! CRYOWRF saltation BC variables
        real    :: upart, conc_col, k_salt, c_salt_loc

        ! Subcycled vertical advection variables
        real    :: flux_v(0:snc_N_loc)  ! Vertical interface fluxes
        real    :: w_iface              ! Interpolated interface velocity
        real    :: max_Cr_col, dt_sub   ! Per-column max CFL and subcycle timestep
        integer :: n_sub, iter_v        ! Subcycle count and iterator
        real, parameter :: CFL_MAX_V = 0.9  ! Max vertical Courant number per subcycle

        ! Thomas algorithm arrays (per column)
        real :: a(snc_N_loc), b(snc_N_loc), c(snc_N_loc), d_rhs(snc_N_loc)
        real :: Vq(snc_N_loc), Vn(snc_N_loc)
        real :: Vq_iface_down(0:snc_N_loc), Vq_iface_up(0:snc_N_loc)
        real :: Vn_iface_down(0:snc_N_loc), Vn_iface_up(0:snc_N_loc)
        real :: V_net_q, V_net_n, dt_thomas

        ! Advection work arrays

        real, allocatable :: rho_fm(:,:,:), div(:,:,:)
        real, allocatable :: U_m_fm(:,:,:), V_m_fm(:,:,:), denom_fm(:,:,:)
        real, allocatable :: qs_fm_old(:,:,:), ns_fm_old(:,:,:)
        real, allocatable :: jaco_fm(:,:,:), jaco_u_fm(:,:,:), jaco_v_fm(:,:,:)
        real    :: t_fac
        integer :: flux_corr, max_iters
        logical :: using_snowpack, saltation_doorschot

        associate( &
            density    => domain%vars_3d(domain%var_indx(kVARS%density)%v)%data_3d, &
            jacobian   => domain%vars_3d(domain%var_indx(kVARS%jacobian)%v)%data_3d, &
            jacobian_u => domain%vars_3d(domain%var_indx(kVARS%jacobian_u)%v)%data_3d, &
            jacobian_v => domain%vars_3d(domain%var_indx(kVARS%jacobian_v)%v)%data_3d, &
            swe        => domain%vars_2d(domain%var_indx(kVARS%snow_water_equivalent)%v)%data_2d, &
            skin_temp  => domain%vars_2d(domain%var_indx(kVARS%skin_temperature)%v)%data_2d, &
            bs_swe_exch => domain%vars_2d(domain%var_indx(kVARS%bs_swe_exchange)%v)%data_2d, &
            z0         => domain%vars_2d(domain%var_indx(kVARS%roughness_z0)%v)%data_2d, &
            psfc       => domain%vars_2d(domain%var_indx(kVARS%surface_pressure)%v)%data_2d, &
            ustar_2d   => domain%vars_2d(domain%var_indx(kVARS%ustar)%v)%data_2d, &
            bs_ustar_t => domain%vars_2d(domain%var_indx(kVARS%bs_threshold_ustar)%v)%data_2d, &
            bs_salt_mass => domain%vars_2d(domain%var_indx(kVARS%bs_saltation_flux)%v)%data_2d, &
            bs_salt_height => domain%vars_2d(domain%var_indx(kVARS%bs_saltation_height)%v)%data_2d, &
            bs_salt_conc => domain%vars_2d(domain%var_indx(kVARS%bs_saltation_concentration)%v)%data_2d, &
            bs_susp_flux => domain%vars_2d(domain%var_indx(kVARS%bs_suspension_flux)%v)%data_2d, &
            bs_subl    => domain%vars_2d(domain%var_indx(kVARS%bs_sublimation_flux)%v)%data_2d, &
            bs_swe_salt => domain%vars_2d(domain%var_indx(kVARS%bs_drift_swe_salt)%v)%data_2d, &
            bs_swe_susp => domain%vars_2d(domain%var_indx(kVARS%bs_drift_swe_susp)%v)%data_2d, &
            bs_swe_subl => domain%vars_2d(domain%var_indx(kVARS%bs_drift_swe_subl)%v)%data_2d,  &
            qs_fm => domain%vars_3d(domain%var_indx(kVARS%qs_fm)%v)%data_3d, &
            ns_fm => domain%vars_3d(domain%var_indx(kVARS%ns_fm)%v)%data_3d, &
            qs_atm => domain%vars_3d(domain%var_indx(kVARS%snow_mass)%v)%data_3d, &
            ns_atm => domain%vars_3d(domain%var_indx(kVARS%snow_number)%v)%data_3d, &
            dz_atm => domain%vars_3d(domain%var_indx(kVARS%dz_interface)%v)%data_3d &
        )


        ! ===================================================================
        ! Operator A — 3D advection on fine mesh (Eq. 21c) via RK3
        ! ===================================================================

        ! Build density field on fine mesh (broadcast surface density to all levels)
        allocate(rho_fm(ims:ime, 1:snc_N_loc, jms:jme))
        allocate(qs_fm_old(ims:ime, snc_N_loc, jms:jme))
        allocate(ns_fm_old(ims:ime, snc_N_loc, jms:jme))
        allocate(div(ims:ime, kms:kme, jms:jme))
        allocate(jaco_fm  (ims:ime, 1:snc_N_loc, jms:jme))
        allocate(jaco_u_fm(ims:ime, 1:snc_N_loc, jms:jme))
        allocate(jaco_v_fm(ims:ime, 1:snc_N_loc, jms:jme))

        ! Halo exchange before saving initial state for RK3 advection
        call exch_fine_mesh_3d(domain)

        using_snowpack = (options%physics%snowmodel == kSM_SNOWPACK)
        saltation_doorschot = (options%sm%saltation_model == kSALTATION_DOORSCHOT)
        !$acc data create(rho_fm, qs_fm_old, ns_fm_old, div, jaco_fm, jaco_u_fm, jaco_v_fm)

        ! Zero sublimation accumulator for this step
        !$acc parallel loop gang vector collapse(2) default(present)
        do j = jms, jme
            do i = ims, ime
                subl_mass_2d(i,j) = 0.0
                bs_subl(i,j) = 0.0
            enddo
        enddo
        !$acc parallel loop gang vector collapse(3) default(present)
        do j = jms, jme
            do k = 1, snc_N_loc
                do i = ims, ime
                    rho_fm(i,k,j) = density(i,kts,j)
                    jaco_fm(i,k,j)   = 1.0!jacobian(i,kts,j)
                    jaco_u_fm(i,k,j) = 1.0!jacobian_u(i,kts,j)
                    jaco_v_fm(i,k,j) = 1.0!jacobian_v(i,kts,j)
                    ! Save initial states for RK3
                    qs_fm_old(i,k,j) = qs_fm(i,k,j)
                    ns_fm_old(i,k,j) = ns_fm(i,k,j)
                enddo
            enddo
        enddo

        call calc_divergence(div, domain, horz_only=.True., use_dqdt=.False., advect_density=.False.)

        ! Allocate fine-mesh wind/denom arrays here (not inside adv_std_compute_wind_2d_fm)
        ! and manage their device lifecycle directly on these local variables.
        ! This avoids NVHPC OpenACC issues with !$acc enter/exit data on allocatable
        ! dummy arguments, which can leave stale device mappings across calls and
        ! cause "partially present" errors when nest contexts switch.
        allocate(U_m_fm  (its-2:ite+3, 1:snc_N_loc, jts-2:jte+3))
        allocate(V_m_fm  (its-2:ite+3, 1:snc_N_loc, jts-2:jte+3))
        allocate(denom_fm(ims:ime,     1:snc_N_loc, jms:jme))
        !$acc enter data create(U_m_fm, V_m_fm, denom_fm)

        ! Compute 3D fine-mesh wind Courant numbers (U_m_fm, V_m_fm, denom_fm)
        call adv_std_compute_wind_2d_fm(u_fm, v_fm, rho_fm, &
            jaco_fm, jaco_u_fm, jaco_v_fm, dx, dt, &
            U_m_fm, V_m_fm, denom_fm, &
            ims, ime, 1, snc_N_loc, jms, jme, its, ite, jts, jte)

        if (options%adv%flux_corr == kFLUXCOR_MONO) then
            ! Initialize flux correction arrays and sign arrays (used on 3rd substep)
            call init_fluxcorr_fm(ims, ime, 1, snc_N_loc, jms, jme, its, ite, jts, jte)
            call set_sign_arrays_fm(U_m_fm, V_m_fm, &
                ims, ime, 1, snc_N_loc, jms, jme, its, ite, jts, jte)
        endif

        if (options%time%RK3) then
            max_iters = 3
        else
            max_iters = 1
        endif
        ! ===================================================================
        ! IMEX-RK3: Explicit horizontal advection + implicit vertical transport
        ! Each RK3 substep: horizontal advection (3rd-order) followed by
        ! implicit Thomas solver (diffusion + effective velocity + ghost cell)
        ! ===================================================================

        do iter = 1, max_iters
            if (options%time%RK3) then
                select case(iter)
                case (1)
                    t_fac = 1.0/3.0
                    flux_corr = 0
                case (2)
                    t_fac = 0.5
                    flux_corr = 0
                case (3)
                    t_fac = 1.0
                    flux_corr = kFLUXCOR_MONO
                end select
            else
                t_fac = 1.0
                flux_corr = 0
            endif

            ! Halo exchange before each substep
            if (iter > 1) call exch_fine_mesh_3d(domain)

            ! --- Explicit horizontal advection for q_bs ---
            call flux_2d_fm(qs_fm, U_m_fm, V_m_fm, t_fac, &
                ims, ime, 1, snc_N_loc, jms, jme, its, ite, jts, jte)

            if (flux_corr == kFLUXCOR_MONO) then
                call WRF_flux_corr_fm(qs_fm_old, U_m_fm, V_m_fm, &
                    flux_x_fm, flux_y_fm, denom_fm, &
                    ims, ime, 1, snc_N_loc, jms, jme, its, ite, jts, jte)
            endif

            call sum_kernel_2d_fm(qs_fm_old, qs_fm, denom_fm, &
                ims, ime, 1, snc_N_loc, jms, jme, its, ite, jts, jte)

            ! --- Explicit horizontal advection for N_bs ---
            call flux_2d_fm(ns_fm, U_m_fm, V_m_fm, t_fac, &
                ims, ime, 1, snc_N_loc, jms, jme, its, ite, jts, jte)

            if (flux_corr == kFLUXCOR_MONO) then
                call WRF_flux_corr_fm(ns_fm_old, U_m_fm, V_m_fm, &
                    flux_x_fm, flux_y_fm, denom_fm, &
                    ims, ime, 1, snc_N_loc, jms, jme, its, ite, jts, jte)
            endif

            call sum_kernel_2d_fm(ns_fm_old, ns_fm, denom_fm, &
                ims, ime, 1, snc_N_loc, jms, jme, its, ite, jts, jte)

            ! --- Implicit vertical transport (Thomas solver) within each RK substep ---
            dt_thomas = dt * t_fac

            !$acc parallel loop gang vector collapse(2) private(a, b, c, d_rhs, Vq, Vn, &
            !$acc& Vq_iface_down, Vq_iface_up, Vn_iface_down, Vn_iface_up) default(present)
            do j = jts, jte
                do i = its, ite
                    rho_air = density(i,kts,j)

                    ! --- Saltation ghost-cell values from SNOWPACK Doorschot-Lehning ---
                    ! SNOWPACK writes bs_salt_conc (kg/m^3) and bs_salt_height (m) directly;
                    ! convert concentration to mixing ratio via / rho_air for use as the
                    ! bottom ghost-cell BC of the suspension diffusion equation.
                    q_salt_val = 0.0
                    n_salt_val = 0.0
                    if (using_snowpack .and. saltation_doorschot) then
                        if (ustar_2d(i,j) > bs_ustar_t(i,j) .and. bs_salt_conc(i,j) > 0.0 .and. &
                            swe(i,j) >= 0.5 .and. skin_temp(i,j) < 272.15) then
                            q_salt_val = bs_salt_conc(i,j) / rho_air
                            n_salt_val = 0.45 * q_salt_val &
                                        * (6.0 / (PI_CONST * D_MEAN_SALT**3 * RHO_ICE))
                        endif
                    else
                        bs_salt_height(i,j) = 0.15   ! [m] fixed saltation reference height

                        ! Safety gates (CRYOWRF Coupler.cpp lines 969-975)
                        if (ustar_2d(i,j) > bs_ustar_t(i,j) .and. bs_salt_mass(i,j) > 0.0 .and. swe(i,j) >= 0.5 &
                            .and. skin_temp(i,j) < 272.15) then
                            ! Particle speed: 75% of log-law wind at h_salt (Nishimura 2014)
                            upart = 0.75 * ustar_2d(i,j) / 0.4 * log(bs_salt_height(i,j) / 0.002)
                            if (upart > 0.0) then
                                conc_col   = bs_salt_mass(i,j) / upart  ! [kg/m^2] column mass
                                ! Capped inverse scale height (Melo et al. 2024 Fig. 4)
                                k_salt     = max(0.45 * gravity / max(ustar_2d(i,j)**2, 1.0e-4), &
                                                1.0 / 0.04)
                                bs_salt_conc(i,j) = conc_col * k_salt * exp(-k_salt * bs_salt_height(i,j))  ! [kg/m^3]
                                q_salt_val = bs_salt_conc(i,j) / rho_air  ! [kg/kg]
                                ! Number conc with fixed d=140μm (CRYOWRF Coupler.cpp line 990)
                                n_salt_val = 0.45 * q_salt_val &
                                        * (6.0 / (PI_CONST * 0.00014**3 * RHO_ICE))
                            endif
                        endif
                    endif

                    ! --- Compute settling velocities and effective interface velocities ---
                    !$acc loop
                    do k = 1, snc_N_loc
                        lambda = ((max(ns_fm(i,k,j),1e-10) * RHO_ICE * PI_CONST * gamma(alpha + 3.0)) &
                                / (max(qs_fm(i,k,j), 1e-10) * 6.0 * gamma(alpha)) )**(1.0/3.0)
                        eta = rho_air * NU_AIR
                        phi_best = (4.0/3.0) * rho_air * (RHO_ICE - rho_air) * gravity / (eta * eta)
                        X_best = phi_best * (alpha/lambda)**3

                        if (X_best <= MITCHELL_X1) then
                            am_loc = MITCHELL_AM_S;  bm_loc = MITCHELL_BM_S
                        else if (X_best <= MITCHELL_X2) then
                            am_loc = MITCHELL_AM_I;  bm_loc = MITCHELL_BM_I
                        else
                            am_loc = MITCHELL_AM_T;  bm_loc = MITCHELL_BM_T
                        endif

                        V_coef = am_loc * eta * phi_best**(bm_loc) / rho_air
                        Vq(k) = V_coef * gamma(3*bm_loc + 2 + alpha) / (gamma(alpha) * gamma(alpha + 3)) * lambda**(1-3*bm_loc)
                        Vn(k) = V_coef * gamma(3*bm_loc - 1 + alpha) / (gamma(alpha)) * lambda**(1-3*bm_loc)
                        Vq(k) = max(SETTLE_MIN, min(Vq(k), 2.0))
                        Vn(k) = max(SETTLE_MIN, min(Vn(k), 2.0))
                    enddo

                    ! Effective velocity at each interface: V_net = Vq_avg - w_fm (positive = downward)
                    Vq_iface_down(0) = max(Vq(1), 0.0); Vq_iface_up(0) = 0.0
                    Vn_iface_down(0) = max(Vn(1), 0.0); Vn_iface_up(0) = 0.0

                    !$acc loop seq
                    do k = 1, snc_N_loc - 1
                        w_iface = 0.5 * (w_fm(i,k,j) + w_fm(i,k+1,j))
                        V_net_q = 0.5*(Vq(k) + Vq(k+1)) - w_iface
                        V_net_n = 0.5*(Vn(k) + Vn(k+1)) - w_iface
                        Vq_iface_down(k) = max(V_net_q, 0.0)
                        Vq_iface_up(k)   = max(-V_net_q, 0.0)
                        Vn_iface_down(k) = max(V_net_n, 0.0)
                        Vn_iface_up(k)   = max(-V_net_n, 0.0)
                    enddo

                    Vq_iface_down(snc_N_loc) = 0.0; Vq_iface_up(snc_N_loc) = 0.0
                    Vn_iface_down(snc_N_loc) = 0.0; Vn_iface_up(snc_N_loc) = 0.0

                    ! ============ Thomas solve for q_bs ============
                    a = 0.0; b = 0.0; c = 0.0; d_rhs = 0.0

                    ! k=1: Ghost-cell bottom BC from saltation
                    dz_salt_interface = 0.5 * snc_dz(i,1,j)
                    ghost_coeff = dt_thomas * Kh_fm(i,1,j) / (dz_salt_interface * snc_dz(i,1,j))
                    if (ghost_coeff * q_salt_val * rho_air * snc_dz(i,1,j) > swe(i,j) &
                        .and. q_salt_val > 0.0) then
                        q_ghost_limited = swe(i,j) / (ghost_coeff * rho_air * snc_dz(i,1,j))
                    else
                        q_ghost_limited = q_salt_val
                    endif

                    b(1) = 1.0 + ghost_coeff
                    ! Effective velocity at interface above cell 1
                    b(1) = b(1) + dt_thomas * Vq_iface_down(0) / snc_dz(i,1,j)   ! loss below (ghost)
                    b(1) = b(1) + dt_thomas * Vq_iface_up(1)   / snc_dz(i,1,j)   ! loss above
                    if (snc_N_loc > 1) then
                        dz_above = 0.5 * (snc_dz(i,1,j) + snc_dz(i,2,j))
                        Kh_half  = 0.5 * (Kh_fm(i,2,j) + Kh_fm(i,1,j))
                        c(1)     = -dt_thomas * Kh_half / (dz_above * snc_dz(i,1,j))
                        b(1)     = b(1) - c(1)
                        c(1)     = c(1) - dt_thomas * Vq_iface_down(1) / snc_dz(i,1,j) ! gain from k+1
                    endif
                    d_rhs(1) = qs_fm(i,1,j) + ghost_coeff * q_ghost_limited

                    ! Interior levels
                    !$acc loop
                    do k = 2, snc_N_loc - 1
                        dz_below = 0.5 * (snc_dz(i,k-1,j) + snc_dz(i,k,j))
                        dz_above = 0.5 * (snc_dz(i,k,j)   + snc_dz(i,k+1,j))

                        Kh_half = 0.5 * (Kh_fm(i,k,j) + Kh_fm(i,k-1,j))
                        a(k) = -dt_thomas * Kh_half / (dz_below * snc_dz(i,k,j))

                        Kh_half = 0.5 * (Kh_fm(i,k+1,j) + Kh_fm(i,k,j))
                        c(k) = -dt_thomas * Kh_half / (dz_above * snc_dz(i,k,j))

                        b(k) = 1.0 - a(k) - c(k)

                        ! Effective velocity upwind terms
                        c(k) = c(k) - dt_thomas * Vq_iface_down(k)   / snc_dz(i,k,j)  ! gain from above
                        b(k) = b(k) + dt_thomas * Vq_iface_up(k)     / snc_dz(i,k,j)  ! loss upward
                        b(k) = b(k) + dt_thomas * Vq_iface_down(k-1) / snc_dz(i,k,j)  ! loss downward
                        a(k) = a(k) - dt_thomas * Vq_iface_up(k-1)   / snc_dz(i,k,j)  ! gain from below

                        d_rhs(k) = qs_fm(i,k,j)
                    enddo

                    ! k=N: Top boundary (zero-gradient)
                    if (snc_N_loc > 1) then
                        dz_below = 0.5 * (snc_dz(i,snc_N_loc-1,j) + snc_dz(i,snc_N_loc,j))
                        Kh_half = 0.5 * (Kh_fm(i,snc_N_loc,j) + Kh_fm(i,snc_N_loc-1,j))
                        a(snc_N_loc) = -dt_thomas * Kh_half / (dz_below * snc_dz(i,snc_N_loc,j))
                        c(snc_N_loc) = 0.0
                        b(snc_N_loc) = 1.0 - a(snc_N_loc)
                        b(snc_N_loc) = b(snc_N_loc) + dt_thomas * Vq_iface_down(snc_N_loc-1) / snc_dz(i,snc_N_loc,j)
                        a(snc_N_loc) = a(snc_N_loc) - dt_thomas * Vq_iface_up(snc_N_loc-1) / snc_dz(i,snc_N_loc,j)
                        d_rhs(snc_N_loc) = qs_fm(i,snc_N_loc,j)
                    endif

                    call thomas_solve(snc_N_loc, a, b, c, d_rhs)
                    !$acc loop
                    do k = 1, snc_N_loc
                        qs_fm(i,k,j) = max(BS_QS_MIN, d_rhs(k))
                    enddo

                    ! Mass tracking: ghost-cell entrainment and settling deposition
                    ! (SWE is NOT modified here — the net mass budget is accumulated in
                    ! bs_swe_susp and applied to SNOWPACK elements at the next SNOWPACK step)
                    dq_entrain = ghost_coeff * (q_ghost_limited - qs_fm(i,1,j)) * rho_air * snc_dz(i,1,j)
                    dep_mass_susp = rho_air * Vq_iface_down(0) * qs_fm(i,1,j) * dt_thomas
                    bs_susp_flux(i,j) = (- dq_entrain + dep_mass_susp)/dt

                    ! ============ Thomas solve for N_bs ============
                    a = 0.0; b = 0.0; c = 0.0; d_rhs = 0.0

                    if (q_salt_val > 0.0) then
                        n_ghost_limited = q_ghost_limited * (n_salt_val / q_salt_val)
                    else
                        n_ghost_limited = 0.0
                    endif

                    b(1) = 1.0 + ghost_coeff
                    b(1) = b(1) + dt_thomas * Vn_iface_down(0) / snc_dz(i,1,j)
                    b(1) = b(1) + dt_thomas * Vn_iface_up(1)   / snc_dz(i,1,j)
                    if (snc_N_loc > 1) then
                        dz_above = 0.5 * (snc_dz(i,1,j) + snc_dz(i,2,j))
                        Kh_half  = 0.5 * (Kh_fm(i,2,j) + Kh_fm(i,1,j))
                        c(1)     = -dt_thomas * Kh_half / (dz_above * snc_dz(i,1,j))
                        b(1)     = b(1) - c(1)
                        c(1)     = c(1) - dt_thomas * Vn_iface_down(1) / snc_dz(i,1,j)
                    endif
                    d_rhs(1) = ns_fm(i,1,j) + ghost_coeff * n_ghost_limited

                    !$acc loop
                    do k = 2, snc_N_loc - 1
                        dz_below = 0.5 * (snc_dz(i,k-1,j) + snc_dz(i,k,j))
                        dz_above = 0.5 * (snc_dz(i,k,j)   + snc_dz(i,k+1,j))

                        Kh_half = 0.5 * (Kh_fm(i,k,j) + Kh_fm(i,k-1,j))
                        a(k) = -dt_thomas * Kh_half / (dz_below * snc_dz(i,k,j))

                        Kh_half = 0.5 * (Kh_fm(i,k+1,j) + Kh_fm(i,k,j))
                        c(k) = -dt_thomas * Kh_half / (dz_above * snc_dz(i,k,j))

                        b(k) = 1.0 - a(k) - c(k)

                        c(k) = c(k) - dt_thomas * Vn_iface_down(k)   / snc_dz(i,k,j)
                        b(k) = b(k) + dt_thomas * Vn_iface_up(k)     / snc_dz(i,k,j)
                        b(k) = b(k) + dt_thomas * Vn_iface_down(k-1) / snc_dz(i,k,j)
                        a(k) = a(k) - dt_thomas * Vn_iface_up(k-1)   / snc_dz(i,k,j)

                        d_rhs(k) = ns_fm(i,k,j)
                    enddo

                    if (snc_N_loc > 1) then
                        dz_below = 0.5 * (snc_dz(i,snc_N_loc-1,j) + snc_dz(i,snc_N_loc,j))
                        Kh_half = 0.5 * (Kh_fm(i,snc_N_loc,j) + Kh_fm(i,snc_N_loc-1,j))
                        a(snc_N_loc) = -dt_thomas * Kh_half / (dz_below * snc_dz(i,snc_N_loc,j))
                        c(snc_N_loc) = 0.0
                        b(snc_N_loc) = 1.0 - a(snc_N_loc)
                        b(snc_N_loc) = b(snc_N_loc) + dt_thomas * Vn_iface_down(snc_N_loc-1) / snc_dz(i,snc_N_loc,j)
                        a(snc_N_loc) = a(snc_N_loc) - dt_thomas * Vn_iface_up(snc_N_loc-1) / snc_dz(i,snc_N_loc,j)
                        d_rhs(snc_N_loc) = ns_fm(i,snc_N_loc,j)
                    endif

                    call thomas_solve(snc_N_loc, a, b, c, d_rhs)
                    !$acc loop
                    do k = 1, snc_N_loc
                        ns_fm(i,k,j) = max(BS_NS_MIN, d_rhs(k))
                    enddo

                enddo  ! i
            enddo  ! j

        enddo  ! RK3 iter



        ! ===================================================================
        ! Post-RK3: Sublimation + Surface mass balance
        ! (Thomas solver now runs inside RK3 substeps above)
        ! ===================================================================
        !$acc parallel loop gang vector collapse(2) default(present)
        do j = jts, jte
            do i = its, ite
                rho_air = density(i,kts,j)

                ! ===================================================================
                ! Sublimation (Eq. 21b / Eq. 13, explicit)
                ! Sharma et al. (2023) Eq. 13: integrated over gamma distribution
                ! ===================================================================

                !$acc loop seq
                do k = 1, snc_N_loc
                    T_air  = t_fm(i,k,j)
                    qv_air = qv_fm(i,k,j)

                    ! Saturation vapor pressure over ice (Buck 1981)
                    e_sat = 611.15 * exp(22.452 * (T_air - 273.15) / (T_air - 0.61))
                    rh = relative_humidity(T_air, qv_air, psfc(i,j))

                    ! ! Saturation deficit (negative when subsaturated)
                    sigma_i = rh - 1.0

                    if (qs_fm(i,k,j) > 0.0 .and. ns_fm(i,k,j) > 0.0) then
                        ! Mean volume diameter
                        D_m = max( alpha/( (max(ns_fm(i,k,j), 1.0e-10) * RHO_ICE * PI_CONST * gamma(alpha + 3.0)) / (max(qs_fm(i,k,j), 1.0e-10) * 6.0 * gamma(alpha)))**(1.0/3.0), 20.0E-6)

                        ! Gamma distribution slope parameter (alpha=3)
                        lambda = alpha / D_m

                        ! Thermodynamic coefficient A (Eq. 12a)
                        A_thermo = 2.0 * PI_CONST / &
                                   (XLS/(KA_AIR * T_air) * (XLS/(RV * T_air) - 1.0) + &
                                    RV * T_air / (DV_AIR * e_sat))

                        ! Best number parameter phi (Mitchell 1996)
                        eta = rho_air * NU_AIR
                        phi_best = (4.0/3.0) * rho_air * (RHO_ICE - rho_air) * gravity / (eta * eta)

                        ! Mean Best number for regime selection
                        X_best = phi_best * D_m**3

                        ! Select Mitchell (1996) regime coefficients
                        if (X_best <= MITCHELL_X1) then
                            am_loc = MITCHELL_AM_S;  bm_loc = MITCHELL_BM_S;  gr_loc = GAMMA_RATIO_S
                        else if (X_best <= MITCHELL_X2) then
                            am_loc = MITCHELL_AM_I;  bm_loc = MITCHELL_BM_I;  gr_loc = GAMMA_RATIO_I
                        else
                            am_loc = MITCHELL_AM_T;  bm_loc = MITCHELL_BM_T;  gr_loc = GAMMA_RATIO_T
                        endif

                        ! Eq. 13: S_q integrated over gamma distribution
                        S_q = ns_fm(i,k,j) * (sigma_i * &
                              (0.78 * alpha * A_thermo / lambda + &
                               0.308 * A_thermo * SC_THIRD * sqrt(am_loc) * phi_best**(bm_loc/2.0) * &
                               gr_loc / lambda**(1.0 + 1.5*bm_loc)) )



                        if (S_q <= 0.0) then
                            ! Limit to available mass
                            S_q = max(S_q, (-qs_fm(i,k,j) / dt))
                            ! Eq. 14: S_N (Morrison & Grabowski 2008)
                            S_N = S_q * (ns_fm(i,k,j) / max(qs_fm(i,k,j), 1.0e-10))
                        else
                            S_q = S_q !min(S_q,DBLE((bs_qv_spec(k)-loc_qsat)*0.999/dt_in))
                            ! Limit to available number
                            S_N = 0.0
                        endif
                        S_N = MAX(S_N, -ns_fm(i,k,j)/dt)

                        ! Update fields
                        qs_fm(i,k,j) = max(BS_QS_MIN, qs_fm(i,k,j) + S_q * dt)
                        ns_fm(i,k,j) = max(BS_NS_MIN, ns_fm(i,k,j) + S_N * dt)

                        ! Accumulate sublimation mass (positive = mass lost from snow)
                        subl_mass_2d(i,j) = subl_mass_2d(i,j) - S_q * rho_air * snc_dz(i,k,j) * dt
                    endif
                enddo

                ! ===================================================================
                ! Surface mass balance (Eq. 22)
                ! Suspension entrainment/deposition already applied to SWE above
                ! (in the flux BC and post-Thomas deposition steps)
                ! ===================================================================

                ! Saltation deposition from horizontal divergence
                ! (SWE is NOT modified here — net budget applied via bs_swe_exchange)
                ! divergence is positive for horizontal spreading (net erosion), negative for convergence (net deposition) so
                ! we multiply by -1 to get the correct sign for deposition/erosion
                dep_mass_salt = (bs_salt_mass(i,j) * (-div(i,1,j)) / wind_fm(i,snc_N_loc,j)) * dt
                bs_swe_salt(i,j) = bs_swe_salt(i,j) + dep_mass_salt

                bs_swe_susp(i,j) = bs_swe_susp(i,j) + bs_susp_flux(i,j) * dt

                ! Sublimation diagnostics
                bs_subl(i,j)     = subl_mass_2d(i,j) * XLS / dt  ! W/m^2
                bs_swe_subl(i,j) = bs_swe_subl(i,j) + subl_mass_2d(i,j)

                ! Accumulate net drift mass exchange for SNOWPACK element modification.
                ! Positive = deposition (mass added to snowpack), negative = erosion.
                ! Applied to SNOWPACK elements at the NEXT snowpack_fortran_step call.
                bs_swe_exch(i,j) = bs_swe_exch(i,j) + bs_susp_flux(i,j) * dt + dep_mass_salt

            enddo
        enddo
        !$acc end data


        end associate
        ! Cleanup advection arrays — handle U_m_fm/V_m_fm/denom_fm directly here
        ! (on the local variables) to avoid NVHPC dummy-argument tracking issues.
        !$acc exit data delete(U_m_fm, V_m_fm, denom_fm)
        deallocate(U_m_fm, V_m_fm, denom_fm)
        ! Module-level flux_x_fm/flux_y_fm cleanup is still in adv_std module
        call adv_std_clean_wind_arrays_fm()
        deallocate(rho_fm, qs_fm_old, ns_fm_old, jaco_fm, jaco_u_fm, jaco_v_fm)

    end subroutine snow_drift_integrate


    !>----------------------------------------------------------
    !! Couple fine mesh top to 3D atmospheric grid
    !! Flux continuity at fine mesh top ↔ first atm level
    !!----------------------------------------------------------
    subroutine couple_to_3d_grid(domain, dt)
        implicit none
        type(domain_t), intent(inout) :: domain
        real,           intent(in)    :: dt

        integer :: i, j
        real    :: rho_air, dz_atm, dz_iface, Kh_iface
        real    :: flux_up_q, flux_up_n, flux_down_q, flux_down_n
        real    :: max_flux_q, max_flux_n
        real    :: lambda, eta, phi_best, X_best
        real    :: am_loc, bm_loc, V_coef, Vq_settle, Vn_settle, Vq_down, Vn_down

        associate( &
            density     => domain%vars_3d(domain%var_indx(kVARS%density)%v)%data_3d, &
            dz          => domain%vars_3d(domain%var_indx(kVARS%dz_interface)%v)%data_3d, &
            qs => domain%vars_3d(domain%var_indx(kVARS%snow_mass)%v)%data_3d, &
            ns  => domain%vars_3d(domain%var_indx(kVARS%snow_number)%v)%data_3d, &
            qs_fm => domain%vars_3d(domain%var_indx(kVARS%qs_fm)%v)%data_3d, &
            ns_fm => domain%vars_3d(domain%var_indx(kVARS%ns_fm)%v)%data_3d &

        )

        !$acc parallel loop gang vector collapse(2)
        do j = jts, jte
            do i = its, ite
                rho_air = density(i,kts,j)
                dz_atm  = dz(i,kts,j)

                ! Interface distance: half fine mesh top cell + half atmospheric cell
                dz_iface = 0.5 * snc_dz(i,snc_N_loc,j) + 0.5 * dz_atm

                ! ============================================================
                ! Upward flux: diffusion from fine mesh top → atmospheric kts
                ! ============================================================
                Kh_iface = Kh_fm(i,snc_N_loc,j)

                ! Diffusive flux (positive = upward, fine mesh → atmosphere)
                flux_up_q = Kh_iface * (qs_fm(i,snc_N_loc,j) - qs(i,kts,j)) / dz_iface
                flux_up_n = Kh_iface * (ns_fm(i,snc_N_loc,j) - ns(i,kts,j))  / dz_iface

                flux_up_q = flux_up_q + max(w_fm(i,snc_N_loc,j), 0.0) * qs_fm(i,snc_N_loc,j)
                flux_up_n = flux_up_n + max(w_fm(i,snc_N_loc,j), 0.0) * ns_fm(i,snc_N_loc,j)
                ! Only allow upward flux (downward diffusion handled by Thomas solver top BC)
                flux_up_q = max(flux_up_q, 0.0)
                flux_up_n = max(flux_up_n, 0.0)

                ! CFL stability limit
                max_flux_q = max(0.0, qs_fm(i,snc_N_loc,j)) * snc_dz(i,snc_N_loc,j) / dt
                max_flux_n = max(0.0, ns_fm(i,snc_N_loc,j)) * snc_dz(i,snc_N_loc,j) / dt
                flux_up_q = max(0.0, min(flux_up_q, max_flux_q))
                flux_up_n = max(0.0, min(flux_up_n, max_flux_n))

                ! ============================================================
                ! Downward flux: settling from atmospheric kts → fine mesh top
                ! ============================================================
                Vq_settle = 0.0
                Vn_settle = 0.0

                if (qs(i,kts,j) > BS_QS_MIN) then
                    ! Lambda from atmospheric mass and number (same as Operator B)
                    lambda = ((ns(i,kts,j) * RHO_ICE * PI_CONST * gamma(alpha + 3.0)) / &
                              (max(qs(i,kts,j), BS_QS_MIN) * 6.0 * gamma(alpha)))**(1.0/3.0)

                    ! Best number parameter (Mitchell 1996)
                    eta = rho_air * NU_AIR
                    phi_best = (4.0/3.0) * rho_air * (RHO_ICE - rho_air) * gravity / (eta * eta)

                    ! Mean Best number for regime selection
                    X_best = phi_best * (alpha/lambda)**3

                    ! Select Mitchell (1996) regime coefficients
                    if (X_best <= MITCHELL_X1) then
                        am_loc = MITCHELL_AM_S;  bm_loc = MITCHELL_BM_S
                    else if (X_best <= MITCHELL_X2) then
                        am_loc = MITCHELL_AM_I;  bm_loc = MITCHELL_BM_I
                    else
                        am_loc = MITCHELL_AM_T;  bm_loc = MITCHELL_BM_T
                    endif

                    V_coef = am_loc * eta * phi_best**(bm_loc) / rho_air

                    Vq_settle = V_coef * gamma(3*bm_loc + 2 + alpha) / (gamma(alpha) * gamma(alpha + 3)) * lambda**(1-3*bm_loc)
                    Vn_settle = V_coef * gamma(3*bm_loc - 1 + alpha) / (gamma(alpha)) * lambda**(1-3*bm_loc)
                    Vq_settle = max(SETTLE_MIN, min(Vq_settle, 2.0))
                    Vn_settle = max(SETTLE_MIN, min(Vn_settle, 2.0))
                endif

                ! Combined downward velocity: settling + downward advection
                Vq_down = Vq_settle + max(0.0, -w_fm(i,snc_N_loc,j))
                Vn_down = Vn_settle + max(0.0, -w_fm(i,snc_N_loc,j))

                ! Single-cell exact settling: fraction = min(V*dt/dz, 1.0)
                ! max(0,...) prevents negative qs from producing reversed fluxes
                flux_down_q = max(0.0, qs(i,kts,j)) * min(Vq_down * dt / dz_atm, 1.0) * dz_atm / dt
                flux_down_n = max(0.0, ns(i,kts,j)) * min(Vn_down * dt / dz_atm, 1.0) * dz_atm / dt

                ! CFL stability limit
                max_flux_q = max(0.0, qs(i,kts,j)) * dz_atm / dt
                max_flux_n = max(0.0, ns(i,kts,j)) * dz_atm / dt
                flux_down_q = max(0.0, min(flux_down_q, max_flux_q))
                flux_down_n = max(0.0, min(flux_down_n, max_flux_n))

                ! ============================================================
                ! Mass-conserving application
                ! ============================================================
                ! Upward flux: remove from fine mesh top, add to atmosphere
                ! Downward flux: remove from atmosphere, add to fine mesh top

                ! qs_flux_cache(i,j) = qs_fm(i,snc_N_loc,j) + max(0.0, (flux_down_q - flux_up_q) * dt)
                ! ns_flux_cache(i,j) = ns_fm(i,snc_N_loc,j) + max(0.0, (flux_down_n - flux_up_n) * dt)

                qs_flux_cache(i,j) = flux_up_q - flux_down_q
                ns_flux_cache(i,j) = flux_up_n - flux_down_n

            enddo
        enddo

        end associate

    end subroutine couple_to_3d_grid


    !>----------------------------------------------------------
    !! Accumulate diagnostic fields
    !!----------------------------------------------------------
    ! subroutine accumulate_diagnostics(domain)
    !     implicit none
    !     type(domain_t), intent(inout) :: domain
    !     integer :: i, j

    !     associate( &
    !         bs_swe_exch => domain%vars_2d(domain%var_indx(kVARS%bs_swe_exchange)%v)%data_2d, &
    !         bs_drift_swe_salt => domain%vars_2d(domain%var_indx(kVARS%bs_drift_swe_salt)%v)%data_2d, &
    !         bs_drift_swe_susp => domain%vars_2d(domain%var_indx(kVARS%bs_drift_swe_susp)%v)%data_2d, &
    !         bs_drift_swe_subl => domain%vars_2d(domain%var_indx(kVARS%bs_drift_swe_subl)%v)%data_2d &
    !     )

    !     ! Sublimation adds moisture and cools the air
    !     ! bs_subl is in W/m^2 (positive = sublimation occurring, energy consumed)
    !     !$acc parallel loop gang vector collapse(2) present(bs_swe_exch, bs_drift_swe_salt, bs_drift_swe_susp, bs_drift_swe_subl)
    !     do j = jts, jte
    !         do i = its, ite
    !             ! Mass exchange = net of saltation + suspension + sublimation
    !             bs_swe_exch(i,j) = &
    !                 bs_drift_swe_salt(i,j) + &
    !                 bs_drift_swe_susp(i,j) !+ &
    !                 ! bs_drift_swe_subl(i,j)
    !         enddo
    !     enddo
    !     end associate

    ! end subroutine accumulate_diagnostics


    !>----------------------------------------------------------
    !! Apply blowing snow sublimation feedback to atmosphere
    !! Runs on GPU with OpenACC, called from integrate_physics_tendencies
    !!----------------------------------------------------------
    subroutine snow_drift_apply_feedback(domain, options, dt)
        implicit none
        type(domain_t),  intent(inout) :: domain
        type(options_t), intent(in)    :: options
        real,            intent(in)    :: dt

        integer :: i, j, its_l, ite_l, jts_l, jte_l, kts_l
        real :: subl_flux

        if (.not.(options%sm%suspension_layer == 1 .and. options%sm%bs_atm_feedback .and. options%physics%snowmodel > 0)) return

        its_l = domain%its ; ite_l = domain%ite
        jts_l = domain%jts ; jte_l = domain%jte
        kts_l = domain%kts

        associate( &
            dz      => domain%vars_3d(domain%var_indx(kVARS%dz_interface)%v)%data_3d, &
            qs      => domain%vars_3d(domain%var_indx(kVARS%snow_mass)%v)%data_3d, &
            ns      => domain%vars_3d(domain%var_indx(kVARS%snow_number)%v)%data_3d, &
            qs_fm => domain%vars_3d(domain%var_indx(kVARS%qs_fm)%v)%data_3d, &
            ns_fm => domain%vars_3d(domain%var_indx(kVARS%ns_fm)%v)%data_3d &
        )

        !$acc parallel loop gang vector collapse(2) present(dz, qs, ns, qs_fm, ns_fm, qs_flux_cache, ns_flux_cache, snc_dz)
        do j = jts_l, jte_l
            do i = its_l, ite_l
                if (abs(qs_flux_cache(i,j)) > 0) then
                    ! Fluxes are positive upward
                    qs_fm(i,snc_N_loc,j) = max(qs_fm(i,snc_N_loc,j) - qs_flux_cache(i,j) * dt / snc_dz(i,snc_N_loc,j),0.0)
                    ns_fm(i,snc_N_loc,j) = max(ns_fm(i,snc_N_loc,j) - ns_flux_cache(i,j) * dt / snc_dz(i,snc_N_loc,j),0.0)

                    qs(i,kts,j) = max(qs(i,kts,j) + qs_flux_cache(i,j) * dt / dz(i,kts,j),0.0)
                    ns(i,kts,j) = max(ns(i,kts,j) + ns_flux_cache(i,j) * dt / dz(i,kts,j),0.0)
                endif
            enddo
        enddo

        end associate

    end subroutine snow_drift_apply_feedback



    !>----------------------------------------------------------
    !! Batch halo exchange for fine mesh variables (qs_fm, ns_fm)
    !! Uses the domain's existing batch exchange infrastructure with NCCL support
    !!----------------------------------------------------------
    subroutine exch_fine_mesh_3d(domain)
        implicit none
        type(domain_t), intent(inout) :: domain

        call domain%halo%halo_3d_send_batch(exch_vars_fm, domain%vars_3d)
        call domain%halo%halo_3d_retrieve_batch(exch_vars_fm, domain%vars_3d)

    end subroutine exch_fine_mesh_3d


    !>----------------------------------------------------------
    !! Thomas algorithm (tridiagonal solver)
    !! Solves: a(k)*x(k-1) + b(k)*x(k) + c(k)*x(k+1) = d(k)
    !!----------------------------------------------------------
    subroutine thomas_solve(n, a, b, c, d)
        !$acc routine seq
        implicit none
        integer, intent(in)    :: n
        real, intent(in)       :: a(n), c(n)
        real, intent(inout)    :: b(n), d(n)

        integer :: k
        real    :: w

        ! Forward sweep
        do k = 2, n
            w = a(k) / b(k-1)
            b(k) = b(k) - w * c(k-1)
            d(k) = d(k) - w * d(k-1)
        enddo

        ! Back substitution
        d(n) = d(n) / b(n)
        do k = n-1, 1, -1
            d(k) = (d(k) - c(k) * d(k+1)) / b(k)
        enddo

    end subroutine thomas_solve

end module snow_drift
