!!
!!----------------------------------------------------------
module module_sm_SNOWPACKdrv
    use icar_constants
    use options_interface,   only : options_t
    use meta_data_interface,  only : meta_data_t
    use domain_interface,    only : domain_t
    use io_routines, only : io_read, io_var_reversed
    use output_metadata,            only : get_varmeta
    use netcdf
    use iso_c_binding
    use meteoio_mod
    use snowpack_mod
    use snowpack_vector_helpers
    use namelist_utils, only : translate_numeric_mapping
    use mod_atm_utilities, only : relative_humidity
    use string, only : str, to_upper
    use omp_lib, only : omp_get_wtime
    use mod_wrf_constants, only: STBOLT

    implicit none

    type(snowpack_config) :: snp_cfg
    type(snow_station), allocatable :: stations(:,:)
    logical, allocatable :: run_snowpack_flag(:,:) 

    integer, parameter :: kMAX_SNOWPACK_ELEMENTS = 100
    logical :: truncation_warned = .false.

    integer :: ims, ime, jms, jme, its, ite, jts, jte
    integer :: i, j, k

    interface
        subroutine snowpack_SHROUD_memory_destructor(cap) &
            bind(C, name="snowpack_SHROUD_memory_destructor")
            import :: snowpack_SHROUD_capsule_data
            type(snowpack_SHROUD_capsule_data), intent(inout) :: cap
        end subroutine
    end interface

    private
    public :: sm_SNOWPACK_init,sm_SNOWPACK
   
    contains

    subroutine sm_SNOWPACK_init(domain, options, context_chng)
        ! Initialize SNOWPACK variables
        type(domain_t), intent(inout) :: domain
        type(options_t), intent(in) :: options
        logical, optional, intent(in) :: context_chng

        integer :: i,j,k, hj, hi, i_s, i_e, j_s, j_e
        logical :: context_change
        type(config) :: cfg
        type(element_data) :: elem
        integer(c_short) :: elem_id

        !!

        if (present(context_chng)) then
            context_change = context_chng
        else
            context_change = .False.
        endif

        ims = domain%grid%ims
        ime = domain%grid%ime
        jms = domain%grid%jms
        jme = domain%grid%jme
        its = domain%grid%its
        ite = domain%grid%ite
        jts = domain%grid%jts
        jte = domain%grid%jte

        ! Define SNOWPACK specific variables
        ! 1. Configuration
        cfg = config()

        ! Snowpack section (from CRYOWRF tests/io.ini)
        ! update_interval is in SECONDS (e.g. 300s = 5 min).
        ! CALCULATION_STEP_LENGTH is in MINUTES for SNOWPACK.
        ! Previously passed raw seconds (300) as minutes → 18000s timestep bug.
        call cfg%add_key("CALCULATION_STEP_LENGTH", "Snowpack", trim(str(int(options%lsm%update_interval / 60.0))))
        call cfg%add_key("METEO_STEP_LENGTH", "Snowpack", trim(str(int(options%lsm%update_interval / 60.0))))
        call cfg%add_key("HEIGHT_OF_WIND_VALUE", "Snowpack", "12.0")
        call cfg%add_key("HEIGHT_OF_METEO_VALUES", "Snowpack", "12.0")
        call cfg%add_key("ENFORCE_MEASURED_SNOW_HEIGHTS", "Snowpack", "false")
        call cfg%add_key("SW_MODE", "Snowpack", "INCOMING")
        call cfg%add_key("FORCING", "Snowpack", "ATMOS")
        call cfg%add_key("ATMOSPHERIC_STABILITY", "Snowpack", trim(to_upper(translate_numeric_mapping("snowpack_atmospheric_stability", options%sm%snowpack_atmospheric_stability))))
        call cfg%add_key("ROUGHNESS_LENGTH", "Snowpack", "0.002")
        call cfg%add_key("CHANGE_BC", "Snowpack", "false")
        call cfg%add_key("THRESH_CHANGE_BC", "Snowpack", "-1.0")
        call cfg%add_key("SNP_SOIL", "Snowpack", "false")
        call cfg%add_key("SOIL_FLUX", "Snowpack", "false")
        call cfg%add_key("GEO_HEAT", "Snowpack", "0.0")
        call cfg%add_key("CANOPY", "Snowpack", "false")
        call cfg%add_key("MEAS_TSS", "Snowpack", "false")

        ! SnowpackAdvanced section (from CRYOWRF tests/io.ini)
        call cfg%add_key("VARIANT", "SnowpackAdvanced", trim(to_upper(translate_numeric_mapping("snowpack_variant", options%sm%snowpack_variant))))
        call cfg%add_key("ALPINE3D", "SnowpackAdvanced", "false")
        call cfg%add_key("SNOW_EROSION", "SnowpackAdvanced", "false")
        call cfg%add_key("SNOW_REDISTRIBUTION", "SnowpackAdvanced", "false")
        call cfg%add_key("COUPLEDPHASECHANGES", "SnowpackAdvanced", "true")
        call cfg%add_key("HN_DENSITY", "SnowpackAdvanced", "PARAMETERIZED")
        call cfg%add_key("HN_DENSITY_PARAMETERIZATION", "SnowpackAdvanced", "LEHNING_NEW")
        call cfg%add_key("HN_DENSITY_FIXEDVALUE", "SnowpackAdvanced", "100.0")
        call cfg%add_key("RIME_INDEX", "SnowpackAdvanced", "false")
        call cfg%add_key("NEWSNOW_LWC", "SnowpackAdvanced", "false")
        call cfg%add_key("READ_DSM", "SnowpackAdvanced", "false")
        call cfg%add_key("SNOW_ALBEDO", "SnowpackAdvanced", "PARAMETERIZED")
        call cfg%add_key("ALBEDO_PARAMETERIZATION", "SnowpackAdvanced", trim(to_upper(translate_numeric_mapping("snowpack_albedo_parameterization", options%sm%snowpack_albedo_parameterization))))
        call cfg%add_key("ALBEDO_AVERAGE_SCHMUCKI", "SnowpackAdvanced", "ALL_DATA")
        call cfg%add_key("ALBEDO_FIXEDVALUE", "SnowpackAdvanced", "0.65")
        call cfg%add_key("ALBEDO_AGING", "SnowpackAdvanced", "false")
        call cfg%add_key("SW_ABSORPTION_SCHEME", "SnowpackAdvanced", "MULTI_BAND")
        call cfg%add_key("WATERTRANSPORTMODEL_SNOW", "SnowpackAdvanced", "BUCKET")
        call cfg%add_key("WATERTRANSPORTMODEL_SOIL", "SnowpackAdvanced", "BUCKET")
        call cfg%add_key("LB_COND_WATERFLUX", "SnowpackAdvanced", "FREEDRAINAGE")
        call cfg%add_key("METAMORPHISM_MODEL", "SnowpackAdvanced", "DEFAULT")
        call cfg%add_key("STRENGTH_MODEL", "SnowpackAdvanced", "DEFAULT")
        call cfg%add_key("VISCOSITY_MODEL", "SnowpackAdvanced", "DEFAULT")
        call cfg%add_key("MINIMUM_L_ELEMENT", "SnowpackAdvanced", "0.001")
        call cfg%add_key("HEIGHT_NEW_ELEM", "SnowpackAdvanced", "0.001")
        call cfg%add_key("COMBINE_ELEMENTS", "SnowpackAdvanced", "true")
        call cfg%add_key("COMB_THRESH_L", "SnowpackAdvanced", "0.01")
        call cfg%add_key("REDUCE_N_ELEMENTS", "SnowpackAdvanced", str(options%sm%snowpack_reduce_n_elements))
        call cfg%add_key("FIXED_POSITIONS", "SnowpackAdvanced", "0.25 0.5 1.0")
        call cfg%add_key("NUMBER_FIXED_RATES", "SnowpackAdvanced", "0")
        call cfg%add_key("MAX_NUMBER_MEAS_TEMPERATURES", "SnowpackAdvanced", "7")
        call cfg%add_key("MIN_DEPTH_SUBSURF", "SnowpackAdvanced", "0.0")
        call cfg%add_key("MEAS_INCOMING_LONGWAVE", "SnowpackAdvanced", "true")
        call cfg%add_key("ALLOW_ADAPTIVE_TIMESTEPPING", "SnowpackAdvanced", "true")
        call cfg%add_key("NUMBER_SLOPES", "SnowpackAdvanced", "1")
        call cfg%add_key("MASS_BALANCE", "SnowpackAdvanced", "false")
        call cfg%add_key("FORCE_ADD_SNOWFALL", "SnowpackAdvanced", "true")
        call cfg%add_key("ADJUST_HEIGHT_OF_WIND_VALUE", "SnowpackAdvanced", "false")
        call cfg%add_key("SOOT_PPMV", "SnowpackAdvanced", "0.0")
        call cfg%add_key("ENABLE_VAPOUR_TRANSPORT", "SnowpackAdvanced", str(options%sm%snowpack_enable_vapour_transport,c_bool=.True.))
        call cfg%add_key("JAM", "SnowpackAdvanced", "false")
        call cfg%add_key("THRESH_DTEMP_AIR_SNOW", "SnowpackAdvanced", "60.0")
        call cfg%add_key("THRESH_RAIN", "SnowpackAdvanced", "1.2")
        call cfg%add_key("THRESH_RH", "SnowpackAdvanced", "0.0")
        call cfg%add_key("WATER_LAYER", "SnowpackAdvanced", "false")
        call cfg%add_key("T_CRAZY_MAX", "SnowpackAdvanced", "340.0")
        call cfg%add_key("T_CRAZY_MIN", "SnowpackAdvanced", "120.0")
        call cfg%add_key("NEW_SNOW_GRAIN_SIZE", "SnowpackAdvanced", "0.2")
        call cfg%add_key("SNOW_PREPARATION", "SnowpackAdvanced", "false")
        call cfg%add_key("ADVECTIVE_HEAT", "SnowpackAdvanced", "false")
        call cfg%add_key("WIND_SCALING_FACTOR", "SnowpackAdvanced", "1.0")
        call cfg%add_key("PERP_TO_SLOPE", "SnowpackAdvanced", "false")
        call cfg%add_key("DETECT_GRASS", "SnowpackAdvanced", "false")
        call cfg%add_key("PREF_FLOW", "SnowpackAdvanced", "false")
        call cfg%add_key("ALLOW_FREEZING_RAIN", "SnowpackAdvanced", "true")

        snp_cfg = snowpack_config(cfg)

        if (allocated(stations)) deallocate(stations)
        if (allocated(run_snowpack_flag)) deallocate(run_snowpack_flag)
        allocate(stations(its:ite,jts:jte))
        allocate(run_snowpack_flag(its:ite,jts:jte))

        run_snowpack_flag = domain%vars_2d(domain%var_indx(kVARS%land_mask)%v)%data_2di(its:ite,jts:jte)/=kLC_WATER

        ! ---------------------------- Setup 2D grid of snow stations ----------------------------
        do j = jts, jte
            do i = its, ite
                
                ! cycle if this is a water pixel
                if (.not.run_snowpack_flag(i,j)) cycle

                stations(i,j) = snow_station(.false., .false.)

                call snowpack_station_resize(stations(i,j)%cxxmem%addr, 0_c_size_t)

                ! Station-level properties (must be set after init)
                call stations(i,j)%set_cos_sl(real(cos(domain%vars_2d(domain%var_indx(kVARS%slope_angle)%v)%data_2d(i,j)),kind=8))       ! cos(slope angle) - flat ground
                call stations(i,j)%set_sector(0_c_size_t)  ! Slope sector
                call stations(i,j)%set_altitude(real(domain%vars_2d(domain%var_indx(kVARS%terrain)%v)%data_2d(i,j), kind=8))
                call stations(i,j)%set_p_albedo(0.2d0)     ! Parameterized albedo (bare soil)
                call stations(i,j)%set_albedo(real(domain%vars_2d(domain%var_indx(kVARS%albedo)%v)%data_2d(i,j), kind=8))       ! Surface albedo
                call stations(i,j)%set_soil_alb(0.2d0)     ! Soil albedo - typical bare soil
                call stations(i,j)%set_soil_emissivity(0.98d0) ! Soil emissivity
                call stations(i,j)%set_bare_soil_z0(0.002d0)   ! Bare soil roughness [m]
                call stations(i,j)%set_soil_node(0_c_size_t)   ! No soil layers

            enddo
        enddo

        ! Read initial SNOWPACK state from domain file if variable names provided
        if (options%domain%snowpack_nlayers_var /= "") then
            call read_snowpack_state(domain, options)
        endif
        
        call update_run_snowpack_flag(domain)

        ! copy in any initial snowpack state if we are supposed to run snowpack already at a given pixel
        call SnowStationsIn(domain, stations)


    end subroutine sm_SNOWPACK_init

    subroutine sm_SNOWPACK(domain,options,lsm_dt,current_rain,current_snow,windspd)
        ! SNOWPACK snow model step
        type(domain_t), intent(inout) :: domain
        type(options_t), intent(in) :: options
        real, intent(in) :: lsm_dt
        real, allocatable, dimension(:,:), intent(in) :: current_rain, current_snow, windspd

        integer :: i, j
        real(c_double) :: cumu_precip
        type(snowpack) :: snp
        type(bound_cond), allocatable :: bdata(:,:)
        type(surface_fluxes), allocatable :: sdata(:,:)
        type(current_meteo), allocatable :: meteo(:,:)

        real(kind=8) :: start_time, end_time, t_setup, t_loop, t_copyout, t_cleanup

        ! Initialize cumulative precipitation
        cumu_precip = 0.0d0 ! dummy input var if forceaddsnowfall is set to true in the config options

        call update_run_snowpack_flag(domain,current_snow)

        ! start timing loop, get wall-clock time
        start_time = omp_get_wtime()

        !! ------------------------------- copy existing snowpack state (TODO: future work, for now just use persistent station array with single nest) -------------------------------
        call MeteoIn(domain, current_rain, current_snow, windspd, meteo)
        call SurfaceFluxesIn(domain, sdata)
        call BoundaryConditionsIn(domain, bdata)
        t_setup = omp_get_wtime()

        ! Loop over horizontal grid points
        ! Create one Snowpack solver per thread (not per cell) to avoid
        ! massive concurrent malloc/free pressure from C++ heap contention
        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, snp) FIRSTPRIVATE(cumu_precip)
        snp = snowpack(snp_cfg)
        !$OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC,1)
        do j = jts, jte
            do i = its, ite
            !! ------------------------------- Run snowpack step -------------------------------
            if (run_snowpack_flag(i,j)) then
                call snp%run_snowpack_model(meteo(i,j), stations(i,j), cumu_precip, bdata(i,j), sdata(i,j))
                call sdata(i,j)%collect_surface_fluxes(bdata(i,j), stations(i,j), meteo(i,j))
            endif
            end do
        end do
        !$OMP END DO
        call snowpack_SHROUD_memory_destructor(snp%cxxmem)
        !$OMP END PARALLEL
        t_loop = omp_get_wtime()

        !! ------------------------------- copy out new snowpack state -------------------------------
        call SurfaceFluxesOut(domain, sdata)
        call BoundaryConditionsOut(domain, bdata)
        call SnowStationsOut(domain, stations)
        t_copyout = omp_get_wtime()

        ! Clean up C++ objects before local arrays go out of scope
        do j = jts, jte
            do i = its, ite
                if (run_snowpack_flag(i,j)) then
                    call snowpack_SHROUD_memory_destructor(meteo(i,j)%cxxmem)
                    call snowpack_SHROUD_memory_destructor(sdata(i,j)%cxxmem)
                    call snowpack_SHROUD_memory_destructor(bdata(i,j)%cxxmem)
                endif
            end do
        end do

        ! end timing
        end_time = omp_get_wtime()
        if (STD_OUT_PE) write(*,'(A,F8.3,A,4(A,F8.3))') &
            "SNOWPACK total: ", end_time - start_time, "s", &
            "  setup=", t_setup - start_time, &
            "  loop=", t_loop - t_setup, &
            "  copyout=", t_copyout - t_loop, &
            "  cleanup=", end_time - t_copyout

    end subroutine sm_SNOWPACK


    !>----------------------------------- MeteoIn ----------------------------------
    ! Initializes and fills meteo data structure for points that will run SNOWPACK
    !>------------------------------------------------------------------------------
    subroutine MeteoIn(domain,current_rain,current_snow,windspd,meteo)
        ! Copy meteorological data from SNOWPACK into HICAR variable structures
        type(domain_t), intent(in) :: domain
        real, allocatable, dimension(:,:), intent(in) :: current_rain, current_snow, windspd
        type(current_meteo), allocatable, intent(out) :: meteo(:,:)

        real(kind=8), dimension(its:ite,jts:jte) :: total_precip, precip_phase, rh
        real(kind=8) :: lw_out

        allocate(meteo(its:ite,jts:jte))

        associate(                                                                                           &
                    temperature => domain%vars_3d(domain%var_indx(kVARS%temperature)%v)%data_3d,             &
                    qv => domain%vars_3d(domain%var_indx(kVARS%water_vapor)%v)%data_3d, &
                    p => domain%vars_3d(domain%var_indx(kVARS%pressure)%v)%data_3d,                       &
                    roughness_z0 => domain%vars_2d(domain%var_indx(kVARS%roughness_z0)%v)%data_2d,       &
                    incoming_shortwave => domain%vars_2d(domain%var_indx(kVARS%shortwave)%v)%data_2d, &
                    ! shortwave_terrain => domain%vars_2d(domain%var_indx(kVARS%shortwave_terrain)%v)%data_2d, &
                    longwave_in => domain%vars_2d(domain%var_indx(kVARS%longwave)%v)%data_2d,   &
                    ground_temp => domain%vars_3d(domain%var_indx(kVARS%soil_temperature)%v)%data_3d,               &
                    skin_temperature => domain%vars_2d(domain%var_indx(kVARS%skin_temperature)%v)%data_2d              &
                 )

        total_precip = current_rain(its:ite,jts:jte) + current_snow(its:ite,jts:jte)
        precip_phase = current_rain(its:ite,jts:jte) / (total_precip+1.0d-6)  ! avoid div by zero
        rh = relative_humidity(temperature(its:ite,1,jts:jte),qv(its:ite,1,jts:jte),p(its:ite,1,jts:jte))

        ! Loop over horizontal grid points
        do j = jts, jte
            do i = its, ite
                if (.not.(run_snowpack_flag(i,j))) cycle
                meteo(i,j) = current_meteo()

                !! ------------------------------- Set meteo variables -------------------------------
                call meteo(i,j)%set_date_julian(real(domain%sim_time%mjd(),kind=8))  ! Date (Julian day 2460310 ≈ Jan 1, 2024)
                call meteo(i,j)%set_ta(real(temperature(i,1,j), kind=8))          ! Air temperature [K] - CRYOWRF Coupler: l_TA
                call meteo(i,j)%set_rh(rh(i,j))            ! Relative humidity [0-1] - CRYOWRF Coupler: l_RH
                call meteo(i,j)%set_rh_avg(rh(i,j))        ! Running mean RH, set same as RH
                call meteo(i,j)%set_vw(real(windspd(i,j), kind=8))             ! Wind velocity [m/s] - CRYOWRF Coupler: l_VW
                ! call meteo(i,j)%set_vw_avg(5.0d0)         ! Running mean wind velocity
                ! call meteo(i,j)%set_vw_max(8.0d0)         ! Max wind velocity - CRYOWRF Coupler: l_VW_MAX
                call meteo(i,j)%set_vw_drift(real(windspd(i,j), kind=8))       ! Drift wind velocity
                ! call meteo(i,j)%set_ustar(0.3d0)          ! Friction velocity [m/s]
                call meteo(i,j)%set_z0(real(roughness_z0(i,j), kind=8))           ! Roughness length [m] - matches ROUGHNESS_LENGTH config -- CURRENTLY OVERWRITTEN IN SNOWPACK
                call meteo(i,j)%set_iswr(real(incoming_shortwave(i,j), kind=8))         ! Incoming shortwave [W/m2] - CRYOWRF Coupler: l_iswr
                ! call meteo(i,j)%set_rswr(real(shortwave_terrain(i,j), kind=8))    ! Reflected shortwave [W/m2] from HICAR terrain radiation

                call meteo(i,j)%set_ea(-999.0d0)            ! Atmospheric emissivity - CRYOWRF Coupler: Mdata.ea, setting to -999 to force recalculation in SNOWPACK

                ! Compute net longwave from incoming LW (HICAR) minus outgoing LW (surface emission)
                ! Positive = net energy into snowpack
                lw_out = 0.98d0 * STBOLT * skin_temperature(i,j)**4
                call meteo(i,j)%set_lw_net(real(longwave_in(i,j) - lw_out, kind=8))
                call meteo(i,j)%set_tss(real(skin_temperature(i,j), kind=8))          ! Snow surface temperature [K]

                ! NOTE: because SOIL_FLUX is set to false above, ts0 has to be the ground surface temperature from NoahMP
                call meteo(i,j)%set_ts0(real(ground_temp(i,1,j), kind=8))          ! Bottom temperature [K] - CRYOWRF Coupler: TSG
                call meteo(i,j)%set_psum(total_precip(i,j))           ! Precipitation sum since previous SNOWPACK call [mm] - CRYOWRF Coupler: l_psum
                call meteo(i,j)%set_psum_ph(precip_phase(i,j))        ! Precip phase (0=solid, 1=liquid)

                
                call meteo(i,j)%set_dw(0.0d0)           ! Wind direction [deg] - CRYOWRF Coupler: l_DW
                call meteo(i,j)%set_dw_drift(0.0d0)     ! Drift wind direction
                call meteo(i,j)%set_psi_s(0.0d0)          ! Stability correction
                call meteo(i,j)%set_hs(0.0d0)             ! Measured snow height [m] -- needed to run, but only used if running on IMIS data, which doesn't happen when coupled to HICAR
                call meteo(i,j)%set_geo_heat(0.0d0)       ! Geothermal heat flux [W/m2] - CRYOWRF io.ini: GEO_HEAT
            end do
        end do

        end associate

    end subroutine MeteoIn

    !>--------------------------------- SurfaceFluxesIn --------------------------------
    ! Initializes surface fluxes structure for points that will run SNOWPACK
    !>---------------------------------------------------------------------------------
    subroutine SurfaceFluxesIn(domain,sfluxes)
        type(domain_t), intent(in) :: domain
        type(surface_fluxes), allocatable, intent(out) :: sfluxes(:,:)

        allocate(sfluxes(its:ite,jts:jte))
        ! Loop over horizontal grid points
        do j = jts, jte
            do i = its, ite
                if (.not.(run_snowpack_flag(i,j))) cycle
                sfluxes(i,j) = surface_fluxes()
            enddo
        enddo

    end subroutine SurfaceFluxesIn

    !>-------------------------------------------------- SurfaceFluxesOut -------------------------------------------------
    ! Copies out data from surface fluxes data structure (summed fluxes since last SNOWPACK call, i.e. mainly mass fluxes)
    !>---------------------------------------------------------------------------------------------------------------------
    subroutine SurfaceFluxesOut(domain,sfluxes)
        type(domain_t), intent(inout) :: domain
        type(surface_fluxes), intent(in) :: sfluxes(its:ite,jts:jte)

        associate(                                                                                           &
                    sublimation => domain%vars_2d(domain%var_indx(kVARS%dSWE_subl)%v)%data_2d,       &
                    meltflux_out => domain%vars_2d(domain%var_indx(kVARS%meltflux_out_tstep)%v)%data_2d &
                 )
        ! Loop over horizontal grid points
        do j = jts, jte
            do i = its, ite
                if (.not.(run_snowpack_flag(i,j))) cycle
                
                ! Negate: SNOWPACK uses positive = into snowpack (deposition),
                ! HICAR uses positive = mass lost from snow (sublimation)
                sublimation(i,j) = -snowpack_get_mass_sublimation(sfluxes(i,j)%cxxmem%addr)

                ! Snowmelt runoff (kg/m²), positive = water leaving snowpack base
                meltflux_out(i,j) = snowpack_get_mass_runoff(sfluxes(i,j)%cxxmem%addr)

            enddo
        enddo

        end associate

    end subroutine SurfaceFluxesOut

    !>--------------------------------- BoundaryConditionsIn --------------------------------
    ! Initializes boundary conditions structure for points that will run SNOWPACK
    !>---------------------------------------------------------------------------------------
    subroutine BoundaryConditionsIn(domain,bconds)
        type(domain_t), intent(in) :: domain
        type(bound_cond), allocatable, intent(out) :: bconds(:,:)

        allocate(bconds(its:ite,jts:jte))

        ! Loop over horizontal grid points
        do j = jts, jte
            do i = its, ite
                if (.not.(run_snowpack_flag(i,j))) cycle
                bconds(i,j) = bound_cond()
            enddo
        enddo

    end subroutine BoundaryConditionsIn

    !>---------------------------------------- BoundaryConditionsOut ---------------------------------------
    ! Copies out state of SNOWPACK boundary conditions (instantaneous fluxes) to the HICAR domain structure
    !>------------------------------------------------------------------------------------------------------
    subroutine BoundaryConditionsOut(domain,bconds)
        type(domain_t), intent(inout) :: domain
        type(bound_cond), intent(inout) :: bconds(its:ite,jts:jte)

        associate(                                                                                           &
                    sensible_heat => domain%vars_2d(domain%var_indx(kVARS%sensible_heat)%v)%data_2d,               &
                    latent_heat => domain%vars_2d(domain%var_indx(kVARS%latent_heat)%v)%data_2d                 &
                 )

        ! Loop over horizontal grid points
        do j = jts, jte
            do i = its, ite
                if (run_snowpack_flag(i,j)) then
                    ! Negate: SNOWPACK uses positive-downward (into snow),
                    ! HICAR uses positive-upward (into atmosphere)
                    sensible_heat(i,j) = -bconds(i,j)%get_qs()
                    latent_heat(i,j) = -bconds(i,j)%get_ql()
                endif
            enddo
        enddo

        end associate

    end subroutine BoundaryConditionsOut


    !>--------------------------------- SnowStationsIn --------------------------------
    ! Copies in snow stations state from HICAR domain for points that will run SNOWPACK
    !>---------------------------------------------------------------------------------
    subroutine SnowStationsIn(domain,stations_in)
        type(domain_t), intent(in) :: domain
        type(snow_station), intent(inout) :: stations_in(its:ite,jts:jte)

        integer(c_size_t) :: n_elem, n_node
        integer :: ne_arr(1:kMAX_SNOWPACK_ELEMENTS)
        integer(c_short) :: sn_k, hicar_k
        type(element_data) :: elem

        ne_arr = 1

        associate(                                                                                           &
                    snow_height => domain%vars_2d(domain%var_indx(kVARS%snow_height)%v)%data_2d,               &
                    snow_water_equivalent => domain%vars_2d(domain%var_indx(kVARS%snow_water_equivalent)%v)%data_2d,                 &
                    albedo => domain%vars_2d(domain%var_indx(kVARS%albedo)%v)%data_2d,                 &
                    n_snow_layers => domain%vars_2d(domain%var_indx(kVARS%snow_nlayers)%v)%data_2di,                 &
                    depositionDate => domain%vars_3d(domain%var_indx(kVARS%depositionDate)%v)%data_3d,       &
                    Layer_Thick => domain%vars_3d(domain%var_indx(kVARS%Ds)%v)%data_3d,             &
                    T_elem => domain%vars_3d(domain%var_indx(kVARS%snow_temperature)%v)%data_3d,                                 &
                    T_node => domain%vars_3d(domain%var_indx(kVARS%snow_temperature_i)%v)%data_3d,                                 &
                    Vol_Frac_I => domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_I)%v)%data_3d,               &
                    Vol_Frac_W => domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_W)%v)%data_3d,               &
                    Vol_Frac_A => domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_A)%v)%data_3d,               &
                    Vol_Frac_S => domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_S)%v)%data_3d,               &
                    Vol_Frac_WP => domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_WP)%v)%data_3d,             &
                    Rg => domain%vars_3d(domain%var_indx(kVARS%Rg)%v)%data_3d,                               &
                    Rb => domain%vars_3d(domain%var_indx(kVARS%Rb)%v)%data_3d,                               &
                    Dd => domain%vars_3d(domain%var_indx(kVARS%Dd)%v)%data_3d,                               &
                    Sp => domain%vars_3d(domain%var_indx(kVARS%Sp)%v)%data_3d,                               &
                    mk => domain%vars_3d(domain%var_indx(kVARS%mk)%v)%data_3d,                               &
                    mass_hoar => domain%vars_3d(domain%var_indx(kVARS%mass_hoar)%v)%data_3d,                 &
                    CDot => domain%vars_3d(domain%var_indx(kVARS%CDot)%v)%data_3d,                           &
                    snow_stress => domain%vars_3d(domain%var_indx(kVARS%snow_stress)%v)%data_3d,                      &
                    N3 => domain%vars_3d(domain%var_indx(kVARS%N3)%v)%data_3d                        &
                 )

        do j = jts, jte
            do i = its, ite
                if (.not.(run_snowpack_flag(i,j))) cycle

                do k = 1, n_snow_layers(i,j)
                    sn_k = k
                    hicar_k = n_snow_layers(i,j) - k + 1
                    elem = element_data(sn_k)               
                    call elem%set_theta([real(Vol_Frac_S(i,hicar_k,j),kind=8), &
                                            real(Vol_Frac_I(i,hicar_k,j),kind=8), &
                                            real(Vol_Frac_W(i,hicar_k,j),kind=8), &
                                            real(Vol_Frac_WP(i,hicar_k,j),kind=8), &
                                            real(Vol_Frac_A(i,hicar_k,j),kind=8)])

                    call elem%set_deposition_date_julian(real(depositionDate(i,hicar_k,j),kind=8))
                    ! Add element to station using the new addElement method
                    call stations_in(i,j)%add_element(elem)
                enddo

                n_elem = snowpack_get_num_elements(stations_in(i,j)%cxxmem%addr)
                n_node = snowpack_get_num_nodes(stations_in(i,j)%cxxmem%addr)

                call snowpack_set_all_node_T(stations_in(i,j)%cxxmem%addr, real(T_node(i,n_node:1:-1,j), kind=8), n_node)
                call snowpack_set_all_element_L(stations_in(i,j)%cxxmem%addr, real(Layer_Thick(i,n_elem:1:-1,j), kind=8), n_elem)
                call snowpack_set_all_element_Te(stations_in(i,j)%cxxmem%addr, real(T_elem(i,n_elem:1:-1,j), kind=8), n_elem)
                call snowpack_set_all_element_rg(stations_in(i,j)%cxxmem%addr, real(Rg(i,n_elem:1:-1,j), kind=8), n_elem)
                call snowpack_set_all_element_rb(stations_in(i,j)%cxxmem%addr, real(Rb(i,n_elem:1:-1,j), kind=8), n_elem)
                call snowpack_set_all_element_dd(stations_in(i,j)%cxxmem%addr, real(Dd(i,n_elem:1:-1,j), kind=8), n_elem)
                call snowpack_set_all_element_sp(stations_in(i,j)%cxxmem%addr, real(Sp(i,n_elem:1:-1,j), kind=8), n_elem)
                call snowpack_set_all_element_mk(stations_in(i,j)%cxxmem%addr, int(mk(i,n_elem:1:-1,j), kind=c_short), n_elem)
                call snowpack_set_all_element_CDot(stations_in(i,j)%cxxmem%addr, real(CDot(i,n_elem:1:-1,j), kind=8), n_elem)
                ! call snowpack_set_all_element_metamo(stations_in(i,j)%cxxmem%addr, real(metamo(i,n_elem:1:-1,j), kind=8), n_elem)
                ! call snowpack_set_all_element_N3(stations_in(i,j)%cxxmem%addr, real(N3(i,n_elem:1:-1,j), kind=8), n_elem)

                call stations_in(i,j)%set_c_h(real(snow_height(i,j), kind=8))
                call stations_in(i,j)%set_albedo(real(albedo(i,j), kind=8))

            enddo
        enddo

        end associate

    end subroutine SnowStationsIn


    !>--------------------------------- SnowStationsOut ------------------------
    ! Copies out state of SNOWPACK snow stations to the HICAR domain structure
    !>--------------------------------------------------------------------------
    subroutine SnowStationsOut(domain,stations_in)
        type(domain_t), intent(inout) :: domain
        type(snow_station), intent(inout) :: stations_in(its:ite,jts:jte)

        integer(c_size_t) :: n_elem, n_node
        integer :: n_save, n_save_nodes
        real(c_double) :: tmp_arr(kMAX_SNOWPACK_ELEMENTS), tmp_node_arr(kMAX_SNOWPACK_ELEMENTS+1)
        integer(c_short) :: tmp_mk(kMAX_SNOWPACK_ELEMENTS)

        associate(                      &
                    snow_height => domain%vars_2d(domain%var_indx(kVARS%snow_height)%v)%data_2d,               &
                    snow_water_equivalent => domain%vars_2d(domain%var_indx(kVARS%snow_water_equivalent)%v)%data_2d,                 &
                    land_emissivity => domain%vars_2d(domain%var_indx(kVARS%land_emissivity)%v)%data_2d,                 &
                    roughness_z0 => domain%vars_2d(domain%var_indx(kVARS%roughness_z0)%v)%data_2d,                 &
                    skin_temperature => domain%vars_2d(domain%var_indx(kVARS%skin_temperature)%v)%data_2d,                 &
                    ground_temp => domain%vars_3d(domain%var_indx(kVARS%soil_temperature)%v)%data_3d,                 &
                    albedo => domain%vars_2d(domain%var_indx(kVARS%albedo)%v)%data_2d,                 &
                    n_snow_layers => domain%vars_2d(domain%var_indx(kVARS%snow_nlayers)%v)%data_2di,                 &
                    depositionDate => domain%vars_3d(domain%var_indx(kVARS%depositionDate)%v)%data_3d,       &
                    Layer_Thick => domain%vars_3d(domain%var_indx(kVARS%Ds)%v)%data_3d,             &
                    T_elem => domain%vars_3d(domain%var_indx(kVARS%snow_temperature)%v)%data_3d,                                 &
                    T_node => domain%vars_3d(domain%var_indx(kVARS%snow_temperature_i)%v)%data_3d,                                 &
                    Vol_Frac_I => domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_I)%v)%data_3d,               &
                    Vol_Frac_W => domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_W)%v)%data_3d,               &
                    Vol_Frac_A => domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_A)%v)%data_3d,               &
                    Vol_Frac_S => domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_S)%v)%data_3d,               &
                    Vol_Frac_WP => domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_WP)%v)%data_3d,             &
                    Rg => domain%vars_3d(domain%var_indx(kVARS%Rg)%v)%data_3d,                               &
                    Rb => domain%vars_3d(domain%var_indx(kVARS%Rb)%v)%data_3d,                               &
                    Dd => domain%vars_3d(domain%var_indx(kVARS%Dd)%v)%data_3d,                               &
                    Sp => domain%vars_3d(domain%var_indx(kVARS%Sp)%v)%data_3d,                               &
                    mk => domain%vars_3d(domain%var_indx(kVARS%mk)%v)%data_3d,                               &
                    mass_hoar => domain%vars_3d(domain%var_indx(kVARS%mass_hoar)%v)%data_3d,                 &
                    CDot => domain%vars_3d(domain%var_indx(kVARS%CDot)%v)%data_3d,                           &
                    snow_stress => domain%vars_3d(domain%var_indx(kVARS%snow_stress)%v)%data_3d,                      &
                    N3 => domain%vars_3d(domain%var_indx(kVARS%N3)%v)%data_3d                               &
            )
        do j = jts, jte
            do i = its, ite
                if (run_snowpack_flag(i,j)) then
                
                    ! First set bulk station properties
                    snow_height(i,j) = stations_in(i,j)%get_c_h()
                    snow_water_equivalent(i,j) = stations_in(i,j)%get_swe()

                    if (snow_height(i,j) > 0.0) then
                        albedo(i,j) = stations_in(i,j)%get_albedo()
                        n_elem = snowpack_get_num_elements(stations_in(i,j)%cxxmem%addr)
                        n_node = snowpack_get_num_nodes(stations_in(i,j)%cxxmem%addr)
                        land_emissivity(i,j) = 0.98 
                        roughness_z0(i,j) = 0.002 !m
                    else
                        albedo(i,j) = stations_in(i,j)%get_soil_alb()
                        n_elem = 0
                        n_node = 0
                        skin_temperature(i,j) = 273.15
                        ground_temp(i,1,j) = 273.15
                    endif

                    ! Safety check: SNOWPACK element count must fit in tmp buffer
                    if (n_elem > kMAX_SNOWPACK_ELEMENTS) then
                        write(*,'(A,I0,A,I0,A,I0,A,I0,A)') "ERROR: SNOWPACK element count (", n_elem, &
                            ") exceeds buffer size (", kMAX_SNOWPACK_ELEMENTS, ") at (", i, ",", j, ")"
                        stop "SNOWPACK element count exceeds kMAX_SNOWPACK_ELEMENTS"
                    endif

                    ! Clamp to domain array size — keep top (surface) layers, discard deeper layers
                    ! (following CRYOWRF approach: SNOWPACK C++ objects retain full resolution,
                    !  domain arrays store only the top sm_nsnow_max layers for output/restart)
                    n_save = min(int(n_elem), kSNOW_GRID_Z)
                    n_save_nodes = n_save + 1

                    if (n_elem > kSNOW_GRID_Z .and. .not. truncation_warned) then
                        write(*,'(A,I0,A,I0,A)') "WARNING: SNOWPACK has ", n_elem, &
                            " layers but sm_nsnow_max = ", kSNOW_GRID_Z, &
                            ". Truncating deeper layers from domain arrays."
                        truncation_warned = .true.
                    endif

                    ! Store clamped layer count in domain (consistent with domain array sizes)
                    n_snow_layers(i,j) = n_save

                    ! Copy layer data from SNOWPACK to domain arrays.
                    ! SNOWPACK order: tmp_arr(1)=ground, tmp_arr(n_elem)=surface
                    ! HICAR order: layer 1=surface, layer n=ground (reversed)
                    ! With truncation: keep top n_save from SNOWPACK = indices n_elem down to n_elem-n_save+1

                    call snowpack_get_all_node_T(stations_in(i,j)%cxxmem%addr, tmp_node_arr, n_node)
                    if (n_save_nodes > 0) then
                        T_node(i,1:n_save_nodes,j) = real(tmp_node_arr(n_node:n_node-n_save_nodes+1:-1))
                        skin_temperature(i,j) = T_node(i,1,j)
                        ground_temp(i,1,j) = T_node(i,n_save_nodes,j)
                    endif

                    call snowpack_get_all_element_deposition_julian(stations_in(i,j)%cxxmem%addr, tmp_arr, n_elem)
                    if (n_save > 0) depositionDate(i,1:n_save,j) = real(tmp_arr(n_elem:n_elem-n_save+1:-1))

                    call snowpack_get_all_element_L(stations_in(i,j)%cxxmem%addr, tmp_arr, n_elem)
                    if (n_save > 0) Layer_Thick(i,1:n_save,j) = real(tmp_arr(n_elem:n_elem-n_save+1:-1))

                    call snowpack_get_all_element_Te(stations_in(i,j)%cxxmem%addr, tmp_arr, n_elem)
                    if (n_save > 0) T_elem(i,1:n_save,j) = real(tmp_arr(n_elem:n_elem-n_save+1:-1))

                    call snowpack_get_all_element_theta_ice(stations_in(i,j)%cxxmem%addr, tmp_arr, n_elem)
                    if (n_save > 0) Vol_Frac_I(i,1:n_save,j) = real(tmp_arr(n_elem:n_elem-n_save+1:-1))

                    call snowpack_get_all_element_theta_water(stations_in(i,j)%cxxmem%addr, tmp_arr, n_elem)
                    if (n_save > 0) Vol_Frac_W(i,1:n_save,j) = real(tmp_arr(n_elem:n_elem-n_save+1:-1))

                    call snowpack_get_all_element_theta_air(stations_in(i,j)%cxxmem%addr, tmp_arr, n_elem)
                    if (n_save > 0) Vol_Frac_A(i,1:n_save,j) = real(tmp_arr(n_elem:n_elem-n_save+1:-1))

                    call snowpack_get_all_element_theta_soil(stations_in(i,j)%cxxmem%addr, tmp_arr, n_elem)
                    if (n_save > 0) Vol_Frac_S(i,1:n_save,j) = real(tmp_arr(n_elem:n_elem-n_save+1:-1))

                    call snowpack_get_all_element_theta_water_pref(stations_in(i,j)%cxxmem%addr, tmp_arr, n_elem)
                    if (n_save > 0) Vol_Frac_WP(i,1:n_save,j) = real(tmp_arr(n_elem:n_elem-n_save+1:-1))

                    call snowpack_get_all_element_rg(stations_in(i,j)%cxxmem%addr, tmp_arr, n_elem)
                    if (n_save > 0) Rg(i,1:n_save,j) = real(tmp_arr(n_elem:n_elem-n_save+1:-1))

                    call snowpack_get_all_element_rb(stations_in(i,j)%cxxmem%addr, tmp_arr, n_elem)
                    if (n_save > 0) Rb(i,1:n_save,j) = real(tmp_arr(n_elem:n_elem-n_save+1:-1))

                    call snowpack_get_all_element_dd(stations_in(i,j)%cxxmem%addr, tmp_arr, n_elem)
                    if (n_save > 0) Dd(i,1:n_save,j) = real(tmp_arr(n_elem:n_elem-n_save+1:-1))

                    call snowpack_get_all_element_sp(stations_in(i,j)%cxxmem%addr, tmp_arr, n_elem)
                    if (n_save > 0) Sp(i,1:n_save,j) = real(tmp_arr(n_elem:n_elem-n_save+1:-1))

                    call snowpack_get_all_element_mk(stations_in(i,j)%cxxmem%addr, tmp_mk, n_elem)
                    if (n_save > 0) mk(i,1:n_save,j) = real(tmp_mk(n_elem:n_elem-n_save+1:-1))

                    ! TODO: FIX FOR LATER ! call snowpack_get_all_element_mass_hoar(stations_in(i,j)%cxxmem%addr, tmp_arr, n_elem)
                    ! TODO: FIX FOR LATER ! mass_hoar(i,1:n_save,j) = real(tmp_arr(n_elem:n_elem-n_save+1:-1))

                    call snowpack_get_all_element_CDot(stations_in(i,j)%cxxmem%addr, tmp_arr, n_elem)
                    if (n_save > 0) CDot(i,1:n_save,j) = real(tmp_arr(n_elem:n_elem-n_save+1:-1))

                    ! call snowpack_get_all_element_metamo(stations_in(i,j)%cxxmem%addr, tmp_arr, n_elem)
                    ! if (n_save > 0) metamo(i,1:n_save,j) = real(tmp_arr(n_elem:n_elem-n_save+1:-1))

                    call snowpack_get_all_element_N3(stations_in(i,j)%cxxmem%addr, tmp_arr, n_elem)
                    if (n_save > 0) N3(i,1:n_save,j) = real(tmp_arr(n_elem:n_elem-n_save+1:-1))
                else
                    n_save = 0
                    n_save_nodes = 0
                endif
                ! Zero out unused layers above n_save
                T_node(i,n_save_nodes+1:,j) = 0.0
                depositionDate(i,n_save+1:,j) = 0.0
                T_elem(i,n_save+1:,j) = 0.0
                Vol_Frac_I(i,n_save+1:,j) = 0.0
                Vol_Frac_W(i,n_save+1:,j) = 0.0
                Vol_Frac_A(i,n_save+1:,j) = 0.0
                Vol_Frac_S(i,n_save+1:,j) = 0.0
                Vol_Frac_WP(i,n_save+1:,j) = 0.0
                Rg(i,n_save+1:,j) = 0.0
                Rb(i,n_save+1:,j) = 0.0
                Dd(i,n_save+1:,j) = 0.0
                Sp(i,n_save+1:,j) = 0.0
                mk(i,n_save+1:,j) = 0.0
                CDot(i,n_save+1:,j) = 0.0
                snow_stress(i,n_save+1:,j) = 0.0
                N3(i,n_save+1:,j) = 0.0
            enddo
        enddo

        end associate

    end subroutine SnowStationsOut

    subroutine update_run_snowpack_flag(domain,current_snow)
        implicit none
        type(domain_t), intent(in) :: domain
        real, optional, allocatable, intent(in) :: current_snow(:, :)

        if (.not.allocated(run_snowpack_flag)) allocate(run_snowpack_flag(its:ite,jts:jte))

        ! Update module-level flag for if we want to run SNOWPACK. Condition is: 
        ! - 1. is there already snow on the ground (snow_height > 0)
        ! - 2. is there new snowfall (current_snow > 0)
        ! - 3. is it not a water pixel (land_surface_type != water)

        if (present(current_snow)) then
            run_snowpack_flag(its:ite,jts:jte) = ((domain%vars_2d(domain%var_indx(kVARS%snow_height)%v)%data_2d(its:ite,jts:jte) > 0.0) .or. &
                                (current_snow(its:ite,jts:jte) > 0.0)) .and. &
                                (domain%vars_2d(domain%var_indx(kVARS%land_mask)%v)%data_2di(its:ite,jts:jte)/=kLC_WATER)
        else
            run_snowpack_flag(its:ite,jts:jte) = (domain%vars_2d(domain%var_indx(kVARS%snow_height)%v)%data_2d(its:ite,jts:jte) > 0.0) .and. &
                                (domain%vars_2d(domain%var_indx(kVARS%land_mask)%v)%data_2di(its:ite,jts:jte)/=kLC_WATER)
        endif

    end subroutine

    !>--------------------------------- read_snowpack_state --------------------------------
    ! Reads SNOWPACK state variables from the domain file into domain arrays.
    ! Variables not found in the file are silently skipped (domain keeps defaults).
    ! Called before SnowStationsIn so that stations get initialized from file data.
    !>--------------------------------------------------------------------------------------
    subroutine read_snowpack_state(domain, options)
        implicit none
        type(domain_t), intent(inout) :: domain
        type(options_t), intent(in)   :: options

        character(len=kMAX_FILE_LENGTH) :: filename
        integer :: nz
        
        filename = options%domain%init_conditions_file
        if (STD_OUT_PE) write(*,*) "Reading SNOWPACK initial state from: ", trim(filename)

        ! 2D: snow_nlayers (always read — gate was checked by caller)
        call try_read_snp_var(filename, options%domain%snowpack_nlayers_var, domain, kVARS%snow_nlayers)



        call try_read_snp_var(filename, options%domain%snowpack_vfi_var, domain, kVARS%Vol_Frac_I)
        call try_read_snp_var(filename, options%domain%snowpack_vfw_var, domain, kVARS%Vol_Frac_W)
        call try_read_snp_var(filename, options%domain%snowpack_vfa_var, domain, kVARS%Vol_Frac_A)
        call try_read_snp_var(filename, options%domain%snowpack_vfs_var, domain, kVARS%Vol_Frac_S)
        call try_read_snp_var(filename, options%domain%snowpack_vfwp_var, domain, kVARS%Vol_Frac_WP)
        call try_read_snp_var(filename, options%domain%snowpack_deposition_var, domain, kVARS%depositionDate)
        call try_read_snp_var(filename, options%domain%snowpack_ds_var, domain, kVARS%Ds)
        call try_read_snp_var(filename, options%domain%snowpack_tsnow_var, domain, kVARS%snow_temperature)
        call try_read_snp_var(filename, options%domain%snowpack_tsnow_i_var, domain, kVARS%snow_temperature_i)
        call try_read_snp_var(filename, options%domain%snowpack_rg_var, domain, kVARS%Rg)
        call try_read_snp_var(filename, options%domain%snowpack_rb_var, domain, kVARS%Rb)
        call try_read_snp_var(filename, options%domain%snowpack_dd_var, domain, kVARS%Dd)
        call try_read_snp_var(filename, options%domain%snowpack_sp_var, domain, kVARS%Sp)
        call try_read_snp_var(filename, options%domain%snowpack_mk_var, domain, kVARS%mk)
        call try_read_snp_var(filename, options%domain%snowpack_cdot_var, domain, kVARS%CDot)
        ! call try_read_snp_var(filename, options%domain%snowpack_snow_stress_var, domain, kVARS%snow_stress)

        domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_I)%v)%data_3d = max(min(domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_I)%v)%data_3d,1.0),0.0)
        domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_W)%v)%data_3d = max(min(domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_W)%v)%data_3d,1.0),0.0)
        domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_A)%v)%data_3d = max(min(domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_A)%v)%data_3d,1.0),0.0)
        domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_S)%v)%data_3d = max(min(domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_S)%v)%data_3d,1.0),0.0)
        domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_WP)%v)%data_3d = max(min(domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_WP)%v)%data_3d,1.0),0.0)

        if (options%domain%snowpack_tsnow_i_var == "" .and. options%domain%snowpack_nlayers_var /= "") then
            if (STD_OUT_PE) write(*,*) "Snow Temperature on interface not provided, filling in with snow temperature at elements..."
            do j = jms, jme
                do i = ims, ime
                    nz = domain%vars_2d(domain%var_indx(kVARS%snow_nlayers)%v)%data_2di(i,j)
                    domain%vars_3d(domain%var_indx(kVARS%snow_temperature_i)%v)%data_3d(i,2:nz,j) = 0.5 * (domain%vars_3d(domain%var_indx(kVARS%snow_temperature)%v)%data_3d(i,1:nz-1,j) + &
                                                                                                            domain%vars_3d(domain%var_indx(kVARS%snow_temperature)%v)%data_3d(i,2:nz,j))
                    domain%vars_3d(domain%var_indx(kVARS%snow_temperature_i)%v)%data_3d(i,nz+1,j) = domain%vars_3d(domain%var_indx(kVARS%snow_temperature)%v)%data_3d(i,nz,j) 
                    domain%vars_3d(domain%var_indx(kVARS%snow_temperature_i)%v)%data_3d(i,1,j) = domain%vars_3d(domain%var_indx(kVARS%snow_temperature)%v)%data_3d(i,1,j) 
                enddo
            enddo
        endif

        do j = jms, jme
            do i = ims, ime
                if (domain%vars_2d(domain%var_indx(kVARS%snow_nlayers)%v)%data_2di(i,j) > 0) then
                    domain%vars_2d(domain%var_indx(kVARS%snow_height)%v)%data_2d(i,j) = &
                        sum(domain%vars_3d(domain%var_indx(kVARS%Ds)%v)%data_3d(i,1:domain%vars_2d(domain%var_indx(kVARS%snow_nlayers)%v)%data_2di(i,j),j))
                else
                    domain%vars_2d(domain%var_indx(kVARS%snow_height)%v)%data_2d(i,j) = 0.0
                endif
            enddo
        enddo

        if (STD_OUT_PE) write(*,*) "Done reading SNOWPACK initial state"

    end subroutine read_snowpack_state


    subroutine try_read_snp_var(filename, varname, domain, kvar)
        implicit none
        character(len=*), intent(in) :: filename, varname
        type(domain_t), intent(inout) :: domain
        integer, intent(in) :: kvar

        type(meta_data_t) :: var

        var = get_varmeta(kvar)

        if (varname /= "") then
            if (var%two_d) then
                if (var%dtype==kINTEGER) then
                    call try_read_snp_2d_int(filename, varname, domain, kvar)
                elseif (var%dtype==kREAL) then
                    ! not yet used
                endif
            elseif (var%three_d) then
                call try_read_snp_3d(filename, varname, domain, kvar)
            endif
        else
            if (STD_OUT_PE) write(*,*) "-------------------------------------------------------------------------"
            if (STD_OUT_PE) write(*,*) "WARNING: SNOWPACK static data partial, data for snowpack variable: ",trim(var%name)
            if (STD_OUT_PE) write(*,*) "WARNING: not provided. Behavior undefined. Setting to 0 and continuing"
            if (STD_OUT_PE) write(*,*) "-------------------------------------------------------------------------"
            if (var%two_d) then
                if (var%dtype==kINTEGER) then
                    domain%vars_2d(domain%var_indx(kvar)%v)%data_2di = 0
                elseif (var%dtype==kREAL) then
                    domain%vars_2d(domain%var_indx(kvar)%v)%data_2d = 0.0
                endif
            elseif (var%three_d) then
                domain%vars_3d(domain%var_indx(kvar)%v)%data_3d = 0.0
            endif
        endif   

    end subroutine try_read_snp_var

    !> Try to read a 2D integer variable from file into domain.
    !! Silently skips if the variable is not in the file or not allocated in domain.
    subroutine try_read_snp_2d_int(filename, varname, domain, kvar)
        implicit none
        character(len=*), intent(in) :: filename, varname
        type(domain_t), intent(inout) :: domain
        integer, intent(in) :: kvar

        real, allocatable :: temp_2d(:,:)
        integer :: ncid, varid, err

        if (domain%var_indx(kvar)%v <= 0) return

        err = nf90_open(filename, NF90_NOWRITE, ncid)
        if (err /= NF90_NOERR) return
        err = nf90_inq_varid(ncid, varname, varid)
        err =  nf90_close(ncid)

        if (err /= NF90_NOERR) then
            if (STD_OUT_PE) write(*,*) "  SNOWPACK init: variable '", trim(varname), "' not found, skipping"
            return
        endif

        call io_read(filename, varname, temp_2d)
        domain%vars_2d(domain%var_indx(kvar)%v)%data_2di(ims:ime,jms:jme) = &
            int(temp_2d(ims:ime,jms:jme))
        if (STD_OUT_PE) write(*,*) "  SNOWPACK init: read '", trim(varname), "'"

    end subroutine try_read_snp_2d_int


    !> Try to read a 3D variable from file into domain.
    !! Silently skips if the variable is not in the file or not allocated in domain.
    subroutine try_read_snp_3d(filename, varname, domain, kvar)
        implicit none
        character(len=*), intent(in) :: filename, varname
        type(domain_t), intent(inout) :: domain
        integer, intent(in) :: kvar

        real, allocatable :: temp_3d(:,:,:), data3d(:,:,:)
        integer :: ncid, varid, err, nx_file, ny_file, nz_file, nz_dom
        integer :: i, j, n_layer, no_attribute
        logical :: layer_1_top

        if (domain%var_indx(kvar)%v <= 0) return

        err = nf90_open(filename, NF90_NOWRITE, ncid)
        if (err /= NF90_NOERR) return
        err = nf90_inq_varid(ncid, varname, varid)
        err = nf90_close(ncid)

        if (err /= NF90_NOERR) then
            if (STD_OUT_PE) write(*,*) "  SNOWPACK init: variable '", trim(varname), "' not found, skipping"
            return
        endif

        call io_read(filename, varname, temp_3d)

        ! Map file data (x, z, y) into domain array (x, z, y)
        nz_dom = size(domain%vars_3d(domain%var_indx(kvar)%v)%data_3d, 2)
        nx_file = size(temp_3d,1)
        ny_file = size(temp_3d,2)
        nz_file = size(temp_3d,3)

        temp_3d = reshape(temp_3d, shape=[nx_file,nz_file,ny_file], order=[1,3,2])

        domain%vars_3d(domain%var_indx(kvar)%v)%data_3d(ims:ime, 1:min(nz_file,nz_dom), jms:jme) = &
            temp_3d(ims:ime, 1:min(nz_file,nz_dom), jms:jme)

        layer_1_top = .True.
        layer_1_top = io_var_reversed(filename, varname, no_attribute)

        if ((.not. layer_1_top) .and. (no_attribute == 0)) then
            allocate(data3d(ims:ime,1:min(nz_file,nz_dom),jms:jme))
            data3d = domain%vars_3d(domain%var_indx(kvar)%v)%data_3d
            domain%vars_3d(domain%var_indx(kvar)%v)%data_3d = -999
            do j = jms, jme
            do i = ims, ime
                n_layer = domain%vars_2d(domain%var_indx(kVARS%snow_nlayers)%v)%data_2di(i,j)
                domain%vars_3d(domain%var_indx(kvar)%v)%data_3d(i,n_layer:1:-1,j) = data3d(i,1:n_layer,j)
            enddo
            enddo
            deallocate(data3d)
        endif

        if (STD_OUT_PE) write(*,*) "  SNOWPACK init: read '", trim(varname), "'"

    end subroutine try_read_snp_3d

end module module_sm_SNOWPACKdrv