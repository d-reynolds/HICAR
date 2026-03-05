!!
!!----------------------------------------------------------
module module_sm_SNOWPACKdrv
    use time_object,         only : Time_type
    use mod_wrf_constants,   only : piconst, XLS
    use icar_constants
    use iso_fortran_env,       only : output_unit
    use options_interface,   only : options_t
    use variable_interface,  only : variable_t
    use domain_interface,    only : domain_t
    use iso_c_binding
    use meteoio_mod
    use snowpack_mod
    use snowpack_vector_helpers
    use namelist_utils, only : translate_numeric_mapping
    use mod_atm_utilities, only : relative_humidity
    use string, only : str, to_upper

    implicit none

    type(snowpack_config) :: snp_cfg
    type(snow_station), allocatable :: stations(:,:)
    logical, allocatable :: run_snowpack_flag(:,:) 

    integer :: ims, ime, jms, jme, its, ite, jts, jte
    integer :: i, j, k

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
        call cfg%add_key("CALCULATION_STEP_LENGTH", "Snowpack", trim(str(options%lsm%update_interval)))
        call cfg%add_key("METEO_STEP_LENGTH", "Snowpack", trim(str(options%lsm%update_interval)))
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
        call cfg%add_key("REDUCE_N_ELEMENTS", "SnowpackAdvanced", str(options%sm%snowpack_reduce_n_elements,c_bool=.True.))
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
            enddo
        enddo

        ! TODO: temporary auto-init snowpack to 15cm everywhere for testing
        if (.not.(options%restart%restart .or. context_change)) then
            ! ---------------------------- Setup 2D grid of snow stations ----------------------------
            do j = jts, jte
                do i = its, ite
                    
                    ! cycle if this is a water pixel
                    if (.not.run_snowpack_flag(i,j)) cycle
                    ! Station-level properties (must be set after init)
                    call stations(i,j)%set_cos_sl(real(cos(domain%vars_2d(domain%var_indx(kVARS%slope_angle)%v)%data_2d(i,j)),kind=8))       ! cos(slope angle) - flat ground
                    call stations(i,j)%set_sector(0_c_size_t)  ! Slope sector
                    call stations(i,j)%set_p_albedo(0.2d0)     ! Parameterized albedo (bare soil)
                    call stations(i,j)%set_albedo(0.2d0)       ! Surface albedo
                    call stations(i,j)%set_soil_alb(0.2d0)     ! Soil albedo - typical bare soil
                    call stations(i,j)%set_soil_emissivity(0.98d0) ! Soil emissivity
                    call stations(i,j)%set_bare_soil_z0(0.002d0)   ! Bare soil roughness [m]
                    call stations(i,j)%set_soil_node(0_c_size_t)   ! No soil layers

                    ! 7. ElementData - create element and add to station
                    do k = 1, kSNOW_GRID_Z
                        elem_id = k
                        elem = element_data(elem_id)
                        call elem%set_l(0.05d0)
                        call elem%set_te(268.0d0)
                        call elem%set_rg(0.5d0)
                        call elem%set_rb(0.25d0)
                        call elem%set_dd(0.3d0)
                        call elem%set_sp(0.7d0)

                        ! Set theta values: [SOIL, ICE, WATER, WATER_PREF, AIR]
                        call elem%set_theta([0.0d0, 0.5d0, 0.1d0, 0.0d0, 0.4d0])

                        ! Set deposition date (Julian day 2460310 ≈ Jan 1, 2024)
                        call elem%set_deposition_date_julian(2460310.0d0)
                        
                        ! Add element to station using the new addElement method
                        call stations(i,j)%add_element(elem)
                    enddo

                    call stations(i,j)%set_c_h(kSNOW_GRID_Z*0.05d0 )
                enddo
            enddo

            call SnowStationsOut(domain, stations)
        else
            call update_run_snowpack_flag(domain)
                
            call SnowStationsInitRestart(domain, stations)

            call SnowStationsIn(domain, stations)

        endif

    end subroutine sm_SNOWPACK_init

    subroutine sm_SNOWPACK(domain,options,lsm_dt,current_rain,current_snow,windspd)
        ! SNOWPACK snow model step
        type(domain_t), intent(inout) :: domain
        type(options_t), intent(in) :: options
        real, intent(in) :: lsm_dt
        real, allocatable, dimension(:,:), intent(in) :: current_rain, current_snow, windspd

        real(c_double) :: cumu_precip
        type(snowpack) :: snp
        type(bound_cond), allocatable :: bdata(:,:)
        type(surface_fluxes), allocatable :: sdata(:,:)
        type(current_meteo), allocatable :: meteo(:,:)

        real :: start_time, end_time, elapsed_time
        real(c_double) :: Te_arr(100), Rho_arr(100), L_arr(100)
        real(c_double) :: theta_ice_arr(100), theta_water_arr(100), theta_air_arr(100)
        real(c_double) :: T_arr(100), z_arr(100)

        ! Initialize cumulative precipitation
        cumu_precip = 0.0d0 ! dummy input var if forceaddsnowfall is set to true in the config options

        call update_run_snowpack_flag(domain,current_snow)

        ! start timing loop, get clock time
        call cpu_time(start_time)

        !! ------------------------------- copy existing snowpack state (TODO: future work, for now just use persistent station array with single nest) -------------------------------
        call MeteoIn(domain, current_rain, current_snow, windspd, meteo)
        call SurfaceFluxesIn(domain, sdata)
        call BoundaryConditionsIn(domain, bdata)
        ! call SnowStationsIn(domain, stations)

        ! Loop over horizontal grid points
        do j = jts, jte
            do i = its, ite
            !! ------------------------------- Run snowpack step -------------------------------
            if (run_snowpack_flag(i,j)) then
                snp = snowpack(snp_cfg)
                call snp%run_snowpack_model(meteo(i,j), stations(i,j), cumu_precip, bdata(i,j), sdata(i,j))
                call sdata(i,j)%collect_surface_fluxes(bdata(i,j), stations(i,j), meteo(i,j))
            endif
            end do
        end do

        !! ------------------------------- copy out new snowpack state -------------------------------
        call SurfaceFluxesOut(domain, sdata)
        call BoundaryConditionsOut(domain, bdata)
        call SnowStationsOut(domain, stations)

        ! end timing loop, get clock time
        call cpu_time(end_time)
        elapsed_time = end_time - start_time
        if (STD_OUT_PE) write(*,'(A,F8.3,A)') "SNOWPACK time elapsed: ", elapsed_time, " seconds"

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

        allocate(meteo(its:ite,jts:jte))

        associate(                                                                                           &
                    temperature => domain%vars_3d(domain%var_indx(kVARS%temperature)%v)%data_3d,             &
                    qv => domain%vars_3d(domain%var_indx(kVARS%water_vapor)%v)%data_3d, &
                    p => domain%vars_3d(domain%var_indx(kVARS%pressure)%v)%data_3d,                       &
                    roughness_z0 => domain%vars_2d(domain%var_indx(kVARS%roughness_z0)%v)%data_2d,       &
                    incoming_shortwave => domain%vars_2d(domain%var_indx(kVARS%shortwave)%v)%data_2d, &
                    net_longwave => domain%vars_2d(domain%var_indx(kVARS%longwave)%v)%data_2d,   &
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
                call meteo(i,j)%set_date_julian(2460310.0d0)  ! Date (Julian day 2460310 ≈ Jan 1, 2024)
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

                call meteo(i,j)%set_ea(-999.0d0)            ! Atmospheric emissivity - CRYOWRF Coupler: Mdata.ea, setting to -999 to force recalculation in SNOWPACK

                ! SET NET LW TO NET_LW FROM driving model...
                call meteo(i,j)%set_lw_net(real(net_longwave(i,j), kind=8))       ! Net longwave [W/m2] - computed from ilwr-outgoing SET to incoming longwave to 0.98*sb*temp^4
                ! ... SET TSS to the value used in HICAR to compute outgoing LW
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
                    sublimation => domain%vars_2d(domain%var_indx(kVARS%dSWE_subl)%v)%data_2d        &
                 )
        ! Loop over horizontal grid points
        do j = jts, jte
            do i = its, ite
                if (.not.(run_snowpack_flag(i,j))) cycle
                
                sublimation(i,j) = snowpack_get_mass_sublimation(sfluxes(i,j)%cxxmem%addr)

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
                    sensible_heat(i,j) = bconds(i,j)%get_qs()
                    latent_heat(i,j) = bconds(i,j)%get_ql()
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
        integer :: ne_arr(1:kSNOW_GRID_Z)

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
                    metamo => domain%vars_3d(domain%var_indx(kVARS%metamo)%v)%data_3d,                      &
                    N3 => domain%vars_3d(domain%var_indx(kVARS%N3)%v)%data_3d                        &
                 )

        do j = jts, jte
            do i = its, ite
                if (.not.(run_snowpack_flag(i,j))) cycle
                n_elem = snowpack_get_num_elements(stations_in(i,j)%cxxmem%addr)
                n_node = snowpack_get_num_nodes(stations_in(i,j)%cxxmem%addr)

                call snowpack_set_all_node_T(stations_in(i,j)%cxxmem%addr, real(T_node(i,1:n_node,j), kind=8), n_node)
                call snowpack_set_all_element_deposition_julian(stations_in(i,j)%cxxmem%addr, real(depositionDate(i,1:n_elem,j), kind=8), n_elem)
                call snowpack_set_all_element_L(stations_in(i,j)%cxxmem%addr, real(Layer_Thick(i,1:n_elem,j), kind=8), n_elem)
                call snowpack_set_all_element_Te(stations_in(i,j)%cxxmem%addr, real(T_elem(i,1:n_elem,j), kind=8), n_elem)
                call snowpack_set_all_element_theta_ice(stations_in(i,j)%cxxmem%addr, real(Vol_Frac_I(i,1:n_elem,j), kind=8), n_elem)
                call snowpack_set_all_element_theta_water(stations_in(i,j)%cxxmem%addr, real(Vol_Frac_W(i,1:n_elem,j), kind=8), n_elem)
                call snowpack_set_all_element_theta_air(stations_in(i,j)%cxxmem%addr, real(Vol_Frac_A(i,1:n_elem,j), kind=8), n_elem)
                call snowpack_set_all_element_theta_soil(stations_in(i,j)%cxxmem%addr, real(Vol_Frac_S(i,1:n_elem,j), kind=8), n_elem)
                call snowpack_set_all_element_theta_water_pref(stations_in(i,j)%cxxmem%addr, real(Vol_Frac_WP(i,1:n_elem,j), kind=8), n_elem)
                call snowpack_set_all_element_rg(stations_in(i,j)%cxxmem%addr, real(Rg(i,1:n_elem,j), kind=8), n_elem)
                call snowpack_set_all_element_rb(stations_in(i,j)%cxxmem%addr, real(Rb(i,1:n_elem,j), kind=8), n_elem)
                call snowpack_set_all_element_dd(stations_in(i,j)%cxxmem%addr, real(Dd(i,1:n_elem,j), kind=8), n_elem)
                call snowpack_set_all_element_sp(stations_in(i,j)%cxxmem%addr, real(Sp(i,1:n_elem,j), kind=8), n_elem)
                call snowpack_set_all_element_mk(stations_in(i,j)%cxxmem%addr, int(mk(i,1:n_elem,j), kind=c_short), n_elem)
                ! TODO: FIX FOR LATER ! call snowpack_set_all_element_mass_hoar(stations_in(i,j)%cxxmem%addr, real(mass_hoar(i,1:n_elem,j), kind=8), n_elem)
                ! call snowpack_set_all_element_ne(stations_in(i,j)%cxxmem%addr, ne_arr(1:n_elem), n_elem)
                call snowpack_set_all_element_CDot(stations_in(i,j)%cxxmem%addr, real(CDot(i,1:n_elem,j), kind=8), n_elem)
                call snowpack_set_all_element_metamo(stations_in(i,j)%cxxmem%addr, real(metamo(i,1:n_elem,j), kind=8), n_elem)
                call snowpack_set_all_element_N3(stations_in(i,j)%cxxmem%addr, real(N3(i,1:n_elem,j), kind=8), n_elem)

                call stations_in(i,j)%set_c_h(real(snow_height(i,j), kind=8))
                ! ! call stations_in(i,j)%set_swe(real(snow_water_equivalent(i,j), kind=8))
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
        real(c_double) :: tmp_arr(100)
        integer(c_short) :: tmp_mk(100)

        associate(                      &
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
                    metamo => domain%vars_3d(domain%var_indx(kVARS%metamo)%v)%data_3d,                      &
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
                    else
                        albedo(i,j) = stations_in(i,j)%get_soil_alb()
                        n_elem = 0
                        n_node = 0
                    endif

                    ! Now handle element(layer) and node(interface) properties
                    n_snow_layers(i,j) = n_elem


                    call snowpack_get_all_node_T(stations_in(i,j)%cxxmem%addr, tmp_arr, n_node)
                    T_node(i,1:n_node,j) = real(tmp_arr(1:n_node))
                    T_node(i,n_node+1:,j) = 0.0

                    call snowpack_get_all_element_deposition_julian(stations_in(i,j)%cxxmem%addr, tmp_arr, n_elem)
                    depositionDate(i,1:n_elem,j) = real(tmp_arr(1:n_elem))
                    depositionDate(i,n_elem+1:,j) = 0.0

                    call snowpack_get_all_element_L(stations_in(i,j)%cxxmem%addr, tmp_arr, n_elem)
                    Layer_Thick(i,1:n_elem,j) = real(tmp_arr(1:n_elem))
                    Layer_Thick(i,n_elem+1:,j) = 0.0

                    call snowpack_get_all_element_Te(stations_in(i,j)%cxxmem%addr, tmp_arr, n_elem)
                    T_elem(i,1:n_elem,j) = real(tmp_arr(1:n_elem))
                    T_elem(i,n_elem+1:,j) = 0.0

                    call snowpack_get_all_element_theta_ice(stations_in(i,j)%cxxmem%addr, tmp_arr, n_elem)
                    Vol_Frac_I(i,1:n_elem,j) = real(tmp_arr(1:n_elem))
                    Vol_Frac_I(i,n_elem+1:,j) = 0.0

                    call snowpack_get_all_element_theta_water(stations_in(i,j)%cxxmem%addr, tmp_arr, n_elem)
                    Vol_Frac_W(i,1:n_elem,j) = real(tmp_arr(1:n_elem))
                    Vol_Frac_W(i,n_elem+1:,j) = 0.0

                    call snowpack_get_all_element_theta_air(stations_in(i,j)%cxxmem%addr, tmp_arr, n_elem)
                    Vol_Frac_A(i,1:n_elem,j) = real(tmp_arr(1:n_elem))
                    Vol_Frac_A(i,n_elem+1:,j) = 0.0

                    call snowpack_get_all_element_theta_soil(stations_in(i,j)%cxxmem%addr, tmp_arr, n_elem)
                    Vol_Frac_S(i,1:n_elem,j) = real(tmp_arr(1:n_elem))
                    Vol_Frac_S(i,n_elem+1:,j) = 0.0

                    call snowpack_get_all_element_theta_water_pref(stations_in(i,j)%cxxmem%addr, tmp_arr, n_elem)
                    Vol_Frac_WP(i,1:n_elem,j) = real(tmp_arr(1:n_elem))
                    Vol_Frac_WP(i,n_elem+1:,j) = 0.0

                    call snowpack_get_all_element_rg(stations_in(i,j)%cxxmem%addr, tmp_arr, n_elem)
                    Rg(i,1:n_elem,j) = real(tmp_arr(1:n_elem))
                    Rg(i,n_elem+1:,j) = 0.0

                    call snowpack_get_all_element_rb(stations_in(i,j)%cxxmem%addr, tmp_arr, n_elem)
                    Rb(i,1:n_elem,j) = real(tmp_arr(1:n_elem))
                    Rb(i,n_elem+1:,j) = 0.0

                    call snowpack_get_all_element_dd(stations_in(i,j)%cxxmem%addr, tmp_arr, n_elem)
                    Dd(i,1:n_elem,j) = real(tmp_arr(1:n_elem))
                    Dd(i,n_elem+1:,j) = 0.0

                    call snowpack_get_all_element_sp(stations_in(i,j)%cxxmem%addr, tmp_arr, n_elem)
                    Sp(i,1:n_elem,j) = real(tmp_arr(1:n_elem))
                    Sp(i,n_elem+1:,j) = 0.0

                    call snowpack_get_all_element_mk(stations_in(i,j)%cxxmem%addr, tmp_mk, n_elem)
                    mk(i,1:n_elem,j) = real(tmp_mk(1:n_elem))
                    mk(i,n_elem+1:,j) = 0.0

                    ! TODO: FIX FOR LATER ! call snowpack_get_all_element_mass_hoar(stations_in(i,j)%cxxmem%addr, tmp_arr, n_elem)
                    ! TODO: FIX FOR LATER ! mass_hoar(i,1:n_elem,j) = real(tmp_arr(1:n_elem))

                    call snowpack_get_all_element_CDot(stations_in(i,j)%cxxmem%addr, tmp_arr, n_elem)
                    CDot(i,1:n_elem,j) = real(tmp_arr(1:n_elem))
                    CDot(i,n_elem+1:,j) = 0.0

                    call snowpack_get_all_element_metamo(stations_in(i,j)%cxxmem%addr, tmp_arr, n_elem)
                    metamo(i,1:n_elem,j) = real(tmp_arr(1:n_elem))
                    metamo(i,n_elem+1:,j) = 0.0

                    call snowpack_get_all_element_N3(stations_in(i,j)%cxxmem%addr, tmp_arr, n_elem)
                    N3(i,1:n_elem,j) = real(tmp_arr(1:n_elem))
                    N3(i,n_elem+1:,j) = 0.0

                endif
            enddo
        enddo

        end associate

    end subroutine SnowStationsOut

    subroutine SnowStationsInitRestart(domain, stations_in)
        implicit none
        type(domain_t), intent(in) :: domain
        type(snow_station), intent(inout) :: stations_in(its:ite,jts:jte)

        integer(c_short) :: elem_id
        type(element_data) :: elem

        associate(                      &
                    snow_height => domain%vars_2d(domain%var_indx(kVARS%snow_height)%v)%data_2d,               &
                    snow_water_equivalent => domain%vars_2d(domain%var_indx(kVARS%snow_water_equivalent)%v)%data_2d,                 &
                    albedo => domain%vars_2d(domain%var_indx(kVARS%albedo)%v)%data_2d,                 &
                    n_snow_layers => domain%vars_2d(domain%var_indx(kVARS%snow_nlayers)%v)%data_2di,                 &
                    depositionDate => domain%vars_3d(domain%var_indx(kVARS%depositionDate)%v)%data_3d       &
            )

            do j = jts, jte
                do i = its, ite
                    
                    ! cycle if this pixel is not being run
                    if (.not.run_snowpack_flag(i,j)) cycle

                    ! Station-level properties (must be set after init)
                    call stations_in(i,j)%set_cos_sl(real(cos(domain%vars_2d(domain%var_indx(kVARS%slope_angle)%v)%data_2d(i,j)), kind=8))       ! cos(slope angle) - flat ground
                    call stations_in(i,j)%set_sector(0_c_size_t)  ! Slope sector
                    call stations_in(i,j)%set_p_albedo(0.2d0)     ! Parameterized albedo (bare soil)
                    call stations_in(i,j)%set_albedo(real(albedo(i,j), kind=8))       ! Surface albedo
                    call stations_in(i,j)%set_soil_alb(0.2d0)     ! Soil albedo - typical bare soil
                    call stations_in(i,j)%set_soil_emissivity(0.98d0) ! Soil emissivity
                    call stations_in(i,j)%set_bare_soil_z0(0.002d0)   ! Bare soil roughness [m]
                    call stations_in(i,j)%set_soil_node(0_c_size_t)   ! No soil layers

                    do k = 1, n_snow_layers(i,j)
                        elem_id = k
                        elem = element_data(elem_id)               

                        ! Add element to station using the new addElement method
                        call stations_in(i,j)%add_element(elem)
                    enddo

                enddo
            enddo

            end associate

    end subroutine SnowStationsInitRestart

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

end module module_sm_SNOWPACKdrv