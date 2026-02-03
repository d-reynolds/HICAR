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

        allocate(stations(its:ite,jts:jte))
        ! ---------------------------- Setup 2D grid of snow stations ----------------------------
        do j = jts, jte
            do i = its, ite
                ! 3. SnowStation - bare ground initialization
                ! Following pattern from SnowStation::initialize() in DataClasses.cc
                !                if_use_canopy, if_use_soil
                stations(i,j) = snow_station(.false., .false.)

                call snowpack_station_resize(stations(i,j)%cxxmem%addr, 0_c_size_t)

                ! Station-level properties (must be set after init)
                call stations(i,j)%set_cos_sl(1.0d0)       ! cos(slope angle) - flat ground
                call stations(i,j)%set_sector(0_c_size_t)  ! Slope sector
                call stations(i,j)%set_p_albedo(0.2d0)     ! Parameterized albedo (bare soil)
                call stations(i,j)%set_albedo(0.2d0)       ! Surface albedo
                call stations(i,j)%set_soil_alb(0.2d0)     ! Soil albedo - typical bare soil
                call stations(i,j)%set_soil_emissivity(0.98d0) ! Soil emissivity
                call stations(i,j)%set_bare_soil_z0(0.002d0)   ! Bare soil roughness [m]
                call stations(i,j)%set_soil_node(0_c_size_t)   ! No soil layers

                ! 7. ElementData - create element and add to station
                do k = 1, 5
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
            enddo
        enddo

    end subroutine sm_SNOWPACK_init

    subroutine sm_SNOWPACK(domain,options,lsm_dt,current_rain,current_snow,windspd)
        ! SNOWPACK snow model step
        type(domain_t), intent(inout) :: domain
        type(options_t), intent(in) :: options
        real, intent(in) :: lsm_dt
        real, allocatable, dimension(:,:), intent(in) :: current_rain, current_snow, windspd

        real(c_double) :: cumu_precip
        type(snowpack) :: snp
        type(bound_cond) :: bdata
        type(surface_fluxes) :: sdata
        type(current_meteo) :: meteo

        real(kind=8), dimension(its:ite,jts:jte) :: total_precip, precip_phase, rh
        integer(c_size_t) :: n_elems, n_nodes
        real :: start_time, end_time, elapsed_time
        real(c_double) :: Te_arr(100), Rho_arr(100), L_arr(100)
        real(c_double) :: theta_ice_arr(100), theta_water_arr(100), theta_air_arr(100)
        real(c_double) :: T_arr(100), z_arr(100)

        associate(                                                                                           &
                    temperature => domain%vars_3d(domain%var_indx(kVARS%temperature)%v)%data_3d,             &
                    qv => domain%vars_3d(domain%var_indx(kVARS%water_vapor)%v)%data_3d, &
                    p => domain%vars_3d(domain%var_indx(kVARS%pressure)%v)%data_3d,                       &
                    roughness_z0 => domain%vars_2d(domain%var_indx(kVARS%roughness_z0)%v)%data_2d,       &
                    incoming_shortwave => domain%vars_2d(domain%var_indx(kVARS%shortwave)%v)%data_2d, &
                    net_longwave => domain%vars_2d(domain%var_indx(kVARS%longwave)%v)%data_2d,   &
                    ground_temp => domain%vars_3d(domain%var_indx(kVARS%soil_temperature)%v)%data_3d,               &
                    skin_temperature => domain%vars_2d(domain%var_indx(kVARS%skin_temperature)%v)%data_2d,              &
                    sensible_heat => domain%vars_2d(domain%var_indx(kVARS%sensible_heat)%v)%data_2d,               &
                    latent_heat => domain%vars_2d(domain%var_indx(kVARS%latent_heat)%v)%data_2d                 &
                 )

        total_precip = current_rain(its:ite,jts:jte) + current_snow(its:ite,jts:jte)
        precip_phase = current_rain(its:ite,jts:jte) / (total_precip+1.0d-6)  ! avoid div by zero
        rh = relative_humidity(temperature(its:ite,1,jts:jte),qv(its:ite,1,jts:jte),p(its:ite,1,jts:jte))

        ! Initialize cumulative precipitation
        cumu_precip = 0.0d0 ! dummy input var if forceaddsnowfall is set to true in the config options

        ! start timing loop, get clock time
        call cpu_time(start_time)

        ! Loop over horizontal grid points
        do j = jts, jte
            do i = its, ite
            !! ------------------------------- Initialize new in/out structures (should be deprecated to initialize one outside of loop) -------------------------------
            bdata = bound_cond()
            sdata = surface_fluxes()
            meteo = current_meteo()

            !! ------------------------------- Set meteo variables -------------------------------
            call meteo%set_date_julian(2460310.0d0)  ! Date (Julian day 2460310 ≈ Jan 1, 2024)
            call meteo%set_ta(real(temperature(i,1,j), kind=8))          ! Air temperature [K] - CRYOWRF Coupler: l_TA
            call meteo%set_rh(rh(i,j))            ! Relative humidity [0-1] - CRYOWRF Coupler: l_RH
            call meteo%set_rh_avg(rh(i,j))        ! Running mean RH, set same as RH
            call meteo%set_vw(real(windspd(i,j), kind=8))             ! Wind velocity [m/s] - CRYOWRF Coupler: l_VW
            ! call meteo%set_vw_avg(5.0d0)         ! Running mean wind velocity
            ! call meteo%set_vw_max(8.0d0)         ! Max wind velocity - CRYOWRF Coupler: l_VW_MAX
            call meteo%set_vw_drift(real(windspd(i,j), kind=8))       ! Drift wind velocity
            ! call meteo%set_ustar(0.3d0)          ! Friction velocity [m/s]
            call meteo%set_z0(real(roughness_z0(i,j), kind=8))           ! Roughness length [m] - matches ROUGHNESS_LENGTH config -- CURRENTLY OVERWRITTEN IN SNOWPACK
            call meteo%set_iswr(real(incoming_shortwave(i,j), kind=8))         ! Incoming shortwave [W/m2] - CRYOWRF Coupler: l_iswr

            call meteo%set_ea(-999.0d0)            ! Atmospheric emissivity - CRYOWRF Coupler: Mdata.ea, setting to -999 to force recalculation in SNOWPACK

            ! SET NET LW TO NET_LW FROM driving model...
            call meteo%set_lw_net(real(net_longwave(i,j), kind=8))       ! Net longwave [W/m2] - computed from ilwr-outgoing SET to incoming longwave to 0.98*sb*temp^4
            ! ... SET TSS to the value used in HICAR to compute outgoing LW
            call meteo%set_tss(real(skin_temperature(i,j), kind=8))          ! Snow surface temperature [K]

            ! NOTE: because SOIL_FLUX is set to false above, ts0 has to be the ground surface temperature from NoahMP
            call meteo%set_ts0(real(ground_temp(i,1,j), kind=8))          ! Bottom temperature [K] - CRYOWRF Coupler: TSG
            call meteo%set_psum(total_precip(i,j))           ! Precipitation sum since previous SNOWPACK call [mm] - CRYOWRF Coupler: l_psum
            call meteo%set_psum_ph(precip_phase(i,j))        ! Precip phase (0=solid, 1=liquid)

            
            call meteo%set_dw(0.0d0)           ! Wind direction [deg] - CRYOWRF Coupler: l_DW
            call meteo%set_dw_drift(0.0d0)     ! Drift wind direction
            call meteo%set_psi_s(0.0d0)          ! Stability correction
            call meteo%set_hs(0.0d0)             ! Measured snow height [m] -- needed to run, but only used if running on IMIS data, which doesn't happen when coupled to HICAR
            call meteo%set_geo_heat(0.0d0)       ! Geothermal heat flux [W/m2] - CRYOWRF io.ini: GEO_HEAT
            ! if (STD_OUT_PE) write(*,'(A,F6.2,A,F4.1,A)') "CurrentMeteo:  Ta=", meteo%get_ta(), " K, VW=", meteo%get_vw(), " m/s"


            !! ------------------------------- copy existing snowpack state (TODO: future work, for now just use persistent station array with single nest) -------------------------------



            !! ------------------------------- Run snowpack step -------------------------------
            snp = snowpack(snp_cfg)
            if (i==its .and. j==jts .and. STD_OUT_PE) write(*,'(A)') "Calling run_snowpack_model..."
            call snp%run_snowpack_model(meteo, stations(i,j), cumu_precip, bdata, sdata)

            !! ------------------------------- copy out new snowpack state -------------------------------
            if (i==its .and. j==jts .and. STD_OUT_PE) write(*,'(A)') "Model step completed!"
            if (i==its .and. j==jts .and. STD_OUT_PE) write(*,'(A,F8.3)') "  cumu_precip = ", cumu_precip
            if (i==its .and. j==jts .and. STD_OUT_PE) write(*,'(A,F8.2)') "  qs from bdata (sensible) = ", bdata%get_qs()
            if (i==its .and. j==jts .and. STD_OUT_PE) write(*,'(A,F8.2)') "  ql (latent)   = ", bdata%get_ql()
            if (i==its .and. j==jts .and. STD_OUT_PE) write(*,'(A,F8.2)') "  lw_out = ", bdata%get_lw_out()
            if (i==its .and. j==jts .and. STD_OUT_PE) write(*,'(A,F8.2)') "  lw_out = ", stations(i,j)%get_c_h()
            sensible_heat(i,j) = bdata%get_qs()
            latent_heat(i,j) = bdata%get_ql()

            ! ! Get bulk element arrays
            ! call snowpack_get_all_element_Te(stations(i,j)%cxxmem%addr, Te_arr, n_elems)
            ! call snowpack_get_all_element_Rho(stations(i,j)%cxxmem%addr, Rho_arr, n_elems)
            ! call snowpack_get_all_element_L(stations(i,j)%cxxmem%addr, L_arr, n_elems)
            ! call snowpack_get_all_element_theta_ice(stations(i,j)%cxxmem%addr, theta_ice_arr, n_elems)
            ! call snowpack_get_all_element_theta_water(stations(i,j)%cxxmem%addr, theta_water_arr, n_elems)
            ! call snowpack_get_all_element_theta_air(stations(i,j)%cxxmem%addr, theta_air_arr, n_elems)

            ! if (STD_OUT_PE) write(*,'(A)') "Elements:"
            ! if (STD_OUT_PE) write(*,'(A)') "  Elem    L(cm)    Te(K)   Rho(kg/m3)  theta_i  theta_w  theta_a"
            ! if (STD_OUT_PE) write(*,'(A)') "  ----  -------  -------  ----------  -------  -------  -------"
            ! do i = 1, n_elems
            !     if (STD_OUT_PE) write(*,'(I6,F9.2,F9.2,F12.1,F9.3,F9.3,F9.3)') i, &
            !         L_arr(i)*100, Te_arr(i), Rho_arr(i), &
            !         theta_ice_arr(i), theta_water_arr(i), theta_air_arr(i)
            ! end do

            ! ! Get bulk node arrays
            ! n_nodes = snowpack_get_num_nodes(stations(i,j)%cxxmem%addr)
            ! call snowpack_get_all_node_T(stations(i,j)%cxxmem%addr, T_arr, n_nodes)
            ! call snowpack_get_all_node_z(stations(i,j)%cxxmem%addr, z_arr, n_nodes)

            ! if (STD_OUT_PE) write(*,'(A)') ""
            ! if (STD_OUT_PE) write(*,'(A)') "Nodes:"
            ! if (STD_OUT_PE) write(*,'(A)') "  Node     z(m)      T(K)"
            ! if (STD_OUT_PE) write(*,'(A)') "  ----  --------  --------"
            ! do i = 1, n_nodes
            !     if (STD_OUT_PE) write(*,'(I6,F10.4,F10.2)') i, z_arr(i), T_arr(i)
            ! end do
            end do
        end do

        end associate

        ! end timing loop, get clock time
        call cpu_time(end_time)
        elapsed_time = end_time - start_time
        if (STD_OUT_PE) write(*,'(A,F8.3,A)') "SNOWPACK time elapsed: ", elapsed_time, " seconds"

    end subroutine sm_SNOWPACK

end module module_sm_SNOWPACKdrv