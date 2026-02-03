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

    implicit none

    type(config) :: cfg
    type(snowpack_config) :: snp_cfg
    type(current_meteo) :: meteo
    type(snow_station) :: station
    type(bound_cond) :: bdata
    type(surface_fluxes) :: sdata
    type(node_data) :: node
    type(element_data) :: elem
    integer(c_short) :: elem_id
    integer(c_size_t) :: n_elems, n_nodes

    integer :: ims, ime, jms, jme

    private
    public :: sm_SNOWPACK_init,sm_SNOWPACK
   
    contains

    subroutine sm_SNOWPACK_init(domain,options, context_chng)
        ! Initialize SNOWPACK variables
        type(domain_t), intent(inout) :: domain
        type(options_t), intent(in) :: options
        logical, optional, intent(in) :: context_chng
        integer :: i,j,k, hj, hi, i_s, i_e, j_s, j_e
        logical :: context_change
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

        ! Define SNOWPACK specific variables
        ! 1. Configuration
        cfg = config()

        ! Snowpack section
        call cfg%add_key("CALCULATION_STEP_LENGTH", "Snowpack", "900")
        call cfg%add_key("HEIGHT_OF_WIND_VALUE", "Snowpack", "10.0")
        call cfg%add_key("HEIGHT_OF_METEO_VALUES", "Snowpack", "2.0")
        call cfg%add_key("ENFORCE_MEASURED_SNOW_HEIGHTS", "Snowpack", "false")
        call cfg%add_key("SW_MODE", "Snowpack", "INCOMING")
        call cfg%add_key("ATMOSPHERIC_STABILITY", "Snowpack", "MO_MICHLMAYR")
        call cfg%add_key("ROUGHNESS_LENGTH", "Snowpack", "0.002")
        call cfg%add_key("CHANGE_BC", "Snowpack", "false")
        call cfg%add_key("THRESH_CHANGE_BC", "Snowpack", "-1.0")
        call cfg%add_key("SNP_SOIL", "Snowpack", "false")
        call cfg%add_key("SOIL_FLUX", "Snowpack", "false")
        call cfg%add_key("GEO_HEAT", "Snowpack", "0.06")
        call cfg%add_key("CANOPY", "Snowpack", "false")
        call cfg%add_key("MEAS_TSS", "Snowpack", "false")

        ! SnowpackAdvanced section
        call cfg%add_key("VARIANT", "SnowpackAdvanced", "DEFAULT")
        call cfg%add_key("ALPINE3D", "SnowpackAdvanced", "false")
        call cfg%add_key("SNOW_EROSION", "SnowpackAdvanced", "false")
        call cfg%add_key("SNOW_REDISTRIBUTION", "SnowpackAdvanced", "false")
        call cfg%add_key("COUPLEDPHASECHANGES", "SnowpackAdvanced", "true")
        call cfg%add_key("HN_DENSITY", "SnowpackAdvanced", "PARAMETERIZED")
        call cfg%add_key("HN_DENSITY_PARAMETERIZATION", "SnowpackAdvanced", "LEHNING_NEW")
        call cfg%add_key("HN_DENSITY_FIXEDVALUE", "SnowpackAdvanced", "100.0")
        call cfg%add_key("RIME_INDEX", "SnowpackAdvanced", "false")
        call cfg%add_key("NEWSNOW_LWC", "SnowpackAdvanced", "0.0")
        call cfg%add_key("READ_DSM", "SnowpackAdvanced", "false")
        call cfg%add_key("SNOW_ALBEDO", "SnowpackAdvanced", "PARAMETERIZED")
        call cfg%add_key("ALBEDO_PARAMETERIZATION", "SnowpackAdvanced", "LEHNING_2")
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
        call cfg%add_key("MINIMUM_L_ELEMENT", "SnowpackAdvanced", "0.0025")
        call cfg%add_key("HEIGHT_NEW_ELEM", "SnowpackAdvanced", "0.02")
        call cfg%add_key("COMBINE_ELEMENTS", "SnowpackAdvanced", "true")
        call cfg%add_key("REDUCE_N_ELEMENTS", "SnowpackAdvanced", "0")
        call cfg%add_key("FIXED_POSITIONS", "SnowpackAdvanced", "0.25 0.5 1.0")
        call cfg%add_key("NUMBER_FIXED_RATES", "SnowpackAdvanced", "0")
        call cfg%add_key("MAX_NUMBER_MEAS_TEMPERATURES", "SnowpackAdvanced", "5")
        call cfg%add_key("MIN_DEPTH_SUBSURF", "SnowpackAdvanced", "0.0")

        snp_cfg = snowpack_config(cfg)


    end subroutine sm_SNOWPACK_init

    subroutine sm_SNOWPACK(domain,options,lsm_dt)
        ! SNOWPACK snow model step
        type(domain_t), intent(inout) :: domain
        type(options_t), intent(in) :: options
        real, intent(in) :: lsm_dt

        ! Local variables
        integer :: i,j,k,nz

        ! 2. CurrentMeteo - meteorological forcing
        meteo = current_meteo()
        call meteo%set_ta(268.15d0)
        call meteo%set_rh(0.75d0)
        call meteo%set_vw(5.0d0)
        call meteo%set_dw(225.0d0)
        call meteo%set_iswr(150.0d0)
        call meteo%set_psum(0.5d0)
        call meteo%set_psum_ph(0.0d0)
        call meteo%set_tss(271.0d0)
        write(*,'(A,F6.2,A,F4.1,A)') "CurrentMeteo:  Ta=", meteo%get_ta(), " K, VW=", meteo%get_vw(), " m/s"

        ! 3. SnowStation - snowpack state
        station = snow_station(.false., .false.)
        call station%set_cos_sl(1.0d0)
        call station%set_albedo(0.82d0)
        call station%set_c_h(0.5d0)
        call station%set_m_h(0.48d0)
        n_elems = station%get_number_of_elements()
        n_nodes = station%get_number_of_nodes()
        write(*,'(A,F5.2,A,I0)') "SnowStation:   cH=", station%get_c_h(), " m, elements=", n_elems

        ! 4. BoundCond - boundary conditions
        bdata = bound_cond()
        call bdata%set_lw_out(315.0d0)
        call bdata%set_qs(25.0d0)
        call bdata%set_ql(10.0d0)
        write(*,'(A,F5.1,A)') "BoundCond:     qs=", bdata%get_qs(), " W/m2"

        ! 5. SurfaceFluxes - output fluxes
        sdata = surface_fluxes()
        call sdata%set_sw_in(150.0d0)
        call sdata%set_sw_out(120.0d0)
        call sdata%set_lw_in(265.0d0)
        call sdata%set_lw_out(315.0d0)
        write(*,'(A,F5.1,A)') "SurfaceFluxes: SW_in=", sdata%get_sw_in(), " W/m2"

        ! 6. NodeData - standalone layer node
        node = node_data()
        call node%set_z(0.25d0)
        call node%set_t(270.0d0)
        call node%set_s_n(1.5d0)
        write(*,'(A,F5.2,A,F6.1,A,F4.1)') "NodeData:      z=", node%get_z(), " m, T=", node%get_t(), " K, S_n=", node%get_s_n()

        ! 7. ElementData - standalone layer element
        elem_id = 1
        elem = element_data(elem_id)
        call elem%set_l(0.05d0)
        call elem%set_te(268.0d0)
        call elem%set_rho(250.0d0)
        call elem%set_rg(0.5d0)
        call elem%set_dd(0.3d0)
        call elem%set_sp(0.7d0)
        write(*,'(A,F4.1,A,F5.1,A,F3.1,A,F3.1)') "ElementData:   L=", elem%get_l()*100, &
            " cm, Rho=", elem%get_rho(), " kg/m3, dd=", elem%get_dd(), ", sp=", elem%get_sp()


        ! ! Loop over horizontal grid points
        ! do j = domain%j_start, domain%j_end
        !     do i = domain%i_start, domain%i_end

        !         ! Retrieve SNOWPACK variables
        !     end do
        ! end do

    end subroutine sm_SNOWPACK

end module module_sm_SNOWPACKdrv