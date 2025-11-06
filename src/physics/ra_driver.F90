!>----------------------------------------------------------
!! This module provides a wrapper to call various radiation models
!! It sets up variables specific to the physics package to be used
!!
!! The main entry point to the code is rad(domain,options,dt)
!!
!! <pre>
!! Call tree graph :
!!  radiation_init->[ external initialization routines]
!!  rad->[  external radiation routines]
!!
!! High level routine descriptions / purpose
!!   radiation_init     - initializes physics package
!!   rad                - sets up and calls main physics package
!!
!! Inputs: domain, options, dt
!!      domain,options  = as defined in data_structures
!!      dt              = time step (seconds)
!! </pre>
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!----------------------------------------------------------
module radiation
    use module_ra_simple, only: ra_simple
    use module_ra_rrtmg_lw, only: rrtmg_lwinit, rrtmg_lwrad
    use module_ra_rrtmg_sw, only: rrtmg_swinit, rrtmg_swrad
    use options_interface,  only : options_t
    use domain_interface,   only : domain_t
    use iso_fortran_env, only: real64
    use time_object,        only : Time_type
    use icar_constants, only : kVARS, kRA_BASIC, kRA_SIMPLE, kRA_RRTMG, kRA_RRTMGP, STD_OUT_PE, kMP_THOMP_AER, kMAX_NESTS
    use mod_wrf_constants, only : cp, R_d, gravity, DEGRAD, DPD, piconst
    use mod_atm_utilities, only : cal_cldfra3, calc_solar_elevation, calc_solar_date
    use mo_rte_kind,           only: wp, i8, wl
    use mo_rte_lw,             only: rte_lw
    use mo_rte_sw,             only: rte_sw
    use mo_aerosol_optics_rrtmgp_merra ! Includes aerosol type integers
    use mo_rte_config,         only: rte_config_checks
    use mo_optical_props,      only: ty_optical_props, &
                                   ty_optical_props_arry, ty_optical_props_1scl, ty_optical_props_2str
    use mo_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp
    use mo_cloud_optics_rrtmgp,only: ty_cloud_optics_rrtmgp
    use mo_gas_concentrations, only: ty_gas_concs
    use mo_source_functions,   only: ty_source_func_lw
    use mo_fluxes,             only: ty_fluxes_broadband
    use mo_load_coefficients,  only: load_and_init
    use mo_load_cloud_coefficients, &
                                only: load_cld_lutcoeff
    use mo_load_aerosol_coefficients, &
                                only: load_aero_lutcoeff
    implicit none

    !! following constants taken from icon source code on October 2nd 2025, release 2025.04
    !! ----------------------------------------------
    ! ICON
    !
    ! ---------------------------------------------------------------
    ! Copyright (C) 2004-2025, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
    ! Contact information: icon-model.org
    !
    ! See AUTHORS.TXT for a list of authors
    ! See LICENSES/ for license information
    ! SPDX-License-Identifier: BSD-3-Clause
    ! ---------------------------------------------------------------
    ! This module provides physical constants for the ICON general circulation models.
    !
    ! Physical constants are grouped as follows:
    ! - Natural constants
    ! - Molar weights
    ! - Earth and Earth orbit constants
    ! - Thermodynamic constants for the dry and moist atmosphere
    ! - Constants used for the computation of lookup tables of the saturation
    !    mixing ratio over liquid water (*c_les*) or ice(*c_ies*)
    !    (to be shifted to the module that computes the lookup tables)
    !> Molar weights
    !! -------------
    !!
    !! Pure species
    REAL(wp), PARAMETER :: amco2 = 44.011_wp        !>[g/mol] CO2
    REAL(wp), PARAMETER :: amch4 = 16.043_wp        !! [g/mol] CH4
    REAL(wp), PARAMETER :: amo3  = 47.9982_wp       !! [g/mol] O3
    REAL(wp), PARAMETER :: amo2  = 31.9988_wp       !! [g/mol] O2
    REAL(wp), PARAMETER :: amn2o = 44.013_wp        !! [g/mol] N2O
    REAL(wp), PARAMETER :: amc11 =137.3686_wp       !! [g/mol] CFC11
    REAL(wp), PARAMETER :: amc12 =120.9140_wp       !! [g/mol] CFC12
    REAL(wp), PARAMETER :: amw   = 18.0154_wp       !! [g/mol] H2O
    REAL(wp), PARAMETER :: amo   = 15.9994_wp       !! [g/mol] O
    REAL(wp), PARAMETER :: amno  = 30.0061398_wp    !! [g/mol] NO
    REAL(wp), PARAMETER :: amn2  = 28.0134_wp       !! [g/mol] N2
    REAL(wp), PARAMETER :: amso4 = 96.0626_wp       !! [g/mol] SO4
    REAL(wp), PARAMETER :: ams   = 32.06_wp         !! [g/mol] S
    !
    !> Mixed species
    REAL(wp), PARAMETER :: amd   = 28.970_wp        !> [g/mol] dry air



    integer :: update_interval
    real*8  :: last_model_time(kMAX_NESTS)
    real    :: solar_constant
    real    :: p_top = 100000.0
    
    !! MJ added to aggregate radiation over output interval
    real, allocatable, dimension(:,:) :: shortwave_cached, cos_project_angle, solar_elevation_store, solar_azimuth_store
    real*8 :: counter
    real*8  :: Delta_t !! MJ added to detect the time for outputting 
    integer :: ims, ime, jms, jme, kms, kme
    integer :: its, ite, jts, jte, kts, kte
    integer :: ids, ide, jds, jde, kds, kde

    type(ty_gas_optics_rrtmgp)   :: k_dist_lw, k_dist_sw
    type(ty_cloud_optics_rrtmgp) :: cloud_optics_lw, cloud_optics_sw
    type(ty_aerosol_optics_rrtmgp_merra)   &
                                    :: aerosol_optics_lw, aerosol_optics_sw
    logical :: do_aerosols = .False.
    integer, parameter :: ngas = 8
    integer :: nblocks = 1

    character(len=3), dimension(ngas), parameter :: &
                gas_names = ['h2o', 'o3', 'co2', 'n2o', 'co ', 'ch4', 'o2 ', 'n2 ']
    type(ty_gas_concs)           :: gas_concs

    private
    public :: radiation_init, ra_var_request, rad_apply_dtheta, rad
    
contains

    subroutine radiation_init(domain,options, context_chng)
        type(domain_t), intent(inout) :: domain
        type(options_t),intent(in) :: options
        logical, optional, intent(in) :: context_chng

        logical :: context_change
        integer :: i
        character(len=512) :: k_dist_sw_file = 'rrtmgp_support/rrtmgp-gas-sw-g112.nc'
        character(len=512) :: k_dist_lw_file = 'rrtmgp_support/rrtmgp-gas-lw-g128.nc'
        character(len=512) :: cloud_optics_sw_file = 'rrtmgp_support/rrtmgp-clouds-sw-bnd.nc'
        character(len=512) :: cloud_optics_lw_file = 'rrtmgp_support/rrtmgp-clouds-lw-bnd.nc'
        character(len=512) :: aerosol_optics_sw_file = ''
        character(len=512) :: aerosol_optics_lw_file = ''

        if (present(context_chng)) then
            context_change = context_chng
        else
            context_change = .false.
        endif
        
        ims = domain%grid%ims
        ime = domain%grid%ime
        jms = domain%grid%jms
        jme = domain%grid%jme
        kms = domain%grid%kms
        kme = domain%grid%kme
        its = domain%grid%its
        ite = domain%grid%ite
        jts = domain%grid%jts
        jte = domain%grid%jte
        kts = domain%grid%kts
        kte = domain%grid%kte

        ids = domain%grid%ids
        ide = domain%grid%ide
        jds = domain%grid%jds
        jde = domain%grid%jde
        kds = domain%grid%kds
        kde = domain%grid%kde
        
        if (options%physics%radiation == 0) return

        update_interval=options%rad%update_interval_rad ! 30 min, 1800 s   600 ! 10 min (600 s)

        !Saftey bound, in case update_interval is 0, or very small
        if (.not.(context_change)) then
            if (update_interval<=10) then
                last_model_time(domain%nest_indx) = domain%sim_time%seconds()-10
            else
                last_model_time(domain%nest_indx) = domain%sim_time%seconds()-update_interval
            endif
        endif
        
        ! if (options%physics%radiation_downScaling==1) then
            if (allocated(cos_project_angle)) then
                !$acc exit data delete(cos_project_angle)
                deallocate(cos_project_angle)
            endif

            allocate(cos_project_angle(ims:ime,jms:jme)) !! MJ added


            if (allocated(solar_elevation_store)) then
                !$acc exit data delete(solar_elevation_store)
                deallocate(solar_elevation_store)
            endif

            allocate(solar_elevation_store(ims:ime,jms:jme)) !! MJ added


            if (allocated(solar_azimuth_store)) then
                !$acc exit data delete(solar_azimuth_store)
                deallocate(solar_azimuth_store)
            endif

            allocate(solar_azimuth_store(ims:ime,jms:jme)) !! MJ added

            if (allocated(shortwave_cached)) then
                !$acc exit data delete(shortwave_cached)
                deallocate(shortwave_cached)
            endif

            allocate(shortwave_cached(ims:ime,jms:jme))

            !$acc enter data create(cos_project_angle, solar_elevation_store, solar_azimuth_store, shortwave_cached)
        ! endif
        
        ! If we are just changing nest contexts, we don't need to reinitialize the radiation modules
        ! if (context_change) return

        if (STD_OUT_PE .and. .not.context_change) write(*,*) "Initializing Radiation"

        if (options%physics%radiation==kRA_BASIC) then
            if (STD_OUT_PE .and. .not.context_change) write(*,*) "    Basic Radiation"
        endif
        if (options%physics%radiation==kRA_SIMPLE .or. (options%physics%landsurface>0 .and. options%physics%radiation==0)) then
            if (STD_OUT_PE .and. .not.context_change) write(*,*) "    Simple Radiation"
            !call ra_simple_init(domain)
        endif!! MJ added to detect the time for outputting 

        if (options%physics%radiation==kRA_RRTMG) then
            if (STD_OUT_PE .and. .not.context_change) write(*,*) "    RRTMG"

            if (options%physics%microphysics .ne. kMP_THOMP_AER) then
               if (STD_OUT_PE .and. .not.context_change)write(*,*) '    NOTE: When running RRTMG, microphysics option 5 works best.'
            endif

            !$acc update host(domain%vars_3d(domain%var_indx(kVARS%pressure_interface)%v)%data_3d)
            ! This will capture the highest pressure level of all nests in this simulation
            p_top = min(p_top, minval(domain%vars_3d(domain%var_indx(kVARS%pressure_interface)%v)%data_3d(domain%ims:domain%ime,domain%kme+1,domain%jms:domain%jme)))

            call rrtmg_lwinit(                           &
                p_top=p_top,     allowed_to_read=.TRUE. ,                &
                ids=domain%ids, ide=domain%ide, jds=domain%jds, jde=domain%jde, kds=domain%kds, kde=domain%kde,                &
                ims=domain%ims, ime=domain%ime, jms=domain%jms, jme=domain%jme, kms=domain%kms, kme=domain%kme,                &
                its=domain%its, ite=domain%ite, jts=domain%jts, jte=domain%jte, kts=domain%kts, kte=domain%kte                 )

            call rrtmg_swinit(                           &
                allowed_to_read=.TRUE.,                     &
                ids=domain%ids, ide=domain%ide, jds=domain%jds, jde=domain%jde, kds=domain%kds, kde=domain%kde,                &
                ims=domain%ims, ime=domain%ime, jms=domain%jms, jme=domain%jme, kms=domain%kms, kme=domain%kme,                &
                its=domain%its, ite=domain%ite, jts=domain%jts, jte=domain%jte, kts=domain%kts, kte=domain%kte                 )
                domain%tend%th_swrad = 0
                domain%tend%th_lwrad = 0
        else if (options%physics%radiation == kRA_RRTMGP) then
            if (STD_OUT_PE .and. .not.context_change) write(*,*) "    RRTMGP"    
            ! if (.not.(k_dist_sw%is_loaded())) then

                !determine number of blocks to use for RRTMGP. This is done to reduce
                ! the GPU memory requirements. A good base is 200x200 horizontal points per block
                nblocks = MAX(1, ((ite - its + 1)*(jte - jts + 1)) / 40000)
    
                if (k_dist_sw%is_loaded() .or. k_dist_lw%is_loaded()) then
                    call k_dist_sw%finalize()
                    call k_dist_lw%finalize()
                    call cloud_optics_sw%finalize()
                    call cloud_optics_lw%finalize()
                endif
                call rte_config_checks(logical(.false., wl))

                call set_atmo_gas_conc(domain)

                ! ----------------------------------------------------------------------------
                ! load data into classes
                call load_and_init(k_dist_sw, k_dist_sw_file, gas_concs)
                call load_and_init(k_dist_lw, k_dist_lw_file, gas_concs)

                !
                call load_cld_lutcoeff(cloud_optics_lw, cloud_optics_lw_file)
                call load_cld_lutcoeff(cloud_optics_sw, cloud_optics_sw_file)
                call stop_on_err(cloud_optics_lw%set_ice_roughness(2))
                call stop_on_err(cloud_optics_sw%set_ice_roughness(2))

                if (do_aerosols) then
                    ! Load aerosol optics coefficients from lookup tables
                    call load_aero_lutcoeff (aerosol_optics_lw, aerosol_optics_lw_file)
                    call load_aero_lutcoeff (aerosol_optics_sw, aerosol_optics_sw_file)
                end if
            ! endif
        endif

    end subroutine radiation_init


    subroutine ra_var_request(options)
        implicit none
        type(options_t), intent(inout) :: options

        if (options%physics%radiation == kRA_SIMPLE) then
            call ra_simple_var_request(options)
        endif

        if (options%physics%radiation == kRA_RRTMG .or. options%physics%radiation == kRA_RRTMGP) then
            call ra_rrtmg_var_request(options)
        endif
        
        !! MJ added: the vars requested if we have radiation_downScaling  
        if (options%physics%radiation_downScaling==1) then        
            call options%alloc_vars( [kVARS%slope, kVARS%slope_angle, kVARS%aspect_angle, kVARS%svf, kVARS%hlm, kVARS%shortwave_direct, &
                                      kVARS%shortwave_diffuse, kVARS%shortwave_direct_above]) 
        endif

    end subroutine ra_var_request


    !> ----------------------------------------------
    !! Communicate to the master process requesting the variables requred to be allocated, used for restart files, and advected
    !!
    !! ----------------------------------------------
    subroutine ra_simple_var_request(options)
        implicit none
        type(options_t), intent(inout) :: options


        ! List the variables that are required to be allocated for the simple radiation code
        call options%alloc_vars( &
                     [kVARS%pressure,    kVARS%potential_temperature,   kVARS%exner,        kVARS%cloud_fraction,   &
                      kVARS%shortwave,   kVARS%longwave, kVARS%cosine_zenith_angle])

        ! List the variables that are required to be advected for the simple radiation code
        call options%advect_vars( &
                      [kVARS%potential_temperature] )

        ! List the variables that are required when restarting for the simple radiation code
        call options%restart_vars( &
                       [kVARS%pressure,     kVARS%potential_temperature, kVARS%shortwave,   kVARS%longwave, kVARS%cloud_fraction, kVARS%cosine_zenith_angle] )

    end subroutine ra_simple_var_request


    !> ----------------------------------------------
    !! Communicate to the master process requesting the variables requred to be allocated, used for restart files, and advected
    !!
    !! ----------------------------------------------
    subroutine ra_rrtmg_var_request(options)
        implicit none
        type(options_t), intent(inout) :: options

        ! List the variables that are required to be allocated for the simple radiation code
        call options%alloc_vars( &
                     [kVARS%pressure,     kVARS%pressure_interface,    kVARS%potential_temperature,   kVARS%exner,            &
                      kVARS%shortwave, kVARS%shortwave_direct, kVARS%shortwave_diffuse,   kVARS%longwave,                                                                     &
                      kVARS%out_longwave_rad, &
                      kVARS%land_mask,    kVARS%snow_water_equivalent,                                                        &
                      kVARS%dz_interface, kVARS%skin_temperature,      kVARS%temperature,             kVARS%density,          &
                      kVARS%longwave_cloud_forcing,                    kVARS%land_emissivity,         kVARS%temperature_interface,  &
                      kVARS%cosine_zenith_angle,                       kVARS%shortwave_cloud_forcing, kVARS%tend_swrad,           &
                      kVARS%tend_th_lwrad, kVARS%tend_th_swrad, kVARS%cloud_fraction, kVARS%albedo])


        ! List the variables that are required when restarting for the simple radiation code
        call options%restart_vars( &
                     [kVARS%pressure,     kVARS%pressure_interface,    kVARS%potential_temperature,   kVARS%exner,            &
                      kVARS%water_vapor,  kVARS%shortwave, kVARS%shortwave_direct, kVARS%shortwave_diffuse,    kVARS%longwave,                                                 &          
                      kVARS%out_longwave_rad, &
                      kVARS%snow_water_equivalent,                                                                            &
                      kVARS%dz_interface, kVARS%skin_temperature,      kVARS%temperature,             kVARS%density,          &
                      kVARS%longwave_cloud_forcing,                    kVARS%land_emissivity, kVARS%temperature_interface,    &
                      kVARS%cosine_zenith_angle,                       kVARS%shortwave_cloud_forcing, kVars%tend_swrad] )

    end subroutine ra_rrtmg_var_request


    subroutine rad(domain, options, dt, halo, subset)
        implicit none

        type(domain_t), intent(inout) :: domain
        type(options_t),intent(in)    :: options
        real,           intent(in)    :: dt
        integer,        intent(in), optional :: halo, subset


        real, dimension(:,:,:,:), pointer :: tauaer_sw=>null(), ssaaer_sw=>null(), asyaer_sw=>null()
        real, allocatable:: albedo(:,:),gsw(:,:)
        real, allocatable:: t_1d(:), p_1d(:), Dz_1d(:), qv_1d(:), qc_1d(:), qi_1d(:), qs_1d(:), cf_1d(:)
        real, allocatable :: qi(:,:,:), qc(:,:,:), qs(:,:,:), cldfra(:,:,:)

        real :: gridkm, ra_dt, declin, hour_frac, air_mass_lay, cld_frc
        real :: relmax, relmin, reimax, reimin
        real(real64) :: date_seconds, julian_day
        integer :: i, k, j, col_indx,sim_month

        logical :: f_qr, f_qc, f_qi, F_QI2, F_QI3, f_qs, f_qg, f_qv, f_qndrop
        integer :: mp_options, F_REC, F_REI, F_RES


        real(wp), dimension(:,:),   allocatable :: p_lay, t_lay, p_lev, t_lev ! t_lev is only needed for LW
        real(wp), dimension(:,:),   allocatable :: q, o3

        integer  :: nbnd, ngpt

        !   Longwave
        type(ty_source_func_lw)               :: lw_sources
        !   Shortwave
        real(wp), dimension(:,:), allocatable :: toa_flux

        !
        ! Clouds
        !
        real(wp), allocatable, dimension(:,:) :: lwp, iwp, swp, rel, rei, res, zwp
        !
        ! Aerosols
        !
        ! logical :: cell_has_aerosols
        ! integer,  dimension(:,:), allocatable :: aero_type
        !                                         ! MERRA2/GOCART aerosol type
        ! real(wp), dimension(:,:), allocatable :: aero_size
        !                                         ! Aerosol size for dust and sea salt
        ! real(wp), dimension(:,:), allocatable :: aero_mass
        !                                         ! Aerosol mass column (kg/m2)
        ! real(wp), dimension(:,:), allocatable :: relhum
        !                                         ! Relative humidity (fraction)
        ! logical, dimension(:,:), allocatable  :: aero_mask
        !                                         ! Aerosol mask

        real(wp), dimension(:),     allocatable :: t_sfc, mu0
        real(wp), dimension(:,:),   allocatable :: emis_sfc, sfc_alb_dir, sfc_alb_dif ! First dimension is band

        character(len=8) :: char_input
        integer :: ncol, nlay, block, jb_s, jb_e
        real(real64) :: sun_declin_deg, eq_of_time_minutes
        !
        ! Timing variables
        !
        type(ty_optical_props_1scl)   :: atmos_lw, clouds_lw, aerosols_lw, snow_lw
        type(ty_optical_props_2str)   :: atmos_sw, clouds_sw, aerosols_sw, snow_sw
        type(ty_fluxes_broadband)    :: fluxes_lw, fluxes_sw

        REAL (wp), dimension(:,:), target, allocatable ::      &
            flx_up, & !< upward SW flux profile, all sky
            flx_dn, & !< downward SW flux profile, all sky
            flx_dnsw_dir !< downward SW flux profile, all sky


        !! MJ added
        real :: trans_atm, trans_atm_dir, max_dir_1, max_dir_2, max_dir, elev_th, ratio_dif, tzone
        integer :: zdx, zdx_max
        
        if (options%physics%radiation == 0) return
        
        !We only need to calculate these variables if we are using terrain shading, otherwise only call on each radiation update
        if (options%physics%radiation_downScaling == 1 .or. &
            ((domain%sim_time%seconds() - last_model_time(domain%nest_indx)) >= update_interval)) then
            tzone = options%rad%tzone
            date_seconds = domain%sim_time%seconds()
            sim_month = domain%sim_time%month
            julian_day = domain%sim_time%date_to_jd()
            hour_frac = domain%sim_time%TOD_hours()

            call calc_solar_date(date_seconds, julian_day, sun_declin_deg, eq_of_time_minutes)

            associate(lat => domain%vars_2d(domain%var_indx(kVARS%latitude)%v)%data_2d, &
                      lon => domain%vars_2d(domain%var_indx(kVARS%longitude)%v)%data_2d, &
                      cosine_zenith_angle => domain%vars_2d(domain%var_indx(kVARS%cosine_zenith_angle)%v)%data_2d)

            !$acc parallel loop gang async(1)
            do j = jms,jme
               !! MJ used corr version, as other does not work in Erupe
                call calc_solar_elevation(solar_elevation=solar_elevation_store(:,j), hour_frac=hour_frac, sun_declin_deg=sun_declin_deg, eq_of_time_minutes=eq_of_time_minutes, tzone=tzone, &
                    lon=lon, lat=lat, j=j, &
                    ims=ims,ime=ime,jms=jms,jme=jme,its=its,ite=ite, solar_azimuth=solar_azimuth_store(:,j))
            enddo

            if (options%physics%radiation_downScaling == 1) then
                associate(slope => domain%vars_2d(domain%var_indx(kVARS%slope_angle)%v)%data_2d, &
                          aspect => domain%vars_2d(domain%var_indx(kVARS%aspect_angle)%v)%data_2d)
                !$acc parallel loop gang vector collapse(2) async(1)
                do j = jms,jme
                    do i = ims,ime
                        cosine_zenith_angle(i,j)=sin(solar_elevation_store(i,j))

                        cos_project_angle(i,j)= cos(slope(i,j))*sin(solar_elevation_store(i,j)) + &
                                                sin(slope(i,j))*cos(solar_elevation_store(i,j))   &
                                                *cos(solar_azimuth_store(i,j)-aspect(i,j))
                    enddo
                enddo
                end associate
            else
                !$acc parallel loop gang vector collapse(2) async(1)
                do j = jms,jme
                    do i = ims,ime
                        cosine_zenith_angle(i,j)=sin(solar_elevation_store(i,j))
                    enddo
                enddo
            endif
            end associate
        endif

        !If we are not over the update interval, don't run any of this, since it contains allocations, etc...
        if ((domain%sim_time%seconds() - last_model_time(domain%nest_indx)) >= update_interval) then

            associate(albedo_dom => domain%vars_3d(domain%var_indx(kVARS%albedo)%v)%data_3d, &
                      shortwave => domain%vars_2d(domain%var_indx(kVARS%shortwave)%v)%data_2d, &
                      shortwave_diffuse => domain%vars_2d(domain%var_indx(kVARS%shortwave_diffuse)%v)%data_2d, &
                      shortwave_direct => domain%vars_2d(domain%var_indx(kVARS%shortwave_direct)%v)%data_2d, &
                      longwave => domain%vars_2d(domain%var_indx(kVARS%longwave)%v)%data_2d, &
                      out_longwave_rad => domain%vars_2d(domain%var_indx(kVARS%out_longwave_rad)%v)%data_2d, &
                      cloud_fraction => domain%vars_2d(domain%var_indx(kVARS%cloud_fraction)%v)%data_2d, &
                      land_mask => domain%vars_2d(domain%var_indx(kVARS%land_mask)%v)%data_2di, &
                      dz_interface => domain%vars_3d(domain%var_indx(kVARS%dz_interface)%v)%data_3d, &
                      temperature => domain%vars_3d(domain%var_indx(kVARS%temperature)%v)%data_3d, &
                      pressure => domain%vars_3d(domain%var_indx(kVARS%pressure)%v)%data_3d, &
                      potential_temperature => domain%vars_3d(domain%var_indx(kVARS%potential_temperature)%v)%data_3d, &
                      temperature_i => domain%vars_3d(domain%var_indx(kVARS%temperature_interface)%v)%data_3d, &
                      pressure_i => domain%vars_3d(domain%var_indx(kVARS%pressure_interface)%v)%data_3d, &
                      exner => domain%vars_3d(domain%var_indx(kVARS%exner)%v)%data_3d, &
                      water_vapor => domain%vars_3d(domain%var_indx(kVARS%water_vapor)%v)%data_3d, &
                      density => domain%vars_3d(domain%var_indx(kVARS%density)%v)%data_3d, &
                      re_cloud => domain%vars_3d(domain%var_indx(kVARS%re_cloud)%v)%data_3d, &
                      re_ice => domain%vars_3d(domain%var_indx(kVARS%re_ice)%v)%data_3d, &
                      re_snow => domain%vars_3d(domain%var_indx(kVARS%re_snow)%v)%data_3d, &
                      land_emissivity => domain%vars_2d(domain%var_indx(kVARS%land_emissivity)%v)%data_2d, &
                      cosine_zenith_angle => domain%vars_2d(domain%var_indx(kVARS%cosine_zenith_angle)%v)%data_2d, &
                      qv_dom => domain%vars_3d(domain%var_indx(kVARS%water_vapor)%v)%data_3d, &
                      qc_dom => domain%vars_3d(domain%var_indx(kVARS%cloud_water_mass)%v)%data_3d, &
                      qi_dom => domain%vars_3d(domain%var_indx(kVARS%ice_mass)%v)%data_3d, &
                      i2_dom => domain%vars_3d(domain%var_indx(kVARS%ice2_mass)%v)%data_3d, &
                      i3_dom => domain%vars_3d(domain%var_indx(kVARS%ice3_mass)%v)%data_3d, &
                      qs_dom => domain%vars_3d(domain%var_indx(kVARS%snow_mass)%v)%data_3d)
            ra_dt = domain%sim_time%seconds() - last_model_time(domain%nest_indx)
            last_model_time(domain%nest_indx) = domain%sim_time%seconds()

            F_QI=.false.
            F_QI2 = .false.
            F_QI3 = .false.
            F_QC=.false.
            F_QR=.false.
            F_QS=.false.
            F_QG=.false.
            f_qndrop=.false.
            F_QV=.false.

            F_REC=0
            F_REI=0
            F_RES=0

            F_QI=(domain%var_indx(kVARS%ice_mass)%v > 0)
            F_QC=(domain%var_indx(kVARS%cloud_water_mass)%v > 0)
            F_QR=(domain%var_indx(kVARS%rain_mass)%v > 0)
            F_QS=(domain%var_indx(kVARS%snow_mass)%v > 0)
            F_QV=(domain%var_indx(kVARS%water_vapor)%v > 0)
            F_QG=(domain%var_indx(kVARS%graupel_mass)%v > 0)
            F_QNDROP=(domain%var_indx(kVARS%cloud_number)%v > 0)
            F_QI2=(domain%var_indx(kVARS%ice2_mass)%v > 0)
            F_QI3=(domain%var_indx(kVARS%ice3_mass)%v > 0)

            if (domain%var_indx(kVARS%re_cloud)%v > 0) F_REC = 1
            if (domain%var_indx(kVARS%re_ice)%v > 0) F_REI = 1
            if (domain%var_indx(kVARS%re_snow)%v > 0) F_RES = 1

            allocate(t_1d(kms:kme))
            allocate(p_1d(kms:kme))
            allocate(Dz_1d(kms:kme))
            allocate(qv_1d(kms:kme))
            allocate(qc_1d(kms:kme))
            allocate(qi_1d(kms:kme))
            allocate(qs_1d(kms:kme))
            allocate(cf_1d(kms:kme))

            allocate(qi(ims:ime,kms:kme,jms:jme))
            allocate(qc(ims:ime,kms:kme,jms:jme))
            allocate(qs(ims:ime,kms:kme,jms:jme))

            allocate(cldfra(ims:ime,kms:kme,jms:jme))

            allocate(albedo(ims:ime,jms:jme))
            allocate(gsw(ims:ime,jms:jme))

            !$acc data create(t_1d, p_1d, Dz_1d, qv_1d, qc_1d, qi_1d, qs_1d, cf_1d, qi, qc, qs, cldfra, albedo, gsw)

            !$acc kernels
            qi = 0
            qc = 0
            qs = 0
            cldfra=0
            !$acc end kernels

            if (F_QI) then
                !$acc parallel loop gang vector collapse(3) present(qi_dom, qi)
                do j = jms,jme
                    do k = kms,kme
                        do i = ims,ime
                        qi(i,k,j) = qi_dom(i,k,j)
                        enddo
                    enddo
                enddo
            endif
            if (F_QC) then
                !$acc parallel loop gang vector collapse(3) present(qc_dom, qc)
                do j = jms,jme
                    do k = kms,kme
                        do i = ims,ime
                        qc(i,k,j) = qc_dom(i,k,j)
                        enddo
                    enddo
                enddo
            endif
            if (F_QS) then
                !$acc parallel loop gang vector collapse(3) present(qs_dom, qs)
                do j = jms,jme
                    do k = kms,kme
                        do i = ims,ime
                        qs(i,k,j) = qs_dom(i,k,j)
                        enddo
                    enddo
                enddo
            endif
            if (F_QI2) then
                !$acc parallel loop gang vector collapse(3) present(i2_dom, qi)
                do j = jms,jme
                    do k = kms,kme
                        do i = ims,ime
                        qi(i,k,j) = qi(i,k,j) + i2_dom(i,k,j)
                        enddo
                    enddo
                enddo
            endif
            if (F_QI3) then
                !$acc parallel loop gang vector collapse(3) present(i3_dom, qi)
                do j = jms,jme
                    do k = kms,kme
                        do i = ims,ime
                        qi(i,k,j) = qi(i,k,j) + i3_dom(i,k,j)
                        enddo
                    enddo
                enddo
            endif


            mp_options=0

            !Calculate solar constant
            call radconst(domain%sim_time%day_of_year(), declin, solar_constant)
            !$acc wait(1)
            if (options%physics%radiation==kRA_SIMPLE) then
                call ra_simple(theta = domain%vars_3d(domain%var_indx(kVARS%potential_temperature)%v)%data_3d,         &
                               pii= domain%vars_3d(domain%var_indx(kVARS%exner)%v)%data_3d,                            &
                               qv = domain%vars_3d(domain%var_indx(kVARS%water_vapor)%v)%data_3d,                      &
                               qc = domain%vars_3d(domain%var_indx(kVARS%cloud_water_mass)%v)%data_3d,                 &
                               qs = domain%vars_3d(domain%var_indx(kVARS%snow_mass)%v)%data_3d + domain%vars_3d(domain%var_indx(kVARS%ice_mass)%v)%data_3d + domain%vars_3d(domain%var_indx(kVARS%graupel_mass)%v)%data_3d,                                    &
                               qr = domain%vars_3d(domain%var_indx(kVARS%rain_mass)%v)%data_3d,                        &
                               p =  domain%vars_3d(domain%var_indx(kVARS%pressure)%v)%data_3d,                         &
                               swdown =  domain%vars_2d(domain%var_indx(kVARS%shortwave)%v)%data_2d,                   &
                               lwdown =  domain%vars_2d(domain%var_indx(kVARS%longwave)%v)%data_2d,                    &
                               cloud_cover =  domain%vars_2d(domain%var_indx(kVARS%cloud_fraction)%v)%data_2d,         &
                               lat = domain%vars_2d(domain%var_indx(kVARS%latitude)%v)%data_2d,                        &
                               lon = domain%vars_2d(domain%var_indx(kVARS%longitude)%v)%data_2d,                       &
                               cosine_zenith_angle = domain%vars_2d(domain%var_indx(kVARS%cosine_zenith_angle)%v)%data_2d,      &
                               date = domain%sim_time,                             &
                               options = options,                                    &
                               dt = ra_dt,                                           &
                               ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, &
                               its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte, F_runlw=.True.)
            else if (options%physics%radiation==kRA_RRTMG .or. options%physics%radiation==kRA_RRTMGP) then

                if (options%lsm%monthly_albedo) then
                    !$acc parallel loop gang vector collapse(2) present(albedo_dom, albedo)
                    do j = jms,jme
                        do i = ims,ime
                            ALBEDO(i,j) = albedo_dom(i, sim_month, j)
                        enddo
                    enddo
                else
                    !$acc parallel loop gang vector collapse(2) present(albedo_dom, albedo)
                    do j = jms,jme
                        do i = ims,ime
                            ALBEDO(i,j) = albedo_dom(i, 1, j)
                        enddo
                    enddo
                endif

                ! domain%tend%th_swrad = 0
                ! domain%vars_2d(domain%var_indx(kVARS%shortwave)%v)%data_2d = 0
                ! Calculate cloud fraction
                If (options%rad%icloud == 3) THEN
                    IF ( F_QC .AND. F_QI ) THEN
                        gridkm = domain%dx/1000
                        !$acc parallel loop gang vector collapse(2) present(cloud_fraction)
                        do j = jts,jte
                            do i = its,ite
                                cloud_fraction(i,j) = 0
                            enddo
                        enddo
                        !$acc parallel loop gang(static: 1) vector collapse(2) present(qv_dom, pressure, temperature, &
                        !$acc                    dz_interface, land_mask, cloud_fraction, qi, qc, qs, cldfra, p_1d, t_1d, qv_1d, qc_1d, qi_1d, qs_1d, Dz_1d, cf_1d) async(1)
                        DO j = jts,jte
                            DO i = its,ite
                                !$acc loop
                                DO k = kms,kme
                                    cf_1d(k) = cldfra(i,k,j)
                                    qv_1d(k) = qv_dom(i,k,j)
                                    qc_1d(k) = qc(i,k,j)
                                    qi_1d(k) = qi(i,k,j)
                                    qs_1d(k) = qs(i,k,j)
                                    p_1d(k) = pressure(i,k,j)
                                    t_1d(k) = temperature(i,k,j)
                                    Dz_1d(k) = dz_interface(i,k,j)
                                ENDDO
                                CALL cal_cldfra3(cf_1d, qv_1d, &
                                              qc_1d, qi_1d, qs_1d, &
                                              Dz_1d, &
                                              p_1d, &
                                              t_1d, &
                                              real(land_mask(i,j)), &
                                              gridkm, 1.5, kms, kme,        &
                                              modify_qvapor=.false., use_multilayer=.False.)
                                !$acc loop
                                DO k = kts,kte
                                    cldfra(i,k,j) = cf_1d(k)
                                    cloud_fraction(i,j) = max(cloud_fraction(i,j), cldfra(i,k,j))
                                ENDDO
                            ENDDO
                        ENDDO
                    END IF
                END IF

                if (options%physics%radiation==kRA_RRTMG) then
                    call RRTMG_SWRAD(rthratensw=domain%tend%th_swrad,         &
    !                swupt, swuptc, swuptcln, swdnt, swdntc, swdntcln, &
    !                swupb, swupbc, swupbcln, swdnb, swdnbc, swdnbcln, &
    !                      swupflx, swupflxc, swdnflx, swdnflxc,      &
                        swdnb = domain%vars_2d(domain%var_indx(kVARS%shortwave)%v)%data_2d,                     &
                        swcf = domain%vars_2d(domain%var_indx(kVARS%shortwave_cloud_forcing)%v)%data_2d,        &
                        gsw = gsw,                                            &
                        xtime = 0., gmt = 0.,                                 &  ! not used
                        xlat = domain%vars_2d(domain%var_indx(kVARS%latitude)%v)%data_2d,                       &  ! not used
                        xlong = domain%vars_2d(domain%var_indx(kVARS%longitude)%v)%data_2d,                     &  ! not used
                        radt = 0., degrad = 0., declin = 0.,                  &  ! not used
                        coszr = domain%vars_2d(domain%var_indx(kVARS%cosine_zenith_angle)%v)%data_2d,           &
                        julday = 0,                                           &  ! not used
                        solcon = solar_constant,                              &
                        albedo = albedo,                                      &
                        t3d = domain%vars_3d(domain%var_indx(kVARS%temperature)%v)%data_3d,                     &
                        t8w = domain%vars_3d(domain%var_indx(kVARS%temperature_interface)%v)%data_3d,           &
                        tsk = domain%vars_2d(domain%var_indx(kVARS%skin_temperature)%v)%data_2d,                &
                        p3d = domain%vars_3d(domain%var_indx(kVARS%pressure)%v)%data_3d,                        &
                        p8w = domain%vars_3d(domain%var_indx(kVARS%pressure_interface)%v)%data_3d,              &
                        pi3d = domain%vars_3d(domain%var_indx(kVARS%exner)%v)%data_3d,                          &
                        rho3d = domain%vars_3d(domain%var_indx(kVARS%density)%v)%data_3d,                       &
                        dz8w = domain%vars_3d(domain%var_indx(kVARS%dz_interface)%v)%data_3d,                   &
                        cldfra3d=cldfra,                                      &
                        !, lradius, iradius,                                  &
                        is_cammgmp_used = .False.,                            &
                        r = R_d,                                               &
                        g = gravity,                                          &
                        re_cloud = domain%vars_3d(domain%var_indx(kVARS%re_cloud)%v)%data_3d,                   &
                        re_ice   = domain%vars_3d(domain%var_indx(kVARS%re_ice)%v)%data_3d,                     &
                        re_snow  = domain%vars_3d(domain%var_indx(kVARS%re_snow)%v)%data_3d,                    &
                        has_reqc=F_REC,                                           & ! use with icloud > 0
                        has_reqi=F_REI,                                           & ! use with icloud > 0
                        has_reqs=F_RES,                                           & ! use with icloud > 0 ! G. Thompson
                        icloud = options%rad%icloud,                  & ! set to nonzero if effective radius is available from microphysics
                        warm_rain = .False.,                                  & ! when a dding WSM3scheme, add option for .True.
                        cldovrlp=options%rad%cldovrlp,                & ! J. Henderson AER: cldovrlp namelist value
                        !f_ice_phy, f_rain_phy,                               &
                        xland=real(domain%vars_2d(domain%var_indx(kVARS%land_mask)%v)%data_2di),                         &
                        xice=real(domain%vars_2d(domain%var_indx(kVARS%land_mask)%v)%data_2di)*0,                        & ! should add a variable for sea ice fraction
                        snow=domain%vars_2d(domain%var_indx(kVARS%snow_water_equivalent)%v)%data_2d,            &
                        qv3d=domain%vars_3d(domain%var_indx(kVARS%water_vapor)%v)%data_3d,                      &
                        qc3d=qc,                                              &
                        qr3d=domain%vars_3d(domain%var_indx(kVARS%rain_mass)%v)%data_3d,                                              &
                        qi3d=qi,                                              &
                        qs3d=qs,                                              &
                        qg3d=domain%vars_3d(domain%var_indx(kVARS%graupel_mass)%v)%data_3d,                                              &
                        !o3input, o33d,                                       &
                        aer_opt=0,                                            &
                        !aerod,                                               &
                        no_src = 1,                                           &
    !                   alswvisdir, alswvisdif,                               &  !Zhenxin ssib alb comp (06/20/2011)
    !                   alswnirdir, alswnirdif,                               &  !Zhenxin ssib alb comp (06/20/2011)
    !                   swvisdir, swvisdif,                                   &  !Zhenxin ssib swr comp (06/20/2011)
    !                   swnirdir, swnirdif,                                   &  !Zhenxin ssib swi comp (06/20/2011)
                        sf_surface_physics=1,                                 &  !Zhenxin
                        f_qv=f_qv, f_qc=f_qc, f_qr=f_qr,                      &
                        f_qi=f_qi, f_qs=f_qs, f_qg=f_qg,                      &
                        !tauaer300,tauaer400,tauaer600,tauaer999,             & ! czhao
                        !gaer300,gaer400,gaer600,gaer999,                     & ! czhao
                        !waer300,waer400,waer600,waer999,                     & ! czhao
    !                   aer_ra_feedback,                                      &
    !jdfcz              progn,prescribe,                                      &
                        calc_clean_atm_diag=0,                                &
    !                    qndrop3d=domain%vars_3d(domain%var_indx(kVARS%cloud_number)%v)%data_3d,                 &
                        f_qndrop=f_qndrop,                                    & !czhao
                        mp_physics=0,                                         & !wang 2014/12
                        ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde, &
                        ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, &
                        its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte-1, &
                        !swupflx, swupflxc,                                   &
                        !swdnflx, swdnflxc,                                   &
                        tauaer3d_sw=tauaer_sw,                                & ! jararias 2013/11
                        ssaaer3d_sw=ssaaer_sw,                                & ! jararias 2013/11
                        asyaer3d_sw=asyaer_sw,                                &
                        swddir = domain%vars_2d(domain%var_indx(kVARS%shortwave_direct)%v)%data_2d,             &
    !                   swddni = domain%vars_2d(domain%var_indx(kVARS%skin_temperature)%v)%data_2d,             &
                        swddif = domain%vars_2d(domain%var_indx(kVARS%shortwave_diffuse)%v)%data_2d,             & ! jararias 2013/08
    !                   swdownc = domain%vars_2d(domain%var_indx(kVARS%skin_temperature)%v)%data_2d,            &
    !                   swddnic = domain%vars_2d(domain%var_indx(kVARS%skin_temperature)%v)%data_2d,            &
    !                   swddirc = domain%vars_2d(domain%var_indx(kVARS%skin_temperature)%v)%data_2d,            &   ! PAJ
                        xcoszen = domain%vars_2d(domain%var_indx(kVARS%cosine_zenith_angle)%v)%data_2d,         &  ! NEED TO CALCULATE THIS.
                        yr=domain%sim_time%year,                            &
                        julian=domain%sim_time%day_of_year(),               &
                        mp_options=mp_options                               )

                    
                    ! DR added December 2023 -- this is necesarry because HICAR does not use a pressure coordinate, so
                    ! p_top is not fixed throughout the simulation. p_top must be reset for each call of RRTMG_LW, 
                    ! which is what RRTMG_LWINIT does. Pass allowed_to_read=.False.
                    ! to avoid re-reading look up tables (already done in init).
                    CALL RRTMG_LWINIT(minval(domain%vars_3d(domain%var_indx(kVARS%pressure_interface)%v)%data_3d(:,domain%kme+1,:)), &
                                    allowed_to_read=.FALSE.,         &
                                    ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde, &
                                    ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, &
                                    its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte)
                    
                    call RRTMG_LWRAD(rthratenlw=domain%tend%th_lwrad,                 &
    !                           lwupt, lwuptc, lwuptcln, lwdnt, lwdntc, lwdntcln,     &        !if lwupt defined, all MUST be defined
    !                           lwupb, lwupbc, lwupbcln, lwdnb, lwdnbc, lwdnbcln,     &
                                glw = domain%vars_2d(domain%var_indx(kVARS%longwave)%v)%data_2d,                        &
                                olr = domain%vars_2d(domain%var_indx(kVARS%out_longwave_rad)%v)%data_2d,                &
                                lwcf = domain%vars_2d(domain%var_indx(kVARS%longwave_cloud_forcing)%v)%data_2d,         &
                                emiss = domain%vars_2d(domain%var_indx(kVARS%land_emissivity)%v)%data_2d,               &
                                p8w = domain%vars_3d(domain%var_indx(kVARS%pressure_interface)%v)%data_3d,              &
                                p3d = domain%vars_3d(domain%var_indx(kVARS%pressure)%v)%data_3d,                        &
                                pi3d = domain%vars_3d(domain%var_indx(kVARS%exner)%v)%data_3d,                          &
                                dz8w = domain%vars_3d(domain%var_indx(kVARS%dz_interface)%v)%data_3d,                   &
                                tsk = domain%vars_2d(domain%var_indx(kVARS%skin_temperature)%v)%data_2d,                &
                                t3d = domain%vars_3d(domain%var_indx(kVARS%temperature)%v)%data_3d,                     &
                                t8w = domain%vars_3d(domain%var_indx(kVARS%temperature_interface)%v)%data_3d,           &     ! temperature interface
                                rho3d = domain%vars_3d(domain%var_indx(kVARS%density)%v)%data_3d,                       &
                                r = R_d,                                               &
                                g = gravity,                                          &
                                icloud = options%rad%icloud,                  & ! set to nonzero if effective radius is available from microphysics
                                warm_rain = .False.,                                  & ! when a dding WSM3scheme, add option for .True.
                                cldfra3d = cldfra,                                    &
                                cldovrlp=options%rad%cldovrlp,                & ! J. Henderson AER: cldovrlp namelist value
    !                            lradius,iradius,                                     & !goes with CAMMGMP (Morrison Gettelman CAM mp)
                                is_cammgmp_used = .False.,                            & !goes with CAMMGMP (Morrison Gettelman CAM mp)
    !                            f_ice_phy, f_rain_phy,                               & !goes with MP option 5 (Ferrier)
                                xland=real(domain%vars_2d(domain%var_indx(kVARS%land_mask)%v)%data_2di),                         &
                                xice=real(domain%vars_2d(domain%var_indx(kVARS%land_mask)%v)%data_2di)*0,                        & ! should add a variable for sea ice fraction
                                snow=domain%vars_2d(domain%var_indx(kVARS%snow_water_equivalent)%v)%data_2d,            &
                                qv3d=domain%vars_3d(domain%var_indx(kVARS%water_vapor)%v)%data_3d,                      &
                                qc3d=qc,                                              &
                                qr3d=domain%vars_3d(domain%var_indx(kVARS%rain_mass)%v)%data_3d,                                              &
                                qi3d=qi,                                              &
                                qs3d=qs,                                              &
                                qg3d=domain%vars_3d(domain%var_indx(kVARS%graupel_mass)%v)%data_3d,                                              &
    !                           o3input, o33d,                                        &
                                f_qv=f_qv, f_qc=f_qc, f_qr=f_qr,                      &
                                f_qi=f_qi, f_qs=f_qs, f_qg=f_qg,                      &
                                re_cloud = domain%vars_3d(domain%var_indx(kVARS%re_cloud)%v)%data_3d,                   &
                                re_ice   = domain%vars_3d(domain%var_indx(kVARS%re_ice)%v)%data_3d,                     &
                                re_snow  = domain%vars_3d(domain%var_indx(kVARS%re_snow)%v)%data_3d,                    &
                                has_reqc=F_REC,                                       & ! use with icloud > 0
                                has_reqi=F_REI,                                       & ! use with icloud > 0
                                has_reqs=F_RES,                                       & ! use with icloud > 0 ! G. Thompson
    !                           tauaerlw1,tauaerlw2,tauaerlw3,tauaerlw4,              & ! czhao
    !                           tauaerlw5,tauaerlw6,tauaerlw7,tauaerlw8,              & ! czhao
    !                           tauaerlw9,tauaerlw10,tauaerlw11,tauaerlw12,           & ! czhao
    !                           tauaerlw13,tauaerlw14,tauaerlw15,tauaerlw16,          & ! czhao
    !                           aer_ra_feedback,                                      & !czhao
    !                    !jdfcz progn,prescribe,                                      & !czhao
                                calc_clean_atm_diag=0,                                & ! used with wrf_chem !czhao
    !                            qndrop3d=domain%vars_3d(domain%var_indx(kVARS%cloud_number)%v)%data_3d,                 & ! used with icould > 0
                                f_qndrop=f_qndrop,                                    & ! if icloud > 0, use this
                            !ccc added for time varying gases.
                                yr=domain%sim_time%year,                             &
                                julian=domain%sim_time%day_of_year(),                &
                            !ccc
                                mp_physics=0,                                          &
                                ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde,  &
                                ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme,  &
                                its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte-1,  &
    !                           lwupflx, lwupflxc, lwdnflx, lwdnflxc,                  &
                                read_ghg=options%rad%read_ghg                  &
                                )
                                                    
                    !$acc update device(domain%vars_3d, domain%vars_2d, domain%tend)

                else if (options%physics%radiation==kRA_RRTMGP) then

                    !cut RRTMGP up into blocks to reduce memory usage
                    do block = 1, nblocks
                    if (block == 1 .and. nblocks==1) then
                        jb_s = jts
                        jb_e = jte
                    elseif (block == 1) then
                        jb_s = jts
                        jb_e = MIN(jts + block*(jte - jts + 1)/nblocks, jte)
                    elseif (block == nblocks) then
                        jb_s = jb_e + 1
                        jb_e = jte
                    else
                        ngpt = jb_e - jb_s + 1
                        jb_s = jb_e + 1
                        jb_e = MIN(jb_s + ngpt - 1, jte)
                    end if
                    ! ----------------------------------------------------------------------------
                    ! 
                    ! ------------------------------ BEGIN INIT  ---------------------------------
                    !
                    ! ----------------------------------------------------------------------------

                    ncol = (ite - its + 1)*(jb_e - jb_s + 1)
                    nlay = (kte - kts + 1)

                    reimin = MAX(cloud_optics_lw%get_min_radius_ice(), cloud_optics_sw%get_min_radius_ice()) 
                    reimax = MIN(cloud_optics_lw%get_max_radius_ice(), cloud_optics_sw%get_max_radius_ice()) 

                    relmin = MAX(cloud_optics_lw%get_min_radius_liq(), cloud_optics_sw%get_min_radius_liq())
                    relmax = MIN(cloud_optics_lw%get_max_radius_liq(), cloud_optics_sw%get_max_radius_liq())


                    allocate(p_lay(ncol, nlay), t_lay(ncol, nlay), p_lev(ncol, nlay+1), t_lev(ncol, nlay+1), q(ncol, nlay), &
                            rel(ncol, nlay), rei(ncol, nlay), lwp(ncol, nlay), zwp(ncol,nlay), iwp(ncol, nlay), swp(ncol, nlay), res(ncol, nlay))
                    !$acc data create(p_lay, t_lay, p_lev, t_lev, q, rel, rei, lwp, zwp, iwp, swp, res)

                    !$acc parallel loop gang vector collapse(2) private(col_indx) present(density, dz_interface, pressure, temperature, &
                    !$acc            qv_dom, re_cloud, re_ice, re_snow, pressure_i, temperature_i, p_lay, t_lay, q, rel, rei, lwp, zwp, iwp, swp, res, p_lev, t_lev, qi) wait(1)
                    do j = jb_s, jb_e
                        do i = its, ite
                            col_indx = (i-its+1) + (j - jb_s)*(ite - its + 1)
                            !$acc loop
                            do k = kts, kte
                                air_mass_lay = density(i,k,j) * dz_interface(i,k,j)
                                p_lay(col_indx,k) = pressure(i,k,j)
                                t_lay(col_indx,k) = temperature(i,k,j)
                                q(col_indx,k)   = qv_dom(i,k,j)*(amd/amw)
                                ! o3(col_indx,k)   = domain%vars_3d(domain%var_indx(kVARS%ozone)%v)%data_3d(i,k,j)


                                zwp(col_indx,k) = 0._wp
                                swp(col_indx,k) = 1000.0*qs(i,k,j) * air_mass_lay

                                if (cldfra(i,k,j) > 0.0_wp) then
                                    cld_frc = MAX(EPSILON(1.0_wp),cldfra(i,k,j))
                                    
                                    lwp(col_indx,k) = 1000.0*qc(i,k,j) * air_mass_lay / cld_frc
                                    iwp(col_indx,k) = 1000.0*qi(i,k,j) * air_mass_lay / cld_frc
                                    rel(col_indx,k) = max(relmin, min(relmax, re_cloud(i,k,j)))
                                    rei(col_indx,k) = max(reimin, min(reimax, re_ice(i,k,j)))
                                    res(col_indx,k) = max(reimin, min(reimax, re_snow(i,k,j)))
                                else
                                    lwp(col_indx,k) = 0.0_wp
                                    iwp(col_indx,k) = 0.0_wp
                                    rel(col_indx,k) = relmin
                                    rei(col_indx,k) = reimin
                                    res(col_indx,k) = reimin
                                end if
                            enddo

                            !$acc loop
                            do k = kts, kte+1
                                p_lev(col_indx,k) = pressure_i(i,k,j)
                                t_lev(col_indx,k) = temperature_i(i,k,j)
                            enddo
                        enddo
                    enddo

                    !$acc update host(p_lay)

                    call set_atmo_gas_conc(domain, block)
                    CALL stop_on_err(gas_concs%set_vmr('h2o',   q))

                    ! ----------------------------------------------------------------------------
                    !  Boundary conditions depending on whether the k-distribution being supplied
                    !   is LW or SW

                    ! ----------------------------------------------------------------------------
                    ! 
                    ! ------------------------------   END INIT  ---------------------------------
                    !
                    ! ----------------------------------------------------------------------------

                    ! ----------------------------------------------------------------------------
                    ! 
                    ! ----------------------------   PERFORM LW  ---------------------------------
                    !
                    ! ----------------------------------------------------------------------------
                    !
                    ! Problem sizes
                    !
                    nbnd = k_dist_lw%get_nband()
                    ngpt = k_dist_lw%get_ngpt()
                    ! LW calculations neglect scattering; SW calculations use the 2-stream approximation

                    allocate(t_sfc(ncol), emis_sfc(nbnd, ncol))
                    allocate(flx_up(ncol,nlay+1), flx_dn(ncol,nlay+1))
                    call stop_on_err(atmos_lw%alloc_1scl(ncol, nlay, k_dist_lw))
                    call stop_on_err(lw_sources%alloc(ncol, nlay, k_dist_lw))
                    CALL stop_on_err(clouds_lw%alloc_1scl(ncol, nlay, k_dist_lw%get_band_lims_wavenumber()))
                    CALL stop_on_err(snow_lw%alloc_1scl(ncol, nlay, k_dist_lw%get_band_lims_wavenumber()))
                    if (do_aerosols) CALL stop_on_err(aerosols_lw%alloc_1scl(ncol, nlay, &
                                             k_dist_lw%get_band_lims_wavenumber()))




                    !$acc data create(t_sfc, emis_sfc, flx_up, flx_dn) &
                    !$ACC   CREATE(lw_sources, atmos_lw, clouds_lw, snow_lw) &
                    !$ACC   CREATE(lw_sources%lay_source, lw_sources%lev_source) &
                    !$ACC   CREATE(lw_sources%sfc_source, lw_sources%sfc_source_Jac) &
                    !$acc   create(atmos_lw%tau, clouds_lw%tau, snow_lw%tau)
                    ! !$acc data create(aerosols_lw, aerosols_lw%tau)


                    ! lw_sources is threadprivate
                    ! Surface temperature

                    !$acc parallel loop gang vector collapse(2) present(t_sfc, emis_sfc, t_lay, land_emissivity)
                    do j = jb_s, jb_e
                        do i = its, ite
                            col_indx = (i-its+1) + (j - jb_s)*(ite - its + 1)
                            t_sfc(col_indx) = t_lay(col_indx,1)
                            !$acc loop
                            do k = 1, nbnd
                                emis_sfc(k,col_indx) = land_emissivity(i,j)
                            enddo
                        enddo
                    enddo

                    call stop_on_err(k_dist_lw%gas_optics(p_lay, p_lev, &
                                                        t_lay, t_sfc, &
                                                        gas_concs,    &
                                                        atmos_lw,        &
                                                        lw_sources, tlev=t_lev))

                    !
                    ! Cloud optics
                    !
                    call stop_on_err(cloud_optics_lw%cloud_optics(lwp, iwp, rel, rei, clouds_lw))
                    call stop_on_err(clouds_lw%increment(atmos_lw))

                    call stop_on_err(cloud_optics_lw%cloud_optics(zwp, swp, res, res, snow_lw))
                    call stop_on_err(snow_lw%increment(atmos_lw))
                    !
                    ! Aerosol optics
                    !
                    ! if(do_aerosols) then
                    !     call stop_on_err(aerosol_optics_lw%aerosol_optics(aero_type, aero_size,  &
                    !                                                 aero_mass, relhum, aerosols_lw))
                    !     call stop_on_err(aerosols_lw%increment(atmos_lw))
                    !     call aerosols_lw%finalize()
                    ! end if

                    fluxes_lw%flux_up => flx_up
                    fluxes_lw%flux_dn => flx_dn
                    call stop_on_err(rte_lw(atmos_lw,      &
                                            lw_sources, &
                                            emis_sfc,   &
                                            fluxes_lw))

                    !$acc parallel loop gang vector collapse(2) present(out_longwave_rad, longwave, flx_up, flx_dn) copyin(kVARS) 
                    do j = jb_s,jb_e 
                        do i = its,ite
                            out_longwave_rad(i,j) = flx_up((i - its + 1) + (j - jb_s)*(ite-its + 1), 1)
                            longwave(i,j) = flx_dn((i - its + 1) + (j - jb_s)*(ite-its + 1), 1)
                        enddo
                    enddo
                    ! Debug prints
                    !$acc        end data
                    call lw_sources%finalize()
                    call atmos_lw%finalize()
                    call clouds_lw%finalize()
                    call snow_lw%finalize()

                    ! ----------------------------------------------------------------------------
                    ! 
                    ! ----------------------------   PERFORM SW  ---------------------------------
                    !
                    ! ----------------------------------------------------------------------------
                    !
                    ! Problem sizes
                    !
                    nbnd = k_dist_sw%get_nband()
                    ngpt = k_dist_sw%get_ngpt()
                    ! LW calculations neglect scattering; SW calculations use the 2-stream approximation

                    call stop_on_err(atmos_sw%alloc_2str( ncol, nlay, k_dist_sw))
                    CALL stop_on_err(clouds_sw%alloc_2str(ncol, nlay, k_dist_sw%get_band_lims_wavenumber()))
                    CALL stop_on_err(snow_sw%alloc_2str(ncol, nlay, k_dist_sw%get_band_lims_wavenumber()))
                    if (do_aerosols) CALL stop_on_err(aerosols_sw%alloc_2str(ncol, nlay, &
                                             k_dist_sw%get_band_lims_wavenumber()))

                    ! toa_flux is threadprivate
                    allocate(toa_flux(ncol, ngpt))
                    allocate(sfc_alb_dir(nbnd, ncol), sfc_alb_dif(nbnd, ncol), mu0(ncol))
                    allocate(flx_dnsw_dir(ncol,nlay+1))


                    !$acc data create(atmos_sw, toa_flux, sfc_alb_dir, sfc_alb_dif, mu0, flx_up, flx_dn, flx_dnsw_dir, clouds_sw, snow_sw) &
                    !$acc   create(atmos_sw%tau, atmos_sw%ssa, atmos_sw%g) &
                    !$acc   create(clouds_sw%tau, clouds_sw%ssa, clouds_sw%g) &
                    !$acc   create(snow_sw%tau, snow_sw%ssa, snow_sw%g)
                    ! !$acc data create(aerosols_sw, aerosols_sw%tau, aerosols_sw%ssa, aerosols_sw%g)

                    !$acc parallel loop gang vector collapse(2) private(col_indx) present(cosine_zenith_angle, ALBEDO, sfc_alb_dir, sfc_alb_dif, mu0)
                    do j = jb_s, jb_e
                        do i = its, ite
                            col_indx = (i-its+1) + (j - jb_s)*(ite - its + 1)

                            mu0(col_indx) = cosine_zenith_angle(i,j)

                            !$acc loop
                            do k = 1, nbnd
                                sfc_alb_dir(k,col_indx) = ALBEDO(i,j)
                                sfc_alb_dif(k,col_indx) = ALBEDO(i,j)
                            enddo
                        enddo
                    enddo

                    call stop_on_err(k_dist_sw%gas_optics(p_lay, p_lev, &
                                                        t_lay,        &
                                                        gas_concs,    &
                                                        atmos_sw,        &
                                                        toa_flux))

                    !
                    ! Cloud optics
                    !
                    call stop_on_err(cloud_optics_sw%cloud_optics(lwp, iwp, rel, rei, clouds_sw))
                    call stop_on_err(clouds_sw%delta_scale())
                    call stop_on_err(clouds_sw%increment(atmos_sw))

                    call stop_on_err(cloud_optics_sw%cloud_optics(zwp, swp, res, res, snow_sw))
                    call stop_on_err(snow_sw%delta_scale())
                    call stop_on_err(snow_sw%increment(atmos_sw))
                    !
                    ! Aerosol optics
                    !
                    ! if(do_aerosols) then
                    !     call stop_on_err(aerosol_optics_sw%aerosol_optics(aero_type, aero_size,  &
                    !                                                 aero_mass, relhum, aerosols_sw))
                    !     call stop_on_err(aerosols_sw%increment(atmos_sw))
                    !     call aerosols_sw%finalize()
                    ! endif

                    fluxes_sw%flux_up => flx_up
                    fluxes_sw%flux_dn => flx_dn
                    fluxes_sw%flux_dn_dir => flx_dnsw_dir

                    call stop_on_err(rte_sw(atmos_sw, mu0,   toa_flux, &
                                            sfc_alb_dir, sfc_alb_dif, &
                                            fluxes_sw))
                    
                    !$acc parallel loop gang vector collapse(2) present(shortwave, shortwave_direct, shortwave_diffuse, flx_up, flx_dn, flx_dnsw_dir) copyin(kVARS) 
                    do j = jb_s,jb_e
                        do i = its,ite
                            shortwave(i,j) = flx_dn((i - its + 1) + (j - jb_s)*(ite-its + 1), 1)
                            shortwave_direct(i,j) = flx_dnsw_dir((i - its + 1) + (j - jb_s)*(ite-its + 1), 1)
                            shortwave_diffuse(i,j) = shortwave(i,j) - shortwave_direct(i,j)
                        enddo
                    enddo
                    !$acc end data
                    !$acc end data
                    call atmos_sw%finalize()
                    call clouds_sw%finalize()
                    call snow_sw%finalize()

                    deallocate(p_lay, t_lay, p_lev, t_lev, q, rel, rei, lwp, zwp, iwp, swp, res)
                    deallocate(toa_flux, sfc_alb_dir, sfc_alb_dif, mu0, flx_up, flx_dn, flx_dnsw_dir, t_sfc, emis_sfc)

                    enddo
                endif 

                ! If the user has provided sky view fraction, then apply this to the diffuse SW now, 
                ! since svf is time-invariant

                associate(tend_swrad => domain%vars_3d(domain%var_indx(kVARS%tend_swrad)%v)%data_3d, &
                          tend_th_swrad => domain%tend%th_swrad)

                if (domain%var_indx(kVARS%svf)%v > 0) then
                    associate(svf => domain%vars_2d(domain%var_indx(kVARS%svf)%v)%data_2d)
                    !$acc parallel loop gang vector collapse(2) present(shortwave_diffuse, svf)
                    do j = jts,jte
                        do i = its,ite
                            shortwave_diffuse(i,j)=shortwave_diffuse(i,j)*svf(i,j)
                        enddo
                    enddo
                    end associate
                endif

                !$acc parallel loop gang vector collapse(3) present(tend_swrad, tend_th_swrad)
                do j = jts,jte
                    do k = kts,kte
                        do i = its,ite
                            tend_swrad(i,k,j) = tend_th_swrad(i,k,j)
                        enddo
                    enddo
                enddo
                end associate
            endif ! end if rrtmg or rrtmgp
            ! cache shortwave from RRTMG_SWRAD for downscaling.
            ! needed if we are to call the terrain shading routine more frequently than RRTMG_SWRAD
            if (options%physics%radiation_downScaling==1) then
                !$acc kernels present(shortwave, shortwave_cached) copyin(kVARS)
                shortwave_cached = shortwave
                !$acc end kernels
            endif
            !$acc end data
            end associate
        end if
        
        
        !! MJ: note that radiation down scaling works only for simple and rrtmg schemes as they provide the above-topography radiation per horizontal plane
        !! MJ corrected, as calc_solar_elevation has largley understimated the zenith angle in Switzerland
        !! MJ added: this is Tobias Jonas (TJ) scheme based on swr function in metDataWizard/PROCESS_COSMO_DATA_1E2E.m and also https://github.com/Tobias-Jonas-SLF/HPEval
        if (options%physics%radiation_downScaling==1) then            
            !! partitioning the total radiation per horizontal plane into the diffusive and direct ones based on https://www.sciencedirect.com/science/article/pii/S0168192320300058, HPEval
            if (.not.(options%physics%radiation==kRA_RRTMG .or. options%physics%radiation==kRA_RRTMGP)) then
                ratio_dif=0.            
                do j = jts,jte
                    do i = its,ite
                        trans_atm = max(min(domain%vars_2d(domain%var_indx(kVARS%shortwave)%v)%data_2d(i,j)/&
                                ( solar_constant* (sin(solar_elevation_store(i,j)+1.e-4)) ),1.),0.)   ! atmospheric transmissivity
                        if (trans_atm<=0.22) then
                            ratio_dif=1.-0.09*trans_atm  
                        elseif (0.22<trans_atm .and. trans_atm<=0.8) then
                            ratio_dif=0.95-0.16*trans_atm+4.39*trans_atm**2.-16.64*trans_atm**3.+12.34*trans_atm**4.   
                        elseif (trans_atm>0.8) then
                            ratio_dif=0.165
                        endif
                        domain%vars_2d(domain%var_indx(kVARS%shortwave_diffuse)%v)%data_2d(i,j)= &
                                ratio_dif*shortwave_cached(i,j)*domain%vars_2d(domain%var_indx(kVARS%svf)%v)%data_2d(i,j)
                    enddo
                enddo                
            endif
            !!
            zdx_max = ubound(domain%vars_3d(domain%var_indx(kVARS%hlm)%v)%data_3d,1)

            associate(shortwave => domain%vars_2d(domain%var_indx(kVARS%shortwave)%v)%data_2d, &
                      shortwave_direct => domain%vars_2d(domain%var_indx(kVARS%shortwave_direct)%v)%data_2d, &
                      shortwave_diffuse => domain%vars_2d(domain%var_indx(kVARS%shortwave_diffuse)%v)%data_2d, &
                      shortwave_direct_above => domain%vars_2d(domain%var_indx(kVARS%shortwave_direct_above)%v)%data_2d, &
                      hlm => domain%vars_3d(domain%var_indx(kVARS%hlm)%v)%data_3d)
            
            !$acc parallel loop gang vector collapse(2) present(shortwave_cached, shortwave_diffuse, shortwave_direct, shortwave, shortwave_direct_above, solar_elevation_store, solar_azimuth_store, cos_project_angle) copyin(kVARS) wait(1)
            do j = jts,jte
                do i = its,ite
                    shortwave_direct(i,j) = max( shortwave_cached(i,j) - shortwave_diffuse(i,j),0.0)

                    ! determin maximum allowed direct swr
                    trans_atm_dir = max(min(shortwave_direct(i,j)/&
                                    (solar_constant*sin(solar_elevation_store(i,j)+1.e-4)),1.),0.)  ! atmospheric transmissivity for direct sw radiation
                    max_dir_1     = solar_constant*exp(log(1.-0.165)/max(sin(solar_elevation_store(i,j)),1.e-4))            
                    max_dir_2     = solar_constant*trans_atm_dir                          
                    max_dir       = min(max_dir_1,max_dir_2)                     ! applying both above criteria 1 and 2                    
                    
                    !!
                    zdx=floor(solar_azimuth_store(i,j)*(180./piconst)/4.0) !! MJ added= we have 90 by 4 deg for hlm ...zidx is the right index based on solar azimuthal angle

                    zdx = max(min(zdx,zdx_max),1)
                    elev_th=(90.-hlm(i,zdx,j))*DEGRAD !! MJ added: it is the solar elevation threshold above which we see the sun from the pixel  
                    if (solar_elevation_store(i,j)>=elev_th) then
                        shortwave_direct_above(i,j)=min(shortwave_direct(i,j),max_dir)
                        shortwave_direct(i,j) = min(shortwave_direct(i,j)/            &
                                                               max(sin(solar_elevation_store(i,j)),0.01),max_dir) * &
                                                               max(cos_project_angle(i,j),0.)
                    else
                        shortwave_direct_above(i,j)=0.
                        shortwave_direct(i,j)=0.
                    endif
                    shortwave(i,j) = shortwave_diffuse(i,j) + shortwave_direct(i,j)
                enddo
            enddo  
            end associate
        endif
    end subroutine rad
    
    subroutine rad_apply_dtheta(domain, options, dt)
        implicit none

        type(domain_t), intent(inout) :: domain
        type(options_t),intent(in)    :: options
        real,           intent(in)    :: dt

        integer :: i, j, k
        !If using the RRTMG scheme, then apply tendencies here. 
        !This is done to allow for the halo exchange to be done before the radiation tendencies are applied
        !This is now outside of interval loop, so this will be called every phys timestep
        if (options%physics%radiation==kRA_RRTMG .or. options%physics%radiation==kRA_RRTMGP) then
            associate(pot_temp => domain%vars_3d(domain%var_indx(kVARS%potential_temperature)%v)%data_3d,  &
                      tend_th_swrad => domain%tend%th_swrad, &
                      tend_th_lwrad => domain%tend%th_lwrad, &
                      tend_swrad    => domain%vars_3d(domain%var_indx(kVARS%tend_swrad)%v)%data_3d)

            !$acc parallel loop gang vector collapse(3) present(pot_temp, tend_th_lwrad, tend_th_swrad, tend_swrad) async(1)
            do j = jts,jte
                do k = kts,kte
                    do i = its,ite
                        pot_temp(i,k,j) = pot_temp(i,k,j)+tend_th_lwrad(i,k,j)*dt+tend_th_swrad(i,k,j)*dt
                        tend_swrad(i,k,j) = tend_th_swrad(i,k,j)
                    enddo
                enddo
            enddo

            end associate
            !$acc wait(1)
        endif

    end subroutine rad_apply_dtheta

    subroutine set_atmo_gas_conc(domain, block_num)
        implicit none
§
        type(domain_t), intent(in) :: domain
        integer, optional, intent(in) :: block_num

        integer :: i, j, k, col_indx, b_num, jb_s, jb_e

        b_num = 1
        if (present(block_num)) b_num = block_num

        call gas_concs%reset()

        call stop_on_err(gas_concs%init(gas_names))

        DO i = 1, ngas
            CALL stop_on_err(gas_concs%set_vmr(gas_names(i), 0._wp))
        ENDDO

        ! Set perscribed atmospheric composition. Water vapor will be set dynamically at each call. 
        ! This can be changed in the future to read from a file or from different climate change scenarios
        call stop_on_err(gas_concs%set_vmr("co2", 410.e-6_wp))
        call stop_on_err(gas_concs%set_vmr("ch4", 1650.e-9_wp))
        call stop_on_err(gas_concs%set_vmr("n2o", 306.e-9_wp))
        call stop_on_err(gas_concs%set_vmr("n2 ",  0.7808_wp))
        call stop_on_err(gas_concs%set_vmr("o2 ",  0.2095_wp))
        call stop_on_err(gas_concs%set_vmr("co ",  0._wp))

    end subroutine set_atmo_gas_conc
    
! DR 2023: Taken from WRF mod_radiation_driver
!---------------------------------------------------------------------
!BOP
! !IROUTINE: radconst - compute radiation terms
! !INTERFAC:
   SUBROUTINE radconst(JULIAN, DECLIN, SOLCON)
!---------------------------------------------------------------------
   IMPLICIT NONE
!---------------------------------------------------------------------

! !ARGUMENTS:
   REAL, INTENT(IN   )      ::       JULIAN
   REAL, INTENT(OUT  )      ::       DECLIN,SOLCON
   REAL                     ::       OBECL,SINOB,SXLONG,ARG,  &
                                     DECDEG,DJUL,RJUL,ECCFAC
!
! !DESCRIPTION:
! Compute terms used in radiation physics 
!EOP

! for short wave radiation

   DECLIN=0.
   SOLCON=0.

!-----OBECL : OBLIQUITY = 23.5 DEGREE.
        
   OBECL=23.5*DEGRAD
   SINOB=SIN(OBECL)
        
!-----CALCULATE LONGITUDE OF THE SUN FROM VERNAL EQUINOX:
        
   IF(JULIAN.GE.80.)SXLONG=DPD*(JULIAN-80.)
   IF(JULIAN.LT.80.)SXLONG=DPD*(JULIAN+285.)
   SXLONG=SXLONG*DEGRAD
   ARG=SINOB*SIN(SXLONG)
   DECLIN=ASIN(ARG)
   DECDEG=DECLIN/DEGRAD
!----SOLAR CONSTANT ECCENTRICITY FACTOR (PALTRIDGE AND PLATT 1976)
   DJUL=JULIAN*360./365.
   RJUL=DJUL*DEGRAD
   ECCFAC=1.000110+0.034221*COS(RJUL)+0.001280*SIN(RJUL)+0.000719*  &
          COS(2*RJUL)+0.000077*SIN(2*RJUL)
   SOLCON=1370.*ECCFAC
   
   END SUBROUTINE radconst
   
  SUBROUTINE stop_on_err(msg)
    USE iso_fortran_env, ONLY : error_unit
    CHARACTER(len=*), INTENT(IN) :: msg

    IF(msg /= "") THEN
      WRITE(error_unit, *) msg
      STOP
    END IF
  END SUBROUTINE


end module radiation
