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
    use icar_constants, only : kVARS, kRA_BASIC, kRA_SIMPLE, kRA_RRTMG, STD_OUT_PE, kMP_THOMP_AER, kMAX_NESTS
    use mod_wrf_constants, only : cp, R_d, gravity, DEGRAD, DPD, piconst
    use mod_atm_utilities, only : cal_cldfra3, calc_solar_elevation

    implicit none
    integer :: update_interval
    real*8  :: last_model_time(kMAX_NESTS)
    real    :: solar_constant
    real    :: p_top = 100000.0
    
    !! MJ added to aggregate radiation over output interval
    real, allocatable, dimension(:,:) :: shortwave_cached, cos_project_angle, solar_elevation_store, solar_azimuth_store
    real, allocatable                 :: solar_azimuth(:), solar_elevation(:)
    real*8 :: counter
    real*8  :: Delta_t !! MJ added to detect the time for outputting 
    integer :: ims, ime, jms, jme, kms, kme
    integer :: its, ite, jts, jte, kts, kte
    integer :: ids, ide, jds, jde, kds, kde

    private
    public :: radiation_init, ra_var_request, rad_apply_dtheta, rad
    
contains

    subroutine radiation_init(domain,options, context_chng)
        type(domain_t), intent(inout) :: domain
        type(options_t),intent(in) :: options
        logical, optional, intent(in) :: context_chng

        logical :: context_change

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

        update_interval=options%rad%update_interval_rrtmg ! 30 min, 1800 s   600 ! 10 min (600 s)

        !Saftey bound, in case update_interval is 0, or very small
        if (.not.(context_change)) then
            if (update_interval<=10) then
                last_model_time(domain%nest_indx) = domain%sim_time%seconds()-10
            else
                last_model_time(domain%nest_indx) = domain%sim_time%seconds()-update_interval
            endif
        endif
        
        if (options%physics%radiation_downScaling==1) then
            if (allocated(cos_project_angle)) deallocate(cos_project_angle)
            allocate(cos_project_angle(ims:ime,jms:jme)) !! MJ added
            if (allocated(solar_elevation_store)) deallocate(solar_elevation_store)
            allocate(solar_elevation_store(ims:ime,jms:jme)) !! MJ added
            if (allocated(solar_azimuth_store)) deallocate(solar_azimuth_store)
            allocate(solar_azimuth_store(ims:ime,jms:jme)) !! MJ added
            if (allocated(shortwave_cached)) deallocate(shortwave_cached)
            allocate(shortwave_cached(ims:ime,jms:jme))
        endif
        
        if (allocated(solar_elevation)) deallocate(solar_elevation)
        allocate(solar_elevation(ims:ime))
        if (allocated(solar_azimuth)) deallocate(solar_azimuth)
        allocate(solar_azimuth(ims:ime)) !! MJ added

        ! If we are just changing nest contexts, we don't need to reinitialize the radiation modules
        if (context_change) return

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
            if(.not.allocated(domain%tend%th_lwrad)) &
                allocate(domain%tend%th_lwrad(domain%ims:domain%ime,domain%kms:domain%kme,domain%jms:domain%jme))
            if(.not.allocated(domain%tend%th_swrad)) &
                allocate(domain%tend%th_swrad(domain%ims:domain%ime,domain%kms:domain%kme,domain%jms:domain%jme))

            if (options%physics%microphysics .ne. kMP_THOMP_AER) then
               if (STD_OUT_PE .and. .not.context_change)write(*,*) '    NOTE: When running RRTMG, microphysics option 5 works best.'
            endif

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
        endif

    end subroutine radiation_init


    subroutine ra_var_request(options)
        implicit none
        type(options_t), intent(inout) :: options

        if (options%physics%radiation == kRA_SIMPLE) then
            call ra_simple_var_request(options)
        endif

        if (options%physics%radiation == kRA_RRTMG) then
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
                      kVARS%cloud_fraction, kVARS%albedo])


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
        real, allocatable :: qc(:,:,:),qi(:,:,:), qs(:,:,:), qg(:,:,:), qr(:,:,:), cldfra(:,:,:)
        real, allocatable :: re_c(:,:,:),re_i(:,:,:), re_s(:,:,:)

        real, allocatable :: xland(:,:)

        real :: gridkm, ra_dt, declin
        integer :: i, k, j

        logical :: f_qr, f_qc, f_qi, F_QI2, F_QI3, f_qs, f_qg, f_qv, f_qndrop
        integer :: mp_options, F_REC, F_REI, F_RES
        
        
        !! MJ added
        real :: trans_atm, trans_atm_dir, max_dir_1, max_dir_2, max_dir, elev_th, ratio_dif
        integer :: zdx, zdx_max
        
        if (options%physics%radiation == 0) return
        
        !We only need to calculate these variables if we are using terrain shading, otherwise only call on each radiation update
        if (options%physics%radiation_downScaling == 1 .or. &
            ((domain%sim_time%seconds() - last_model_time(domain%nest_indx)) >= update_interval)) then
            do j = jms,jme
               !! MJ used corr version, as other does not work in Erupe
                solar_elevation  = calc_solar_elevation(date=domain%sim_time, tzone=options%rad%tzone, &
                    lon=domain%vars_2d(domain%var_indx(kVARS%longitude)%v)%data_2d, lat=domain%vars_2d(domain%var_indx(kVARS%latitude)%v)%data_2d, j=j, &
                    ims=ims,ime=ime,jms=jms,jme=jme,its=its,ite=ite, solar_azimuth=solar_azimuth)
                domain%vars_2d(domain%var_indx(kVARS%cosine_zenith_angle)%v)%data_2d(its:ite,j)=sin(solar_elevation(its:ite))
                
                !If we are doing terrain shading, we will need these later!
                if (options%physics%radiation_downScaling == 1) then
                    cos_project_angle(its:ite,j)= cos(domain%vars_2d(domain%var_indx(kVARS%slope_angle)%v)%data_2d(its:ite,j))*sin(solar_elevation(its:ite)) + &
                                                  sin(domain%vars_2d(domain%var_indx(kVARS%slope_angle)%v)%data_2d(its:ite,j))*cos(solar_elevation(its:ite))   &
                                                  *cos(solar_azimuth(its:ite)-domain%vars_2d(domain%var_indx(kVARS%aspect_angle)%v)%data_2d(its:ite,j))

                    solar_elevation_store(its:ite,j) = solar_elevation(its:ite)
                    solar_azimuth_store(its:ite,j) = solar_azimuth(its:ite)
                endif
            enddo
        endif
        !If we are not over the update interval, don't run any of this, since it contains allocations, etc...
        if ((domain%sim_time%seconds() - last_model_time(domain%nest_indx)) >= update_interval) then

            ra_dt = domain%sim_time%seconds() - last_model_time(domain%nest_indx)
            last_model_time(domain%nest_indx) = domain%sim_time%seconds()

            allocate(t_1d(kms:kme))
            allocate(p_1d(kms:kme))
            allocate(Dz_1d(kms:kme))
            allocate(qv_1d(kms:kme))
            allocate(qc_1d(kms:kme))
            allocate(qi_1d(kms:kme))
            allocate(qs_1d(kms:kme))
            allocate(cf_1d(kms:kme))

            allocate(qc(ims:ime,kms:kme,jms:jme))
            allocate(qi(ims:ime,kms:kme,jms:jme))
            allocate(qs(ims:ime,kms:kme,jms:jme))
            allocate(qg(ims:ime,kms:kme,jms:jme))
            allocate(qr(ims:ime,kms:kme,jms:jme))

            allocate(re_c(ims:ime,kms:kme,jms:jme))
            allocate(re_i(ims:ime,kms:kme,jms:jme))
            allocate(re_s(ims:ime,kms:kme,jms:jme))

            allocate(cldfra(ims:ime,kms:kme,jms:jme))
            allocate(xland(ims:ime,jms:jme))

            allocate(albedo(ims:ime,jms:jme))
            allocate(gsw(ims:ime,jms:jme))

            ! Note, need to link NoahMP to update albedo

            qc = 0
            qi = 0
            qs = 0
            qg = 0
            qr = 0

            re_c = 0
            re_i = 0
            re_s = 0

            cldfra=0

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


            if (F_QG) qg(:,:,:) = domain%vars_3d(domain%var_indx(kVARS%graupel_mass)%v)%data_3d
            if (F_QC) qc(:,:,:) = domain%vars_3d(domain%var_indx(kVARS%cloud_water_mass)%v)%data_3d
            if (F_QI) qi(:,:,:) = domain%vars_3d(domain%var_indx(kVARS%ice_mass)%v)%data_3d
            if (F_QI2) qi(:,:,:) = qi + domain%vars_3d(domain%var_indx(kVARS%ice2_mass)%v)%data_3d
            if (F_QI3) qi(:,:,:) = qi + domain%vars_3d(domain%var_indx(kVARS%ice3_mass)%v)%data_3d
            if (F_QS) qs(:,:,:) = domain%vars_3d(domain%var_indx(kVARS%snow_mass)%v)%data_3d
            if (F_QR) qr(:,:,:) = domain%vars_3d(domain%var_indx(kVARS%rain_mass)%v)%data_3d

            if (F_REC > 0) re_c(:,:,:) = domain%vars_3d(domain%var_indx(kVARS%re_cloud)%v)%data_3d
            if (F_REI > 0) re_i(:,:,:) = domain%vars_3d(domain%var_indx(kVARS%re_ice)%v)%data_3d
            if (F_RES > 0) re_s(:,:,:) = domain%vars_3d(domain%var_indx(kVARS%re_snow)%v)%data_3d

            mp_options=0

            !Calculate solar constant
            call radconst(domain%sim_time%day_of_year(), declin, solar_constant)
           
            if (options%physics%radiation==kRA_SIMPLE) then
                call ra_simple(theta = domain%vars_3d(domain%var_indx(kVARS%potential_temperature)%v)%data_3d,         &
                               pii= domain%vars_3d(domain%var_indx(kVARS%exner)%v)%data_3d,                            &
                               qv = domain%vars_3d(domain%var_indx(kVARS%water_vapor)%v)%data_3d,                      &
                               qc = qc,                 &
                               qs = qs + qi + qg,                                    &
                               qr = qr,                        &
                               p =  domain%vars_3d(domain%var_indx(kVARS%pressure)%v)%data_3d,                         &
                               swdown =  domain%vars_2d(domain%var_indx(kVARS%shortwave)%v)%data_2d,                   &
                               lwdown =  domain%vars_2d(domain%var_indx(kVARS%longwave)%v)%data_2d,                    &
                               cloud_cover =  domain%vars_2d(domain%var_indx(kVARS%cloud_fraction)%v)%data_2d,         &
                               lat = domain%vars_2d(domain%var_indx(kVARS%latitude)%v)%data_2d,                        &
                               lon = domain%vars_2d(domain%var_indx(kVARS%longitude)%v)%data_2d,                       &
                               date = domain%sim_time,                             &
                               options = options,                                    &
                               dt = ra_dt,                                           &
                               ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, &
                               its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte, F_runlw=.True.)
            endif

            if (options%physics%radiation==kRA_RRTMG) then

                if (options%lsm%monthly_albedo) then
                    ALBEDO = domain%vars_3d(domain%var_indx(kVARS%albedo)%v)%data_3d(:, domain%sim_time%month, :)
                else
                    ALBEDO = domain%vars_3d(domain%var_indx(kVARS%albedo)%v)%data_3d(:, 1, :)
                endif

                domain%tend%th_swrad = 0
                domain%vars_2d(domain%var_indx(kVARS%shortwave)%v)%data_2d = 0
                ! Calculate cloud fraction
                If (options%rad%icloud == 3) THEN
                    IF ( F_QC .AND. F_QI ) THEN
                        gridkm = domain%dx/1000
                        XLAND = domain%vars_2d(domain%var_indx(kVARS%land_mask)%v)%data_2di
                        domain%vars_2d(domain%var_indx(kVARS%cloud_fraction)%v)%data_2d = 0
                        DO j = jts,jte
                            DO i = its,ite
                                DO k = kts,kte
                                    p_1d(k) = domain%vars_3d(domain%var_indx(kVARS%pressure)%v)%data_3d(i,k,j) !p(i,k,j)
                                    t_1d(k) = domain%vars_3d(domain%var_indx(kVARS%temperature)%v)%data_3d(i,k,j)
                                    qv_1d(k) = domain%vars_3d(domain%var_indx(kVARS%water_vapor)%v)%data_3d(i,k,j)
                                    qc_1d(k) = qc(i,k,j)
                                    qi_1d(k) = qi(i,k,j)
                                    qs_1d(k) = qs(i,k,j)
                                    Dz_1d(k) = domain%vars_3d(domain%var_indx(kVARS%dz_interface)%v)%data_3d(i,k,j)
                                    cf_1d(k) = cldfra(i,k,j)
                                ENDDO
                                CALL cal_cldfra3(cf_1d, qv_1d, qc_1d, qi_1d, qs_1d, Dz_1d, &
                 &                              p_1d, t_1d, XLAND(i,j), gridkm,        &
                 &                              .false., 1.5, kms, kme)

                                DO k = kts,kte
                                    ! qc, qi and qs are locally recalculated in cal_cldfra3 base on RH to account for subgrid clouds                                     qc(i,k,j) = qc_1d(k)
                                    qc(i,k,j) = qc_1d(k)
                                    qi(i,k,j) = qi_1d(k)
                                    qs(i,k,j) = qs_1d(k)
                                    cldfra(i,k,j) = cf_1d(k)
                                    domain%vars_2d(domain%var_indx(kVARS%cloud_fraction)%v)%data_2d(i,j) = max(domain%vars_2d(domain%var_indx(kVARS%cloud_fraction)%v)%data_2d(i,j), cf_1d(k))
                                ENDDO
                            ENDDO
                        ENDDO
                    END IF
                END IF

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
                    re_cloud = re_c,                   &
                    re_ice   = re_i,                     &
                    re_snow  = re_s,                    &
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
                    qr3d=qr,                                              &
                    qi3d=qi,                                              &
                    qs3d=qs,                                              &
                    qg3d=qg,                                              &
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

                    ! cache shortwave from RRTMG_SWRAD for downscaling.
                    ! needed if we are to call the terrain shading routine more frequently than RRTMG_SWRAD
                    if (options%physics%radiation_downScaling==1) then
                        shortwave_cached = domain%vars_2d(domain%var_indx(kVARS%shortwave)%v)%data_2d
                    endif
                    
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
                            qr3d=qr,                                              &
                            qi3d=qi,                                              &
                            qs3d=qs,                                              &
                            qg3d=qg,                                              &
!                           o3input, o33d,                                        &
                            f_qv=f_qv, f_qc=f_qc, f_qr=f_qr,                      &
                            f_qi=f_qi, f_qs=f_qs, f_qg=f_qg,                      &
                            re_cloud = re_c,                   &
                            re_ice   = re_i,                     &
                            re_snow  = re_s,                    &
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
                            
                ! If the user has provided sky view fraction, then apply this to the diffuse SW now, 
                ! since svf is time-invariant
                if (domain%var_indx(kVARS%svf)%v > 0) then
                    domain%vars_2d(domain%var_indx(kVARS%shortwave_diffuse)%v)%data_2d=domain%vars_2d(domain%var_indx(kVARS%shortwave_diffuse)%v)%data_2d*domain%vars_2d(domain%var_indx(kVARS%svf)%v)%data_2d
                endif
                
                domain%vars_3d(domain%var_indx(kVARS%tend_swrad)%v)%data_3d = domain%tend%th_swrad
            endif
        endif
        
        
        !! MJ: note that radiation down scaling works only for simple and rrtmg schemes as they provide the above-topography radiation per horizontal plane
        !! MJ corrected, as calc_solar_elevation has largley understimated the zenith angle in Switzerland
        !! MJ added: this is Tobias Jonas (TJ) scheme based on swr function in metDataWizard/PROCESS_COSMO_DATA_1E2E.m and also https://github.com/Tobias-Jonas-SLF/HPEval
        if (options%physics%radiation_downScaling==1) then            
            !! partitioning the total radiation per horizontal plane into the diffusive and direct ones based on https://www.sciencedirect.com/science/article/pii/S0168192320300058, HPEval
            if (.not.(options%physics%radiation==kRA_RRTMG)) then
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
            do j = jts,jte
                do i = its,ite
                    domain%vars_2d(domain%var_indx(kVARS%shortwave_direct)%v)%data_2d(i,j) = max( shortwave_cached(i,j) - &
                                                                    domain%vars_2d(domain%var_indx(kVARS%shortwave_diffuse)%v)%data_2d(i,j),0.0)

                    ! determin maximum allowed direct swr
                    trans_atm_dir = max(min(domain%vars_2d(domain%var_indx(kVARS%shortwave_direct)%v)%data_2d(i,j)/&
                                    (solar_constant*sin(solar_elevation_store(i,j)+1.e-4)),1.),0.)  ! atmospheric transmissivity for direct sw radiation
                    max_dir_1     = solar_constant*exp(log(1.-0.165)/max(sin(solar_elevation_store(i,j)),1.e-4))            
                    max_dir_2     = solar_constant*trans_atm_dir                          
                    max_dir       = min(max_dir_1,max_dir_2)                     ! applying both above criteria 1 and 2                    
                    
                    !!
                    zdx=floor(solar_azimuth_store(i,j)*(180./piconst)/4.0) !! MJ added= we have 90 by 4 deg for hlm ...zidx is the right index based on solar azimuthal angle

                    zdx = max(min(zdx,zdx_max),1)
                    elev_th=(90.-domain%vars_3d(domain%var_indx(kVARS%hlm)%v)%data_3d(i,zdx,j))*DEGRAD !! MJ added: it is the solar elevation threshold above which we see the sun from the pixel  
                    if (solar_elevation_store(i,j)>=elev_th) then
                        domain%vars_2d(domain%var_indx(kVARS%shortwave_direct_above)%v)%data_2d(i,j)=min(domain%vars_2d(domain%var_indx(kVARS%shortwave_direct)%v)%data_2d(i,j),max_dir)
                        domain%vars_2d(domain%var_indx(kVARS%shortwave_direct)%v)%data_2d(i,j) = min(domain%vars_2d(domain%var_indx(kVARS%shortwave_direct)%v)%data_2d(i,j)/            &
                                                               max(sin(solar_elevation_store(i,j)),0.01),max_dir) * &
                                                               max(cos_project_angle(i,j),0.)
                    else
                        domain%vars_2d(domain%var_indx(kVARS%shortwave_direct_above)%v)%data_2d(i,j)=0.
                        domain%vars_2d(domain%var_indx(kVARS%shortwave_direct)%v)%data_2d(i,j)=0.
                    endif
                    domain%vars_2d(domain%var_indx(kVARS%shortwave)%v)%data_2d(i,j) = domain%vars_2d(domain%var_indx(kVARS%shortwave_diffuse)%v)%data_2d(i,j) + &
                                                          domain%vars_2d(domain%var_indx(kVARS%shortwave_direct)%v)%data_2d(i,j)
                enddo
            enddo           
        endif
    end subroutine rad
    
    subroutine rad_apply_dtheta(domain, options, dt)
        implicit none

        type(domain_t), intent(inout) :: domain
        type(options_t),intent(in)    :: options
        real,           intent(in)    :: dt

        !If using the RRTMG scheme, then apply tendencies here. 
        !This is done to allow for the halo exchange to be done before the radiation tendencies are applied
        !This is now outside of interval loop, so this will be called every phys timestep
        if (options%physics%radiation==kRA_RRTMG) then
            domain%vars_3d(domain%var_indx(kVARS%potential_temperature)%v)%data_3d = domain%vars_3d(domain%var_indx(kVARS%potential_temperature)%v)%data_3d+domain%tend%th_lwrad*dt+domain%tend%th_swrad*dt
            ! domain%vars_3d(domain%var_indx(kVARS%temperature)%v)%data_3d = domain%vars_3d(domain%var_indx(kVARS%potential_temperature)%v)%data_3d*domain%vars_3d(domain%var_indx(kVARS%exner)%v)%data_3d
            domain%vars_3d(domain%var_indx(kVARS%tend_swrad)%v)%data_3d = domain%tend%th_swrad
        endif

    end subroutine rad_apply_dtheta

    
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
   
end module radiation
