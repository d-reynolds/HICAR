!>----------------------------------------------------------
!! This module provides a wrapper to call various microphysics models
!! It sets up variables specific to the physics package to be used including
!! history variables not currently stored in the domain level data structure
!!
!! The main entry point to the code is mp(domain,options,dt)
!!
!! <pre>
!! Call tree graph :
!!  mp_init->[ external initialization routines]
!!  mp->[   external microphysics routines]
!!  mp_finish
!!
!! High level routine descriptions / purpose
!!   mp_init            - allocates module data and initializes physics package
!!   mp                 - sets up and calls main physics package
!!   mp_finish          - deallocates module memory, place to do the same for physics
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
module microphysics
    ! use data_structures
    use icar_constants
    use mod_wrf_constants,          only: EP_1, EP_2, cp, cpv, XLS, XLV, XLF, R_d, R_v, gravity, epsilon, cliq, cice, psat, rhowater, rhosnow, rhoair0
    use module_mp_thompson_aer,     only: mp_gt_driver_aer, thompson_aer_init
    use module_mp_thompson,         only: mp_gt_driver, thompson_init
    use MODULE_MP_MORR_TWO_MOMENT,  only: MORR_TWO_MOMENT_INIT, MP_MORR_TWO_MOMENT
    use module_mp_wsm6,             only: wsm6, wsm6init
    use module_mp_wsm3,             only: wsm3, wsm3init
    use module_mp_simple,           only: mp_simple_driver
    use module_mp_jensen_ishmael,   only: mp_jensen_ishmael, jensen_ishmael_init

    use time_object,                only: Time_type
    use options_interface,          only: options_t
    use domain_interface,           only: domain_t
    use wind,                       only: calc_w_real

    implicit none

    ! permit the microphysics to update on a longer time step than the advection
    integer :: update_interval
    real*8 :: last_model_time
    ! temporary variables
    real,allocatable,dimension(:,:) :: SR, last_rain, last_snow, last_graup, refl_10cm


    public :: mp, mp_var_request
contains


    !>----------------------------------------------------------
    !! Initialize microphysical routines
    !!
    !! This routine will call the initialization routines for the specified microphysics packages.
    !! It also initializes any module level variables, e.g. update_interval
    !!
    !! @param   options     ICAR model options to specify required initializations
    !!
    !!----------------------------------------------------------
    subroutine mp_init(domain,options,context_chng)
        implicit none
        type(domain_t), intent(in) :: domain
        type(options_t), intent(in) :: options
        logical, optional, intent(in) :: context_chng

        logical :: context_change
        if (present(context_chng)) then
            context_change = context_chng
        else
            context_change = .false.
        endif

        update_interval = options%mp%update_interval
        last_model_time = -999

        call allocate_module_variables(domain%ims, domain%ime, domain%jms, domain%jme, domain%kms, domain%kme)

        ! If we are just changing the nest context, no need to re-initialize microphysics routines themselves
        if (context_change) return

        if (STD_OUT_PE .and. .not.context_change) write(*,*) ""
        if (STD_OUT_PE .and. .not.context_change) write(*,*) "Initializing Microphysics"
        if (options%physics%microphysics    == kMP_THOMPSON) then
            if (STD_OUT_PE .and. .not.context_change) write(*,*) "    Thompson Microphysics"
            call thompson_init(options%mp)
        elseif (options%physics%microphysics    == kMP_THOMP_AER) then
            if (STD_OUT_PE .and. .not.context_change) write(*,*) "    Thompson Eidhammer Microphysics"
            if (STD_OUT_PE .and. .not.context_change) write(*,*) "    ERROR: not yet implemented (init function needs arguments), exiting..."
            stop
!            call thompson_aer_init()
        elseif (options%physics%microphysics == kMP_SB04) then
            if (STD_OUT_PE .and. .not.context_change) write(*,*) "    Simple Microphysics"
        elseif (options%physics%microphysics==kMP_MORRISON) then
            if (STD_OUT_PE .and. .not.context_change) write(*,*) "    Morrison Microphysics"
            call MORR_TWO_MOMENT_INIT(hail_opt=1)
        elseif (options%physics%microphysics==kMP_ISHMAEL) then
            if (STD_OUT_PE .and. .not.context_change) write(*,*) "    Jensen-Ischmael Microphysics"
            call jensen_ishmael_init()
        elseif (options%physics%microphysics==kMP_WSM6) then
            if (STD_OUT_PE .and. .not.context_change) write(*,*) "    WSM6 Microphysics"
            call wsm6init(rhoair0,rhowater,rhosnow,cliq,cpv)
        elseif (options%physics%microphysics==kMP_WSM3) then
            if (STD_OUT_PE .and. .not.context_change) write(*,*) "    WSM3 Microphysics"
            call wsm3init(rhoair0,rhowater,rhosnow,cliq,cpv, allowed_to_read=.True.)
        endif

    end subroutine mp_init


    subroutine mp_thompson_aer_var_request(options)
        implicit none
        type(options_t), intent(inout) :: options

        ! List the variables that are required to be allocated for the simple microphysics
        call options%alloc_vars( &
                     [kVARS%pressure,    kVARS%potential_temperature,   kVARS%exner,        kVARS%density,      &
                      kVARS%water_vapor, kVARS%cloud_water_mass,             kVARS%rain_mass,  kVARS%rain_number, &
                      kVARS%snow_mass, kVARS%ice_mass,               kVARS%w,            kVARS%ice_number,      &
                      kVARS%snowfall,    kVARS%precipitation,           kVARS%graupel,      kVARS%graupel_mass,     &
                      kVARS%dz, kVARS%re_cloud, kVARS%re_ice, kVARS%re_snow ])

        ! List the variables that are required to be advected for the simple microphysics
        call options%advect_vars( &
                      [kVARS%potential_temperature, kVARS%water_vapor, kVARS%cloud_water_mass,  kVARS%rain_number, &
                       kVARS%snow_mass,           kVARS%ice_mass,   &
                       kVARS%rain_mass,           kVARS%ice_number, kVARS%graupel_mass   ] )

        ! List the variables that are required to be allocated for the simple microphysics
        call options%restart_vars( &
                       [kVARS%pressure,     kVARS%potential_temperature,    kVARS%water_vapor,   &
                        kVARS%cloud_water_mass,  kVARS%rain_mass,              kVARS%snow_mass,   &
                        kVARS%precipitation,kVARS%snowfall,                 kVARS%graupel,       &
                        kVARS%dz,           kVARS%snow_mass,              kVARS%ice_mass,     &
                        kVARS%rain_number, kVARS%rain_mass,                      &
                        kVARS%ice_number,  kVARS%graupel_mass,                   &
                        kVARS%re_cloud, kVARS%re_ice, kVARS%re_snow  ] )


    end subroutine

    subroutine mp_wsm6_var_request(options)
        implicit none
        type(options_t), intent(inout) :: options

        ! List the variables that are required to be allocated for the simple microphysics
        call options%alloc_vars( &
                     [kVARS%pressure,    kVARS%potential_temperature,   kVARS%exner,        kVARS%density,      &
                      kVARS%water_vapor, kVARS%cloud_water_mass,             kVARS%rain_mass,                      &
                      kVARS%snow_mass, kVARS%ice_mass,               kVARS%dz,                               &
                      kVARS%snowfall,    kVARS%precipitation,           kVARS%graupel,   kVARS%graupel_mass])

        ! List the variables that are required to be advected for the simple microphysics
        call options%advect_vars( &
                      [kVARS%potential_temperature, kVARS%water_vapor, kVARS%cloud_water_mass,  &
                       kVARS%snow_mass,           kVARS%ice_mass,   &
                       kVARS%rain_mass,           kVARS%graupel_mass   ] )

        ! List the variables that are required to be allocated for the simple microphysics
        call options%restart_vars( &
                       [kVARS%pressure,     kVARS%potential_temperature,    kVARS%water_vapor,   &
                        kVARS%cloud_water_mass,  kVARS%rain_mass,              kVARS%snow_mass,   &
                        kVARS%precipitation,kVARS%snowfall,                 kVARS%graupel,       &
                        kVARS%dz,           kVARS%snow_mass,              kVARS%ice_mass,     &
                        kVARS%rain_mass,  kVARS%graupel_mass] )


    end subroutine


    subroutine mp_wsm3_var_request(options)
        implicit none
        type(options_t), intent(inout) :: options

        ! List the variables that are required to be allocated for the simple microphysics
        call options%alloc_vars( &
                     [kVARS%pressure,    kVARS%potential_temperature,   kVARS%exner,        kVARS%density,      &
                      kVARS%water_vapor, kVARS%cloud_water_mass,             kVARS%rain_mass,                      &
                      kVARS%dz,          kVARS%snowfall,                kVARS%precipitation ])

        ! List the variables that are required to be advected for the simple microphysics
        call options%advect_vars( &
                      [kVARS%potential_temperature, kVARS%water_vapor, kVARS%cloud_water_mass,  &
                       kVARS%rain_mass ] )

        ! List the variables that are required to be allocated for the simple microphysics
        call options%restart_vars( &
                       [kVARS%pressure,     kVARS%potential_temperature,    kVARS%water_vapor,   &
                        kVARS%cloud_water_mass,  kVARS%rain_mass,                                   &
                        kVARS%precipitation,kVARS%snowfall,                 kVARS%graupel,       &
                        kVARS%dz, kVARS%rain_mass ] )


    end subroutine

    subroutine mp_ishmael_var_request(options)
        implicit none
        type(options_t), intent(inout) :: options

        ! List the variables that are required to be allocated for the simple microphysics
        call options%alloc_vars( &
                     [kVARS%pressure,    kVARS%potential_temperature,   kVARS%exner,        kVARS%density, kVARS%w,     &
                       kVARS%water_vapor, kVARS%cloud_water_mass,  kVARS%rain_number, &
                       kVARS%ice_mass,   kVARS%rain_mass,           kVARS%ice_number, kVARS%ice1_a, &
                       kVARS%ice1_c, kVARS%ice2_mass, kVARS%ice2_number, kVARS%ice2_a, kVARS%ice2_c, &
                       kVARS%ice3_mass, kVARS%ice3_number, kVARS%ice3_a, kVARS%ice3_c, &
                       kVARS%snowfall,    kVARS%precipitation,  kVARS%dz,   kVARS%re_cloud, kVARS%re_ice, &
                       kVARS%ice1_rho, kVARS%ice1_phi, kVARS%ice1_vmi, kVARS%ice2_rho, kVARS%ice2_phi, kVARS%ice2_vmi, &
			kVARS%ice3_rho, kVARS%ice3_phi, kVARS%ice3_vmi  ])

        ! List the variables that are required to be advected for the simple microphysics
        call options%advect_vars( &
                      [kVARS%potential_temperature, kVARS%water_vapor, kVARS%cloud_water_mass,  kVARS%rain_number, &
                       kVARS%ice_mass,   kVARS%rain_mass,           kVARS%ice_number, kVARS%ice1_a, &
                       kVARS%ice1_c, kVARS%ice2_mass, kVARS%ice2_number, kVARS%ice2_a, kVARS%ice2_c, &
                       kVARS%ice3_mass, kVARS%ice3_number, kVARS%ice3_a, kVARS%ice3_c] )

        ! List the variables that are required to be allocated for the simple microphysics
        call options%restart_vars( &
                     [kVARS%pressure,    kVARS%potential_temperature,   kVARS%exner,        kVARS%density, kVARS%w,     &
                       kVARS%water_vapor, kVARS%cloud_water_mass,  kVARS%rain_number, &
                       kVARS%ice_mass,   kVARS%rain_mass,           kVARS%ice_number, kVARS%ice1_a, &
                       kVARS%ice1_c, kVARS%ice2_mass, kVARS%ice2_number, kVARS%ice2_a, kVARS%ice2_c, &
                       kVARS%ice3_mass, kVARS%ice3_number, kVARS%ice3_a, kVARS%ice3_c, &
                       kVARS%snowfall,    kVARS%precipitation,  kVARS%dz,   kVARS%re_cloud, kVARS%re_ice    ])



    end subroutine

    subroutine mp_morr_var_request(options)
        implicit none
        type(options_t), intent(inout) :: options

        ! List the variables that are required to be allocated for the simple microphysics
        call options%alloc_vars( &
                     [kVARS%pressure,    kVARS%potential_temperature,   kVARS%exner,        kVARS%density,      &
                      kVARS%water_vapor, kVARS%cloud_water_mass,             kVARS%rain_mass,  kVARS%rain_number, &
                      kVARS%snow_mass, kVARS%ice_mass,               kVARS%w,            kVARS%ice_number,      &
                      kVARS%snowfall,    kVARS%precipitation,           kVARS%graupel,      kVARS%graupel_mass,     &
                      kVARS%graupel_number, kVARS%snow_number, &
                      kVARS%tend_qr, kVARS%tend_qs, kVARS%tend_qi, kVARS%dz,   &
                      kVARS%re_cloud, kVARS%re_ice, kVARS%re_snow, kVARS%ice1_vmi, kVARS%ice2_vmi    ])

        ! List the variables that are required to be advected for the simple microphysics
        call options%advect_vars( &
                      [kVARS%potential_temperature, kVARS%water_vapor, kVARS%cloud_water_mass,  kVARS%rain_number, &
                       kVARS%snow_mass,           kVARS%ice_mass,   &
                       kVARS%rain_mass,           kVARS%ice_number, kVARS%graupel_mass, &
                       kVARS%graupel_number, kVARS%snow_number] )

        ! List the variables that are required to be allocated for the simple microphysics
        call options%restart_vars( &
                       [kVARS%pressure,     kVARS%potential_temperature,    kVARS%water_vapor,   &
                        kVARS%cloud_water_mass,  kVARS%rain_mass,              kVARS%snow_mass,   &
                        kVARS%precipitation,kVARS%snowfall,                 kVARS%graupel,       &
                        kVARS%dz,           kVARS%snow_mass,              kVARS%ice_mass,     &
                        kVARS%rain_number, kVARS%rain_mass,  &
                        kVARS%ice_number,  kVARS%graupel_mass, &
                        kVARS%snow_number, kVARS%graupel_number, &
                        kVARS%re_cloud, kVARS%re_ice, kVARS%re_snow   ] )



    end subroutine

    subroutine mp_simple_var_request(options)
        implicit none
        type(options_t), intent(inout) :: options


        ! List the variables that are required to be allocated for the simple microphysics
        call options%alloc_vars( &
                     [kVARS%pressure,    kVARS%potential_temperature,   kVARS%exner,        kVARS%density,      &
                      kVARS%water_vapor, kVARS%cloud_water_mass,             kVARS%rain_mass,  kVARS%snow_mass,  &
                      kVARS%precipitation, kVARS%snowfall,              kVARS%dz])

        ! List the variables that are required to be advected for the simple microphysics
        call options%advect_vars( &
                      [kVARS%potential_temperature, kVARS%water_vapor, kVARS%cloud_water_mass,   &
                       kVARS%rain_mass,           kVARS%snow_mass   ] )

        ! List the variables that are required to be allocated for the simple microphysics
        call options%restart_vars( &
                       [kVARS%pressure,     kVARS%potential_temperature,    kVARS%water_vapor,  &
                        kVARS%cloud_water_mass,  kVARS%rain_mass,              kVARS%snow_mass,  &
                        kVARS%precipitation,kVARS%snowfall,                 kVARS%dz] )

    end subroutine mp_simple_var_request


    subroutine mp_var_request(options)
        implicit none
        type(options_t), intent(inout) :: options

        if (options%physics%microphysics    == kMP_THOMPSON) then
            call mp_thompson_aer_var_request(options)
        elseif (options%physics%microphysics    == kMP_THOMP_AER) then
            call mp_thompson_aer_var_request(options)
        elseif (options%physics%microphysics == kMP_SB04) then
            call mp_simple_var_request(options)
        elseif (options%physics%microphysics==kMP_MORRISON) then
            call mp_morr_var_request(options)
        elseif (options%physics%microphysics==kMP_ISHMAEL) then
            call mp_ishmael_var_request(options)
        elseif (options%physics%microphysics==kMP_WSM6) then
            call mp_wsm6_var_request(options)
        elseif (options%physics%microphysics==kMP_WSM3) then
            call mp_wsm3_var_request(options)
        endif
    end subroutine mp_var_request


    subroutine allocate_module_variables(ims,ime,jms,jme,kms,kme)
        implicit none
        integer, intent(in) :: ims,ime,jms,jme,kms,kme

        ! snow rain ratio
        if (allocated(SR)) deallocate(SR)
        allocate(SR(ims:ime,jms:jme))
        SR=0
        ! last snow amount
        if (allocated(last_snow)) deallocate(last_snow)
        allocate(last_snow(ims:ime,jms:jme))
        last_snow=0
        ! last rain amount
        if (allocated(last_rain)) deallocate(last_rain)
        allocate(last_rain(ims:ime,jms:jme))
        last_rain=0
        ! temporary precip amount
        if (allocated(last_graup)) deallocate(last_graup)
        allocate(last_graup(ims:ime,jms:jme))
        last_graup=0

        if (allocated(refl_10cm)) deallocate(refl_10cm)
        allocate(refl_10cm(ims:ime,jms:jme))
        refl_10cm=0
        ! if (.not.allocated(domain%tend%qr)) then
        !     allocate(domain%tend%qr(nx,nz,ny))
        !     domain%tend%qr=0
        ! endif
        ! if (.not.allocated(domain%tend%qs)) then
        !     allocate(domain%tend%qs(nx,nz,ny))
        !     domain%tend%qs=0
        ! endif
        ! if (.not.allocated(domain%tend%qi)) then
        !     allocate(domain%tend%qi(nx,nz,ny))
        !     domain%tend%qi=0
        ! endif


    end subroutine


    subroutine process_subdomain(domain, options, dt,       &
                                its,ite, jts,jte, kts,kte,  &
                                ims,ime, jms,jme, kms,kme,  &
                                ids,ide, jds,jde, kds,kde)
        implicit none
        type(domain_t), intent(inout) :: domain
        type(options_t),intent(in)    :: options
        real,           intent(in)    :: dt
        integer,        intent(in)    :: its,ite, jts,jte, kts,kte
        integer,        intent(in)    :: ims,ime, jms,jme, kms,kme
        integer,        intent(in)    :: ids,ide, jds,jde, kds,kde

        ! run the thompson microphysics
        if (options%physics%microphysics==kMP_THOMPSON) then
            ! call the thompson microphysics
            call mp_gt_driver(qv = domain%vars_3d(domain%var_indx(kVARS%water_vapor)%v)%data_3d,                      &
                              th = domain%vars_3d(domain%var_indx(kVARS%potential_temperature)%v)%data_3d,            &
                              qc = domain%vars_3d(domain%var_indx(kVARS%cloud_water_mass)%v)%data_3d,                 &
                              qi = domain%vars_3d(domain%var_indx(kVARS%ice_mass)%v)%data_3d,                   &
                              ni = domain%vars_3d(domain%var_indx(kVARS%ice_number)%v)%data_3d,                 &
                              qr = domain%vars_3d(domain%var_indx(kVARS%rain_mass)%v)%data_3d,                        &
                              nr = domain%vars_3d(domain%var_indx(kVARS%rain_number)%v)%data_3d,                      &
                              qs = domain%vars_3d(domain%var_indx(kVARS%snow_mass)%v)%data_3d,                        &
                              qg = domain%vars_3d(domain%var_indx(kVARS%graupel_mass)%v)%data_3d,                     &
                              pii= domain%vars_3d(domain%var_indx(kVARS%exner)%v)%data_3d,                            &
                              p =  domain%vars_3d(domain%var_indx(kVARS%pressure)%v)%data_3d,                         &
                              dz = domain%vars_3d(domain%var_indx(kVARS%dz)%v)%data_3d,                          &
                              dt_in = dt,                                           &
                              itimestep = 1,                                        & ! not used in thompson
                              RAINNC = domain%vars_2d(domain%var_indx(kVARS%precipitation)%v)%data_2d,    &
                              RAINNCV = last_rain,                                & ! not used outside thompson (yet)
                              SR = SR,                                              & ! not used outside thompson (yet)
                              SNOWNC = domain%vars_2d(domain%var_indx(kVARS%snowfall)%v)%data_2d,         &
                              GRAUPELNC = domain%vars_2d(domain%var_indx(kVARS%graupel)%v)%data_2d,       &
                              ids = ids, ide = ide,                   & ! domain dims
                              jds = jds, jde = jde,                   &
                              kds = kds, kde = kde,                   &
                              ims = ims, ime = ime,                   & ! memory dims
                              jms = jms, jme = jme,                   &
                              kms = kms, kme = kme,                   &
                              its = its, ite = ite,                   & ! tile dims
                              jts = jts, jte = jte,                   &
                              kts = kts, kte = kte)

        elseif (options%physics%microphysics==kMP_THOMP_AER) then
            ! call the thompson microphysics
            call mp_gt_driver_aer(qv = domain%vars_3d(domain%var_indx(kVARS%water_vapor)%v)%data_3d,                      &
                                  th = domain%vars_3d(domain%var_indx(kVARS%potential_temperature)%v)%data_3d,            &
                                  qc = domain%vars_3d(domain%var_indx(kVARS%cloud_water_mass)%v)%data_3d,                 &
                                  qi = domain%vars_3d(domain%var_indx(kVARS%ice_mass)%v)%data_3d,                   &
                                  ni = domain%vars_3d(domain%var_indx(kVARS%ice_number)%v)%data_3d,                 &
                                  qr = domain%vars_3d(domain%var_indx(kVARS%rain_mass)%v)%data_3d,                        &
                                  nr = domain%vars_3d(domain%var_indx(kVARS%rain_number)%v)%data_3d,                      &
                                  qs = domain%vars_3d(domain%var_indx(kVARS%snow_mass)%v)%data_3d,                        &
                                  qg = domain%vars_3d(domain%var_indx(kVARS%graupel_mass)%v)%data_3d,                     &
                                  pii= domain%vars_3d(domain%var_indx(kVARS%exner)%v)%data_3d,                            &
                                  p =  domain%vars_3d(domain%var_indx(kVARS%pressure)%v)%data_3d,                         &
                                  w =  domain%vars_3d(domain%var_indx(kVARS%w_real)%v)%data_3d,                           &
                                  dz = domain%vars_3d(domain%var_indx(kVARS%dz_interface)%v)%data_3d,                     &
                                  dt_in = dt,                                           &
                                  RAINNC = domain%vars_2d(domain%var_indx(kVARS%precipitation)%v)%data_2d,    &
                                  SNOWNC = domain%vars_2d(domain%var_indx(kVARS%snowfall)%v)%data_2d,         &
                                  GRAUPELNC = domain%vars_2d(domain%var_indx(kVARS%graupel)%v)%data_2d,       &
                                  re_cloud = domain%vars_3d(domain%var_indx(kVARS%re_cloud)%v)%data_3d,                   &
                                  re_ice   = domain%vars_3d(domain%var_indx(kVARS%re_ice)%v)%data_3d,                     &
                                  re_snow  = domain%vars_3d(domain%var_indx(kVARS%re_snow)%v)%data_3d,                    &
                                  has_reqc=1, has_reqi=1, has_reqs=1,                   &
                                  ids = ids, ide = ide,                   & ! domain dims
                                  jds = jds, jde = jde,                   &
                                  kds = kds, kde = kde,                   &
                                  ims = ims, ime = ime,                   & ! memory dims
                                  jms = jms, jme = jme,                   &
                                  kms = kms, kme = kme,                   &
                                  its = its, ite = ite,                   & ! tile dims
                                  jts = jts, jte = jte,                   &
                                  kts = kts, kte = kte)

        elseif (options%physics%microphysics==kMP_SB04) then
            ! call the simple microphysics routine of SB04
            call mp_simple_driver(domain%vars_3d(domain%var_indx(kVARS%pressure)%v)%data_3d,                  &
                                  domain%vars_3d(domain%var_indx(kVARS%potential_temperature)%v)%data_3d,     &
                                  domain%vars_3d(domain%var_indx(kVARS%exner)%v)%data_3d,                     &
                                  domain%vars_3d(domain%var_indx(kVARS%density)%v)%data_3d,                   &
                                  domain%vars_3d(domain%var_indx(kVARS%water_vapor)%v)%data_3d,               &
                                  domain%vars_3d(domain%var_indx(kVARS%cloud_water_mass)%v)%data_3d,          &
                                  domain%vars_3d(domain%var_indx(kVARS%rain_mass)%v)%data_3d,                 &
                                  domain%vars_3d(domain%var_indx(kVARS%snow_mass)%v)%data_3d,                 &
                                  domain%vars_2d(domain%var_indx(kVARS%precipitation)%v)%data_2d, &
                                  domain%vars_2d(domain%var_indx(kVARS%snowfall)%v)%data_2d,      &
                                  dt,                                       &
                                  domain%vars_3d(domain%var_indx(kVARS%dz)%v)%data_3d,                   &
                                  ims = ims, ime = ime,                   & ! memory dims
                                  jms = jms, jme = jme,                   &
                                  kms = kms, kme = kme,                   &
                                  its = its, ite = ite,                   & ! tile dims
                                  jts = jts, jte = jte,                   &
                                  kts = kts, kte = kte)

        elseif (options%physics%microphysics==kMP_MORRISON) then
            call MP_MORR_TWO_MOMENT(ITIMESTEP = 1,                   &
                             TH = domain%vars_3d(domain%var_indx(kVARS%potential_temperature)%v)%data_3d,   &
                             QV = domain%vars_3d(domain%var_indx(kVARS%water_vapor)%v)%data_3d,             &
                             QC = domain%vars_3d(domain%var_indx(kVARS%cloud_water_mass)%v)%data_3d,        &
                             QR = domain%vars_3d(domain%var_indx(kVARS%rain_mass)%v)%data_3d,               &
                             QI = domain%vars_3d(domain%var_indx(kVARS%ice_mass)%v)%data_3d,          &
                             QS = domain%vars_3d(domain%var_indx(kVARS%snow_mass)%v)%data_3d,               &
                             QG = domain%vars_3d(domain%var_indx(kVARS%graupel_mass)%v)%data_3d,            &
                             NI = domain%vars_3d(domain%var_indx(kVARS%ice_number)%v)%data_3d,        &
                             NS = domain%vars_3d(domain%var_indx(kVARS%snow_number)%v)%data_3d,             &
                             NR = domain%vars_3d(domain%var_indx(kVARS%rain_number)%v)%data_3d,             &
                             NG = domain%vars_3d(domain%var_indx(kVARS%graupel_number)%v)%data_3d,          &
                             RHO = domain%vars_3d(domain%var_indx(kVARS%density)%v)%data_3d,                &
                             PII = domain%vars_3d(domain%var_indx(kVARS%exner)%v)%data_3d,                  &
                             P = domain%vars_3d(domain%var_indx(kVARS%pressure)%v)%data_3d,                 &
                             DT_IN = dt, DZ = domain%vars_3d(domain%var_indx(kVARS%dz_interface)%v)%data_3d,     &
                             W = domain%vars_3d(domain%var_indx(kVARS%w_real)%v)%data_3d,                        &
                             RAINNC = domain%vars_2d(domain%var_indx(kVARS%precipitation)%v)%data_2d, &
                             RAINNCV = last_rain, SR=SR,                  &
                             SNOWNC = domain%vars_2d(domain%var_indx(kVARS%snowfall)%v)%data_2d,&
                             SNOWNCV = last_snow,                         &
                             GRAUPELNC = domain%vars_2d(domain%var_indx(kVARS%graupel)%v)%data_2d,          &
                             GRAUPELNCV = last_graup,                    & ! hm added 7/13/13
                             EFFC = domain%vars_3d(domain%var_indx(kVARS%re_cloud)%v)%data_3d,          &
                             EFFI = domain%vars_3d(domain%var_indx(kVARS%re_ice)%v)%data_3d,            &
                             EFFS = domain%vars_3d(domain%var_indx(kVARS%re_snow)%v)%data_3d,           &
                             ISED3D = domain%vars_3d(domain%var_indx(kVARS%ice2_vmi)%v)%data_3d,           &
                             SSED3D = domain%vars_3d(domain%var_indx(kVARS%ice1_vmi)%v)%data_3d,           &
                             refl_10cm = refl_10cm, diagflag = .False.,   &
                             do_radar_ref=0,                              & ! GT added for reflectivity calcs
                             qrcuten=domain%tend%qr,                      &
                             qscuten=domain%tend%qs,                      &
                             qicuten=domain%tend%qi,                      &
                             ids = ids, ide = ide,                   & ! domain dims
                             jds = jds, jde = jde,                   &
                             kds = kds, kde = kde,                   &
                             ims = ims, ime = ime,                   & ! memory dims
                             jms = jms, jme = jme,                   &
                             kms = kms, kme = kme,                   &
                             its = its, ite = ite,                   & ! tile dims
                             jts = jts, jte = jte,                   &
                             kts = kts, kte = kte)
        elseif (options%physics%microphysics==kMP_ISHMAEL) then
            call mp_jensen_ishmael(ITIMESTEP=1,                &  !*                                                                         
                             DT_IN=dt,                           &  !*
                             P=domain%vars_3d(domain%var_indx(kVARS%pressure)%v)%data_3d,                                &  !*
                             DZ=domain%vars_3d(domain%var_indx(kVARS%dz_interface)%v)%data_3d,                            &  !* !
                             TH= domain%vars_3d(domain%var_indx(kVARS%potential_temperature)%v)%data_3d,                              &  !*
                             QV = domain%vars_3d(domain%var_indx(kVARS%water_vapor)%v)%data_3d,             &
                             QC = domain%vars_3d(domain%var_indx(kVARS%cloud_water_mass)%v)%data_3d,        &
                             QR = domain%vars_3d(domain%var_indx(kVARS%rain_mass)%v)%data_3d,               &
                             NR = domain%vars_3d(domain%var_indx(kVARS%rain_number)%v)%data_3d,             &
                             QI1 = domain%vars_3d(domain%var_indx(kVARS%ice_mass)%v)%data_3d,          &
                             NI1 = domain%vars_3d(domain%var_indx(kVARS%ice_number)%v)%data_3d,        &
                             AI1 = domain%vars_3d(domain%var_indx(kVARS%ice1_a)%v)%data_3d,                     &  !*
                             CI1 = domain%vars_3d(domain%var_indx(kVARS%ice1_c)%v)%data_3d,                     &  !*
                             QI2 = domain%vars_3d(domain%var_indx(kVARS%ice2_mass)%v)%data_3d,                     &  !*
                             NI2 = domain%vars_3d(domain%var_indx(kVARS%ice2_number)%v)%data_3d,                     &  !*
                             AI2 = domain%vars_3d(domain%var_indx(kVARS%ice2_a)%v)%data_3d,                     &  !*
                             CI2 = domain%vars_3d(domain%var_indx(kVARS%ice2_c)%v)%data_3d,                     &  !*
                             QI3 = domain%vars_3d(domain%var_indx(kVARS%ice3_mass)%v)%data_3d,                     &  !*
                             NI3 = domain%vars_3d(domain%var_indx(kVARS%ice3_number)%v)%data_3d,                     &  !*
                             AI3 = domain%vars_3d(domain%var_indx(kVARS%ice3_a)%v)%data_3d,                     &  !*
                             CI3 = domain%vars_3d(domain%var_indx(kVARS%ice3_c)%v)%data_3d,                     &  !*
                             IDS=ids,IDE=ide, JDS=jds,JDE=jde, KDS=kds,KDE=kde, &
                             IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=kms,KME=kme, &
                             ITS=its,ITE=ite, JTS=jts,JTE=jte, KTS=kts,KTE=kte, &
                             RAINNC = domain%vars_2d(domain%var_indx(kVARS%precipitation)%v)%data_2d, &
                             RAINNCV = last_rain,                    &
                             SNOWNC = domain%vars_2d(domain%var_indx(kVARS%snowfall)%v)%data_2d,&
                             SNOWNCV = last_snow,                         &
                             diag_effc3d=domain%vars_3d(domain%var_indx(kVARS%re_cloud)%v)%data_3d,               &
                             diag_effi3d=domain%vars_3d(domain%var_indx(kVARS%re_ice)%v)%data_3d,                 &
                             !diag_dbz3d=refl_10cm,               &
                             diag_vmi3d_1=domain%vars_3d(domain%var_indx(kVARS%ice1_vmi)%v)%data_3d,                 &
                             !diag_di3d_1=di3d,                   &
                             diag_rhopo3d_1=domain%vars_3d(domain%var_indx(kVARS%ice1_rho)%v)%data_3d,          &
                             diag_phii3d_1=domain%vars_3d(domain%var_indx(kVARS%ice1_phi)%v)%data_3d,           &
                             diag_vmi3d_2=domain%vars_3d(domain%var_indx(kVARS%ice2_vmi)%v)%data_3d,               &
                             !diag_di3d_2=di3d_2,                 
                             diag_rhopo3d_2=domain%vars_3d(domain%var_indx(kVARS%ice2_rho)%v)%data_3d,          &
                             diag_phii3d_2=domain%vars_3d(domain%var_indx(kVARS%ice2_phi)%v)%data_3d,           &
                             diag_vmi3d_3=domain%vars_3d(domain%var_indx(kVARS%ice3_vmi)%v)%data_3d,               &
                          
                          !diag_di3d_3=di3d_3,                 &
                             diag_rhopo3d_3=domain%vars_3d(domain%var_indx(kVARS%ice3_rho)%v)%data_3d,          &
                             diag_phii3d_3=domain%vars_3d(domain%var_indx(kVARS%ice3_phi)%v)%data_3d            &
                             !diag_itype_1=itype,                 &
                             !diag_itype_2=itype_2,               &
                             !diag_itype_3=itype_3                &
                             )
        elseif (options%physics%microphysics==kMP_WSM6) then
            call wsm6(q = domain%vars_3d(domain%var_indx(kVARS%water_vapor)%v)%data_3d,                      &
                              th = domain%vars_3d(domain%var_indx(kVARS%potential_temperature)%v)%data_3d,            &
                              qc = domain%vars_3d(domain%var_indx(kVARS%cloud_water_mass)%v)%data_3d,                 &
                              qi = domain%vars_3d(domain%var_indx(kVARS%ice_mass)%v)%data_3d,                   &
                              qr = domain%vars_3d(domain%var_indx(kVARS%rain_mass)%v)%data_3d,                        &
                              qs = domain%vars_3d(domain%var_indx(kVARS%snow_mass)%v)%data_3d,                        &
                              qg = domain%vars_3d(domain%var_indx(kVARS%graupel_mass)%v)%data_3d,                     &
                              pii= domain%vars_3d(domain%var_indx(kVARS%exner)%v)%data_3d,                            &
                              p =  domain%vars_3d(domain%var_indx(kVARS%pressure)%v)%data_3d,                         &
                              delz = domain%vars_3d(domain%var_indx(kVARS%dz_interface)%v)%data_3d,                          &
                              den = domain%vars_3d(domain%var_indx(kVARS%density)%v)%data_3d,                   &
                              delt = dt,                                           &
                              g = gravity,                                          &
                              cpd = cp, cpv = cpv, rd = R_d, rv = R_v, t0c = 273.15,          &
                              ep1 = EP_1, ep2 = EP_2, qmin = epsilon,                                &
                              XLS = XLS, XLV0 = XLV, XLF0 = XLF,                    &
                              den0 = rhoair0, denr = rhowater,                  &
                              cliq = cliq, cice = cice, psat = psat,                                   &
                              rain = domain%vars_2d(domain%var_indx(kVARS%precipitation)%v)%data_2d,    &
                              rainncv = last_rain,                                & ! not used outside thompson (yet)
                              sr = SR,                                              & ! not used outside thompson (yet)
                              snow = domain%vars_2d(domain%var_indx(kVARS%snowfall)%v)%data_2d,         &
                              graupel = domain%vars_2d(domain%var_indx(kVARS%graupel)%v)%data_2d,       &
                              ids = ids, ide = ide,                   & ! domain dims
                              jds = jds, jde = jde,                   &
                              kds = kds, kde = kde,                   &
                              ims = ims, ime = ime,                   & ! memory dims
                              jms = jms, jme = jme,                   &
                              kms = kms, kme = kme,                   &
                              its = its, ite = ite,                   & ! tile dims
                              jts = jts, jte = jte,                   &
                              kts = kts, kte = kte)

        elseif (options%physics%microphysics==kMP_WSM3) then

            call wsm3(        q = domain%vars_3d(domain%var_indx(kVARS%water_vapor)%v)%data_3d,                       &
                              th = domain%vars_3d(domain%var_indx(kVARS%potential_temperature)%v)%data_3d,            &
                              qci = domain%vars_3d(domain%var_indx(kVARS%cloud_water_mass)%v)%data_3d,                &
                              qrs = domain%vars_3d(domain%var_indx(kVARS%rain_mass)%v)%data_3d,                       &
                              w = domain%vars_3d(domain%var_indx(kVARS%w_real)%v)%data_3d,                            &
                              pii= domain%vars_3d(domain%var_indx(kVARS%exner)%v)%data_3d,                            &
                              p =  domain%vars_3d(domain%var_indx(kVARS%pressure)%v)%data_3d,                         &
                              delz = domain%vars_3d(domain%var_indx(kVARS%dz)%v)%data_3d,                        &
                              den = domain%vars_3d(domain%var_indx(kVARS%density)%v)%data_3d,                         &
                              delt = dt,                                            &
                              g = gravity,                                          &
                              cpd = cp, cpv = cpv, rd = R_d, rv = R_v, t0c = 273.15,  &
                              ep1 = EP_1, ep2 = EP_2, qmin = epsilon,                 &
                              XLS = XLS, XLV0 = XLV, XLF0 = XLF,                    &
                              den0 = rhoair0, denr = rhowater,                      &
                              cliq = cliq, cice = cice, psat = psat,                &
                              rain = domain%vars_2d(domain%var_indx(kVARS%precipitation)%v)%data_2d,      &
                              rainncv = last_rain,                                & ! not used outside thompson (yet)
                              sr = SR,                                              & ! not used outside thompson (yet)
                              snow = domain%vars_2d(domain%var_indx(kVARS%snowfall)%v)%data_2d,           &
                              snowncv = last_snow,                                  &
                              has_reqc=0, has_reqi=0, has_reqs=0,     &
                              ids = ids, ide = ide,                   & ! domain dims
                              jds = jds, jde = jde,                   &
                              kds = kds, kde = kde,                   &
                              ims = ims, ime = ime,                   & ! memory dims
                              jms = jms, jme = jme,                   &
                              kms = kms, kme = kme,                   &
                              its = its, ite = ite,                   & ! tile dims
                              jts = jts, jte = jte,                   &
                              kts = kts, kte = kte)
        endif

    end subroutine process_subdomain

    subroutine process_halo(domain, options, dt, halo,  &
                            its,ite, jts,jte, kts,kte,  &
                            ims,ime, jms,jme, kms,kme,  &
                            ids,ide, jds,jde, kds,kde)
        implicit none
        type(domain_t), intent(inout) :: domain
        type(options_t),intent(in)    :: options
        real,           intent(in)    :: dt
        integer,        intent(in)    :: halo
        integer,        intent(in)    :: its,ite, jts,jte, kts,kte
        integer,        intent(in)    :: ims,ime, jms,jme, kms,kme
        integer,        intent(in)    :: ids,ide, jds,jde, kds,kde
        integer :: halo_its,halo_ite, halo_jts,halo_jte, halo_kts,halo_kte

        ! process the western halo
        halo_ite = its+halo-1
        call process_subdomain(domain, options, dt, &
                    its,halo_ite, jts,jte, kts,kte, &
                    ims,ime, jms,jme, kms,kme,      &
                    ids,ide, jds,jde, kds,kde)


        ! process the eastern halo
        halo_its = ite-halo+1
        call process_subdomain(domain, options, dt, &
                    halo_its,ite, jts,jte, kts,kte, &
                    ims,ime, jms,jme, kms,kme,      &
                    ids,ide, jds,jde, kds,kde)

        ! for the top and bottom halos, we no longer process the corner elements, so subset i tile
        halo_its = its+halo
        halo_ite = ite-halo

        ! process the southern halo
        halo_jte = jts+halo-1
        call process_subdomain(domain, options, dt,           &
                    halo_its,halo_ite, jts,halo_jte, kts,kte, &
                    ims,ime, jms,jme, kms,kme,                &
                    ids,ide, jds,jde, kds,kde)


        ! process the northern halo
        halo_jts = jte-halo+1
        call process_subdomain(domain, options, dt,           &
                    halo_its,halo_ite, halo_jts,jte, kts,kte, &
                    ims,ime, jms,jme, kms,kme,                &
                    ids,ide, jds,jde, kds,kde)


    end subroutine process_halo


    !>----------------------------------------------------------
    !! Microphysical driver
    !!
    !! This routine handles calling the individual microphysics routine specified, that
    !! includes creating and passing any temporary variables, and checking when to update
    !! the microphysics based on the specified update_interval.
    !!
    !! @param   domain      ICAR model domain structure
    !! @param   options     ICAR model options structure
    !! @param   dt_in       Current driving time step (this is the advection step)
    !!
    !!----------------------------------------------------------
    subroutine mp(domain, options, dt_in, halo, subset)
        implicit none
        type(domain_t), intent(inout) :: domain
        type(options_t),intent(in)    :: options
        real,           intent(in)    :: dt_in
        integer,        intent(in),   optional :: halo, subset

        real :: mp_dt
        integer ::ids,ide,jds,jde,kds,kde, itimestep=1
        integer ::ims,ime,jms,jme,kms,kme
        integer ::its,ite,jts,jte,kts,kte, nx,ny,nz

        if(options%physics%microphysics==0) return
        
        ids = domain%grid%ids;   ims = domain%grid%ims;   its = domain%grid%its
        ide = domain%grid%ide;   ime = domain%grid%ime;   ite = domain%grid%ite
        nx  = domain%grid%nx
        kds = domain%grid%kds;   kms = domain%grid%kms;   kts = domain%grid%kts
        kde = domain%grid%kde;   kme = domain%grid%kme;   kte = domain%grid%kte
        nz  = domain%grid%nz
        jds = domain%grid%jds;   jms = domain%grid%jms;   jts = domain%grid%jts
        jde = domain%grid%jde;   jme = domain%grid%jme;   jte = domain%grid%jte
        ny  = domain%grid%ny

        ! if this is the first time mp is called, set last time such that mp will update
        if (last_model_time==-999) then
            last_model_time = (domain%sim_time%seconds() - max(real(update_interval), dt_in))
        endif


        ! only run the microphysics if the next time step would put it over the update_interval time
        if (((domain%sim_time%seconds() + dt_in)-last_model_time)>=update_interval) then

            ! calculate the actual time step for the microphysics
            mp_dt = domain%sim_time%seconds()-last_model_time


            ! set the current tile to the top layer to process microphysics for
            if (options%mp%top_mp_level>0) then
                kte = min(kte, options%mp%top_mp_level)
            endif

            ! reset the counter so we know that *this* is the last time we've run the microphysics
            ! NOTE, ONLY reset this when running the inner subset... ideally probably need a separate counter for the halo and subset
            !last_model_time = domain%sim_time%seconds()
            
            call calc_w_real(domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d,      &
                         domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d,      &
                         domain%vars_3d(domain%var_indx(kVARS%w)%v)%data_3d,      &
                         domain%vars_3d(domain%var_indx(kVARS%w_real)%v)%data_3d,      &
                         domain%vars_3d(domain%var_indx(kVARS%dzdx_u)%v)%data_3d, domain%vars_3d(domain%var_indx(kVARS%dzdy_v)%v)%data_3d, domain%vars_3d(domain%var_indx(kVARS%dzdx)%v)%data_3d, domain%vars_3d(domain%var_indx(kVARS%dzdy)%v)%data_3d,   &
                         domain%vars_3d(domain%var_indx(kVARS%jacobian)%v)%data_3d)
                             

            if (present(subset)) then
                last_model_time = domain%sim_time%seconds()
                call process_subdomain(domain, options, mp_dt,                 &
                                       its = its + subset, ite = ite - subset, &
                                       jts = jts + subset, jte = jte - subset, &
                                       kts = kts,          kte = kte,          &
                                       ims = ims, ime = ime,                   & ! memory dims
                                       jms = jms, jme = jme,                   &
                                       kms = kms, kme = kme,                   &
                                       ids = ids, ide = ide,                   & ! domain dims
                                       jds = jds, jde = jde,                   &
                                       kds = kds, kde = kde)
            endif

            if (present(halo)) then
                call process_halo(domain, options, mp_dt, halo, &
                                       its = its, ite = ite,    &
                                       jts = jts, jte = jte,    &
                                       kts = kts, kte = kte,    &
                                       ims = ims, ime = ime,    & ! memory dims
                                       jms = jms, jme = jme,    &
                                       kms = kms, kme = kme,    &
                                       ids = ids, ide = ide,    & ! domain dims
                                       jds = jds, jde = jde,    &
                                       kds = kds, kde = kde)

            endif

            if ((.not.present(halo)).and.(.not.present(subset))) then
                last_model_time = domain%sim_time%seconds()                             
                call process_subdomain(domain, options, mp_dt,  &
                                        its = its, ite = ite,    &
                                        jts = jts, jte = jte,    &
                                        kts = kts, kte = kte,    &
                                        ims = ims, ime = ime,    & ! memory dims
                                        jms = jms, jme = jme,    &
                                        kms = kms, kme = kme,    &
                                        ids = ids, ide = ide,    & ! domain dims
                                        jds = jds, jde = jde,    &
                                        kds = kds, kde = kde)

            endif

        endif

    end subroutine mp

    !>----------------------------------------------------------
    !! Finalize microphysical routines
    !!
    !! This routine will call the finalization routines (if any) for the specified microphysics packages.
    !! It also deallocates any module level variables, e.g. SR
    !!
    !! @param   options     ICAR model options to specify required initializations
    !!
    !!----------------------------------------------------------
    subroutine mp_finish(options)
        implicit none
        type(options_t),intent(in)::options

        if (allocated(SR)) then
            deallocate(SR)
        endif
        if (allocated(last_rain)) then
            deallocate(last_rain)
        endif
        if (allocated(last_snow)) then
            deallocate(last_snow)
        endif
        if (allocated(last_graup)) then
            deallocate(last_graup)
        endif

    end subroutine mp_finish
end module microphysics
