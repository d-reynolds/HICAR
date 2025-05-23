!> ----------------------------------------------------------------------------
!!  Driver to call different convection schemes
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!! ----------------------------------------------------------------------------
module convection
    use icar_constants
    use mod_wrf_constants,  only: svpt0, EP_1, EP_2, gravity, XLS, XLV, cp, R_d
    use module_cu_tiedtke,  only: tiedtkeinit, CU_TIEDTKE
    ! use module_cu_kf,       only: kfinit, KFCPS
    use module_cu_nsas,  only: nsasinit, cu_nsas
    use module_cu_bmj,  only: bmjinit, bmjdrv

    use options_interface,   only : options_t
    use domain_interface,    only : domain_t

    implicit none
    private
    public :: init_convection, convect, cu_var_request

    logical,allocatable, dimension(:,:) ::  CU_ACT_FLAG
    real,   allocatable, dimension(:,:) ::  XLAND, RAINCV, PRATEC, NCA
    real,   allocatable, dimension(:,:,:):: W0AVG, w_stochastic
    real,   allocatable, dimension(:,:)::   lowest_convection_layer, highest_convection_layer! , hpbl  ! pbl height I guess? Used in NSAS.
    real,   allocatable, dimension(:,:)::   CLDEFI ! precipitation efficiency (for BMJ scheme) (dimensionless)
    real,   allocatable, dimension(:,:)::   CUBOT,CUTOP, CONVCLD  ! for BMJ scheme
    real,   allocatable, dimension(:,:,:):: CCLDFRA, QCCONV, QICONV ! for BMJ scheme
    integer :: ids,ide,jds,jde,kds,kde ! Domain dimensions
    integer :: ims,ime,jms,jme,kms,kme ! Local Memory dimensions
    integer :: its,ite,jts,jte,kts,kte ! Processing Tile dimensions
    integer, allocatable, dimension(:,:) :: lowlyr   ! for bmj scheme!, kpbl
    logical :: bmj_rad_feedback

contains

    subroutine cu_var_request(options)
        implicit none
        type(options_t),intent(inout) :: options

        if (options%physics%convection == kCU_TIEDTKE) then
            call options%alloc_vars( &
                         [kVARS%water_vapor, kVARS%potential_temperature, kVARS%temperature,                        &
                         kVARS%cloud_water_mass, kVARS%ice_mass, kVARS%precipitation, kVARS%convective_precipitation,   &
                         kVARS%exner, kVARS%dz_interface, kVARS%density, kVARS%pressure_interface, kVARS%pressure,  &
                         kVARS%sensible_heat, kVARS%latent_heat, kVARS%u, kVARS%v, kVARS%w, kVARS%land_mask,        &
                         kVARS%tend_qv, kVARS%tend_th, kVARS%tend_qc, kVARS%tend_qi, kVARS%tend_qs, kVARS%tend_qr,  &
                         kVARS%tend_u, kVARS%tend_v, kVARS%tend_qv_pbl, kVARS%tend_qv_adv, kVARS%QFX])

             call options%advect_vars([kVARS%potential_temperature, kVARS%water_vapor])

             call options%restart_vars( &
                         [kVARS%water_vapor, kVARS%potential_temperature, kVARS%temperature, kVARS%QFX,             &
                         kVARS%cloud_water_mass, kVARS%ice_mass, kVARS%precipitation, kVARS%convective_precipitation,   &
                         kVARS%sensible_heat, kVARS%latent_heat, kVARS%u, kVARS%v, kVARS%pressure])

        else if (options%physics%convection == kCU_NSAS) then  ! Not checked thoroughly yet
            call options%alloc_vars( &
                         [kVARS%water_vapor, kVARS%potential_temperature, kVARS%temperature, kVARS%QFX,             &
                         kVARS%cloud_water_mass, kVARS%ice_mass, kVARS%precipitation, kVARS%convective_precipitation,   &
                         kVARS%exner, kVARS%dz_interface, kVARS%density, kVARS%pressure_interface, kVARS%pressure,  &
                         kVARS%sensible_heat, kVARS%latent_heat, kVARS%u, kVARS%v, kVARS%w, kVARS%land_mask,        &
                         kVARS%tend_qv, kVARS%tend_th, kVARS%tend_qc, kVARS%tend_qi, kVARS%tend_qs, kVARS%tend_qr,  &
                         kVARS%tend_u, kVARS%tend_v, kVARS%hpbl])

             call options%advect_vars([kVARS%potential_temperature, kVARS%water_vapor])

             call options%restart_vars( &
                         [kVARS%water_vapor, kVARS%potential_temperature, kVARS%temperature, kVARS%QFX,             &
                         kVARS%cloud_water_mass, kVARS%ice_mass, kVARS%precipitation, kVARS%convective_precipitation,   &
                         kVARS%sensible_heat, kVARS%latent_heat, kVARS%u, kVARS%v, kVARS%pressure, kVARS%hpbl])

        else if (options%physics%convection == kCU_BMJ) then  ! Not checked thoroughly yet
            call options%alloc_vars( &
                         [kVARS%water_vapor, kVARS%potential_temperature, kVARS%temperature,                        &
                         kVARS%cloud_water_mass, kVARS%ice_mass, kVARS%precipitation, kVARS%convective_precipitation,   &
                         kVARS%exner, kVARS%dz_interface, kVARS%density, kVARS%pressure_interface, kVARS%pressure,  &
                         kVARS%sensible_heat, kVARS%latent_heat, kVARS%u, kVARS%v, kVARS%w, kVARS%land_mask,        &
                         kVARS%tend_qv, kVARS%tend_th, kVARS%tend_qc, kVARS%tend_qi, kVARS%tend_qs, kVARS%tend_qr,  &
                         kVARS%tend_u, kVARS%tend_v, kVARS%kpbl])

             call options%advect_vars([kVARS%potential_temperature, kVARS%water_vapor])

             call options%restart_vars( &
                         [kVARS%water_vapor, kVARS%potential_temperature, kVARS%temperature,                        &
                         kVARS%cloud_water_mass, kVARS%ice_mass, kVARS%precipitation, kVARS%convective_precipitation,   &
                         kVARS%sensible_heat, kVARS%latent_heat, kVARS%u, kVARS%v, kVARS%pressure, kVARS%kpbl])

        endif


    end subroutine cu_var_request


    subroutine init_convection(domain,options, context_chng)
        implicit none
        type(domain_t),  intent(inout) :: domain
        type(options_t), intent(in) :: options
        logical, optional, intent(in) :: context_chng

        logical :: context_change, restart
        integer :: i, j

        if (present(context_chng)) then
            context_change = context_chng
        else
            context_change = .false.
        endif

        restart = context_change .or. options%restart%restart

        if (STD_OUT_PE .and. .not.context_change) write(*,*) "Initializing Cumulus Scheme"

        ! module level variables for easy access... need to think about tiling to permit halo processing separately.
        ids = domain%grid%ids
        ide = domain%grid%ide
        jds = domain%grid%jds
        jde = domain%grid%jde
        kds = domain%grid%kds
        kde = domain%grid%kde
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

        if (options%physics%convection > 0) then
            if (allocated(CU_ACT_FLAG)) deallocate(CU_ACT_FLAG)
            allocate(CU_ACT_FLAG(ims:ime,jms:jme), source=.False.)
            CU_ACT_FLAG(its:ite, jts:jte) = .True.

            if (allocated(RAINCV)) deallocate(RAINCV)
            if (allocated(PRATEC)) deallocate(PRATEC)
            allocate(RAINCV(ims:ime,jms:jme))
            allocate(PRATEC(ims:ime,jms:jme))

            RAINCV = 0
            PRATEC = 0
            
            if (allocated(W0AVG)) deallocate(W0AVG)
            allocate(W0AVG(ims:ime,kms:kme,jms:jme))
            W0AVG = 0

            if (allocated(XLAND)) deallocate(XLAND)
            allocate(XLAND(ims:ime,jms:jme))
            XLAND = domain%vars_2d(domain%var_indx(KVARS%land_mask)%v)%data_2di
            where(domain%vars_2d(domain%var_indx(KVARS%land_mask)%v)%data_2di == 0) XLAND = 2 ! 0 is water if using "LANDMASK" as input

            if (allocated(w_stochastic)) deallocate(w_stochastic)
            allocate(w_stochastic(ims:ime,kms:kme,jms:jme))

        endif

        if (options%physics%convection == kCU_TIEDTKE) then
            if (STD_OUT_PE .and. .not.context_change) write(*,*) "    Tiedtke Cumulus scheme"

            call tiedtkeinit(domain%tend%th,domain%tend%qv,   &
                             domain%tend%qc,domain%tend%qi,   &
                             domain%tend%u, domain%tend%v,    &
                             restart,1,1,0,            &
                             .true.,                          &
                             ids,ide, jds,jde, kds,kde,       &
                             ims,ime, jms,jme, kms,kme,       &
                             its,ite, jts,jte, kts,kte  )

         ! elseif (options%physics%convection==kCU_KAINFR) then
         !     write(*,*) "    Kain-Fritsch Cumulus scheme"
         !
         !     allocate(NCA(ids:ide,jds:jde))
         !     NCA=0
         !     call kfinit(domain%tend%th,domain%tend%qv,   &
         !                 domain%tend%qc,domain%tend%qr,   &
         !                 domain%tend%qi,domain%tend%qs,   &
         !                 NCA,W0AVG,1,1,                   & !...,P_QI,P_QS
         !                 0,.false.,.true.,                & ! P_FIRST_SCALAR, restart, allowed_to_read
         !                 ids, ide, jds, jde, kds, kde,  &
         !                 ids, ide, jds, jde, kds, kde,    &  ! !!!!!!!!!!!! kme iso kde!
         !                 ids, ide, jds, jde, kds, kde-1)


        else if (options%physics%convection==kCU_NSAS) then
            if (STD_OUT_PE .and. .not.context_change) write(*,*) "     NSAS Cumulus scheme"

            if (allocated(lowest_convection_layer)) deallocate(lowest_convection_layer)
            if (allocated(highest_convection_layer)) deallocate(highest_convection_layer)
            allocate(lowest_convection_layer(ims:ime,jms:jme))
            allocate(highest_convection_layer(ims:ime,jms:jme))

            ! NSAS expects a PBL height, But since ICAR does not currently compute this, fix at rough esitmate of 1000m ?
            if (options%physics%boundarylayer/=kPBL_YSU) then
                do i=ims,ime
                    do j=jms,jme
                        domain%vars_2d(domain%var_indx(kVARS%hpbl)%v)%data_2d(i,j)=1000.  ! or take terrain height into account?
                    enddo
                enddo
            endif

            call nsasinit(  rthcuten=domain%tend%th                             &
                            ,rqvcuten=domain%tend%qv                            &
                            ,rqccuten=domain%tend%qc                            &
                            ,rqicuten=domain%tend%qi                            &
                            ,rucuten=domain%tend%u                              &
                            ,rvcuten=domain%tend%v                              &
                            ,restart=restart             &  ! options%restart                            &
                            ,p_qc=1                                             & ! copied from Tiedke; no idea as to what these 3 flags represent.
                            ,p_qi=1                                             &
                            ,p_first_scalar=0                                   &
                            ,allowed_to_read=.true.                             &
                            ,ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde   &
                            ,ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme   &
                            ,its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte &  ! -1 or not?
                            )

        else if (options%physics%convection==kCU_BMJ) then
            if (STD_OUT_PE .and. .not.context_change) write(*,*) "     BMJ Cumulus scheme"

            if (allocated(CLDEFI)) deallocate(CLDEFI)
            if (allocated(CUBOT)) deallocate(CUBOT)
            if (allocated(CUTOP)) deallocate(CUTOP)
            if (allocated(CCLDFRA)) deallocate(CCLDFRA)
            if (allocated(QCCONV)) deallocate(QCCONV)
            if (allocated(QICONV)) deallocate(QICONV)
            if (allocated(CONVCLD)) deallocate(CONVCLD)
            if (allocated(lowlyr)) deallocate(lowlyr)
            allocate( CLDEFI(ims:ime,jms:jme) )
            allocate( CUBOT(IMS:IME,JMS:JME) )
            allocate( CUTOP(IMS:IME,JMS:JME) )  ! INTENT:OUT
            allocate( CCLDFRA(ims:ime,jms:jme,kms:kme) )
            allocate( QCCONV(ims:ime,jms:jme,kms:kme) )
            allocate( QICONV(ims:ime,jms:jme,kms:kme) )
            allocate( CONVCLD(ims:ime,jms:jme) )
            allocate( lowlyr(ims:ime,jms:jme) )
            ! allocate( kpbl(ims:ime,jms:jme) ) ! domain%vars_2d(domain%var_indx(kVARS%kpbl)%v)%data_2d when YSU is selected


            ! ----- the values below are just educated (or not) guesses. Need to refine.  -------
            do i=ims,ime
                do j=jms,jme
                    lowlyr(i,j) = 1 ! ???
                    if (options%physics%boundarylayer/=kPBL_YSU) domain%vars_2d(domain%var_indx(kVARS%kpbl)%v)%data_2d(i,j) = 10 ! without YSU, we assign a stationary value to kpbl
                enddo
            enddo



            bmj_rad_feedback = .true.  !??

            call BMJINIT( rthcuten=domain%tend%th                               &
                         ,rqvcuten=domain%tend%qv                               &
                         ,rqccuten=domain%tend%qc                               &
                         ,rqrcuten=domain%tend%qr                               & !-- RQRCUTEN      Qr tendency due to cumulus scheme precipitation (kg/kg/s)
                         ,cldefi=CLDEFI                                         & !-- CLDEFI        precipitation efficiency (for BMJ scheme) (dimensionless)
                         ,lowlyr=lowlyr                                         & !-- LOWLYR        index of lowest model layer above the ground
                         ,cp=cp                                                 & ! cp  = 1012.0    ! J/kg/K   specific heat capacity of moist STP air?                 CP
                         ,RD=r_d                      & ! r_d (wrf_constatns) = 287., while Rd  = 287.058 (icar_constants)
                         ,RESTART=restart                                      &
                         ,allowed_to_read=.true.                                & !
                         ,ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde  &
                         ,ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme  &
                         ,its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte  &  ! -1 or not?
                         )

         endif

         if ((options%cu%stochastic_cu /= kNO_STOCHASTIC) .and.  (STD_OUT_PE .and. .not.context_change)) then
            write(*,*)"      stochastic W pertubation for convection triggering"  ! to check that it actually turns on/of
         elseif ((options%cu%stochastic_cu == kNO_STOCHASTIC) .and.  (STD_OUT_PE .and. .not.context_change)) then
            write(*,*)"      No stochastic W pertubation for convection triggering"
         endif

    end subroutine init_convection

subroutine convect(domain,options,dt_in)
    implicit none
    type(domain_t),  intent(inout) :: domain
    type(options_t), intent(in)    :: options
    real, intent(in) :: dt_in

    integer :: j,itimestep,STEPCU, dx_factor_nsas
    real :: internal_dt
    real, allocatable :: znu(:)

    if (options%physics%convection==0) return

    itimestep = 1
    STEPCU = 1

    !$omp parallel private(j) &
    !$omp default(shared)
    !$omp do schedule(static)
    do j=jts,jte
        RAINCV(:,j)     = 0
        domain%tend%qv(:,:,j) = 0
        domain%tend%th(:,:,j) = 0
        domain%tend%qc(:,:,j) = 0
        domain%tend%qi(:,:,j) = 0
        domain%tend%qs(:,:,j) = 0
        domain%tend%qr(:,:,j) = 0
        domain%tend%u(:,:,j)  = 0
        domain%tend%v(:,:,j)  = 0
    enddo
    !$omp end do
    !$omp end parallel

    ! Stochastic disturbance for  vertical speeds:
    if (options%cu%stochastic_cu /= kNO_STOCHASTIC) then
        call random_number(w_stochastic)
        block
            integer :: i
            do i=1,size(w_stochastic,2)
                w_stochastic(:,i,:) = domain%vars_2d(domain%var_indx(kVARS%sensible_heat)%v)%data_2d/500 + w_stochastic(:,i,:)
            enddo
        end block
        W0AVG = domain%vars_3d(domain%var_indx(kVARS%w_real)%v)%data_3d+(w_stochastic * options%cu%stochastic_cu - options%cu%stochastic_cu*0.75) ! e.g. * 20 - 15)
    else
        W0AVG = domain%vars_3d(domain%var_indx(kVARS%w_real)%v)%data_3d
    endif


    if (options%physics%convection==kCU_TIEDTKE) then

        call calc_znu(domain,znu)
        call CU_TIEDTKE(                                        &
                 dt_in, itimestep, STEPCU                       &
                ,RAINCV, PRATEC                                 &
                ,domain%vars_2d(domain%var_indx(kVARS%QFX)%v)%data_2d                             &
                ,domain%vars_2d(domain%var_indx(kVARS%sensible_heat)%v)%data_2d                   &
                ,znu                                     &
                ,domain%vars_3d(domain%var_indx(kVARS%u_mass)%v)%data_3d                          &
                ,domain%vars_3d(domain%var_indx(kVARS%v_mass)%v)%data_3d                          &
                ,W0AVG                                          &
                ,domain%vars_3d(domain%var_indx(kVARS%temperature)%v)%data_3d                     &
                ,domain%vars_3d(domain%var_indx(kVARS%water_vapor)%v)%data_3d                     &
                ,domain%vars_3d(domain%var_indx(kVARS%cloud_water_mass)%v)%data_3d                &
                ,domain%vars_3d(domain%var_indx(kVARS%ice_mass)%v)%data_3d                  &
                ,domain%vars_3d(domain%var_indx(kVARS%exner)%v)%data_3d                           &
                ,domain%vars_3d(domain%var_indx(kVARS%density)%v)%data_3d                         &
                ,domain%tend%qv_adv                             &
                ,domain%tend%qv_pbl                             &
                ,domain%vars_3d(domain%var_indx(kVARS%dz_interface)%v)%data_3d                    & ! dz8w
                ,domain%vars_3d(domain%var_indx(kVARS%pressure)%v)%data_3d                        & !
                ,domain%vars_3d(domain%var_indx(kVARS%pressure_interface)%v)%data_3d              &
                ,XLAND, CU_ACT_FLAG                             &
                ,ids,ide, jds,jde, kds,kde                      &
                ,ims,ime, jms,jme, kms,kme                      &
                ,its,ite, jts,jte, kts,kte-1                    &
                ,domain%tend%th, domain%tend%qv, domain%tend%qc &
                ,domain%tend%qi, domain%tend%u,  domain%tend%v  &
                ,.True.,.True.,.True.,.True.,.True.             &
                )

    ! elseif (options%physics%convection==kCU_KAINFR) then
    !     call KFCPS(                                          &
    !            ids,ide, jds,jde, kds,kde                     & ! domain
    !           ,ids,ide, jds,jde, kds,kde                     & ! memory
    !           ,ids+1,ide-1, jds+1,jde-1, kds,kde             & ! tile
    !           ,dt_in,itimestep,domain%dx,dt_in,.false.       & ! dt KTAU, dx, cu_dt, adapt_step
    !           ,domain%rho                                    & ! rho
    !           ,RAINCV,PRATEC,NCA                             & ! convective rain, convective rain rate,
    !           ,domain%Um,domain%Vm,domain%th,domain%t        &
    !           ,domain%w_real                                 &
    !           ,domain%qv,domain%dz_inter,domain%p,domain%pii &
    !           ,W0AVG                                         & ! "average" vertical velocity (computed internally from w)
    !           ,XLV0,XLV1,XLS0,XLS1,cp,Rd,gravity,EP1         & ! physical "constants"
    !           ,EP2,SVP1,SVP2,SVP3,SVPT0                      & ! physical "constants"
    !           ,STEPCU,CU_ACT_FLAG,.false.                    & ! number of steps between CU calls, boolean grid to act on , warm_rain_only
    !         ! optional arguments
    !         ! ,F_QV    ,F_QC    ,F_QR    ,F_QI    ,F_QS      &
    !           ,.True., .True.,  .True.,  .True.,  .True.     &
    !           ,domain%tend%qv, domain%tend%qc                &
    !           ,domain%tend%qr, domain%tend%qi                &
    !           ,domain%tend%qs, domain%tend%th)

    elseif (options%physics%convection==kCU_NSAS) then

        ! set dx_factor_nsas (cu_nsas.f90 line 567):
        if (domain%dx.le.1000) then
            dx_factor_nsas=1  !       if (dx_factor_nsas == 1) then   dx_factor = 250. / delx ! assume 2.5 ms-1 (1km) and 1.125 cms-1 (200km) (cu_nsas.f90 line 567)
        else
            dx_factor_nsas=2
        endif

        call cu_nsas( &
            dt=dt_in               & !-- dt          time step (s)
            ,dx=domain%dx           & !-- dx          grid interval (m)
            ,p3di=domain%vars_3d(domain%var_indx(kVARS%pressure_interface)%v)%data_3d              & !-- p3di        3d pressure (pa) at interface level
            ,p3d=domain%vars_3d(domain%var_indx(kVARS%pressure)%v)%data_3d                         & !-- p3d         3d pressure (pa)
            ,pi3d=domain%vars_3d(domain%var_indx(kVARS%exner)%v)%data_3d                           & !-- pi3d        3d exner function (dimensionless)
            ,qc3d=domain%vars_3d(domain%var_indx(kVARS%cloud_water_mass)%v)%data_3d                & !-- qc3d        cloud water mixing ratio (kg/kg)
            ,qi3d=domain%vars_3d(domain%var_indx(kVARS%ice_mass)%v)%data_3d                  & !-- qi3d        cloud ice mixing ratio (kg/kg)
            ,rho3d=domain%vars_3d(domain%var_indx(kVARS%density)%v)%data_3d                        &
            ,itimestep=itimestep                                 &
            ,stepcu=STEPCU                                       &
            ,hbot=lowest_convection_layer                        & !-- HBOT          index of lowest model layer with convection  intent(out) ???
            ,htop=highest_convection_layer                       & !-- HTOP          index of highest model layer with convection
            ,cu_act_flag=CU_ACT_FLAG                             &
            ,rthcuten=domain%tend%th                             &
            ,rqvcuten=domain%tend%qv                             & !-- RQVCUTEN      Qv tendency due to cumulus scheme precipitation (kg/kg/s)
            ,rqccuten=domain%tend%qc                             & !--   ...etc...
            ,rqicuten=domain%tend%qi                             &
            ,rucuten=domain%tend%u                               &
            ,rvcuten=domain%tend%v                               &
            ,qv3d=domain%vars_3d(domain%var_indx(kVARS%water_vapor)%v)%data_3d                     & !-- qv3d        3d water vapor mixing ratio (kg/kg)
            ,t3d=domain%vars_3d(domain%var_indx(kVARS%temperature)%v)%data_3d                      & !-- t3d         temperature (k)
            ,raincv=RAINCV                                       & !-- raincv      cumulus scheme precipitation (mm)
            ,pratec=PRATEC                                       &
            ,xland=XLAND                                         &
            ,dz8w=domain%vars_3d(domain%var_indx(kVARS%dz_interface)%v)%data_3d                    & !-- dz8w        dz between full levels (m)
            ,w=W0AVG                                             & !-- w           vertical velocity (m/s)
            ,u3d=domain%vars_3d(domain%var_indx(kVARS%u_mass)%v)%data_3d                           & !-- u3d         3d u-velocity interpolated to theta points (m/s)
            ,v3d=domain%vars_3d(domain%var_indx(kVARS%v_mass)%v)%data_3d                           & !-- v3d         3d v-velocity interpolated to theta points (m/s)
            ,hpbl=domain%vars_2d(domain%var_indx(kVARS%hpbl)%v)%data_2d                            &  ! ?? PBL height I assume?
            ,hfx=domain%vars_2d(domain%var_indx(kVARS%sensible_heat)%v)%data_2d                     & !  HFX  - net upward heat flux at the surface (W/m^2)
            ,qfx=domain%vars_2d(domain%var_indx(kVARS%QFX)%v)%data_2d                              & !  QFX  - net upward moisture flux at the surface (kg/m^2/s)
            ,mp_physics=5        & ! - sets ncloud: - integer no_cloud(0),no_ice(1),cloud+ice(2) (see ln 141 cu_nsas.f90)  use options%physics%microphysics?
            ,dx_factor_nsas=dx_factor_nsas                       & !
            ,p_qc=1                                              & ! copied from Tiedke; no idea as to what these 3 flags represent.
            ,p_qi=1                                              &
            ,p_first_scalar=0                                    &
            ! ,pgcon=""                                            & ! pgcon_use  = 0.55  ! Zhang & Wu (2003,JAS) ! 0.55 is a physically-based value used by GFS
            ,cp=cp                                               & ! cp  = 1012.0    ! J/kg/K   specific heat capacity of moist STP air?
            ,cliq=4190.                                          & ! from wrf_constants.f90
            ,cpv=4.*461.6                                        & ! from wrf_constants.f90
            ,g=gravity                                           &
            ,xlv=2.5E6                                           & ! from wrf_constants.f90
            ,r_d=287.                                            & ! from wrf_constants.f90
            ,r_v=461.6                                           & ! from wrf_constants.f90
            ,ep_1=EP_1                                  & ! EP_1=R_v/R_d-1. (wrf_constants.f90) 461.6/287.-1.
            ,ep_2=EP_2                                     & ! EP_2=R_d/R_v  = 287./461.6
            ,cice=2106.                                           & ! from wrf_constants.f90                                          &
            ,xls=2.85E6                                           & ! from wrf_constants.f90
            ,psat=610.78                                           & ! from wrf_constants.f90
            ,f_qi=.true.                                             & ! optional arg
            ,f_qc=.true.                                             & ! optional arg
            ,ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde,  &
            ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme,  &
            its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte-1 )
        ! ! ! Definitions (for coupling - temporary):
        ! - - - - - - - - - - - - - - - - -
        ! microphysics scheme --> ncloud
        !
        !      !! ncloud   - integer no_cloud(0),no_ice(1),cloud+ice(2)
        !    if (mp_physics .eq. 0) then
        !     ncloud = 0
        !   elseif ( mp_physics .eq. 1 .or. mp_physics .eq. 3 ) then
        !     ncloud = 1
        !   else
        !     ncloud = 2
        !   endif
        ! - - - - - - - - - - - - - - - - -


    elseif (options%physics%convection==kCU_BMJ) then
        CALL BMJDRV(                                            &
                TH=domain%vars_3d(domain%var_indx(kVARS%potential_temperature)%v)%data_3d         &
                ,T=domain%vars_3d(domain%var_indx(kVARS%temperature)%v)%data_3d                   &
                ,RAINCV=RAINCV                                  &
                ,PRATEC=PRATEC                                  &
                ,rho=domain%vars_3d(domain%var_indx(kVARS%density)%v)%data_3d                     &
                ,dt=dt_in                                       & !-- dt          time step (s)
                ,ITIMESTEP=itimestep ,STEPCU=stepcu             &
                ,CUTOP=CUTOP, CUBOT=CUBOT                       &
                ,KPBL=domain%vars_2d(domain%var_indx(kVARS%kpbl)%v)%data_2di                                      &!-- KPBL INTENT(IN)          layer index of the PBL?
                ,dz8w=domain%vars_3d(domain%var_indx(kVARS%dz_interface)%v)%data_3d               &!-- dz8w          dz between full levels (m)
                ,PINT=domain%vars_3d(domain%var_indx(kVARS%pressure_interface)%v)%data_3d         &  !-- p8w           pressure at full levels (Pa)
                ,PMID=domain%vars_3d(domain%var_indx(kVARS%pressure)%v)%data_3d                   &!-- Pcps        3D hydrostatic pressure at half levels (Pa)
                ,PI=domain%vars_3d(domain%var_indx(kVARS%exner)%v)%data_3d                        &   ! exner
                ,CP=cp ,R=r_d ,ELWV=xlv ,ELIV=xls ,G=gravity    &
                ,TFRZ=svpt0 ,D608=EP_1 ,CLDEFI=cldefi           &
                ,LOWLYR=lowlyr                                  &
                ,XLAND=XLAND                                    &
                ,CU_ACT_FLAG=CU_ACT_FLAG                        &
                ,QV=domain%vars_3d(domain%var_indx(kVARS%water_vapor)%v)%data_3d                  &
                ,CCLDFRA=CCLDFRA                                &                 !-- CCLDFRA       convective cloud fraction (for BMJ scheme)
                ,CONVCLD=CONVCLD                                &  !-- CONVCLD       Convective cloud (for BMJ scheme) (kg/m^2)
                ,QCCONV=qcconv                                  &! QCCONV     convective cloud mixing ratio (kg/kg)
                ,QICONV=qiconv                                  &! QICONV     convective ice mixing ratio (kg/kg)
                ,BMJ_RAD_FEEDBACK=bmj_rad_feedback                 &
                ,IDS=ids,IDE=ide,JDS=jds,JDE=jde,KDS=kds,KDE=kde   &
                ,IMS=ims,IME=ime,JMS=jms,JME=jme,KMS=kms,KME=kme   &
                ,ITS=its,ITE=ite,JTS=jts,JTE=jte,KTS=kts,KTE=kte-1 &
                ! optionals
                ,rthcuten=domain%tend%th                           &  !-- RQVCUTEN      Qv tendency due to cumulus scheme precipitation (kg/kg/s)
                ,rqvcuten=domain%tend%qv                           &  !-- RTHCUTEN      Theta tendency due to cumulus scheme precipitation (K/s)
                )


    endif

    ! add domain wide tendency terms
    ! use a separate dt to make it easier to apply on a different dt
    internal_dt = dt_in

    ! if (options%physics%convection==kCU_TIEDTKE) then
    if ( &
        (options%physics%convection==kCU_TIEDTKE) .or.  &
        (options%physics%convection==kCU_NSAS) .or.     &
        (options%physics%convection==kCU_BMJ)           &
        ) then
        ! $omp parallel private(j) &
        ! $omp default(shared)
        ! $omp do schedule(static)
        do j=jts,jte
            if (options%cu%tendency_fraction > 0) then
                ! Here we adjust the tendencies calculated by the cumulus scheme based on the namelist setting tendency_fraction (cu_parameters)
                !-- RTHCUTEN      Theta tendency due to cumulus scheme precipitation (K/s)  (domain%tend%qv)
                !-- RQVCUTEN      Qv tendency due to cumulus scheme precipitation (kg/kg/s)  (domain%tend%th)
                ! -- ...etc
                if (options%cu%tend_qv_fraction > 0) domain%vars_3d(domain%var_indx(kVARS%water_vapor)%v)%data_3d(:,:,j)           = domain%vars_3d(domain%var_indx(kVARS%water_vapor)%v)%data_3d(:,:,j)           + domain%tend%qv(:,:,j)*internal_dt * options%cu%tend_qv_fraction
                if (options%cu%tend_qc_fraction > 0) domain%vars_3d(domain%var_indx(kVARS%cloud_water_mass)%v)%data_3d(:,:,j)      = domain%vars_3d(domain%var_indx(kVARS%cloud_water_mass)%v)%data_3d(:,:,j)      + domain%tend%qc(:,:,j)*internal_dt * options%cu%tend_qc_fraction
                if (options%cu%tend_th_fraction > 0) domain%vars_3d(domain%var_indx(kVARS%potential_temperature)%v)%data_3d(:,:,j) = domain%vars_3d(domain%var_indx(kVARS%potential_temperature)%v)%data_3d(:,:,j) + domain%tend%th(:,:,j)*internal_dt * options%cu%tend_th_fraction
                if (options%cu%tend_qi_fraction > 0) domain%vars_3d(domain%var_indx(kVARS%ice_mass)%v)%data_3d(:,:,j)        = domain%vars_3d(domain%var_indx(kVARS%ice_mass)%v)%data_3d(:,:,j)        + domain%tend%qi(:,:,j)*internal_dt * options%cu%tend_qi_fraction
            endif
            ! if (options%physics%convection==kCU_KAINFR) then
            !     domain%qsnow(:,:,j) =domain%qsnow(:,:,j) + domain%tend%Qs(:,:,j)*internal_dt
            !     domain%qrain(:,:,j) =domain%qrain(:,:,j) + domain%tend%Qr(:,:,j)*internal_dt
            ! endif

            domain%vars_2d(domain%var_indx(kVARS%precipitation)%v)%data_2d(:,j)  = domain%vars_2d(domain%var_indx(kVARS%precipitation)%v)%data_2d(:,j) + RAINCV(:,j)
            domain%vars_2d(domain%var_indx(kVARS%convective_precipitation)%v)%data_2d(:,j) = domain%vars_2d(domain%var_indx(kVARS%convective_precipitation)%v)%data_2d(:,j) + RAINCV(:,j)

            ! if (options%physics%convection==kCU_TIEDTKE) then
            !     domain%u_cu(ids+1:ide,:,j) = 0.999*domain%u_cu(ids+1:ide,:,j) + (domain%tend%u(ids:ide-1,:,j)+domain%tend%u(ids+1:ide,:,j))/2 * internal_dt
            !     if (j>jds) then
            !         domain%v_cu(:,:,j) = 0.999*domain%v_cu(:,:,j) + (domain%tend%v(:,:,j)+domain%tend%v(:,:,j-1))/2 * internal_dt
            !     endif
            ! endif
        enddo
        ! $omp end do
        ! $omp end parallel

    endif

end subroutine convect

    !>------------------------------------------------------------
    !! Calculate the ZNU and ZNW variables
    !!
    !! @param domain    Model domain structure
    !!
    !!------------------------------------------------------------
    subroutine calc_znu(domain, znu)
        implicit none
        type(domain_t), intent(in) :: domain
        real, allocatable, intent(out) :: znu(:)

        integer :: i, xpt, ypt
        real    :: ptop
        integer :: kms, kme

        kms = domain%kms
        kme = domain%kme
        allocate(znu(kms:kme))
        ! one grid point into the domain gets a non-boundary point
        xpt = domain%ims + domain%grid%halo_size
        ypt = domain%jms + domain%grid%halo_size

        associate(p     => domain%vars_3d(domain%var_indx(kVARS%pressure)%v)%data_3d,                         &
                psfc  => domain%vars_2d(domain%var_indx(kVARS%surface_pressure)%v)%data_2d(xpt, ypt))

        ptop = p(xpt,kme,ypt) - (p(xpt,kme-1,ypt) - p(xpt,kme,ypt))/2.0 !NOT CORRECT
        ptop = max(ptop,1.0)

        do i=kms, kme
            znu(i) = (p(xpt,i,ypt) - ptop) / (psfc - ptop)
        enddo

        ! if (domain%var_indx(kVARS%znw)%v > 0) then
        !     do i = kms, kme
        !         if (i > kms) then
        !             domain%vars_1d(domain%var_indx(kVARS%znw)%v)%data_1d(i) = ((p(xpt,i,ypt) + p(xpt,i-1,ypt)) / 2 - ptop) / (psfc-ptop)
        !         else
        !             domain%vars_1d(domain%var_indx(kVARS%znw)%v)%data_1d(i) = 1
        !         endif
        !     enddo
        ! endif

        end associate

    end subroutine calc_znu

end module convection
