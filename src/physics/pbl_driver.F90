!>----------------------------------------------------------
!! This module provides a wrapper to call various PBL models
!! It sets up variables specific to the physics package to be used including both
!!
!! The main entry point to the code is pbl(domain,options,dt)
!!
!! <pre>
!! Call tree graph :
!!  pbl_init->[ external initialization routines]
!!  pbl->[  external PBL routines]
!!  pbl_finalize
!!
!! High level routine descriptions / purpose
!!   pbl_init           - initializes physics package
!!   pbl                - sets up and calls main physics package
!!   pbl_finalize       - permits physics package cleanup (close files, deallocate memory)
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
module planetary_boundary_layer
    use domain_interface,   only : domain_t
    use options_interface,  only : options_t
    !use pbl_simple,    only : simple_pbl, finalize_simple_pbl, init_simple_pbl
    !use pbl_diagnostic, only : diagnostic_pbl, finalize_diagnostic_pbl, init_diagnostic_pbl
    !use module_bl_ysu, only : ysuinit, ysu
    use module_bl_ysu, only : ysuinit, ysu
    use module_bl_ysu_gpu, only : ysuinit_gpu, ysu_gpu
    use mod_wrf_constants, only : EOMEG, XLV, XLS, r_v, R_d, KARMAN, gravity, EP_1, EP_2, cp, rcp, rovg
    use icar_constants !, only : karman,stefan_boltzmann
    use ieee_arithmetic ! for debugging
    use array_utilities, only : array_offset_x_3d, array_offset_y_3d


    implicit none
    real,allocatable, dimension(:,:)    ::  windspd, regime, sensible_heat_tmp, qfx_tmp
    ! integer, allocatable, dimension(:,:) :: kpbl2d
    real, allocatable, dimension(:,:,:) :: RTHRATEN
    real, allocatable, dimension(:,:,:) :: trid_a, trid_b, trid_c, trid_rhs  ! work arrays for scalar PBL diffusion

    private
    public :: pbl_var_request, pbl_init, pbl, pbl_finalize, pbl_apply_tend

    integer :: ids, ide, jds, jde, kds, kde,  &
               ims, ime, jms, jme, kms, kme,  &
               its, ite, jts, jte, kts, kte, j, k, i

    logical :: allowed_to_read, flag_qi

contains

    subroutine pbl_var_request(options)
        implicit none
        type(options_t),intent(inout) :: options

        if (options%physics%boundarylayer==kPBL_YSU) then

            call options%alloc_vars( &
                         [kVARS%water_vapor, kVARS%potential_temperature, kVARS%temperature,                &
                         kVARS%exner, kVARS%dz_interface, kVARS%density, kVARS%pressure_interface,          &
                         kVARS%skin_temperature, kVARS%terrain, kVARS%ground_surf_temperature,              &
                         kVARS%sensible_heat, kVARS%latent_heat, kVARS%u_10m, kVARS%v_10m,                  &
                         kVARS%humidity_2m, kVARS%surface_pressure, kVARS%ground_heat_flux,                 &
                         kVARS%roughness_z0, kVARS%ustar, kVARS%ice_mass,                                  &
                         kVARS%tend_th_pbl, kVARS%tend_qc_pbl, kVARS%tend_qi_pbl,  kVARS%temperature_2m,    &
                         kVARS%tend_u, kVARS%tend_v, kVARS%tend_qv_pbl, kVARS%pressure, kVARS%kpbl,         &
                         kVARS%fm, kVARS%fh, kVARS%QFX, kVARS%br,                                          &
                         kVARS%land_mask, kVARS%cloud_water_mass, kVARS%coeff_heat_exchange_3d, kVARS%coeff_momentum_exchange_3d, kVARS%hpbl ]) !kVARS%tend_qv_adv,kVARS%tend_qv, kVARS%tend_qs, kVARS%tend_qr,, kVARS%u_mass, kVARS%v_mass,
!           kVARS%coeff_momentum_drag, ??

             call options%restart_vars( &
                        [kVARS%water_vapor, kVARS%potential_temperature, kVARS%temperature,                &
                        kVARS%exner, kVARS%dz_interface, kVARS%density, kVARS%pressure_interface,          &
                        kVARS%skin_temperature, kVARS%terrain, kVARS%ground_surf_temperature,              &
                        kVARS%sensible_heat, kVARS%latent_heat, kVARS%u_10m, kVARS%v_10m,                  &
                        kVARS%humidity_2m, kVARS%surface_pressure, kVARS%ground_heat_flux,                 &
                        kVARS%roughness_z0, kVARS%ice_mass, kVARS%QFX,       &
                        kVARS%temperature_2m,    &
                        kVARS%pressure,         &
                        kVARS%cloud_water_mass,kVARS%coeff_heat_exchange_3d, kVARS%coeff_momentum_exchange_3d, kVARS%hpbl  ]) !kVARS%u_mass, kVARS%v_mass,
        endif
    end subroutine pbl_var_request


    subroutine pbl_init(domain,options,context_chng)
        implicit none
        type(domain_t),     intent(inout)   :: domain
        type(options_t),    intent(in)      :: options
        logical, optional, intent(in)       :: context_chng

        logical :: context_change, restart

        if (present(context_chng)) then
            context_change = context_chng
        else
            context_change = .False.
        endif
        ids = domain%ids ; ide = domain%ide ; jds = domain%jds ; jde = domain%jde ; kds = domain%kds ; kde = domain%kde
        ims = domain%ims ; ime = domain%ime ; jms = domain%jms ; jme = domain%jme ; kms = domain%kms ; kme = domain%kme
        its = domain%its ; ite = domain%ite ; jts = domain%jts ; jte = domain%jte ; kts = domain%kts ; kte = domain%kte

        allowed_to_read = .True.
        restart = .False.
        flag_qi = .true.

        if (STD_OUT_PE .and. .not.context_change) write(*,*) "Initializing PBL Scheme"

        !if (options%physics%boundarylayer==kPBL_SIMPLE) then
        !    if (STD_OUT_PE .and. .not.context_change) write(*,*) "    Simple PBL"
        !    call init_simple_pbl(domain, options)
        !else if (options%physics%boundarylayer==kPBL_DIAGNOSTIC) then
        !    if (STD_OUT_PE .and. .not.context_change) write(*,*) "    Diagnostic PBL"
        !    call init_diagnostic_pbl(domain, options)
        if (options%physics%boundarylayer==kPBL_YSU) then

            if (STD_OUT_PE .and. .not.context_change) write(*,*) "    YSU PBL"

            ! allocate local vars YSU:
            if (allocated(windspd)) then
                !$acc exit data delete(windspd)
                deallocate(windspd)
            endif
            allocate(windspd(ims:ime, jms:jme))
            ! allocate(hpbl(ims:ime, jms:jme))  ! this should go to domain object for convective modules!!
            !allocate(u10d(ims:ime, jms:jme))
            !allocate(v10d(ims:ime, jms:jme))
            ! allocate(kpbl2d(ims:ime, jms:jme)) ! domain%kpbl now
            ! allocate(CHS(ims:ime,jms:jme))
            ! CHS = 0.01
            !allocate(xland_real(ims:ime,jms:jme))
            !xland_real=real(domain%land_mask)
            if (allocated(regime)) then
                !$acc exit data delete(regime)
                deallocate(regime)
            endif
            allocate(regime(ims:ime,jms:jme))
            if (allocated(RTHRATEN)) then
                !$acc exit data delete(RTHRATEN)
                deallocate(RTHRATEN)
            endif
            allocate(RTHRATEN(ims:ime, kms:kme, jms:jme)) !initialize radiative heating tendencies and set to 0 in case user turns on ysu radiative heating w/o radiations scheme
            RTHRATEN = 0.0

            if (allocated(sensible_heat_tmp)) then
                !$acc exit data delete(sensible_heat_tmp)
                deallocate(sensible_heat_tmp)
            endif
            allocate(sensible_heat_tmp(ims:ime, jms:jme))

            if (allocated(qfx_tmp)) then
                !$acc exit data delete(qfx_tmp)
                deallocate(qfx_tmp)
            endif
            allocate(qfx_tmp(ims:ime, jms:jme))

            if (allocated(trid_a)) then
                !$acc exit data delete(trid_a, trid_b, trid_c, trid_rhs)
                deallocate(trid_a, trid_b, trid_c, trid_rhs)
            endif
            allocate(trid_a(ims:ime, kms:kme, jms:jme))
            allocate(trid_b(ims:ime, kms:kme, jms:jme))
            allocate(trid_c(ims:ime, kms:kme, jms:jme))
            allocate(trid_rhs(ims:ime, kms:kme, jms:jme))

            !$acc enter data copyin(windspd,regime,RTHRATEN,qfx_tmp,sensible_heat_tmp) create(trid_a,trid_b,trid_c,trid_rhs)
            
            ! initialize tendencies (this is done in ysu init but only for tiles, not mem (ie its vs ims))
            ! BK: check if this actually matters ???
            if(.not.context_change)then
                do j = jms,jme
                do k = kms,kme
                do i = ims,ime
                    domain%tend%u(i,k,j) = 0.
                    domain%tend%v(i,k,j) = 0.
                    domain%tend%th_pbl(i,k,j) = 0.
                    domain%tend%qv_pbl(i,k,j) = 0.
                    domain%tend%qc_pbl(i,k,j) = 0.
                    domain%tend%qi_pbl(i,k,j) = 0.
                enddo
                enddo
                enddo
            endif


            call ysuinit_gpu(rublten=domain%tend%u                  &
                        ,rvblten=domain%tend%v                  &
                        ,rthblten=domain%tend%th_pbl            &
                        ,rqvblten=domain%tend%qv_pbl            &
                        ,rqcblten=domain%tend%qc_pbl            &
                        ,rqiblten=domain%tend%qi_pbl            &
                        ,p_qi=1                                 &
                        ,p_first_scalar=1                       &
                        ,restart=context_change                        &
                        ,allowed_to_read= allowed_to_read      &
                        ,ids=ids, ide=ide, jds=jds, jde=jde     &
                        ,kds=kds, kde=kde, ims=ims, ime=ime     &
                        ,jms=jms, jme=jme, kms=kms, kme=kme     &
                        ,its=its, ite=ite, jts=jts, jte=jte     &
                        ,kts=kts, kte=kte-1)
        endif
    end subroutine pbl_init

    subroutine pbl(domain, options, dt_in)
        implicit none
        type(domain_t),  intent(inout)  :: domain
        type(options_t), intent(in)     :: options
        real,            intent(in)     :: dt_in  !  =real(dt%seconds())

        ! if (options%physics%boundarylayer==kPBL_SIMPLE) then
        !     call simple_pbl(domain% potential_temperature %data_3d,     &
        !                     domain% water_vapor           %data_3d,     &
        !                     domain% cloud_water_mass      %data_3d,     &
        !                     domain% cloud_ice_mass        %data_3d,     &
        !                     domain% rain_mass             %data_3d,     &
        !                     domain% snow_mass             %data_3d,     &
        !                     domain% u_mass                %data_3d,     &
        !                     domain% v_mass                %data_3d,     &
        !                     domain% exner                 %data_3d,     &
        !                     domain% density               %data_3d,     &
        !                     domain% z                     %data_3d,     &
        !                     domain% dz_mass               %data_3d,     &
        !                     domain% terrain               %data_2d,     &
        !                     domain% land_mask,                          &
        !                     its, ite, jts, jte, kts, kte,               &
        !                     dt_in)
        !                     ! domain% qv_pbl_tendency     %data_3d)
        ! endif
        ! if (options%physics%boundarylayer==kPBL_DIAGNOSTIC) then
        !     call diagnostic_pbl(domain,dt_in)
        !                     ! domain% qv_pbl_tendency     %data_3d)
        ! endif
        if (options%physics%boundarylayer==kPBL_YSU) then

            ! Reset tendencies before the next pbl call. (not sure if necessary)

            associate(tend_u => domain%tend%u, tend_v => domain%tend%v, tend_th_pbl => domain%tend%th_pbl, &
                      tend_qv_pbl => domain%tend%qv_pbl, tend_qc_pbl => domain%tend%qc_pbl, tend_qi_pbl => domain%tend%qi_pbl, &
                      sensible_heat => domain%vars_2d(domain%var_indx(kVARS%sensible_heat)%v)%data_2d, &
                      qfx => domain%vars_2d(domain%var_indx(kVARS%qfx)%v)%data_2d, &
                      u_10m => domain%vars_2d(domain%var_indx(kVARS%u_10m)%v)%data_2d, &
                      v_10m => domain%vars_2d(domain%var_indx(kVARS%v_10m)%v)%data_2d)
            !$acc parallel present(tend_u, tend_v, tend_th_pbl, tend_qv_pbl, tend_qc_pbl, tend_qi_pbl, &
            !$acc &               u_10m, v_10m, windspd,regime)
            !$acc loop gang vector collapse(3)
            do j = jms,jme
            do k = kms,kme
            do i = ims,ime
                tend_u(i,k,j)       = 0
                tend_v(i,k,j)       = 0
                tend_th_pbl(i,k,j)  = 0
                tend_qv_pbl(i,k,j)  = 0
                tend_qc_pbl(i,k,j)  = 0
                tend_qi_pbl(i,k,j)  = 0
            enddo
            enddo
            enddo

            !$acc loop gang vector collapse(2)
            do j = jms,jme
            do i = ims,ime
                windspd(i,j) = sqrt(u_10m(i,j)**2 + v_10m(i,j)**2) ! as it is done in lsm_driver.
                if (windspd(i,j)==0) windspd(i,j)=1e-5
            enddo
            enddo
            !$acc end parallel

            ! Radiative heating tendencies are only allocated when the
            ! radiation scheme provides them (RRTMG/RRTMGP); referencing
            ! vars_3d(-1) when they are absent is an out-of-bounds access.
            if (domain%var_indx(kVARS%tend_th_lwrad)%v > 0 .and. &
                domain%var_indx(kVARS%tend_th_swrad)%v > 0) then
                associate(tend_th_lwrad => domain%vars_3d(domain%var_indx(kVARS%tend_th_lwrad)%v)%data_3d, &
                          tend_th_swrad => domain%vars_3d(domain%var_indx(kVARS%tend_th_swrad)%v)%data_3d)
                !$acc parallel loop gang vector collapse(3) present(tend_th_lwrad, tend_th_swrad, RTHRATEN)
                do j = jms,jme
                do k = kms,kme
                do i = ims,ime
                    RTHRATEN(i,k,j) = tend_th_lwrad(i,k,j) + tend_th_swrad(i,k,j)
                enddo
                enddo
                enddo
                end associate
            endif


            ! Add blowing snow sublimation feedback to the sensible and latent heat
            ! flux fields. The PBL scheme reads these and
            ! distributes the fluxes vertically through the boundary layer.
            ! bs_subl is in W/m^2 (positive = sublimation occurring, energy consumed)
            if (options%sm%suspension_layer == 1 .and. options%sm%bs_atm_feedback .and. options%physics%snowmodel  > 0) then
                associate(bs_subl => domain%vars_2d(domain%var_indx(kVARS%bs_sublimation_flux)%v)%data_2d)
                !$acc parallel loop gang vector collapse(2) default(present)
                do j = jms, jme
                    do i = ims ,ime
                        ! Sublimation cools air → negative sensible heat contribution
                        sensible_heat_tmp(i,j) = sensible_heat(i,j) - bs_subl(i,j)  ! W/m^2, positive = sublimation

                        ! Sublimation adds moisture → positive latent heat contribution
                        qfx_tmp(i,j) = qfx(i,j) + bs_subl(i,j) / XLS
                    enddo
                enddo
                end associate
            else
                !$acc parallel loop gang vector collapse(2) default(present)
                do j = jms, jme
                    do i = ims, ime
                        sensible_heat_tmp(i,j) = sensible_heat(i,j)
                        qfx_tmp(i,j) = qfx(i,j)
                    enddo
                enddo
            endif
            end associate

            call ysu_gpu(u3d=domain%vars_3d(domain%var_indx(kVARS%u_mass)%v)%data_3d                           & !-- u3d         3d u-velocity interpolated to theta points (m/s)
                    ,v3d=domain%vars_3d(domain%var_indx(kVARS%v_mass)%v)%data_3d                           & !-- v3d         3d v-velocity interpolated to theta points (m/s)
                    ,th3d=domain%vars_3d(domain%var_indx(kVARS%potential_temperature)%v)%data_3d           &
                    ,t3d=domain%vars_3d(domain%var_indx(kVARS%temperature)%v)%data_3d                      &
                    ,qv3d=domain%vars_3d(domain%var_indx(kVARS%water_vapor)%v)%data_3d                     &
                    ,qc3d=domain%vars_3d(domain%var_indx(kVARS%cloud_water_mass)%v)%data_3d                & !-- qc3d        cloud water mixing ratio (kg/kg)
                    ,qi3d=domain%vars_3d(domain%var_indx(kVARS%ice_mass)%v)%data_3d                  & !-- qi3d        cloud ice mixing ratio (kg/kg)
                    ,p3d=domain%vars_3d(domain%var_indx(kVARS%pressure)%v)%data_3d                         & !-- p3d         3d pressure (pa)
                    ,p3di=domain%vars_3d(domain%var_indx(kVARS%pressure_interface)%v)%data_3d              & !-- p3di        3d pressure (pa) at interface level
                    ,pi3d=domain%vars_3d(domain%var_indx(kVARS%exner)%v)%data_3d                           & !-- pi3d        3d exner function (dimensionless)
                    ,rublten=domain%tend%u                  & ! i/o (mass grid; applied to staggered u under kRANS_WINDS)
                    ,rvblten=domain%tend%v                  & ! i/o (mass grid; applied to staggered v under kRANS_WINDS)
                    ,rthblten=domain%tend%th_pbl            & ! i/o
                    ,rqvblten=domain%tend%qv_pbl            & ! i/o
                    ,rqcblten=domain%tend%qc_pbl            & ! i/o
                    ,flag_qi=.True.                         &
                    ,cp=cp                                  &
                    ,g=gravity                              &
                    ,rovcp=rcp                            & ! rovcp = Rd/cp
                    ,rd=R_d                                 &  ! J/(kg K) specific gas constant for dry air
                    ,rovg=rovg                              &
                    ,dz8w=domain%vars_3d(domain%var_indx(kVARS%dz_interface)%v)%data_3d       & !-- dz8w        dz between full levels (m)
                    ,z8w=domain%vars_3d(domain%var_indx(kVARS%z_interface)%v)%data_3d               & !-- z8w         height of full levels (m)
                    ,xlv=XLV                    & !-- xlv         latent heat of vaporization (j/kg)
                    ,rv=r_v                                  &  ! J/(kg K) specific gas constant for wet/moist air
                    ,psfc=domain%vars_2d(domain%var_indx(kVARS%surface_pressure)%v)%data_2d   &
                    ,znt=domain%vars_2d(domain%var_indx(kVARS%roughness_z0)%v)%data_2d       &  ! i/o -- znt		roughness length (m) (input only)
                    ,ust=domain%vars_2d(domain%var_indx(kVARS%ustar)%v)%data_2d                       & ! i/o -- ust		u* in similarity theory (m/s)
                    ,hpbl=domain%vars_2d(domain%var_indx(kVARS%hpbl)%v)%data_2d               & ! i/o -- hpbl	pbl height (m) - intent(inout)
                    ,psim=domain%vars_2d(domain%var_indx(kVARS%psim)%v)%data_2d               & !-- psim        similarity stability function for momentum - intent(in)
                    ,psih=domain%vars_2d(domain%var_indx(kVARS%psih)%v)%data_2d               & !-- psih        similarity stability function for heat- intent(in)
                    ,xland=domain%vars_2d(domain%var_indx(kVARS%land_mask)%v)%data_2di                               &
                    ,hfx=sensible_heat_tmp                     & !  HFX  - net upward heat flux at the surface (W/m^2)
                    ,qfx=qfx_tmp           & !  QFX  - net upward moisture flux at the surface (kg/m^2/s)
                    !,UOCE=uoce,VOCE=voce                                  & !ocean currents -- not currently used
                    !,CTOPO=ctopo,CTOPO2=ctopo2                            & !optional, only applied to momentum tendencies, not currently used
                    ! ,tke_pbl=domain%vars_2d(domain%var_indx(kVARS%tke_pbl)%v)%data_2d               &
                    ,YSU_TOPDOWN_PBLMIX=options%pbl%ysu_topdown_pblmix                &
                    ,wspd=windspd                           & ! i/o -- wspd        wind speed at lowest model level (m/s)
                    ,br=domain%vars_2d(domain%var_indx(kVARS%br)%v)%data_2d                   & !-- br          bulk richardson number in surface layer
                    ,dt=dt_in                               & !-- dt		time step (s)
                    ,kpbl2d=domain%vars_2d(domain%var_indx(kVARS%kpbl)%v)%data_2di                          & ! o --     ?? k layer of pbl top??
                    ,ep1=EP_1                                & !-- ep1         constant for virtual temperature (r_v/r_d - 1) (dimensionless)
                    ,ep2=EP_2                                & !-- ep2         constant for specific humidity calculation
                    ,karman=karman                          & !-- karman      von karman constant
                    ,RTHRATEN=RTHRATEN                                    &
!                    ,WSTAR=wstar,DELTA=delta                              &  !Output variables of YSU which we currently dont use
                    ,exch_h=domain%vars_3d(domain%var_indx(kVARS%coeff_heat_exchange_3d)%v)%data_3d  & ! i/o -- exch_h ! exchange coefficient for heat, K m/s , but 3d??
                    ! ,exch_m=domain%vars_3d(domain%var_indx(kVARS%coeff_momentum_exchange_3d)%v)%data_3d  & ! i/o -- exch_h ! exchange coefficient for heat, K m/s , but 3d??
                    ,u10=domain%vars_2d(domain%var_indx(kVARS%u_10m)%v)%data_2d               &
                    ,v10=domain%vars_2d(domain%var_indx(kVARS%v_10m)%v)%data_2d               &
                    ,ids=ids, ide=ide, jds=jds, jde=jde     &
                    ,kds=kds, kde=kde, ims=ims, ime=ime     &
                    ,jms=jms, jme=jme, kms=kms, kme=kme     &
                    ,its=its, ite=ite, jts=jts, jte=jte     &
                    ,kts=kts, kte=kte-1                     &
                !optional
                    ,rho=domain%vars_3d(domain%var_indx(kVARS%density)%v)%data_3d                           & !-- rho        3d density (kg/m^3)
                    ,rqiblten=domain%tend%qi_pbl            & ! i/o
                    ,regime=regime                          )!  i/o -- regime	flag indicating pbl regime (stable, unstable, etc.) - not used?


                    ! if(STD_OUT_PE .and. options%general%debug) write(*,*) "  pbl height/lev is:", maxval(domain%vars_2d(domain%var_indx(kVARS%hpbl)%v)%data_2d ),"m/", maxval(domain%kpbl)  ! uncomment if you want to see the pbl height.


        endif ! End YSU call

    end subroutine pbl

    subroutine pbl_apply_tend(domain,options,dt)
        implicit none
        type(domain_t),  intent(inout)  :: domain
        type(options_t), intent(in)     :: options
        real,            intent(in)     :: dt

        if (options%physics%boundarylayer==kPBL_YSU) then
            !> ------------  add tendency terms  ------------
            !
            ! Here the tendency terms that were calculated by the ysu routine are added to the domain-wide fields.
            ! Momentum tendencies (u, v) are only applied when running the prognostic RANS solver — the
            ! diagnostic mass-conserving solver overwrites winds every forcing step and would simply discard
            ! a PBL momentum tendency, so applying it would be wasted work (and the surface-layer drag would
            ! then be double-counted on top of the diagnostic balance).
            !
            !> -----------------------------------------------

            if (options%physics%windtype == kRANS_WINDS) then
                ! Apply mass-grid momentum tendencies to the staggered u-/v-grids.
                ! Interpolation is inlined on the device — avoids the host-side reallocate
                ! that array_offset_x_3d would do every call.
                !
                ! Interior u-face at i uses the average of mass-cell tendencies at i-1 and i.
                ! The outermost faces (i = ims at global west boundary, ime+1 at east) are
                ! left untouched: they are owned by the lateral-BC machinery (forcing +
                ! Davies relaxation) and would otherwise out-of-bound `tend_u_mass(i-1)`
                ! at i = its == ims on a global-west boundary rank.
                associate(u_data      => domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d, &
                          v_data      => domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d, &
                          tend_u_mass => domain%tend%u, &
                          tend_v_mass => domain%tend%v)
                !$acc parallel loop gang vector collapse(3) default(present)
                do j = jts, jte
                    do k = kts, kte
                        do i = max(its, ims+1), ite
                            u_data(i,k,j) = u_data(i,k,j) + &
                                0.5 * (tend_u_mass(i-1,k,j) + tend_u_mass(i,k,j)) * dt
                        end do
                    end do
                end do
                !$acc parallel loop gang vector collapse(3) default(present)
                do j = max(jts, jms+1), jte
                    do k = kts, kte
                        do i = its, ite
                            v_data(i,k,j) = v_data(i,k,j) + &
                                0.5 * (tend_v_mass(i,k,j-1) + tend_v_mass(i,k,j)) * dt
                        end do
                    end do
                end do
                end associate
            endif

            ! add mass grid tendencies

            associate(qv => domain%vars_3d(domain%var_indx(kVARS%water_vapor)%v)%data_3d, &
                        qc => domain%vars_3d(domain%var_indx(kVARS%cloud_water_mass)%v)%data_3d, &
                        th => domain%vars_3d(domain%var_indx(kVARS%potential_temperature)%v)%data_3d, &
                        qi => domain%vars_3d(domain%var_indx(kVARS%ice_mass)%v)%data_3d, &
                        qv_tend => domain%tend%qv_pbl, &
                        qc_tend => domain%tend%qc_pbl, &
                        th_tend => domain%tend%th_pbl, &
                        qi_tend => domain%tend%qi_pbl)

            !$acc parallel loop gang vector collapse(3) present(qv,qc,th,qi,qv_tend,qc_tend,th_tend,qi_tend)
            do j = jts, jte
                do k = kts, kte
                    do i = its, ite
                            qv(i,k,j) = qv(i,k,j) + qv_tend(i,k,j) * dt
                            qc(i,k,j) = qc(i,k,j) + qc_tend(i,k,j) * dt
                            th(i,k,j) = th(i,k,j) + th_tend(i,k,j) * dt
                            qi(i,k,j) = qi(i,k,j) + qi_tend(i,k,j) * dt
                    end do
                end do
            end do
            end associate

            ! Implicit vertical diffusion of hydrometeor number concentrations
            ! using the eddy diffusivity (exch_h) computed by YSU.
            ! Follows WRF's scalar_pblmix approach with zero-flux BCs.
            associate(exch_h => domain%vars_3d(domain%var_indx(kVARS%coeff_heat_exchange_3d)%v)%data_3d, &
                      rho    => domain%vars_3d(domain%var_indx(kVARS%density)%v)%data_3d,                &
                      dz     => domain%vars_3d(domain%var_indx(kVARS%dz_interface)%v)%data_3d)

            if (domain%var_indx(kVARS%snow_mass)%v > 0) &
                call pbl_scalar_diff(domain%vars_3d(domain%var_indx(kVARS%snow_mass)%v)%data_3d, exch_h, rho, dz, dt)
            if (domain%var_indx(kVARS%ice_number)%v > 0) &
                call pbl_scalar_diff(domain%vars_3d(domain%var_indx(kVARS%ice_number)%v)%data_3d, exch_h, rho, dz, dt)
            if (domain%var_indx(kVARS%rain_number)%v > 0) &
                call pbl_scalar_diff(domain%vars_3d(domain%var_indx(kVARS%rain_number)%v)%data_3d, exch_h, rho, dz, dt)
            if (domain%var_indx(kVARS%snow_number)%v > 0) &
                call pbl_scalar_diff(domain%vars_3d(domain%var_indx(kVARS%snow_number)%v)%data_3d, exch_h, rho, dz, dt)
            if (domain%var_indx(kVARS%graupel_number)%v > 0) &
                call pbl_scalar_diff(domain%vars_3d(domain%var_indx(kVARS%graupel_number)%v)%data_3d, exch_h, rho, dz, dt)

            end associate

        endif
    end subroutine pbl_apply_tend


    subroutine pbl_scalar_diff(phi, exch_h, rho, dz, dt_in)
        !---------------------------------------------------------------
        ! Apply implicit vertical diffusion to a scalar field using
        ! the eddy diffusivity (exch_h) output by the PBL scheme.
        ! Follows WRF's scalar_pblmix approach (Thomas algorithm).
        ! Zero-flux boundary conditions at surface and model top.
        !
        ! Uses module-level work arrays trid_a, trid_b, trid_c, trid_rhs
        ! and module-level index variables.
        !---------------------------------------------------------------
        implicit none
        real, dimension(ims:ime, kms:kme, jms:jme), intent(inout) :: phi
        real, dimension(ims:ime, kms:kme, jms:jme), intent(in)    :: exch_h, rho, dz
        real, intent(in) :: dt_in

        real    :: rho_int, dz_half, cddz_below, cddz_above, fk
        integer :: i2, j2, k2, kte_pbl

        kte_pbl = kte - 1   ! YSU operates on kts:kte-1

        ! Step 1: Build tridiagonal coefficients
        ! exch_h(i,k,j) = diffusivity at interface between levels k-1 and k [m^2/s]
        !$acc parallel present(trid_a, trid_b, trid_c, trid_rhs, phi, exch_h, rho, dz)
        !$acc loop gang vector collapse(3) private(rho_int, dz_half, cddz_below, cddz_above)
        do j2 = jts, jte
        do k2 = kts, kte_pbl
        do i2 = its, ite
            ! Sub-diagonal: coupling with level below
            if (k2 > kts) then
                rho_int    = 0.5*(rho(i2,k2-1,j2) + rho(i2,k2,j2))
                dz_half    = 0.5*(dz(i2,k2-1,j2)  + dz(i2,k2,j2))
                cddz_below = rho_int * exch_h(i2,k2,j2) / dz_half
                trid_a(i2,k2,j2) = -dt_in * cddz_below / (rho(i2,k2,j2) * dz(i2,k2,j2))
            else
                trid_a(i2,k2,j2) = 0.0   ! zero-flux surface BC
            endif

            ! Super-diagonal: coupling with level above
            if (k2 < kte_pbl) then
                rho_int    = 0.5*(rho(i2,k2,j2) + rho(i2,k2+1,j2))
                dz_half    = 0.5*(dz(i2,k2,j2)  + dz(i2,k2+1,j2))
                cddz_above = rho_int * exch_h(i2,k2+1,j2) / dz_half
                trid_c(i2,k2,j2) = -dt_in * cddz_above / (rho(i2,k2,j2) * dz(i2,k2,j2))
            else
                trid_c(i2,k2,j2) = 0.0   ! zero-flux top BC
            endif

            ! Main diagonal (diagonally dominant)
            trid_b(i2,k2,j2)   = 1.0 - trid_a(i2,k2,j2) - trid_c(i2,k2,j2)
            ! RHS = current scalar value
            trid_rhs(i2,k2,j2) = phi(i2,k2,j2)
        enddo
        enddo
        enddo
        !$acc end parallel

        ! Step 2: Thomas algorithm forward elimination
        ! First level
        !$acc parallel present(trid_a, trid_b, trid_c, trid_rhs)
        !$acc loop gang vector collapse(2) private(fk)
        do j2 = jts, jte
        do i2 = its, ite
            fk = 1.0 / trid_b(i2,kts,j2)
            trid_c(i2,kts,j2)   = fk * trid_c(i2,kts,j2)
            trid_rhs(i2,kts,j2) = fk * trid_rhs(i2,kts,j2)
        enddo
        enddo
        !$acc end parallel

        ! Interior + top levels
        !$acc parallel present(trid_a, trid_b, trid_c, trid_rhs)
        !$acc loop gang vector collapse(2) private(fk)
        do j2 = jts, jte
        do i2 = its, ite
        !$acc loop seq
            do k2 = kts+1, kte_pbl
                fk = 1.0 / (trid_b(i2,k2,j2) - trid_a(i2,k2,j2) * trid_c(i2,k2-1,j2))
                trid_c(i2,k2,j2)   = fk * trid_c(i2,k2,j2)
                trid_rhs(i2,k2,j2) = fk * (trid_rhs(i2,k2,j2) - trid_a(i2,k2,j2) * trid_rhs(i2,k2-1,j2))
            enddo
        enddo
        enddo
        !$acc end parallel

        ! Step 3: Back substitution — write result into phi, clamp non-negative
        !$acc parallel present(phi, trid_c, trid_rhs)
        !$acc loop gang vector collapse(2)
        do j2 = jts, jte
        do i2 = its, ite
            phi(i2,kte_pbl,j2) = max(0.0, trid_rhs(i2,kte_pbl,j2))
        !$acc loop seq
            do k2 = kte_pbl-1, kts, -1
                phi(i2,k2,j2) = max(0.0, trid_rhs(i2,k2,j2) - trid_c(i2,k2,j2) * phi(i2,k2+1,j2))
            enddo
        enddo
        enddo
        !$acc end parallel

    end subroutine pbl_scalar_diff


    subroutine pbl_finalize(options)
        implicit none
        type(options_t), intent(in) :: options

        !if (options%physics%boundarylayer==kPBL_SIMPLE) then
        !    call finalize_simple_pbl()
        !else if (options%physics%boundarylayer==kPBL_DIAGNOSTIC) then
        !    call finalize_diagnostic_pbl()
        !endif


        !$acc exit data delete(windspd,regime,RTHRATEN,trid_a,trid_b,trid_c,trid_rhs)
        if (allocated(trid_a)) deallocate(trid_a, trid_b, trid_c, trid_rhs)

    end subroutine pbl_finalize
end module planetary_boundary_layer
