!> ----------------------------------------------------------------------------
!!  Driver to call different advection schemes
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!! ----------------------------------------------------------------------------
module advection
    use icar_constants
    use adv_std,                    only : adv_std_init, adv_std_var_request, adv_std_advect3d, adv_std_compute_wind, adv_std_clean_wind_arrays
    use adv_fluxcorr,               only : init_fluxcorr, set_sign_arrays, compute_upwind_fluxes_async
    ! use debug_module,               only : domain_fix
    use options_interface,          only: options_t
    use domain_interface,           only: domain_t
    use variable_dict_interface,  only : var_dict_t
    use variable_interface,       only : variable_t
    use timer_interface,          only : timer_t
    use data_structures,          only : index_type

    implicit none
    private
    integer :: ims, ime, kms, kme, jms, jme
    public :: advect, adv_init, adv_var_request
contains

    subroutine adv_init(domain,options, context_chng)
        implicit none
        type(domain_t), intent(inout) :: domain
        type(options_t),intent(in) :: options
        logical, optional, intent(in) :: context_chng

        logical :: context_change

        if (present(context_chng)) then
            context_change = context_chng
        else
            context_change = .false.
        endif

        if (STD_OUT_PE .and. .not.context_change) write(*,*) ""
        if (STD_OUT_PE .and. .not.context_change) write(*,*) "Initializing Advection"
        if (options%physics%advection==kADV_STD) then
            if (STD_OUT_PE .and. .not.context_change) write(*,*) "    Standard"
            call adv_std_init(domain,options)
            if (options%adv%flux_corr > 0) call init_fluxcorr(domain)
        endif
        
        ims = domain%ims; ime = domain%ime
        kms = domain%kms; kme = domain%kme
        jms = domain%jms; jme = domain%jme
        
    end subroutine adv_init

    subroutine adv_var_request(options)
        implicit none
        type(options_t),intent(inout) :: options

        !if (options%physics%advection==kADV_UPWIND) then
        !    call upwind_var_request(options)
        if (options%physics%advection==kADV_STD) then
            call adv_std_var_request(options)
        endif
    end subroutine
    
    subroutine advect(domain, options, dt,flux_time, flux_corr_time, sum_time, adv_wind_time)
        implicit none
        type(domain_t), intent(inout) :: domain
        type(options_t),intent(in)    :: options
        real,intent(in) :: dt
        type(timer_t), intent(inout) :: flux_time, flux_corr_time, sum_time, adv_wind_time
        type(variable_t) :: var_to_advect
        integer :: n

        real, allocatable :: U_m(:,:,:), V_m(:,:,:), W_m(:,:,:), denom(:,:,:)

        if (options%physics%advection==kADV_STD) then

            call adv_wind_time%start()
            call adv_std_compute_wind(domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d, &
                                      domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d, &
                                      domain%vars_3d(domain%var_indx(kVARS%w)%v)%data_3d, &
                                      domain%vars_3d(domain%var_indx(kVARS%density)%v)%data_3d, &
                                      domain%vars_3d(domain%var_indx(kVARS%jacobian)%v)%data_3d, &
                                      domain%vars_3d(domain%var_indx(kVARS%jacobian_u)%v)%data_3d, &
                                      domain%vars_3d(domain%var_indx(kVARS%jacobian_v)%v)%data_3d, &
                                      domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d, &
                                      domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d, &
                                      domain%dx,options,dt, U_m, V_m, W_m, denom)
            if (options%adv%flux_corr==kFLUXCOR_MONO) call set_sign_arrays(U_m,V_m,W_m)
            call adv_wind_time%stop()
        else
            return
        endif

        if (options%time%RK3 .and. options%physics%advection==kADV_STD) then
            call RK3_adv(domain, options, U_m, V_m, W_m, denom, flux_time, flux_corr_time, sum_time, adv_wind_time)
        elseif (options%physics%advection==kADV_STD) then
            do n = 1, size(domain%adv_vars)

                ! if (options%physics%advection==kADV_STD) then
    
                !     call RK3_adv(domain%vars_3d(domain%adv_vars(n)%v)%data_3d, temp, options, U_m, V_m, W_m, denom, &
                !          domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d, flux_time, flux_corr_time, sum_time, adv_wind_time, n)
                ! 
            ! else
                if (options%physics%advection==kADV_STD) then
                    call adv_std_advect3d(domain%vars_3d(domain%adv_vars(n)%v)%data_3d,domain%vars_3d(domain%adv_vars(n)%v)%data_3d, &
                        U_m, V_m, W_m, denom, domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d,flux_time, flux_corr_time, sum_time)
                endif
            enddo
        endif

        call adv_std_clean_wind_arrays(U_m, V_m, W_m, denom)
        ! if (options%adv%flux_corr==kFLUXCOR_MONO) call clear_flux_sign_arrays()
    end subroutine advect

    subroutine RK3_adv(domain, options, U_m, V_m, W_m, denom, flux_time, flux_corr_time, sum_time, adv_wind_time)
        implicit none
        type(domain_t), intent(inout) :: domain
        type(options_t), intent(in) :: options
        real, allocatable, intent(in) :: U_m(:,:,:), V_m(:,:,:), W_m(:,:,:), denom(:,:,:)
        type(timer_t), intent(inout) :: flux_time, flux_corr_time, sum_time, adv_wind_time

        integer :: RK3_step, n, n_adv, flux_corr, q_id
        real :: t_fac
        real, allocatable :: temp_all(:,:,:,:)

        n_adv = size(domain%adv_vars)

        ! Allocate 4D temporary array: one 3D slab per advected variable
        allocate(temp_all(ims:ime, kms:kme, jms:jme, n_adv))
        !$acc data create(temp_all)

        ! Save initial states for all advected variables
        do n = 1, n_adv
            associate( curr_var_data3d => domain%vars_3d(domain%adv_vars(n)%v)%data_3d)
            !$acc kernels present(curr_var_data3d, temp_all) async(n)
            temp_all(:,:,:,n) = curr_var_data3d(:,:,:)
            !$acc end kernels
            end associate
        enddo
        !$acc wait

        ! RK3 loop (outer) over variables (inner), with one batch exchange per step
        do RK3_step = 1, 3

            select case(RK3_step)
            case (1)
                t_fac = 1.0/3.0
                flux_corr = 0
            case (2)
                t_fac = 0.5
                flux_corr = 0
            case (3)
                t_fac = 1.0
                flux_corr = options%adv%flux_corr
            end select

            do n = 1, n_adv
                q_id = n

                if (flux_corr==kFLUXCOR_MONO) then
                    call flux_corr_time%start()
                    call compute_upwind_fluxes_async(temp_all(:,:,:,n), U_m, V_m, W_m, &
                        domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d, denom, (q_id+100))
                    call flux_corr_time%stop()
                endif


                call adv_std_advect3d(domain%vars_3d(domain%adv_vars(n)%v)%data_3d, &
                    temp_all(:,:,:,n), U_m, V_m, W_m, denom, &
                    domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d, flux_time, flux_corr_time, sum_time, &
                    t_factor_in=t_fac, q_id_in=q_id, flux_corr_in=flux_corr)
            enddo

            if (RK3_step < 3) then
                ! Single batch exchange for all advected variables after each RK3 step
                call domain%halo_3d_send()
                call domain%halo_3d_retrieve()
            endif
        enddo

        !$acc end data
        deallocate(temp_all)

    end subroutine RK3_adv

end module advection
