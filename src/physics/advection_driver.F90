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
    use adv_mpdata,                 only : mpdata_init, mpdata_advect3d, mpdata_compute_wind
    use adv_fluxcorr,               only : init_fluxcorr, set_sign_arrays, clear_flux_sign_arrays, compute_upwind_fluxes_async
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
        else if(options%physics%advection==kADV_MPDATA) then
            if (STD_OUT_PE .and. .not.context_change) write(*,*) "    MP-DATA"
            call mpdata_init(domain,options)
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
        if (options%physics%advection==kADV_MPDATA) then
            call adv_std_var_request(options)
        else
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

        real, allocatable :: U_m(:,:,:), V_m(:,:,:), W_m(:,:,:), denom(:,:,:), temp(:,:,:)

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
            if (options%time%RK3) then
                ! Allocate temporary array for RK3
                allocate(temp(ims:ime,kms:kme,jms:jme))
                !$acc enter data create(temp)
            endif
            call adv_wind_time%stop()

        else if(options%physics%advection==kADV_MPDATA) then
            call mpdata_compute_wind(domain,options,dt)
        else
            return
        endif

        do n = 1, size(domain%adv_vars)

            if (options%time%RK3) then
                if (options%physics%advection==kADV_STD) then
    
                    call RK3_adv(domain%vars_3d(domain%adv_vars(n)%v)%data_3d, temp, options, U_m, V_m, W_m, denom, &
                         domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d, flux_time, flux_corr_time, sum_time, adv_wind_time, n)
                    
                else if(options%physics%advection==kADV_MPDATA) then
                    ! Not yet implemented (is it compatable w/ RK3?)
                endif
            else
                if (options%physics%advection==kADV_STD) then
                    call adv_std_advect3d(domain%vars_3d(domain%adv_vars(n)%v)%data_3d,domain%vars_3d(domain%adv_vars(n)%v)%data_3d, &
                        U_m, V_m, W_m, denom, domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d,flux_time, flux_corr_time, sum_time)
                !else if(options%physics%advection==kADV_MPDATA) then                                    
                !    call mpdata_advect3d(var, rho, jaco, dz, options)
                endif
            endif
        enddo

        call adv_std_clean_wind_arrays(U_m, V_m, W_m, denom)
        if (options%adv%flux_corr==kFLUXCOR_MONO) call clear_flux_sign_arrays()
        if (options%physics%advection==kADV_STD .and. options%time%RK3) then
            !$acc exit data delete(temp)
            deallocate(temp)
        endif
    end subroutine advect

    subroutine RK3_adv(var, temp, options, U_m, V_m, W_m, denom, dz, flux_time, flux_corr_time, sum_time, adv_wind_time, q_id)
        implicit none
        real, allocatable, intent(inout) :: var(:,:,:), temp(:,:,:)
        type(options_t), intent(in) :: options
        real, allocatable, intent(in) :: U_m(:,:,:), V_m(:,:,:), W_m(:,:,:), denom(:,:,:), dz(:,:,:)
        type(timer_t), intent(inout) :: flux_time, flux_corr_time, sum_time, adv_wind_time
        integer, intent(in) :: q_id

        integer :: j, k, i

        !$acc kernels present(var, temp) async(q_id)
        temp(:,:,:) = var(:,:,:)
        !$acc end kernels

        !Initial advection-tendency calculations
        if (options%adv%flux_corr==kFLUXCOR_MONO) then
            !$acc wait(q_id) !temp is needed right away
            call compute_upwind_fluxes_async(temp,U_m, V_m, W_m, dz, denom, (q_id+100))
        endif

        ! wait on q_id will be handeled inside adv_std_advect3d
        call adv_std_advect3d(var, temp, U_m, V_m, W_m, denom, dz,flux_time, flux_corr_time, sum_time,t_factor_in=0.333,q_id_in=q_id)
        call adv_std_advect3d(var, temp, U_m, V_m, W_m, denom, dz,flux_time, flux_corr_time, sum_time,t_factor_in=0.5,q_id_in=q_id)
        !final advection call with tendency-fluxes
        call adv_std_advect3d(var, temp, U_m, V_m, W_m, denom, dz,flux_time, flux_corr_time, sum_time,flux_corr_in=options%adv%flux_corr,q_id_in=q_id)
        
    end subroutine RK3_adv

end module advection
