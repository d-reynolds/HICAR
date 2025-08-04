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
    use adv_fluxcorr,               only : init_fluxcorr, set_sign_arrays, clear_flux_sign_arrays
    ! use debug_module,               only : domain_fix
    use options_interface,          only: options_t
    use domain_interface,           only: domain_t
    use variable_dict_interface,  only : var_dict_t
    use variable_interface,       only : variable_t
    use timer_interface,          only : timer_t

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
    
    subroutine advect(domain, options, dt,flux_time, flux_up_time, flux_corr_time, sum_time, adv_wind_time)
        implicit none
        type(domain_t), intent(inout) :: domain
        type(options_t),intent(in)    :: options
        real,intent(in) :: dt
        type(timer_t), intent(inout) :: flux_time, flux_up_time, flux_corr_time, sum_time, adv_wind_time
        type(variable_t) :: var_to_advect
        integer :: j, k, i, n

        real, dimension(ims:ime,kms:kme,jms:jme) :: temp

        real, allocatable :: U_m(:,:,:), V_m(:,:,:), W_m(:,:,:), denom(:,:,:)

        !$acc data present(domain,kVARS) create(temp)

        if (options%physics%advection==kADV_STD) then

            call adv_wind_time%start()
            call adv_std_compute_wind(domain,options,dt, U_m, V_m, W_m, denom)
            if (options%adv%flux_corr==kFLUXCOR_MONO) call set_sign_arrays(U_m,V_m,W_m)
            call adv_wind_time%stop()

        else if(options%physics%advection==kADV_MPDATA) then
            call mpdata_compute_wind(domain,options,dt)
        endif

        do n = 1, size(domain%adv_vars)

            !$acc parallel loop gang vector collapse(3)
            do j = jms,jme
            do k = kms,kme
            do i = ims,ime
                temp(i,k,j) = domain%vars_3d(domain%adv_vars(n)%v)%data_3d(i,k,j)
            enddo
            enddo
            enddo

            if (options%time%RK3) then
                if (options%physics%advection==kADV_STD) then
    
                    !Initial advection-tendency calculations

                    call adv_std_advect3d(domain%vars_3d(domain%adv_vars(n)%v)%data_3d, temp, U_m, V_m, W_m, denom, domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d,flux_time, flux_up_time, flux_corr_time, sum_time,t_factor_in=0.333)
                    call adv_std_advect3d(domain%vars_3d(domain%adv_vars(n)%v)%data_3d, temp, U_m, V_m, W_m, denom, domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d,flux_time, flux_up_time, flux_corr_time, sum_time,t_factor_in=0.5)
                    !final advection call with tendency-fluxes
                    call adv_std_advect3d(domain%vars_3d(domain%adv_vars(n)%v)%data_3d, temp, U_m, V_m, W_m, denom, domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d,flux_time, flux_up_time, flux_corr_time, sum_time,flux_corr_in=options%adv%flux_corr)
                    
                else if(options%physics%advection==kADV_MPDATA) then
                    ! Not yet implemented (is it compatable w/ RK3?)
                endif
            else
                if (options%physics%advection==kADV_STD) then
                    call adv_std_advect3d(domain%vars_3d(domain%adv_vars(n)%v)%data_3d,temp, U_m, V_m, W_m, denom, domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d,flux_time, flux_up_time, flux_corr_time, sum_time)
                !else if(options%physics%advection==kADV_MPDATA) then                                    
                !    call mpdata_advect3d(var, rho, jaco, dz, options)
                endif
            endif
        enddo
        !$acc end data
        call adv_std_clean_wind_arrays(U_m, V_m, W_m, denom)
        if (options%adv%flux_corr==kFLUXCOR_MONO) call clear_flux_sign_arrays()

    end subroutine advect

end module advection
