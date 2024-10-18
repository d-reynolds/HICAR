!> ----------------------------------------------------------------------------
!!  Driver to call different advection schemes
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!! ----------------------------------------------------------------------------
module advection
    use data_structures
    use icar_constants
    use adv_std,                    only : adv_std_init, adv_std_var_request, adv_std_advect3d, adv_fluxcorr_advect3d, adv_std_compute_wind
    use adv_mpdata,                 only : mpdata_init, mpdata_advect3d, mpdata_compute_wind
    use adv_fluxcorr,               only : init_fluxcorr
    ! use debug_module,               only : domain_fix
    use options_interface,          only: options_t
    use domain_interface,           only: domain_t
    use variable_dict_interface,  only : var_dict_t
    use variable_interface,       only : variable_t

    implicit none
    private
    integer :: ims, ime, kms, kme, jms, jme
    public :: advect, adv_init, adv_var_request
contains

    subroutine adv_init(domain,options)
        implicit none
        type(domain_t), intent(inout) :: domain
        type(options_t),intent(in) :: options

        if (STD_OUT_PE) write(*,*) ""
        if (STD_OUT_PE) write(*,*) "Initializing Advection"
        if (options%physics%advection==kADV_STD) then
            if (STD_OUT_PE) write(*,*) "    Standard"
            call adv_std_init(domain,options)
        else if(options%physics%advection==kADV_MPDATA) then
            if (STD_OUT_PE) write(*,*) "    MP-DATA"
            call mpdata_init(domain,options)
        endif
        
        ims = domain%ims; ime = domain%ime
        kms = domain%kms; kme = domain%kme
        jms = domain%jms; jme = domain%jme
        
        if (options%adv%flux_corr > 0) call init_fluxcorr(domain)
        !Allocate storage variable for temp-quantities

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
    
    subroutine advect(domain, options, dt)
        implicit none
        type(domain_t), intent(inout) :: domain
        type(options_t),intent(in)    :: options
        real,intent(in) :: dt
        type(variable_t) :: var_to_advect
        integer :: j, k, i

        real, dimension(ims:ime,kms:kme,jms:jme) :: adv_dz, temp

        real, allocatable :: U_m(:,:,:), V_m(:,:,:), W_m(:,:,:), denom(:,:,:)
        ! integer :: nx, nz, ny
        !
        ! nx=size(domain%p,1)
        ! nz=size(domain%p,2)
        ! ny=size(domain%p,3)
        !
        ! if (.not.allocated(domain%tend%qv_adv)) then
        !     allocate(domain%tend%qv_adv(nx,nz,ny))
        !     domain%tend%qv_adv=0
        ! endif
        

        if (options%physics%advection==kADV_STD) then
            call adv_std_compute_wind(domain,options,dt, U_m, V_m, W_m, denom)
        else if(options%physics%advection==kADV_MPDATA) then
            call mpdata_compute_wind(domain,options,dt)
        endif
        
        adv_dz = domain%advection_dz

        ! !$OMP PARALLEL default(shared) &
        ! !$OMP private(i,j,k) &
        ! !$OMP firstprivate(adv_dz, options,U_m,V_m,W_m,denom)
        ! !$OMP DO
        do i = 1, domain%adv_vars%n_vars

            if (options%time%RK3) then
                if (options%physics%advection==kADV_STD) then
    
                    !Initial advection-tendency calculations
                    temp  = domain%adv_vars%var_list(i)%var%data_3d

                    call adv_std_advect3d(temp,domain%adv_vars%var_list(i)%var%data_3d, U_m, V_m, W_m, denom, adv_dz,t_factor_in=0.333)
                    call adv_std_advect3d(temp,domain%adv_vars%var_list(i)%var%data_3d, U_m, V_m, W_m, denom, adv_dz,t_factor_in=0.5)
    
                    !final advection call with tendency-fluxes
                    call adv_fluxcorr_advect3d(temp,domain%adv_vars%var_list(i)%var%data_3d, U_m, V_m, W_m, denom, adv_dz)
    
                    domain%adv_vars%var_list(i)%var%data_3d = temp
                else if(options%physics%advection==kADV_MPDATA) then
                    ! Not yet implemented (is it compatable w/ RK3?)
                endif
            else
                if (options%physics%advection==kADV_STD) then
                    call adv_std_advect3d(domain%adv_vars%var_list(i)%var%data_3d,domain%adv_vars%var_list(i)%var%data_3d, U_m, V_m, W_m, denom, adv_dz)
                !else if(options%physics%advection==kADV_MPDATA) then                                    
                !    call mpdata_advect3d(var, rho, jaco, dz, options)
                endif
            endif
        enddo
        ! !$OMP END DO
        ! !$OMP END PARALLEL
    end subroutine advect

end module advection
