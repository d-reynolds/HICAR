!>------------------------------------------------------------
!! Module to manage the ICAR wind field, including calls to linear winds
!! importantly it also rotates the wind field into the ICAR grid and
!! balances the U, V, and W fields for "mass" conservation
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!------------------------------------------------------------
module wind    

    use linear_theory_winds, only : linear_perturb, setup_linwinds
    ! use wind_iterative,      only : calc_iter_winds, init_iter_winds
#ifdef USE_AMGX
    use wind_iterative_amgx,     only : calc_iter_winds_amgx, init_iter_winds_amgx
#elif defined USE_PETSC
    use wind_iterative_petsc,      only : calc_iter_winds_petsc, init_iter_winds_petsc
#endif
    use iso_fortran_env, only : output_unit
    use icar_constants
    use domain_interface,  only : domain_t
    use options_interface, only : options_t
    use wind_surf, only         : calc_Sx, apply_Sx
    use wind_thermal, only      : apply_thermal_winds, init_thermal_winds
    use mod_atm_utilities,   only : calc_froude, calc_Ri, calc_dry_stability
    use array_utilities,      only : smooth_array
    use debug_module,     only : domain_check_winds

    implicit none
    private
    public:: balance_uvw, update_winds, init_winds, calc_w_real, wind_var_request, update_wind_dqdt

    integer :: ids, ide, jds, jde, kds, kde,  &
               ims, ime, jms, jme, kms, kme,  &
               its, ite, jts, jte, kts, kte,  &
               i_s, i_e, j_s, j_e

    logical :: first_wind=.True.
    real, parameter::deg2rad=0.017453293 !2*pi/360
    real, parameter :: rad2deg=57.2957779371
    real, parameter :: DEFAULT_FR_L = 1000.0
contains


    subroutine wind_linear_var_request(options)
        implicit none
        type(options_t), intent(inout) :: options

        ! List the variables that are required to be allocated for the linear wind solution
        call options%alloc_vars( &
                        [kVARS%nsquared,    kVARS%potential_temperature,   kVARS%exner,            &
                            kVARS%water_vapor, kVARS%cloud_water_mass,             kVARS%rain_mass,      &
                            kVARS%u,           kVARS%v,                       kVARS%w,                &
                            kVARS%dz, kVARS%global_terrain])

        ! List the variables that are required to be advected
        ! call options%advect_vars( &
        !               [, &
        !                , &
        !                ] )

        ! List the variables that are required for restarts with the linear wind solution
        call options%restart_vars( &
                        [kVARS%nsquared,    kVARS%potential_temperature,                           &
                            kVARS%water_vapor, kVARS%cloud_water_mass,             kVARS%rain_mass,      &
                            kVARS%u,           kVARS%v,                       kVARS%w,                &
                            kVARS%dz ])

    end subroutine

    subroutine wind_var_request(options)
        implicit none
        type(options_t), intent(inout) :: options

        if (options%wind%linear_theory) then
            call wind_linear_var_request(options)
        endif

        call options%alloc_vars([kVARS%blk_ri, kVARS%froude])

        if (options%physics%windtype == kITERATIVE_WINDS) then
            call options%alloc_vars([kVARS%wind_alpha])

            if (options%wind%alpha_const<0) call options%alloc_vars([kVARS%froude_terrain])

            call options%restart_vars([kVARS%w_real])
        endif
        
        
        if (options%wind%thermal) then
            call options%alloc_vars([kVARS%potential_temperature, kVARS%skin_temperature])
            call options%exch_vars([kVARS%skin_temperature])
            
            call options%restart_vars([kVARS%potential_temperature, kVARS%skin_temperature])
        endif

        if (options%wind%Sx) then
            call options%alloc_vars([kVARS%Sx, kVARS%TPI])
        endif
    end subroutine wind_var_request




    !------------------------------------------------------------------------------
    ! subroutine balance_uvw
    !
    ! Purpose:
    !   This subroutine balances the u, v, and w wind components in the domain
    !   by calculating the divergence of the wind field and adjusting the
    !   w component to ensure mass conservation.
    !
    ! Input:
    !   domain   - Derived data type containing the domain information
    !   options  - Derived data type containing various options
    !   update_in (optional) - Logical variable indicating which variable data array to update
    !
    ! Output:
    !   domain   - Derived data type with updated w component
    !
    ! Method:
    !   1. Calculate the divergence of the wind field
    !   2. Adjust the w component to balance the divergence
    !   3. (Optional) Perform the same for the convective wind field
    !
    ! Note:
    !   The convective wind field balancing is currently commented out.
    !
    !------------------------------------------------------------------------------    
    subroutine balance_uvw(domain, adv_den)
        ! This subroutine balances the u, v, and w wind components in the domain
        
        implicit none
        
        ! domain: a derived data type containing the domain information
        type(domain_t), intent(inout) :: domain
        
        ! options: a derived data type containing various options
        logical, intent(in) :: adv_den
                
        ! divergence: a 3D array to store the divergence of the wind field
        real, dimension(ims:ime, kms:kme, jms:jme) :: divergence
        
        ! Associate various variables from the domain data structure for easier access

        !$acc data create(divergence)
        associate(dx => domain%dx, &
                    rho => domain%vars_3d(domain%var_indx(kVARS%density)%v)%data_3d, &
                    dz => domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d, &
                    jaco_u => domain%vars_3d(domain%var_indx(kVARS%jacobian_u)%v)%data_3d, &
                    jaco_v => domain%vars_3d(domain%var_indx(kVARS%jacobian_v)%v)%data_3d, &
                    jaco_w => domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d)
        
        call calc_divergence(divergence, domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d, domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d, domain%vars_3d(domain%var_indx(kVARS%w)%v)%data_3d, &
                domain%vars_3d(domain%var_indx(kVARS%jacobian)%v)%data_3d, jaco_u, jaco_v, jaco_w, dz, dx, rho, adv_den, horz_only=.True.)
        call calc_w(domain%vars_3d(domain%var_indx(kVARS%w)%v)%data_3d, divergence, dz, jaco_w, rho, adv_den)
        
        end associate
        !$acc end data
        !$acc wait (10)
        !------------------------------------------------------------
        ! Now do the same for the convective wind field if needed
        !------------------------------------------------------------
        
        ! if (options%physics%convection > 0) then
        ! ! calculate horizontal divergence 
        ! dv = domain%v_cu(2:nx-1,i,3:ny) - domain%v_cu(2:nx-1,i,2:ny-1)
        ! du = domain%u_cu(3:nx,i,2:ny-1) - domain%u_cu(2:nx-1,i,2:ny-1)
        ! divergence = du + dv
        ! ! Then calculate w to balance
        ! if (i==1) then
        ! ! if this is the first model level start from 0 at the ground
        ! domain%w_cu(2:nx-1,i,2:ny-1) = 0 - divergence
        ! else
        ! ! else calculate w as a change from w at the level below
        ! domain%w_cu(2:nx-1,i,2:ny-1) = domain%w_cu(2:nx-1,i-1,2:ny-1) - divergence
        ! endif
        ! endif
        
    end subroutine balance_uvw

    subroutine calc_w(w,div,dz,jaco_w,rho,adv_den)
        real,    intent(inout)                                   :: w(ims:ime,kms:kme,jms:jme)
        real,    dimension(ims:ime,kms:kme,jms:jme), intent(in)  :: div, dz, jaco_w, rho
        logical, intent(in)    :: adv_den
        
        real, dimension(ims:ime,kms:kme,jms:jme) :: rho_i
        integer :: i, j, k

        !$acc data present(w, div, dz, jaco_w, rho) create(rho_i)
        !$acc parallel loop gang vector collapse(3) async(1)
        do j = jms, jme
            do k = kms, kme-1
            do i = ims, ime
                rho_i(i,k,j) = ( rho(i,k,j)*dz(i,k+1,j) + rho(i,k+1,j)*dz(i,k,j) ) / (dz(i,k,j)+dz(i,k+1,j))
            enddo
            enddo
        enddo
        
        !$acc parallel loop gang vector collapse(2) async(1)
        do j = jms, jme
            do i = ims, ime
            rho_i(i,kme,j) = rho(i,kme,j)
            enddo
        enddo
        
        !$acc parallel loop gang collapse(2) wait(1) async(10)
        do j = j_s,j_e
        do i = i_s,i_e
        !$acc loop seq
        do k = kms,kme
            if (adv_den) then
                if (k==kms) then
                    w(i,k,j) = 0 - div(i,k,j) * dz(i,k,j) &
                                                / (jaco_w(i,k,j) * rho_i(i,k,j) )
                else
                    w(i,k,j) = ( (w(i,k-1,j) * rho_i(i,k-1,j) &
                                                * jaco_w(i,k-1,j)) - div(i,k,j) * &
                                                dz(i,k,j)) / (jaco_w(i,k,j) *  rho_i(i,k,j))
                endif
            else
                if (k==kms) then
                    w(i,k,j) = (0 - div(i,k,j) * dz(i,k,j)) / (jaco_w(i,k,j) )
                else 
                    w(i,k,j) = (w(i,k-1,j) * jaco_w(i,k-1,j) - &
                                                div(i,k,j) * dz(i,k,j))/ (jaco_w(i,k,j) )
                end if
            end if
        end do
        end do
        end do

        !$acc end data
    end subroutine

    subroutine calc_divergence(div, u, v, w, jaco, jaco_u, jaco_v, jaco_w, dz, dx, rho, advect_density, horz_only)
        implicit none
        real,           intent(inout) :: div(ims:ime,kms:kme,jms:jme)
        real, dimension(ims:ime,kms:kme,jms:jme),   intent(in)    :: w, dz, jaco_w, rho, jaco
        real, dimension(ims:ime+1,kms:kme,jms:jme), intent(in)    :: u, jaco_u
        real, dimension(ims:ime,kms:kme,jms:jme+1), intent(in)    :: v, jaco_v
        real,           intent(in)    :: dx
        logical,intent(in)    :: advect_density
        logical, optional, intent(in) :: horz_only
        
        real, dimension(ims:ime,kms:kme,jms:jme) :: w_met
        real, dimension(ims:ime+1,kms:kme,jms:jme) :: u_met
        real, dimension(ims:ime,kms:kme,jms:jme+1) :: v_met
        real, dimension(ims:ime,kms:kme-1,jms:jme) :: rho_i
        logical :: horz
        integer :: i, j, k

        horz = .False.
        if (present(horz_only)) horz=horz_only

        !$acc data present(div, u, v, w, dz, jaco, jaco_u, jaco_v, jaco_w, rho, dx) create(rho_i, u_met, v_met, w_met)

        !Multiplication of U/V by metric terms, converting jacobian to staggered-grid where possible, otherwise making assumption of
        !Constant jacobian at edges
        

        if (advect_density) then
            !$acc parallel async(0)
            !$acc loop gang vector collapse(3)
            do j = jms, jme
                do k = kms, kme
                do i = ims+1, ime
                    u_met(i,k,j) = u(i,k,j) * jaco_u(i,k,j) * (rho(i-1,k,j) + rho(i,k,j))/2
                enddo
                enddo
            enddo
            !$acc loop gang vector collapse(3)
            do j = jms+1, jme
                do k = kms, kme
                do i = ims, ime
                    v_met(i,k,j) = v(i,k,j) * jaco_v(i,k,j) * (rho(i,k,j-1) + rho(i,k,j))/2
                enddo
                enddo
            enddo
            !Handle edges assuming constant density gradient
            !$acc loop gang vector collapse(2)
            do j = jms, jme
            do k = kms, kme
                u_met(ims,k,j) = u(ims,k,j) * jaco_u(ims,k,j) * (1.5*rho(ims,k,j) - 0.5*rho(ims+1,k,j))
                u_met(ime+1,k,j) = u(ime+1,k,j) * jaco_u(ime+1,k,j) * (1.5*rho(ime,k,j) - 0.5*rho(ime-1,k,j))
            enddo
            enddo
            !$acc loop gang vector collapse(2)
            do k = kms, kme
            do i = ims, ime
                v_met(i,k,jms) = v(i,k,jms) * jaco_v(i,k,jms) * (1.5*rho(i,k,jms) - 0.5*rho(i,k,jms+1))
                v_met(i,k,jme+1) = v(i,k,jme+1) * jaco_v(i,k,jme+1) * (1.5*rho(i,k,jme) - 0.5*rho(i,k,jme-1))
            enddo
            enddo

            !$acc loop gang vector collapse(3)
            do j = jms, jme
                do k = kms, kme-1
                do i = ims, ime
                    !Interpolate density to w grid
                    rho_i(i,k,j) = ( rho(i,k,j)*dz(i,k+1,j) + rho(i,k+1,j)*dz(i,k,j) ) / (dz(i,k,j)+dz(i,k+1,j))
                enddo
                enddo
            enddo
            !$acc end parallel
        else
            !$acc parallel async(0)
            !$acc loop gang vector collapse(3)
            do j = jms, jme
                do k = kms, kme
                do i = ims, ime+1
                    u_met(i,k,j) = u(i,k,j) * jaco_u(i,k,j)
                enddo
                enddo
            enddo
            !$acc loop gang vector collapse(3)
            do j = jms, jme+1
                do k = kms, kme
                do i = ims, ime
                    v_met(i,k,j) = v(i,k,j) * jaco_v(i,k,j)
                enddo
                enddo
            enddo
            !$acc end parallel
        end if

        !$acc parallel loop gang vector collapse(3) async(1) wait(0)
        do j = jms, jme
            do k = kms, kme
            do i = ims, ime
                div(i,k,j) = (u_met(i+1, k, j) - u_met(i, k, j) + &
                              v_met(i, k, j+1) - v_met(i, k, j)) / dx

            enddo
            enddo
        enddo

        if (.NOT.(horz)) then
            if (advect_density) then
                !$acc parallel async(1)
                !$acc loop gang vector collapse(3)
                do j = jms, jme
                    do k = kms, kme-1
                    do i = ims, ime
                        !Interpolate density to w grid
                        w_met(i,k,j) = w(i,k,j) * jaco_w(i,k,j) * rho_i(i,k,j)
                    enddo
                    enddo
                enddo
                !$acc loop gang vector collapse(2)
                do j = jms,jme
                do i = ims,ime
                    w_met(i,kme,j) = w(i,kme,j) * jaco_w(i,kme,j) * rho(i,kme,j)
                enddo
                enddo
                !$acc end parallel
            else
                !$acc parallel loop gang vector collapse(3) async(1)
                do j = jms, jme
                    do k = kms, kme
                        do i = ims, ime
                            w_met(i,k,j) = w(i,k,j) * jaco_w(i,k,j)
                        enddo
                    enddo
                enddo
            end if

            !$acc parallel loop gang vector collapse(3) async(1)
            do j = jms,jme
            do k = kms,kme
            do i = ims,ime
                if (k == kms) then
                    div(i, k, j) = div(i, k, j) + w_met(i, k, j)/(dz(i, k, j))
                else
                    div(i, k, j) = div(i, k, j) + &
                                   (w_met(i,k,j)-w_met(i,k-1,j))/(dz(i,k,j))
                endif
                div(i, k, j) = div(i, k, j) / jaco(i,k,j)

            enddo
            enddo
            enddo
        endif

        !$acc end data
        !$acc wait(1)
    end subroutine calc_divergence
    

    !>------------------------------------------------------------
    !! Correct for a grid that is locally rotated with respect to EW,NS
    !!
    !! Assumes forcing winds are EW, NS relative, not grid relative.
    !!
    !!------------------------------------------------------------
    subroutine make_winds_grid_relative(u, v, sintheta, costheta)
        real, intent(inout)             :: u(ims:ime+1,kms:kme,jms:jme), v(ims:ime,kms:kme,jms:jme+1)
        real, intent(in)    :: sintheta(ims:ime,jms:jme), costheta(ims:ime,jms:jme)
        
        real, dimension(ims+1:ime,kms:kme,jms:jme) :: v_ustag
        real, dimension(ims:ime,kms:kme,jms+1:jme) :: u_vstag
        
        real, dimension(ims+1:ime,jms:jme) :: costheta_ustag, sintheta_ustag
        real, dimension(ims:ime,jms+1:jme) :: costheta_vstag, sintheta_vstag
        integer :: i, j, k

        !$acc data present(u, v, sintheta, costheta) create(v_ustag, u_vstag, costheta_ustag, sintheta_ustag, costheta_vstag, sintheta_vstag)

        !$acc kernels
        v_ustag = (v(ims:ime-1,:,jms:jme)+v(ims+1:ime,:,jms:jme)+v(ims:ime-1,:,jms+1:jme+1)+v(ims+1:ime,:,jms+1:jme+1))/4
        u_vstag = (u(ims:ime,:,jms:jme-1)+u(ims:ime,:,jms+1:jme)+u(ims+1:ime+1,:,jms:jme-1)+u(ims+1:ime+1,:,jms+1:jme))/4
        
        costheta_ustag = (costheta(ims+1:ime,jms:jme)+costheta(ims:ime-1,jms:jme))/2
        sintheta_ustag = (sintheta(ims+1:ime,jms:jme)+sintheta(ims:ime-1,jms:jme))/2
        
        costheta_vstag = (costheta(ims:ime,jms+1:jme)+costheta(ims:ime,jms:jme-1))/2
        sintheta_vstag = (sintheta(ims:ime,jms+1:jme)+sintheta(ims:ime,jms:jme-1))/2

        do k = kms,kme
            u(ims,k,:)       = u(ims,k,:) * costheta_ustag(ims+1,:) + v_ustag(ims+1,k,:) * sintheta_ustag(ims+1,:)
            u(ime+1,k,:)     = u(ime+1,k,:) * costheta_ustag(ime,:) + v_ustag(ime,k,:) * sintheta_ustag(ime,:)
        
            v(:,k,jms)       = v(:,k,jms) * costheta_vstag(:,jms+1) + u_vstag(:,k,jms+1) * sintheta_vstag(:,jms+1)
            v(:,k,jme+1)     = v(:,k,jme+1) * costheta_vstag(:,jme) + u_vstag(:,k,jme) * sintheta_vstag(:,jme)
            
            u(ims+1:ime,k,:) = u(ims+1:ime,k,:) * costheta_ustag - v_ustag(:,k,:) * sintheta_ustag
            v(:,k,jms+1:jme) = v(:,k,jms+1:jme) * costheta_vstag + u_vstag(:,k,:) * sintheta_vstag
        enddo
        !$acc end kernels
        !$acc end data
    end subroutine

    subroutine apply_base_from_forcing(domain, w_var_given)
        implicit none
        type(domain_t), intent(inout) :: domain
        logical, intent(in) :: w_var_given

        integer :: i, j, k

        associate(u => domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d, &
                  v => domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d, &
                  w => domain%vars_3d(domain%var_indx(kVARS%w_real)%v)%data_3d, &
                  fu => domain%forcing_hi(domain%forcing_var_indx(kVARS%u)%v)%data_3d, &
                  fv => domain%forcing_hi(domain%forcing_var_indx(kVARS%v)%v)%data_3d, &
                  fw => domain%forcing_hi(domain%forcing_var_indx(kVARS%w_real)%v)%data_3d)
        
        !$acc data present(u, v, w, fu, fv)

        !$acc parallel
        !Compute the forcing wind field at the next update step, assuming a linear interpolation through time

        !$acc loop gang vector collapse(3)
        do j = jms, jme
            do k = kms, kme
                do i = ims,ime+1
                    u(i,k,j) = fu(i,k,j)
                enddo
            enddo
        enddo
        !$acc loop gang vector collapse(3)
        do j = jms, jme+1
            do k = kms, kme
                do i = ims,ime
                    v(i,k,j) = fv(i,k,j)
                enddo
            enddo
        enddo
        !$acc end parallel

        if (w_var_given) then
            !$acc parallel loop gang vector collapse(3) present(fw)
            do j = jms, jme
                do k = kms, kme
                    do i = ims,ime
                        w(i,k,j) = fw(i,k,j)
                    enddo
                enddo
            enddo
        else
            !If we have not read in W_real from forcing, set target w_real to 0.0. This minimizes vertical motion in solution
            !$acc parallel loop gang vector collapse(3)
            do j = jms, jme
                do k = kms, kme
                    do i = ims,ime
                        w(i,k,j) = 0.0
                    enddo
                enddo
            enddo
        endif
        !$acc end data
        end associate

        end subroutine apply_base_from_forcing

    !>------------------------------------------------------------
    !! Apply wind field physics and adjustments
    !!
    !! This will call the linear wind module if necessary, otherwise it just updates for
    !! This should ONLY be called once for each forcing step, otherwise effects will be additive.
    !!
    !!------------------------------------------------------------
    recursive subroutine update_winds(domain, options, new_input)
        implicit none
        type(domain_t), intent(inout) :: domain
        type(options_t),intent(in)    :: options
        logical, optional, intent(in)           :: new_input

        real, allocatable, dimension(:,:,:) :: div
        integer :: nx, ny, nz, it
        integer :: i, j, k
        logical :: newinput, w_var_given
        real :: wind_dt_seconds, alpha_const_val

        w_var_given = (options%forcing%wvar/="")
        wind_dt_seconds = options%wind%update_dt%seconds()
        alpha_const_val = options%wind%alpha_const

        newinput = .False.
        if (present(new_input)) newinput = new_input
        
        ! if this is a restart run, we have already read in the winds, so do not do anything
        if (first_wind .and. options%restart%restart) then
            first_wind = .False.
            return
        endif
        associate(u => domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d, &
                  u_dqdt_3d => domain%vars_3d(domain%var_indx(kVARS%u)%v)%dqdt_3d, &
                  v => domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d, &
                  v_dqdt_3d => domain%vars_3d(domain%var_indx(kVARS%v)%v)%dqdt_3d, &
                  fu => domain%forcing_hi(domain%forcing_var_indx(kVARS%u)%v)%data_3d, &
                  fu_dqdt_3d => domain%forcing_hi(domain%forcing_var_indx(kVARS%u)%v)%dqdt_3d, &
                  fv => domain%forcing_hi(domain%forcing_var_indx(kVARS%v)%v)%data_3d, &
                  fv_dqdt_3d => domain%forcing_hi(domain%forcing_var_indx(kVARS%v)%v)%dqdt_3d, &
                  fw => domain%forcing_hi(domain%forcing_var_indx(kVARS%w_real)%v)%data_3d, &
                  fw_dqdt_3d => domain%forcing_hi(domain%forcing_var_indx(kVARS%w_real)%v)%dqdt_3d, &
                  wind_alpha => domain%vars_3d(domain%var_indx(kVARS%wind_alpha)%v)%data_3d)
        !$acc data present(u, v, u_dqdt_3d, v_dqdt_3d, fu, fv, fu_dqdt_3d, fv_dqdt_3d, wind_alpha) copyin(kVARS)
        if (first_wind) then            
            first_wind = .False.

            !do this now, so that we will have some values in data_3d when calling update_stability
            call apply_base_from_forcing(domain, w_var_given)
        endif

        if (( (options%wind%alpha_const<=0 .and. (options%physics%windtype==kITERATIVE_WINDS)) .or. options%wind%Sx) ) then
            call update_stability(domain, options)
        endif

        if (newinput) then

            !$acc parallel
            !Set the state of the forcing winds to be at the end of the next wind update step
            !$acc loop gang vector collapse(3)
            do j = jms, jme
                do k = kms, kme
                    do i = ims, ime+1
                        fu(i,k,j) = fu(i,k,j) + fu_dqdt_3d(i,k,j)*wind_dt_seconds
                    enddo
                enddo
            enddo
            !$acc loop gang vector collapse(3)
            do j = jms, jme+1
                do k = kms, kme
                    do i = ims, ime
                        fv(i,k,j) = fv(i,k,j) + fv_dqdt_3d(i,k,j)*wind_dt_seconds
                    enddo
                enddo
            enddo
            !$acc end parallel
            if (w_var_given) then
                !$acc parallel loop gang vector collapse(3) present(fw, fw_dqdt_3d)
                do j = jms, jme
                    do k = kms, kme
                        do i = ims, ime
                            fw(i,k,j) = fw(i,k,j) + fw_dqdt_3d(i,k,j)*wind_dt_seconds
                        enddo
                    enddo
                enddo
            endif

            call update_winds(domain, options)

            !$acc parallel
            !Set the state of the forcing winds back to be at the current time
            !$acc loop gang vector collapse(3)
            do j = jms, jme
                do k = kms, kme
                    do i = ims, ime+1
                        fu(i,k,j) = fu(i,k,j)  - fu_dqdt_3d(i,k,j)*wind_dt_seconds

                        u_dqdt_3d(i,k,j) = u(i,k,j)
                    enddo
                enddo
            enddo
            !$acc loop gang vector collapse(3)
            do j = jms, jme+1
                do k = kms, kme
                    do i = ims, ime
                        fv(i,k,j) = fv(i,k,j)  - fv_dqdt_3d(i,k,j)*wind_dt_seconds

                        v_dqdt_3d(i,k,j) = v(i,k,j)
                    enddo
                enddo
            enddo
            !$acc end parallel
            if (w_var_given) then
                !$acc parallel loop gang vector collapse(3) present(fw, fw_dqdt_3d)
                do j = jms, jme
                    do k = kms, kme
                        do i = ims, ime
                            fw(i,k,j) = fw(i,k,j)  - fw_dqdt_3d(i,k,j)*wind_dt_seconds
                        enddo
                    enddo
                enddo
            endif

            ! call a diagnostic update to overwrite the value for density
            ! which was earlier calculated using the future pressure field
            call domain%diagnostic_update()

        endif
        if (options%general%debug) call domain_check_winds(domain, "Pre update_winds::apply_base_from_forcing")
        call apply_base_from_forcing(domain, w_var_given)
        if (options%general%debug) call domain_check_winds(domain, "Post update_winds::apply_base_from_forcing")

        ! rotate winds from cardinal directions to grid orientation (e.g. u is grid relative not truly E-W)
        call make_winds_grid_relative(domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d, domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d, domain%vars_2d(domain%var_indx(kVARS%sintheta)%v)%data_2d, domain%vars_2d(domain%var_indx(kVARS%costheta)%v)%data_2d)
        if (options%general%debug) call domain_check_winds(domain, "Post update_winds::make_winds_grid_relative")

        call domain%halo%exch_var(domain%vars_3d(domain%var_indx(kVARS%u)%v),corners=.True.)
        call domain%halo%exch_var(domain%vars_3d(domain%var_indx(kVARS%v)%v),corners=.True.)
        if (options%general%debug) call domain_check_winds(domain, "Post update_winds::make_winds_grid_relative_exch")

        if (options%wind%Sx) then
            call apply_Sx(domain%vars_4d(domain%var_indx(kVARS%Sx)%v)%data_4d,domain%vars_2d(domain%var_indx(kVARS%TPI)%v)%data_2d, &
                    domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d,domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d, &
                    domain%vars_3d(domain%var_indx(kVARS%blk_ri)%v)%data_3d,domain%vars_3d(domain%var_indx(kVARS%dzdx)%v)%data_3d, &
                    domain%vars_3d(domain%var_indx(kVARS%dzdy)%v)%data_3d, ims, ime, kms, kme, jms, jme, its, ite, jts, jte)
            call domain%halo%exch_var(domain%vars_3d(domain%var_indx(kVARS%u)%v),corners=.True.)
            call domain%halo%exch_var(domain%vars_3d(domain%var_indx(kVARS%v)%v),corners=.True.)
            if (options%general%debug) call domain_check_winds(domain, "Post update_winds::apply_Sx")
        endif 

        if (options%wind%thermal) then
            !Since this is an update call and the sensible heat fluxes can now be quite variable/patch, exchange sensible heat so that corrections are consistent
            call domain%halo_2d_send()
            call domain%halo_2d_retrieve()

            ! If model is running with a pbl scheme that supplies a 3D K_h, pass that here
            if (options%physics%boundarylayer == kPBL_YSU) then
                call apply_thermal_winds(domain%vars_2d(domain%var_indx(kVARS%skin_temperature)%v)%data_2d,domain%vars_3d(domain%var_indx(kVARS%exner)%v)%data_3d,domain%vars_3d(domain%var_indx(kVARS%potential_temperature)%v)%data_3d,  &
                                     domain%vars_3d(domain%var_indx(kVARS%z)%v)%data_3d, domain%vars_3d(domain%var_indx(kVARS%dz)%v)%data_3d,domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d, domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d,&
                                     domain%vars_3d(domain%var_indx(kVARS%dzdx)%v)%data_3d,domain%vars_3d(domain%var_indx(kVARS%dzdy)%v)%data_3d, domain%vars_3d(domain%var_indx(kVARS%coeff_heat_exchange_3d)%v)%data_3d)
            else
                call apply_thermal_winds(domain%vars_2d(domain%var_indx(kVARS%skin_temperature)%v)%data_2d,domain%vars_3d(domain%var_indx(kVARS%exner)%v)%data_3d,domain%vars_3d(domain%var_indx(kVARS%potential_temperature)%v)%data_3d,  &
                                     domain%vars_3d(domain%var_indx(kVARS%z)%v)%data_3d, domain%vars_3d(domain%var_indx(kVARS%dz)%v)%data_3d,                                                       &
                                     domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d, domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d,domain%vars_3d(domain%var_indx(kVARS%dzdx)%v)%data_3d,domain%vars_3d(domain%var_indx(kVARS%dzdy)%v)%data_3d)
            endif
            call domain%halo%exch_var(domain%vars_3d(domain%var_indx(kVARS%u)%v),corners=.True.)
            call domain%halo%exch_var(domain%vars_3d(domain%var_indx(kVARS%v)%v),corners=.True.)
            if (options%general%debug) call domain_check_winds(domain, "Post update_winds::apply_thermal_winds")
        endif 

        ! linear winds
        if (options%wind%linear_theory) then
            call linear_perturb(domain,options,options%lt%vert_smooth,.False.,options%adv%advect_density, update=.False.)
            call domain%halo%exch_var(domain%vars_3d(domain%var_indx(kVARS%u)%v),corners=.True.)
            call domain%halo%exch_var(domain%vars_3d(domain%var_indx(kVARS%v)%v),corners=.True.)
            if (options%general%debug) call domain_check_winds(domain, "Post update_winds::linear_perturb")
        endif
            
        if (options%physics%windtype==kITERATIVE_WINDS) then
            allocate(div(ims:ime,kms:kme,jms:jme))
            !$acc data create(div)
            if (alpha_const_val>0) then
                !$acc parallel loop gang vector collapse(3)
                do j = jms,jme
                    do k = kms,kme
                        do i = ims,ime
                            wind_alpha(i,k,j) = alpha_const_val
                        enddo
                    enddo
                enddo
            else
                call calc_alpha(wind_alpha, domain%vars_3d(domain%var_indx(kVARS%froude)%v)%data_3d)
            endif

            do i = 1, options%wind%wind_iterations
                call calc_idealized_wgrid(domain)

                call calc_divergence(div,domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d,domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d,domain%vars_3d(domain%var_indx(kVARS%w)%v)%data_3d, &
                                domain%vars_3d(domain%var_indx(kVARS%jacobian)%v)%data_3d, domain%vars_3d(domain%var_indx(kVARS%jacobian_u)%v)%data_3d, domain%vars_3d(domain%var_indx(kVARS%jacobian_v)%v)%data_3d,domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d,domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d,domain%dx, &
                                domain%vars_3d(domain%var_indx(kVARS%density)%v)%data_3d,options%adv%advect_density,horz_only=.False.)

#ifdef USE_AMGX                
                call calc_iter_winds_amgx(domain,wind_alpha,div,options%adv%advect_density)
#elif defined USE_PETSC
                call calc_iter_winds_petsc(domain,wind_alpha,div,options%adv%advect_density)
#endif                
                !Exchange u and v, since the outer points are not updated in above function
                call domain%halo%exch_var(domain%vars_3d(domain%var_indx(kVARS%u)%v),corners=.True.)
                call domain%halo%exch_var(domain%vars_3d(domain%var_indx(kVARS%v)%v),corners=.True.)
            end do
            !$acc end data

            if (options%general%debug) call domain_check_winds(domain, "Post update_winds::iterative_winds")
        endif
    ! elseif (options%physics%windtype==kOBRIEN_WINDS) then
    !     call Obrien_winds(domain, options, update_in=.True.)
    ! elseif (options%physics%windtype==kLINEAR_OBRIEN_WINDS) then
    !     call linear_perturb(domain,options,options%lt%vert_smooth,.False.,options%adv%advect_density, update=.False.)
    !     call Obrien_winds(domain, options, update_in=.True.)

        call balance_uvw(domain,options%adv%advect_density)
        
        !reset w_real back to the original forcing field
        call calc_w_real(domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d,  &
                domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d,      &
                domain%vars_3d(domain%var_indx(kVARS%w)%v)%data_3d,      &
                domain%vars_3d(domain%var_indx(kVARS%w_real)%v)%data_3d, &
                domain%vars_3d(domain%var_indx(kVARS%dzdx_u)%v)%data_3d, &
                domain%vars_3d(domain%var_indx(kVARS%dzdy_v)%v)%data_3d, &
                domain%vars_3d(domain%var_indx(kVARS%dzdx)%v)%data_3d,   &
                domain%vars_3d(domain%var_indx(kVARS%dzdy)%v)%data_3d,   &
                domain%vars_3d(domain%var_indx(kVARS%jacobian)%v)%data_3d)

        !$acc end data

        end associate

        if (options%general%debug) call domain_check_winds(domain, "Post update_winds::balance_uvw")

    end subroutine update_winds
    
    subroutine update_wind_dqdt(domain, dt)
        implicit none
        type(domain_t), intent(inout)  :: domain
        real, intent(in)                :: dt

        integer :: i, j, k

        associate(u => domain%vars_3d(domain%var_indx(kVARS%u)%v), &
                  v => domain%vars_3d(domain%var_indx(kVARS%v)%v))

        !$acc parallel present(u%data_3d, u%dqdt_3d, v%data_3d, v%dqdt_3d)
        !$acc loop gang vector collapse(3)
        do j = jms, jme
            do k = kms, kme
                do i = ims, ime+1
                    u%dqdt_3d(i,k,j) = (u%dqdt_3d(i,k,j)-u%data_3d(i,k,j))/dt
                enddo
            enddo
        enddo
        
        !$acc loop gang vector collapse(3)
        do j = jms, jme+1
            do k = kms, kme
                do i = ims, ime
                    v%dqdt_3d(i,k,j) = (v%dqdt_3d(i,k,j)-v%data_3d(i,k,j))/dt
                enddo
            enddo
        enddo
        !$acc end parallel
        end associate

        !If we are not using advect density, then balance_uvw will not be called every physics step, so compute a tendancy here
        ! if (.not.(options%adv%advect_density)) then
        !     domain%vars_3d(domain%var_indx(kVARS%w)%v)%dqdt_3d = (domain%vars_3d(domain%var_indx(kVARS%w)%v)%dqdt_3d-domain%vars_3d(domain%var_indx(kVARS%w)%v)%data_3d)/options%wind%update_dt%seconds()
        ! endif

    end subroutine

    !>------------------------------------------------------------
    !! Helper function to calculate the w_grid we should expect
    !! given some perscribed w_real field, which is assumed to be
    !! stored in the domain%vars_3d(domain%var_indx(kVARS%w_real)%v)%data_3d
    !! field at the time of this function call.
    !!
    !!------------------------------------------------------------
    subroutine calc_idealized_wgrid(domain)
        implicit none
        type(domain_t), intent(inout) :: domain

        real, allocatable, dimension(:,:,:) :: zero_arr

        integer :: i, j, k

        allocate(zero_arr(ims:ime,kms:kme,jms:jme))
        zero_arr = 0.0
        !$acc data copyin(zero_arr)
        !Call this, passing 0 for w_grid, to get vertical components of vertical motion
        call calc_w_real(domain%vars_3d(domain%var_indx(kVARS%u)%v) %data_3d,      &
                        domain%vars_3d(domain%var_indx(kVARS%v)%v) %data_3d,      &
                        zero_arr,      &
                        domain%vars_3d(domain%var_indx(kVARS%w)%v)%data_3d,      &
                        domain%vars_3d(domain%var_indx(kVARS%dzdx_u)%v)%data_3d, domain%vars_3d(domain%var_indx(kVARS%dzdy_v)%v)%data_3d, domain%vars_3d(domain%var_indx(kVARS%dzdx)%v)%data_3d, domain%vars_3d(domain%var_indx(kVARS%dzdy)%v)%data_3d,   &
                        domain%vars_3d(domain%var_indx(kVARS%jacobian)%v)%data_3d)
        !$acc end data


        associate(w => domain%vars_3d(domain%var_indx(kVARS%w)%v)%data_3d, &
                  w_real => domain%vars_3d(domain%var_indx(kVARS%w_real)%v)%data_3d, &
                  advection_dz => domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d, &
                  jacobian_w => domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d)

        !$acc data present(w, w_real, advection_dz, jacobian_w)
        !apply any w_real
        !$acc parallel loop gang vector collapse(3)
        do j = jms, jme
            do k = kms, kme
                do i = ims,ime
                    w(i,k,j) = (w_real(i,k,j)-w(i,k,j))
                enddo
            enddo
        enddo
        !stagger w, which was just calculated at the mass points, to the vertical k-levels, so that we can calculate divergence with it
        !$acc parallel loop gang collapse(2)
        do i = ims,ime
            do j = jms, jme
                !$acc loop seq
                do k = kms, kme-1
                    w(i,k,j) = (w(i,k,j)*advection_dz(i,k+1,j) + w(i,k+1,j)*advection_dz(i,k,j))/ &
                                        (advection_dz(i,k,j)+advection_dz(i,k+1,j))
                    w(i,k,j) = w(i,k,j) / jacobian_w(i,k,j)
                enddo
            enddo
        enddo

        !$acc parallel loop gang vector collapse(2)
        do i = ims,ime
            do j = jms, jme
                w(i,kme,j) = 0.0
            enddo
        enddo
        !$acc end data
        end associate

    end subroutine calc_idealized_wgrid
    
    subroutine calc_alpha(alpha, froude)
        implicit none
        real,    intent(in)    :: froude(ims:ime,kms:kme,jms:jme)
        real,    intent(inout) :: alpha(ims:ime,kms:kme,jms:jme)

        real :: alpha_min, alpha_max
        integer :: i, j, k

        alpha_min = 0.2
        alpha_max = 2.0
        
        !$acc data present(alpha, froude)

        !Following Moussiopoulos, et al. (1988). Bounding low Fr to avoid /0 error and negative Fr

        !$acc kernels
        alpha = 1.0 - 0.5*max(froude**4,0.00001)*(sqrt(1.0+4.0/max(froude**4,0.00001)) - 1.0) 
        alpha = sqrt(max(alpha,0.00001))
        alpha = min(max(alpha,alpha_min),alpha_max)

        !set alpha at top of domain to 0.2 to limit flux accross upper boundary
        alpha(:,kme,:) = alpha_min

        !$acc end kernels

        ! Ensure that there are no sharp transitions in alpha at boundary, 
        ! which can leak boundary effects into model (very high w_grid values result)
        if (jms==jds) then
            !$acc parallel loop gang vector collapse(3)
            do j = jms, jts-1
            do k = kms, kme
            do i = ims, ime
                alpha(i,k,j) = alpha(i,k,jts)
            enddo
            enddo
            enddo
        end if
        if (ims==ids) then
            !$acc parallel loop gang vector collapse(3)
            do j = jms, jme
            do k = kms, kme
            do i = ims, its-1
                alpha(i,k,j) = alpha(its,k,j)
            enddo
            enddo
            enddo
        end if
        if (jme==jde) then
            !$acc parallel loop gang vector collapse(3)
            do j = jte+1, jme
            do k = kms, kme
            do i = ims, ime
                alpha(i,k,j) = alpha(i,k,jte)
            enddo
            enddo
            enddo
        end if
        if (ime==ide) then
            !$acc parallel loop gang vector collapse(3)
            do j = jms, jme
            do k = kms, kme
            do i = ite+1, ime
                alpha(i,k,j) = alpha(ite,k,j)
            enddo
            enddo
            enddo
        end if

        !$acc end data

        !$acc update host(alpha)
        !smooth alpha to avoid sharp transitions
        do k = 1,3
            call smooth_array(alpha,2,ydim=3)
        enddo
        !$acc update device(alpha)

    end subroutine calc_alpha
    
    subroutine calc_w_real(u,v,w_grid,w_real,dzdx_u,dzdy_v,dzdx,dzdy,jaco)

        implicit none
        real, intent(in), dimension(ims:ime,kms:kme,jms:jme)    :: w_grid, jaco, dzdx, dzdy
        real, intent(in), dimension(ims:ime+1,kms:kme,jms:jme)  :: u, dzdx_u
        real, intent(in), dimension(ims:ime,kms:kme,jms:jme+1)  :: v, dzdy_v
        real, intent(inout), dimension(ims:ime,kms:kme,jms:jme) :: w_real
        
        integer :: k, i, j
                
        real, dimension(ims:ime,kms:kme+1,jms:jme)   :: w_grid_ext
        real, dimension(ims:ime,jms:jme)   :: currw
        real, dimension(ims:ime+1,jms:jme) :: uw
        real, dimension(ims:ime,jms:jme+1) :: vw


        !$acc data present(u, v, w_grid, w_real, dzdx_u, dzdy_v, dzdx, dzdy, jaco) create(w_grid_ext)

        !calculate the real vertical motions (including U*dzdx + V*dzdy)
        !$acc kernels
        w_grid_ext(:,kms,:) = 0
        w_grid_ext(:,kms+1:kme+1,:) = w_grid
        !$acc end kernels

        !$acc parallel loop gang vector collapse(3)
        do j = jms, jme
        do k = kms, kme
        do i = ims, ime

            ! compute the U * dz/dx component of vertical motion
            ! uw    =   u(ims:ime+1,k,jms:jme) !* dzdx_u(ims:ime+1,k,jms:jme) *

            ! ! compute the V * dz/dy component of vertical motion
            ! vw    =   v(ims:ime,k,jms:jme+1) !* dzdy_v(ims:ime,k,jms:jme+1)

            ! ! the W grid relative motion
            ! currw = w_grid(ims:ime, k, jms:jme) !* jaco_w(ims:ime, k, jms:jme)

            ! if (options%physics%convection>0) then
            !     currw = currw + domain%w_cu(2:nx-1,z,2:ny-1) * domain%dz_inter(2:nx-1,z,2:ny-1) / domain%dx
            ! endif
            
            ! compute the real vertical velocity of air by combining the different components onto the mass grid
            ! includes vertical interpolation between w_z-1/2 and w_z+1/2
            w_real(i, k, j) = (u(i,k,j) + u(i+1,k,j))*0.5*dzdx(i,k,j) &
                                                 +(v(i,k,j) + v(i,k,j+1))*0.5*dzdy(i,k,j) &
                                                 +(w_grid_ext(i, k, j) + w_grid_ext(i, k+1, j))*jaco(i,k,j) * 0.5
        end do
        end do
        end do
        !$acc end data

    end subroutine calc_w_real
    
    ! >------------------------------------------------------------
    !! O'Brien wind adjustment method
    !!
    !! SLATED FOR DEPRECATION - use iterative winds instead
    !!------------------------------------------------------------
    ! subroutine Obrien_winds(domain, options, update_in)
    !     implicit none
    !     type(domain_t), intent(inout) :: domain
    !     type(options_t),intent(in)    :: options
    !     logical, optional, intent(in) :: update_in

    !     ! interal parameters
    !     real, dimension(ims:ime,kms:kme,jms:jme)   :: div, ADJ,ADJ_coef, U_cor, V_cor, current_w
    !     real, dimension(ims:ime+1,kms:kme,jms:jme) :: current_u
    !     real, dimension(ims:ime,kms:kme,jms:jme+1) :: current_v
    !     real    :: corr_factor
    !     integer :: it, wind_k
    !     logical :: update

    !     update=.False.
    !     if (present(update_in)) update=update_in

    !     !If we are doing an update, we need to swap meta data into data_3d fields so it can be exchanged while balancing
    !     !First, we save a copy of the current data_3d so that we can substitute it back in later
    !     if (update) then
    !          current_u = domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d
    !          current_v = domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d
    !          current_w = domain%vars_3d(domain%var_indx(kVARS%w)%v)%data_3d

    !          domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d = domain%vars_3d(domain%var_indx(kVARS%u)%v)%dqdt_3d
    !          domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d = domain%vars_3d(domain%var_indx(kVARS%v)%v)%dqdt_3d
    !          domain%vars_3d(domain%var_indx(kVARS%w)%v)%data_3d = domain%vars_3d(domain%var_indx(kVARS%w)%v)%dqdt_3d
    !     endif

    !     !Do an initial exchange to make sure the U and V grids are similar for calculating w
    !     call domain%halo%exch_var(domain%vars_3d(domain%var_indx(kVARS%u)%v),do_dqdt=.False.,corners=.True.)
    !     call domain%halo%exch_var(domain%vars_3d(domain%var_indx(kVARS%v)%v),do_dqdt=.False.,corners=.True.)

    !     !First call bal_uvw to generate an initial-guess for vertical winds
    !     call balance_uvw(domain, options%adv%advect_density)

    !     ! Calculate and apply correction to w-winds
    !     wind_k = kme
    !     ! previously this code was solving for 0 vertical motion at the flat z height instead of the top boundary.
    !     ! left in for now as it could be useful to implement something similar in the future.
    !     ! however, this was also creating weird artifacts above the flat z height that need to be fixed if re-implementing.
    !     ! do k = kms,kme
    !     !     if (sum(domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(ims,1:k,jms)) > domain%smooth_height) then
    !     !         wind_k = k
    !     !         exit
    !     !     endif
    !     ! enddo
    !     ! domain%smooth_height = sum(domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(ims,:,jms))
    !     !Compute relative correction factors for U and V based on input speeds
    !     U_cor = ABS(domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d(ims:ime,:,jms:jme))/ &
    !             (ABS(domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d(ims:ime,:,jms:jme))+ABS(domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d(ims:ime,:,jms:jme)))

    !     do i = ims,ime
    !         do j = jms,jme
    !             domain%smooth_height = sum(domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i,:,j)) !
    !             do k = kms,kme
    !                 corr_factor = ((sum(domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i,1:k,j)))/domain%smooth_height)
    !                 corr_factor = min(corr_factor,1.0)
    !                 domain%vars_3d(domain%var_indx(kVARS%w)%v)%data_3d(i,k,j) = domain%vars_3d(domain%var_indx(kVARS%w)%v)%data_3d(i,k,j) - corr_factor * (domain%vars_3d(domain%var_indx(kVARS%w)%v)%data_3d(i,wind_k,j))

    !                 !if ( (domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d(i,k,j)+domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d(i,k,j)) == 0) U_cor(i,k,j) = 0.5
    !             enddo
    !         enddo
    !     enddo

    !     do k = kms,kme
    !         ! Compute this now, since it wont change in the loop
    !         ADJ_coef(:,k,:) = -2/domain%dx
    !     enddo

    !     !V_cor = 1 - U_cor


    !     U_cor = 0.5
    !     V_cor = 0.5

    !     ! Now, fixing w-winds, iterate over U/V to reduce divergence with new w-winds
    !     do it = 0,options%wind%wind_iterations
    !         !Compute divergence in new wind field
    !         call calc_divergence(div, domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d, domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d, domain%vars_3d(domain%var_indx(kVARS%w)%v)%data_3d, &
    !                             domain%vars_3d(domain%var_indx(kVARS%jacobian)%v)%data_3d, domain%vars_3d(domain%var_indx(kVARS%jacobian_u)%v)%data_3d, domain%vars_3d(domain%var_indx(kVARS%jacobian_v)%v)%data_3d, domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d,    &
    !                             domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d, domain%dx, domain%vars_3d(domain%var_indx(kVARS%density)%v)%data_3d, options%adv%advect_density)
    !         !Compute adjustment based on divergence
    !         ADJ = div/ADJ_coef

    !         !Distribute divergence among the U and V fields
    !         domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d(its+1:ite+1,:,jts:jte) = domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d(its+1:ite+1,:,jts:jte) + &
    !                                                     (ADJ(its:ite,:,jts:jte) * U_cor(its+1:ite+1,:,jts:jte))

    !         domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d(its+1:ite+1,:,jts:jte) = domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d(its+1:ite+1,:,jts:jte) - &
    !                                                     (ADJ(its+1:ite+1,:,jts:jte) * U_cor(its+1:ite+1,:,jts:jte))

    !         domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d(its:ite,:,jts+1:jte+1) = domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d(its:ite,:,jts+1:jte+1) + &
    !                                                     (ADJ(its:ite,:,jts:jte) * V_cor(its:ite,:,jts+1:jte+1))

    !         domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d(its:ite,:,jts+1:jte+1) = domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d(its:ite,:,jts+1:jte+1) - &
    !                                                     (ADJ(its:ite,:,jts+1:jte+1) * V_cor(its:ite,:,jts+1:jte+1))
    !         call domain%halo%exch_var(domain%vars_3d(domain%var_indx(kVARS%u)%v),do_dqdt=.False.,corners=.True.)
    !         call domain%halo%exch_var(domain%vars_3d(domain%var_indx(kVARS%v)%v),do_dqdt=.False.,corners=.True.)

    !     enddo

    !     !If an update loop, swap dqdt and data_3d fields back
    !     if (update) then
    !         domain%vars_3d(domain%var_indx(kVARS%u)%v)%dqdt_3d = domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d
    !         domain%vars_3d(domain%var_indx(kVARS%v)%v)%dqdt_3d = domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d
    !         domain%vars_3d(domain%var_indx(kVARS%w)%v)%dqdt_3d = domain%vars_3d(domain%var_indx(kVARS%w)%v)%data_3d

    !         domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d = current_u
    !         domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d = current_v
    !         domain%vars_3d(domain%var_indx(kVARS%w)%v)%data_3d = current_w
    !     endif

    ! end subroutine Obrien_winds

    !>------------------------------------------------------------
    !! Setup initial fields (i.e. grid relative rotation fields)
    !!
    !!------------------------------------------------------------
    subroutine init_winds(domain,options,context_chng)
        type(domain_t),  intent(inout) :: domain
        type(options_t), intent(in)    :: options
        logical, optional, intent(in) :: context_chng

        logical :: context_change

        if (present(context_chng)) then
            context_change = context_chng
        else
            context_change = .false.
        endif

        if (.not.(context_change)) first_wind = .True.

        call set_module_indices(domain)
        
        if (.not.(context_change) .and. first_wind) call allocate_winds(domain, options)

        if (options%wind%linear_theory) then
            call setup_linwinds(domain, options, .False., options%adv%advect_density)
        endif
        if (options%physics%windtype==kITERATIVE_WINDS) then
#ifdef USE_AMGX
            call init_iter_winds_amgx(domain)
#elif defined USE_PETSC
            call init_iter_winds_petsc(domain,options)
#else
            write(*,*) "ERROR: iterative winds selected but model not compiled with AMGX or PETSc libraries."
            stop
#endif
        endif

        if (options%wind%thermal) call init_thermal_winds(domain, options)

    end subroutine init_winds

    subroutine set_module_indices(domain)
        type(domain_t), intent(in) :: domain

        ids = domain%ids ; ide = domain%ide ; jds = domain%jds ; jde = domain%jde ; kds = domain%kds ; kde = domain%kde
        ims = domain%ims ; ime = domain%ime ; jms = domain%jms ; jme = domain%jme ; kms = domain%kms ; kme = domain%kme
        its = domain%its ; ite = domain%ite ; jts = domain%jts ; jte = domain%jte ; kts = domain%kts ; kte = domain%kte

        i_s = its-1
        i_e = ite+1
        j_s = jts-1
        j_e = jte+1
        
        if (ims==ids) i_s = ims
        if (ime==ide) i_e = ime
        if (jms==jds) j_s = jms
        if (jme==jde) j_e = jme

    end subroutine set_module_indices

    !>------------------------------------------------------------
    !! Allocate memory used in various wind related routines
    !!
    !!------------------------------------------------------------
    subroutine allocate_winds(domain, options)
        type(domain_t), intent(inout) :: domain
        type(options_t), intent(in) :: options
        
        if (options%wind%Sx .and. (domain%var_indx(kVARS%Sx)%v > 0)) then
            if (STD_OUT_PE) write(*,*) "    Calculating Sx and TPI for wind modification"
            call calc_Sx(domain, options)
        endif
        
        if (options%wind%alpha_const<0 .and. (options%physics%windtype==kITERATIVE_WINDS)) then
            call compute_terrain_blocking_heights(domain)
        endif

    end subroutine allocate_winds
    
    subroutine update_stability(domain, options)
        implicit none
        type(domain_t), intent(inout) :: domain
        type(options_t), intent(in) :: options

        real, dimension(ims:ime,kms:kme,jms:jme) :: wind_speed, temp_froude, u_m, v_m, u_shear, v_shear, winddir, stability
        integer, dimension(ims:ime,kms:kme,jms:jme) :: dir_indices
        
        integer :: n, ob_k, Ri_k_max
        integer :: i, j, k
        real :: z_top, z_bot, z_mean, th_top, th_bot, obstacle_height, RI_Z_MAX
        integer :: ymin, ymax, xmin, xmax, n_smoothing_passes, nsmooth_gridcells, ubound_terrain
        
        RI_Z_MAX = 100.0
        Ri_k_max = 0
        n_smoothing_passes = 5
        nsmooth_gridcells = 20 !int(500 / domain%dx)
        
        do k = kms,kme
            z_mean =SUM(domain%vars_3d(domain%var_indx(kVARS%z)%v)%data_3d(ims:ime,k,jms:jme))/SIZE(domain%vars_3d(domain%var_indx(kVARS%z)%v)%data_3d(ims:ime,k,jms:jme))
            if (z_mean > RI_Z_MAX .and. Ri_k_max==0) Ri_k_max = max(2,k-1)
        enddo

        associate(u => domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d, &
                  v => domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d, &
                  froude => domain%vars_3d(domain%var_indx(kVARS%froude)%v)%data_3d, &
                  blk_ri => domain%vars_3d(domain%var_indx(kVARS%blk_ri)%v)%data_3d, &
                  z => domain%vars_3d(domain%var_indx(kVARS%z)%v)%data_3d, &
                  potential_temperature => domain%vars_3d(domain%var_indx(kVARS%potential_temperature)%v)%data_3d, &
                  froude_terrain => domain%vars_4d(domain%var_indx(kVARS%froude_terrain)%v)%data_4d)
        !$acc data present(u, v, froude, blk_ri, z, potential_temperature) create(u_m,v_m,u_shear,v_shear,winddir,wind_speed,stability)

        !$acc kernels
        u_m(ims:ime,kms:kme,jms:jme) = (u(ims:ime,kms:kme,jms:jme) + &
                                        u(ims+1:ime+1,kms:kme,jms:jme))/2
        v_m(ims:ime,kms:kme,jms:jme) = (v(ims:ime,kms:kme,jms:jme) + &
                                        v(ims:ime,kms:kme,jms+1:jme+1))/2
        
        u_shear(:,kms,:) = u_m(:,kms+4,:)
        u_shear(:,kms+1:kme,:) = u_m(:,kms+1:kme,:) - u_m(:,kms:kme-1,:)
        v_shear(:,kms,:) = v_m(:,kms+4,:)
        v_shear(:,kms+1:kme,:) = v_m(:,kms+1:kme,:) - v_m(:,kms:kme-1,:)
        !$acc end kernels
        
        !Since we will loop up to nz-1, we set all Fr to 0.1, which will leave the upper layer as very stable
        !Since we will loop up to nz-1, we set all Ri here to 10

        !$acc parallel loop gang vector collapse(3)
        do j = jms, jme
            do k = kms, kme
                do i = ims, ime
                    wind_speed(i,k,j) = sqrt( (u_m(i,k,j))**2 + (v_m(i,k,j))**2 )
                    froude(i,k,j) = 0.1
                    blk_ri(i,k,j) = 10.0
                enddo
            enddo
        enddo

        !$acc parallel loop gang vector collapse(3)
        do j = jms, jme
            do k=kms, kme-1
                do i = ims, ime
                    th_bot = potential_temperature(i,kms,j)
                    th_top = potential_temperature(i,Ri_k_max,j)
                    z_bot  = z(i,kms,j)
                    z_top  = z(i,Ri_k_max,j)
                    stability(i,k,j) = calc_dry_stability(th_top, th_bot, z_top, z_bot) 
                    
                    blk_ri(i,k,j) =  calc_Ri(stability(i,k,j), u_shear(i,kms,j), v_shear(i,kms,j), (z_top-z_bot))
                enddo
            enddo
        enddo


        if (options%wind%alpha_const<0 .and. (options%physics%windtype==kITERATIVE_WINDS)) then


            ubound_terrain = ubound(froude_terrain,4)
            !Compute wind direction for each cell on mass grid
            !$acc parallel loop gang vector collapse(3) copyout(dir_indices)
            do j = jms, jme
                do k=kms, kme
                    do i = ims, ime
                        winddir(i,k,j) = ATAN(-u_m(i,k,j),-v_m(i,k,j))*rad2deg
                        if(winddir(i,k,j) <= 0.0) winddir(i,k,j) = winddir(i,k,j)+360
                        if(winddir(i,k,j) == 360.0) winddir(i,k,j) = 0.0
                        dir_indices(i,k,j) = max(min(int(winddir(i,k,j)/5)+1,ubound_terrain),1)                 
                    enddo
                end do
            end do
            
            !Build grid of Sx values based on wind direction at that cell
            do j = jms, jme
                do k = kms, kme
                    do i = ims, ime
                        temp_froude(i,k,j) = froude_terrain(i,k,j,dir_indices(i,k,j))
                    enddo
                end do
            end do

            !$acc parallel loop gang vector collapse(3) copyin(temp_froude)
            do j = jms, jme
                do k=kms, kme-1
                    do i = ims, ime

                        th_bot = potential_temperature(i,k,j)
                        th_top = potential_temperature(i,k+1,j)
                        z_bot  = z(i,k,j)
                        z_top  = z(i,k+1,j)

                        ! If we have an upwind obstacle, use the obstacle height to calculate a bulk Froude Number over the column
                        ! If there is nothing blocking, we calculate a local bulk Froude Number, using the local th and z indexed 
                        ! above
                        !if (.not.(temp_froude(i,k,j) == DEFAULT_FR_L)) then
                        !The height of the obstacle is calculated from the blocking terrain height (z_obst-z_loc+1000)
                        !    obstacle_height = temp_froude(i,k,j)-DEFAULT_FR_L+z_bot
                        !
                        !    do ob_k = k+1,kme
                        !        if (domain%vars_3d(domain%var_indx(kVARS%z)%v)%data_3d(i,ob_k,j) > obstacle_height) exit
                        !    enddo
                        !    ob_k = min(ob_k,kme)
                        !    th_top = domain%vars_3d(domain%var_indx(kVARS%potential_temperature)%v)%data_3d(i,ob_k,j)
                        !    z_top  = domain%vars_3d(domain%var_indx(kVARS%z)%v)%data_3d(i,ob_k,j)
                        !endif
                        stability(i,k,j) = calc_dry_stability(th_top, th_bot, z_top, z_bot, potential_temperature(i,kms,j)) 

                        !Above function calculates N^2, but froude number wants just N
                        stability(i,k,j) = sqrt(max(stability(i,k,j), 0.0))
                        froude(i,k,j) = calc_froude(stability(i,k,j), temp_froude(i,k,j), wind_speed(i,k,j))
                    enddo
                enddo
            enddo
        endif
        !$acc end data
        end associate
    end subroutine update_stability

    !>-----------------------------------------
    !> Compute a smoothed terrain varience field for use in Froude number calculation
    !>
    !------------------------------------------
    subroutine compute_terrain_blocking_heights(domain) !froude_terrain, terrain)
        implicit none
        type(domain_t), intent(inout) :: domain
        real, dimension(1:72,ims:ime,kms:kme,jms:jme)   :: temp_ft_array
        integer           :: ang, is, js
        integer           :: i, j, k
        integer           :: rear_ang, fore_ang, test_ang, rear_ang_diff, fore_ang_diff, ang_diff, k_max, window_rear, window_fore, window_width
        integer :: x, y, azm_index
        integer :: xs,xe, ys,ye, n, np
        real              :: pt_height, temp_ft, maxFTVal, azm
                
        temp_ft_array = -100000.0

        ! then compute the range of terrain (max-min) in a given window
        do i=ims, ime
            do j=jms, jme
                do k=kms,kme
                    if (k == kms) then
                        pt_height = domain%vars_2d(domain%var_indx(kVARS%neighbor_terrain)%v)%data_2d(i,j)
                    else if (k > kms) then
                        pt_height = pt_height + domain%vars_3d(domain%var_indx(kVARS%dz_interface)%v)%data_3d(i,k,j)
                    end if
                    
                    do is = domain%ihs,domain%ihe
                        do js = domain%jhs,domain%jhe
                        
                            !Compute azimuth ind of point
                            azm = atan2(1.0*(is-i),1.0*(js-j))*rad2deg
                            if(azm < 0) then
                                azm = 360+azm
                            else if(azm >= 360.0) then
                                azm=0.0
                            endif
                            azm_index = int(azm/5)+1
                        
                            temp_ft = domain%vars_2d(domain%var_indx(kVARS%neighbor_terrain)%v)%data_2d(is,js) - pt_height
                            
                            if (temp_ft > temp_ft_array(azm_index,i,k,j)) then
                                                        
                                !Only save scale length if it is greater than the default -- otherwise copy that over
                                if (temp_ft > DEFAULT_FR_L) then
                                    temp_ft_array(azm_index,i,k,j) = temp_ft
                                else
                                    temp_ft_array(azm_index,i,k,j) = DEFAULT_FR_L
                                end if
                            end if
                        enddo
                    enddo

                    !After finding Fr-Terrain in each absolute direction around grid cell, 
                    !Pick max for each 20º window and perform interpolation to other directions if necesarry
                    
                    rear_ang = 1 
                    fore_ang = 1
                    
                    if (.not.( all((temp_ft_array(:,i,k,j) <= -100000.0)) )) then
                    
                        !Perform 20º window max search
                        window_width = 2
                        do ang = 1, 72
                            window_rear = ang-window_width
                            window_fore = ang+window_width
                        
                            if (ang <= window_width) then
                                window_rear = 72-(window_width-ang)
                                
                                maxFTVal = maxval(temp_ft_array(window_rear:72,i,k,j))

                                if (maxval(temp_ft_array(1:window_fore,i,k,j)) > maxFTVal) then
                                    maxFTVal = maxval(temp_ft_array(1:window_fore,i,k,j))
                                end if
                                
                            else if ( ang >= (72-(window_width-1)) ) then
                                window_fore = window_width-(72-ang)
                                
                                maxFTVal = maxval(temp_ft_array(window_rear:72,i,k,j))

                                if (maxval(temp_ft_array(1:window_fore,i,k,j)) > maxFTVal) then
                                    maxFTVal = maxval(temp_ft_array(1:window_fore,i,k,j))
                                end if
                            else
                                maxFTVal = maxval(temp_ft_array(window_rear:window_fore,i,k,j))
                            end if
                            domain%vars_4d(domain%var_indx(kVARS%froude_terrain)%v)%data_4d(i,k,j,ang) = maxFTVal
                        end do                    
                    
                        do ang = 1, 72
                            !Determine indices for interpolation
                            if ( (ang==fore_ang) ) then
                                !Update indices for interpolated Fr-Terrain's
                                rear_ang = ang
                            
                                fore_ang = ang+1
                                if (fore_ang > 72) fore_ang = 1
                                
                                do while (domain%vars_4d(domain%var_indx(kVARS%froude_terrain)%v)%data_4d(i,k,j,fore_ang) <= -100000.0)
                                    fore_ang = fore_ang+1
                                    if (fore_ang > 72) fore_ang = 1
                                end do
                            
                            end if
                            
                            if (ang==1) then
                                rear_ang = 72
                                do while(domain%vars_4d(domain%var_indx(kVARS%froude_terrain)%v)%data_4d(i,k,j,rear_ang) <= -100000.0)
                                    rear_ang = rear_ang-1
                                end do
                            end if
                    
                            !If we did not calculate Fr-Terrain for a given direction
                            if (domain%vars_4d(domain%var_indx(kVARS%froude_terrain)%v)%data_4d(i,k,j,ang) == -100000.0) then
                                !Weight the two surrounding Fr-Terrain values based on our angular-distance to them
                                rear_ang_diff = ang-rear_ang
                                fore_ang_diff = fore_ang-ang
                                ang_diff = fore_ang-rear_ang
                        
                                !Handle wrap-around case
                                if (ang > fore_ang) then
                                    fore_ang_diff = fore_ang+(72-ang)
                                    ang_diff = fore_ang+(72-rear_ang)
                                end if
                        
                                !Interpolation, linearly-weighted by angular-distance from values
                                domain%vars_4d(domain%var_indx(kVARS%froude_terrain)%v)%data_4d(i,k,j,ang) = (domain%vars_4d(domain%var_indx(kVARS%froude_terrain)%v)%data_4d(i,k,j,rear_ang)*fore_ang_diff + &
                                                    domain%vars_4d(domain%var_indx(kVARS%froude_terrain)%v)%data_4d(i,k,j,fore_ang)*rear_ang_diff)/ang_diff

                            end if
                        end do

                    else
                        !IF we only have -100000 for all entries, set to dz
                        domain%vars_4d(domain%var_indx(kVARS%froude_terrain)%v)%data_4d(i,k,j,:) = domain%vars_3d(domain%var_indx(kVARS%dz_interface)%v)%data_3d(i,k,j)
                    end if
                enddo

            enddo
        enddo
                                                               
        if (domain%jms==(domain%jds)) domain%vars_4d(domain%var_indx(kVARS%froude_terrain)%v)%data_4d(:,:,jms,:) = &
                                        domain%vars_4d(domain%var_indx(kVARS%froude_terrain)%v)%data_4d(:,:,jms+1,:)
                        
        if (domain%ims==(domain%ids)) domain%vars_4d(domain%var_indx(kVARS%froude_terrain)%v)%data_4d(ims,:,:,:) = &
                                        domain%vars_4d(domain%var_indx(kVARS%froude_terrain)%v)%data_4d(ims+1,:,:,:)

        if (domain%jme==(domain%jde)) domain%vars_4d(domain%var_indx(kVARS%froude_terrain)%v)%data_4d(:,:,jme,:) = &
                                        domain%vars_4d(domain%var_indx(kVARS%froude_terrain)%v)%data_4d(:,:,jme-1,:)

        if (domain%ime==(domain%ide)) domain%vars_4d(domain%var_indx(kVARS%froude_terrain)%v)%data_4d(ime,:,:,:) = &
                                        domain%vars_4d(domain%var_indx(kVARS%froude_terrain)%v)%data_4d(ime-1,:,:,:)
                                 

    end subroutine compute_terrain_blocking_heights


    !>------------------------------------------------------------
    !! Provides a routine to deallocate memory allocated in allocate_winds
    !!
    !!------------------------------------------------------------
    ! subroutine finalize_winds(domain)
    !     type(domain_t), intent(inout) :: domain
    !
    !     if (allocated(domain%sintheta)) then
    !         deallocate(domain%sintheta)
    !     endif
    !     if (allocated(domain%costheta)) then
    !         deallocate(domain%costheta)
    !     endif
    !     if (allocated(domain%dzdx)) then
    !         deallocate(domain%dzdx)
    !     endif
    !     if (allocated(domain%dzdy)) then
    !         deallocate(domain%dzdy)
    !     endif
    !
    ! end subroutine finalize_winds
end module wind
