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
    use wind_iterative,      only : calc_iter_winds, init_iter_winds

    !use mod_blocking,        only : update_froude_number, initialize_blocking
    use variable_interface,       only : variable_t
    use iso_fortran_env
    use icar_constants
    use domain_interface,  only : domain_t
    use boundary_interface,only : boundary_t
    use options_interface, only : options_t
    use grid_interface,    only : grid_t
    use wind_surf, only         : calc_Sx, apply_Sx
    use wind_thermal, only      : apply_thermal_winds, init_thermal_winds
    use io_routines, only : io_read, io_write
    use mod_atm_utilities,   only : calc_froude, calc_Ri, calc_dry_stability
    use array_utilities,      only : smooth_array
    use time_delta_object,  only : time_delta_t

    implicit none
    private
    public:: balance_uvw, update_winds, init_winds, calc_w_real, wind_var_request, update_wind_dqdt

    integer :: ids, ide, jds, jde, kds, kde,  &
               ims, ime, jms, jme, kms, kme,  &
               its, ite, jts, jte, kts, kte, j, k, i, &
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

        if (options%physics%windtype == kWIND_LINEAR .or. &
            options%physics%windtype == kLINEAR_OBRIEN_WINDS .or. &
            options%physics%windtype == kLINEAR_ITERATIVE_WINDS) then
            call wind_linear_var_request(options)
        endif

        call options%alloc_vars([kVARS%blk_ri, kVARS%froude])

        if (options%physics%windtype == kITERATIVE_WINDS .or. options%physics%windtype == kLINEAR_ITERATIVE_WINDS) then
            call options%alloc_vars([kVARS%wind_alpha, kVARS%froude_terrain])

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
    subroutine balance_uvw(domain, options)
        ! This subroutine balances the u, v, and w wind components in the domain
        
        implicit none
        
        ! domain: a derived data type containing the domain information
        type(domain_t), intent(inout) :: domain
        
        ! options: a derived data type containing various options
        type(options_t), intent(in) :: options
                
        ! divergence: a 3D array to store the divergence of the wind field
        real, dimension(ims:ime, kms:kme, jms:jme) :: divergence
        
        ! Associate various variables from the domain data structure for easier access
        associate(dx => domain%dx, &
                    rho => domain%vars_3d(domain%var_indx(kVARS%density)%v)%data_3d, &
                    dz => domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d, &
                    jaco_u => domain%vars_3d(domain%var_indx(kVARS%jacobian_u)%v)%data_3d, &
                    jaco_v => domain%vars_3d(domain%var_indx(kVARS%jacobian_v)%v)%data_3d, &
                    jaco_w => domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d)
        
        call calc_divergence(divergence, domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d, domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d, domain%vars_3d(domain%var_indx(kVARS%w)%v)%data_3d, &
                domain%vars_3d(domain%var_indx(kVARS%jacobian)%v)%data_3d, jaco_u, jaco_v, jaco_w, dz, dx, rho, options, horz_only=.True.)
        call calc_w(domain%vars_3d(domain%var_indx(kVARS%w)%v)%data_3d, divergence, dz, jaco_w, rho, options%adv%advect_density)
        
        end associate
                
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
        
        real, dimension(ims:ime,kms:kme-1,jms:jme) :: rho_i
        
        rho_i(:,kms:kme-1,:) = ( rho(:,kms:kme-1,:)*dz(:,kms+1:kme,:) + rho(:,kms+1:kme,:)*dz(:,kms:kme-1,:) ) / (dz(:,kms:kme-1,:)+dz(:,kms+1:kme,:))
        
        w = 0
        
        do k = kms,kme
            if (adv_den) then
                if (k==kms) then
                    w(i_s:i_e,k,j_s:j_e) = 0 - div(i_s:i_e,k,j_s:j_e) * dz(i_s:i_e,k,j_s:j_e) &
                                                / (jaco_w(i_s:i_e,k,j_s:j_e) * rho_i(i_s:i_e,k,j_s:j_e) )
                elseif (k==kme) then
                    w(i_s:i_e,k,j_s:j_e) = ((w(i_s:i_e,k-1,j_s:j_e) * rho_i(i_s:i_e,k-1,j_s:j_e) &
                                                * jaco_w(i_s:i_e,k-1,j_s:j_e)) - div(i_s:i_e,k,j_s:j_e) * &
                                                dz(i_s:i_e,k,j_s:j_e)) / (jaco_w(i_s:i_e,k,j_s:j_e) * rho(i_s:i_e,k,j_s:j_e))
                else
                    w(i_s:i_e,k,j_s:j_e) = ( (w(i_s:i_e,k-1,j_s:j_e) * rho_i(i_s:i_e,k-1,j_s:j_e) * &
                                jaco_w(i_s:i_e,k-1,j_s:j_e)) - div(i_s:i_e,k,j_s:j_e) * dz(i_s:i_e,k,j_s:j_e)) / &
                                (jaco_w(i_s:i_e,k,j_s:j_e) *  rho_i(i_s:i_e,k,j_s:j_e) )
                endif
            else
                if (k==kms) then
                    w(i_s:i_e,k,j_s:j_e) = (0 - div(i_s:i_e,k,j_s:j_e) * dz(i_s:i_e,k,j_s:j_e)) / (jaco_w(i_s:i_e,k,j_s:j_e) )
                else 
                    w(i_s:i_e,k,j_s:j_e) = (w(i_s:i_e,k-1,j_s:j_e) * jaco_w(i_s:i_e,k-1,j_s:j_e) - &
                                                div(i_s:i_e,k,j_s:j_e) * dz(i_s:i_e,k,j_s:j_e))/ (jaco_w(i_s:i_e,k,j_s:j_e) )
                end if
            end if
        end do

    end subroutine

    subroutine calc_divergence(div, u, v, w, jaco, jaco_u, jaco_v, jaco_w, dz, dx, rho, options, horz_only)
        implicit none
        real,           intent(inout) :: div(ims:ime,kms:kme,jms:jme)
        real, dimension(ims:ime,kms:kme,jms:jme),   intent(in)    :: w, dz, jaco_w, rho, jaco
        real, dimension(ims:ime+1,kms:kme,jms:jme), intent(in)    :: u, jaco_u
        real, dimension(ims:ime,kms:kme,jms:jme+1), intent(in)    :: v, jaco_v
        real,           intent(in)    :: dx
        type(options_t),intent(in)    :: options
        logical, optional, intent(in) :: horz_only
        
        real, dimension(ims:ime,kms:kme,jms:jme) :: diff_U, diff_V, w_met
        real, dimension(ims:ime+1,kms:kme,jms:jme) :: u_met
        real, dimension(ims:ime,kms:kme,jms:jme+1) :: v_met
        real, dimension(ims:ime,kms:kme-1,jms:jme) :: rho_i
        logical :: horz

        horz = .False.
        if (present(horz_only)) horz=horz_only


        !Multiplication of U/V by metric terms, converting jacobian to staggered-grid where possible, otherwise making assumption of
        !Constant jacobian at edges
        
        if (options%adv%advect_density) then
            u_met(ims+1:ime,:,:) = u(ims+1:ime,:,:) * jaco_u(ims+1:ime,:,:) * (rho(ims:ime-1,:,:) + rho(ims+1:ime,:,:))/2
            u_met(ims,:,:) = u(ims,:,:) * jaco_u(ims,:,:) * (1.5*rho(ims,:,:) - 0.5*rho(ims+1,:,:))
            u_met(ime+1,:,:) = u(ime+1,:,:) * jaco_u(ime+1,:,:) * (1.5*rho(ime,:,:) - 0.5*rho(ime-1,:,:))

            v_met(:,:,jms+1:jme) = v(:,:,jms+1:jme) * jaco_v(:,:,jms+1:jme) * (rho(:,:,jms:jme-1) + rho(:,:,jms+1:jme))/2
            v_met(:,:,jms) = v(:,:,jms) * jaco_v(:,:,jms) * (1.5*rho(:,:,jms) - 0.5*rho(:,:,jms+1))
            v_met(:,:,jme+1) = v(:,:,jme+1) * jaco_v(:,:,jme+1) * (1.5*rho(:,:,jme) - 0.5*rho(:,:,jme-1))

            rho_i(:,kms:kme-1,:) = ( rho(:,kms:kme-1,:)*dz(:,kms+1:kme,:) + rho(:,kms+1:kme,:)*dz(:,kms:kme-1,:) ) / (dz(:,kms:kme-1,:)+dz(:,kms+1:kme,:))
        else
            u_met = u * jaco_u
            v_met = v * jaco_v
        end if

        diff_U = u_met(ims+1:ime+1, :, jms:jme) - u_met(ims:ime, :, jms:jme)
        diff_V = v_met(ims:ime, :, jms+1:jme+1) - v_met(ims:ime, :, jms:jme)

        div(ims:ime,kms:kme,jms:jme) = (diff_U+diff_V) /(dx)

        if (.NOT.(horz)) then
            if (options%adv%advect_density) then
                w_met(:,kme,:) = w(:,kme,:) * jaco_w(:,kme,:) * rho(:,kme,:)
                w_met(:,kms:kme-1,:) = w(:,kms:kme-1,:) * jaco_w(:,kms:kme-1,:) * rho_i
            else
                w_met = w*jaco_w
            end if

            do k = kms,kme
                if (k == kms) then
                    div(ims:ime, k, jms:jme) = div(ims:ime, k, jms:jme) + w_met(ims:ime, k, jms:jme)/(dz(ims:ime, k, jms:jme))
                else
                    div(ims:ime, k, jms:jme) = div(ims:ime, k, jms:jme) + &
                                   (w_met(ims:ime,k,jms:jme)-w_met(ims:ime,k-1,jms:jme))/(dz(ims:ime,k,jms:jme))
                endif
            enddo

            div = div / jaco
        endif

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
        
    end subroutine

    subroutine apply_base_from_forcing(domain, options)
        implicit none
        type(domain_t), intent(inout) :: domain
        type(options_t), intent(in)    :: options

        !Compute the forcing wind field at the next update step, assuming a linear interpolation through time
        do i = ims,ime+1
            do k = kms, kme
                do j = jms, jme
                    domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d(i,k,j) = domain%forcing_hi(domain%forcing_var_indx(kVARS%u)%v)%data_3d(i,k,j)
                enddo
            enddo
        enddo
        
        do i = ims,ime
            do k = kms, kme
                do j = jms, jme+1
                    domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d(i,k,j) = domain%forcing_hi(domain%forcing_var_indx(kVARS%v)%v)%data_3d(i,k,j)
                enddo
            enddo
        enddo
        
        if (.not.(options%forcing%wvar=="")) then
            do i = ims,ime
                do k = kms, kme
                    do j = jms, jme
                        domain%vars_3d(domain%var_indx(kVARS%w_real)%v)%data_3d(i,k,j) = domain%forcing_hi(domain%forcing_var_indx(kVARS%w_real)%v)%data_3d(i,k,j)
                    enddo
                enddo
            enddo
        else
            !If we have not read in W_real from forcing, set target w_real to 0.0. This minimizes vertical motion in solution
            domain%vars_3d(domain%var_indx(kVARS%w_real)%v)%data_3d = 0.0
        endif
            
        end subroutine apply_base_from_forcing

    !>------------------------------------------------------------
    !! Apply wind field physics and adjustments
    !!
    !! This will call the linear wind module if necessary, otherwise it just updates for
    !! This should ONLY be called once for each forcing step, otherwise effects will be additive.
    !!
    !!------------------------------------------------------------
    subroutine update_winds(domain, options, new_input)
        implicit none
        type(domain_t), intent(inout) :: domain
        type(options_t),intent(in)    :: options
        logical, optional, intent(in)           :: new_input

        real, allocatable, dimension(:,:,:) :: div
        integer :: nx, ny, nz, it
        logical :: newinput

        newinput = .False.
        if (present(new_input)) newinput = new_input
        
        if (first_wind) then            
            first_wind = .False.
            ! do nothing if this is a restart run -- then we have read in the already calculated winds
            if (options%restart%restart) return

            !do this now, so that we will have some values in data_3d when calling update_stability
            call apply_base_from_forcing(domain, options)
        endif

        if (( (options%wind%alpha_const<=0 .and. (options%physics%windtype==kLINEAR_ITERATIVE_WINDS &
            .or. options%physics%windtype==kITERATIVE_WINDS)) .or. options%wind%Sx) ) then
            call update_stability(domain, options)
        endif

        if (newinput) then

            !Set the state of the forcing winds to be at the end of the next wind update step
            domain%forcing_hi(domain%forcing_var_indx(kVARS%u)%v)%data_3d = domain%forcing_hi(domain%forcing_var_indx(kVARS%u)%v)%data_3d + &
                        domain%forcing_hi(domain%forcing_var_indx(kVARS%u)%v)%dqdt_3d*options%wind%update_dt%seconds()
            domain%forcing_hi(domain%forcing_var_indx(kVARS%v)%v)%data_3d = domain%forcing_hi(domain%forcing_var_indx(kVARS%v)%v)%data_3d + &
                        domain%forcing_hi(domain%forcing_var_indx(kVARS%v)%v)%dqdt_3d*options%wind%update_dt%seconds()

            if (.not.(options%forcing%wvar=="")) domain%forcing_hi(domain%forcing_var_indx(kVARS%w_real)%v)%data_3d = domain%forcing_hi(domain%forcing_var_indx(kVARS%w_real)%v)%data_3d + domain%forcing_hi(domain%forcing_var_indx(kVARS%w_real)%v)%dqdt_3d*options%wind%update_dt%seconds()
            
            call update_winds(domain, options)

            !Set the state of the forcing winds back to be at the current time
            domain%forcing_hi(domain%forcing_var_indx(kVARS%u)%v)%data_3d = domain%forcing_hi(domain%forcing_var_indx(kVARS%u)%v)%data_3d - &
                        domain%forcing_hi(domain%forcing_var_indx(kVARS%u)%v)%dqdt_3d*options%wind%update_dt%seconds()
            domain%forcing_hi(domain%forcing_var_indx(kVARS%v)%v)%data_3d = domain%forcing_hi(domain%forcing_var_indx(kVARS%v)%v)%data_3d - &
                        domain%forcing_hi(domain%forcing_var_indx(kVARS%v)%v)%dqdt_3d*options%wind%update_dt%seconds()

            if (.not.(options%forcing%wvar=="")) domain%forcing_hi(domain%forcing_var_indx(kVARS%w_real)%v)%data_3d = domain%forcing_hi(domain%forcing_var_indx(kVARS%w_real)%v)%data_3d - domain%forcing_hi(domain%forcing_var_indx(kVARS%w_real)%v)%dqdt_3d*options%wind%update_dt%seconds()

            domain%vars_3d(domain%var_indx(kVARS%u)%v)%dqdt_3d = domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d
            domain%vars_3d(domain%var_indx(kVARS%v)%v)%dqdt_3d = domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d

            ! call a diagnostic update to overwrite the value for density
            ! which was earlier calculated using the future pressure field
            call domain%diagnostic_update()

        endif

        call apply_base_from_forcing(domain, options)

        ! rotate winds from cardinal directions to grid orientation (e.g. u is grid relative not truly E-W)
        call make_winds_grid_relative(domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d, domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d, domain%vars_2d(domain%var_indx(kVARS%sintheta)%v)%data_2d, domain%vars_2d(domain%var_indx(kVARS%costheta)%v)%data_2d)

        call domain%halo%exch_var(domain%vars_3d(domain%var_indx(kVARS%u)%v),corners=.True.)
        call domain%halo%exch_var(domain%vars_3d(domain%var_indx(kVARS%v)%v),corners=.True.)
        if (options%wind%Sx) then
            call apply_Sx(domain%vars_4d(domain%var_indx(kVARS%Sx)%v)%data_4d,domain%vars_2d(domain%var_indx(kVARS%TPI)%v)%data_2d, &
                    domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d,domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d, &
                    domain%vars_3d(domain%var_indx(kVARS%blk_ri)%v)%data_3d,domain%vars_3d(domain%var_indx(kVARS%dzdx)%v)%data_3d, &
                    domain%vars_3d(domain%var_indx(kVARS%dzdy)%v)%data_3d, ims, ime, kms, kme, jms, jme, its, ite, jts, jte)
            call domain%halo%exch_var(domain%vars_3d(domain%var_indx(kVARS%u)%v),corners=.True.)
            call domain%halo%exch_var(domain%vars_3d(domain%var_indx(kVARS%v)%v),corners=.True.)
        endif 

        if (options%wind%thermal) then
            !Since this is an update call and the sensible heat fluxes can now be quite variable/patch, exchange sensible heat so that corrections are consistent
            call domain%halo_2d_send_batch()
            call domain%halo_2d_retrieve_batch()

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
        endif 

        ! linear winds
        if (options%physics%windtype==kWIND_LINEAR .or. options%physics%windtype==kLINEAR_ITERATIVE_WINDS .or. options%physics%windtype==kITERATIVE_WINDS) then
            if (options%physics%windtype==kWIND_LINEAR .or. options%physics%windtype==kLINEAR_ITERATIVE_WINDS) then
                call linear_perturb(domain,options,options%lt%vert_smooth,.False.,options%adv%advect_density, update=.True.)
            endif
            
            if (options%physics%windtype==kLINEAR_ITERATIVE_WINDS .or. options%physics%windtype==kITERATIVE_WINDS) then
                allocate(div(ims:ime,kms:kme,jms:jme))

                if (options%wind%alpha_const>0) then
                    domain%vars_3d(domain%var_indx(kVARS%wind_alpha)%v)%data_3d = options%wind%alpha_const
                else
                    call calc_alpha(domain%vars_3d(domain%var_indx(kVARS%wind_alpha)%v)%data_3d, domain%vars_3d(domain%var_indx(kVARS%froude)%v)%data_3d)
                endif


                call calc_idealized_wgrid(domain)

                call calc_divergence(div,domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d,domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d,domain%vars_3d(domain%var_indx(kVARS%w)%v)%data_3d, &
                                domain%vars_3d(domain%var_indx(kVARS%jacobian)%v)%data_3d, domain%vars_3d(domain%var_indx(kVARS%jacobian_u)%v)%data_3d, domain%vars_3d(domain%var_indx(kVARS%jacobian_v)%v)%data_3d,domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d,domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d,domain%dx, &
                                domain%vars_3d(domain%var_indx(kVARS%density)%v)%data_3d,options,horz_only=.False.)

                call calc_iter_winds(domain,domain%vars_3d(domain%var_indx(kVARS%wind_alpha)%v)%data_3d,div,options%adv%advect_density)

                ! -- Do this all a second time to minimize discretization artifacts
                call calc_idealized_wgrid(domain)

                call calc_divergence(div,domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d,domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d,domain%vars_3d(domain%var_indx(kVARS%w)%v)%data_3d, &
                        domain%vars_3d(domain%var_indx(kVARS%jacobian)%v)%data_3d, domain%vars_3d(domain%var_indx(kVARS%jacobian_u)%v)%data_3d, domain%vars_3d(domain%var_indx(kVARS%jacobian_v)%v)%data_3d,domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d,domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d,domain%dx, &
                        domain%vars_3d(domain%var_indx(kVARS%density)%v)%data_3d,options,horz_only=.False.)

                call calc_iter_winds(domain,domain%vars_3d(domain%var_indx(kVARS%wind_alpha)%v)%data_3d,div,options%adv%advect_density)
            endif
        elseif (options%physics%windtype==kOBRIEN_WINDS) then
            call Obrien_winds(domain, options, update_in=.True.)
        elseif (options%physics%windtype==kLINEAR_OBRIEN_WINDS) then
            call linear_perturb(domain,options,options%lt%vert_smooth,.False.,options%adv%advect_density, update=.True.)
            call Obrien_winds(domain, options, update_in=.True.)
        endif

        call balance_uvw(domain,options)
        
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

        
    end subroutine update_winds
    
    subroutine update_wind_dqdt(domain, options)
        implicit none
        type(domain_t), intent(inout)  :: domain
        type(options_t), intent(in)    :: options

        domain%vars_3d(domain%var_indx(kVARS%u)%v)%dqdt_3d = (domain%vars_3d(domain%var_indx(kVARS%u)%v)%dqdt_3d-domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d)/options%wind%update_dt%seconds()
        domain%vars_3d(domain%var_indx(kVARS%v)%v)%dqdt_3d = (domain%vars_3d(domain%var_indx(kVARS%v)%v)%dqdt_3d-domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d)/options%wind%update_dt%seconds()
        
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

        !Call this, passing 0 for w_grid, to get vertical components of vertical motion
        call calc_w_real(domain%vars_3d(domain%var_indx(kVARS%u)%v) %data_3d,      &
                        domain%vars_3d(domain%var_indx(kVARS%v)%v) %data_3d,      &
                        domain%vars_3d(domain%var_indx(kVARS%w)%v)%data_3d*0.0,      &
                        domain%vars_3d(domain%var_indx(kVARS%w)%v)%data_3d,      &
                        domain%vars_3d(domain%var_indx(kVARS%dzdx_u)%v)%data_3d, domain%vars_3d(domain%var_indx(kVARS%dzdy_v)%v)%data_3d, domain%vars_3d(domain%var_indx(kVARS%dzdx)%v)%data_3d, domain%vars_3d(domain%var_indx(kVARS%dzdy)%v)%data_3d,   &
                        domain%vars_3d(domain%var_indx(kVARS%jacobian)%v)%data_3d)

        !apply any w_real
        domain%vars_3d(domain%var_indx(kVARS%w)%v)%data_3d = (domain%vars_3d(domain%var_indx(kVARS%w_real)%v)%data_3d-domain%vars_3d(domain%var_indx(kVARS%w)%v)%data_3d)

        !stagger w, which was just calculated at the mass points, to the vertical k-levels, so that we can calculate divergence with it
        do i = ims,ime
            do j = jms, jme
                do k = kms, kme-1
                    domain%vars_3d(domain%var_indx(kVARS%w)%v)%data_3d(i,k,j) = (domain%vars_3d(domain%var_indx(kVARS%w)%v)%data_3d(i,k,j)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i,k+1,j) + &
                                                            domain%vars_3d(domain%var_indx(kVARS%w)%v)%data_3d(i,k+1,j)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i,k,j))/ &
                                        (domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i,k,j)+domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i,k+1,j))
                enddo
            enddo
        enddo
        domain%vars_3d(domain%var_indx(kVARS%w)%v)%data_3d(:,kme,:) = 0.0
        domain%vars_3d(domain%var_indx(kVARS%w)%v)%data_3d = domain%vars_3d(domain%var_indx(kVARS%w)%v)%data_3d / domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d

    end subroutine calc_idealized_wgrid
    
    subroutine calc_alpha(alpha, froude)
        implicit none
        real,    intent(in)    :: froude(ims:ime,kms:kme,jms:jme)
        real,    intent(inout) :: alpha(ims:ime,kms:kme,jms:jme)
        integer :: k

        real :: alpha_min, alpha_max

        alpha_min = 0.2
        alpha_max = 2.0
        
        !Following Moussiopoulos, et al. (1988). Bounding low Fr to avoid /0 error and negative Fr
        alpha = 1.0 - 0.5*max(froude**4,0.00001)*(sqrt(1.0+4.0/max(froude**4,0.00001)) - 1.0) 
        alpha = sqrt(alpha)
        alpha = min(max(alpha,alpha_min),alpha_max)

        ! Ensure that there are no sharp transitions in alpha at boundary, 
        ! which can leak boundary effects into model (very high w_grid values result)
        if (jms==jds) then
            do j = jms, jts-1
                alpha(:,:,j) = alpha(:,:,jts)
            end do
        end if
        if (ims==ids) then
            do i = ims, its-1
                alpha(i,:,:) = alpha(its,:,:)
            end do
        end if
        if (jme==jde) then
            do j = jte+1, jme
                alpha(:,:,j) = alpha(:,:,jte)
            end do
        end if
        if (ime==ide) then
            do i = ite+1, ime
                alpha(i,:,:) = alpha(ite,:,:)
            end do
        end if

        !smooth alpha to avoid sharp transitions
        do k = 1,3
            call smooth_array(alpha,2,ydim=3)
        enddo

        !set alpha at top of domain to 0.2 to limit flux accross upper boundary
        alpha(:,kme,:) = alpha_min

    end subroutine calc_alpha
    
    subroutine calc_w_real(u,v,w_grid,w_real,dzdx_u,dzdy_v,dzdx,dzdy,jaco)

        implicit none
        real, intent(in), dimension(ims:ime,kms:kme,jms:jme)    :: w_grid, jaco, dzdx, dzdy
        real, intent(in), dimension(ims:ime+1,kms:kme,jms:jme)  :: u, dzdx_u
        real, intent(in), dimension(ims:ime,kms:kme,jms:jme+1)  :: v, dzdy_v
        real, intent(inout), dimension(ims:ime,kms:kme,jms:jme) :: w_real
        
        integer :: z
                
        real, dimension(ims:ime,jms:jme)   :: lastw
        real, dimension(ims:ime,jms:jme)   :: currw
        real, dimension(ims:ime+1,jms:jme) :: uw
        real, dimension(ims:ime,jms:jme+1) :: vw

        !calculate the real vertical motions (including U*dzdx + V*dzdy)
        lastw = 0!w_grid(i_s:i_e, kms, j_s:j_e) * jaco_w(i_s:i_e, kms, j_s:j_e)
        do z = kms, kme

            ! compute the U * dz/dx component of vertical motion
            uw    =   u(ims:ime+1,z,jms:jme) !* dzdx_u(ims:ime+1,z,jms:jme) *

            ! compute the V * dz/dy component of vertical motion
            vw    =   v(ims:ime,z,jms:jme+1) !* dzdy_v(ims:ime,z,jms:jme+1)

            ! the W grid relative motion
            currw = w_grid(ims:ime, z, jms:jme) !* jaco_w(ims:ime, z, jms:jme)

            ! if (options%physics%convection>0) then
            !     currw = currw + domain%w_cu(2:nx-1,z,2:ny-1) * domain%dz_inter(2:nx-1,z,2:ny-1) / domain%dx
            ! endif
            
            ! compute the real vertical velocity of air by combining the different components onto the mass grid
            ! includes vertical interpolation between w_z-1/2 and w_z+1/2
            w_real(ims:ime, z, jms:jme) = (uw(ims:ime,:) + uw(ims+1:ime+1,:))*0.5*dzdx(:,z,:) &
                                                 +(vw(:,jms:jme) + vw(:,jms+1:jme+1))*0.5*dzdy(:,z,:) &
                                                 +(lastw + currw)*jaco(:,z,:) * 0.5
            lastw = currw ! could avoid this memcopy cost using pointers or a single manual loop unroll
        end do
                
    end subroutine calc_w_real
    
    subroutine Obrien_winds(domain, options, update_in)
        implicit none
        type(domain_t), intent(inout) :: domain
        type(options_t),intent(in)    :: options
        logical, optional, intent(in) :: update_in

        ! interal parameters
        real, dimension(ims:ime,kms:kme,jms:jme)   :: div, ADJ,ADJ_coef, U_cor, V_cor, current_w
        real, dimension(ims:ime+1,kms:kme,jms:jme) :: current_u
        real, dimension(ims:ime,kms:kme,jms:jme+1) :: current_v
        real    :: corr_factor
        integer :: it, wind_k
        logical :: update

        update=.False.
        if (present(update_in)) update=update_in

        !If we are doing an update, we need to swap meta data into data_3d fields so it can be exchanged while balancing
        !First, we save a copy of the current data_3d so that we can substitute it back in later
        if (update) then
             current_u = domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d
             current_v = domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d
             current_w = domain%vars_3d(domain%var_indx(kVARS%w)%v)%data_3d

             domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d = domain%vars_3d(domain%var_indx(kVARS%u)%v)%dqdt_3d
             domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d = domain%vars_3d(domain%var_indx(kVARS%v)%v)%dqdt_3d
             domain%vars_3d(domain%var_indx(kVARS%w)%v)%data_3d = domain%vars_3d(domain%var_indx(kVARS%w)%v)%dqdt_3d
        endif

        !Do an initial exchange to make sure the U and V grids are similar for calculating w
        call domain%halo%exch_var(domain%vars_3d(domain%var_indx(kVARS%u)%v),do_dqdt=.False.,corners=.True.)
        call domain%halo%exch_var(domain%vars_3d(domain%var_indx(kVARS%v)%v),do_dqdt=.False.,corners=.True.)

        !First call bal_uvw to generate an initial-guess for vertical winds
        call balance_uvw(domain, options)

        ! Calculate and apply correction to w-winds
        wind_k = kme
        ! previously this code was solving for 0 vertical motion at the flat z height instead of the top boundary.
        ! left in for now as it could be useful to implement something similar in the future.
        ! however, this was also creating weird artifacts above the flat z height that need to be fixed if re-implementing.
        ! do k = kms,kme
        !     if (sum(domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(ims,1:k,jms)) > domain%smooth_height) then
        !         wind_k = k
        !         exit
        !     endif
        ! enddo
        ! domain%smooth_height = sum(domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(ims,:,jms))
        !Compute relative correction factors for U and V based on input speeds
        U_cor = ABS(domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d(ims:ime,:,jms:jme))/ &
                (ABS(domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d(ims:ime,:,jms:jme))+ABS(domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d(ims:ime,:,jms:jme)))

        do i = ims,ime
            do j = jms,jme
                domain%smooth_height = sum(domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i,:,j)) !
                do k = kms,kme
                    corr_factor = ((sum(domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i,1:k,j)))/domain%smooth_height)
                    corr_factor = min(corr_factor,1.0)
                    domain%vars_3d(domain%var_indx(kVARS%w)%v)%data_3d(i,k,j) = domain%vars_3d(domain%var_indx(kVARS%w)%v)%data_3d(i,k,j) - corr_factor * (domain%vars_3d(domain%var_indx(kVARS%w)%v)%data_3d(i,wind_k,j))

                    !if ( (domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d(i,k,j)+domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d(i,k,j)) == 0) U_cor(i,k,j) = 0.5
                enddo
            enddo
        enddo

        do k = kms,kme
            ! Compute this now, since it wont change in the loop
            ADJ_coef(:,k,:) = -2/domain%dx
        enddo

        !V_cor = 1 - U_cor


        U_cor = 0.5
        V_cor = 0.5

        ! Now, fixing w-winds, iterate over U/V to reduce divergence with new w-winds
        do it = 0,options%wind%wind_iterations
            !Compute divergence in new wind field
            call calc_divergence(div, domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d, domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d, domain%vars_3d(domain%var_indx(kVARS%w)%v)%data_3d, &
                                domain%vars_3d(domain%var_indx(kVARS%jacobian)%v)%data_3d, domain%vars_3d(domain%var_indx(kVARS%jacobian_u)%v)%data_3d, domain%vars_3d(domain%var_indx(kVARS%jacobian_v)%v)%data_3d, domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d,    &
                                domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d, domain%dx, domain%vars_3d(domain%var_indx(kVARS%density)%v)%data_3d, options)
            !Compute adjustment based on divergence
            ADJ = div/ADJ_coef

            !Distribute divergence among the U and V fields
            domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d(its+1:ite+1,:,jts:jte) = domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d(its+1:ite+1,:,jts:jte) + &
                                                        (ADJ(its:ite,:,jts:jte) * U_cor(its+1:ite+1,:,jts:jte))

            domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d(its+1:ite+1,:,jts:jte) = domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d(its+1:ite+1,:,jts:jte) - &
                                                        (ADJ(its+1:ite+1,:,jts:jte) * U_cor(its+1:ite+1,:,jts:jte))

            domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d(its:ite,:,jts+1:jte+1) = domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d(its:ite,:,jts+1:jte+1) + &
                                                        (ADJ(its:ite,:,jts:jte) * V_cor(its:ite,:,jts+1:jte+1))

            domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d(its:ite,:,jts+1:jte+1) = domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d(its:ite,:,jts+1:jte+1) - &
                                                        (ADJ(its:ite,:,jts+1:jte+1) * V_cor(its:ite,:,jts+1:jte+1))
            call domain%halo%exch_var(domain%vars_3d(domain%var_indx(kVARS%u)%v),do_dqdt=.False.,corners=.True.)
            call domain%halo%exch_var(domain%vars_3d(domain%var_indx(kVARS%v)%v),do_dqdt=.False.,corners=.True.)

        enddo

        !If an update loop, swap dqdt and data_3d fields back
        if (update) then
            domain%vars_3d(domain%var_indx(kVARS%u)%v)%dqdt_3d = domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d
            domain%vars_3d(domain%var_indx(kVARS%v)%v)%dqdt_3d = domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d
            domain%vars_3d(domain%var_indx(kVARS%w)%v)%dqdt_3d = domain%vars_3d(domain%var_indx(kVARS%w)%v)%data_3d

            domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d = current_u
            domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d = current_v
            domain%vars_3d(domain%var_indx(kVARS%w)%v)%data_3d = current_w
        endif

    end subroutine Obrien_winds

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

        if (options%physics%windtype==kWIND_LINEAR .or. &
                 options%physics%windtype==kLINEAR_OBRIEN_WINDS .or. &
                 options%physics%windtype==kLINEAR_ITERATIVE_WINDS) then
            call setup_linwinds(domain, options, .False., options%adv%advect_density)
        endif
        if (options%physics%windtype==kITERATIVE_WINDS .or. options%physics%windtype==kLINEAR_ITERATIVE_WINDS) then
            call init_iter_winds(domain,options)
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
        ! note w is special cased because it does not have a forcing variable, so it is not necessarily allocated automatically
        if ((domain%var_indx(kVARS%w)%v > 0)) then
            if (.not.(allocated(domain%vars_3d(domain%var_indx(kVARS%w)%v)%dqdt_3d))) then
                allocate(domain%vars_3d(domain%var_indx(kVARS%w)%v)%dqdt_3d(ims:ime,kms:kme,jms:jme))
                domain%vars_3d(domain%var_indx(kVARS%w)%v)%dqdt_3d = 0
            endif
        endif
        
        if (options%wind%Sx .and. (domain%var_indx(kVARS%Sx)%v > 0)) then
            if (STD_OUT_PE) write(*,*) "    Calculating Sx and TPI for wind modification"
            call calc_Sx(domain, options)
        endif
        
        if (options%wind%alpha_const<0 .and. (options%physics%windtype==kLINEAR_ITERATIVE_WINDS .or. options%physics%windtype==kITERATIVE_WINDS)) then
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
        real :: z_top, z_bot, z_mean, th_top, th_bot, obstacle_height, RI_Z_MAX
        integer :: ymin, ymax, xmin, xmax, n_smoothing_passes, nsmooth_gridcells
        
        RI_Z_MAX = 100.0
        Ri_k_max = 0
        n_smoothing_passes = 5
        nsmooth_gridcells = 20 !int(500 / domain%dx)
        
        u_m(ims:ime,kms:kme,jms:jme) = (domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d(ims:ime,kms:kme,jms:jme) + &
                                        domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d(ims+1:ime+1,kms:kme,jms:jme))/2
        v_m(ims:ime,kms:kme,jms:jme) = (domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d(ims:ime,kms:kme,jms:jme) + &
                                        domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d(ims:ime,kms:kme,jms+1:jme+1))/2
        
        u_shear(:,kms,:) = u_m(:,kms+4,:)
        u_shear(:,kms+1:kme,:) = u_m(:,kms+1:kme,:) - u_m(:,kms:kme-1,:)
        v_shear(:,kms,:) = v_m(:,kms+4,:)
        v_shear(:,kms+1:kme,:) = v_m(:,kms+1:kme,:) - v_m(:,kms:kme-1,:)

        ! If we want variable alpha, then calculate froude terrrain
        if (options%wind%alpha_const<0 .and. (options%physics%windtype==kLINEAR_ITERATIVE_WINDS .or. options%physics%windtype==kITERATIVE_WINDS)) then
                        
            
            !Compute wind direction for each cell on mass grid
            do i = ims, ime
                do j = jms, jme
                    do k=kms, kme
                        winddir(i,k,j) = ATAN(-u_m(i,k,j),-v_m(i,k,j))*rad2deg
                        if(winddir(i,k,j) <= 0.0) winddir(i,k,j) = winddir(i,k,j)+360
                        if(winddir(i,k,j) == 360.0) winddir(i,k,j) = 0.0
                        dir_indices(i,k,j) = int(winddir(i,k,j)/5)+1                    
                    enddo
                end do
            end do
            
            !Build grid of Sx values based on wind direction at that cell
            do i = ims, ime
                do j = jms, jme
                    do k=kms, kme
                        temp_froude(i,k,j) = domain%vars_4d(domain%var_indx(kVARS%froude_terrain)%v)%data_4d(i,k,j,dir_indices(i,k,j))
                    enddo
                end do
            end do
        endif

        wind_speed = sqrt( (u_m)**2 + (v_m)**2 )
        
        !Since we will loop up to nz-1, we set all Fr to 0.1, which will leave the upper layer as very stable
        !Since we will loop up to nz-1, we set all Ri here to 10

        do i = ims,ime
            do k = kms, kme
                do j = jms,jme
                    domain%vars_3d(domain%var_indx(kVARS%froude)%v)%data_3d(i,k,j) = 0.1
                    domain%vars_3d(domain%var_indx(kVARS%blk_ri)%v)%data_3d(i,k,j) = 10.0
                enddo
            enddo
        enddo


        do k = kms,kme
            z_mean =SUM(domain%vars_3d(domain%var_indx(kVARS%z)%v)%data_3d(ims:ime,k,jms:jme))/SIZE(domain%vars_3d(domain%var_indx(kVARS%z)%v)%data_3d(ims:ime,k,jms:jme))
            if (z_mean > RI_Z_MAX .and. Ri_k_max==0) Ri_k_max = max(2,k-1)
        enddo


        do i = ims,ime
            do j = jms,jme
                do k = kms,kme-1
                    th_bot = domain%vars_3d(domain%var_indx(kVARS%potential_temperature)%v)%data_3d(i,kms,j)
                    th_top = domain%vars_3d(domain%var_indx(kVARS%potential_temperature)%v)%data_3d(i,Ri_k_max,j)
                    z_bot  = domain%vars_3d(domain%var_indx(kVARS%z)%v)%data_3d(i,kms,j)
                    z_top  = domain%vars_3d(domain%var_indx(kVARS%z)%v)%data_3d(i,Ri_k_max,j)
                    stability(i,k,j) = calc_dry_stability(th_top, th_bot, z_top, z_bot) 
                    
                    domain%vars_3d(domain%var_indx(kVARS%blk_ri)%v)%data_3d(i,k,j) =  calc_Ri(stability(i,k,j), u_shear(i,kms,j), v_shear(i,kms,j), (z_top-z_bot))
                enddo
            enddo
        enddo


        if (options%wind%alpha_const<0 .and. (options%physics%windtype==kLINEAR_ITERATIVE_WINDS .or. options%physics%windtype==kITERATIVE_WINDS)) then
            do i = ims,ime
                do j = jms,jme
                    do k = kms,kme-1
                
                        th_bot = domain%vars_3d(domain%var_indx(kVARS%potential_temperature)%v)%data_3d(i,k,j)
                        th_top = domain%vars_3d(domain%var_indx(kVARS%potential_temperature)%v)%data_3d(i,k+1,j)
                        z_bot  = domain%vars_3d(domain%var_indx(kVARS%z)%v)%data_3d(i,k,j)
                        z_top  = domain%vars_3d(domain%var_indx(kVARS%z)%v)%data_3d(i,k+1,j)
                        
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
                        stability(i,k,j) = calc_dry_stability(th_top, th_bot, z_top, z_bot, domain%vars_3d(domain%var_indx(kVARS%potential_temperature)%v)%data_3d(i,kms,j)) 

                        !Above function calculates N^2, but froude number wants just N
                        stability(i,k,j) = sqrt(max(stability(i,k,j), 0.0))
                        domain%vars_3d(domain%var_indx(kVARS%froude)%v)%data_3d(i,k,j) = calc_froude(stability(i,k,j), temp_froude(i,k,j), wind_speed(i,k,j))
                    enddo
                enddo
            enddo
        endif

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
        integer           :: rear_ang, fore_ang, test_ang, rear_ang_diff, fore_ang_diff, ang_diff, k_max, window_rear, window_fore, window_width
        integer :: x, y, azm_index
        integer :: xs,xe, ys,ye, n, np
        real, allocatable :: temp_terrain(:,:)
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
