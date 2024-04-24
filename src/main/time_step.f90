!> ----------------------------------------------------------------------------
!!  Main time stepping module.
!!  Calculates a stable time step (dt) and loops over physics calls
!!  Also updates boundaries every time step.
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!! ----------------------------------------------------------------------------
module time_step
    use iso_fortran_env
    use mpi_f08
    use data_structures             ! *_type  types and kCONSTANTS
    use microphysics,               only : mp
    use advection,                  only : advect
    use mod_atm_utilities,          only : exner_function, compute_ivt, compute_iq
    use convection,                 only : convect
    use land_surface,               only : lsm
    use surface_layer,              only : sfc
    use planetary_boundary_layer,   only : pbl
    use radiation,                  only : rad
    use wind,                       only : balance_uvw, update_winds, update_wind_dqdt
    use domain_interface,           only : domain_t
    use boundary_interface,         only : boundary_t
    use options_interface,          only : options_t
    use debug_module,               only : domain_check
    use string,                     only : str
    use timer_interface,            only : timer_t

    implicit none
    private
    double precision, parameter  :: DT_BIG = 36000.0
    double precision  :: future_dt_seconds = DT_BIG
    public :: step

contains

    !>------------------------------------------------------------
    !!  Calculate the maximum stable time step given some CFL criteria
    !!
    !!  For each grid cell, find the mean of the wind speeds from each
    !!  direction * sqrt(3) for the 3D advection CFL limited time step
    !!  Also find the maximum wind speed anywhere in the domain to check
    !!  against a 1D advection limit.
    !!
    !! @param dx  [ scalar ]        horizontal grid cell width  [m]
    !! @param u   [nx+1 x nz x ny]  east west wind speeds       [m/s]
    !! @param v   [nx x nz x ny+1]  North South wind speed      [m/s]
    !! @param w   [nx x nz x ny]    vertical wind speed         [m/s]
    !! @param CFL [ scalar ]        CFL limit to use (e.g. 1.0)
    !! @return dt [ scalar ]        Maximum stable time step    [s]
    !!
    !!------------------------------------------------------------
    function compute_dt(dx, u, v, w, rho, dz, ims, ime, kms, kme, jms, jme, its, ite, jts, jte, CFL, cfl_strictness, use_density) result(dt)
        real,       intent(in)                   :: dx
        real,       intent(in), dimension(ims:ime+1,kms:kme,jms:jme) :: u 
        real,       intent(in), dimension(ims:ime,kms:kme,jms:jme+1) :: v
        real,       intent(in), dimension(ims:ime,kms:kme,jms:jme)   :: w, rho
        real,       intent(in), dimension(kms:kme)     :: dz
        integer,    intent(in)                   :: ims, ime, kms, kme, jms, jme, its, ite, jts, jte
        real,       intent(in)                   :: CFL
        integer,    intent(in)                   :: cfl_strictness
        logical,    intent(in)                   :: use_density
        
        ! output value
        real :: dt
        ! locals
        real :: three_d_cfl = 0.577350269 ! = sqrt(3)/3
        integer :: i, j, k, zoffset
        real :: maxwind3d, maxwind1d, current_wind, sqrt3

        sqrt3 = sqrt(3.0) * 1.001 ! with a safety factor

        maxwind1d = 0
        maxwind3d = 0

        if (cfl_strictness==1) then
            ! to ensure we are stable for 1D advection:
            if (use_density) then
                !maxwind1d = max( maxval(abs(u(2:,:,:) / (rho*dz*dx) )), maxval(abs(v(:,:,2:) / (rho*dz*dx))) )
                !maxwind1d = max( maxwind1d, maxval(abs(w/(rho*dz*dx))) )
            else
                maxwind1d = max( maxval(abs(u)), maxval(abs(v)))
                maxwind1d = max( maxwind1d, maxval(abs(w)))
            endif

            maxwind3d = maxwind1d * sqrt3
        else if (cfl_strictness==5) then

            if (use_density) then
                !maxwind1d = maxval(abs(u(2:,:,:) / (rho*dz*dx) )) &
                !          + maxval(abs(v(:,:,2:) / (rho*dz*dx) )) &
                !          + maxval(abs(w(:,:, :) / (rho*dz*dx) ))
            else
                maxwind3d = maxval(abs(u)) + maxval(abs(v)) + maxval(abs(w))
            endif

        else
            ! to ensure we are stable for 3D advection we'll use the average "max" wind speed
            ! but that average has to be divided by sqrt(3) for stability in 3 dimensional advection
            do j=jts,jte
                do k=kms,kme
                    if (k==kms) then
                        zoffset = 0
                    else
                        zoffset = -1
                    endif

                    do i=its,ite
                        ! just compute the sum of the wind speeds, but take the max over the two
                        ! faces of the grid cell (e.g. east and west sides)
                        ! this will be divided by 3 later by three_d_cfl
                        if (use_density) then
                            !current_wind = (max(abs(u(i,k,j)), abs(u(i+1,k,j))) &
                            !              + max(abs(v(i,k,j)), abs(v(i,k,j+1))) &
                            !              + max(abs(w(i,k,j)), abs(w(i,k+zoffset,j))) ) &
                            !              / (rho(i,k,j) * dz(i,k,j) * dx)
                        else
                            !current_wind = max(abs(u(i,k,j)), abs(u(i+1,k,j))) / dx &
                            !              +max(abs(v(i,k,j)), abs(v(i,k,j+1))) / dx &
                            !              +max(abs(w(i,k,j)), abs(w(i,k+zoffset,j))) / dz(k)
                                          
                            current_wind = max(( max( abs(u(i,k,j)), abs(u(i+1,k,j)) ) / dx), &
                                               ( max( abs(v(i,k,j)), abs(v(i,k,j+1)) ) / dx), &
                                               ( max( abs(w(i,k,j)), abs(w(i,k+zoffset,j)) ) / dz(k) ))
                        endif
                        maxwind3d = max(maxwind3d, current_wind)
                    ENDDO
                ENDDO
            ENDDO

            if (cfl_strictness==2) then
                ! effectively divides by 3 to take the mean and multiplies by the sqrt(3) for the 3D advection limit
                maxwind3d = maxwind3d * three_d_cfl

                ! to ensure we are stable for 1D advection:
                if (use_density) then
                    !maxwind1d = max( maxval(abs(u(2:,:,:) / (rho*dz*dx) )), maxval(abs(v(:,:,2:) / (rho*dz*dx))) )
                    !maxwind1d = max( maxwind1d, maxval(abs(w/(rho*dz*dx))) )
                else
                    maxwind1d = max( maxval(abs(u)), maxval(abs(v)))
                    maxwind1d = max( maxwind1d, maxval(abs(w)))
                endif
                ! also insure stability for 1D advection
                maxwind3d = max(maxwind1d,maxwind3d)

            ! else if (cfl_strictness==3) then
            !   leave maxwind3d as the sum of the max winds
            ! This should be the default, does it need to be multiplied by sqrt(3)?
            elseif (cfl_strictness==4) then
                maxwind3d = maxwind3d * sqrt3
            endif

        endif

        !TESTING: Do we need to multiply maxwind3d by sqrt3 as the comment above suggests?
        ! maxwind3d = maxwind3d * sqrt3

        dt = CFL / maxwind3d

        ! If we have too small a time step throw an error
        ! something is probably wrong in the physics or input data
        if (dt<1e-1) then
            write(*,*) "dt   = ", dt
            write(*,*) "Umax = ", maxval(abs(u))
            write(*,*) "Vmax = ", maxval(abs(v))
            write(*,*) "Wmax = ", maxval(abs(w))
            error stop "ERROR time step too small"
        endif

    end function compute_dt


    !>------------------------------------------------------------
    !!  Prints progress to the terminal if requested
    !!
    !! @param current_time  the current state of the model time
    !! @param end_time      the end of the current full time step (when step will return)
    !! @param time_step     length of full time step to be integrated by step
    !! @param dt            numerical timestep to print for information
    !!
    !!------------------------------------------------------------
    subroutine print_progress(current_time, end_time, time_step, dt, last_time)
        implicit none
        type(Time_type),    intent(in)    :: current_time,    end_time
        type(time_delta_t), intent(in)    :: time_step,       dt
        real,               intent(inout) :: last_time

        type(time_delta_t) :: progress_dt
        real :: time_percent

        ! first compute the current time until reaching the end
        progress_dt  = (end_time - current_time)

        ! convert that to a percentage of the total time required
        time_percent = 100 - progress_dt%seconds() / time_step%seconds()  * 100

        ! finally if it has been at least 5% of the time since the last time we printed output, print output
        if (time_percent > (last_time + 5.0)) then
            last_time = last_time + 5.0
            ! this used to just use the nice $ (or advance="NO") trick, but at least with some mpi implementations, it buffers this output until it crashes
            write(*,"(A,f5.1,A,A)") char(13), max(0.0, time_percent)," %  dt=",trim(dt%as_string())
        endif

    end subroutine print_progress

    !>------------------------------------------------------------
    !! Update the numerical timestep to use
    !!
    !! @param dt            numerical timestep to use
    !! @param options       set options for how to update the time step
    !! @param domain        the full domain structure (need winds to compute dt)
    !! @param end_time      the end of the current full time step (when step will return)
    !!
    !!------------------------------------------------------------
    subroutine update_dt(dt, options, domain, update)
        implicit none
        type(time_delta_t), intent(inout) :: dt
        type(options_t),    intent(in)    :: options
        type(domain_t),     intent(in)    :: domain
        logical, optional,  intent(in)    :: update

        logical :: update_in
        double precision                  :: present_dt_seconds, seconds_out
        ! compute internal timestep dt to maintain stability
        ! courant condition for 3D advection. 
                
        update_in = .False.
        if (present(update)) update_in = update
        
        ! If this is the first step (future_dt_seconds has not yet been set)
        if (future_dt_seconds == DT_BIG) then
            present_dt_seconds = compute_dt(domain%dx, domain%u%data_3d, domain%v%data_3d, &
                            domain%w%data_3d, domain%density%data_3d, options%parameters%dz_levels, &
                            domain%ims, domain%ime, domain%kms, domain%kme, domain%jms, domain%jme, &
                            domain%its, domain%ite, domain%jts, domain%jte, &
                            options%time_options%cfl_reduction_factor, &
                            cfl_strictness=options%time_options%cfl_strictness, use_density=.false.)
        else
            present_dt_seconds = future_dt_seconds
        endif
        
        future_dt_seconds = compute_dt(domain%dx, domain%u%dqdt_3d, domain%v%dqdt_3d, &
                        domain%w%dqdt_3d, domain%density%data_3d, options%parameters%dz_levels, &
                        domain%ims, domain%ime, domain%kms, domain%kme, domain%jms, domain%jme, &
                        domain%its, domain%ite, domain%jts, domain%jte, &
                        options%time_options%cfl_reduction_factor, &
                        cfl_strictness=options%time_options%cfl_strictness, use_density=.false.)
                
        !Minimum dt is min(present_dt_seconds, future_dt_seconds). Then reduce this accross all compute processes
        call MPI_Allreduce(min(present_dt_seconds, future_dt_seconds), seconds_out, 1, MPI_DOUBLE, MPI_MIN, domain%IO_comms)
        
        ! Set dt to the outcome of reduce
        call dt%set(seconds=seconds_out)
        if (this_image()==1) write(*,*) 'time_step: ',dt%seconds()

    end subroutine update_dt
    
    !>------------------------------------------------------------
    !!  Step forward one IO time step.
    !!
    !!  Calculated the internal model time step to satisfy the CFL criteria,
    !!  then updates all forcing update increments for that dt and loops through
    !!  time calling physics modules.
    !!  Also checks to see if it is time to write a model output file.
    !!
    !! @param domain    domain data structure containing model state
    !! @param options   model options structure
    !! @param bc        model boundary conditions data structure
    !! @param next_output   Next time to write an output file (in "model_time")
    !!
    !!------------------------------------------------------------
    subroutine step(domain, forcing, end_time, options, mp_timer, adv_timer, rad_timer, lsm_timer, pbl_timer, exch_timer, send_timer, ret_timer, wait_timer, forcing_timer, diagnostic_timer, wind_bal_timer, wind_timer)
        implicit none
        type(domain_t),     intent(inout)   :: domain
        type(boundary_t),   intent(inout)   :: forcing
        type(Time_type),    intent(in)      :: end_time
        type(options_t),    intent(in)      :: options

        type(timer_t),      intent(inout)   :: mp_timer, adv_timer, rad_timer, lsm_timer, pbl_timer, exch_timer, send_timer, ret_timer, wait_timer, forcing_timer, diagnostic_timer, wind_bal_timer, wind_timer

        real            :: last_print_time
        real, save      :: last_wind_update

        type(time_delta_t) :: time_step_size, dt_saver
        type(time_delta_t), save      :: dt

        last_print_time = 0.0
        call dt_saver%set(seconds=0.0)
        
        time_step_size = end_time - domain%model_time
        
        ! Initialize to just over update_dt to force an update on first loop
        if (domain%model_time==options%parameters%start_time) last_wind_update = options%wind%update_dt%seconds() + 1

        ! now just loop over internal timesteps computing all physics in order (operator splitting...)
        do while (domain%model_time < end_time)
            
            call send_timer%start()
            call domain%halo%halo_3d_send_batch(exch_vars=domain%exch_vars, adv_vars=domain%adv_vars)
            call send_timer%stop()

            !Determine dt
            if (last_wind_update >= options%wind%update_dt%seconds()) then

                call domain%diagnostic_update(options)

                call wind_timer%start()
                call update_winds(domain, forcing, options)
                call wind_timer%stop()
                !Now that new winds have been calculated, get new time step in seconds, and see if they require adapting the time step
                ! Note that there will currently be some discrepancy between using the current density and whatever density will be at 
                ! the next time step, but we assume that it is negligable
                ! and that using a CFL criterion < 1.0 will cover this
                call update_dt(dt, options, domain)
                call update_wind_dqdt(domain, options)

                last_wind_update = 0.0
            endif
            !call update_dt(dt, options, domain, end_time)
            ! Make sure we don't over step the forcing or output period
            if ((domain%model_time + dt) > end_time) then
                dt_saver = dt
                dt = end_time - domain%model_time
            endif
            
            ! ! apply/update boundary conditions including internal wind and pressure changes.
            call forcing_timer%start()
            call domain%apply_forcing(forcing,options,real(dt%seconds()))
            call forcing_timer%stop()

            call diagnostic_timer%start()
            call domain%diagnostic_update(options)
            call diagnostic_timer%stop()

            if (options%parameters%advect_density) then
                ! if using advect_density winds need to be balanced at each update
                call wind_bal_timer%start()
                call balance_uvw(domain,options)
                call wind_bal_timer%stop()
            endif

            ! if an interactive run was requested than print status updates everytime at least 5% of the progress has been made
            if (options%parameters%interactive .and. (this_image()==1)) then
                call print_progress(domain%model_time, end_time, time_step_size, dt, last_print_time)
            endif
            ! this if is to avoid round off errors causing an additional physics call that won't really do anything

            if (real(dt%seconds()) > 1e-3) then

                if (options%parameters%debug) call domain_check(domain, "img: "//trim(str(this_image()))//" init", fix=.True.)

                ! first process the halo section of the domain (currently hard coded at 1 should come from domain?)
                call rad_timer%start()
                call rad(domain, options, real(dt%seconds()))
                if (options%parameters%debug) call domain_check(domain, "img: "//trim(str(this_image()))//" rad(domain", fix=.True.)
                call rad_timer%stop()

                call lsm_timer%start()
                call sfc(domain, options, real(dt%seconds()))!, halo=1)
                call lsm(domain, options, real(dt%seconds()))!, halo=1)
                if (options%parameters%debug) call domain_check(domain, "img: "//trim(str(this_image()))//" lsm")
                call lsm_timer%stop()

                call pbl_timer%start()
                call pbl(domain, options, real(dt%seconds()))!, halo=1)
                call pbl_timer%stop()
                ! balance u/v and re-calculate dt after winds have been modified by pbl:
                ! if (options%physics%boundarylayer==kPBL_YSU) then
                !     call balance_uvw(   domain%u%data_3d,   domain%v%data_3d,   domain%w%data_3d,       &
                !                         domain%jacobian_u,  domain%jacobian_v,  domain%jacobian_w,      &
                !                         domain%advection_dz, domain%dx, domain%jacobian, options    )
                !
                !     call update_dt(dt, options, domain, end_time)
                !
                !     if ((domain%model_time + dt) > end_time) then
                !         dt = end_time - domain%model_time
                !     endif
                ! endif
                if (options%parameters%debug) call domain_check(domain, "img: "//trim(str(this_image()))//" pbl")

                call convect(domain, options, real(dt%seconds()))!, halo=1)
                if (options%parameters%debug) call domain_check(domain, "img: "//trim(str(this_image()))//" convect")

                
                call ret_timer%start()
                call domain%halo%halo_3d_retrieve_batch(exch_vars=domain%exch_vars, adv_vars=domain%adv_vars)
                !call domain%halo%batch_exch(exch_vars=domain%exch_vars, adv_vars=domain%adv_vars)
                call ret_timer%stop()


                call adv_timer%start()
                call advect(domain, options, real(dt%seconds()))
                if (options%parameters%debug) call domain_check(domain, "img: "//trim(str(this_image()))//" advect(domain", fix=.True.)
                call adv_timer%stop()
                                
                call mp_timer%start()

                call mp(domain, options, real(dt%seconds()))
                if (options%parameters%debug) call domain_check(domain, "img: "//trim(str(this_image()))//" mp_halo", fix=.True.)
                call mp_timer%stop()
                
                if (options%parameters%debug) call domain_check(domain, "img: "//trim(str(this_image()))//" domain%halo_send", fix=.True.)


                !If we are in the last ~10 updates of a time step and a variable drops below 0, we have probably over-shot a value of 0. Force back to 0
                if ((end_time%seconds() - domain%model_time%seconds()) < (dt%seconds()*10)) then

                    call domain%enforce_limits()
                endif


                if (options%parameters%debug) call domain_check(domain, "img: "//trim(str(this_image()))//" domain%apply_forcing", fix=.True.)

            endif
            ! step model_time forward
            domain%model_time = domain%model_time + dt
            
            last_wind_update = last_wind_update + dt%seconds()
        enddo

        !If we overwrote dt to edge up to the next output/input step, revert here
        if (dt_saver%seconds()>0.0) dt = dt_saver
        
    end subroutine step
end module time_step
