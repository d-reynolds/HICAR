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
    use microphysics,               only : mp
    use advection,                  only : advect
    use mod_atm_utilities,          only : exner_function, compute_ivt, compute_iq
    use convection,                 only : convect
    use land_surface,               only : lsm, lsm_apply_fluxes
    use surface_layer,              only : sfc
    use planetary_boundary_layer,   only : pbl, pbl_apply_tend
    use radiation,                  only : rad, rad_apply_dtheta
    use wind,                       only : balance_uvw, update_winds, update_wind_dqdt
    use domain_interface,           only : domain_t
    use boundary_interface,         only : boundary_t
    use options_interface,          only : options_t
    use debug_module,               only : domain_check
    use string,                     only : str
    use time_object,                only : Time_type
    use time_delta_object,          only : time_delta_t
    use icar_constants,             only : STD_OUT_PE
    implicit none
    private
    double precision, parameter  :: DT_BIG = 36000.0
    double precision  :: future_dt_seconds = DT_BIG
    integer :: max_i, max_j, max_k

    public :: step, compute_dt

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
    function compute_dt(dx, u, v, w, rho, dz, ims, ime, kms, kme, jms, jme, its, ite, jts, jte, CFL, use_density) result(dt)
        real,       intent(in)                   :: dx
        real,       intent(in), dimension(ims:ime+1,kms:kme,jms:jme) :: u 
        real,       intent(in), dimension(ims:ime,kms:kme,jms:jme+1) :: v
        real,       intent(in), dimension(ims:ime,kms:kme,jms:jme)   :: w, rho
        real,       intent(in), dimension(kms:kme)     :: dz
        integer,    intent(in)                   :: ims, ime, kms, kme, jms, jme, its, ite, jts, jte
        real,       intent(in)                   :: CFL
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
        max_i = 1
        max_j = 1
        max_k = 1

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
                    if (current_wind > maxwind3d) then
                        max_i = i
                        max_j = j
                        max_k = k
                    endif
                    maxwind3d = max(maxwind3d, current_wind)
                ENDDO
            ENDDO
        ENDDO
                
        maxwind3d = maxwind3d !* sqrt3

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
            write(*,"(A2,f5.1,A,A)") "  ", max(0.0, time_percent)," %  dt=",trim(dt%as_string())
            flush(output_unit)
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
    subroutine update_dt(dt, options, domain)
        implicit none
        type(time_delta_t), intent(inout) :: dt
        type(options_t),    intent(in)    :: options
        type(domain_t),     intent(in)    :: domain

        double precision                  :: present_dt_seconds, seconds_out
        ! compute internal timestep dt to maintain stability
        ! courant condition for 3D advection. 
                
        ! If this is the first step (future_dt_seconds has not yet been set)
        !if (future_dt_seconds == DT_BIG) then
        present_dt_seconds = compute_dt(domain%dx, domain%u%data_3d, domain%v%data_3d, &
                        domain%w%data_3d, domain%density%data_3d, options%domain%dz_levels, &
                        domain%ims, domain%ime, domain%kms, domain%kme, domain%jms, domain%jme, &
                        domain%its, domain%ite, domain%jts, domain%jte, &
                        options%time%cfl_reduction_factor, &
                        use_density=.false.)
        !else
        !    present_dt_seconds = future_dt_seconds
        !endif
        
        future_dt_seconds = compute_dt(domain%dx, domain%u%dqdt_3d, domain%v%dqdt_3d, &
                        domain%w%dqdt_3d, domain%density%data_3d, options%domain%dz_levels, &
                        domain%ims, domain%ime, domain%kms, domain%kme, domain%jms, domain%jme, &
                        domain%its, domain%ite, domain%jts, domain%jte, &
                        options%time%cfl_reduction_factor, &
                        use_density=.false.)
                
        !Minimum dt is min(present_dt_seconds, future_dt_seconds). Then reduce this accross all compute processes
        call MPI_Allreduce(min(present_dt_seconds, future_dt_seconds), seconds_out, 1, MPI_DOUBLE, MPI_MIN, domain%compute_comms)
        
        if (min(present_dt_seconds, future_dt_seconds)==seconds_out) then
            write(*,*) 'time_step determining i:      ',max_i
            write(*,*) 'time_step determining j:      ',max_j
            write(*,*) 'time_step determining k:      ',max_k
            if (future_dt_seconds<present_dt_seconds) then
                write(*,*) 'time_step determining w_grid: ',domain%w%dqdt_3d(max_i,max_k,max_j)
                write(*,*) 'time_step determining u: ',domain%u%dqdt_3d(max_i,max_k,max_j)
                write(*,*) 'time_step determining v: ',domain%v%dqdt_3d(max_i,max_k,max_j)
            else
                write(*,*) 'time_step determining w_grid: ',domain%w%data_3d(max_i,max_k,max_j)
                write(*,*) 'time_step determining u: ',domain%u%data_3d(max_i,max_k,max_j)
                write(*,*) 'time_step determining v: ',domain%v%data_3d(max_i,max_k,max_j)
            endif
            write(*,*) 'time_step determining w_real: ',domain%w_real%data_3d(max_i,max_k,max_j)
            write(*,*) 'time_step determining jaco: ',domain%jacobian(max_i,max_k,max_j)
            write(*,*) 'time_step determining dzdx: ',domain%dzdx%data_3d(max_i,max_k,max_j)
            write(*,*) 'time_step determining dzdy: ',domain%dzdy%data_3d(max_i,max_k,max_j)
            flush(output_unit)
        endif
        ! Set dt to the outcome of reduce
        call dt%set(seconds=seconds_out)
        if (STD_OUT_PE) write(*,*) 'time_step: ',trim(dt%as_string())

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
    subroutine step(domain, forcing, end_time, options)
        implicit none
        type(domain_t),     intent(inout)   :: domain
        type(boundary_t),   intent(inout)   :: forcing
        type(Time_type),    intent(in)      :: end_time
        type(options_t),    intent(in)      :: options

        real            :: last_print_time
        real, save      :: last_wind_update
        logical         :: last_loop

        type(time_delta_t) :: time_step_size, dt_saver, max_dt
        type(Time_type) :: next_input_tmp
        type(time_delta_t), save      :: dt

        last_print_time = 0.0
        call dt_saver%set(seconds=0.0)
        
        time_step_size = end_time - domain%sim_time

        next_input_tmp = domain%next_input - domain%input_dt
        call max_dt%set(seconds=1.0)
        ! Initialize to just over update_dt to force an update on first loop after input ingestion
        if  ( domain%sim_time%equals(next_input_tmp, precision=max_dt) ) then
            last_wind_update = options%wind%update_dt%seconds() + 1
        endif
        
        last_loop = .False.
        ! now just loop over internal timesteps computing all physics in order (operator splitting...)
        do while (domain%sim_time < end_time .and. .not.(last_loop))
            
            !Determine dt
            if (last_wind_update >= options%wind%update_dt%seconds() .or. options%wind%wind_only) then
                call domain%wind_timer%start()
                call update_winds(domain, forcing, options)
                call domain%wind_timer%stop()

                !Now that new winds have been calculated, get new time step in seconds, and see if they require adapting the time step
                ! Note that there will currently be some discrepancy between using the current density and whatever density will be at 
                ! the next time step, but we assume that it is negligable
                ! and that using a CFL criterion < 1.0 will cover this
                call update_dt(dt, options, domain)
                call update_wind_dqdt(domain, options)

                last_wind_update = 0.0

                if (options%wind%wind_only) then
                    domain%sim_time = end_time
                    return
                endif
            endif
            !call update_dt(dt, options, domain, end_time)
            ! Make sure we don't over step the forcing or output period
            if ((domain%sim_time + dt) > end_time) then
                dt_saver = dt
                dt = end_time - domain%sim_time

                ! Sometimes, due to very small time differences, in the inequality controling the loop,
                ! the physics loop can run again. Stop that from happening here
                last_loop = .True.
            endif
            
            ! ! apply/update boundary conditions including internal wind and pressure changes.
            call domain%forcing_timer%start()
            call domain%apply_forcing(forcing,options,real(dt%seconds()))
            call domain%forcing_timer%stop()

            call domain%diagnostic_timer%start()
            call domain%diagnostic_update()
            call domain%diagnostic_timer%stop()

            if (options%adv%advect_density) then
                ! if using advect_density winds need to be balanced at each update
                call domain%wind_bal_timer%start()
                call balance_uvw(domain,options)
                call domain%wind_bal_timer%stop()
            endif

            ! if an interactive run was requested than print status updates everytime at least 5% of the progress has been made
            if (options%general%interactive .and. (STD_OUT_PE)) then
                call print_progress(domain%sim_time, end_time, time_step_size, dt, last_print_time)
            endif
            ! this if is to avoid round off errors causing an additional physics call that won't really do anything

            if (real(dt%seconds()) > 1e-3) then

                if (options%general%debug) call domain_check(domain, "init", fix=.True.)

                call domain%send_timer%start()
                call domain%halo%halo_3d_send_batch(exch_vars=domain%exch_vars, adv_vars=domain%adv_vars)
                call domain%halo%halo_2d_send_batch(exch_vars=domain%exch_vars, adv_vars=domain%adv_vars)
                call domain%send_timer%stop()
    
                call domain%rad_timer%start()
                call rad(domain, options, real(dt%seconds()))
                if (options%general%debug) call domain_check(domain, "rad(domain", fix=.True.)
                call domain%rad_timer%stop()


                call domain%lsm_timer%start()
                call sfc(domain, options, real(dt%seconds()))!, halo=1)
                call lsm(domain, options, real(dt%seconds()))!, halo=1)
                if (options%general%debug) call domain_check(domain, "lsm")
                call domain%lsm_timer%stop()

                call domain%pbl_timer%start()
                call pbl(domain, options, real(dt%seconds()))!, halo=1)
                call domain%pbl_timer%stop()

                call domain%ret_timer%start()
                call domain%halo%halo_3d_retrieve_batch(exch_vars=domain%exch_vars, adv_vars=domain%adv_vars)
                call domain%halo%halo_2d_retrieve_batch(exch_vars=domain%exch_vars, adv_vars=domain%adv_vars)
                call integrate_physics_tendencies(domain, options, real(dt%seconds()))
                !call domain%halo%batch_exch(exch_vars=domain%exch_vars, adv_vars=domain%adv_vars)
                call domain%ret_timer%stop()

                if (options%general%debug) call domain_check(domain, "pbl")

                call convect(domain, options, real(dt%seconds()))!, halo=1)
                if (options%general%debug) call domain_check(domain, "convect")

                

                call domain%adv_timer%start()
                call advect(domain, options, real(dt%seconds()),domain%flux_timer, domain%flux_up_timer, domain%flux_corr_timer, domain%sum_timer, domain%adv_wind_timer)
                if (options%general%debug) call domain_check(domain, "advect(domain", fix=.True.)
                call domain%adv_timer%stop()

                
                                
                call domain%mp_timer%start()

                call mp(domain, options, real(dt%seconds()))
                if (options%general%debug) call domain_check(domain, "mp_halo", fix=.True.)
                call domain%mp_timer%stop()
                
                if (options%general%debug) call domain_check(domain, "domain%halo_send", fix=.True.)


                !If we are in the last ~10 updates of a time step and a variable drops below 0, we have probably over-shot a value of 0. Force back to 0
                if ((end_time%seconds() - domain%sim_time%seconds()) < (dt%seconds()*10)) then
                    call domain%enforce_limits()
                endif


                if (options%general%debug) call domain_check(domain, "domain%apply_forcing", fix=.True.)

            endif
            ! step model_time forward
            call domain%increment_sim_time(dt)
            
            last_wind_update = last_wind_update + dt%seconds()
        enddo

        !If we overwrote dt to edge up to the next output/input step, revert here
        if (dt_saver%seconds()>0.0) dt = dt_saver
        
    end subroutine step


    subroutine integrate_physics_tendencies(domain, options, dt)
        implicit none
        type(domain_t),     intent(inout)   :: domain
        type(options_t),    intent(in)      :: options
        real,               intent(in)      :: dt

        call rad_apply_dtheta(domain, options, dt)
        call lsm_apply_fluxes(domain,options,dt)
        call pbl_apply_tend(domain,options,dt)

    end subroutine integrate_physics_tendencies


end module time_step
