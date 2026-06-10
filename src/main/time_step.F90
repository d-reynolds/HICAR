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
    use iso_fortran_env, only : output_unit
    use mpi_f08, only: MPI_Allreduce, MPI_REAL, MPI_MIN
    use string,                     only : as_string
    use microphysics,               only : mp
    use advection,                  only : advect
    use convection,                 only : convect
    use land_surface,               only : lsm, lsm_apply_fluxes
    use surface_layer,              only : sfc
    use planetary_boundary_layer,   only : pbl, pbl_apply_tend
    use radiation,                  only : rad, rad_apply_dtheta
    use snow_drift,                 only : snow_drift_apply_feedback
    use wind,                       only : balance_uvw, update_winds, update_wind_dqdt
    use domain_interface,           only : domain_t
    use options_interface,          only : options_t
    use debug_module,               only : domain_check
    use time_object,                only : Time_type
    use time_delta_object,          only : time_delta_t
    use icar_constants,             only : STD_OUT_PE, kVARS, kRANS_WINDS
    use snow_drift,                 only : snow_drift_step

    implicit none

    private
    real, parameter  :: DT_BIG = 36000.0
    real  :: future_dt_seconds = DT_BIG
    integer :: max_i, max_j, max_k
    real :: max_u, max_v, max_w
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
    function compute_dt(dx, u, v, w, rho, dz, ims, ime, kms, kme, jms, jme, its, ite, jts, jte, CFL, err_msg) result(dt)
        real,       intent(in)                   :: dx
        real,       intent(in), dimension(ims:ime+1,kms:kme,jms:jme) :: u 
        real,       intent(in), dimension(ims:ime,kms:kme,jms:jme+1) :: v
        real,       intent(in), dimension(ims:ime,kms:kme,jms:jme)   :: w, rho, dz
        integer,    intent(in)                   :: ims, ime, kms, kme, jms, jme, its, ite, jts, jte
        real,       intent(in)                   :: CFL
        character(len=*), intent(in), optional  :: err_msg
        
        ! output value
        real :: dt

        ! locals
        integer :: i, j, k, indx(3)
        real :: maxwind3d, maxwind1d, maxwind, sqrt3
        real, dimension(its:ite, kms:kme, jts:jte) :: CFL_wind3d, CFL_wind1d

        sqrt3 = sqrt(3.0)

        CFL_wind3d = 0
        CFL_wind1d = 0

        max_i = 1
        max_j = 1
        max_k = 1

        ! to ensure we are stable for 3D advection we'll use the average "max" wind speed
        ! but that average has to be divided by sqrt(3) for stability in 3 dimensional advection

        !$acc parallel loop gang vector collapse(3) present(u, v, w, rho, dz, dx) copy(CFL_wind3d, CFL_wind1d)
        do j=jts,jte
            do k=kms,kme
                do i=its,ite
                    CFL_wind3d(i,k,j) = max(abs(u(i,k,j)), abs(u(i+1,k,j))) / dx &
                                    +max(abs(v(i,k,j)), abs(v(i,k,j+1))) / dx &
                                    +max(abs(w(i,k,j)), abs(w(i,max(1,k-1),j))) / dz(i,k,j)
                                    
                    CFL_wind1d(i,k,j) = max(( max( abs(u(i,k,j)), abs(u(i+1,k,j)) ) / dx), &
                                         ( max( abs(v(i,k,j)), abs(v(i,k,j+1)) ) / dx), &
                                         ( max( abs(w(i,k,j)), abs(w(i,max(1,k-1),j)) ) / dz(i,k,j) ))
                ENDDO
            ENDDO
        ENDDO

        maxwind3d = maxval(CFL_wind3d)
        maxwind1d = maxval(CFL_wind1d)

        ! According to the WRF technical documentation, the CFL condition for 3D advection is
        ! obtained by taking the CFL criterion for the 1D case and dividing by sqrt(3).
        ! A more rigorous restriction is that CFL_x + CFL_y + CFL_z < CFL (Baldauf 2008), as calculated above.
        ! This is more restrictive for diagonal winds, but allows the restriction to 
        ! relax to the 1D case for one-dimensional winds. Since the CFL constraint on vertical
        ! winds is often much less than 1, the criterion has a maximum of roughly 2x the 1D case.

        ! In practical experience, the WRF case
        ! is often sufficiently stable. So, to allow for faster time steps, we make this the
        ! maximum limiting wind speed condition.
        maxwind = maxwind3d!min(maxwind3d,maxwind1d*sqrt3)

        ! if (maxwind == maxwind1d*sqrt3) then
        !     !get index of maxval in maxwind1d, breaking into i, j, and k indices
        !     indx = findloc(CFL_wind1d,maxwind1d)
        ! else
        !get index of maxval in maxwind3d, breaking into i, j, and k indices
        indx = findloc(CFL_wind3d,maxwind3d)
        ! endif

        max_i = indx(1)+its-1
        max_k = indx(2)
        max_j = indx(3)+jts-1

        !$acc kernels present(u, v, w) copyin(max_i, max_j, max_k) copy(max_u, max_v, max_w)
        max_u = max(abs(u(max_i,max_k,max_j)), abs(u(max_i+1,max_k,max_j)))
        max_v = max(abs(v(max_i,max_k,max_j)), abs(v(max_i,max_k,max_j+1)))
        max_w = max(abs(w(max_i,max_k,max_j)), abs(w(max_i,max(1,max_k-1),max_j)))
        !$acc end kernels

        dt = CFL / maxwind

        ! If we have too small a time step throw an error
        ! something is probably wrong in the physics or input data
        if (dt<1e-1) then
            if (STD_OUT_PE) then 
                if (present(err_msg)) write(*,*) trim(err_msg)
                write(*,*) "dt   = ", dt
                !$acc update host(u, v, w)
                write(*,*) "Umax = ", maxval(abs(u))
                write(*,*) "Vmax = ", maxval(abs(v))
                write(*,*) "Wmax = ", maxval(abs(w))
            endif
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
            write(*,"(A2,f5.1,A,A)") "  ", max(0.0, time_percent)," %  dt=",trim(as_string(dt))
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
    subroutine update_dt(dt, options, domain, quiet)
        implicit none
        type(time_delta_t), intent(inout) :: dt
        type(options_t),    intent(in)    :: options
        type(domain_t),     intent(in)    :: domain
        logical, optional,  intent(in)    :: quiet

        real                  :: present_dt_seconds, seconds_out
        integer :: ierr
        integer :: max_i_pres, max_j_pres, max_k_pres
        real :: max_u_pres, max_v_pres, max_w_pres
        logical :: verbose

        verbose = .True.
        if (present(quiet)) verbose = .not. quiet
        ! compute internal timestep dt to maintain stability
        ! courant condition for 3D advection. 
                
        ! If this is the first step (future_dt_seconds has not yet been set)
        !if (future_dt_seconds == DT_BIG) then

        present_dt_seconds = compute_dt(domain%dx, domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d, domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d, &
                        domain%vars_3d(domain%var_indx(kVARS%w)%v)%data_3d, domain%vars_3d(domain%var_indx(kVARS%density)%v)%data_3d, &
                        domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d, &
                        domain%ims, domain%ime, domain%kms, domain%kme, domain%jms, domain%jme, domain%its, domain%ite, domain%jts, domain%jte, &
                        options%time%cfl_reduction_factor, &
                        err_msg="error computing dt for winds at the current time step")
        max_u_pres = max_u
        max_v_pres = max_v
        max_w_pres = max_w
        max_i_pres = max_i
        max_j_pres = max_j
        max_k_pres = max_k

        if (options%physics%windtype == kRANS_WINDS) then
            ! Under RANS the dqdt arrays are projection workspace, not the
            ! winds at the next forcing time — only the present prognostic
            ! field constrains dt (it is re-evaluated every step).
            future_dt_seconds = present_dt_seconds
        else
            future_dt_seconds = compute_dt(domain%dx, domain%vars_3d(domain%var_indx(kVARS%u)%v)%dqdt_3d, domain%vars_3d(domain%var_indx(kVARS%v)%v)%dqdt_3d, &
                        domain%vars_3d(domain%var_indx(kVARS%w)%v)%dqdt_3d, domain%vars_3d(domain%var_indx(kVARS%density)%v)%data_3d, &
                        domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d, &
                        domain%ims, domain%ime, domain%kms, domain%kme, domain%jms, domain%jme, domain%its, domain%ite, domain%jts, domain%jte, &
                        options%time%cfl_reduction_factor, &
                        err_msg="error computing dt for winds at the future input timestep")
        endif

        !Minimum dt is min(present_dt_seconds, future_dt_seconds). Then reduce this accross all compute processes
        call MPI_Allreduce(min(present_dt_seconds, future_dt_seconds), seconds_out, 1, MPI_REAL, MPI_MIN, domain%compute_comms)
        
        if (min(present_dt_seconds, future_dt_seconds)==seconds_out .and. verbose) then

            if (future_dt_seconds>present_dt_seconds) then
                max_u = max_u_pres
                max_v = max_v_pres
                max_w = max_w_pres
                max_i = max_i_pres
                max_j = max_j_pres
                max_k = max_k_pres
            endif
            write(*,*) 'time_step determining i:      ',max_i
            write(*,*) 'time_step determining j:      ',max_j
            write(*,*) 'time_step determining k:      ',max_k

            write(*,*) 'time_step determining w_grid: ',max_w
            write(*,*) 'time_step determining u: ',max_u
            write(*,*) 'time_step determining v: ',max_v

            write(*,*) 'time_step determining w_real: ',domain%vars_3d(domain%var_indx(kVARS%w_real)%v)%data_3d(max_i,max_k,max_j)
            write(*,*) 'time_step determining jaco: ',domain%vars_3d(domain%var_indx(kVARS%jacobian)%v)%data_3d(max_i,max_k,max_j)
            write(*,*) 'time_step determining dzdx: ',domain%vars_3d(domain%var_indx(kVARS%dzdx)%v)%data_3d(max_i,max_k,max_j)
            write(*,*) 'time_step determining dzdy: ',domain%vars_3d(domain%var_indx(kVARS%dzdy)%v)%data_3d(max_i,max_k,max_j)
            flush(output_unit)
        endif
        ! Set dt to the outcome of reduce
        call dt%set(seconds=seconds_out)
        if (STD_OUT_PE .and. verbose) write(*,*) 'time_step: ',trim(as_string(dt))
        if (STD_OUT_PE .and. verbose) flush(output_unit)

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
    subroutine step(domain, end_time, options)
        implicit none
        type(domain_t),     intent(inout)   :: domain
        type(Time_type),    intent(in)      :: end_time
        type(options_t),    intent(in)      :: options

        real            :: last_print_time
        real, save      :: last_wind_update
        real            :: remaining_step_seconds
        logical         :: last_loop, force_update_winds

        type(time_delta_t) :: time_step_size, max_dt
        type(Time_type) :: tmp_time
        type(time_delta_t), save      :: dt

        last_print_time = 0.0
        call max_dt%set(seconds=1.0)

        time_step_size = end_time - domain%sim_time

        force_update_winds = domain%sim_time%equals(options%general%start_time, precision=max_dt)
        if (options%restart%restart) force_update_winds = domain%sim_time%equals(options%restart%restart_time, precision=max_dt)

        call tmp_time%set(domain%next_input%mjd() - domain%input_dt%days())

        force_update_winds = (force_update_winds .or. domain%sim_time%equals(tmp_time, precision=max_dt) )
        
        last_loop = .False.

        ! now just loop over internal timesteps computing all physics in order (operator splitting...)
        do while (domain%sim_time < end_time .and. .not.(last_loop))
            
            !Determine dt
            remaining_step_seconds = real(end_time%seconds() - domain%sim_time%seconds())

            if (force_update_winds .or. &
                (last_wind_update >= options%wind%update_dt%seconds() .and. &
                 remaining_step_seconds > domain%dt) .or. &
                options%wind%wind_only) then

                call domain%wind_timer%start()
                call domain%halo%exch_var(domain%vars_3d(domain%var_indx(kVARS%density)%v),corners=.True.)
                call update_winds(domain, options)
                call domain%wind_timer%stop()

                !Now that new winds have been calculated, get new time step in seconds, and see if they require adapting the time step
                ! Note that there will currently be some discrepancy between using the current density and whatever density will be at 
                ! the next time step, but we assume that it is negligable
                ! and that using a CFL criterion < 1.0 will cover this
                if (domain%restart_dt > 0.0) then
                    call dt%set(seconds=domain%restart_dt)
                    domain%restart_dt = 0.0
                else
                    call update_dt(dt, options, domain)
                endif
                domain%dt = real(dt%seconds())

                ! Under RANS the dqdt arrays are projection workspace and the
                ! boundary nudging reads the forcing-side arrays directly, so
                ! no wind tendency conversion is needed (or meaningful).
                if (options%physics%windtype /= kRANS_WINDS) then
                    call update_wind_dqdt(domain, real(options%wind%update_dt%seconds()) - domain%forcing_elapsed)
                endif

                ! Reset forcing_elapsed after first wind update so subsequent wind updates
                ! within the same step interval use the full wind_update_dt
                domain%forcing_elapsed = 0.0
                last_wind_update = 0.0
                force_update_winds = .False.

                if (options%wind%wind_only) then
                    domain%sim_time = end_time
                    exit
                endif
            endif

            ! Under RANS the winds evolve every physics step, so the CFL
            ! limit must be re-evaluated each step, not only at wind updates.
            if (options%physics%windtype == kRANS_WINDS) then
                call update_dt(dt, options, domain, quiet=.True.)
                domain%dt = real(dt%seconds())
            endif

            ! Make sure we don't over step the forcing or output period
            call tmp_time%set(domain%sim_time%mjd() + dt%days())

            if (tmp_time > end_time) then
                dt = end_time - domain%sim_time

                ! Sometimes, due to very small time differences, in the inequality controling the loop,
                ! the physics loop can run again. Stop that from happening here
                last_loop = .True.
            endif
            
            ! ! apply/update boundary conditions including internal wind and pressure changes.
            call domain%forcing_timer%start()
            call domain%apply_forcing(options,real(dt%seconds()))
            call domain%forcing_timer%stop()

            call domain%diagnostic_timer%start()
            call domain%diagnostic_update()
            call domain%diagnostic_timer%stop()


            ! if an interactive run was requested than print status updates everytime at least 5% of the progress has been made
            if (options%general%interactive .and. (STD_OUT_PE)) then
                call print_progress(domain%sim_time, end_time, time_step_size, dt, last_print_time)
            endif
            ! this if is to avoid round off errors causing an additional physics call that won't really do anything

            if (real(dt%seconds()) > 1e-3) then

                ! if (options%general%debug) call domain_check(domain, "init", fix=.True.)

                call domain%send_timer%start()
                call domain%halo_3d_send()
                call domain%halo_2d_send()
                call domain%send_timer%stop()

                call domain%rad_timer%start()
                call rad(domain, options, real(dt%seconds()))
                if (options%general%debug) call domain_check(domain, "rad")
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
                call domain%halo_3d_retrieve()
                call domain%halo_2d_retrieve()
                call domain%ret_timer%stop()

                ! balance_uvw is the diagnostic mass-conservation enforcer; under
                ! kRANS_WINDS the projection step inside advect() already enforces
                ! ∇·(ρ₀ u) = 0, and calling balance_uvw on top would discard that work
                ! by re-dividing w to compensate for an already-zero divergence.
                if (options%adv%advect_density .and. options%physics%windtype /= kRANS_WINDS) then
                ! if using advect_density winds need to be balanced at each update
                call domain%wind_bal_timer%start()
                call balance_uvw(domain,options%adv%advect_density)
                call domain%wind_bal_timer%stop()
                endif

                if (options%general%debug) call domain_check(domain, "pbl")
                call convect(domain, options, real(dt%seconds()))!, halo=1)
                if (options%general%debug) call domain_check(domain, "convect")
                

                call domain%adv_timer%start()
                call advect(domain, options, real(dt%seconds()),domain%flux_timer, domain%flux_corr_timer, domain%sum_timer, domain%adv_wind_timer)
                !call domain%enforce_limits()
                if (options%general%debug) call domain_check(domain, "advect")
                call domain%adv_timer%stop()

                call domain%lsm_timer%start()
                call snow_drift_step(domain, options, real(dt%seconds()))!, halo=1)
                if (options%general%debug) call domain_check(domain, "suspension")
                call domain%lsm_timer%stop()

                call integrate_physics_tendencies(domain, options, real(dt%seconds()))

                                
                call domain%mp_timer%start()
                call mp(domain, options, real(dt%seconds()))
                if (options%general%debug) call domain_check(domain, "mp_halo")
                call domain%mp_timer%stop()
                

            endif
            ! step model_time forward
            call domain%increment_sim_time(dt)
            last_wind_update = last_wind_update + dt%seconds()
        enddo

        !If we overwrote dt to edge up to the next output/input step, revert here
        call dt%set(seconds=domain%dt)
        
    end subroutine step


    subroutine integrate_physics_tendencies(domain, options, dt)
        implicit none
        type(domain_t),     intent(inout)   :: domain
        type(options_t),    intent(in)      :: options
        real,               intent(in)      :: dt

        call rad_apply_dtheta(domain, options, dt)
        call lsm_apply_fluxes(domain,options,dt)
        call pbl_apply_tend(domain,options,dt)
        call snow_drift_apply_feedback(domain, options, dt)

    end subroutine integrate_physics_tendencies


end module time_step
