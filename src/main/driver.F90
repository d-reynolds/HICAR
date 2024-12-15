!>-----------------------------------------
!! Main Program
!!
!! Initialize options and memory in init_model
!! Read initial conditions in bc_init (from a restart file if requested)
!! initialize physics packages in init_physics (e.g. tiedke and thompson if used)
!! If this run is a restart run, then set start to the restart timestep
!!      in otherwords, ntimesteps is the number of BC updates from the beginning of the entire model
!!      run, not just from the begining of this restart run
!! calculate model time in seconds based on the time between BC updates (in_dt)
!! Calculate the next model output time from current model time + output time delta (out_dt)
!!
!! Finally, loop until ntimesteps are reached updating boundary conditions and stepping the model forward
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!-----------------------------------------
program icar
    use iso_fortran_env
    use mpi_f08
    use options_interface,  only : options_t
    use domain_interface,   only : domain_t
    use boundary_interface, only : boundary_t
    use output_interface,   only : output_t
    use time_step,          only : step                ! Advance the model forward in time
    use initialization,     only : split_processes, welcome_message, init_options, init_model, init_physics, init_model_state
    use nest_manager,       only : start_nest_context, end_nest_context
    use timer_interface,    only : timer_t
    use time_object,        only : Time_type
    use time_delta_object,  only : time_delta_t
    use icar_constants
    use wind_iterative,     only : finalize_iter_winds
    use wind_iterative_old, only : finalize_iter_winds_old
    use ioserver_interface, only : ioserver_t
    use ioclient_interface, only : ioclient_t
    use io_routines,        only : io_write
    use namelist_utils,     only : get_nml_var_default
    use land_surface,               only : lsm_init

    implicit none

    type(options_t), allocatable :: options(:)
    type(domain_t), allocatable  :: domain(:)
    type(boundary_t), allocatable :: boundary(:)
    type(ioserver_t), allocatable  :: ioserver(:)
    type(ioclient_t), allocatable  :: ioclient(:)

    type(timer_t), allocatable, dimension(:) :: initialization_timer, total_timer, input_timer, &
                                                output_timer, physics_timer, wind_timer, mp_timer, &
                                                adv_timer, rad_timer, lsm_timer, pbl_timer, exch_timer, &
                                                send_timer, ret_timer, wait_timer, forcing_timer, diagnostic_timer, wind_bal_timer

    type(Time_type), allocatable, dimension(:) :: next_output, next_input
    type(Time_type) :: end_of_nest_loop
    type(time_delta_t) :: small_time_delta
    
    integer :: i, ierr, exec_team, n_nests, n, old_nest
    real :: t_val, t_val2, t_val3
    logical :: init_flag, new_input, info_only, gen_nml, only_namelist_check
    logical :: start_time_match = .False.
    character(len=kMAX_FILE_LENGTH) :: namelist_file


    !Initialize MPI if needed
    init_flag = .False.
    call MPI_initialized(init_flag, ierr)
    if (.not.(init_flag)) then
        call MPI_INIT(ierr)
        init_flag = .True.
    endif

    call MPI_Comm_Rank(MPI_COMM_WORLD,PE_RANK_GLOBAL)
    STD_OUT_PE = (PE_RANK_GLOBAL==0)

    call small_time_delta%set(1)

    !-----------------------------------------
    !  Model Initialization
    
    ! Read command line options to determine what kind of run this is
    call read_co(namelist_file, info_only, gen_nml, only_namelist_check)

    if (STD_OUT_PE .and. .not.(gen_nml .or. only_namelist_check .or. info_only)) then
        call welcome_message()
        flush(output_unit)
        write(*,'(/ A)') "--------------------------------------------------------"
        write(*,'(A)')   "Initializing Options"
        write(*,'(A /)') "--------------------------------------------------------"

    endif

    ! Reads user supplied model options
    call init_options(options, namelist_file, info_only=info_only, gen_nml=gen_nml, only_namelist_check=only_namelist_check)
    if (STD_OUT_PE) flush(output_unit)
    if (STD_OUT_PE) write(*,'(/ A)') "--------------------------------------------------------"
    if (STD_OUT_PE) write(*,'(A)')   "Finished reading options, beginning processor assignment"
    if (STD_OUT_PE) write(*,'(A /)') "--------------------------------------------------------"

    n_nests = options(1)%general%nests
    ! Allocate the multiple domains, boundarys, and timers
    allocate(domain(n_nests))
    allocate(boundary(n_nests))
    allocate(ioclient(n_nests))
    allocate(ioserver(n_nests))
    allocate(next_input(n_nests))
    allocate(next_output(n_nests))
    allocate(total_timer(n_nests), initialization_timer(n_nests), input_timer(n_nests), &
             output_timer(n_nests), physics_timer(n_nests), wind_timer(n_nests), &
             mp_timer(n_nests), adv_timer(n_nests), rad_timer(n_nests), &
             lsm_timer(n_nests), pbl_timer(n_nests), exch_timer(n_nests), &
             send_timer(n_nests), ret_timer(n_nests), wait_timer(n_nests), &
             forcing_timer(n_nests), diagnostic_timer(n_nests), wind_bal_timer(n_nests))

    !Determine split of processes which will become I/O servers and which will be compute tasks
    !Also sets constants for the program to keep track of this splitting
    call split_processes(exec_team, domain, ioserver, ioclient)
    if (STD_OUT_PE) write(*,'(/ A)') "--------------------------------------------------------------"
    if (STD_OUT_PE) write(*,'(A)')   "Finished processor assignment, beginning domain initialization"
    if (STD_OUT_PE) write(*,'(A)')   "--------------------------------------------------------------"

    select case(exec_team)
    case(kCOMPUTE_TEAM)
        do i = 1, n_nests
            call total_timer(i)%start()
            call initialization_timer(i)%start()    

            if (STD_OUT_PE) write(*,"(/ A22,I2,A1,A,A16)") "-------------- Domain ",i," (",trim(options(i)%domain%init_conditions_file),") --------------"
            call init_model(options, domain(i), boundary(i), ioclient(i), i)

            call total_timer(i)%stop()
            call initialization_timer(i)%stop()
        enddo
        if (STD_OUT_PE) write(*,'(/ A)') "----------------------------------------------------------------"
        if (STD_OUT_PE) write(*,'(A)')   "Finished domain initialization, beginning physics initialization"
        if (STD_OUT_PE) write(*,'(A)')   "----------------------------------------------------------------"
        old_nest = 1

        ! Need to break the loop here to ensure that the boundary object is first initilaized for all nests
        do i = 1, n_nests
            if (options(i)%general%parent_nest > 0) then
                start_time_match = options(i)%general%start_time == options(options(i)%general%parent_nest)%general%start_time
            endif
            if (options(i)%general%parent_nest == 0 .or. options(i)%restart%restart .or. start_time_match) then
                call total_timer(i)%start()
                call initialization_timer(i)%start()    

                if (old_nest /= i) then
                    call end_nest_context(domain(old_nest), options(old_nest))
                    old_nest = i
                endif

                if (STD_OUT_PE) write(*,"(/ A22,I2,A2,A,A16)") "-------------- Domain ",i," (",trim(options(i)%domain%init_conditions_file),") --------------"
                call init_model_state(options(i), domain(i), boundary(i), ioclient(i))
                next_output(i) = options(i)%general%start_time + options(i)%output%output_dt
                next_input(i) = options(i)%general%start_time !+ options(1)%forcing%input_dt
                old_nest = i

                call total_timer(i)%stop()
                call initialization_timer(i)%stop() 
            end if
        end do

        !  End Compute-Side Initialization
        !--------------------------------------------------------
        if (STD_OUT_PE) write(*,'(/ A)') "------------------------------------------------------"
        if (STD_OUT_PE) write(*,'(A)')   "Initialization complete, beginning physics integration"
        if (STD_OUT_PE) write(*,'(A)')   "------------------------------------------------------"

        do while (ANY(domain%ended .eqv. .False.))
            do i = 1, n_nests

                ! -----------------------------------------------------------
                !  Initialization of child nests, or cycling a completed nest
                ! -----------------------------------------------------------
                ! If we are a child nest, and the domain has not yet been initialized with data...
                if (.not.(domain(i)%started)) then
                    ! ... and if we will be running on the next iteration...
                    if (options(i)%general%start_time - small_time_delta <= domain(options(i)%general%parent_nest)%model_time) then
                        ! ...then initialize ourselves here
                        call total_timer(i)%start()
                        call initialization_timer(i)%start()

                        call end_nest_context(domain(old_nest), options(old_nest))
                        old_nest = i
    
                        if (STD_OUT_PE) write(*,"(/ A22,I2,A2,A,A16)") "-------------- Domain ",i," (",trim(options(i)%domain%init_conditions_file),") --------------"
                        if (STD_OUT_PE) write(*,'(/ A)') "---------------------"
                        if (STD_OUT_PE) write(*,'(A)')   "Waking Domain"
                        if (STD_OUT_PE) write(*,'(A)')   "---------------------"
                
                
                        call init_model_state(options(i), domain(i), boundary(i), ioclient(i))
                        next_output(i) = options(i)%general%start_time + options(i)%output%output_dt
                        next_input(i) = options(i)%general%start_time !+ options(1)%forcing%input_dt

                        call total_timer(i)%stop()
                        call initialization_timer(i)%stop()
                        cycle
                    endif
                endif

                if (domain(i)%ended .or. (domain(i)%started .eqv. .False.)) cycle

                call total_timer(i)%start()

                !Switch nest contexts if needed
                if (old_nest /= i) then
                    call end_nest_context(domain(old_nest), options(old_nest))
                    call start_nest_context(domain(i), options(i))
                    old_nest = i
                endif


                ! -----------------------------------------------------
                !
                !  Read input data
                !
                ! -----------------------------------------------------
                if (STD_OUT_PE) write(*,"(/ A22,I2,A2,A,A16)") "-------------- Domain ",i," (",trim(options(i)%domain%init_conditions_file),") --------------"
                if (STD_OUT_PE) write(*,*) "Updating Boundary conditions"
                next_input(i) = next_input(i) + options(i)%forcing%input_dt
                new_input = .True.
                call input_timer(i)%start()

                call ioclient(i)%receive(boundary(i))
                
                ! after reading all variables that can be read, not compute any remaining variables (e.g. z from p+ps)
                call boundary(i)%update_computed_vars(options(i), update=.True.)
                call domain(i)%interpolate_forcing(boundary(i), update=.True.)

                ! Make the boundary condition dXdt values into units of [X]/s
                call boundary(i)%update_delta_fields(next_input(i) - domain(i)%model_time)
                call domain(i)%update_delta_fields(next_input(i) - domain(i)%model_time)

                call input_timer(i)%stop()

                end_of_nest_loop = step_end(next_input(i), options(i)%general%end_time)

                do while (domain(i)%model_time + small_time_delta < end_of_nest_loop)
                    ! -----------------------------------------------------
                    !
                    !  Integrate physics forward in time
                    !
                    ! -----------------------------------------------------
                    if (STD_OUT_PE) write(*,*) "Running Physics"
                    if (STD_OUT_PE) write(*,*) "  Model time = ", trim(domain(i)%model_time%as_string())
                    if (STD_OUT_PE) write(*,*) "   End  time = ", trim(options(i)%general%end_time%as_string())
                    if (STD_OUT_PE) write(*,*) "  Next Input = ", trim(next_input(i)%as_string())
                    if (STD_OUT_PE) write(*,*) "  Next Output= ", trim(next_output(i)%as_string())
                    if (STD_OUT_PE) flush(output_unit)
                    
                    ! this is the meat of the model physics, run all the physics for the current time step looping over internal timesteps
                    if (.not.(options(i)%wind%wind_only)) then
                        call physics_timer(i)%start()
                        call step(domain(i), boundary(i), step_end(end_of_nest_loop, next_output(i)), new_input, options(i),           &
                                        mp_timer(i), adv_timer(i), rad_timer(i), lsm_timer(i), pbl_timer(i), exch_timer(i), &
                                        send_timer(i), ret_timer(i), wait_timer(i), forcing_timer(i), diagnostic_timer(i), wind_bal_timer(i), wind_timer(i))
                        call physics_timer(i)%stop()
                    elseif (options(i)%wind%wind_only) then
                        call domain(i)%apply_forcing(boundary(i), options(i), real(options(i)%output%output_dt%seconds()))
                        domain(i)%model_time = next_output(i)
                    endif
                    new_input = .False.
                    ! -----------------------------------------------------
                    !  Write output data if it is time
                    ! -----------------------------------------------------
                    if ((domain(i)%model_time + small_time_delta) >= next_output(i)) then
                        if (STD_OUT_PE) write(*,*) "Writing output file"
                        call output_timer(i)%start()
                        call ioclient(i)%push(domain(i))
                        next_output(i) = next_output(i) + options(i)%output%output_dt
                        call output_timer(i)%stop()
                    endif
                enddo
                ! If we have children, who are not yet done, then we need to push our data to them
                if (size(options(i)%general%child_nests) > 0 .and. &
                    ANY(domain(options(i)%general%child_nests)%ended .eqv. .False.)) then
                        call ioclient(i)%update_nest(domain(i))
                endif

                call total_timer(i)%stop()

                if(domain(i)%model_time + small_time_delta > options(i)%general%end_time) then
                    domain(i)%ended = .True.
                    if (STD_OUT_PE) write(*,*) "Domain ",i," has reached the end of its run time."
                endif
            enddo
        end do
        
        if (ANY(options(:)%physics%windtype==kITERATIVE_WINDS) .or. ANY(options(:)%physics%windtype==kLINEAR_ITERATIVE_WINDS)) call finalize_iter_winds() 
        t_val = 0
        do i = 1,n_nests
            t_val = t_val + timer_mean(total_timer(i), domain(1)%compute_comms)
        enddo

        !
        !-----------------------------------------
        if (STD_OUT_PE) then
            call MPI_Comm_Size(MPI_COMM_WORLD,i)


            if (STD_OUT_PE) write(*,'(/ A)') "------------------------------------------------------"
            if (STD_OUT_PE) write(*,'(A)')   "Simulation completed successfully!"
            if (STD_OUT_PE) write(*,'(A /)') "------------------------------------------------------"
            write(*,*) "Model run from : ",trim(options(1)%general%start_time%as_string())
            write(*,*) "           to  : ",trim(options(1)%general%end_time%as_string())
            write(*,*) "Number of images:",i
            write(*,*) ""
            write(*,*) "Timing across all compute images:"
            write(*,*) ""
            write(*,*) "  Total time: ",t_val
            write(*,*) ""
        endif
        do i = 1, n_nests

            if (STD_OUT_PE) write(*,"(A22,I2,A2,A,A16)") "-------------- Domain ",i," (",trim(options(i)%domain%init_conditions_file),") --------------"
            if (STD_OUT_PE) write(*,'(A31 A10 A3 A10 A3 A10)') " ", "mean", " | ", "min", " | ", "max"
            t_val = timer_mean(total_timer(i), domain(1)%compute_comms)
            t_val2 = timer_min(total_timer(i), domain(1)%compute_comms)
            t_val3 = timer_max(total_timer(i), domain(1)%compute_comms)
            if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "total", ":", t_val, " | ", t_val2, " | ", t_val3
            t_val = timer_mean(initialization_timer(i), domain(1)%compute_comms)
            t_val2 = timer_min(initialization_timer(i), domain(1)%compute_comms)
            t_val3 = timer_max(initialization_timer(i), domain(1)%compute_comms)
            if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "init", ":", t_val, " | ", t_val2, " | ", t_val3
            t_val = timer_mean(input_timer(i), domain(1)%compute_comms)
            t_val2 = timer_min(input_timer(i), domain(1)%compute_comms)
            t_val3 = timer_max(input_timer(i), domain(1)%compute_comms)
            if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "input", ":", t_val, " | ", t_val2, " | ", t_val3
            t_val = timer_mean(output_timer(i), domain(1)%compute_comms)
            t_val2 = timer_min(output_timer(i), domain(1)%compute_comms)
            t_val3 = timer_max(output_timer(i), domain(1)%compute_comms)
            if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "output", ":", t_val, " | ", t_val2, " | ", t_val3
            t_val = timer_mean(physics_timer(i), domain(1)%compute_comms)
            t_val2 = timer_min(physics_timer(i), domain(1)%compute_comms)
            t_val3 = timer_max(physics_timer(i), domain(1)%compute_comms)
            if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "physics", ":", t_val, " | ", t_val2, " | ", t_val3
            t_val = timer_mean(mp_timer(i), domain(1)%compute_comms)
            t_val2 = timer_min(mp_timer(i), domain(1)%compute_comms)
            t_val3 = timer_max(mp_timer(i), domain(1)%compute_comms)
            if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "microphysics", ":", t_val, " | ", t_val2, " | ", t_val3
            t_val = timer_mean(adv_timer(i), domain(1)%compute_comms)
            t_val2 = timer_min(adv_timer(i), domain(1)%compute_comms)
            t_val3 = timer_max(adv_timer(i), domain(1)%compute_comms)
            if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "advection", ":", t_val, " | ", t_val2, " | ", t_val3
            t_val = timer_mean(rad_timer(i), domain(1)%compute_comms)
            t_val2 = timer_min(rad_timer(i), domain(1)%compute_comms)
            t_val3 = timer_max(rad_timer(i), domain(1)%compute_comms)
            if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "radiation", ":", t_val, " | ", t_val2, " | ", t_val3
            t_val = timer_mean(lsm_timer(i), domain(1)%compute_comms)
            t_val2 = timer_min(lsm_timer(i), domain(1)%compute_comms)
            t_val3 = timer_max(lsm_timer(i), domain(1)%compute_comms)
            if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "LSM", ":", t_val, " | ", t_val2, " | ", t_val3
            t_val = timer_mean(pbl_timer(i), domain(1)%compute_comms)
            t_val2 = timer_min(pbl_timer(i), domain(1)%compute_comms)
            t_val3 = timer_max(pbl_timer(i), domain(1)%compute_comms)
            if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "PBL", ":", t_val, " | ", t_val2, " | ", t_val3
            t_val = timer_mean(forcing_timer(i), domain(1)%compute_comms)
            t_val2 = timer_min(forcing_timer(i), domain(1)%compute_comms)
            t_val3 = timer_max(forcing_timer(i), domain(1)%compute_comms)
            if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "forcing", ":", t_val, " | ", t_val2, " | ", t_val3
            t_val = timer_mean(wind_bal_timer(i), domain(1)%compute_comms)
            t_val2 = timer_min(wind_bal_timer(i), domain(1)%compute_comms)
            t_val3 = timer_max(wind_bal_timer(i), domain(1)%compute_comms)
            if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "wind bal", ":", t_val, " | ", t_val2, " | ", t_val3
            t_val = timer_mean(diagnostic_timer(i), domain(1)%compute_comms)
            t_val2 = timer_min(diagnostic_timer(i), domain(1)%compute_comms)
            t_val3 = timer_max(diagnostic_timer(i), domain(1)%compute_comms)
            if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "diagnostic", ":", t_val, " | ", t_val2, " | ", t_val3
            t_val = timer_mean(send_timer(i), domain(1)%compute_comms)
            t_val2 = timer_min(send_timer(i), domain(1)%compute_comms)
            t_val3 = timer_max(send_timer(i), domain(1)%compute_comms)
            if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "halo-exchange(send)", ":", t_val, " | ", t_val2, " | ", t_val3
            t_val = timer_mean(ret_timer(i), domain(1)%compute_comms)
            t_val2 = timer_min(ret_timer(i), domain(1)%compute_comms)
            t_val3 = timer_max(ret_timer(i), domain(1)%compute_comms)
            if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "halo-exchange(retrieve)", ":", t_val, " | ", t_val2, " | ", t_val3
            t_val = timer_mean(wait_timer(i), domain(1)%compute_comms)
            t_val2 = timer_min(wait_timer(i), domain(1)%compute_comms)
            t_val3 = timer_max(wait_timer(i), domain(1)%compute_comms)
            if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "halo-exchange(wait)", ":", t_val, " | ", t_val2, " | ", t_val3
            t_val = timer_mean(wind_timer(i), domain(1)%compute_comms)
            t_val2 = timer_min(wind_timer(i), domain(1)%compute_comms)
            t_val3 = timer_max(wind_timer(i), domain(1)%compute_comms)
            if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "winds", ":", t_val, " | ", t_val2, " | ", t_val3
        enddo
    case(kIO_TEAM)
    
        !This is, unfortunately, stupidly, needed to allow coarrays to work. 
        !Perhaps on later releases with better team support it can be removed?
#ifdef  CRAY_PE
        call domain(1)%init(options(1))
#endif
        do i = 1, n_nests
            call ioserver(i)%init(options,i)
        enddo

        do i = 1, n_nests
            if (options(i)%general%parent_nest > 0) then
                start_time_match = (options(i)%general%start_time == options(options(i)%general%parent_nest)%general%start_time)
            endif
            if (options(i)%general%parent_nest == 0 .or. options(i)%restart%restart .or. start_time_match) then

                !Get initial conditions
                if (options(i)%general%parent_nest == 0) call ioserver(i)%read_file()

                if (options(i)%restart%restart) then
                    call ioserver(i)%read_restart_file(options(i))
                else
                    call ioserver(i)%write_file(ioserver(i)%io_time)
                endif

                if ((size(options(i)%general%child_nests) > 0)) then
                    ! This call will gather the model state of the forcing fields from the nest parent
                    call ioserver(i)%gather_forcing(ioserver(options(i)%general%child_nests))
                endif

                ioserver(i)%started = .True.
            end if
        end do

        ! Batch off another round of file reads so that we are always one step ahead of the compute team
        ! This is needed to ensure that the compute team has the data it needs when it needs it
        do i = 1, n_nests
            if (options(i)%general%parent_nest == 0) call ioserver(i)%read_file()
        end do

        !  End IO-Side Initialization
        !--------------------------------------------------------

        do while (ANY(ioserver(:)%ended .eqv. .False.))
            do i = 1, n_nests

                ! If we are a child nest, not at the end of our run time...
                if (ioserver(i)%started .eqv. .False.) then
                    ! ... and if we will be running on the next iteration...
                    if (ioserver(i)%io_time - small_time_delta <= ioserver(options(i)%general%parent_nest)%io_time) then
                        ! ...then initialize ourselves here
                        call ioserver(i)%write_file(ioserver(i)%io_time)
        
                        if ((size(options(i)%general%child_nests) > 0)) then
                            ! This call will gather the model state of the forcing fields from the nest parent
                            call ioserver(i)%gather_forcing(ioserver(options(i)%general%child_nests))
                        endif
        
                        ioserver(i)%started = .True.
                        cycle
                    endif
                endif

                if (ioserver(i)%ended .or. (ioserver(i)%started .eqv. .False.)) cycle

                next_input(i) = ioserver(i)%io_time + options(i)%forcing%input_dt
                next_output(i) = ioserver(i)%io_time + options(i)%output%output_dt

                end_of_nest_loop = step_end(next_input(i), options(i)%general%end_time)

                do while (ioserver(i)%io_time + small_time_delta < end_of_nest_loop)

                    ioserver(i)%io_time = step_end(next_input(i),next_output(i))

                    !See if we even have files to read
                    if (ioserver(i)%files_to_read) then
                        !See of it is time to read.
                        if (ioserver(i)%io_time + options(i)%forcing%input_dt + small_time_delta >= next_input(i)) then
                            call ioserver(i)%read_file()
                        endif
                    endif

                    ! If the next event is a change is nest scope, then we need to wait on this, as the compute team will wait for us
                    if ( (size(options(i)%general%child_nests) > 0) .and. ioserver(i)%io_time + small_time_delta >= end_of_nest_loop) then
                        if (ANY(ioserver(options(i)%general%child_nests)%ended .eqv. .False.)) then
                            ! This call will gather the model state of the forcing fields from the nest parent
                            call ioserver(i)%gather_forcing(ioserver(options(i)%general%child_nests))
                        endif
                    endif

                    ! If we aren't yet done, then wait for an output
                    ! If we are at an output step, do it now
                    if (ioserver(i)%io_time+small_time_delta >= next_output(i) .and. ioserver(i)%io_time <= options(i)%general%end_time) then
                        call ioserver(i)%write_file(ioserver(i)%io_time)
                        next_output(i) = ioserver(i)%io_time + options(i)%output%output_dt
                    endif
                end do

                if (ioserver(i)%io_time + small_time_delta >= options(i)%general%end_time) then
                    ioserver(i)%ended = .True.
                endif
            enddo
            !if (.not.(ioserver(1)%files_to_read) .and. ioserver(n_nests)%io_time + small_time_delta >= options(1)%general%end_time) io_loop = .False.
        enddo
        !If we are done with the program
        call ioserver(1)%close_files()
    end select

    call MPI_Barrier(MPI_COMM_WORLD)
    CALL MPI_Finalize()

contains

    function timer_mean(timer,comms) result(mean_t)
        implicit none
        type(timer_t), intent(inout) :: timer
        type(MPI_Comm), intent(in) :: comms

        real :: mean_t, t_sum
        integer :: ierr
            
        t_sum = timer%get_time()
        call MPI_Allreduce(MPI_IN_PLACE,t_sum,1,MPI_REAL,MPI_SUM,comms,ierr)
        mean_t = t_sum/kNUM_COMPUTE
    
    end function

    function timer_max(timer,comms) result(t_sum)
        implicit none
        type(timer_t), intent(inout) :: timer
        type(MPI_Comm), intent(in) :: comms

        real :: mean_t, t_sum
        integer :: ierr
            
        t_sum = timer%get_time()
        call MPI_Allreduce(MPI_IN_PLACE,t_sum,1,MPI_REAL,MPI_MAX,comms,ierr)
    
    end function

    function timer_min(timer,comms) result(t_sum)
        implicit none
        type(timer_t), intent(inout) :: timer
        type(MPI_Comm), intent(in) :: comms

        real :: mean_t, t_sum
        integer :: ierr
            
        t_sum = timer%get_time()
        call MPI_Allreduce(MPI_IN_PLACE,t_sum,1,MPI_REAL,MPI_MIN,comms,ierr)
    
    end function


    function step_end(time1, time2, end_time) result(min_time)
        implicit none
        type(Time_type), intent(in) :: time1, time2
        type(Time_type), optional, intent(in) :: end_time
        type(Time_type) :: min_time

        if (time1 <= time2 ) then
            min_time = time1
        else
            min_time = time2
        endif

        if (present(end_time)) then
            if (end_time < min_time) then
                min_time = end_time
            endif
        endif
    end function
    
    subroutine read_co(nml_file, info, gen_nml, only_namelist_check)
        implicit none
        character(len=kMAX_FILE_LENGTH), intent(out) :: nml_file
        logical, intent(out) :: info, gen_nml, only_namelist_check

        integer :: cnt, p
        character(len=kMAX_FILE_LENGTH) :: first_arg, arg, default
        logical :: file_exists

        nml_file = ""
        info = .False.
        gen_nml = .False.
        only_namelist_check = .False.

        cnt = command_argument_count()

        ! If there are no command line arguments, throw error
        if (cnt == 0 .and. STD_OUT_PE) then
            write(*,*) "Usage: ./HICAR [-v [variable_name ...|--all]] [--check-nml] [--gen-nml] namelist_file"
            write(*,*) "    -v [variable_name ...|--all]: Print information about the namelist variable(s) variable_name, ... "
            write(*,*) "                                  --all prints out information for all namelist variables."
            write(*,*) "    --check-nml:                  Check the namelist file for errors without running the model."
            write(*,*) "    --gen-nml:                    Generate a namelist file with default values."
            write(*,*) "    namelist_file:                The name of the namelist file to use."
            write(*,*)
            write(*,*) "    Example to generate a namelist with default values:  ./HICAR --gen-nml namelist_file.nml"
            write(*,*) "    Example to check namelist:                           ./HICAR --check-nml namelist_file.nml"
            write(*,*) "    Example to run model:                                ./HICAR namelist_file.nml"
            write(*,*) "    Example to learn about a namelist variable:          ./HICAR -v mp"
            write(*,*) "    Example to generate namelist variable documentation: ./HICAR -v --all > namelist_doc.txt"

            write(*,*)


            stop
        endif

        ! get first command line argument
        call get_command_argument(1, first_arg)

        ! test if argument is a '-v' type argument, indicating that we should print namelist info for this variable
        if (first_arg == '-v') then
            ! if there is no second argument, throw error
            if (cnt >= 2) then
                call get_command_argument(2, arg)
                if (arg == '--all') then
                    info = .True.
                else
                    do p = 2, cnt
                        call get_command_argument(p, arg)
                        default = get_nml_var_default(arg, info=.True.)
                    end do
                    stop
                endif
            elseif (cnt < 2) then
                if (STD_OUT_PE) write(*,*) "ERROR: No variable name provided with -v flag."
                stop
            endif
        elseif (first_arg == '--check-nml') then
            only_namelist_check = .True.
            if (cnt >= 2) then
                call get_command_argument(2, nml_file)
            elseif (cnt < 2) then
                if (STD_OUT_PE) write(*,*) "ERROR: No namelist provided with the --check-nml flag."
                stop
            endif
        elseif (first_arg == '--gen-nml') then
            gen_nml = .True.
            if (cnt >= 2) then
                call get_command_argument(2, nml_file)
            elseif (cnt < 2) then
                if (STD_OUT_PE) write(*,*) "ERROR: No namelist provided with the --gen-nml flag."
                stop
            endif
        else
            nml_file = first_arg
        endif

        if (.not.first_arg=='-v' .and. .not.first_arg=='--gen-nml') then
            ! Check that the options file actually exists
            INQUIRE(file=trim(nml_file), exist=file_exists)

            ! if options file does not exist, print an error and quit
            if (.not.file_exists) then
                if (STD_OUT_PE) write(*,*) "Using options file = ", trim(nml_file)
                stop "Options file does not exist. "
            endif
        endif
    end subroutine

end program

! This is the Doxygen mainpage documentation.  This should be moved to another file at some point.

!>------------------------------------------
!!  @mainpage
!!
!!  @section Introduction
!!  ICAR is a simplified atmospheric model designed primarily for climate downscaling, atmospheric sensitivity tests,
!!  and hopefully educational uses. At this early stage, the model is still undergoing rapid development, and users
!!  are encouraged to get updates frequently.
!!
!!  @section Running_ICAR
!!  To run the model 3D time-varying atmospheric data are required, though an ideal test case can be generated for
!!  simple simulations as well. There are some sample python scripts to help make input forcing files, but the WRF
!!  pre-processing system can also be used. Low-resolution WRF output files can be used directly, various reanalysis
!!  and GCM output files can be used with minimal pre-processing (just get all the variables in the same netcdf file.)
!!  In addition, a high-resolution netCDF topography file is required. This will define the grid that ICAR will run on.
!!  Finally an ICAR options file is used to specify various parameters for the model. A sample options file is provided
!!  in the run/ directory.
!!
!!  @section Developing
!!  This document provides the primary API and code structure documentation. The code is based on github.com/NCAR/icar
!!  Developers are encouraged to fork the main git repository and maintain their own git repository from which to
!!  issue pull requests.
!!
!!  @section Reference
!!  Gutmann, E. D., I. Barstad, M. P. Clark, J. R. Arnold, and R. M. Rasmussen (2016),
!!  The Intermediate Complexity Atmospheric Research Model, J. Hydrometeor, doi:<a href="http://dx.doi.org/10.1175/JHM-D-15-0155.1">10.1175/JHM-D-15-0155.1</a>.
!!
!!------------------------------------------
