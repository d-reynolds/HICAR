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
    use flow_object_interface, only : flow_obj_t
    use domain_interface,   only : domain_t
    use boundary_interface, only : boundary_t
    use output_interface,   only : output_t
    use time_step,          only : step                ! Advance the model forward in time
    use initialization,     only : split_processes, welcome_message, init_options, init_model, init_physics
    use nest_manager,      only : nest_next_up, should_update_nests, can_update_child_nest, all_nests_not_done, &
                                end_nest_context, wake_nest, switch_nest_context
    use timer_interface,    only : timer_t
    use time_object,        only : Time_type
    use time_delta_object,  only : time_delta_t
    use icar_constants
    use ioserver_interface, only : ioserver_t
    use ioclient_interface, only : ioclient_t
    use io_routines,        only : io_write
    use namelist_utils,     only : get_nml_var_default
    use land_surface,               only : lsm_init
    use output_metadata,    only : list_output_vars
    implicit none

    type(options_t), allocatable :: options(:)
    ! type(domain_t) :: domain(kMAX_NESTS) ! Currently hard-coded, could be dynamic, but compile time on 
    !                                      ! the Cray compiler is very slow for dynamically allocating large, derrived type arrays
    type(boundary_t), allocatable :: boundary(:)
    class(flow_obj_t), allocatable  :: components(:)
    type(ioclient_t), allocatable  :: ioclient(:)
    
    integer :: i, ierr, exec_team, n_nests, n, old_nest, PE_RANK_GLOBAL
    real :: t_val, t_val2, t_val3
    logical :: init_flag, info_only, gen_nml, only_namelist_check
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

    !-----------------------------------------
    !  Model Initialization
    
    ! Read command line options to determine what kind of run this is
    call read_co(namelist_file, info_only, gen_nml, only_namelist_check)

    if (STD_OUT_PE .and. .not.(gen_nml .or. only_namelist_check .or. info_only)) then
        call welcome_message()
        write(*,'(/ A)') "--------------------------------------------------------"
        write(*,'(A)')   "Initializing Options"
        write(*,'(A /)') "--------------------------------------------------------"
        flush(output_unit)

    endif

    ! Reads user supplied model options
    call init_options(options, namelist_file, info_only=info_only, gen_nml=gen_nml, only_namelist_check=only_namelist_check)

    if (STD_OUT_PE) flush(output_unit)
    if (STD_OUT_PE) write(*,'(/ A)') "--------------------------------------------------------"
    if (STD_OUT_PE) write(*,'(A)')   "Finished reading options, beginning processor assignment"
    if (STD_OUT_PE) write(*,'(A /)') "--------------------------------------------------------"

    n_nests = options(1)%general%nests

    ! !Allocate the multiple domains, boundarys
    allocate(boundary(n_nests))
    allocate(ioclient(n_nests))

    !Determine split of processes which will become I/O servers and which will be compute tasks
    !Also sets constants for the program to keep track of this splitting
    call split_processes(components, ioclient, n_nests)
    if (STD_OUT_PE) write(*,'(/ A)') "--------------------------------------------------------------"
    if (STD_OUT_PE) write(*,'(A)')   "Finished processor assignment, beginning domain initialization"
    if (STD_OUT_PE) write(*,'(A)')   "--------------------------------------------------------------"

    do i = 1, n_nests
        call init_component(components(i), options, boundary(i), ioclient(i), i)
    enddo

    !  End Compute-Side Initialization
    !--------------------------------------------------------
    if (STD_OUT_PE) write(*,'(/ A)') "------------------------------------------------------"
    if (STD_OUT_PE) write(*,'(A)')   "Initialization complete, beginning physics integration"
    if (STD_OUT_PE) write(*,'(A)')   "------------------------------------------------------"

    call component_loop(components(1:n_nests), options(1:n_nests), boundary(1:n_nests), ioclient(1:n_nests))

    call component_program_end(components(1:n_nests), options(1:n_nests))

    CALL MPI_Finalize()

contains

    subroutine init_component(component, options, boundary, ioclient, nest_index)
        implicit none
        class(flow_obj_t), intent(inout) :: component
        type(options_t), intent(inout) :: options(:)
        type(boundary_t), intent(inout) :: boundary
        type(ioclient_t), intent(inout) :: ioclient
        integer, intent(in) :: nest_index


        ! check if the type of component is a domain or an ioserver
        select type (component)
            type is (domain_t)
                call component%total_timer%start()
                call component%initialization_timer%start()    

                if (STD_OUT_PE) write(*,"(/ A22,I2,A2,A,A16)") "-------------- Domain ",nest_index," (",trim(options(nest_index)%domain%init_conditions_file),") --------------"
                call init_model(options, component, boundary, ioclient, nest_index)

                call component%total_timer%stop()
                call component%initialization_timer%stop()
            type is (ioserver_t)
                call component%init(options,nest_index)
            type is (flow_obj_t)
                write(*,*) "normal, flow object behavior"
        end select

    end subroutine init_component

    subroutine wake_component(component, options, boundary, ioclient)
        implicit none
        class(flow_obj_t), intent(inout) :: component(:)
        type(options_t), intent(inout) :: options
        type(boundary_t), intent(inout) :: boundary
        type(ioclient_t), intent(inout) :: ioclient

        ! check if the type of component is a domain or an ioserver
        select type (component)
            type is (domain_t)
                call component(options%nest_indx)%total_timer%start()
                call component(options%nest_indx)%initialization_timer%start()    

                if (STD_OUT_PE) write(*,"(/ A22,I2,A2,A,A16)") "-------------- Domain ",options%nest_indx," (",trim(options%domain%init_conditions_file),") --------------"
                call wake_nest(options, component(options%nest_indx), boundary, ioclient)

                call component(options%nest_indx)%total_timer%stop()
                call component(options%nest_indx)%initialization_timer%stop() 
            type is (ioserver_t)
                !Get initial conditions
                if (options%general%parent_nest == 0) call component(options%nest_indx)%read_file()

                if (options%restart%restart) then
                    call component(options%nest_indx)%read_restart_file(options)
                else
                    call component(options%nest_indx)%write_file()
                endif

                !See if any nests needs updating
                if (should_update_nests(component,options)) then
                    call update_component_nest(component,options,ioclient)
                endif

                ! Batch off another round of file reads so that we are always one step ahead of the compute team
                ! This is needed to ensure that the compute team does not wait unnescecarily on us
                if (options%general%parent_nest == 0) call component(options%nest_indx)%read_file()

            type is (flow_obj_t)
                write(*,*) "normal, flow object behavior"
        end select

    end subroutine wake_component

    subroutine update_component_nest(component,options,ioclient)
        implicit none
        class(flow_obj_t), intent(inout) :: component(:)
        type(options_t), intent(in) :: options
        type(ioclient_t), intent(inout) :: ioclient

        integer :: n
        ! check if the type of component is a domain or an ioserver
        select type (component)
            type is (domain_t)
                call ioclient%update_nest(component(options%nest_indx))
            type is (ioserver_t)
                ! This call will gather the model state of the forcing fields from the nest parent
                call component(options%nest_indx)%gather_forcing()

                ! now loop over all child nests
                do n = 1, size(options%general%child_nests)
                    !Test if we can update the child nest
                    if ( can_update_child_nest(component(options%nest_indx),component(options%general%child_nests(n))) ) then
                        ! This call will distribute the model state of the forcing fields to the child nest
                        call component(options%nest_indx)%distribute_forcing(component(options%general%child_nests(n)), n)
                    endif
                enddo
        type is (flow_obj_t)
                write(*,*) "normal, flow object behavior"
        end select

    end subroutine update_component_nest

    subroutine component_write(component, options, ioclient)
        implicit none
        class(flow_obj_t), intent(inout) :: component
        type(options_t), intent(inout) :: options
        type(ioclient_t), intent(inout) :: ioclient

        ! check if the type of component is a domain or an ioserver
        select type (component)
            type is (domain_t)
                if (STD_OUT_PE) write(*,*) "Writing output file"
                call component%output_timer%start()
                call ioclient%push(component)
                call component%output_timer%stop()
            type is (ioserver_t)
                call component%write_file()
            type is (flow_obj_t)
                write(*,*) "normal, flow object behavior"
        end select

    end subroutine component_write

    subroutine component_read(component, options, boundary, ioclient)
        implicit none
        class(flow_obj_t), intent(inout) :: component
        type(options_t), intent(inout) :: options
        type(boundary_t), intent(inout) :: boundary
        type(ioclient_t), intent(inout) :: ioclient

        ! check if the type of component is a domain or an ioserver
        select type (component)
            type is (domain_t)
                call component%total_timer%start()

                !Switch nest contexts if needed
                call switch_nest_context(component,options)

                ! -----------------------------------------------------
                !
                !  Read input data
                !
                ! -----------------------------------------------------
                if (STD_OUT_PE) write(*,"(/ A22,I2,A2,A,A16)") "-------------- Domain ",component%nest_indx," (",trim(options%domain%init_conditions_file),") --------------"
                if (STD_OUT_PE) write(*,*) "Updating Boundary conditions"

                call component%input_timer%start()

                call ioclient%receive(boundary, component)
                
                ! after reading all variables that can be read, not compute any remaining variables (e.g. z from p+ps)
                call boundary%update_computed_vars(options, update=.True.)
                call component%interpolate_forcing(boundary, update=.True.)

                ! Make the boundary condition dXdt values into units of [X]/s
                call boundary%update_delta_fields(component%next_input - component%sim_time)
                call component%update_delta_fields(component%next_input - component%sim_time)

                call component%input_timer%stop()
            type is (ioserver_t)
                !See if we even have files to read
                if (component%files_to_read) then
                    call component%read_file()
                endif
            type is (flow_obj_t)
                write(*,*) "normal, flow object behavior"
        end select

    end subroutine component_read

    subroutine component_main_loop(component, options, boundary)
        implicit none
        class(flow_obj_t), intent(inout) :: component
        type(options_t), intent(in) :: options
        type(boundary_t), intent(inout) :: boundary

        ! check if the type of component is a domain or an ioserver
        select type (component)
            type is (domain_t)
                ! -----------------------------------------------------
                !
                !  Integrate physics forward in time
                !
                ! -----------------------------------------------------
                if (STD_OUT_PE) write(*,*) "Running Physics"
                if (STD_OUT_PE) write(*,*) "  Model time = ", trim(component%sim_time%as_string())
                if (STD_OUT_PE) write(*,*) "   End  time = ", trim(component%end_time%as_string())
                if (STD_OUT_PE) write(*,*) "  Next Input = ", trim(component%next_input%as_string())
                if (STD_OUT_PE) write(*,*) "  Next Output= ", trim(component%next_output%as_string())
                if (STD_OUT_PE) flush(output_unit)
                
                ! this is the meat of the model physics, run all the physics for the current time step looping over internal timesteps
                call component%physics_timer%start()
                call step(component, boundary, component%next_flow_event(), options)
                call component%physics_timer%stop()
            type is (ioserver_t)
                call component%set_sim_time(component%next_flow_event())

            type is (flow_obj_t)
                write(*,*) "normal, flow object behavior"
        end select

    end subroutine component_main_loop

    subroutine component_loop(components, options, boundary, ioclient)
        implicit none
        class(flow_obj_t), intent(inout) :: components(:)
        type(options_t), intent(inout) :: options(:)
        type(boundary_t), intent(inout):: boundary(:)
        type(ioclient_t), intent(inout):: ioclient(:)

        integer :: i, n_nests

        n_nests = size(components)

        do while (all_nests_not_done(components))
            do i = 1, n_nests
                if (nest_next_up(components,options(i))) then
                    call wake_component(components, options(i), boundary(i), ioclient(i))
                    cycle
                endif

                if (components(i)%dead_or_asleep()) cycle

                call component_read(components(i), options(i), boundary(i), ioclient(i))

                do while ( .not.(components(i)%time_for_input()) .and. .not.(components(i)%ended) )
                    call component_main_loop(components(i), options(i), boundary(i))

                    ! If it is time for an output, do. But, if we are about to exit this loop, then 
                    ! skip ahead to updating the nest, since this will be a bottleneck for execution.
                    if (components(i)%time_for_output()) then
                        if (.not.(components(i)%time_for_input() )) then
                            call component_write(components(i), options(i), ioclient(i))
                        endif
                    endif
                end do

                if (should_update_nests(components,options(i))) then
                    call update_component_nest(components,options(i),ioclient(i))
                endif

                !Do a write here which may have been skipped above
                if (components(i)%time_for_output()) then
                    call component_write(components(i), options(i), ioclient(i))
                endif

                call component_end_of_nest_loop(components(i), options(i))

            enddo
        enddo

    end subroutine component_loop


    subroutine component_end_of_nest_loop(component,options)
        implicit none
        class(flow_obj_t), intent(inout) :: component
        type(options_t), intent(in) :: options

        ! check if the type of component is a domain or an ioserver
        select type (component)
            type is (domain_t)
                call component%total_timer%stop()

                if (component%ended) then
                    call end_nest_context()
                    call component%release()
                    if (STD_OUT_PE) write(*,*) "Domain ",i," has reached the end of its run time."
                endif
            type is (ioserver_t)
            type is (flow_obj_t)
                write(*,*) "normal, flow object behavior"
        end select

    end subroutine component_end_of_nest_loop

    subroutine component_program_end(component, options)
        implicit none
        class(flow_obj_t), intent(inout) :: component(:)
        type(options_t), intent(in) :: options(:)

        integer :: i, n_nests
        real    :: t_val, t_val2, t_val3

        n_nests = size(component)

        select type (component)
            type is (domain_t)
                t_val = 0
                do i = 1,n_nests
                    t_val = t_val + component(i)%total_timer%mean(component(1)%compute_comms)
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
                    t_val = component(i)%total_timer%mean(component(1)%compute_comms)
                    t_val2 = component(i)%total_timer%min(component(1)%compute_comms)
                    t_val3 = component(i)%total_timer%max(component(1)%compute_comms)
                    if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "total", ":", t_val, " | ", t_val2, " | ", t_val3
                    t_val = component(i)%initialization_timer%mean(component(1)%compute_comms)
                    t_val2 = component(i)%initialization_timer%min(component(1)%compute_comms)
                    t_val3 = component(i)%initialization_timer%max(component(1)%compute_comms)
                    if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "init", ":", t_val, " | ", t_val2, " | ", t_val3
                    t_val = component(i)%input_timer%mean(component(1)%compute_comms)
                    t_val2 = component(i)%input_timer%min(component(1)%compute_comms)
                    t_val3 = component(i)%input_timer%max(component(1)%compute_comms)
                    if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "input", ":", t_val, " | ", t_val2, " | ", t_val3
                    t_val = component(i)%output_timer%mean(component(1)%compute_comms)
                    t_val2 = component(i)%output_timer%min(component(1)%compute_comms)
                    t_val3 = component(i)%output_timer%max(component(1)%compute_comms)
                    if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "output", ":", t_val, " | ", t_val2, " | ", t_val3
                    t_val = component(i)%physics_timer%mean(component(1)%compute_comms)
                    t_val2 = component(i)%physics_timer%min(component(1)%compute_comms)
                    t_val3 = component(i)%physics_timer%max(component(1)%compute_comms)
                    if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "physics", ":", t_val, " | ", t_val2, " | ", t_val3
                    t_val = component(i)%mp_timer%mean(component(1)%compute_comms)
                    t_val2 = component(i)%mp_timer%min(component(1)%compute_comms)
                    t_val3 = component(i)%mp_timer%max(component(1)%compute_comms)
                    if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "microphysics", ":", t_val, " | ", t_val2, " | ", t_val3
                    t_val = component(i)%adv_timer%mean(component(1)%compute_comms)
                    t_val2 = component(i)%adv_timer%min(component(1)%compute_comms)
                    t_val3 = component(i)%adv_timer%max(component(1)%compute_comms)
                    if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "advection", ":", t_val, " | ", t_val2, " | ", t_val3
                    t_val = component(i)%sum_timer%mean(component(1)%compute_comms)
                    t_val2 = component(i)%sum_timer%min(component(1)%compute_comms)
                    t_val3 = component(i)%sum_timer%max(component(1)%compute_comms)
                    if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "advection_sum", ":", t_val, " | ", t_val2, " | ", t_val3
                    t_val = component(i)%flux_timer%mean(component(1)%compute_comms)
                    t_val2 = component(i)%flux_timer%min(component(1)%compute_comms)
                    t_val3 = component(i)%flux_timer%max(component(1)%compute_comms)
                    if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "advection_flux", ":", t_val, " | ", t_val2, " | ", t_val3
                    t_val = component(i)%flux_up_timer%mean(component(1)%compute_comms)
                    t_val2 = component(i)%flux_up_timer%min(component(1)%compute_comms)
                    t_val3 = component(i)%flux_up_timer%max(component(1)%compute_comms)
                    if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "advection_flux_up", ":", t_val, " | ", t_val2, " | ", t_val3
                    t_val = component(i)%flux_corr_timer%mean(component(1)%compute_comms)
                    t_val2 = component(i)%flux_corr_timer%min(component(1)%compute_comms)
                    t_val3 = component(i)%flux_corr_timer%max(component(1)%compute_comms)
                    if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "advection_flux_corr", ":", t_val, " | ", t_val2, " | ", t_val3
                    t_val = component(i)%adv_wind_timer%mean(component(1)%compute_comms)
                    t_val2 = component(i)%adv_wind_timer%min(component(1)%compute_comms)
                    t_val3 = component(i)%adv_wind_timer%max(component(1)%compute_comms)
                    if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "advection_wind", ":", t_val, " | ", t_val2, " | ", t_val3
                    t_val = component(i)%rad_timer%mean(component(1)%compute_comms)
                    t_val2 = component(i)%rad_timer%min(component(1)%compute_comms)
                    t_val3 = component(i)%rad_timer%max(component(1)%compute_comms)
                    if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "radiation", ":", t_val, " | ", t_val2, " | ", t_val3
                    t_val = component(i)%lsm_timer%mean(component(1)%compute_comms)
                    t_val2 = component(i)%lsm_timer%min(component(1)%compute_comms)
                    t_val3 = component(i)%lsm_timer%max(component(1)%compute_comms)
                    if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "LSM", ":", t_val, " | ", t_val2, " | ", t_val3
                    t_val = component(i)%pbl_timer%mean(component(1)%compute_comms)
                    t_val2 = component(i)%pbl_timer%min(component(1)%compute_comms)
                    t_val3 = component(i)%pbl_timer%max(component(1)%compute_comms)
                    if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "PBL", ":", t_val, " | ", t_val2, " | ", t_val3
                    t_val = component(i)%forcing_timer%mean(component(1)%compute_comms)
                    t_val2 = component(i)%forcing_timer%min(component(1)%compute_comms)
                    t_val3 = component(i)%forcing_timer%max(component(1)%compute_comms)
                    if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "forcing", ":", t_val, " | ", t_val2, " | ", t_val3
                    t_val = component(i)%wind_bal_timer%mean(component(1)%compute_comms)
                    t_val2 = component(i)%wind_bal_timer%min(component(1)%compute_comms)
                    t_val3 = component(i)%wind_bal_timer%max(component(1)%compute_comms)
                    if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "wind bal", ":", t_val, " | ", t_val2, " | ", t_val3
                    t_val = component(i)%diagnostic_timer%mean(component(1)%compute_comms)
                    t_val2 = component(i)%diagnostic_timer%min(component(1)%compute_comms)
                    t_val3 = component(i)%diagnostic_timer%max(component(1)%compute_comms)
                    if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "diagnostic", ":", t_val, " | ", t_val2, " | ", t_val3
                    t_val = component(i)%send_timer%mean(component(1)%compute_comms)
                    t_val2 = component(i)%send_timer%min(component(1)%compute_comms)
                    t_val3 = component(i)%send_timer%max(component(1)%compute_comms)
                    if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "halo-exchange(send)", ":", t_val, " | ", t_val2, " | ", t_val3
                    t_val = component(i)%ret_timer%mean(component(1)%compute_comms)
                    t_val2 = component(i)%ret_timer%min(component(1)%compute_comms)
                    t_val3 = component(i)%ret_timer%max(component(1)%compute_comms)
                    if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "halo-exchange(retrieve)", ":", t_val, " | ", t_val2, " | ", t_val3
                    t_val = component(i)%wait_timer%mean(component(1)%compute_comms)
                    t_val2 = component(i)%wait_timer%min(component(1)%compute_comms)
                    t_val3 = component(i)%wait_timer%max(component(1)%compute_comms)
                    if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "halo-exchange(wait)", ":", t_val, " | ", t_val2, " | ", t_val3
                    t_val = component(i)%wind_timer%mean(component(1)%compute_comms)
                    t_val2 = component(i)%wind_timer%min(component(1)%compute_comms)
                    t_val3 = component(i)%wind_timer%max(component(1)%compute_comms)
                    if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "winds", ":", t_val, " | ", t_val2, " | ", t_val3
                enddo
    
            type is (ioserver_t)
                call component(1)%close_files()
        end select

    end subroutine component_program_end
    
    subroutine read_co(nml_file, info, gen_nml, only_namelist_check)
        implicit none
        character(len=kMAX_FILE_LENGTH), intent(out) :: nml_file
        logical, intent(out) :: info, gen_nml, only_namelist_check

        integer :: cnt, p
        character(len=kMAX_FILE_LENGTH) :: first_arg, arg, default
        character(len=kMAX_FILE_LENGTH), allocatable :: keywords(:)
        logical :: file_exists

        nml_file = ""
        info = .False.
        gen_nml = .False.
        only_namelist_check = .False.
        first_arg = ""

        cnt = command_argument_count()

        if (cnt > 0) then
            ! get first command line argument
            call get_command_argument(1,first_arg)
        endif

        ! If there are no command line arguments, throw error
        if ( (cnt == 0 .or. first_arg=='-h' .or. first_arg=='--help') .and. STD_OUT_PE) then
            write(*,*) "Usage: ./HICAR [-v [variable_name ...|--all]] [--check-nml] [--gen-nml] namelist_file"
            write(*,*) "    -v [variable_name ...|--all]: Print information about the namelist variable(s) variable_name, ... "
            write(*,*) "                                  --all prints out information for all namelist variables."
            write(*,*) "    --check-nml:                  Check the namelist file for errors without running the model."
            write(*,*) "    --gen-nml:                    Generate a namelist file with default values."
            write(*,*) "    --out-vars [keywords]:        List all output variables which are related to the space-separated list of keywords."
            write(*,*)
            write(*,*) "    namelist_file:                The name of the namelist file to use."
            write(*,*)
            write(*,*) "    Example to generate a namelist with default values:                ./HICAR --gen-nml namelist_file.nml"
            write(*,*) "    Example to check namelist:                                         ./HICAR --check-nml namelist_file.nml"
            write(*,*) "    Example to run model:                                              ./HICAR namelist_file.nml"
            write(*,*) "    Example to list all output variables related to wind or snow:      ./HICAR --out-vars wind snow"
            write(*,*) "    Example to learn about a namelist variable:                        ./HICAR -v mp"
            write(*,*) "    Example to generate namelist variable documentation:               ./HICAR -v --all > namelist_doc.txt"
            write(*,*)


            stop
        endif

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
        elseif (first_arg == '--out-vars') then
            if (cnt >= 2) then
                allocate(keywords(cnt-1))
                do p = 2, cnt
                    call get_command_argument(p, arg)
                    keywords(p-1) = arg
                end do
                call list_output_vars(keywords)
                stop
            elseif (cnt < 2) then
                if (STD_OUT_PE) write(*,*) "ERROR: No keyword provided with the --out-vars flag."
                stop
            endif
        else
            nml_file = first_arg
        endif

        if (.not.first_arg=='-v' .and. .not.first_arg=='--gen-nml' .and. .not.first_arg=='--out-vars') then
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
