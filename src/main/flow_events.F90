module flow_events
    use options_interface,  only : options_t
    use domain_interface,   only : domain_t
    use boundary_interface, only : boundary_t
    use flow_object_interface, only : flow_obj_t
    use icar_constants!,             only : kITERATIVE_WINDS, kWIND_LINEAR
    use ioserver_interface,         only : ioserver_t
    use ioclient_interface,         only : ioclient_t
    use nest_manager, only: all_nests_not_done, nest_next_up, should_update_nests, &
        can_update_child_nest, end_nest_context, switch_nest_context, wake_nest
    use time_step, only: step
    use initialization, only: init_model
    use time_object, only: Time_type
    use iso_fortran_env
    use mpi_f08

    implicit none
    private
    public:: component_init, component_loop, component_program_end

contains

subroutine component_init(component, options, boundary, ioclient, nest_index)
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
        class default
    end select

end subroutine component_init

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
            if (options%general%parent_nest == 0) call component_read(component(options%nest_indx),options,boundary,ioclient)

            if (options%restart%restart) then
                call component(options%nest_indx)%read_restart_file(options)
            else
                call component_write(component(options%nest_indx),ioclient)
            endif

            !See if any nests needs updating
            if (should_update_nests(component,options)) then
                call update_component_nest(component,options,ioclient)
            endif

            ! Batch off another round of file reads so that we are always one step ahead of the compute team
            ! This is needed to ensure that the compute team does not wait unnescecarily on us
            if (options%general%parent_nest == 0) call component_read(component(options%nest_indx),options,boundary,ioclient)

        class default
            if (.not.(ioclient%parent_comms==MPI_COMM_NULL)) then
                call MPI_Barrier(MPI_COMM_WORLD)
                if (.not.(options%restart%restart)) then
                    call component(options%nest_indx)%increment_output_time()
                    call MPI_Barrier(MPI_COMM_WORLD)
                endif
                ! "update nest"
                if (size(options%general%child_nests) > 0) then
                    call update_component_nest(component,options,ioclient)
                endif
            else
                if (options%general%parent_nest == 0) then
                    call MPI_Barrier(MPI_COMM_WORLD)
                endif
                if (.not.(options%restart%restart)) then
                    call component(options%nest_indx)%increment_output_time()
                    call MPI_Barrier(MPI_COMM_WORLD)
                endif

                !See if any nests needs updating
                if (should_update_nests(component,options)) then
                    call update_component_nest(component,options,ioclient)
                endif

                if (options%general%parent_nest == 0) then
                    call MPI_Barrier(MPI_COMM_WORLD)
                endif
            endif
    end select

    component(options%nest_indx)%started = .true.

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
            if (STD_OUT_PE_IO) write(*,"(/ A22,I2,A16)") "-------------- IOserver",component(options%nest_indx)%nest_indx," --------------"
            if (STD_OUT_PE_IO) write(*,*) "Gathering forcings"
            if (STD_OUT_PE_IO) write(*,*) "  Model time = ", trim(component(options%nest_indx)%sim_time%as_string())

            call component(options%nest_indx)%gather_forcing()
            if (STD_OUT_PE_IO) write(*,*) "Completed gathering forcings"
    
            ! now loop over all child nests
            do n = 1, size(options%general%child_nests)
                !Test if we can update the child nest
                if ( can_update_child_nest(component(options%nest_indx),component(options%general%child_nests(n))) ) then
                    ! This call will distribute the model state of the forcing fields to the child nest
                    if (STD_OUT_PE_IO) write(*,*) "Distributing forcings to child nest: ",options%general%child_nests(n)
                    call component(options%nest_indx)%distribute_forcing(component(options%general%child_nests(n)), n)
                    if (STD_OUT_PE_IO) write(*,*) "Distributed forcings to child nest: ",options%general%child_nests(n)

                endif
            enddo
        class default
            ! if ioclient is null, then we are acting as an ioserver
            if (.not.(ioclient%parent_comms==MPI_COMM_NULL)) then
                call MPI_Barrier(MPI_COMM_WORLD)
            else
                ! now loop over all child nests
                call MPI_Barrier(MPI_COMM_WORLD)
                do n = 1, size(options%general%child_nests)
                    !Test if we can update the child nest
                    if ( can_update_child_nest(component(options%nest_indx),component(options%general%child_nests(n))) ) then
                        ! This call will distribute the model state of the forcing fields to the child nest
                        call MPI_Barrier(MPI_COMM_WORLD)
                    endif
                enddo
            endif
    end select

end subroutine update_component_nest

subroutine component_write(component, ioclient)
    implicit none
    class(flow_obj_t), intent(inout) :: component
    type(ioclient_t), intent(inout) :: ioclient

    ! check if the type of component is a domain or an ioserver
    select type (component)
        type is (domain_t)
            if (STD_OUT_PE) write(*,*) "Writing output file"
            call component%output_timer%start()
            call ioclient%push(component)
            call component%output_timer%stop()
        type is (ioserver_t)
            if (STD_OUT_PE_IO) write(*,"(/ A22,I2,A16)") "-------------- IOserver",component%nest_indx," --------------"
            if (STD_OUT_PE_IO) write(*,*) "Writing out domain"
            if (STD_OUT_PE_IO) write(*,*) "  Model time = ", trim(component%sim_time%as_string())
            if (STD_OUT_PE_IO) write(*,*) "  Next Output = ", trim(component%next_output%as_string())
            call component%write_file()
            if (STD_OUT_PE_IO) write(*,*) "Wrote out domain"
        class default
            call component%increment_output_time()
            call MPI_Barrier(MPI_COMM_WORLD)

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
            call component%update_delta_fields()

            call component%input_timer%stop()
        type is (ioserver_t)
            !See if we even have files to read
            if (STD_OUT_PE_IO) write(*,"(/ A22,I2,A16)") "-------------- IOserver",component%nest_indx," --------------"
            if (STD_OUT_PE_IO) write(*,*) "Reading in Boundary conditions"
            if (STD_OUT_PE_IO) write(*,*) "  Model time = ", trim(component%sim_time%as_string())
            if (STD_OUT_PE_IO) write(*,*) "  Next Input = ", trim(component%next_input%as_string())
            call component%read_file()
            if (STD_OUT_PE_IO) write(*,*) "Read in Boundary conditions"

        class default
            if (.not.(ioclient%parent_comms==MPI_COMM_NULL)) then
                call MPI_Barrier(MPI_COMM_WORLD)
            else
                ! if ioclient comms is null, then we are acting as an ioserver
                if (options%general%parent_nest == 0) then
                    if (component%sim_time < (component%end_time-component%input_dt-component%input_dt)) then
                        call MPI_Barrier(MPI_COMM_WORLD)
                    endif
                endif
            endif
    end select

end subroutine component_read

subroutine component_main_loop(component, options)
    implicit none
    class(flow_obj_t), intent(inout) :: component
    type(options_t), intent(in) :: options

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
            call step(component, component%next_flow_event(), options)
            call component%physics_timer%stop()
        type is (ioserver_t)
            call component%set_sim_time(component%next_flow_event())
        class default
            if (STD_OUT_PE) write(*,*) "  Model time = ", trim(component%sim_time%as_string())
            if (STD_OUT_PE) write(*,*) "   End  time = ", trim(component%end_time%as_string())
            if (STD_OUT_PE) write(*,*) "  Next Input = ", trim(component%next_input%as_string())
            if (STD_OUT_PE) write(*,*) "  Next Output= ", trim(component%next_output%as_string())
            if (STD_OUT_PE) flush(output_unit)

            call component%set_sim_time(component%next_flow_event())
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

            call components(i)%increment_input_time()
            call component_read(components(i), options(i), boundary(i), ioclient(i))

            do while ( .not.(components(i)%time_for_input()) .and. .not.(components(i)%ended) )
                call component_main_loop(components(i), options(i))

                if (should_update_nests(components,options(i))) then
                    call update_component_nest(components,options(i),ioclient(i))
                endif
    
                if (components(i)%time_for_output()) then
                    call component_write(components(i), ioclient(i))
                endif
            end do

            call component_end_of_nest_loop(components,i)
        enddo
    enddo

end subroutine component_loop


subroutine component_end_of_nest_loop(component,nest_indx)
    implicit none
    class(flow_obj_t), intent(inout) :: component(:)
    integer, intent(in) :: nest_indx

    ! check if the type of component is a domain or an ioserver
    select type (component)
        type is (domain_t)
            call component(nest_indx)%total_timer%stop()

            if (component(nest_indx)%ended) then
                call end_nest_context()
                call component(nest_indx)%release()
                if (STD_OUT_PE) write(*,*) "Domain ",nest_indx," has reached the end of its run time."
            endif
        class default
            ! Handle generic flow_obj_t case
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
        class default

    end select

end subroutine component_program_end

end module flow_events
