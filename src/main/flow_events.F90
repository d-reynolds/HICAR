module flow_events
    use iso_fortran_env
    use mpi_f08
    use options_interface,  only : options_t
    use domain_interface,   only : domain_t
    use boundary_interface, only : boundary_t
    use time_object,        only : Time_type
    use flow_object_interface, only : flow_obj_t, comp_arr_t
    use icar_constants!,             only : kITERATIVE_WINDS, kWIND_LINEAR
    use ioserver_interface,         only : ioserver_t
    use ioclient_interface,         only : ioclient_t
    use nest_manager, only: any_nests_not_done, nest_next_up, should_update_nests, &
        can_update_child_nest, end_nest_context, switch_nest_context, wake_nest
    use time_step, only: step
    use initialization, only: init_model, init_model_state
    use wind_iterative_old, only: finalize_petsc
    use string, only : as_string

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
            if (STD_OUT_PE) flush(output_unit)
            call init_model(options, component, boundary, ioclient, nest_index)

            call component%total_timer%stop()
            call component%initialization_timer%stop()
        type is (ioserver_t)
            call component%init(options,nest_index)
        class default
    end select

end subroutine component_init

subroutine wake_component(comp_arr, options, boundary, ioclient)
    implicit none
    type(comp_arr_t), intent(inout) :: comp_arr(:)
    type(options_t), intent(inout) :: options
    type(boundary_t), intent(inout) :: boundary
    type(ioclient_t), intent(inout) :: ioclient

    ! check if the type of component is a domain or an ioserver
    associate (comp => comp_arr(options%nest_indx)%comp)
    select type (comp)
        type is (domain_t)
            call comp%total_timer%start()
            call comp%initialization_timer%start()    

            if (STD_OUT_PE) write(*,"(/ A22,I2,A2,A,A16)") "-------------- Domain ",options%nest_indx," (",trim(options%domain%init_conditions_file),") --------------"
            if (STD_OUT_PE) flush(output_unit)
            
            call wake_nest(comp)
            call init_model_state(options, comp, boundary, ioclient)

            call comp%total_timer%stop()
            call comp%initialization_timer%stop() 
        type is (ioserver_t)

            !Get initial conditions
            if (options%general%parent_nest == 0) call component_read(comp,options,boundary,ioclient)

            if (options%restart%restart) then
                call comp%read_restart_file(options)
            else
                call component_write(comp,ioclient)
            endif

            !See if any nests needs updating
            if (should_update_nests(comp_arr,options)) then
                call update_component_nest(comp_arr,options,ioclient)
            endif

            ! Batch off another round of file reads so that we are always one step ahead of the compute team
            ! This is needed to ensure that the compute team does not wait unnescecarily on us
            if (options%general%parent_nest == 0) call component_read(comp,options,boundary,ioclient)

        class default
            if (.not.(ioclient%parent_comms==MPI_COMM_NULL)) then
                call component_read(comp,options,boundary,ioclient)
                if (.not.(options%restart%restart)) then
                    call component_write(comp,ioclient)
                endif
                ! "update nest"
                if (size(options%general%child_nests) > 0) then
                    call update_component_nest(comp_arr,options,ioclient)
                endif
            else
                if (options%general%parent_nest == 0) then
                    call component_read(comp,options,boundary,ioclient)
                endif
                if (.not.(options%restart%restart)) then
                    call component_write(comp,ioclient)
                endif

                !See if any nests needs updating
                if (should_update_nests(comp_arr,options)) then
                    call update_component_nest(comp_arr,options,ioclient)
                endif

                if (options%general%parent_nest == 0) then
                    call component_read(comp,options,boundary,ioclient)
                endif
            endif
    end select

    comp%started = .true.
    end associate
end subroutine wake_component

subroutine update_component_nest(comp_arr,options,ioclient)
    implicit none
    type(comp_arr_t), intent(inout) :: comp_arr(:)
    type(options_t), intent(in) :: options
    type(ioclient_t), intent(inout) :: ioclient

    class(flow_obj_t), allocatable :: child_nest
    integer :: n
    type(Time_type) :: sim_time_safety_under
    ! check if the type of component is a domain or an ioserver

    associate (comp => comp_arr(options%nest_indx)%comp)
    select type (comp)
        type is (domain_t)
            call ioclient%update_nest(comp)
        type is (ioserver_t)
            ! This call will gather the model state of the forcing fields from the nest parent
            if (STD_OUT_PE_IO) write(*,"(/ A23,I2,A16)") "-------------- IOserver",comp%nest_indx," --------------"
            if (STD_OUT_PE_IO) write(*,*) "Gathering forcings"
            if (STD_OUT_PE_IO) write(*,*) "  Model time = ", trim(as_string(comp%sim_time))
            if (STD_OUT_PE_IO) write(*,"(A23,I2,A16 /)") "-------------- IOserver",comp%nest_indx," --------------"
            if (STD_OUT_PE_IO) flush(output_unit)

            call comp%gather_forcing()
            if (STD_OUT_PE_IO) write(*,"(/ A23,I2,A16)") "-------------- IOserver",comp%nest_indx," --------------"
            if (STD_OUT_PE_IO) write(*,*) "Completed gathering forcings"
            if (STD_OUT_PE_IO) write(*,"(A23,I2,A16 /)") "-------------- IOserver",comp%nest_indx," --------------"
            if (STD_OUT_PE_IO) flush(output_unit)

            ! now loop over all child nests
            do n = 1, size(options%general%child_nests)
                !Test if we can update the child nest
                associate (child_nest => comp_arr(options%general%child_nests(n))%comp)
                select type(child_nest)
                type is (ioserver_t)
                    call sim_time_safety_under%set(child_nest%sim_time%mjd() - comp%small_time_delta%days())

                    if ( (comp%sim_time >= sim_time_safety_under .and. .not.(child_nest%ended)) )then
                        ! This call will distribute the model state of the forcing fields to the child nest
                        ! if (STD_OUT_PE_IO) write(*,*) "Distributing forcings to child nest: ",options%general%child_nests(n)
                        call comp%distribute_forcing(child_nest, n)
                        if (STD_OUT_PE_IO) write(*,"(/ A23,I2,A16)") "-------------- IOserver",comp%nest_indx," --------------"
                        if (STD_OUT_PE_IO) write(*,*) "Distributed forcings to child nest: ",options%general%child_nests(n)
                        if (STD_OUT_PE_IO) write(*,"(A23,I2,A16 /)") "-------------- IOserver",comp%nest_indx," --------------"
                        if (STD_OUT_PE_IO) flush(output_unit)

                    endif
                end select
                end associate
            enddo
        class default
            ! if ioclient is null, then we are acting as an ioserver
            if (.not.(ioclient%parent_comms==MPI_COMM_NULL)) then
                    write(*,*) "Client Barr 3 -- gather on nest ", options%nest_indx
                call MPI_Barrier(MPI_COMM_WORLD)
            else
                ! now loop over all child nests
                write(*,*) "Server Barr 3 -- gather on nest ", options%nest_indx
                call MPI_Barrier(MPI_COMM_WORLD)

                do n = 1, size(options%general%child_nests)
                    !Test if we can update the child nest
                    child_nest = comp_arr(options%general%child_nests(n))%comp
                    call sim_time_safety_under%set(child_nest%sim_time%mjd() - comp%small_time_delta%days())

                    if ( (comp%sim_time >= sim_time_safety_under .and. .not.(child_nest%ended)) )then
                            ! This call will distribute the model state of the forcing fields to the child nest
                        write(*,*) "Server Barr 4 -- scatter to nest ", child_nest%nest_indx
                        call MPI_Barrier(MPI_COMM_WORLD)
                    endif
                enddo
            endif
    end select
    end associate

end subroutine update_component_nest

subroutine component_write(component, ioclient)
    implicit none
    class(flow_obj_t), intent(inout) :: component
    type(ioclient_t), intent(inout) :: ioclient

    ! check if the type of component is a domain or an ioserver
    select type (component)
        type is (domain_t)
            if (STD_OUT_PE) write(*,*) "Writing output file"
            if (STD_OUT_PE) flush(output_unit)
            call component%output_timer%start()
            call ioclient%push(component)
            call component%output_timer%stop()
        type is (ioserver_t)
            if (STD_OUT_PE_IO) write(*,"(/ A23,I2,A16)") "-------------- IOserver",component%nest_indx," --------------"
            if (STD_OUT_PE_IO) write(*,*) "Writing out domain"
            if (STD_OUT_PE_IO) write(*,*) "  Model time = ", trim(as_string(component%sim_time))
            if (STD_OUT_PE_IO) write(*,*) "  Next Output = ", trim(as_string(component%next_output))
            if (STD_OUT_PE_IO) write(*,"(A23,I2,A16 /)") "-------------- IOserver",component%nest_indx," --------------"
            if (STD_OUT_PE_IO) flush(output_unit)
            call component%write_file()
            if (STD_OUT_PE_IO) write(*,"(/ A23,I2,A16)") "-------------- IOserver",component%nest_indx," --------------"
            if (STD_OUT_PE_IO) write(*,*) "Wrote out domain"
            if (STD_OUT_PE_IO) write(*,"(A23,I2,A16 /)") "-------------- IOserver",component%nest_indx," --------------"
            if (STD_OUT_PE_IO) flush(output_unit)
        class default
            call component%increment_output_time()
            write(*,*) "Client/Server Barr 5 -- output on nest ", component%nest_indx
            call MPI_Barrier(MPI_COMM_WORLD)

    end select

end subroutine component_write

subroutine component_read(component, options, boundary, ioclient)
    implicit none
    class(flow_obj_t), intent(inout) :: component
    type(options_t), intent(inout) :: options
    type(boundary_t), intent(inout) :: boundary
    type(ioclient_t), intent(inout) :: ioclient
    type(Time_type) :: end_time_safety_under
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
            if (STD_OUT_PE) flush(output_unit)

            call component%input_timer%start()

            call ioclient%receive(boundary, component)
            
            ! after reading all variables that can be read, not compute any remaining variables (e.g. z from p+ps)
            call boundary%update_computed_vars(options)
            call component%interpolate_forcing(boundary, update=.True.)

            ! Make the boundary condition dXdt values into units of [X]/s
            call component%update_delta_fields()

            call component%input_timer%stop()
        type is (ioserver_t)
            !See if we even have files to read
            if (STD_OUT_PE_IO) write(*,"(/ A23,I2,A16)") "-------------- IOserver",component%nest_indx," --------------"
            if (STD_OUT_PE_IO) write(*,*) "Reading in Boundary conditions"
            if (STD_OUT_PE_IO) write(*,*) "  Model time = ", trim(as_string(component%sim_time))
            if (STD_OUT_PE_IO) write(*,*) "  Next Input = ", trim(as_string(component%next_input))
            if (STD_OUT_PE_IO) write(*,"(A23,I2,A16 /)") "-------------- IOserver",component%nest_indx," --------------"
            if (STD_OUT_PE_IO) flush(output_unit)
            call component%read_file()
            if (STD_OUT_PE_IO) write(*,"(/ A23,I2,A16)") "-------------- IOserver",component%nest_indx," --------------"
            if (STD_OUT_PE_IO) write(*,*) "Read in Boundary conditions"
            if (STD_OUT_PE_IO) write(*,"(A23,I2,A16 /)") "-------------- IOserver",component%nest_indx," --------------"
            if (STD_OUT_PE_IO) flush(output_unit)

        class default
            if (.not.(ioclient%parent_comms==MPI_COMM_NULL)) then
                write(*,*) "Client Barr 5 -- read on nest ", options%nest_indx
                call MPI_Barrier(MPI_COMM_WORLD)
            else
                ! if ioclient comms is null, then we are acting as an ioserver
                if (options%general%parent_nest == 0) then
                    
                    call end_time_safety_under%set(component%end_time%mjd() - component%input_dt%days() - component%small_time_delta%days())

                    if ( component%dead_or_asleep() .or. (component%sim_time < end_time_safety_under) ) then
                        write(*,*) "Server Barr 5 -- read on nest ", options%nest_indx
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
            if (STD_OUT_PE) write(*,*) "  Model time = ", trim(as_string(component%sim_time))
            if (STD_OUT_PE) write(*,*) "   End  time = ", trim(as_string(component%end_time))
            if (STD_OUT_PE) write(*,*) "  Next Input = ", trim(as_string(component%next_input))
            if (STD_OUT_PE) write(*,*) "  Next Output= ", trim(as_string(component%next_output))
            if (STD_OUT_PE) flush(output_unit)
            
            ! this is the meat of the model physics, run all the physics for the current time step looping over internal timesteps
            call component%physics_timer%start()
            call step(component, component%next_flow_event(), options)
            call component%physics_timer%stop()
        type is (ioserver_t)
            call component%set_sim_time(component%next_flow_event())
        class default
            if (STD_OUT_PE) write(*,*) "  Model time = ", trim(as_string(component%sim_time))
            if (STD_OUT_PE) write(*,*) "   End  time = ", trim(as_string(component%end_time))
            if (STD_OUT_PE) write(*,*) "  Next Input = ", trim(as_string(component%next_input))
            if (STD_OUT_PE) write(*,*) "  Next Output= ", trim(as_string(component%next_output))
            if (STD_OUT_PE) flush(output_unit)

            call component%set_sim_time(component%next_flow_event())
    end select

end subroutine component_main_loop

subroutine component_loop(components, options, boundary, ioclient)
    implicit none
    type(comp_arr_t), intent(inout) :: components(:)
    type(options_t), intent(inout) :: options(:)
    type(boundary_t), intent(inout):: boundary(:)
    type(ioclient_t), intent(inout):: ioclient(:)

    integer :: i, n_nests

    n_nests = options(1)%general%nests

    do while (any_nests_not_done(components))
        do i = 1, n_nests
            if (nest_next_up(components,options(i))) then
                call wake_component(components, options(i), boundary(i), ioclient(i))
                cycle
            endif

            if (components(i)%comp%dead_or_asleep()) cycle

            call components(i)%comp%increment_input_time()
            call component_read(components(i)%comp, options(i), boundary(i), ioclient(i))

            do while ( .not.(components(i)%comp%time_for_input()) .and. .not.(components(i)%comp%ended) )
                call component_main_loop(components(i)%comp, options(i))

                if (should_update_nests(components,options(i))) then
                    call update_component_nest(components,options(i),ioclient(i))
                endif

                if (components(i)%comp%time_for_output()) then
                    call component_write(components(i)%comp, ioclient(i))
                endif
            end do

            call component_end_of_nest_loop(components(i)%comp,i)
        enddo
    enddo

end subroutine component_loop


subroutine component_end_of_nest_loop(component,nest_indx)
    implicit none
    class(flow_obj_t), intent(inout) :: component
    integer, intent(in) :: nest_indx

    ! check if the type of component is a domain or an ioserver
    select type (component)
        type is (domain_t)
            call component%total_timer%stop()

            if (component%ended) then
                call end_nest_context()
                call component%release()
                if (STD_OUT_PE) write(*,*) "Domain ",nest_indx," has reached the end of its run time."
                if (STD_OUT_PE) flush(output_unit)
            endif
        class default
            ! Handle generic flow_obj_t case
    end select

end subroutine component_end_of_nest_loop

subroutine component_program_end(component, options)
    implicit none
    type(comp_arr_t), intent(inout) :: component(:)
    type(options_t), intent(in) :: options(:)

    integer :: i, n_nests
    real    :: t_val, t_val2, t_val3

    n_nests = options(1)%general%nests

    t_val = 0

    do i = 1, n_nests
        associate (comp => component(i)%comp)
        select type (comp)
            type is (domain_t)        
                t_val = t_val + comp%total_timer%mean(comp%compute_comms)            
            end select
        end associate
    enddo

    !
    !-----------------------------------------
    if (STD_OUT_PE) then
        call MPI_Comm_Size(MPI_COMM_WORLD,i)


        write(*,'(/ A)') "------------------------------------------------------"
        write(*,'(A)')   "Simulation completed successfully!"
        write(*,'(A /)') "------------------------------------------------------"
        write(*,*) "Model run from : ",trim(as_string(options(1)%general%start_time))
        write(*,*) "           to  : ",trim(as_string(options(1)%general%end_time))
        write(*,*) "Number of images:",i
        write(*,*) ""
        write(*,*) "Timing across all compute images:"
        write(*,*) ""
        write(*,*) "  Total time: ",t_val
        write(*,*) ""
        flush(output_unit)

    endif

    do i = 1, n_nests
        associate (comp => component(i)%comp)
        select type (comp)
            type is (domain_t)
                if (STD_OUT_PE) write(*,"(A22,I2,A2,A,A16)") "-------------- Domain ",i," (",trim(options(i)%domain%init_conditions_file),") --------------"
                if (STD_OUT_PE) write(*,'(A31 A10 A3 A10 A3 A10)') " ", "mean", " | ", "min", " | ", "max"
                t_val = comp%total_timer%mean(comp%compute_comms)
                t_val2 = comp%total_timer%min(comp%compute_comms)
                t_val3 = comp%total_timer%max(comp%compute_comms)
                if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "total", ":", t_val, " | ", t_val2, " | ", t_val3
                t_val = comp%initialization_timer%mean(comp%compute_comms)
                t_val2 = comp%initialization_timer%min(comp%compute_comms)
                t_val3 = comp%initialization_timer%max(comp%compute_comms)
                if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "init", ":", t_val, " | ", t_val2, " | ", t_val3
                t_val = comp%input_timer%mean(comp%compute_comms)
                t_val2 = comp%input_timer%min(comp%compute_comms)
                t_val3 = comp%input_timer%max(comp%compute_comms)
                if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "input", ":", t_val, " | ", t_val2, " | ", t_val3
                t_val = comp%output_timer%mean(comp%compute_comms)
                t_val2 = comp%output_timer%min(comp%compute_comms)
                t_val3 = comp%output_timer%max(comp%compute_comms)
                if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "output", ":", t_val, " | ", t_val2, " | ", t_val3
                t_val = comp%physics_timer%mean(comp%compute_comms)
                t_val2 = comp%physics_timer%min(comp%compute_comms)
                t_val3 = comp%physics_timer%max(comp%compute_comms)
                if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "physics", ":", t_val, " | ", t_val2, " | ", t_val3
                t_val = comp%mp_timer%mean(comp%compute_comms)
                t_val2 = comp%mp_timer%min(comp%compute_comms)
                t_val3 = comp%mp_timer%max(comp%compute_comms)
                if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "microphysics", ":", t_val, " | ", t_val2, " | ", t_val3
                t_val = comp%adv_timer%mean(comp%compute_comms)
                t_val2 = comp%adv_timer%min(comp%compute_comms)
                t_val3 = comp%adv_timer%max(comp%compute_comms)
                if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "advection", ":", t_val, " | ", t_val2, " | ", t_val3
                t_val = comp%sum_timer%mean(comp%compute_comms)
                t_val2 = comp%sum_timer%min(comp%compute_comms)
                t_val3 = comp%sum_timer%max(comp%compute_comms)
                if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "advection_sum", ":", t_val, " | ", t_val2, " | ", t_val3
                t_val = comp%flux_timer%mean(comp%compute_comms)
                t_val2 = comp%flux_timer%min(comp%compute_comms)
                t_val3 = comp%flux_timer%max(comp%compute_comms)
                if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "advection_flux", ":", t_val, " | ", t_val2, " | ", t_val3
                t_val = comp%flux_up_timer%mean(comp%compute_comms)
                t_val2 = comp%flux_up_timer%min(comp%compute_comms)
                t_val3 = comp%flux_up_timer%max(comp%compute_comms)
                if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "advection_flux_up", ":", t_val, " | ", t_val2, " | ", t_val3
                t_val = comp%flux_corr_timer%mean(comp%compute_comms)
                t_val2 = comp%flux_corr_timer%min(comp%compute_comms)
                t_val3 = comp%flux_corr_timer%max(comp%compute_comms)
                if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "advection_flux_corr", ":", t_val, " | ", t_val2, " | ", t_val3
                t_val = comp%adv_wind_timer%mean(comp%compute_comms)
                t_val2 = comp%adv_wind_timer%min(comp%compute_comms)
                t_val3 = comp%adv_wind_timer%max(comp%compute_comms)
                if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "advection_wind", ":", t_val, " | ", t_val2, " | ", t_val3
                t_val = comp%rad_timer%mean(comp%compute_comms)
                t_val2 = comp%rad_timer%min(comp%compute_comms)
                t_val3 = comp%rad_timer%max(comp%compute_comms)
                if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "radiation", ":", t_val, " | ", t_val2, " | ", t_val3
                t_val = comp%lsm_timer%mean(comp%compute_comms)
                t_val2 = comp%lsm_timer%min(comp%compute_comms)
                t_val3 = comp%lsm_timer%max(comp%compute_comms)
                if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "LSM", ":", t_val, " | ", t_val2, " | ", t_val3
                t_val = comp%pbl_timer%mean(comp%compute_comms)
                t_val2 = comp%pbl_timer%min(comp%compute_comms)
                t_val3 = comp%pbl_timer%max(comp%compute_comms)
                if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "PBL", ":", t_val, " | ", t_val2, " | ", t_val3
                t_val = comp%forcing_timer%mean(comp%compute_comms)
                t_val2 = comp%forcing_timer%min(comp%compute_comms)
                t_val3 = comp%forcing_timer%max(comp%compute_comms)
                if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "forcing", ":", t_val, " | ", t_val2, " | ", t_val3
                t_val = comp%wind_bal_timer%mean(comp%compute_comms)
                t_val2 = comp%wind_bal_timer%min(comp%compute_comms)
                t_val3 = comp%wind_bal_timer%max(comp%compute_comms)
                if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "wind bal", ":", t_val, " | ", t_val2, " | ", t_val3
                t_val = comp%diagnostic_timer%mean(comp%compute_comms)
                t_val2 = comp%diagnostic_timer%min(comp%compute_comms)
                t_val3 = comp%diagnostic_timer%max(comp%compute_comms)
                if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "diagnostic", ":", t_val, " | ", t_val2, " | ", t_val3
                t_val = comp%send_timer%mean(comp%compute_comms)
                t_val2 = comp%send_timer%min(comp%compute_comms)
                t_val3 = comp%send_timer%max(comp%compute_comms)
                if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "halo-exchange(send)", ":", t_val, " | ", t_val2, " | ", t_val3
                t_val = comp%ret_timer%mean(comp%compute_comms)
                t_val2 = comp%ret_timer%min(comp%compute_comms)
                t_val3 = comp%ret_timer%max(comp%compute_comms)
                if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "halo-exchange(retrieve)", ":", t_val, " | ", t_val2, " | ", t_val3
                t_val = comp%wait_timer%mean(comp%compute_comms)
                t_val2 = comp%wait_timer%min(comp%compute_comms)
                t_val3 = comp%wait_timer%max(comp%compute_comms)
                if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "halo-exchange(wait)", ":", t_val, " | ", t_val2, " | ", t_val3
                t_val = comp%wind_timer%mean(comp%compute_comms)
                t_val2 = comp%wind_timer%min(comp%compute_comms)
                t_val3 = comp%wind_timer%max(comp%compute_comms)
                if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "winds", ":", t_val, " | ", t_val2, " | ", t_val3
#ifdef _OPENACC
                t_val = comp%cpu_gpu_timer%mean(comp%compute_comms)
                t_val2 = comp%cpu_gpu_timer%min(comp%compute_comms)
                t_val3 = comp%cpu_gpu_timer%max(comp%compute_comms)
                if (STD_OUT_PE) write(*,'(A30 A1 F10.3 A3 F10.3 A3 F10.3)') "CPU-GPU transfers", ":", t_val, " | ", t_val2, " | ", t_val3
#endif
                if (STD_OUT_PE) flush(output_unit)
                if ( ANY(options%physics%windtype == kITERATIVE_WINDS) .or. ANY(options%physics%windtype == kLINEAR_ITERATIVE_WINDS)) then
                    call finalize_petsc()
                endif
            type is (ioserver_t)
                if (i==1) call comp%close_files()
            class default
        end select
        end associate
    end do

end subroutine component_program_end

end module flow_events
