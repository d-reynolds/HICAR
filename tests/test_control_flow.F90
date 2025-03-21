!!! test_control_flow.F90
! Purpose: Test the control loop of the main driver program.
! 
! Relies on a single IO process and a single compute process, which essentially ping-pong
! between the two processes, seeing if we can get them to stall.
module test_control_flow

    use mpi_f08
    use options_interface, only: options_t
    use time_object, only: Time_type
    use time_delta_object, only: time_delta_t
    use flow_object_interface, only: flow_obj_t
    use icar_constants
    use nest_manager, only: all_nests_not_done, nest_next_up, should_update_nests, &
        can_update_child_nest
    use testdrive, only : new_unittest, unittest_type, error_type, check, test_failed

    implicit none
    private
    
    public :: collect_control_flow_suite
    
    contains

    !> Collect all exported unit tests
    subroutine collect_control_flow_suite(testsuite)
        !> Collection of tests
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
      
        testsuite = [ &
            new_unittest("control_flow", test_ctrl_flow) &
          ]
      
    end subroutine collect_control_flow_suite

    !> Test the control flow of the main driver program
    subroutine test_ctrl_flow(error)
        type(error_type), allocatable, intent(out) :: error

        ! Split the processes into two teams based on even or odd rank
        type(MPI_COMM) :: exec_team, io_team
        integer :: num_procs, my_rank
        integer :: ierr, n_nests, i, my_io_rank, my_exec_rank
        integer :: count, n
        type(flow_obj_t), allocatable :: flow_obj(:)
        type(options_t), allocatable :: options(:)

        ! Get the number of processes and the rank of this process
        call MPI_Comm_size(MPI_COMM_WORLD, num_procs, ierr)
        call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
        ! Create the IO team
        call MPI_Comm_split(MPI_COMM_WORLD, my_rank / 2, my_rank, io_team, ierr)
        ! Create the compute team
        call MPI_Comm_split(MPI_COMM_WORLD, my_rank - (my_rank / 2) * 2, my_rank, exec_team, ierr)
        ! Get the rank of this process in the IO team
        call MPI_Comm_rank(io_team, my_io_rank, ierr)
        ! Get the rank of this process in the compute team
        call MPI_Comm_rank(exec_team, my_exec_rank, ierr)

        n_nests = 2

        allocate(flow_obj(n_nests))
        allocate(options(n_nests))
        ! initialize the objects
        do i = 1, n_nests
            call options(i)%init()
            options(i)%general%nests = n_nests
            ! options(i)%general%nest_indx = i
            call options(i)%general%start_time%init(options(i)%general%calendar)
            call options(i)%general%start_time%set("2020-01-01 00:00:00")
            call options(i)%general%end_time%init(options(i)%general%calendar)
            call options(i)%general%end_time%set("2020-01-02 00:00:00")

            call options(i)%forcing%input_dt%set(seconds=3600)
            call options(i)%output%output_dt%set(seconds=3600)

            if (i < n_nests) then
                options(i)%general%parent_nest = i - 1
                allocate(options(i)%general%child_nests(0))
            else
                options(i)%general%parent_nest = i - 1
                allocate(options(i)%general%child_nests(1))
                options(i)%general%child_nests = [i]
            end if

    
            call flow_obj(i)%init_flow_obj(options(i), i)
        end do

        if (my_exec_rank >= 0) then

            do i = 1, n_nests
                call flow_obj(i)%increment_input_time()
                call flow_obj(i)%increment_output_time()
            enddo
            !> This is the compute team, so we will run the compute loop
            do while (all_nests_not_done(flow_obj(1:n_nests)))
                do i = 1, n_nests
                    ! -----------------------------------------------------------
                    !  Initialization of child nests, or cycling a completed nest
                    ! -----------------------------------------------------------
                    ! If we are a child nest, and the domain has not yet been initialized with data...
                    if (nest_next_up(flow_obj,options(i), i)) then
                        ! ...then initialize ourselves here
                        call MPI_Barrier(MPI_COMM_WORLD, ierr)
                        call flow_obj(i)%increment_output_time()
                    endif
    
                    if (flow_obj(i)%dead_or_asleep()) cycle
        
                    call MPI_Barrier(MPI_COMM_WORLD, ierr)
                    call flow_obj(i)%increment_input_time()

                    do while ( .not.(flow_obj(i)%time_for_input()) .and. .not.(flow_obj(i)%ended) )
    
                        call flow_obj(i)%set_sim_time(flow_obj(i)%next_flow_event())

                        if (flow_obj(i)%time_for_output()) then
                            call MPI_Barrier(MPI_COMM_WORLD, ierr)
                            call flow_obj(i)%increment_output_time()
                        endif
                    enddo
    
                    ! If we have children, who are not yet done, then we need to push our data to them
                    if (should_update_nests(flow_obj, options(i), i)) then
                        call MPI_Barrier(MPI_COMM_WORLD, ierr)
                    endif
                enddo
            enddo
        else if (my_io_rank >= 0) then

            do i = 1, n_nests
                call flow_obj(i)%increment_input_time()
                call flow_obj(i)%increment_output_time()
            enddo
            !> This is the IO team, so we will run the IO loop
            do i = 1, n_nests
                if (options(i)%general%parent_nest == 0) call flow_obj(i)%increment_input_time()
            end do
            do while (all_nests_not_done(flow_obj(1:n_nests)))
                do i = 1, n_nests
                    if (nest_next_up(flow_obj,options(i),i)) then
                        ! ...then initialize ourselves here
                        call MPI_Barrier(MPI_COMM_WORLD, ierr)
                        call flow_obj(i)%increment_output_time()
                    endif
    
                    if (flow_obj(i)%dead_or_asleep()) cycle

                    call MPI_Barrier(MPI_COMM_WORLD, ierr)
                    call flow_obj(i)%increment_input_time()

                    do while ( .not.(flow_obj(i)%time_for_input()) .and. .not.(flow_obj(i)%ended) )

                        call flow_obj(i)%set_sim_time(flow_obj(i)%next_flow_event())
    
                        ! now loop over all child nests
                        if (should_update_nests(flow_obj,options(i),i)) then
                            call MPI_Barrier(MPI_COMM_WORLD, ierr)
                            do n = 1, size(options(i)%general%child_nests)
                                !Test if we can update the child nest
                                if ( can_update_child_nest(flow_obj(i),flow_obj(options(i)%general%child_nests(n))) ) then
                                    ! This call will distribute the model state of the forcing fields to the child nest
                                    call MPI_Barrier(io_team, ierr)
                                endif
                            enddo
                        endif
                        if (flow_obj(i)%time_for_output()) then
                            call MPI_Barrier(MPI_COMM_WORLD, ierr)
                            call flow_obj(i)%increment_output_time()
                        endif
                    enddo
                enddo
            enddo
        end if

        write(*,*) "All nests are done."

    end subroutine test_ctrl_flow

end module test_control_flow
