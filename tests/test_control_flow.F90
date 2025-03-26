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
    use ioclient_interface, only: ioclient_t
    use boundary_interface, only: boundary_t
    use flow_events,        only: component_loop
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
            new_unittest("standard", test_standard_prgm), &
            new_unittest("standard_restart", test_standard_restart_prgm), &
            new_unittest("nested", test_nested_prgm), &
            new_unittest("nested_restart", test_nested_restart_prgm) &
          ]
      
    end subroutine collect_control_flow_suite

    !> Test the control flow of the main driver program
    subroutine test_standard_prgm(error)
        type(error_type), allocatable, intent(out) :: error

        type(options_t), allocatable :: in_options(:)
        integer :: i, n_nests

        n_nests = 1
        allocate(in_options(n_nests))

        do i = 1, n_nests
            in_options(i)%general%calendar = "gregorian"
            call in_options(i)%general%start_time%init(in_options(i)%general%calendar)
            call in_options(i)%general%start_time%set("2000-01-01 00:00:00")
            call in_options(i)%general%end_time%init(in_options(i)%general%calendar)
            call in_options(i)%general%end_time%set("2000-01-02 00:00:00")
            call in_options(i)%forcing%input_dt%set(seconds=3600)
            call in_options(i)%output%output_dt%set(seconds=3600)
            in_options(i)%general%parent_nest = i - 1
            if (i<n_nests) then
                allocate(in_options(i)%general%child_nests(1))
                in_options(i)%general%child_nests = [i+1]
            else
                allocate(in_options(i)%general%child_nests(0))
            end if
        enddo

        call test_main_prg(in_options)

    end subroutine test_standard_prgm

    !> Test the control flow of the main driver program
    subroutine test_standard_restart_prgm(error)
        type(error_type), allocatable, intent(out) :: error

        type(options_t), allocatable :: in_options(:)
        integer :: i, n_nests

        n_nests = 1
        allocate(in_options(n_nests))

        do i = 1, n_nests
            in_options(i)%general%calendar = "gregorian"
            call in_options(i)%general%start_time%init(in_options(i)%general%calendar)
            call in_options(i)%general%start_time%set("2000-01-01 00:00:00")
            call in_options(i)%general%end_time%init(in_options(i)%general%calendar)
            call in_options(i)%general%end_time%set("2000-01-02 00:00:00")
            in_options(i)%restart%restart = .true.
            call in_options(i)%restart%restart_time%init(in_options(i)%general%calendar)
            call in_options(i)%restart%restart_time%set("2000-01-01 03:30:00")
            call in_options(i)%forcing%input_dt%set(seconds=3600)
            call in_options(i)%output%output_dt%set(seconds=600)
            in_options(i)%general%parent_nest = i - 1
            if (i<n_nests) then
                allocate(in_options(i)%general%child_nests(1))
                in_options(i)%general%child_nests = [i+1]
            else
                allocate(in_options(i)%general%child_nests(0))
            end if
        enddo

        call test_main_prg(in_options)

    end subroutine test_standard_restart_prgm

    !> Test the control flow of the main driver program
    subroutine test_nested_prgm(error)
        type(error_type), allocatable, intent(out) :: error

        type(options_t), allocatable :: in_options(:)
        integer :: i, n_nests

        n_nests = 3
        allocate(in_options(n_nests))

        do i = 1, n_nests
            in_options(i)%general%calendar = "gregorian"
            call in_options(i)%general%start_time%init(in_options(i)%general%calendar)
            call in_options(i)%general%start_time%set("2000-01-01 00:00:00")
            call in_options(i)%general%end_time%init(in_options(i)%general%calendar)
            call in_options(i)%general%end_time%set("2000-01-01 00:30:00")
            call in_options(i)%forcing%input_dt%set(seconds=3600)
            call in_options(i)%output%output_dt%set(seconds=600)
            in_options(i)%general%parent_nest = i - 1
            if (i<n_nests) then
                allocate(in_options(i)%general%child_nests(1))
                in_options(i)%general%child_nests = [i+1]
            else
                allocate(in_options(i)%general%child_nests(0))
            end if
        enddo

        call test_main_prg(in_options)

    end subroutine test_nested_prgm

    !> Test the control flow of the main driver program
    subroutine test_nested_restart_prgm(error)
        type(error_type), allocatable, intent(out) :: error

        type(options_t), allocatable :: in_options(:)
        integer :: i, n_nests

        n_nests = 3
        allocate(in_options(n_nests))

        do i = 1, n_nests
            in_options(i)%general%calendar = "gregorian"
            call in_options(i)%general%start_time%init(in_options(i)%general%calendar)
            call in_options(i)%general%start_time%set("2000-01-01 00:00:00")
            call in_options(i)%general%end_time%init(in_options(i)%general%calendar)
            call in_options(i)%general%end_time%set("2000-01-02 00:00:00")
            in_options(i)%restart%restart = .true.
            call in_options(i)%restart%restart_time%init(in_options(i)%general%calendar)
            call in_options(i)%restart%restart_time%set("2000-01-01 05:00:00")
            call in_options(i)%forcing%input_dt%set(seconds=3600)
            call in_options(i)%output%output_dt%set(seconds=3600)
            in_options(i)%general%parent_nest = i - 1
            if (i<n_nests) then
                allocate(in_options(i)%general%child_nests(1))
                in_options(i)%general%child_nests = [i+1]
            else
                allocate(in_options(i)%general%child_nests(0))
            end if
        enddo

        call test_main_prg(in_options)

    end subroutine test_nested_restart_prgm

    subroutine test_main_prg(in_options)
        type(options_t), intent(in) :: in_options(:)

        integer :: n_nests, i
        logical :: is_io, is_exec
        integer :: my_rank, num_procs, my_io_rank, my_exec_rank
        integer :: ierr
        type(flow_obj_t), allocatable :: flow_obj(:)
        type(options_t), allocatable :: options(:)
        type(boundary_t), allocatable :: boundary(:)
        type(ioclient_t), allocatable :: ioclient(:)
        type(MPI_Comm) :: io_team, exec_team

        n_nests = size(in_options)

        allocate(flow_obj(n_nests))
        allocate(options(n_nests))
        allocate(boundary(n_nests))
        allocate(ioclient(n_nests))

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

        if (my_io_rank > 0) then
            ! This is an IO process
            ioclient%parent_comms = MPI_COMM_NULL
        else
            ! This is a compute process
            ioclient%parent_comms = MPI_COMM_WORLD
        end if

        write (*,*) "entering main prog"

        ! initialize the objects
        do i = 1, n_nests
            call options(i)%init()
            write (*,*) "options init"

            options(i)%general%nests = n_nests
            options(i)%nest_indx = i
            call options(i)%general%start_time%init(in_options(i)%general%calendar)
            options(i)%general%start_time = in_options(i)%general%start_time
            call options(i)%general%end_time%init(in_options(i)%general%calendar)
            options(i)%general%end_time = in_options(i)%general%end_time

            options(i)%restart%restart = in_options(i)%restart%restart
            options(i)%restart%restart_time = in_options(i)%restart%restart_time
            call options(i)%forcing%input_dt%set(seconds=in_options(i)%forcing%input_dt%seconds())
            call options(i)%output%output_dt%set(seconds=in_options(i)%output%output_dt%seconds())


            options(i)%general%parent_nest = in_options(i)%general%parent_nest

            allocate(options(i)%general%child_nests(size(in_options(i)%general%child_nests)))

            if (size(in_options(i)%general%child_nests) > 0) then
                options(i)%general%child_nests = in_options(i)%general%child_nests
            end if
            write (*,*) "initializing flow obj"

            call flow_obj(i)%init_flow_obj(options(i), i)
        end do
        write (*,*) "Flow objects initialized"
        call component_loop(flow_obj(1:n_nests), options(1:n_nests), boundary(1:n_nests), ioclient(1:n_nests))
        write (*,*) "Finished test prg for rank ", my_rank

        ! Ensure that all processes have really finished the main loop, and sync is completed
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
    end subroutine test_main_prg


end module test_control_flow
