! test_advect.F90
! This program initializes a grid and sets up MPI communication for a simple advection test.
! It tests all possible configurations of advection, of varying numerical orders, and use of the monotonic flux limiter

module test_advect
          
    use variable_dict_interface, only: var_dict_t
    use variable_interface,      only: variable_t
    use mpi_f08
    use grid_interface, only: grid_t
    use halo_interface, only: halo_t
    use icar_constants
    use testdrive, only : new_unittest, unittest_type, error_type, check, test_failed

    implicit none
    private
    
    ! Hard coded time step factor to be used with RK3 time stepping for 5th and 3rd order advection
    real :: time_step_factor(2) = (/ 1.6, 1.4/)
    integer :: h_order(3) = (/ 1, 3, 5/)
    integer :: v_order(3) = (/ 1, 3, 5/)
    integer :: i, j, k, l

    public :: collect_advect_suite
    
    contains
    


    !> Collect all exported unit tests
    subroutine collect_advect_suite(testsuite)
        !> Collection of tests
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
        
        testsuite = [ &
            new_unittest("adv_h1_v1", test_adv_h1_v1) &
            ! new_unittest("adv_h3_v3", test_adv_h3_v3), &
            ! new_unittest("adv_h3_v3", test_adv_h5_v3), &
            ! new_unittest("adv_h3_v3", test_adv_h5_v5), &
            ]
    
    end subroutine collect_advect_suite
    
    subroutine test_adv_h1_v1(error)
        type(error_type), allocatable, intent(out) :: error

        integer :: my_index, ierr
        type(grid_t) :: grid
        type(MPI_Comm) :: comms

        call MPI_Comm_Rank(MPI_COMM_WORLD,my_index,ierr)
        my_index = my_index + 1

        !set the comm MPI_com type to be MPI_COMM_WORLD
        CALL MPI_Comm_dup( MPI_COMM_WORLD, comms, ierr )

        call advect_standard(h_order=1,v_order=1,error=error)

    end subroutine test_adv_h1_v1
    
    subroutine advect_standard(h_order,v_order,error)
        integer, intent(in) :: h_order,v_order
        type(error_type), allocatable, intent(out) :: error

        integer :: my_index, ierr
        type(grid_t) :: grid
        type(MPI_Comm) :: comms

        call MPI_Comm_Rank(MPI_COMM_WORLD,my_index,ierr)
        my_index = my_index + 1

        !set the comm MPI_com type to be MPI_COMM_WORLD
        CALL MPI_Comm_dup( MPI_COMM_WORLD, comms, ierr )

        !initialize grids
        call grid%set_grid_dimensions( 20, 20, 20, image=my_index, comms=comms, adv_order=3)

    end subroutine advect_standard


end module test_advect