!!! test_halo_exch.F90
! Purpose: Test the halo exchange functionality
! This program initializes a grid and variables, performs halo exchanges, and verifies the results.
! It checks the correctness of the halo exchange by comparing the exchanged values with expected values.
! 
module test_halo_exch

    use variable_dict_interface, only: var_dict_t
    use variable_interface,      only: variable_t
    use mpi_f08
    use grid_interface, only: grid_t
    use halo_interface, only: halo_t
    use icar_constants
    use testdrive, only : new_unittest, unittest_type, error_type, check, test_failed

    implicit none
    private
    
    public :: collect_halo_exch_suite
    
    contains

    !> Collect all exported unit tests
    subroutine collect_halo_exch_suite(testsuite)
        !> Collection of tests
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
      
        testsuite = [ &
          new_unittest("batch_exch", test_batch_exch), &
            new_unittest("var_exch", test_var_exch), &
            new_unittest("var_u_exch", test_var_u_exch), &
            new_unittest("var_v_exch", test_var_v_exch) &
          ]
      
    end subroutine collect_halo_exch_suite


    
    subroutine test_batch_exch(error)
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

        call halo_exch_standard(grid,error,batch_in=.True.,test_str_in="batch exchange")

    end subroutine test_batch_exch

    subroutine test_var_exch(error)
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

        call halo_exch_standard(grid,error,test_str_in="var exchange")
        if (allocated(error)) return
        call halo_exch_standard(grid,error,corners_in=.True.,test_str_in="var exchange, corners")

    end subroutine test_var_exch

    subroutine test_var_u_exch(error)
        type(error_type), allocatable, intent(out) :: error

        integer :: my_index, ierr
        type(grid_t) :: grid
        type(MPI_Comm) :: comms

        call MPI_Comm_Rank(MPI_COMM_WORLD,my_index,ierr)
        my_index = my_index + 1

        !set the comm MPI_com type to be MPI_COMM_WORLD
        CALL MPI_Comm_dup( MPI_COMM_WORLD, comms, ierr )

        !initialize grids
        call grid%set_grid_dimensions( 20, 20, 20, image=my_index, comms=comms, adv_order=3, nx_extra = 1)

        call halo_exch_standard(grid,error,test_str_in="var staggered on u grid")
        if (allocated(error)) return
        call halo_exch_standard(grid,error,corners_in=.True.,test_str_in="var staggered on u grid, corners")

    end subroutine test_var_u_exch

    subroutine test_var_v_exch(error)
        type(error_type), allocatable, intent(out) :: error

        integer :: my_index, ierr
        type(grid_t) :: grid
        type(MPI_Comm) :: comms

        call MPI_Comm_Rank(MPI_COMM_WORLD,my_index,ierr)
        my_index = my_index + 1

        !set the comm MPI_com type to be MPI_COMM_WORLD
        CALL MPI_Comm_dup( MPI_COMM_WORLD, comms, ierr )

        !initialize grids
        call grid%set_grid_dimensions( 20, 20, 20, image=my_index, comms=comms, adv_order=3, ny_extra = 1)

        call halo_exch_standard(grid,error,test_str_in="var staggered on v grid")
        if (allocated(error)) return
        call halo_exch_standard(grid,error,corners_in=.True.,test_str_in="var staggered on v grid, corners")

    end subroutine test_var_v_exch

    subroutine halo_exch_standard(grid,error,batch_in,corners_in,test_str_in)
        implicit none
        type(grid_t), intent(in) :: grid
        type(error_type), allocatable, intent(out) :: error
        logical, optional, intent(in) :: batch_in
        logical, optional, intent(in) :: corners_in
        character(len=*), optional, intent(in) :: test_str_in

        type(halo_t) :: halo
        type(var_dict_t) :: exch_vars, adv_vars
        type(variable_t) :: var
        integer :: my_index
        integer :: ierr
        type(MPI_Comm) :: comms
        logical :: batch, corners
        character(len=100) :: test_str
    
        batch = .false.
        corners = .false.
        test_str = ""
        if (present(batch_in)) batch = batch_in
        if (present(corners_in)) corners = corners_in
        if (present(test_str_in)) test_str = test_str_in

        call MPI_Comm_Rank(MPI_COMM_WORLD,my_index,ierr)
        my_index = my_index + 1

        !set the comm MPI_com type to be MPI_COMM_WORLD
        CALL MPI_Comm_dup( MPI_COMM_WORLD, comms, ierr )

        !initialize variables to exchange
        call var%initialize(grid)

        !populate adv_vars with two test variables
        call adv_vars%add_var('var', var) 
        !call adv_vars%add_var('temperature', temperature) 

        call halo%init(exch_vars, adv_vars, grid, comms)
        
        ! Initialize fields with my_index values
        var%data_3d = my_index

        if (batch) then
            call halo%batch_exch(exch_vars, adv_vars)
        else
            call halo%exch_var(var, corners=corners)
        endif

        ! Verify exchange worked
        ! check if this image is not on the eastern boundary
        if (.not.(grid%ximg == grid%ximages)) then
            !check that the eastern halo is filled with the value of my_index for the eastern neighbor
            if (.not.(ALL(var%data_3d((grid%ite+1):grid%ime,:,grid%jts:grid%jte) == my_index+1))) then
                call test_failed(error, "Halo exch failed", "Failed for eastern halo exchange, "//trim(test_str))
                return
            endif
        endif

        ! check if this image is not on the western boundary
        if (.not.(grid%ximg == 1)) then
            !check that the western halo is filled with the value of my_index for the western neighbor
            if (.not.(ALL(var%data_3d(grid%ims:grid%its-1,:,grid%jts:grid%jte) == my_index-1))) then
                call test_failed(error, "Halo exch failed", "Failed for western halo exchange, "//trim(test_str))
                return
            endif
        endif

        ! check if this image is not on the southern boundary
        if (.not.(grid%yimg == 1)) then
            !check that the southern halo is filled with the value of my_index for the southern neighbor
            if (.not.(ALL(var%data_3d(grid%its:grid%ite,:,grid%jms:grid%jts-1) == my_index-grid%ximages))) then
                call test_failed(error, "Halo exch failed", "Failed for southern halo exchange, "//trim(test_str))
                return
            endif
        endif

        ! check if this image is not on the northern boundary
        if (.not.(grid%yimg == grid%yimages)) then
            !check that the northern halo is filled with the value of my_index for the northern neighbor
            if (.not.(ALL(var%data_3d(grid%its:grid%ite,:,grid%jte+1:grid%jme) == my_index+grid%ximages))) then
                call test_failed(error, "Halo exch failed", "Failed for northern halo exchange, "//trim(test_str))
                return
            endif
        endif

        if (.not.(corners)) return

        ! if this image is in the north eastern corner
        if ((grid%yimg == 1 .and. grid%ximg == 1)) then
            !check that the north eastern corner halo is filled with the value of my_index for the north eastern neighbor
            if (.not.(ALL(var%data_3d(grid%ite+1:grid%ime,:,grid%jte+1:grid%jme) == my_index+grid%ximages+1))) then
                call test_failed(error, "Halo exch failed", "Failed for north eastern halo exchange, "//trim(test_str))
                return
            endif
        endif

        ! if this image is in the north western corner
        if ((grid%yimg == 1 .and. grid%ximg == grid%ximages)) then
            !check that the north western corner halo is filled with the value of my_index for the north western neighbor
            if (.not.(ALL(var%data_3d(grid%ims:grid%its-1,:,grid%jte+1:grid%jme) == my_index+grid%ximages-1))) then
                call test_failed(error, "Halo exch failed", "Failed for north western halo exchange, "//trim(test_str))
                return
            endif
        endif

        ! if this image is in the south eastern corner
        if ((grid%yimg == grid%yimages .and. grid%ximg == 1)) then
            !check that the south eastern corner halo is filled with the value of my_index for the south eastern neighbor
            if (.not.(ALL(var%data_3d(grid%ite+1:grid%ime,:,grid%jms:grid%jts-1) == my_index-grid%ximages+1))) then
                call test_failed(error, "Halo exch failed", "Failed for south eastern halo exchange, "//trim(test_str))
                return
                stop
            endif
        endif

        ! if this image is in the south western corner
        if ((grid%yimg == grid%yimages .and. grid%ximg == grid%ximages)) then
            !check that the south western corner halo is filled with the value of my_index for the south western neighbor
            if (.not.(ALL(var%data_3d(grid%ims:grid%its-1,:,grid%jms:grid%jts-1) == my_index-grid%ximages-1))) then
                call test_failed(error, "Halo exch failed", "Failed for south western halo exchange, "//trim(test_str))
                return
            endif
        endif
    end subroutine halo_exch_standard

end module test_halo_exch
