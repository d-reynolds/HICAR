!!! test_halo_exch.F90
! Purpose: Test the halo exchange functionality
! This program initializes a grid and variables, performs halo exchanges, and verifies the results.
! It checks the correctness of the halo exchange by comparing the exchanged values with expected values.
! 
module test_halo_exch

    use variable_interface,      only: variable_t
    use mpi_f08
    use grid_interface, only: grid_t
    use halo_interface, only: halo_t
    use icar_constants
    use data_structures, only: index_type
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
            new_unittest("var_exch_2d", test_var_exch_2d), &
            new_unittest("var_u_exch", test_var_u_exch), &
            new_unittest("var_v_exch", test_var_v_exch) &
          ]
      
    end subroutine collect_halo_exch_suite


    
    subroutine test_batch_exch(error)
        type(error_type), allocatable, intent(out) :: error

        integer :: my_index, ierr
        type(grid_t) :: grid, grid_2d
        type(MPI_Comm) :: comms

        call MPI_Comm_Rank(MPI_COMM_WORLD,my_index,ierr)
        my_index = my_index + 1

        !set the comm MPI_com type to be MPI_COMM_WORLD
        CALL MPI_Comm_dup( MPI_COMM_WORLD, comms, ierr )

        !initialize grids
        call grid%set_grid_dimensions( 171, 517, 20, image=my_index, comms=comms, adv_order=3)

        call halo_exch_standard(grid,error,batch_in=.True.,test_str_in="batch exchange")

        call grid_2d%set_grid_dimensions( 171, 517, 0, image=my_index, global_nz=20, comms=comms, adv_order=3)

        call halo_exch_standard(grid_2d,error,batch_in=.True.,test_str_in="batch exchange 2d")

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
        call grid%set_grid_dimensions( 171, 517, 20, image=my_index, comms=comms, adv_order=3)

        call halo_exch_standard(grid,error,test_str_in="var exchange")
        if (allocated(error)) return
        call halo_exch_standard(grid,error,corners_in=.True.,test_str_in="var exchange, corners")
        call halo_exch_standard(grid,error,corners_in=.True.,do_dqdt=.True.,test_str_in="var exchange, corners, dqdt")

    end subroutine test_var_exch

    subroutine test_var_exch_2d(error)
        type(error_type), allocatable, intent(out) :: error

        integer :: my_index, ierr
        type(grid_t) :: grid
        type(MPI_Comm) :: comms

        call MPI_Comm_Rank(MPI_COMM_WORLD,my_index,ierr)
        my_index = my_index + 1

        !set the comm MPI_com type to be MPI_COMM_WORLD
        CALL MPI_Comm_dup( MPI_COMM_WORLD, comms, ierr )

        !initialize grids
        call grid%set_grid_dimensions( 171, 517, 0, image=my_index, global_nz=20, comms=comms, adv_order=3)

        call halo_exch_standard(grid,error,test_str_in="var exchange")
        if (allocated(error)) return
        call halo_exch_standard(grid,error,corners_in=.True.,test_str_in="var exchange, corners")

    end subroutine test_var_exch_2d


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
        call grid%set_grid_dimensions( 171, 517, 20, image=my_index, comms=comms, adv_order=3, nx_extra = 1)

        call halo_exch_standard(grid,error,test_str_in="var staggered on u grid")
        if (allocated(error)) return
        call halo_exch_standard(grid,error,corners_in=.True.,test_str_in="var staggered on u grid, corners")
        call halo_exch_standard(grid,error,corners_in=.True.,do_dqdt=.True.,test_str_in="var staggered on u grid, corners, dqdt")

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
        call grid%set_grid_dimensions( 171, 517, 20, image=my_index, comms=comms, adv_order=3, ny_extra = 1)

        call halo_exch_standard(grid,error,test_str_in="var staggered on v grid")
        if (allocated(error)) return
        call halo_exch_standard(grid,error,corners_in=.True.,test_str_in="var staggered on v grid, corners")
        call halo_exch_standard(grid,error,corners_in=.True.,do_dqdt=.True.,test_str_in="var staggered on v grid, corners, dqdt")

    end subroutine test_var_v_exch

    subroutine halo_exch_standard(grid,error,batch_in,corners_in,do_dqdt,test_str_in)
        implicit none
        type(grid_t), intent(in) :: grid
        type(error_type), allocatable, intent(out) :: error
        logical, optional, intent(in) :: batch_in, corners_in, do_dqdt
        character(len=*), optional, intent(in) :: test_str_in

        type(halo_t) :: halo
        type(index_type), allocatable :: exch_vars(:)
        type(variable_t) :: var, var_data(1), exch_var
        integer :: my_index
        integer :: ierr
        type(MPI_Comm) :: comms
        type(grid_t) :: grid_3d
        logical :: batch, corners, interior, dqdt
        logical :: north, south, east, west
        logical :: northeast, northwest, southeast, southwest
        integer :: i, j, k
        real    :: val
        character(len=100) :: test_str
    
        batch = .false.
        corners = .false.
        dqdt = .False.
        test_str = ""
        if (present(batch_in)) batch = batch_in
        if (present(corners_in)) corners = corners_in
        if (present(do_dqdt)) dqdt = do_dqdt
        if (present(test_str_in)) test_str = test_str_in

        call MPI_Comm_Rank(MPI_COMM_WORLD,my_index,ierr)
        my_index = my_index + 1

        !set the comm MPI_com type to be MPI_COMM_WORLD
        CALL MPI_Comm_dup( MPI_COMM_WORLD, comms, ierr )

        if (grid%is2d) then

            call grid_3d%set_grid_dimensions( grid%nx_global, grid%ny_global, 20, image=my_index, comms=comms, adv_order=3)

            !initialize variables to exchange
            call var%initialize(kVARS%water_vapor,grid_3d,forcing_var=.True.)
            call exch_var%initialize(kVARS%snow_height,grid,forcing_var=.True.)

            allocate(exch_vars(2))

            exch_var%id = kVARS%snow_height
            exch_vars(1)%v = 1
            exch_vars(1)%id = kVARS%snow_height
            exch_vars(2)%v = 2
            exch_vars(2)%id = kVARS%water_vapor

        else

            !in case input grid has staggered dimensions, create a clean grid here with no stagger. This will be used to initialize the halo object
            call grid_3d%set_grid_dimensions( grid%nx_global-grid%nx_e, grid%ny_global-grid%ny_e, 20, image=my_index, comms=comms, adv_order=3)

            !initialize variables to exchange
            call var%initialize(kVARS%water_vapor,grid,forcing_var=dqdt)

            allocate(exch_vars(1))
            exch_vars(1)%v = 1
            exch_vars(1)%id = kVARS%water_vapor

        endif

        if (batch) then
            !initialize variables to exchange
            call var_data(1)%initialize(exch_vars(1)%id,grid,forcing_var=dqdt)
        endif

        call halo%init(exch_vars, grid_3d, comms)

        do i = grid%its, grid%ite
            do j = grid%jts, grid%jte
                do k = 1, max(1,grid%kte)
                    ! Set tile cells to global index values
                    if (grid%is3d) then
                        var%data_3d(i,k,j) = i+(j-1)*grid%nx_global+(k-1)*grid%nx_global*grid%ny_global
                    else
                        exch_var%data_2d(i,j) = i+(j-1)*grid%nx_global
                    endif
                end do
            end do
        end do
        ! Initialize boundary halo cells to their global index on sides where
        ! no neighbor exists to exchange with (global domain boundaries only).
        if (grid%ximg == 1) then
            do i = grid%ims, grid%its-1
                do j = grid%jts, grid%jte
                    do k = 1, max(1,grid%kte)
                        if (grid%is3d) then
                            var%data_3d(i,k,j) = i+(j-1)*grid%nx_global+(k-1)*grid%nx_global*grid%ny_global
                        else
                            exch_var%data_2d(i,j) = i+(j-1)*grid%nx_global
                        endif
                    end do
                end do
            end do
        endif
        if (grid%ximg == grid%ximages) then
            do i = grid%ite+1, grid%ime
                do j = grid%jts, grid%jte
                    do k = 1, max(1,grid%kte)
                        if (grid%is3d) then
                            var%data_3d(i,k,j) = i+(j-1)*grid%nx_global+(k-1)*grid%nx_global*grid%ny_global
                        else
                            exch_var%data_2d(i,j) = i+(j-1)*grid%nx_global
                        endif
                    end do
                end do
            end do
        endif
        if (grid%yimg == 1) then
            do i = grid%its, grid%ite
                do j = grid%jms, grid%jts-1
                    do k = 1, max(1,grid%kte)
                        if (grid%is3d) then
                            var%data_3d(i,k,j) = i+(j-1)*grid%nx_global+(k-1)*grid%nx_global*grid%ny_global
                        else
                            exch_var%data_2d(i,j) = i+(j-1)*grid%nx_global
                        endif
                    end do
                end do
            end do
        endif
        if (grid%yimg == grid%yimages) then
            do i = grid%its, grid%ite
                do j = grid%jte+1, grid%jme
                    do k = 1, max(1,grid%kte)
                        if (grid%is3d) then
                            var%data_3d(i,k,j) = i+(j-1)*grid%nx_global+(k-1)*grid%nx_global*grid%ny_global
                        else
                            exch_var%data_2d(i,j) = i+(j-1)*grid%nx_global
                        endif
                    end do
                end do
            end do
        endif
        ! Initialize corner halo cells on processes at two domain boundaries
        ! (southwest corner: ximg==1, yimg==1, etc.)
        if (grid%ximg == 1 .and. grid%yimg == 1) then
            do i = grid%ims, grid%its-1
                do j = grid%jms, grid%jts-1
                    do k = 1, max(1,grid%kte)
                        if (grid%is3d) then
                            var%data_3d(i,k,j) = i+(j-1)*grid%nx_global+(k-1)*grid%nx_global*grid%ny_global
                        else
                            exch_var%data_2d(i,j) = i+(j-1)*grid%nx_global
                        endif
                    end do
                end do
            end do
        endif
        if (grid%ximg == 1 .and. grid%yimg == grid%yimages) then
            do i = grid%ims, grid%its-1
                do j = grid%jte+1, grid%jme
                    do k = 1, max(1,grid%kte)
                        if (grid%is3d) then
                            var%data_3d(i,k,j) = i+(j-1)*grid%nx_global+(k-1)*grid%nx_global*grid%ny_global
                        else
                            exch_var%data_2d(i,j) = i+(j-1)*grid%nx_global
                        endif
                    end do
                end do
            end do
        endif
        if (grid%ximg == grid%ximages .and. grid%yimg == 1) then
            do i = grid%ite+1, grid%ime
                do j = grid%jms, grid%jts-1
                    do k = 1, max(1,grid%kte)
                        if (grid%is3d) then
                            var%data_3d(i,k,j) = i+(j-1)*grid%nx_global+(k-1)*grid%nx_global*grid%ny_global
                        else
                            exch_var%data_2d(i,j) = i+(j-1)*grid%nx_global
                        endif
                    end do
                end do
            end do
        endif
        if (grid%ximg == grid%ximages .and. grid%yimg == grid%yimages) then
            do i = grid%ite+1, grid%ime
                do j = grid%jte+1, grid%jme
                    do k = 1, max(1,grid%kte)
                        if (grid%is3d) then
                            var%data_3d(i,k,j) = i+(j-1)*grid%nx_global+(k-1)*grid%nx_global*grid%ny_global
                        else
                            exch_var%data_2d(i,j) = i+(j-1)*grid%nx_global
                        endif
                    end do
                end do
            end do
        endif
        if (batch .and. grid%is2d) var_data(1)%data_2d = exch_var%data_2d
        if (batch .and. grid%is3d) var_data(1)%data_3d = var%data_3d
        if (dqdt) var%dqdt_3d = var%data_3d

        !$acc data copy(halo, exch_vars)

        if (batch) then
            if (grid%is3d) then
                !$acc update device(var_data(1)%data_3d)
                call halo%halo_3d_send_batch(exch_vars, var_data)
                call halo%halo_3d_retrieve_batch(exch_vars, var_data)
                !$acc update host(var_data(1)%data_3d)
            else
                !$acc update device(var_data(1)%data_2d)
                call halo%halo_2d_send_batch(exch_vars, var_data)
                call halo%halo_2d_retrieve_batch(exch_vars, var_data)
                !$acc update host(var_data(1)%data_2d)
            endif
        else
            if (grid%is3d) then
                if (dqdt) then
                    !$acc update device(var%dqdt_3d)
                endif
                !$acc update device(var%data_3d)
                call halo%exch_var(var, do_dqdt=dqdt, corners=corners)
                !$acc update host(var%data_3d)
                if (dqdt) then
                    !$acc update host(var%dqdt_3d)
                endif
            else
                !$acc update device(exch_var%data_2d)
                call halo%exch_var(exch_var, corners=corners)
                !$acc update host(exch_var%data_2d)
            endif
        endif
        !$acc end data

        if (batch .and. grid%is3d) var%data_3d = var_data(1)%data_3d
        if (batch .and. grid%is2d) exch_var%data_2d = var_data(1)%data_2d
        
        ! now loop through all memory indices and check that the value in var%data_3d
        ! is equal to the expected value. If the value is not equal to the expected value
        ! check where we are (if we are in the middle, edge, or corner) and set the appropriate
        ! flag to false.

        interior = .True.
        north = .True.
        south = .True.
        east = .True.
        west = .True.
        northeast = .True.
        northwest = .True.
        southeast = .True.
        southwest = .True.

        do i = grid%ims, grid%ime
            do j = grid%jms, grid%jme
                do k = 1, max(1,grid%kte)
                    if (grid%is3d) val = var%data_3d(i,k,j)
                    if (grid%is2d) val = exch_var%data_2d(i,j)
                    if (dqdt) val = var%dqdt_3d(i,k,j)
                    if (val /= i+(j-1)*grid%nx_global+(k-1)*grid%nx_global*grid%ny_global) then
                        if (i < grid%its) then
                            if (j < grid%jts) then
                                southwest = .False.
                            elseif (j > grid%jte) then
                                northwest = .False.
                            else
                                west = .False.
                            endif
                        elseif (i > grid%ite) then
                            if (j < grid%jts) then
                                southeast = .False.
                            elseif (j > grid%jte) then
                                northeast = .False.
                            else
                                east = .False.
                            endif
                        elseif (j < grid%jts) then
                            south = .False.
                        elseif (j > grid%jte) then
                            north = .False.
                        else
                            interior = .False.
                        endif
                    endif
                end do
            end do
        end do

        call halo%finalize()

        ! Verify exchange worked
        ! All cells should match the global index formula:
        ! - Tile cells: set during initialization, should be unchanged
        ! - Exchanged halos: filled by neighbor with matching global indices
        ! - Boundary halos: set to global index before exchange, should be unchanged

        ! check that all of the interior values remained unchanged
        if (.not.(interior)) then
            call test_failed(error, "Halo exch failed", "variable interior overwritten, "//trim(test_str))
            if (grid%is3d) write(*,*) "maxval of interior data is: ", maxval(var%data_3d(grid%its:grid%ite,1,grid%jts:grid%jte))
            if (grid%is3d) write(*,*) "minval of interior data is: ", minval(var%data_3d(grid%its:grid%ite,1,grid%jts:grid%jte))
            if (grid%is2d) write(*,*) "maxval of interior data is: ", maxval(exch_var%data_2d(grid%its:grid%ite,grid%jts:grid%jte))
            if (grid%is2d) write(*,*) "minval of interior data is: ", minval(exch_var%data_2d(grid%its:grid%ite,grid%jts:grid%jte))
            return
        endif

        ! check eastern halo (exchanged by neighbor, or boundary-preserved)
        if (.not.(east)) then
            call test_failed(error, "Halo exch failed", "Failed for eastern halo, "//trim(test_str))
            if (grid%is3d) write(*,*) "east data is: ", var%data_3d(grid%ite+1:grid%ime,1,grid%jts:grid%jte)
            if (grid%is2d) write(*,*) "east data is: ", exch_var%data_2d(grid%ite+1:grid%ime,grid%jts:grid%jte)
            return
        endif

        ! check western halo
        if (.not.(west)) then
            call test_failed(error, "Halo exch failed", "Failed for western halo, "//trim(test_str))
            if (grid%is3d) write(*,*) "west data is: ", var%data_3d(grid%ims:grid%its-1,1,grid%jts:grid%jte)
            if (grid%is2d) write(*,*) "west data is: ", exch_var%data_2d(grid%ims:grid%its-1,grid%jts:grid%jte)
            return
        endif

        ! check southern halo
        if (.not.(south)) then
            call test_failed(error, "Halo exch failed", "Failed for southern halo, "//trim(test_str))
            if (grid%is3d) write(*,*) "south data is: ", var%data_3d(grid%its:grid%ite,1,grid%jms:grid%jts-1)
            if (grid%is2d) write(*,*) "south data is: ", exch_var%data_2d(grid%its:grid%ite,grid%jms:grid%jts-1)
            return
        endif

        ! check northern halo
        if (.not.(north)) then
            call test_failed(error, "Halo exch failed", "Failed for northern halo, "//trim(test_str))
            if (grid%is3d) write(*,*) "north data is: ", var%data_3d(grid%its:grid%ite,1,grid%jte+1:grid%jme)
            if (grid%is2d) write(*,*) "north data is: ", exch_var%data_2d(grid%its:grid%ite,grid%jte+1:grid%jme)
            return
        endif

        ! Check corner halos only when corners are exchanged (corners=True or batch).
        ! When corners are not exchanged, corner cells on interior processes are
        ! undefined and cannot be verified.
        if (.not.(corners .or. batch) .or. (grid%yimages == 1 .and. grid%ximages == 1)) return
        ! corner exchange not yet implemented for 2d batch exchange, so skip corner checks in that case
        if (batch .and. grid%is2d) return
        
        if (.not.(northeast)) then
            call test_failed(error, "Halo exch failed", "Failed for northeast corner halo, "//trim(test_str))
            if (grid%is3d) write(*,*) "northeast data is: ", var%data_3d(grid%ite+1:grid%ime,1,grid%jte+1:grid%jme)
            return
        endif

        if (.not.(northwest)) then
            call test_failed(error, "Halo exch failed", "Failed for northwest corner halo, "//trim(test_str))
            if (grid%is3d) write(*,*) "northwest data is: ", var%data_3d(grid%ims:grid%its-1,1,grid%jte+1:grid%jme)
            return
        endif

        if (.not.(southeast)) then
            call test_failed(error, "Halo exch failed", "Failed for southeast corner halo, "//trim(test_str))
            if (grid%is3d) write(*,*) "southeast data is: ", var%data_3d(grid%ite+1:grid%ime,1,grid%jms:grid%jts-1)
            return
        endif

        if (.not.(southwest)) then
            call test_failed(error, "Halo exch failed", "Failed for southwest corner halo, "//trim(test_str))
            if (grid%is3d) write(*,*) "southwest data is: ", var%data_3d(grid%ims:grid%its-1,1,grid%jms:grid%jts-1)
            return
        endif
    end subroutine halo_exch_standard

end module test_halo_exch
