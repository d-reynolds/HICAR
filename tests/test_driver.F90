! @file test_driver.F90
! @brief This program manages all the tests for the HICAR code. Mostly entirely taken from the test-drive library: https://github.com/fortran-lang/test-drive/tree/main
! @author Dylan Reynolds
! @date 2025-03-05
! @version 1.0
! @details This program initializes MPI, creates test suites, and runs the tests. It also handles output redirection for non-root processes.
program test_driver
    use testdrive, only : run_testsuite, new_testsuite, testsuite_type, &
                            select_suite, run_selected, get_argument
    use test_halo_exch, only : collect_halo_exch_suite
    use test_advect, only : collect_advect_suite
    use mpi_f08
    use iso_fortran_env
    use icar_constants, only: STD_OUT_PE
    use string, only: to_lower
    implicit none
    integer :: stat, is, ierr, my_index, null_unit
    character(len=:), allocatable :: suite_name, test_name
    type(testsuite_type), allocatable :: testsuites(:)
    character(len=*), parameter :: fmt = '("#", *(1x, a))'
    logical :: init_flag
    logical :: no_test_run = .True.
    character(len=9) :: file

    stat = 0

    !Initialize MPI if needed
    init_flag = .False.
    call MPI_initialized(init_flag, ierr)
    if (.not.(init_flag)) then
        call MPI_INIT(ierr)
        init_flag = .True.
    endif


    call MPI_Comm_Rank(MPI_COMM_WORLD,my_index,ierr)
    my_index = my_index + 1

    if (my_index /= 1) then
        ! Redirect stdout and stderr to /dev/null for non-root processes
        file = '.tmp_dupplicate_test_output'
        open(newunit=null_unit, file=file, status='replace')
        ! Redirect standard output
        close(output_unit)
        open(output_unit, file=file, status='replace')
        ! Redirect error output
        close(error_unit)
        open(error_unit, file=file, status='replace')
    else
        STD_OUT_PE = .True.
    end if

  


    ! -------------------------------------------------------------------------------------------------
    ! -------------------------------------------------------------------------------------------------
    ! Create the test suites, only code in the program body which 
    ! the developer must touch when adding new test suites
    ! -------------------------------------------------------------------------------------------------
    ! -------------------------------------------------------------------------------------------------
    testsuites = [ &
        new_testsuite("halo_exch", collect_halo_exch_suite), &
        new_testsuite("advection", collect_advect_suite) &
        ]
    ! -------------------------------------------------------------------------------------------------
    ! -------------------------------------------------------------------------------------------------
    ! -------------------------------------------------------------------------------------------------







    call get_argument(1, suite_name)
    call get_argument(2, test_name)

    if (allocated(suite_name)) then
        if (.not.(to_lower(suite_name) == 'all')) then
            is = select_suite(testsuites, suite_name)
            if (is > 0 .and. is <= size(testsuites)) then
                if (allocated(test_name)) then
                    write(error_unit, fmt) "Suite:", testsuites(is)%name
                    call run_selected(testsuites(is)%collect, test_name, error_unit, stat)
                    if (stat < 0) then
                    error stop 1
                    end if
                else
                    write(error_unit, fmt) "Testing:", testsuites(is)%name
                    call run_testsuite(testsuites(is)%collect, error_unit, stat)
                end if
            else
                write(error_unit, fmt) "Available testsuites"
                do is = 1, size(testsuites)
                    write(error_unit, fmt) "-", testsuites(is)%name
                end do
                error stop 1
            end if

            no_test_run = .False.
        endif
    endif
    if (no_test_run) then
        do is = 1, size(testsuites)
        write(error_unit, fmt) "Testing:", testsuites(is)%name
        call run_testsuite(testsuites(is)%collect, error_unit, stat)
        end do
    end if

    if (stat > 0) then
      write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
    end if

    call MPI_Finalize(ierr)
  
  end program test_driver