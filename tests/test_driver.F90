program test_driver
    use, intrinsic :: iso_fortran_env, only : error_unit
    use testdrive, only : run_testsuite, new_testsuite, testsuite_type
    use test_halo_exch, only : collect_halo_exch_suite
    use mpi_f08
    use iso_fortran_env
    implicit none
    integer :: stat, is, ierr, my_index, null_unit
    type(testsuite_type), allocatable :: testsuites(:)
    character(len=*), parameter :: fmt = '("#", *(1x, a))'
    logical :: init_flag
    character(len=9) :: file

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
    end if

    stat = 0
  
    testsuites = [ &
      new_testsuite("halo_exch", collect_halo_exch_suite) &
      ]
  
    do is = 1, size(testsuites)
      if (my_index == 1) write(error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
    end do
  
    if (stat > 0) then
      write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
    end if

    call MPI_Finalize(ierr)
  
  end program test_driver