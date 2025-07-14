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
#ifdef _OPENACC
    use openacc
#endif
    use options_interface,  only : options_t
    use flow_object_interface, only : comp_arr_t
    use boundary_interface, only : boundary_t
    use ioclient_interface, only : ioclient_t

    use initialization,     only : split_processes, init_options
    use icar_constants
    use namelist_utils,     only : get_nml_var_default
    use output_metadata,    only : list_output_vars
    use flow_events,        only : component_init, component_loop, component_program_end
    implicit none

    type(options_t), allocatable :: options(:)
    type(boundary_t), allocatable :: boundary(:)
    type(comp_arr_t)  :: components(kMAX_NESTS)
    type(ioclient_t), allocatable  :: ioclient(:)
    
    integer :: i, ierr, exec_team, n_nests, n, old_nest, PE_RANK_GLOBAL
    real :: t_val, t_val2, t_val3
    logical :: init_flag, info_only, gen_nml, only_namelist_check
    logical :: start_time_match = .False.
    character(len=kMAX_FILE_LENGTH) :: namelist_file
#ifdef _OPENACC
    integer :: dev, devNum, local_rank, comm_size
    type(MPI_Comm) :: local_comm
    integer(acc_device_kind) :: devtype
#endif

    ! Read command line options to determine what kind of run this is
    call read_co(namelist_file, info_only, gen_nml, only_namelist_check)

    !Initialize MPI if needed
    init_flag = .False.
    call MPI_initialized(init_flag)
    if (.not.(init_flag)) then
        call MPI_INIT()
        init_flag = .True.
    endif

#ifdef _OPENACC
!
! ****** Set the Accelerator device number based on local rank
!
     call MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, &
          MPI_INFO_NULL, local_comm)
     call MPI_Comm_rank(local_comm, local_rank)
     call MPI_Comm_size(local_comm, comm_size)

     devtype = acc_get_device_type()
     devNum = acc_get_num_devices(devtype)
     dev = mod(local_rank,devNum)
     if ((local_rank + 1) < comm_size) then
        call acc_set_device_num(dev, devtype)
        call acc_init(devtype)
     endif
# endif

    call MPI_Comm_Rank(MPI_COMM_WORLD,PE_RANK_GLOBAL)
    STD_OUT_PE = (PE_RANK_GLOBAL==0)

    !-----------------------------------------
    !  Model Initialization

    ! Reads user supplied model options
    call init_options(options, namelist_file, info_only=info_only, gen_nml=gen_nml, only_namelist_check=only_namelist_check)

    if (STD_OUT_PE) write(*,'(/ A)') "--------------------------------------------------------"
    if (STD_OUT_PE) write(*,'(A)')   "Finished reading options, beginning processor assignment"
    if (STD_OUT_PE) write(*,'(A /)') "--------------------------------------------------------"
    if (STD_OUT_PE) flush(output_unit)

    n_nests = options(1)%general%nests

    ! !Allocate the multiple domains, boundarys
    allocate(boundary(n_nests))
    allocate(ioclient(n_nests))

    !Determine split of processes which will become I/O servers and which will be compute tasks
    !Also sets constants for the program to keep track of this splitting
    call split_processes(components, ioclient, n_nests, options(1))
    if (STD_OUT_PE) write(*,'(/ A)') "--------------------------------------------------------------"
    if (STD_OUT_PE) write(*,'(A)')   "Finished processor assignment, beginning domain initialization"
    if (STD_OUT_PE) write(*,'(A)')   "--------------------------------------------------------------"
    if (STD_OUT_PE) flush(output_unit)

    do i = 1, n_nests
        call component_init(components(i)%comp, options, boundary(i), ioclient(i), i)
    enddo

    !  End Compute-Side Initialization
    !--------------------------------------------------------
    if (STD_OUT_PE) write(*,'(/ A)') "------------------------------------------------------"
    if (STD_OUT_PE) write(*,'(A)')   "Initialization complete, beginning physics integration"
    if (STD_OUT_PE) write(*,'(A)')   "------------------------------------------------------"

    call component_loop(components(1:n_nests), options, boundary, ioclient)

    call component_program_end(components(1:n_nests), options)

    CALL MPI_Finalize()

contains
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

        !Turn this on to enable output by default -- likely that we are only running as a single process
        !If we do any output from this subroutine
        STD_OUT_PE = .True.

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
        STD_OUT_PE = .False.
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
