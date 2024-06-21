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
    use options_interface,  only : options_t
    use domain_interface,   only : domain_t
    use boundary_interface, only : boundary_t
    use output_interface,   only : output_t
    use time_step,          only : step                ! Advance the model forward in time
    use initialization,     only : split_processes, welcome_message, init_options, init_model, init_physics, init_model_state
    use timer_interface,    only : timer_t
    use time_object,        only : Time_type
    use time_delta_object,  only : time_delta_t
    use icar_constants
    use wind_iterative,     only : finalize_iter_winds
    use ioserver_interface, only : ioserver_t
    use ioclient_interface, only : ioclient_t
    use io_routines,        only : io_write
    use namelist_utils,     only : get_nml_var_default
    use land_surface,               only : lsm_init

    implicit none

    type(options_t) :: options
    type(domain_t)  :: domain
    type(boundary_t):: boundary
    !type(event_type)  :: written_ev[*], write_ev[*], read_ev[*], child_read_ev[*], end_ev[*]
    type(output_t)  :: restart_dataset
    type(output_t)  :: output_dataset
    type(ioserver_t)  :: ioserver
    type(ioclient_t)  :: ioclient

    type(timer_t)   :: initialization_timer, total_timer, input_timer, output_timer, physics_timer, wind_timer, mp_timer, adv_timer, rad_timer, lsm_timer, pbl_timer, exch_timer, send_timer, ret_timer, wait_timer, forcing_timer, diagnostic_timer, wind_bal_timer
    type(Time_type) :: next_output, next_input
    type(time_delta_t) :: small_time_delta
    
    integer :: i, ierr, exec_team
    real :: t_val
    logical :: init_flag, io_loop, info_only, gen_nml, only_namelist_check
    character(len=MAXFILELENGTH) :: namelist_file


    !Initialize MPI if needed
    init_flag = .False.
    call MPI_initialized(init_flag, ierr)
    if (.not.(init_flag)) then
        call MPI_INIT(ierr)
        init_flag = .True.
    endif

    call small_time_delta%set(1)

    call MPI_Comm_Rank(MPI_COMM_WORLD,PE_RANK_GLOBAL)
    STD_OUT_PE = (PE_RANK_GLOBAL==0)

    !-----------------------------------------
    !  Model Initialization
    
    ! Read command line options to determine what kind of run this is
    call read_co(namelist_file, info_only, gen_nml, only_namelist_check)

    if (STD_OUT_PE .and. .not.(gen_nml .or. only_namelist_check .or. info_only)) then
        call welcome_message()
        flush(output_unit)
        if (STD_OUT_PE) write(*,*) "Initializing Options"
    endif
    ! Reads user supplied model options
    call init_options(options, namelist_file, info_only=info_only, gen_nml=gen_nml, only_namelist_check=only_namelist_check)
    if (STD_OUT_PE) flush(output_unit)

    !Determine split of processes which will become I/O servers and which will be compute tasks
    !Also sets constants for the program to keep track of this splitting
    call split_processes(exec_team, domain, ioserver, ioclient)

    select case(exec_team)
    case(kCOMPUTE_TEAM)
        call total_timer%start()
        call initialization_timer%start()    

        call init_model(options,domain,boundary)
        call ioclient%init(domain, boundary, options)

        if (STD_OUT_PE) write(*,*) "Receiving initial data"
        if (STD_OUT_PE) flush(output_unit)
        call ioclient%receive(boundary)

        if (STD_OUT_PE) write(*,*) "Populating boundary object"
        if (STD_OUT_PE) flush(output_unit)
        call boundary%update_computed_vars(options, update=options%forcing%time_varying_z)

        if (STD_OUT_PE) write(*,*) "Inizializing model state"
        if (STD_OUT_PE) flush(output_unit)
        call init_model_state(options, domain, boundary) ! added boundary structure

        if (options%restart%restart) then
            if (STD_OUT_PE) write(*,*) "Reading restart data"
            call ioclient%receive_rst(domain, options)
        endif

        ! physics drivers need to be initialized after restart data are potentially read in.
        call init_physics(options, domain, boundary)
        if (STD_OUT_PE) flush(output_unit)

        call output_timer%start()
        call ioclient%push(domain)
        call output_timer%stop()

        call initialization_timer%stop() 
        !  End Compute-Side Initialization
        !--------------------------------------------------------

        next_output = options%general%start_time + options%output%output_dt
        next_input = options%general%start_time

        if (STD_OUT_PE) write(*,*) "Initialization complete, beginning physics integration."
        do while (domain%model_time + small_time_delta < options%general%end_time)
            ! -----------------------------------------------------
            !
            !  Read input data if necessary
            !
            ! -----------------------------------------------------
            if ((domain%model_time + small_time_delta) >= next_input) then
                if (STD_OUT_PE) write(*,*) ""
                if (STD_OUT_PE) write(*,*) " ----------------------------------------------------------------------"
                if (STD_OUT_PE) write(*,*) "Updating Boundary conditions"
                next_input = next_input + options%forcing%input_dt

                call input_timer%start()

                call ioclient%receive(boundary)
                
                ! after reading all variables that can be read, not compute any remaining variables (e.g. z from p+ps)
                call boundary%update_computed_vars(options, update=.True.)

                call domain%interpolate_forcing(boundary, update=.True.)
                
                ! Make the boundary condition dXdt values into units of [X]/s
                call boundary%update_delta_fields(next_input - domain%model_time)
                call domain%update_delta_fields(next_input - domain%model_time)
                
                call input_timer%stop()
            endif

            ! -----------------------------------------------------
            !
            !  Integrate physics forward in time
            !
            ! -----------------------------------------------------
            if (STD_OUT_PE) write(*,*) "Running Physics"
            if (STD_OUT_PE) write(*,*) "  Model time = ", trim(domain%model_time%as_string())
            if (STD_OUT_PE) write(*,*) "   End  time = ", trim(options%general%end_time%as_string())
            if (STD_OUT_PE) write(*,*) "  Next Input = ", trim(next_input%as_string())
            if (STD_OUT_PE) write(*,*) "  Next Output= ", trim(next_output%as_string())
            if (STD_OUT_PE) flush(output_unit)
            
            ! this is the meat of the model physics, run all the physics for the current time step looping over internal timesteps
            if (.not.(options%wind%wind_only)) then
                call physics_timer%start()
                call step(domain, boundary, step_end(next_input, next_output, options%general%end_time), options,           &
                                mp_timer, adv_timer, rad_timer, lsm_timer, pbl_timer, exch_timer, &
                                send_timer, ret_timer, wait_timer, forcing_timer, diagnostic_timer, wind_bal_timer, wind_timer)
                call physics_timer%stop()
            elseif (options%wind%wind_only) then
                call domain%apply_forcing(boundary, options, real(options%output%output_dt%seconds()))
                domain%model_time = next_output
            endif

            ! -----------------------------------------------------
            !  Write output data if it is time
            ! -----------------------------------------------------
            if ((domain%model_time + small_time_delta) >= next_output) then
                if (STD_OUT_PE) write(*,*) "Writing output file"
                call output_timer%start()
                call ioclient%push(domain)
                next_output = next_output + options%output%output_dt
                call output_timer%stop()
            endif
        end do
        
        if (options%physics%windtype==kITERATIVE_WINDS .or. options%physics%windtype==kLINEAR_ITERATIVE_WINDS) call finalize_iter_winds() 
    !
    !-----------------------------------------
        call total_timer%stop()
        if (STD_OUT_PE) then
            call MPI_Comm_Size(MPI_COMM_WORLD,i)
            write(*,*) ""
            write(*,*) "Model run from : ",trim(options%general%start_time%as_string())
            write(*,*) "           to  : ",trim(options%general%end_time%as_string())
            write(*,*) "Domain : ",trim(options%domain%init_conditions_file)
            write(*,*) "Number of images:",i
            write(*,*) ""
            write(*,*) "Average timing across compute images:"
        endif
        t_val = timer_mean(total_timer, domain%compute_comms)
        if (STD_OUT_PE) write(*,*) "total          : ", t_val 
        t_val = timer_mean(initialization_timer, domain%compute_comms)
        if (STD_OUT_PE) write(*,*) "init           : ", t_val
        t_val = timer_mean(input_timer, domain%compute_comms)
        if (STD_OUT_PE) write(*,*) "input          : ", t_val
        t_val = timer_mean(output_timer, domain%compute_comms)
        if (STD_OUT_PE) write(*,*) "output         : ", t_val
        t_val = timer_mean(physics_timer, domain%compute_comms)
        if (STD_OUT_PE) write(*,*) "physics        : ", t_val
        t_val = timer_mean(mp_timer, domain%compute_comms)
        if (STD_OUT_PE) write(*,*) "microphysics   : ", t_val
        t_val = timer_mean(adv_timer, domain%compute_comms)
        if (STD_OUT_PE) write(*,*) "advection      : ", t_val
        t_val = timer_mean(rad_timer, domain%compute_comms)
        if (STD_OUT_PE) write(*,*) "radiation      : ", t_val
        t_val = timer_mean(lsm_timer, domain%compute_comms)
        if (STD_OUT_PE) write(*,*) "LSM            : ", t_val
        t_val = timer_mean(pbl_timer, domain%compute_comms)
        if (STD_OUT_PE) write(*,*) "PBL            : ", t_val
        t_val = timer_mean(forcing_timer, domain%compute_comms)
        if (STD_OUT_PE) write(*,*) "forcing        : ", t_val 
        t_val = timer_mean(wind_bal_timer, domain%compute_comms)
        if (STD_OUT_PE) write(*,*) "wind bal       : ", t_val
        t_val = timer_mean(diagnostic_timer, domain%compute_comms)
        if (STD_OUT_PE) write(*,*) "diagnostic     : ", t_val
        t_val = timer_mean(send_timer, domain%compute_comms)
        if (STD_OUT_PE) write(*,*) "halo-exchange(send)  : ", t_val 
        t_val = timer_mean(ret_timer, domain%compute_comms)
        if (STD_OUT_PE) write(*,*) "halo-exchange(retrieve)  : ", t_val
        t_val = timer_mean(wind_timer, domain%compute_comms)
        if (STD_OUT_PE) write(*,*) "winds          : ", t_val     

    case(kIO_TEAM)
    
        !This is, unfortunately, stupidly, needed to allow coarrays to work. 
        !Perhaps on later releases with better team support it can be removed?
#ifdef  CRAY_PE
        call domain%init(options)
#endif
        call ioserver%init(options)

        !Get initial conditions
        call ioserver%read_file()
        
        if (options%restart%restart) then
            call ioserver%read_restart_file(options)
        endif

        !Get first boundary conditions
        call ioserver%read_file()

        !  End IO-Side Initialization
        !--------------------------------------------------------

        next_output = options%general%start_time
        next_input = options%general%start_time

        io_loop = .True.

        do while (io_loop)
            
            ! If we aren't yet done, then wait for an output
            ! If we are at an output step, do it now
            if (ioserver%io_time+small_time_delta >= next_output .and. ioserver%io_time <= options%general%end_time) then
                call ioserver%write_file(ioserver%io_time)
                next_output = next_output + options%output%output_dt
            endif

            !See if we even have files to read
            if (ioserver%files_to_read) then
                !See of it is time to read.
                if (ioserver%io_time+small_time_delta >= next_input) then
                    call ioserver%read_file()
                    next_input = next_input + options%forcing%input_dt
                endif

                !If we still have files to read, set io_time to be the soonest of either input or output time
                ioserver%io_time = step_end(next_input,next_output)
            else
                !If we have no more files to read, the only thing we could otherwise do is write output
                ioserver%io_time = next_output
            endif

            if (.not.(ioserver%files_to_read) .and. ioserver%io_time > options%general%end_time) io_loop = .False.

        enddo
        !If we are done with the program
        call ioserver%close_files()
    end select
    
    CALL MPI_Finalize()
contains

    function timer_mean(timer,comms) result(mean_t)
        implicit none
        type(timer_t), intent(inout) :: timer
        type(MPI_Comm), intent(in) :: comms

        real :: mean_t, t_sum
        integer :: ierr
            
        t_sum = timer%get_time()
        call MPI_Allreduce(MPI_IN_PLACE,t_sum,1,MPI_REAL,MPI_SUM,comms,ierr)
        mean_t = t_sum/kNUM_COMPUTE
    
    end function

    function step_end(time1, time2, end_time) result(min_time)
        implicit none
        type(Time_type), intent(in) :: time1, time2
        type(Time_type), optional, intent(in) :: end_time
        type(Time_type) :: min_time

        if (time1 <= time2 ) then
            min_time = time1
        else
            min_time = time2
        endif

        if (present(end_time)) then
            if (end_time < min_time) then
                min_time = end_time
            endif
        endif
    end function
    
    subroutine read_co(nml_file, info, gen_nml, only_namelist_check)
        implicit none
        character(len=MAXFILELENGTH), intent(out) :: nml_file
        logical, intent(out) :: info, gen_nml, only_namelist_check

        integer :: cnt, p
        character(len=MAXFILELENGTH) :: first_arg, arg, default
        logical :: file_exists

        nml_file = ""
        info = .False.
        gen_nml = .False.
        only_namelist_check = .False.

        cnt = command_argument_count()

        ! If there are no command line arguments, throw error
        if (cnt == 0 .and. STD_OUT_PE) then
            write(*,*) "Usage: ./HICAR [-v [variable_name ...|--all]] [--check-nml] [--gen-nml] namelist_file"
            write(*,*) "    -v [variable_name ...|--all]: Print information about the namelist variable(s) variable_name, ... "
            write(*,*) "                                  --all prints out information for all namelist variables."
            write(*,*) "    --check-nml:                  Check the namelist file for errors without running the model."
            write(*,*) "    --gen-nml:                    Generate a namelist file with default values."
            write(*,*) "    namelist_file:                The name of the namelist file to use."
            write(*,*)
            write(*,*) "    Example to check namelist:                           ./HICAR --check-nml namelist_file.nml"
            write(*,*) "    Example to run model:                                ./HICAR namelist_file.nml"
            write(*,*) "    Example to check namelist variable:                  ./HICAR -v mp"
            write(*,*) "    Example to generate namelist variable documentation: ./HICAR -v --all > namelist_doc.txt"

            write(*,*)


            stop
        endif

        ! get first command line argument
        call get_command_argument(1, first_arg)

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
        else
            nml_file = first_arg
        endif

        if (.not.first_arg=='-v' .and. .not.first_arg=='--gen-nml') then
            ! Check that the options file actually exists
            INQUIRE(file=trim(nml_file), exist=file_exists)

            ! if options file does not exist, print an error and quit
            if (.not.file_exists) then
                if (STD_OUT_PE) write(*,*) "Using options file = ", trim(nml_file)
                stop "Options file does not exist. "
            endif
        endif
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
