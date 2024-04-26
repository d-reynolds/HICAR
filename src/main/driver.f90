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
    use, intrinsic iso_fortran_env
    use mpi_f08
    use options_interface,  only : options_t
    use domain_interface,   only : domain_t
    use boundary_interface, only : boundary_t
    use output_interface,   only : output_t
    use time_step,          only : step                ! Advance the model forward in time
    use initialization,     only : init_model, init_physics, init_model_state
    use timer_interface,    only : timer_t
    use time_object,        only : Time_type
    use time_delta_object,  only : time_delta_t
    use icar_constants
    use wind_iterative,     only : init_petsc_comms, finalize_iter_winds
    use ioserver_interface, only : ioserver_t
    use ioclient_interface, only : ioclient_t
    use io_routines,        only : io_write

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
    logical :: init_flag, io_loop

    !if (STD_OUT_PE) write(*,*) 'Before MPI init'
    !if (STD_OUT_PE) flush(output_unit)

    !Initialize MPI if needed
    init_flag = .False.
    call MPI_initialized(init_flag, ierr)
    if (.not.(init_flag)) then
        call MPI_INIT(ierr)
        init_flag = .True.
    endif

    !Determine split of processes which will become I/O servers and which will be compute tasks
    !Also sets constants for the program to keep track of this splitting

    !We should first read in options here

    call split_processes(exec_team, domain, ioserver, ioclient)

    call small_time_delta%set(1)
        
    !-----------------------------------------
    !  Model Initialization
    !
    ! Reads config options and initializes domain and boundary conditions
    ! First the compute domain has to be established
    
    ! Model structure must be initialized first. 
    ! The domain dimensions will be used when establishing I/O client-server relations

    if (STD_OUT_PE) write(*,*) 'Before init'
    if (STD_OUT_PE) flush(output_unit)

    call init_model(options, domain, boundary)
    if (STD_OUT_PE) flush(output_unit)
    !call init_IO(exec_team, domain, boundary, options, ioclient, ioserver)        
    !if (STD_OUT_PE) flush(output_unit)

    select case(exec_team)
    case(kCOMPUTE_TEAM)
        call total_timer%start()
        call initialization_timer%start()    

        call ioclient%init(domain, boundary, options)

        if (STD_OUT_PE) write(*,*) "Receiving initial data"
        if (STD_OUT_PE) flush(output_unit)

        call ioclient%receive(boundary)

        if (STD_OUT_PE) write(*,*) "Populating boundary object"
        if (STD_OUT_PE) flush(output_unit)

        call boundary%update_computed_vars(options, update=options%parameters%time_varying_z)

        if (STD_OUT_PE) write(*,*) "Inizializing model state"
        if (STD_OUT_PE) flush(output_unit)

        call init_model_state(options, domain, boundary) ! added boundary structure

        if (options%parameters%restart) then
            if (STD_OUT_PE) write(*,*) "Reading restart data"
            call ioclient%receive_rst(domain)
        endif
        if (STD_OUT_PE) flush(output_unit)
        ! physics drivers need to be initialized after restart data are potentially read in.
        call init_physics(options, domain, boundary)
        if (STD_OUT_PE) flush(output_unit)

        call output_timer%start()
        call ioclient%push(domain)
        call output_timer%stop()

        call initialization_timer%stop() 

        next_output = options%parameters%start_time + options%io_options%output_dt
        next_input = options%parameters%start_time + options%io_options%input_dt


        if (STD_OUT_PE) write(*,*) "Initialization complete, beginning physics integration."
        do while (domain%model_time + small_time_delta < options%parameters%end_time)

            ! -----------------------------------------------------
            !
            !  Read input data if necessary
            !
            ! -----------------------------------------------------
            if ((domain%model_time + small_time_delta + options%io_options%input_dt) >= next_input) then
                if (STD_OUT_PE) write(*,*) ""
                if (STD_OUT_PE) write(*,*) " ----------------------------------------------------------------------"
                if (STD_OUT_PE) write(*,*) "Updating Boundary conditions"
                
                call input_timer%start()

                call ioclient%receive(boundary)
                
                ! after reading all variables that can be read, not compute any remaining variables (e.g. z from p+ps)
                call boundary%update_computed_vars(options, update=.True.)

                ! if the vertical levels of the forcing data change over time, they need to be interpolated to the original levels here.
                if (options%parameters%time_varying_z) then
                    call boundary%interpolate_original_levels(options)
                endif
                call domain%interpolate_forcing(boundary, update=.True.)
                
                ! Make the boundary condition dXdt values into units of [X]/s
                call boundary%update_delta_fields(next_input - domain%model_time)
                call domain%update_delta_fields(next_input - domain%model_time)
                
                next_input = next_input + options%io_options%input_dt
                call input_timer%stop()
            endif


            ! -----------------------------------------------------
            !
            !  Integrate physics forward in time
            !
            ! -----------------------------------------------------
            if (STD_OUT_PE) write(*,*) "Running Physics"
            if (STD_OUT_PE) write(*,*) "  Model time = ", trim(domain%model_time%as_string())
            if (STD_OUT_PE) write(*,*) "   End  time = ", trim(options%parameters%end_time%as_string())
            if (STD_OUT_PE) write(*,*) "  Next Input = ", trim(next_input%as_string())
            if (STD_OUT_PE) write(*,*) "  Next Output= ", trim(next_output%as_string())
            if (STD_OUT_PE) flush(output_unit)
            
            ! this is the meat of the model physics, run all the physics for the current time step looping over internal timesteps
            if (.not.(options%wind%wind_only)) then
                call physics_timer%start()
                call step(domain, boundary, step_end(next_input, next_output, options%parameters%end_time), options,           &
                                mp_timer, adv_timer, rad_timer, lsm_timer, pbl_timer, exch_timer, &
                                send_timer, ret_timer, wait_timer, forcing_timer, diagnostic_timer, wind_bal_timer, wind_timer)
                call physics_timer%stop()
            elseif (options%wind%wind_only) then
                call domain%apply_forcing(boundary, options, real(options%io_options%output_dt%seconds()))
                domain%model_time = next_output
            endif

            ! -----------------------------------------------------
            !  Write output data if it is time
            ! -----------------------------------------------------
            if ((domain%model_time + small_time_delta) >= next_output) then
                if (STD_OUT_PE) write(*,*) "Writing output file"
                call output_timer%start()
                call ioclient%push(domain)
                next_output = next_output + options%io_options%output_dt
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
            write(*,*) "Model run from : ",trim(options%parameters%start_time%as_string())
            write(*,*) "           to  : ",trim(options%parameters%end_time%as_string())
            write(*,*) "Domain : ",trim(options%parameters%init_conditions_file)
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
    
        call ioserver%init(domain, options)

        call ioserver%read_file()
        
        if (options%parameters%restart) then
            !sync all !Make sure we don't overwrite the file previously read in before the compute task has used it
            call ioserver%read_restart_file(options)
        endif

        !We do not do output here so that we can continue straight to reading the next input file and not block the compute processes
        next_output = options%parameters%start_time
        next_input = options%parameters%start_time !+ options%io_options%input_dt

        io_loop = .True.

        do while (io_loop)
            !See if we even have files to read
            if (ioserver%files_to_read) then
                !See of it is time to read
                if (ioserver%io_time+small_time_delta >= next_input) then
                    call ioserver%read_file()
                    next_input = next_input + options%io_options%input_dt
                endif
            endif
            
            ! If we aren't yet done, then wait for an output
            if (ioserver%io_time <= options%parameters%end_time) then
                call ioserver%write_file(ioserver%io_time)
                ioserver%io_time = ioserver%io_time + options%io_options%output_dt
            endif
            
            if (.not.(ioserver%files_to_read) .and. ioserver%io_time > options%parameters%end_time) io_loop = .False.
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
        type(Time_type), intent(in) :: time1, time2, end_time
        type(Time_type) :: min_time

        if (time1 <= time2 ) then
            min_time = time1
        else
            min_time = time2
        endif

        if (end_time < min_time) then
            min_time = end_time
        endif
    end function


    subroutine split_processes(exec_team, domain, ioserver, ioclient)
        implicit none
        integer, intent(inout) :: exec_team
        type(domain_t), intent(inout) :: domain
        type(ioserver_t), intent(inout) :: ioserver
        type(ioclient_t), intent(inout) :: ioclient

        integer :: n, k, name_len, color, ierr, node_name_i, num_PE
        character(len=MPI_MAX_PROCESSOR_NAME) :: node_name
        integer, allocatable :: node_names(:) 

        type(MPI_Comm) :: globalComm, splitComm
 
        call MPI_Comm_Size(MPI_COMM_WORLD,num_PE)

        allocate(node_names(num_PE))

        node_names = 0
        node_name_i = 0
       
        call MPI_Comm_Rank(MPI_COMM_WORLD,PE_RANK_GLOBAL)
        STD_OUT_PE = (PE_RANK_GLOBAL==0)

        !First, determine information about node/CPU configuration
        call MPI_Get_processor_name(node_name, name_len, ierr)
        do n = 1,name_len
            node_name_i = node_name_i + ichar(node_name(n:n))*n*10
        enddo

        node_names(PE_RANK_GLOBAL+1) = node_name_i

        !Get list of node names on all processes        
        call MPI_Allreduce(MPI_IN_PLACE,node_names,num_PE,MPI_INT,MPI_MAX,MPI_COMM_WORLD,ierr)
        
        !Set global constants related to resource distribution
        kNUM_PROC_PER_NODE = count(node_names==node_names(1))

        !Assign one io process per node, this results in best co-array transfer times
        kNUM_SERVERS = ceiling(num_PE*1.0/kNUM_PROC_PER_NODE)
        kNUM_COMPUTE = num_PE-kNUM_SERVERS
        
        if ((mod(kNUM_COMPUTE,2) /= 0) .and. STD_OUT_PE) then
            write(*,*) 'WARNING: number of compute processes is odd-numbered.' 
            write(*,*) 'One process per node is used for I/O.'
            write(*,*) 'If the total number of compute processes is odd-numbered,'
            write(*,*) 'this may lead to errors with domain decomposition'
        endif

        k = 1
        allocate(DOM_IMG_INDX(kNUM_COMPUTE))
        do n = 1,kNUM_COMPUTE
            if (mod(k,(num_PE/kNUM_SERVERS)) == 0) k = k+1
            DOM_IMG_INDX(n) = k
            k = k+1
        enddo


        !-----------------------------------------
        ! Assign Compute and IO processes
        !-----------------------------------------
        if (mod((PE_RANK_GLOBAL+1),(num_PE/kNUM_SERVERS)) == 0) then
            exec_team = kIO_TEAM
            color = 1
        else
            exec_team = kCOMPUTE_TEAM
            color = 0
        endif
       
        CALL MPI_Comm_dup( MPI_COMM_WORLD, globalComm, ierr )
        ! Assign all images in the IO team to the IO_comms MPI communicator. Use image indexing within initial team to get indexing of global MPI ranks
        CALL MPI_Comm_split( globalComm, color, PE_RANK_GLOBAL, splitComm, ierr )

        select case (exec_team)
        case (kCOMPUTE_TEAM)
            CALL MPI_Comm_dup( splitComm, domain%compute_comms, ierr )
            ioserver%IO_comms = MPI_COMM_NULL
        case (kIO_TEAM)
            CALL MPI_Comm_dup( splitComm, ioserver%IO_comms, ierr )
            domain%compute_comms = MPI_COMM_NULL
        end select

        !-----------------------------------------
        ! Group server and client processes
        !-----------------------------------------
        CALL MPI_Comm_dup( MPI_COMM_WORLD, globalComm, ierr )
        ! Group IO clients with their related server process. This is basically just grouping processes by node
        CALL MPI_Comm_split( globalComm, node_names(PE_RANK_GLOBAL+1), PE_RANK_GLOBAL, splitComm, ierr )

        select case (exec_team)
        case (kCOMPUTE_TEAM)
            CALL MPI_Comm_dup( splitComm, ioclient%parent_comms, ierr )
            ioserver%client_comms = MPI_COMM_NULL
        case (kIO_TEAM)
            CALL MPI_Comm_dup( splitComm, ioserver%client_comms, ierr )
            ioclient%parent_comms = MPI_COMM_NULL
        end select

    end subroutine split_processes
    
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
