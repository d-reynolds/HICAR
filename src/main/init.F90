!> ----------------------------------------------------------------------------
!!  Model Initialization includes allocating memory for boundary and domain
!!      data structures.  It reads all of the options from the namelist
!!      file (or files).  It also reads in Lat/Lon and Terrain data.  This module
!!      also sets up geographic (and vertical) look uptables for the forcing data
!!      Finally, there is a driver routine to initialize all model physics packages
!!
!!   The module has been updated to allow arbitrary named variables
!!       this allows the use of e.g. ERAi, but still is not as flexible as it could be
!!
!!   The use of various python wrapper scripts in helpers/ makes it easy to add new
!!       datasets, and make them conform to the expectations of the current system.
!!      For now there are no plans to near term plans to substantially modify this.
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!! ----------------------------------------------------------------------------
module initialization
    !use data_structures
    use options_interface,  only : options_t
    use domain_interface,   only : domain_t
    use boundary_interface, only : boundary_t
    use microphysics,               only : mp_init
    use advection,                  only : adv_init
    use radiation,                  only : radiation_init
    use convection,                 only : init_convection
    use planetary_boundary_layer,   only : pbl_init
    use land_surface,               only : lsm_init
    use surface_layer,              only : sfc_init
    use io_routines,                only : io_read, io_write
    use mod_atm_utilities,          only : init_atm_utilities
    use wind,                       only : update_winds, init_winds

    use icar_constants!,             only : kITERATIVE_WINDS, kWIND_LINEAR
    use ioserver_interface,         only : ioserver_t
    use ioclient_interface,         only : ioclient_t
    use iso_fortran_env
    use mpi_f08


    ! use io_routines,                only : io_read, &
    !                                        io_write3d,io_write3di, io_write
    ! use geo,                        only : geo_LUT, geo_interp, geo_interp2d, standardize_coordinates
    ! use vertical_interpolation,     only : vLUT, vinterp
    ! use wind,                       only : init_winds
    ! use initialize_options,         only : init_options
    ! use string,                     only : str


    implicit none
    private
    public::split_processes, welcome_message, init_options, init_model, init_physics, init_model_state

contains


    subroutine split_processes(exec_team, domain, ioserver, ioclient)
        implicit none
        integer, intent(inout) :: exec_team
        type(domain_t), intent(inout) :: domain
        type(ioserver_t), intent(inout) :: ioserver
        type(ioclient_t), intent(inout) :: ioclient

        integer :: n, k, name_len, color, ierr, node_name_i, num_PE
        integer :: num_threads

        character(len=MPI_MAX_PROCESSOR_NAME) :: node_name
        integer, allocatable :: node_names(:) 
        type(MPI_Comm) :: globalComm, splitComm

#if defined(_OPENMP)
        num_threads = omp_get_max_threads()
#else
        num_threads = 1
#endif    

        call MPI_Comm_Size(MPI_COMM_WORLD,num_PE)

        if (STD_OUT_PE) then
            write(*,*) "  Number of processing elements:",num_PE
            write(*,*) "  Max number of OpenMP Threads:",num_threads
        endif
    
        allocate(node_names(num_PE))

        node_names = 0
        node_name_i = 0

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

    subroutine init_options(options, namelist_file, info_only, gen_nml, only_namelist_check)
        implicit none
        type(options_t), intent(inout) :: options
        character(len=*), intent(in) :: namelist_file
        logical, optional, intent(in) :: info_only, gen_nml, only_namelist_check

        ! read in options file
        if (STD_OUT_PE) write(*,*) "Initializing Options"
        if (STD_OUT_PE) flush(output_unit)
        call options%init(namelist_file, info_only=info_only, gen_nml=gen_nml, only_namelist_check=only_namelist_check)

    end subroutine init_options

    subroutine init_model(options,domain,boundary)
        implicit none
        type(options_t), intent(inout) :: options
        type(domain_t),  intent(inout) :: domain
        type(boundary_t),intent(inout) :: boundary ! forcing file for init conditions
        
        if (STD_OUT_PE) write(*,*) "Initializing Domain"
        if (STD_OUT_PE) flush(output_unit)
        call domain%init(options)

        if (STD_OUT_PE) write(*,*) "Initializing boundary condition data structure"
        if (STD_OUT_PE) flush(output_unit)
        call boundary%init(options,domain%latitude%data_2d,domain%longitude%data_2d,domain%variables_to_force)

        if (STD_OUT_PE) write(*,*) "Initializing atmospheric utilities"
        if (STD_OUT_PE) flush(output_unit)
        ! initialize the atmospheric helper utilities
        call init_atm_utilities(options)

        if (STD_OUT_PE) write(*,'(/ A)') "Finished basic initialization"
        if (STD_OUT_PE) write(*,'(A /)') "---------------------------------------"

    end subroutine init_model

    subroutine init_physics(options, domain, forcing)
        implicit none
        type(options_t), intent(inout) :: options
        type(domain_t),  intent(inout) :: domain
        type(boundary_t),intent(in)    :: forcing


        if (STD_OUT_PE) write(*,*) "Init initial winds"
        if (STD_OUT_PE) flush(output_unit)
        call init_winds(domain,options)

        if (STD_OUT_PE) write(*,*) "Updating initial winds"
        if (STD_OUT_PE) flush(output_unit)
        call update_winds(domain, forcing, options)

        ! initialize microphysics code (e.g. compute look up tables in Thompson et al)
        call mp_init(options) !this could easily be moved to init_model...
        if (STD_OUT_PE) flush(output_unit)
        call init_convection(domain,options)
        if (STD_OUT_PE) flush(output_unit)

        call pbl_init(domain,options)
        if (STD_OUT_PE) flush(output_unit)

        call radiation_init(domain,options)
        if (STD_OUT_PE) flush(output_unit)

        call lsm_init(domain,options)
        if (STD_OUT_PE) flush(output_unit)

        call sfc_init(domain,options)
        if (STD_OUT_PE) flush(output_unit)

        call adv_init(domain,options)
        if (STD_OUT_PE) flush(output_unit)

    end subroutine init_physics

    subroutine init_model_state(options,domain,boundary)
        implicit none
        type(options_t), intent(inout) :: options
        type(domain_t),  intent(inout) :: domain
        type(boundary_t),intent(inout) :: boundary ! forcing file for init conditions

        if (STD_OUT_PE) write(*,*) "Reading Initial conditions from boundary dataset"
        call domain%get_initial_conditions(boundary, options)

    end subroutine init_model_state


    subroutine welcome_message()
        implicit none

        write(*,*) ""
        write(*,*) "============================================================"
        write(*,*) "|                                                          |"
        write(*,*) "|  The Intermediate Complexity Atmospheric Research Model  |"
        write(*,*) "|                          (ICAR)                          |"
        write(*,*) "|                                                          |"
        write(*,*) "|   Developed at NCAR:                                     |"
        write(*,*) "|     The National Center for Atmospheric Research         |"
        write(*,*) "|     NCAR is sponsored by the National Science Foundation |"
        write(*,*) "|                                                          |"
        write(*,*) "|   Version: ",kVERSION_STRING,"                                         |"
        write(*,*) "|                                                          |"
        write(*,*) "============================================================"
        write(*,*) ""

    end subroutine welcome_message

end module
