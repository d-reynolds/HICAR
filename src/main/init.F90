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
    use options_interface,  only : options_t, general_namelist, inter_nest_options_check
    use options_types,      only : general_options_type
    use domain_interface,   only : domain_t
    use boundary_interface, only : boundary_t
    use flow_object_interface, only : flow_obj_t
    use microphysics,               only : mp_init, mp_var_request
    use advection,                  only : adv_init, adv_var_request
    use radiation,                  only : radiation_init, ra_var_request
    use convection,                 only : init_convection, cu_var_request
    use planetary_boundary_layer,   only : pbl_init, pbl_var_request
    use land_surface,               only : lsm_init, lsm_var_request
    use surface_layer,              only : sfc_init, sfc_var_request
    use io_routines,                only : io_read, io_write, io_newunit
    use wind,                       only : update_winds, init_winds, wind_var_request
    use omp_lib,                    only : omp_get_max_threads
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


    subroutine split_processes(components, ioclient, n_nests)
        implicit none
        class(flow_obj_t), allocatable, intent(inout) :: components(:)
        type(ioclient_t), intent(inout) :: ioclient(:)
        integer, intent(in) :: n_nests

        integer :: n, k, name_len, color, ierr, node_name_i, num_PE
        integer :: num_threads, found, PE_RANK_GLOBAL, NUM_SERVERS, NUM_COMPUTE, NUM_IO_PER_NODE, NUM_PROC_PER_NODE

        character(len=MPI_MAX_PROCESSOR_NAME) :: node_name, ENV_IO_PER_NODE
        integer, allocatable :: node_names(:) 
        type(MPI_Comm) :: globalComm, splitComm, shared_comm

#if defined(_OPENMP)
        num_threads = omp_get_max_threads()
#else
        num_threads = 1
#endif    

        ! See if the environment variable HICAR_IO_PER_NODE is set
        ! If it is, then we will use that to determine the number of IO processes per node
        call get_environment_variable('HICAR_IO_PER_NODE',ENV_IO_PER_NODE,found)

        ! If the environment variable is not set, then we will use 1
        if (found == 0) then
            NUM_IO_PER_NODE = 1
        else
            ! If the number of IO processes per node is set, then we will use that
            read(ENV_IO_PER_NODE,*) NUM_IO_PER_NODE
        endif

        call MPI_Comm_Size(MPI_COMM_WORLD,num_PE)
        call MPI_Comm_Rank(MPI_COMM_WORLD,PE_RANK_GLOBAL)

        !Discover the number of processors per node
        call MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, shared_comm, ierr)
        call MPI_Comm_Size(shared_comm,NUM_PROC_PER_NODE)

        !Assign one io process per node, this results in best co-array transfer times
        NUM_SERVERS = ceiling(num_PE*NUM_IO_PER_NODE*1.0/NUM_PROC_PER_NODE)
        NUM_COMPUTE = num_PE-NUM_SERVERS
        
        if ((mod(NUM_COMPUTE,2) /= 0) .and. STD_OUT_PE) then
            write(*,*) 'WARNING: number of compute processes is odd-numbered.' 
            write(*,*) 'One process per node is used for I/O.'
            write(*,*) 'If the total number of compute processes is odd-numbered,'
            write(*,*) 'this may lead to errors with domain decomposition'
        endif

        !-----------------------------------------
        ! Assign Compute and IO processes
        !-----------------------------------------
        if (mod((PE_RANK_GLOBAL+1),(num_PE/NUM_SERVERS)) == 0) then
            color = 1
            allocate(ioserver_t::components(n_nests))
        else
            color = 0
            allocate(domain_t::components(n_nests))
        endif
    
        CALL MPI_Comm_dup( MPI_COMM_WORLD, globalComm, ierr )
        ! Assign all images in the IO team to the IO_comms MPI communicator. Use image indexing within initial team to get indexing of global MPI ranks
        CALL MPI_Comm_split( globalComm, color, PE_RANK_GLOBAL, splitComm, ierr )

        do n = 1, n_nests
            associate(current_component => components(n))
                select type (current_component)
                    type is (domain_t)
                        CALL MPI_Comm_dup( splitComm, current_component%compute_comms, ierr )
                    type is (ioserver_t)
                        CALL MPI_Comm_dup( splitComm, current_component%IO_comms, ierr )
                end select
            end associate
        enddo
        !-----------------------------------------
        ! Group server and client processes
        !-----------------------------------------

        ! Assign IO clients to the same communicator as their related server process
        color = (PE_RANK_GLOBAL) / (num_PE/NUM_SERVERS)

        CALL MPI_Comm_dup( MPI_COMM_WORLD, globalComm, ierr )
        ! Group IO clients with their related server process. This is basically just grouping processes by node
        CALL MPI_Comm_split( globalComm, color, PE_RANK_GLOBAL, splitComm, ierr )

        do n = 1, n_nests
            associate(current_component => components(n))
                select type (current_component)
                    type is (domain_t)
                        CALL MPI_Comm_dup( splitComm, ioclient(n)%parent_comms, ierr )
                    type is (ioserver_t)
                        CALL MPI_Comm_dup( splitComm, current_component%client_comms, ierr )
                end select
            end associate
        enddo
        if (STD_OUT_PE) then
            write(*,*) "  Number of processing elements:          ",num_PE
            write(*,*) "  Number of compute elements:             ",NUM_COMPUTE
            write(*,*) "  Number of IO elements:                  ",NUM_SERVERS
            write(*,*) "  Number of processing elements per node: ",NUM_PROC_PER_NODE
            write(*,*) "  Number of IO processes per node:        ",NUM_IO_PER_NODE

#if defined(_OPENMP)
            write(*,*) "  Max number of OpenMP Threads:           ",num_threads
#endif
        endif


    end subroutine split_processes

    subroutine init_options(options, namelist_file, info_only, gen_nml, only_namelist_check)
        implicit none
        type(options_t), allocatable, intent(inout) :: options(:)
        character(len=*), intent(in) :: namelist_file
        logical, optional, intent(in) :: info_only, gen_nml, only_namelist_check

        type(general_options_type) :: dummy_general_options
        integer :: nests, i, name_unit, rc

        namelist /general/  nests

        ! We need at least one domain
        nests = 1

        ! If we are generating or printing namelist option information, we don't want to read anything, and only need to run the init routine once
        if ( .not.(info_only .or. gen_nml) )then
            ! First thing, read from options file how many nests we have
            call general_namelist(namelist_file,   dummy_general_options, 1, info_only=info_only, gen_nml=gen_nml)
            nests = dummy_general_options%nests
        endif

        allocate(options(nests))

        ! read in options file
        do i = 1, nests
            call options(i)%init(namelist_file, i, info_only=info_only, gen_nml=gen_nml)
            call collect_physics_requests(options(i))
            call options(i)%check()

            if (options(i)%general%parent_nest > 0) then
                ! Setup options%forcing to expect synthetic forcing from parent domain
                call options(i)%setup_synthetic_forcing()
            endif
        enddo

        call inter_nest_options_check(options)

        ! If this run was just done to check the namelist options, stop now
        if (only_namelist_check) then
            if (STD_OUT_PE) write(*,*) 'Namelist options check complete.'
            stop
        endif

    end subroutine init_options

    subroutine init_model(options,domain,boundary,ioclient,nest_indx)
        implicit none
        type(options_t), intent(inout) :: options(:)
        type(domain_t),  intent(inout) :: domain
        type(boundary_t),intent(inout) :: boundary ! forcing file for init conditions
        type(ioclient_t),intent(inout) :: ioclient
        integer, intent(in) :: nest_indx
        
        if (STD_OUT_PE .and. nest_indx == 1) write(*,*) "Initializing Domain"
        if (STD_OUT_PE .and. nest_indx == 1) flush(output_unit)
        call domain%init(options(nest_indx), nest_indx)

        if (options(nest_indx)%general%parent_nest == 0) then
            if (STD_OUT_PE) write(*,*) "Initializing boundary condition data structure"
            if (STD_OUT_PE) flush(output_unit)    
            call boundary%init(options(nest_indx), domain%vars_2d(domain%var_indx(kVARS%latitude)%v)%data_2d, domain%vars_2d(domain%var_indx(kVARS%longitude)%v)%data_2d)
        else
            call boundary%init(options(nest_indx), domain%vars_2d(domain%var_indx(kVARS%latitude)%v)%data_2d, domain%vars_2d(domain%var_indx(kVARS%longitude)%v)%data_2d, &
                                parent_options=options(options(nest_indx)%general%parent_nest))
        endif

        call ioclient%init(domain, boundary, options, nest_indx)

    end subroutine init_model

    subroutine init_model_state(options, domain, boundary, ioclient)
        implicit none
        type(options_t), intent(inout) :: options
        type(domain_t),  intent(inout) :: domain
        type(boundary_t),intent(inout) :: boundary ! forcing file for init conditions
        type(ioclient_t),intent(inout) :: ioclient

        if (STD_OUT_PE) write(*,*) "Receiving initial data"
        if (STD_OUT_PE) flush(output_unit)
        call ioclient%receive(boundary, domain)

        if (STD_OUT_PE) write(*,*) "Populating boundary object"
        if (STD_OUT_PE) flush(output_unit)
        call boundary%update_computed_vars(options, update=options%forcing%time_varying_z)

        if (STD_OUT_PE) write(*,*) "Initializing forcing interpolation"
        call domain%get_initial_conditions(boundary, options)    

        if (options%restart%restart) then
            if (STD_OUT_PE) write(*,*) "Reading restart data"
            call ioclient%receive_rst(domain, options)
        endif

        if (STD_OUT_PE) write(*,*) "Initializing physics"
        ! physics drivers need to be initialized after restart data are potentially read in.
        call init_physics(options, domain, boundary)
        if (STD_OUT_PE) flush(output_unit)

        if (.not.(options%restart%restart)) call ioclient%push(domain)

        if (size(options%general%child_nests) > 0) call ioclient%update_nest(domain)

    end subroutine init_model_state


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
        call mp_init(domain, options) !this could easily be moved to init_model...
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


    subroutine welcome_message()
        implicit none

        write(*,*) ""
        write(*,*) "============================================================"
        write(*,*) "|                                                          |"
        write(*,*) "|  The Intermediate Complexity Atmospheric Research Model  |"
        write(*,*) "|                       (ICAR/HICAR)                       |"
        write(*,*) "|                                                          |"
        write(*,*) "|   ICAR Developed at NCAR:                                |"
        write(*,*) "|     The National Center for Atmospheric Research         |"
        write(*,*) "|     NCAR is sponsored by the National Science Foundation |"
        write(*,*) "|   HICAR Developed at SLF:                                |"
        write(*,*) "|     The Swiss Institute for Snow and Avalanche Research  |"
        write(*,*) "|                                                          |"
        write(*,*) "|   Version: ",kVERSION_STRING,"                                         |"
        write(*,*) "|                                                          |"
        write(*,*) "============================================================"
        write(*,*) ""

    end subroutine welcome_message

    !> -------------------------------
    !! Call all physics driver var_request routines
    !!
    !! var_request routines allow physics modules to requested
    !! which variables they need to have allocated, advected, and written in restart files
    !!
    !! -------------------------------
    subroutine collect_physics_requests(options)
        type(options_t), intent(inout) :: options

        call ra_var_request(options)
        call lsm_var_request(options)
        call sfc_var_request(options)
        call pbl_var_request(options)
        call cu_var_request(options)
        call mp_var_request(options)
        call adv_var_request(options)
        call wind_var_request(options)

    end subroutine collect_physics_requests

end module
