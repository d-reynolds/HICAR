submodule(output_interface) output_implementation
  use output_metadata,          only : get_varmeta, get_varindx, get_varname
  use debug_module,             only : check_ncdf
  use iso_fortran_env,          only : output_unit
  use time_io,                  only : get_output_time, find_timestep_in_filelist
  use string,                   only : str, split_str, as_string
  use io_routines,              only : io_newunit
  implicit none

contains


    module subroutine init_output(this, options, its, ite, kts, kte, jts, jte, ide, kde, jde)
        implicit none
        class(output_t),  intent(inout)  :: this
        type(options_t),  intent(in)     :: options
        integer,          intent(in)     :: its, ite, kts, kte, jts, jte, ide, kde, jde
        
        integer :: i
        character(len=kMAX_STRING_LENGTH), allocatable :: tokens(:)
        character(len=kMAX_STRING_LENGTH) :: tmp_str

        allocate(this%variables(kINITIAL_VAR_SIZE))
        allocate(this%var_meta(kINITIAL_VAR_SIZE))

        this%n_vars = 0
        this%n_dims      = 0
        this%is_initialized = .True.

        this%output_counter = 1
        
        this%its = its; this%ite = ite; this%kts = kts; this%kte = kte; this%jts = jts; this%jte = jte
        this%file_dim_len = (/ide, jde, kde/)

        call set_attrs(this, options)
        
        ! Set the output filename to be the location of the output folder, and the name of the initial conditions file

        ! Get the name of the base initial conditions file. Could contain a path, so we need to strip it off
        ! Using occurarnces of "/" in the string to split the string. The last token is the filename
        tokens = split_str(trim(options%domain%init_conditions_file), "/")

        ! store the base file name temporarily here.
        tmp_str = trim(tokens(size(tokens)))

        ! Set the output filename to be the location of the output folder, and the name of the initial conditions file, with ".nc" removed
        this%base_out_file_name = trim(options%output%output_folder) // '/' // tmp_str(1:(len_trim(tmp_str)-3)) // "_"
        this%output_count = options%output%frames_per_outfile
        this%base_rst_file_name = trim(options%restart%restart_folder) // '/' // tmp_str(1:(len_trim(tmp_str)-3)) // "_"
        
        write(this%output_fn, '(A,A,".nc")')    &
                trim(this%base_out_file_name),   &
                trim(as_string(options%general%start_time,this%file_date_format))

        call this%add_variables(options%output%vars_for_output + options%vars_for_restart)
        
    end subroutine init_output
    
    !If this is a restart run, set output counter and filename to pick up
    !where we left off
    module subroutine init_restart(this, options, par_comms, out_var_indices)
        implicit none
        class(output_t),  intent(inout)  :: this
        type(options_t),  intent(in)     :: options
        type(MPI_Comm),   intent(in)     :: par_comms
        integer,          intent(in)     :: out_var_indices(:)
        
        integer                       :: error
        character(len=kMAX_FILE_LENGTH), allocatable :: file_list(:)

        !Expect the worse
        error = 1

        call get_outputfiles(this,file_list)
        
        if (size(file_list) > 0) then
            !Find output file and step in output filelist
            this%output_counter = find_timestep_in_filelist(file_list, 'time', options%restart%restart_time,this%output_fn,error=error)
        endif
        
        !If there was no error getting output step, then a file already exists which contains our time step
        if (error == 0) then
            !Open output file, setting out_ncfile_id
            call open_file(this, this%output_fn, options%restart%restart_time, par_comms, out_var_indices)
            this%out_ncfile_id = this%active_nc_id

            ! if this is a restart run, then there will be no initial call to save_out_file, meaning
            ! that the output counter should be auto-incremented to 2 (1 + 1) here 
            this%output_counter = this%output_counter + 1
        ! if there was an error getting the output step, then we assume that there is no file for the current output time step
        else
            write(*,*) "No existing output file found for restart run. Technically, this should never happen."
            write(*,*) "If you wanted to restart on purpose without any existing output files, just be aware."
            write(*,*) "We will set the output counter to 1 and start producing a new output file."

            ! set back to 1, in case this was set when calling find_timestep_in_filelist
            this%output_counter = 1
            write(this%output_fn, '(A,A,".nc")')    &
                trim(this%base_out_file_name),   &
                trim(as_string(options%restart%restart_time,this%file_date_format))
        endif

    end subroutine init_restart
    
    subroutine get_outputfiles(this,file_list)
        implicit none
        class(output_t),  intent(in)   :: this
        character(len=*), allocatable, intent(inout):: file_list(:)
        
        integer :: i, error, nfiles, unit, PE_RANK_GLOBAL
        character(len=kMAX_FILE_LENGTH)  :: file
        character(len=kMAX_FILE_LENGTH), allocatable :: temp_list(:)
        character(len=kMAX_FILE_LENGTH) :: temp_file
        character(len=kMAX_FILE_LENGTH) :: cmd_str
        
        allocate(temp_list(MAX_NUMBER_FILES))
        call MPI_Comm_Rank(MPI_COMM_WORLD,PE_RANK_GLOBAL)

        temp_file = '.tmp_outfiles'//trim(str(PE_RANK_GLOBAL))//'.txt'
        cmd_str = 'ls '//trim(this%base_out_file_name)//'*.nc > '//trim(temp_file)

        call EXECUTE_COMMAND_LINE( cmd_str, EXITSTAT=error )

        if (error /= 0) then
            return
        endif
        
        open( io_newunit(unit), file = trim(temp_file) )
        i=0
        error=0
        do while (error==0)
            read(unit, '(A)', iostat=error) file
            if (error==0) then
                i=i+1
                temp_list(i) = trim(file)
            endif
        enddo
        close(unit)
        
        cmd_str = 'rm '//trim(temp_file)
        call system( cmd_str )
        
        nfiles = i
        allocate(file_list(nfiles))
        file_list(1:nfiles) = temp_list(1:nfiles)
        deallocate(temp_list)
    end subroutine

    !> -------------------------------
    !! Populare the metadata structure
    !!
    !! -------------------------------
    module subroutine set_attrs(this, options)
        class(output_t),  intent(inout)  :: this
        type(options_t), intent(in)    :: options
        character*60 :: a_string

        call this%add_attribute("comment",options%general%comment)
        call this%add_attribute("source","ICAR version:"//trim(options%general%version))

        ! Add info on grid setting:
        write(a_string,*) options%domain%sleve
        call this%add_attribute("sleve",a_string)
        if (options%domain%sleve) then
          write(a_string,*) options%domain%terrain_smooth_windowsize
          call this%add_attribute("terrain_smooth_windowsize",a_string )
          write(a_string,*) options%domain%terrain_smooth_cycles
          call this%add_attribute("terrain_smooth_cycles",a_string )
          write(a_string,*) options%domain%decay_rate_L_topo
          call this%add_attribute("decay_rate_L_topo",a_string )
          write(a_string,*) options%domain%decay_rate_s_topo
          call this%add_attribute("decay_rate_S_topo",a_string )
          write(a_string,*) options%domain%sleve_n
          call this%add_attribute("sleve_n",a_string )
        endif
        ! Add some more info on physics settings:
        write(a_string,*) options%physics%boundarylayer
        call this%add_attribute("pbl", a_string )
        write(a_string,*) options%physics%landsurface
        call this%add_attribute("lsm", a_string )
        write(a_string,*) options%physics%watersurface
        call this%add_attribute("water", a_string )
        write(a_string,*) options%physics%microphysics
        call this%add_attribute("mp", a_string )
        write(a_string,*) options%physics%radiation
        call this%add_attribute("rad", a_string )
        write(a_string,*) options%physics%convection
        call this%add_attribute("conv", a_string )
        write(a_string,*) options%physics%advection
        call this%add_attribute("adv", a_string )
        write(a_string,*) options%physics%windtype
        call this%add_attribute("wind", a_string )


        call this%add_attribute("ids",str(1))
        call this%add_attribute("ide",str(this%file_dim_len(1)))
        call this%add_attribute("jds",str(1))
        call this%add_attribute("jde",str(this%file_dim_len(2)))
        call this%add_attribute("kds",str(1))
        call this%add_attribute("kde",str(this%file_dim_len(3)))

        !call this%add_attribute("ims",str(this%grid%ims))
        !call this%add_attribute("ime",str(this%grid%ime))
        !call this%add_attribute("jms",str(this%grid%jms))
        !call this%add_attribute("jme",str(this%grid%jme))
        !call this%add_attribute("kms",str(this%grid%kms))
        !call this%add_attribute("kme",str(this%grid%kme))

        !call this%add_attribute("its",str(this%grid%its))
        !call this%add_attribute("ite",str(this%grid%ite))
        !call this%add_attribute("jts",str(this%grid%jts))
        !call this%add_attribute("jte",str(this%grid%jte))
        !call this%add_attribute("kts",str(this%grid%kts))
        !call this%add_attribute("kte",str(this%grid%kte))

    end subroutine

    module subroutine add_to_output(this, in_variable)
        class(output_t),   intent(inout) :: this
        type(meta_data_t),  intent(in) :: in_variable
        
        type(variable_t) :: variable
        type(meta_data_t) :: var_meta
        integer :: n, nx, ny, nz
        
        nx = this%ite - this%its + 1 + in_variable%xstag
        ny = this%jte - this%jts + 1 + in_variable%ystag
        nz = this%kte - this%kts ! no "+1" since the writer nz bounds are domain%kts to domain%kte+1 by default
        var_meta = in_variable

        if (var_meta%three_d) then
            if (var_meta%dimensions(2)== "level_i") nz = nz+1

            if (var_meta%dim_len(2)>0) nz = var_meta%dim_len(2)

            call variable%initialize(var_meta%id, (/nx, nz, ny/))
            variable%data_3d = kEMPT_BUFF
            var_meta%dim_len = (/ nx, nz, ny /)
        else
            call variable%initialize(var_meta%id, (/nx, ny/))
            variable%data_2d = kEMPT_BUFF
            var_meta%dim_len = (/ nx, ny /)
        endif
        if (this%n_vars == size(this%variables)) call this%increase_var_capacity()

        this%n_vars = this%n_vars + 1
        this%variables(this%n_vars) = variable
        this%var_meta(this%n_vars) = var_meta

    end subroutine

    module subroutine save_out_file(this, time, par_comms, out_var_indices)
        class(output_t),  intent(inout) :: this
        type(Time_type),  intent(in)  :: time
        type(MPI_Comm),   intent(in)     :: par_comms
        integer,          intent(in)  :: out_var_indices(:)


        !Check if we should change the file
        if (this%output_counter > this%output_count) then
            write(this%output_fn, '(A,A,".nc")')    &
                trim(this%base_out_file_name),   &
                trim(as_string(time,this%file_date_format))

            this%output_counter = 1
            if (this%out_ncfile_id > 0) then
                call check_ncdf(nf90_close(this%out_ncfile_id), "Closing output file ")
                this%out_ncfile_id = -1
            endif
        endif

        this%active_nc_id = this%out_ncfile_id
        if (this%active_nc_id < 0) then
            call open_file(this, this%output_fn, time, par_comms,out_var_indices)
            this%out_ncfile_id = this%active_nc_id
        endif
        call flush(output_unit)

        if (.not.this%block_checked) call block_hunter(this)
        ! store output
        call save_data(this, this%output_counter, time, out_var_indices)

        !In case we had creating set to true, set to false
        this%creating = .false.
        
        this%active_nc_id = -1
                
        this%output_counter  = this%output_counter + 1

        if (this%out_ncfile_id > 0) then
            call check_ncdf(nf90_close(this%out_ncfile_id), "Closing output file ")
            this%out_ncfile_id = -1
        endif

    end subroutine
    
    module subroutine save_rst_file(this, time, par_comms, rst_var_indices)
        class(output_t),  intent(inout) :: this
        type(Time_type),  intent(in)  :: time
        type(MPI_Comm),   intent(in)     :: par_comms
        integer,          intent(in)  :: rst_var_indices(:)

        this%active_nc_id = this%rst_ncfile_id
        
        if (this%active_nc_id < 0) then
            write(this%restart_fn, '(A,A,".nc")')    &
                trim(this%base_rst_file_name),   &
                trim(as_string(time,this%file_date_format))

            call open_file(this, this%restart_fn, time, par_comms,rst_var_indices)
            this%rst_ncfile_id = this%active_nc_id
        endif

        ! store output
        call save_data(this, 1, time, rst_var_indices)
        
        !In case we had creating set to true, set to false
        this%creating = .false.

        !Close the file
        call check_ncdf(nf90_close(this%rst_ncfile_id), "Closing restart file ")

        this%rst_ncfile_id = -1
        this%active_nc_id = -1

    end subroutine

    
    !See if this outputer has 'blocks', or a piece of it which it should not write
    !This can occur due to using only images on a given node, which results in stepped
    !patterns if the output is not appropriately treated
    !
    !The goal of the routine is to set all start and cnt object fields so that either
    !(1) if there is no block, the standard fields contain the whole bounds of the object
    !(2) if there is a block, the standard fields contain the continuous block of the object
    ! and the 'b' fields contain the discontinuous block
    subroutine block_hunter(this)
        class(output_t),  intent(inout) :: this

        integer :: i_s, i_e, k_s, k_e, j_s, j_e, nx, ny, nz
        integer :: i_s_b, i_e_b, j_s_b, j_e_b
        integer :: i_s_b2, i_e_b2, j_s_b2, j_e_b2
        logical :: found_var = .False.
        integer :: i
        real, allocatable, dimension(:,:) :: datas
        
        this%block_checked = .True.
        this%is_blocked = .False.
        this%blocked_LL = .False.
        this%blocked_UR = .False.

        i_s_b = 2; i_e_b = 1; j_s_b = 2; j_e_b = 1;
        i_s_b2 = 2; i_e_b2 = 1; j_s_b2 = 2; j_e_b2 = 1;

        i_s = this%its
        i_e = this%ite
        j_s = this%jts
        j_e = this%jte
        k_s = this%kts
        k_e = this%kte

        nx = i_e - i_s + 1
        ny = j_e - j_s + 1
        allocate(datas(nx,ny))

        do i = 1,size(this%variables)
            if ( (this%variables(i)%xstag==0) .and. (this%variables(i)%ystag==0) ) then
                if (this%variables(i)%three_d) then
                    if (ALL(this%variables(i)%data_3d(1:nx,1,1:ny) == kEMPT_BUFF)) cycle
                    datas = this%variables(i)%data_3d(1:nx,1,1:ny)
                    found_var = .True.
                    exit
                else if (this%variables(i)%two_d) then
                    if (ALL(this%variables(i)%data_2d == kEMPT_BUFF)) cycle
                    datas = this%variables(i)%data_2d(1:nx,1:ny)
                    found_var = .True.
                    exit
                endif
            endif
        enddo

        if (found_var .eqv. .False.) then
            write(*,*) "No output variable found which is not on a staggered grid."
            write(*,*) "Cannot determine block. This should never happen"
            flush(output_unit)
            stop "block_hunter function in output_obj failed to find a variable"
        endif

        !Check each corner to see where this starts
        !LL
        if (datas(1,1) == kEMPT_BUFF) then
            i_s_b = findloc(datas(:,1),kEMPT_BUFF,dim=1,back=.True.) + 1
            i_e_b = i_e
            j_e_b = findloc(datas(1,:),kEMPT_BUFF,dim=1,back=.True.) + j_s - 1
            j_s_b = j_s
            this%blocked_LL = .True.
        endif
        !UR
        if (datas(nx,ny) == kEMPT_BUFF) then
            i_s_b2 = i_s 
            i_e_b2 = findloc(datas(:,ny),kEMPT_BUFF,dim=1,back=.False.) - 1
            j_e_b2 = j_e 
            j_s_b2 = findloc(datas(nx,:),kEMPT_BUFF,dim=1,back=.False.) + j_s - 1
            this%blocked_UR = .True.
        endif
        
        if (this%blocked_LL) j_s = j_e_b+1
        if (this%blocked_UR) j_e = j_s_b2-1
        
        this%start_3d = (/ i_s, j_s, k_s /)
        nx = max((i_e-i_s+1),0)
        ny = max((j_e-j_s+1),0)
        nz = (k_e-k_s+1)
        this%cnt_3d = (/ nx, ny, nz  /)
        this%cnt_2d = (/ nx, ny /)

        if (ANY( (/nx, ny, nz/) == 0)) then
            this%cnt_3d = 0
            this%cnt_2d = 0
            this%start_3d = (/ 1, 1, 1 /)
        endif


        !Compute block start and cnts accordingly
        if (this%blocked_LL) then
            this%start_3d_b = (/ i_s_b, j_s_b, k_s /)
            this%cnt_3d_b = (/ (i_e_b-i_s_b+1), (j_e_b-j_s_b+1), (k_e-k_s+1)  /)
            this%cnt_2d_b = (/ (i_e_b-i_s_b+1), (j_e_b-j_s_b+1) /)
        else
            this%start_3d_b = (/ 1, 1, 1 /)
            this%cnt_3d_b = (/ 0, 0, 0  /)
            this%cnt_2d_b = (/ 0, 0 /)
        endif 
        if (this%blocked_UR) then
            this%start_3d_b2 = (/ i_s_b2, j_s_b2, k_s /)
            this%cnt_3d_b2 = (/ (i_e_b2-i_s_b2+1), (j_e_b2-j_s_b2+1), (k_e-k_s+1)  /)
            this%cnt_2d_b2 = (/ (i_e_b2-i_s_b2+1), (j_e_b2-j_s_b2+1) /)
        else
            this%start_3d_b2 = (/ 1, 1, 1 /)
            this%cnt_3d_b2 = (/ 0, 0, 0  /)
            this%cnt_2d_b2 = (/ 0, 0 /)
        endif
    end subroutine block_hunter
    
    subroutine open_file(this, filename, time, par_comms, var_indx_list)
        class(output_t),                 intent(inout) :: this
        character(len=kMAX_FILE_LENGTH), intent(in)    :: filename

        type(Time_type),                 intent(in)    :: time
        type(MPI_Comm),                  intent(in)    :: par_comms
        integer,                         intent(in)    :: var_indx_list(:)
        type(MPI_Info) :: par_comm_info
        integer :: err
        
        ! open file
        call MPI_Comm_get_info(par_comms, par_comm_info)
        err = nf90_open(filename, IOR(NF90_WRITE,NF90_NETCDF4), this%active_nc_id, &
                comm = par_comms%MPI_VAL, info = par_comm_info%MPI_VAL)
        if (err /= NF90_NOERR) then
            call check_ncdf( nf90_create(filename, IOR(NF90_CLOBBER,NF90_NETCDF4), this%active_nc_id, &
                    comm = par_comms%MPI_VAL, info = par_comm_info%MPI_VAL), "Opening:"//trim(filename))
        else
            ! in case we need to add a new variable when setting up variables
            call check_ncdf(nf90_redef(this%active_nc_id), "Setting redefine mode for: "//trim(filename))
        endif

        !This command cuts down on write time and must be done once the file is opened
        call check_ncdf( nf90_set_fill(this%active_nc_id, nf90_nofill, err), "Setting fill mode to none")

        ! define variables or find variable IDs (and dimensions)
        call setup_variables(this, time, var_indx_list)

        if (this%creating) then
            ! add global attributes such as the image number, domain dimension, creation time
            call add_global_attributes(this)
        endif

        !Set creatint to true, so that we will output any static variables on the first pass.
        this%creating=.True.

        ! End define mode. This tells netCDF we are done defining metadata.
        call check_ncdf( nf90_enddef(this%active_nc_id), "end define mode" )
    

    end subroutine

    module subroutine add_variables(this, vars_to_out)
        class(output_t),  intent(inout)  :: this
        integer, dimension(:), intent(in):: vars_to_out

        integer :: i

        !Loop through domain vars_to_out, get the var index for the given variable name, and add that var's meta data to local list
        do i = 1,kMAX_STORAGE_VARS
            if (vars_to_out(i)>0) then
                !i should equal the kVARS index of the variable we want to output
                call this%add_to_output(get_varmeta(i))
            endif
        enddo
    end subroutine

    subroutine add_global_attributes(this)
        implicit none
        class(output_t), intent(inout)  :: this
        integer :: i

        character(len=19)       :: todays_date_time
        integer,dimension(8)    :: date_time
        character(len=49)       :: date_format
        character(len=5)        :: UTCoffset
        character(len=64)       :: err
        integer                 :: ncid

        ncid = this%active_nc_id

        err="Creating global attributes"
        call check_ncdf( nf90_put_att(ncid,NF90_GLOBAL,"Conventions","CF-1.6"), trim(err))
        call check_ncdf( nf90_put_att(ncid,NF90_GLOBAL,"title","Intermediate Complexity Atmospheric Research (ICAR) model output"), trim(err))
        call check_ncdf( nf90_put_att(ncid,NF90_GLOBAL,"institution","National Center for Atmospheric Research"), trim(err))
        call check_ncdf( nf90_put_att(ncid,NF90_GLOBAL,"history","Created:"//todays_date_time//UTCoffset), trim(err))
        call check_ncdf( nf90_put_att(ncid,NF90_GLOBAL,"references", &
                    "Gutmann et al. 2016: The Intermediate Complexity Atmospheric Model (ICAR). J.Hydrometeor. doi:10.1175/JHM-D-15-0155.1, 2016."), trim(err))
        call check_ncdf( nf90_put_att(ncid,NF90_GLOBAL,"contact","Ethan Gutmann : gutmann@ucar.edu"), trim(err))
        call check_ncdf( nf90_put_att(ncid,NF90_GLOBAL,"git",VERSION), trim(err))

        if (this%n_attrs > 0) then
            do i=1,this%n_attrs
                call check_ncdf( nf90_put_att(   ncid,             &
                                            NF90_GLOBAL,                &
                                            trim(this%attributes(i)%name),    &
                                            trim(this%attributes(i)%value)),  &
                                            "global attr:"//trim(this%attributes(i)%name))
            enddo
        endif

        call date_and_time(values=date_time,zone=UTCoffset)
        date_format='(I4,"/",I2.2,"/",I2.2," ",I2.2,":",I2.2,":",I2.2)'
        write(todays_date_time,date_format) date_time(1:3),date_time(5:7)

        call check_ncdf(nf90_put_att(ncid, NF90_GLOBAL,"history","Created:"//todays_date_time//UTCoffset), "global attr")
        !call check_ncdf(nf90_put_att(ncid, NF90_GLOBAL, "image", this_image()))

    end subroutine add_global_attributes

    subroutine setup_variables(this, time, var_indx_list)
        implicit none
        class(output_t), intent(inout) :: this
        type(Time_type), intent(in)    :: time
        integer, intent(in)            :: var_indx_list(:)
        integer :: i

        ! iterate through variables creating or setting up variable IDs if they exist, also dimensions
        do i=1,size(var_indx_list)            
            call setup_dims_for_var(this, this%var_meta(var_indx_list(i)))
            call setup_variable(this, this%var_meta(var_indx_list(i)))
        end do

        call setup_time_variable(this, time)

    end subroutine setup_variables

    subroutine setup_time_variable(this, time)
        implicit none
        class(output_t), intent(inout) :: this
        type(Time_type), intent(in)    :: time
        integer :: err
        character(len=kMAX_NAME_LENGTH) :: calendar

        associate(var => this%time)
        var%name = "time"
        var%dimensions = [ "time" ]

        select case (time%calendar)
            case(GREGORIAN)
                calendar = "proleptic_gregorian"
            case(NOLEAP)
                calendar = "noleap"
            case(THREESIXTY)
                calendar = "360_day"
            case default
                calendar = "standard"
        end select


        err = nf90_inq_varid(this%active_nc_id, var%name, var%file_var_id)

        ! if the variable was not found in the netcdf file then we will define it.
        if (err /= NF90_NOERR) then

            if (allocated(var%dim_ids)) deallocate(var%dim_ids)
            allocate(var%dim_ids(1))

            ! Try to find the dimension ID if it exists already.
            err = nf90_inq_dimid(this%active_nc_id, trim(var%dimensions(1)), var%dim_ids(1))

            ! if the dimension doesn't exist in the file, create it.
            if (err/=NF90_NOERR) then
                call check_ncdf( nf90_def_dim(this%active_nc_id, trim(var%dimensions(1)), NF90_UNLIMITED, &
                            var%dim_ids(1) ), "def_dim"//var%dimensions(1) )
            endif

            write(this%time_units,kOUTPUT_FMT) time%year_zero,time%month_zero,time%day_zero,time%hour_zero

            call check_ncdf( nf90_def_var(this%active_nc_id, var%name, NF90_DOUBLE, var%dim_ids(1), var%file_var_id), "Defining time" )
            call check_ncdf( nf90_put_att(this%active_nc_id, var%file_var_id,"standard_name","time"))
            call check_ncdf( nf90_put_att(this%active_nc_id, var%file_var_id,"calendar",trim(calendar)))
            call check_ncdf( nf90_put_att(this%active_nc_id, var%file_var_id,"units",this%time_units))
            call check_ncdf( nf90_put_att(this%active_nc_id, var%file_var_id,"UTCoffset","0"))

        else
            err = nf90_get_att(this%active_nc_id, var%file_var_id, "units", this%time_units)
        endif
        end associate

    end subroutine setup_time_variable


    subroutine save_data(this, current_step, time, var_indx_list)
        implicit none
        class(output_t), intent(in) :: this
        integer,         intent(in) :: current_step
        type(Time_type), intent(in) :: time
        integer, intent(in)         :: var_indx_list(:)
        
        type(Time_type) :: output_time
        real, allocatable :: var_3d(:,:,:)
        integer :: i, k_s, k_e
        
        integer :: start_three_D_t(4), start_three_D_t_b(4), start_three_D_t_b2(4)
        integer :: start_two_D_t(3), start_two_D_t_b(3), start_two_D_t_b2(3)
        integer :: cnt_3d(3), cnt_3d_b(3), cnt_3d_b2(3)
        integer :: cnt_2d(2)
        integer :: v_i_s, v_i_e, v_j_s, v_j_e
        integer :: v_i_s_b, v_i_e_b, v_j_s_b, v_j_e_b
        integer :: v_i_s_b2, v_i_e_b2, v_j_s_b2, v_j_e_b2
    
        do i=1,size(var_indx_list)
            associate(var => this%var_meta(var_indx_list(i)))                
                k_s = this%kts
                k_e = var%dim_len(2)
                
                start_three_D_t = (/ this%start_3d(1), this%start_3d(2), this%start_3d(3), current_step /)
                cnt_3d = (/ this%cnt_3d(1), this%cnt_3d(2), (k_e-k_s+1) /)
                    
                start_three_D_t_b = (/ this%start_3d_b(1), this%start_3d_b(2), this%start_3d_b(3), current_step /)
                cnt_3d_b = (/ this%cnt_3d_b(1), this%cnt_3d_b(2), (k_e-k_s+1) /)
                        
                start_three_D_t_b2 = (/ this%start_3d_b2(1), this%start_3d_b2(2), this%start_3d_b2(3), current_step /)
                cnt_3d_b2 = (/ this%cnt_3d_b2(1), this%cnt_3d_b2(2), (k_e-k_s+1) /)
                    
                start_two_D_t = (/ this%start_3d(1), this%start_3d(2), current_step /)
                start_two_D_t_b = (/ this%start_3d_b(1), this%start_3d_b(2), current_step /)
                start_two_D_t_b2 = (/ this%start_3d_b2(1), this%start_3d_b2(2), current_step /)

                if (this%ite == this%file_dim_len(1)  .and. cnt_3d(1) > 0) cnt_3d(1) = cnt_3d(1) + var%xstag
                if ((this%start_3d_b(1) - this%its + cnt_3d_b(1)) == this%file_dim_len(1)) cnt_3d_b(1) = cnt_3d_b(1) + var%xstag
                if ((this%start_3d_b2(1) - this%its + cnt_3d_b2(1)) == this%file_dim_len(1)) cnt_3d_b2(1) = cnt_3d_b2(1) + var%xstag
                
                if (this%jte == this%file_dim_len(2) .and. cnt_3d(2) > 0) cnt_3d(2) = cnt_3d(2) + var%ystag
                if ((this%start_3d_b(2) - this%its + cnt_3d_b(2)) == this%file_dim_len(2)) cnt_3d_b(2) = cnt_3d_b(2) + var%ystag
                if ((this%start_3d_b2(2) - this%jts + cnt_3d_b2(2)) == this%file_dim_len(2)) cnt_3d_b2(2) = cnt_3d_b2(2) + var%ystag

                v_i_s = this%start_3d(1) - this%its + 1
                v_i_e = v_i_s + cnt_3d(1) - 1
        
                v_j_s = this%start_3d(2) - this%jts + 1
                v_j_e = v_j_s + cnt_3d(2) - 1

                v_i_s_b = this%start_3d_b(1) - this%its + 1
                v_i_e_b = v_i_s_b + cnt_3d_b(1) - 1
        
                v_j_s_b = this%start_3d_b(2) - this%jts + 1
                v_j_e_b = v_j_s_b + cnt_3d_b(2) - 1
        
                v_i_s_b2 = this%start_3d_b2(1) - this%its + 1
                v_i_e_b2 = v_i_s_b2 + cnt_3d_b2(1) - 1
        
                v_j_s_b2 = this%start_3d_b2(2) - this%jts + 1
                v_j_e_b2 = v_j_s_b2 + cnt_3d_b2(2) - 1
                
                !Safety checks so that we don't index outside of the lower array bound (can happen if no block)
                v_i_s = max(v_i_s,1)
                v_i_e = max(v_i_e,1)

                v_j_s = max(v_j_s,1)
                v_j_e = max(v_j_e,1)

                v_i_s_b = max(v_i_s_b,1)
                v_i_e_b = max(v_i_e_b,1)

                v_j_s_b = max(v_j_s_b,1)
                v_j_e_b = max(v_j_e_b,1)
                
                v_i_s_b2 = max(v_i_s_b2,1)
                v_i_e_b2 = max(v_i_e_b2,1)

                v_j_s_b2 = max(v_j_s_b2,1)
                v_j_e_b2 = max(v_j_e_b2,1)
                call check_ncdf( nf90_var_par_access(this%active_nc_id, var%file_var_id, nf90_collective))

                if (var%three_d) then
                    var_3d = reshape(this%variables(var_indx_list(i))%data_3d, &
                        shape=(/ ubound(this%variables(var_indx_list(i))%data_3d,1), ubound(this%variables(var_indx_list(i))%data_3d,3), ubound(this%variables(var_indx_list(i))%data_3d,2) /), order=[1,3,2])
                    if (var%unlimited_dim) then
                        call check_ncdf( nf90_put_var(this%active_nc_id, var%file_var_id, &
                            var_3d(v_i_s:v_i_e,v_j_s:v_j_e,k_s:k_e), start_three_D_t, count=(/cnt_3d(1), cnt_3d(2), cnt_3d(3), 1/)), "saving:"//trim(var%name) )
                        call check_ncdf( nf90_put_var(this%active_nc_id, var%file_var_id, &
                            var_3d(v_i_s_b:v_i_e_b,v_j_s_b:v_j_e_b,k_s:k_e), start_three_D_t_b, count=(/cnt_3d_b(1), cnt_3d_b(2), cnt_3d_b(3), 1/)), "saving LL block:"//trim(var%name) )
                        call check_ncdf( nf90_put_var(this%active_nc_id, var%file_var_id, &
                            var_3d(v_i_s_b2:v_i_e_b2,v_j_s_b2:v_j_e_b2,k_s:k_e), start_three_D_t_b2, count=(/cnt_3d_b2(1), cnt_3d_b2(2), cnt_3d_b2(3), 1/)), "saving UR block:"//trim(var%name) )
                    elseif (this%creating) then
                        call check_ncdf( nf90_put_var(this%active_nc_id, var%file_var_id, var_3d(v_i_s:v_i_e,v_j_s:v_j_e,k_s:k_e), &
                            start=this%start_3d,count=cnt_3d ), "saving:"//trim(var%name) )
                        call check_ncdf( nf90_put_var(this%active_nc_id, var%file_var_id, var_3d(v_i_s_b:v_i_e_b,v_j_s_b:v_j_e_b,k_s:k_e), &
                            start=this%start_3d_b,count=cnt_3d_b ), "saving LL block:"//trim(var%name) )
                        call check_ncdf( nf90_put_var(this%active_nc_id, var%file_var_id, &
                            var_3d(v_i_s_b2:v_i_e_b2,v_j_s_b2:v_j_e_b2,k_s:k_e), start=this%start_3d_b2,&
                            count=cnt_3d_b2 ), "saving UR block:"//trim(var%name) )
                    endif
                elseif (var%two_d) then
                    if (var%dtype == kREAL) then
                        if (var%unlimited_dim) then
                            call check_ncdf( nf90_put_var(this%active_nc_id, var%file_var_id,  this%variables(var_indx_list(i))%data_2d(v_i_s:v_i_e,v_j_s:v_j_e), &
                                    start_two_D_t,count=(/ cnt_3d(1), cnt_3d(2), 1/)), "saving:"//trim(var%name) )
                            call check_ncdf( nf90_put_var(this%active_nc_id, var%file_var_id,  this%variables(var_indx_list(i))%data_2d(v_i_s_b:v_i_e_b,v_j_s_b:v_j_e_b), &
                                start_two_D_t_b,count=(/ cnt_3d_b(1), cnt_3d_b(2), 1/)), "saving LL block:"//trim(var%name) )
                            call check_ncdf( nf90_put_var(this%active_nc_id, var%file_var_id,  this%variables(var_indx_list(i))%data_2d(v_i_s_b2:v_i_e_b2,v_j_s_b2:v_j_e_b2), &
                                start_two_D_t_b2,count=(/ cnt_3d_b2(1), cnt_3d_b2(2), 1/)), "saving UR block:"//trim(var%name) )
                        elseif (this%creating) then
                            call check_ncdf( nf90_put_var(this%active_nc_id, var%file_var_id,  this%variables(var_indx_list(i))%data_2d(v_i_s:v_i_e,v_j_s:v_j_e), &
                                        start=(/ this%start_3d(1), this%start_3d(2) /), &
                                        count=(/ cnt_3d(1), cnt_3d(2) /)), "saving:"//trim(var%name) )
                            call check_ncdf( nf90_put_var(this%active_nc_id, var%file_var_id,  this%variables(var_indx_list(i))%data_2d(v_i_s_b:v_i_e_b,v_j_s_b:v_j_e_b), &
                                    start=(/ this%start_3d_b(1), this%start_3d_b(2) /), &
                                    count=(/ cnt_3d_b(1), cnt_3d_b(2) /)), "saving LL block:"//trim(var%name) )
                            call check_ncdf( nf90_put_var(this%active_nc_id, var%file_var_id,  this%variables(var_indx_list(i))%data_2d(v_i_s_b2:v_i_e_b2,v_j_s_b2:v_j_e_b2), &
                                    start=(/ this%start_3d_b2(1), this%start_3d_b2(2) /), &
                                    count=(/ cnt_3d_b2(1), cnt_3d_b2(2) /)), "saving UR block:"//trim(var%name) )
                        endif    
                    else if (var%dtype == kINTEGER) then
                        if (var%unlimited_dim) then
                            call check_ncdf( nf90_put_var(this%active_nc_id, var%file_var_id,  INT(this%variables(var_indx_list(i))%data_2d(v_i_s:v_i_e,v_j_s:v_j_e)), &
                                    start_two_D_t,count=(/ cnt_3d(1), cnt_3d(2), 1/)), "saving:"//trim(var%name) )
                            call check_ncdf( nf90_put_var(this%active_nc_id, var%file_var_id,  INT(this%variables(var_indx_list(i))%data_2d(v_i_s_b:v_i_e_b,v_j_s_b:v_j_e_b)), &
                                start_two_D_t_b,count=(/ cnt_3d_b(1), cnt_3d_b(2), 1/)), "saving LL block:"//trim(var%name) )
                            call check_ncdf( nf90_put_var(this%active_nc_id, var%file_var_id,  INT(this%variables(var_indx_list(i))%data_2d(v_i_s_b2:v_i_e_b2,v_j_s_b2:v_j_e_b2)), &
                                start_two_D_t_b2,count=(/ cnt_3d_b2(1), cnt_3d_b2(2), 1/)), "saving UR block:"//trim(var%name) )
                        elseif (this%creating) then
                            call check_ncdf( nf90_put_var(this%active_nc_id, var%file_var_id,  INT(this%variables(var_indx_list(i))%data_2d(v_i_s:v_i_e,v_j_s:v_j_e)), &
                                        start=(/ this%start_3d(1), this%start_3d(2) /), &
                                        count=(/ cnt_3d(1), cnt_3d(2) /)), "saving:"//trim(var%name) )
                            call check_ncdf( nf90_put_var(this%active_nc_id, var%file_var_id,  INT(this%variables(var_indx_list(i))%data_2d(v_i_s_b:v_i_e_b,v_j_s_b:v_j_e_b)), &
                                    start=(/ this%start_3d_b(1), this%start_3d_b(2) /), &
                                    count=(/ cnt_3d_b(1), cnt_3d_b(2) /)), "saving LL block:"//trim(var%name) )
                            call check_ncdf( nf90_put_var(this%active_nc_id, var%file_var_id,  INT(this%variables(var_indx_list(i))%data_2d(v_i_s_b2:v_i_e_b2,v_j_s_b2:v_j_e_b2)), &
                                    start=(/ this%start_3d_b2(1), this%start_3d_b2(2) /), &
                                    count=(/ cnt_3d_b2(1), cnt_3d_b2(2) /)), "saving UR block:"//trim(var%name) )
                        endif    
    
                    ! elseif (var%dtype == kDOUBLE) then
                    !     if (var%unlimited_dim) then
                    !         call check_ncdf( nf90_put_var(this%active_nc_id, var%file_var_id,  var%data_2dd(v_i_s:v_i_e,v_j_s:v_j_e), &
                    !                 start_two_D_t,count=(/ cnt_3d(1), cnt_3d(2), 1/)), "saving:"//trim(var%name) )
                    !         if (this%blocked_UR .or. this%blocked_LL) then
                    !             call check_ncdf( nf90_put_var(this%active_nc_id, var%file_var_id,  var%data_2dd(v_i_s_b:v_i_e_b,v_j_s_b:v_j_e_b), &
                    !                 start_two_D_t_b,count=(/ cnt_3d_b(1), cnt_3d_b(2), 1/)), "saving LL block:"//trim(var%name) )
                    !             call check_ncdf( nf90_put_var(this%active_nc_id, var%file_var_id,  var%data_2dd(v_i_s_b2:v_i_e_b2,v_j_s_b2:v_j_e_b2), &
                    !                 start_two_D_t_b2,count=(/ cnt_3d_b2(1), cnt_3d_b2(2), 1/)), "saving UR block:"//trim(var%name) )
                    !         endif
                    !     elseif (this%creating) then
                    !         call check_ncdf( nf90_put_var(this%active_nc_id, var%file_var_id,  var%data_2dd(v_i_s:v_i_e,v_j_s:v_j_e), &
                    !                     start=(/ this%start_3d(1), this%start_3d(2) /), &
                    !                     count=(/ cnt_3d(1), cnt_3d(2) /)), "saving:"//trim(var%name) )
                    !         if (this%blocked_UR .or. this%blocked_LL) then
                    !             call check_ncdf( nf90_put_var(this%active_nc_id, var%file_var_id,  var%data_2dd(v_i_s_b:v_i_e_b,v_j_s_b:v_j_e_b), &
                    !                     start=(/ this%start_3d_b(1), this%start_3d_b(2) /), &
                    !                     count=(/ cnt_3d_b(1), cnt_3d_b(2) /)), "saving LL block:"//trim(var%name) )
                    !             call check_ncdf( nf90_put_var(this%active_nc_id, var%file_var_id,  var%data_2dd(v_i_s_b2:v_i_e_b2,v_j_s_b2:v_j_e_b2), &
                    !                     start=(/ this%start_3d_b2(1), this%start_3d_b2(2) /), &
                    !                     count=(/ cnt_3d_b2(1), cnt_3d_b2(2) /)), "saving UR block:"//trim(var%name) )
                    !         endif
                    !     endif
                    endif
                endif
            end associate
        end do
        
        call check_ncdf( nf90_var_par_access(this%active_nc_id, this%time%file_var_id, nf90_collective))

        output_time = get_output_time(time, units=this%time_units, round_seconds=.False.)

        call check_ncdf( nf90_put_var(this%active_nc_id, this%time%file_var_id, dble(output_time%mjd()), [current_step] ),   &
                    "saving:"//trim(this%time%name) )

    end subroutine save_data

    subroutine setup_dims_for_var(this, var_meta)
        implicit none
        class(output_t),    intent(inout) :: this
        type(meta_data_t),   intent(inout) :: var_meta
        integer :: i, err, dim_len
        character(len=kMAX_DIM_LENGTH), allocatable :: dimensions(:)
        character(len=kMAX_DIM_LENGTH) :: tmp_dim

        if (allocated(var_meta%dim_ids)) deallocate(var_meta%dim_ids)

        allocate(var_meta%dim_ids(size(var_meta%dimensions)))

        dimensions = var_meta%dimensions
        if (var_meta%three_d .or. var_meta%four_d) then
            tmp_dim = dimensions(2)
            dimensions(2) = dimensions(3)
            dimensions(3) = tmp_dim
        endif
        do i = 1, size(var_meta%dim_ids)

            ! Try to find the dimension ID if it exists already.
            err = nf90_inq_dimid(this%active_nc_id, trim(dimensions(i)), var_meta%dim_ids(i))

            ! probably the dimension doesn't exist in the file, so we will create it.
            if (err/=NF90_NOERR) then
                ! assume that the last dimension should be the unlimited dimension (generally a good idea...)
                if (var_meta%unlimited_dim .and. (i==size(var_meta%dim_ids)) .and. trim(dimensions(i))=="time") then
                    call check_ncdf( nf90_def_dim(this%active_nc_id, trim(dimensions(i)), NF90_UNLIMITED, &
                                var_meta%dim_ids(i) ), "def_dim"//dimensions(i) )
                else
                    if (i == 1) dim_len = this%file_dim_len(1)+var_meta%xstag
                    if (i == 2) dim_len = this%file_dim_len(2)+var_meta%ystag
                    if (i == 3) dim_len = var_meta%dim_len(2)

                    call check_ncdf( nf90_def_dim(this%active_nc_id, dimensions(i), dim_len,       &
                                var_meta%dim_ids(i) ), "def_dim"//dimensions(i) )
                endif
            endif
        end do

    end subroutine setup_dims_for_var

    subroutine setup_variable(this, var_meta)
        implicit none
        class(output_t),   intent(inout) :: this
        type(meta_data_t),  intent(inout) :: var_meta
        integer :: i, n, err

        err = nf90_inq_varid(this%active_nc_id, var_meta%name, var_meta%file_var_id)

        ! if the variable was not found in the netcdf file then we will define it.
        if (err /= NF90_NOERR) then
        
            if (var_meta%dtype == kREAL) then
                call check_ncdf( nf90_def_var(this%active_nc_id, var_meta%name, NF90_REAL, var_meta%dim_ids, var_meta%file_var_id), &
                        "Defining variable:"//trim(var_meta%name) )
            elseif (var_meta%dtype == kINTEGER) then
                call check_ncdf( nf90_def_var(this%active_nc_id, var_meta%name, NF90_INT, var_meta%dim_ids, var_meta%file_var_id), &
                        "Defining variable:"//trim(var_meta%name) )
            elseif (var_meta%dtype == kDOUBLE) then
                call check_ncdf( nf90_def_var(this%active_nc_id, var_meta%name, NF90_DOUBLE, var_meta%dim_ids, var_meta%file_var_id), &
                        "Defining variable:"//trim(var_meta%name) )
            endif

            ! setup attributes
            do i=1,size(var_meta%attributes)
                call check_ncdf( nf90_put_att(this%active_nc_id,                &
                                         var_meta%file_var_id,                    &
                                         trim(var_meta%attributes(i)%name),        &
                                         trim(var_meta%attributes(i)%value)),      &
                            "saving attribute "//trim(var_meta%attributes(i)%name)//" for "//trim(var_meta%name))
            enddo
        endif
        

    end subroutine setup_variable

    module subroutine increase_var_capacity(this)
        implicit none
        class(output_t),   intent(inout)  :: this
        type(variable_t),  allocatable :: new_variables(:)
        type(meta_data_t), allocatable :: new_meta_data(:)

        ! assert allocated(this%variables)
        allocate(new_variables, source=this%variables)
        allocate(new_meta_data, source=this%var_meta)
        ! new_variables = this%variables

        deallocate(this%variables)
        deallocate(this%var_meta)

        allocate(this%variables(size(new_variables)*2))
        this%variables(:size(new_variables)) = new_variables

        allocate(this%var_meta(size(new_meta_data)*2))
        this%var_meta(:size(new_meta_data)) = new_meta_data

        deallocate(new_variables)
        deallocate(new_meta_data)

    end subroutine

    module subroutine close_output_files(this)
        implicit none
        class(output_t),   intent(inout)  :: this

        if (this%out_ncfile_id > 0) then
            call check_ncdf(nf90_close(this%out_ncfile_id), "Closing output file ")
            this%out_ncfile_id = -1
        endif
        if (this%rst_ncfile_id > 0) then
            call check_ncdf(nf90_close(this%rst_ncfile_id), "Closing restart file ")
            this%rst_ncfile_id = -1
        endif

    end subroutine close_output_files

end submodule
