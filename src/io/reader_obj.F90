
!>----------------------------------------------------------
!!  Define the interface for the output object
!!
!!  Output objects store all of the data and references to data necessary to write
!!  an output file.  This includes primarily internal netcdf related IDs.
!!  Output objects also store an array of variables to output.
!!  These variables maintain pointers to the data to be output as well as
!!  Metadata (e.g. dimension names, units, other attributes)
!!
!!  @author
!!  Dylan Reynolds (dylan.reynolds@slf.ch)
!!  Ethan Gutmann (gutmann@ucar.edu)
!! 
!!----------------------------------------------------------
submodule(reader_interface) reader_implementation
  use debug_module,       only : check_ncdf
  use timer_interface,    only : timer_t
  use string,             only : as_string
  use time_io,            only : read_times, find_timestep_in_filelist, var_has_time_dim
  use array_utilities,    only : interpolate_in_z
  use io_routines,        only : io_read, io_getdims, io_var_reversed, can_file_parallel
  implicit none

contains


    module subroutine init_reader(this, its, ite, kts, kte, jts, jte, options)
        implicit none
        class(reader_t), intent(inout) :: this
        integer, intent(in) :: its, ite, kts, kte, jts, jte
        type(options_t), intent(in) :: options

        
        type(dim_arrays_type), allocatable           :: var_dimensions(:)
        character(len=kMAX_NAME_LENGTH), allocatable :: vars_to_read(:)
        type(meta_data_t) :: var_meta
        integer :: i, nx, ny, nz, dims(3)
        real, dimension(:,:,:), allocatable :: tmp_read_var
        

        this%file_list = options%forcing%boundary_files
        this%time_var  = options%forcing%time_var
        this%lat_var  = options%forcing%latvar

        this%model_end_time = options%general%end_time
        this%input_dt   = options%forcing%input_dt

        this%ncfile_id = -1

        if (options%restart%restart) then
            call set_curfile_curstep(this, options%restart%restart_time, this%file_list, this%time_var, options%restart%restart)
        else
            call set_curfile_curstep(this, options%general%start_time, this%file_list, this%time_var, options%restart%restart)
        endif
        !Now setup dimensions of reader_object so we know what part of input file that we should read
        this%its = its; this%ite = ite
        this%kts = kts; this%kte = kte
        this%jts = jts; this%jte = jte

        !get global input file dimensions
        call io_read(this%file_list(this%curfile), options%forcing%qvvar,   tmp_read_var,   this%curstep)

        this%ids = 1; this%ide = ubound(tmp_read_var, 1)
        this%jds = 1; this%jde = ubound(tmp_read_var, 2)

        ! the parameters option type can't contain allocatable arrays because it is a coarray
        ! so we need to allocate the vars_to_read and var_dimensions outside of the options type
        call setup_variable_lists(options%forcing%vars_to_read, options%forcing%dim_list, vars_to_read, var_dimensions)

        this%n_vars = size(vars_to_read)
        allocate(this%var_meta(this%n_vars))

        do i=1, this%n_vars
            nx = (this%ite-this%its+1)
            nz = (this%kte-this%kts+1)
            ny = (this%jte-this%jts+1)

            if (var_dimensions(i)%dims(1) == this%ide+1) nx = nx+1
            if (var_dimensions(i)%dims(2) == this%jde+1) ny = ny+1

            var_meta%name = trim(vars_to_read(i))
            if (allocated(var_meta%dim_len)) deallocate(var_meta%dim_len)

            if (var_dimensions(i)%num_dims==2) then
                allocate(var_meta%dim_len(2))
                var_meta%dim_len(1) = nx
                var_meta%dim_len(2) = ny
                var_meta%three_d = .False.
                var_meta%two_d = .True.
            elseif (var_dimensions(i)%num_dims==3) then
                allocate(var_meta%dim_len(3))
                var_meta%dim_len(1) = nx
                var_meta%dim_len(2) = nz
                var_meta%dim_len(3) = ny
                var_meta%three_d = .True.
                var_meta%two_d = .False.
            endif
            this%var_meta(i) = var_meta
        end do
        
    end subroutine init_reader

    !>------------------------------------------------------------
    !! Reads a new set of forcing data for the next time step
    !!
    !!------------------------------------------------------------
    module subroutine read_next_step(this, buffer, par_comms)
        class(reader_t), intent(inout) :: this
        real, allocatable, intent(inout) :: buffer(:,:,:,:)
        type(MPI_Comm), intent(in)              :: par_comms
        type(MPI_Info) :: par_comm_info

        real, allocatable :: data4d(:,:,:,:), data3d(:,:,:), data2d(:,:), data1d(:)
        type(meta_data_t)  :: var
        integer :: nx, ny, nz, data_nz, err, varid, n, ndims, i, k, j
        integer, dimension(4) :: start_3d_t, cnt_3d_t
        integer, dimension(3) :: start_2d_t, cnt_2d_t, start_3d, cnt_3d
        integer, dimension(2) :: start_2d, cnt_2d, start_1d_t, cnt_1d_t
        integer, dimension(1) :: start_1d, cnt_1d
        integer :: var_id, y_s, mid_y_local, mid_y_global
        integer, allocatable :: dims(:)
        character(len=kMAX_DIM_LENGTH) :: time_dim_name, attr_val
        logical :: has_time, flip_y, is_reversed, is_file_parallel

        if (allocated(buffer)) deallocate(buffer)
        allocate(buffer(this%n_vars,this%its:this%ite+1,this%kts:this%kte,this%jts:this%jte+1))

        !See if we must open the file

        is_file_parallel = .False.!can_file_parallel(this%file_list(this%curfile))
        if (this%ncfile_id < 0) then
            if (.not.(is_file_parallel)) then
                call check_ncdf( nf90_open(this%file_list(this%curfile), IOR(NF90_NOWRITE,NF90_NETCDF4), this%ncfile_id), " Opening file "//trim(this%file_list(this%curfile)))
            else
                par_comm_info = MPI_INFO_NULL
                call MPI_Comm_get_info(par_comms, par_comm_info)
                call check_ncdf( nf90_open(this%file_list(this%curfile), IOR(NF90_NOWRITE,NF90_NETCDF4), this%ncfile_id, &
                        comm = par_comms%MPI_VAL, info = par_comm_info%MPI_VAL), " Opening file "//trim(this%file_list(this%curfile)))
            endif
        endif
        
        flip_y = io_var_reversed(this%file_list(this%curfile), this%lat_var)

        !if flip_y, then the y index for starting needs to be flipped about the midpoint of the y dimension
        if (flip_y) then
            !first get midpoint of this y read 
            ny = this%jte - this%jts + 1
            mid_y_local = ny/2 + this%jts - 1
            y_s = this%jde - mid_y_local - ny/2
            ! mid_y_global = (this%jde - this%jds+1)/2
            ! if (mid_y_local <= mid_y_global) y_s = mid_y_global + (mid_y_global-mid_y_local) - ny/2
            ! if (mid_y_local > mid_y_global) y_s = mid_y_global - (mid_y_local-mid_y_global) - ny/2
        else
            y_s = this%jts
        endif
        ! setup start/count arrays accordingly
        start_3d_t = (/ this%its,y_s,this%kts,this%curstep /)
        start_2d_t = (/ this%its,y_s,this%curstep /)
        start_1d_t = (/ this%kts,this%curstep /)

        start_3d = (/ this%its,y_s,this%kts /)
        start_2d = (/ this%its,y_s /)
        start_1d = (/ this%kts /)

        cnt_3d_t = (/ (this%ite-this%its+1),(this%jte-this%jts+1),(this%kte-this%kts+1),1 /)
        cnt_2d_t = (/ (this%ite-this%its+1),(this%jte-this%jts+1),1 /)
        cnt_1d_t = (/ (this%kte-this%kts+1),1 /)

        cnt_3d = (/ (this%ite-this%its+1),(this%jte-this%jts+1),(this%kte-this%kts+1) /)
        cnt_2d = (/ (this%ite-this%its+1),(this%jte-this%jts+1) /)
        cnt_1d = (/ (this%kte-this%kts+1) /)

        ! determine if we need to flip the y dimension. We will determine this by looking for the
        ! latitude variable and checking for the attributes stored_direction = "decreasing" or "down" or
        ! positive = "down"


        ! loop through the list of variables that need to be read in
        do n = 1, size(this%var_meta)
            ! get the next variable in the structure
            var = this%var_meta(n)
            if (var%file_var_id < 0) call check_ncdf( nf90_inq_varid(this%ncfile_id, var%name, var%file_var_id), " Getting var ID for "//trim(var%name))
            if (is_file_parallel) call check_ncdf( nf90_var_par_access(this%ncfile_id, var%file_var_id, nf90_collective))

            !get number of dimensions
            call check_ncdf( nf90_inquire_variable(this%ncfile_id, var%file_var_id, ndims = ndims), " Getting dim length for "//trim(var%name))

            !see if one of the dimensions is unlimited (time)
            has_time = var_has_time_dim(this%file_list(this%curfile), var%name, this%time_var)

            if (var%three_d) then
                nx = var%dim_len(1)
                ny = var%dim_len(3)
                nz = var%dim_len(2)

                if (ndims == 4) then
                    if (allocated(data4d)) deallocate(data4d)
                    if (has_time) then

                        !get dimension length of 3rd dimension to see if it is staggered in vertical
                        call io_getdims(this%file_list(this%curfile), var%name, dims)
                        data_nz = dims(3)

                        cnt_3d_t = (/ nx, ny, data_nz,1 /)
                        allocate(data4d(nx,ny,data_nz,1))

                        call check_ncdf( nf90_get_var(this%ncfile_id, var%file_var_id, data4d, start=start_3d_t, count=cnt_3d_t), " Getting 3D var with time dim: "//trim(var%name))

                        !check if variable is staggered in vertical
                        ! if so, we need to interpolate it to the mass levels
                        if (data_nz == nz+1) then
                            if (allocated(data3d)) deallocate(data3d)
                            allocate(data3d(nx,ny,data_nz))
                            data3d = data4d(:,:,:,1)
                            call interpolate_in_z(data3d,zdim=3)
                            deallocate(data4d)
                            allocate(data4d(nx,ny,nz,1))
                            data4d(:,:,:,1) = data3d(:,:,:)
                            deallocate(data3d)
                        endif

                        buffer(n,this%its:this%its+nx-1,this%kts:this%kts+nz-1,this%jts:this%jts+ny-1) = reshape(data4d(:,:,:,1), shape=[nx,nz,ny], order=[1,3,2])
                    else
                        write(*,*) "Error: ", trim(var%name), " is spatially 3D, but forcing data is spatially 4D"
                        stop "Internal variable is spatially 3D, but forcing data is spatially 4D"
                    endif
                elseif (ndims==3) then
                    if (allocated(data3d)) deallocate(data3d)
                    if (has_time) then

                        cnt_2d_t = (/ nx, ny, 1 /)

                        !if 2D spatial data has been provided for a 3D variable, assume that we are given a 2D plane and should replicate it in the vertical
                        allocate(data3d(nx,ny,1))
                        call check_ncdf( nf90_get_var(this%ncfile_id, var%file_var_id, data3d, start=start_2d_t, count=cnt_2d_t), " Getting 2D var with time dim for 3D var: "//trim(var%name))
                        ! Broadcast the 2D data across all vertical levels
                        do k = this%kts, this%kte
                            buffer(n,this%its:this%its+nx-1,k,this%jts:this%jts+ny-1) = data3d(:,:,1)
                        enddo
                    else

                        !get dimension length of 3rd dimension to see if it is staggered in vertical
                        call io_getdims(this%file_list(this%curfile), var%name, dims)
                        data_nz = dims(3)

                        cnt_3d = (/ nx, ny, data_nz /)
                        allocate(data3d(nx,ny,data_nz))
                        
                        call check_ncdf( nf90_get_var(this%ncfile_id, var%file_var_id, data3d, start=start_3d, count=cnt_3d), " Getting 3D var: "//trim(var%name))

                        !check if variable is staggered in vertical
                        ! if so, we need to interpolate it to the mass levels
                        if (data_nz == nz+1) then
                            call interpolate_in_z(data3d,zdim=3)
                        endif

                        buffer(n,this%its:this%its+nx-1,this%kts:this%kts+nz-1,this%jts:this%jts+ny-1) = reshape(data3d(:,:,:), shape=[nx,nz,ny], order=[1,3,2])
                    endif
                elseif (ndims==2) then
                    if (allocated(data2d)) deallocate(data2d)
                    if (has_time) then

                        !get dimension length of 1st dimension to see if it is staggered in vertical
                        call io_getdims(this%file_list(this%curfile), var%name, dims)
                        data_nz = dims(1)

                        cnt_1d_t = (/ data_nz, 1 /)
                        allocate(data2d(data_nz,1))

                        !if 1D spatial data has been provided for a 3D variable, assume that we are given a 1D line and should fill each z level with it
                        call check_ncdf( nf90_get_var(this%ncfile_id, var%file_var_id, data2d, start=start_1d_t, count=cnt_1d_t), " Getting 1D var with time dim for 3D var: "//trim(var%name))

                        !check if variable is staggered in vertical
                        ! if so, we need to interpolate it to the mass levels
                        if (data_nz == nz+1) then
                            if (allocated(data1d)) deallocate(data1d)
                            allocate(data1d(nz))
                            data1d = (data2d(2:nz+1,1)+data2d(1:nz,1))/2.0
                            deallocate(data2d)
                            allocate(data2d(nz,1))
                            data2d(:,1) = data1d
                        endif

                        !reverse the order of the 1D data if needed
                        is_reversed = io_var_reversed(this%file_list(this%curfile), var%name)
                        if (is_reversed) data2d(:,1) = data2d(size(data2d,1):1:-1,1)

                        ! Broadcast the 1D data across all horizontal levels
                        do j = this%jts, this%jts+ny-1
                            do i = this%its, this%its+nx-1
                                buffer(n,i,this%kts:this%kte,j) = data2d(:,1)
                            enddo
                        enddo
                    else

                        cnt_2d = (/ nx, ny /)

                        !if 2D spatial data has been provided for a 3D variable, assume that we are given a 2D plane and should replicate it in the vertical
                        allocate(data2d(nx,ny))
                        call check_ncdf( nf90_get_var(this%ncfile_id, var%file_var_id, data2d, start=start_2d, count=cnt_2d), " Getting 2D var for 3D var: "//trim(var%name))
                        ! Broadcast the 2D data across all vertical levels
                        do k = this%kts, this%kte
                            buffer(n,this%its:this%its+nx-1,k,this%jts:this%jts+ny-1) = data2d(:,:)
                        enddo
                    endif
                elseif (ndims==1) then
                    if (allocated(data1d)) deallocate(data1d)
                    if (.not.(has_time)) then

                        !get dimension length of 1st dimension to see if it is staggered in vertical
                        call io_getdims(this%file_list(this%curfile), var%name, dims)
                        data_nz = dims(1)

                        cnt_1d = (/ data_nz /)
                        allocate(data1d(data_nz))

                        !if 1D spatial data has been provided for a 3D variable, assume that we are given a 1D line and should fill each z level with it
                        call check_ncdf( nf90_get_var(this%ncfile_id, var%file_var_id, data1d, start=start_1d, count=cnt_1d), " Getting 1D var for 3D var: "//trim(var%name))
                        
                        !check if variable is staggered in vertical
                        ! if so, we need to interpolate it to the mass levels
                        if (data_nz == nz+1) then
                            if (allocated(data2d)) deallocate(data2d)
                            allocate(data2d(nz,1))
                            data2d(:,1) = (data1d(2:nz+1)+data1d(1:nz))/2.0
                            deallocate(data1d)
                            allocate(data1d(nz))
                            data1d = data2d(:,1)
                        endif

                        !reverse the order of the 1D data if needed
                        is_reversed = io_var_reversed(this%file_list(this%curfile), var%name)
                        if (is_reversed) data1d = data1d(size(data1d):1:-1)

                        ! Broadcast the 1D data across all horizontal levels
                        do j = this%jts, this%jts+ny-1
                            do i = this%its, this%its+nx-1
                                buffer(n,i,this%kts:this%kte,j) = data1d(this%kts:this%kte)
                            enddo
                        enddo
                    endif
                endif

                if (flip_y) then
                    !invert buffer about the y axis using a temporary array to avoid overlap
                    if (allocated(data3d)) deallocate(data3d)
                    allocate(data3d(nx,nz,ny))
                    data3d = buffer(n,this%its:this%its+nx-1,this%kts:this%kte,this%jts:this%jts+ny-1)
                    buffer(n,this%its:this%its+nx-1,this%kts:this%kte,this%jts:this%jts+ny-1) = data3d(:,:,ny:1:-1)
                    deallocate(data3d)
                endif
            else if (var%two_d) then
                nx = var%dim_len(1)
                ny = var%dim_len(2)

                cnt_2d_t = (/ nx, ny, 1 /)
                cnt_2d = (/ nx, ny /)

                if (ndims == 3) then
                    call check_ncdf( nf90_get_var(this%ncfile_id, var%file_var_id, buffer(n,this%its:this%its+nx-1,1,this%jts:this%jts+ny-1), start=start_2d_t, count=cnt_2d_t), " Getting 2D var with time dim: "//trim(var%name))
                elseif (ndims==2) then
                    call check_ncdf( nf90_get_var(this%ncfile_id, var%file_var_id, buffer(n,this%its:this%its+nx-1,1,this%jts:this%jts+ny-1), start=start_2d, count=cnt_2d), " Getting 2D var: "//trim(var%name))
                endif
                if (flip_y) then
                    !invert buffer about the y axis using a temporary array to avoid overlap
                    if (allocated(data2d)) deallocate(data2d)
                    allocate(data2d(nx,ny))
                    data2d = buffer(n,this%its:this%its+nx-1,1,this%jts:this%jts+ny-1)
                    buffer(n,this%its:this%its+nx-1,1,this%jts:this%jts+ny-1) = data2d(:,ny:1:-1)
                    deallocate(data2d)
                endif
            endif

        end do
        
        call update_forcing_step(this)

    end subroutine
    
    
    !>------------------------------------------------------------
    !! Update the curstep and curfile (increments curstep and curfile if necessary)
    !!
    !!------------------------------------------------------------
    subroutine update_forcing_step(this)
        implicit none
        type(reader_t),   intent(inout) :: this

        integer :: steps_in_file
        type(Time_type), allocatable :: times_in_file(:)
        type(Time_type) :: time_tmp


        this%curstep = this%curstep + 1 ! this may be all we have to do most of the time
        ! check that we haven't stepped passed the end of the current file
        steps_in_file = get_n_timesteps(this, this%time_var, 0)

        call this%input_time%set(this%input_time%mjd() + this%input_dt%days())

        if (steps_in_file < this%curstep) then
            ! close current file
            call check_ncdf(nf90_close(this%ncfile_id), "Closing file "//trim(this%file_list(this%curfile)))
            this%ncfile_id = -1
            
            ! if we have, use the next file
            this%curfile = this%curfile + 1
            ! and the first timestep in the next file
            this%curstep = 1



            ! if we have run out of input files, stop with an error message
            if (this%curfile > size(this%file_list)) then
                this%eof = .True.
            endif

        endif

        ! Check if the next file to read is beyond the model end time.
        if (.not.(this%eof)) then
            !Get time step of the next input step
            call read_times(this%file_list(this%curfile), this%time_var, times_in_file)
            
            call time_tmp%set(times_in_file(this%curstep)%mjd() - this%input_dt%days())

            if (time_tmp >= this%model_end_time) then
                this%eof = .True.
            !Check if the next time step in the file is equal to the input time 
            else if (.not.(times_in_file(this%curstep) == this%input_time)) then
                ! warn the user and exit
                write(*,*) "Warning: The next time step in the file is not equal to the input time.  The next time step in the file is: ", as_string(times_in_file(this%curstep))
                write(*,*) "The input time is: ", as_string(this%input_time)
                write(*,*) "The current file is: ", this%file_list(this%curfile)
                write(*,*) "The current step is: ", this%curstep
                stop
            endif
        endif

    end subroutine

    !>------------------------------------------------------------
    !! Figure out how many time steps are in a file based on a specified variable
    !!
    !! By default assumes that the variable has three dimensions.  If not, var_space_dims must be set
    !!------------------------------------------------------------
    function get_n_timesteps(this, varname, var_space_dims) result(steps_in_file)
        implicit none
        type(reader_t),   intent(in) :: this
        character(len=*), intent(in) :: varname
        integer,          intent(in), optional :: var_space_dims
        integer :: steps_in_file, i

        integer :: space_dims, varid, ndims, numDims, dimlen, dims(100), dimIds(100)

        space_dims=3
        if (present(var_space_dims)) space_dims = var_space_dims

        ! Get the varid of the variable, based on its name.
        call check_ncdf(nf90_inq_varid(this%ncfile_id, varname, varid),varname)
        ! find the number of dimensions
        call check_ncdf(nf90_inquire_variable(this%ncfile_id, varid, ndims = numDims),varname)
        ! find the dimension IDs
        call check_ncdf(nf90_inquire_variable(this%ncfile_id, varid, dimids = dimIds(:numDims)),varname)
        dims(1)=numDims
        ! finally, find the length of each dimension
        do i=1,numDims
            call check_ncdf(nf90_inquire_dimension(this%ncfile_id, dimIds(i), len = dimlen))
            dims(i+1)=dimlen
        end do

        if (dims(1) == space_dims) then
            steps_in_file = 1
        else
            steps_in_file = dims(dims(1)+1)
        endif

    end function


    !>------------------------------------------------------------
    !! Set the boundary data structure to the correct time step / file in the list of files
    !!
    !! Reads the time_var from each file successively until it finds a timestep that matches time
    !!------------------------------------------------------------
    subroutine set_curfile_curstep(this, time, file_list, time_var, restart)
        implicit none
        class(reader_t), intent(inout) :: this
        type(Time_type),    intent(in) :: time
        character(len=*),   intent(in) :: file_list(:)
        character(len=*),   intent(in) :: time_var
        logical,            intent(in) :: restart

        type(Time_type), allocatable :: times_in_file(:)
        character(len=kMAX_FILE_LENGTH) :: filename
        integer          :: error, n

        ! if this is a restart run, it is acceptable to find a non-exact first file time, 
        ! in which case we take the forward time (assuming restart was written between input steps)
        if (restart) then
            this%curstep = find_timestep_in_filelist(file_list, time_var, time, filename, forward=.False., error=error)
        else
            this%curstep = find_timestep_in_filelist(file_list, time_var, time, filename, error=error)
        endif

        if (error==1) then
            stop "Ran out of files to process while searching for matching time variable!"
        endif
        
        do n = 1,size(file_list)
            if (trim(file_list(n))==trim(filename)) this%curfile = n
        enddo

        call read_times(file_list(this%curfile), time_var, times_in_file)
        this%input_time = times_in_file(this%curstep)

        this%eof = .False.
    end subroutine

    !>------------------------------------------------------------
    !! Setup the vars_to_read and var_dimensions arrays given a master set of variables
    !!
    !! Count the number of variables specified, then allocate and store those variables in a list just their size.
    !! The master list will have all variables, but not all will be set
    !!------------------------------------------------------------
    subroutine setup_variable_lists(master_var_list, master_dim_list, vars_to_read, var_dimensions)
        implicit none
        character(len=kMAX_NAME_LENGTH), intent(in)                 :: master_var_list(:)
        type(dim_arrays_type),           intent(in)                 :: master_dim_list(:)
        character(len=kMAX_NAME_LENGTH), intent(inout), allocatable :: vars_to_read(:)
        type(dim_arrays_type),           intent(inout), allocatable :: var_dimensions(:)

        integer :: n_valid_vars
        integer :: i, curvar, err

        n_valid_vars = 0
        do i=1, size(master_var_list)
            if (trim(master_var_list(i)) /= '') then
                !check if "_computed" is in the variable name, if so skip it
                if (index(trim(master_var_list(i)), '_computed') > 0) cycle
                n_valid_vars = n_valid_vars + 1
            endif
        enddo

        allocate(vars_to_read(  n_valid_vars), stat=err)
        if (err /= 0) stop "vars_to_read: Allocation request denied"

        allocate(var_dimensions(  n_valid_vars), stat=err)
        if (err /= 0) stop "var_dimensions: Allocation request denied"

        curvar = 1
        do i=1, size(master_var_list)
            if (trim(master_var_list(i)) /= '') then
                !check if "_computed" is in the variable name, if so skip it
                if (index(trim(master_var_list(i)), '_computed') > 0) cycle
                vars_to_read(curvar) = master_var_list(i)
                var_dimensions(curvar)%num_dims = master_dim_list(i)%num_dims
                var_dimensions(curvar)%dims(1:master_dim_list(i)%num_dims) = master_dim_list(i)%dims(1:master_dim_list(i)%num_dims)
                ! if (STD_OUT_PE) print *, "in variable list: ", vars_to_read(curvar)
                curvar = curvar + 1
            endif
        enddo
    end subroutine

    module subroutine close_file(this)
        implicit none
        class(reader_t),   intent(inout)  :: this

        if (this%ncfile_id > 0) then
            call check_ncdf(nf90_close(this%ncfile_id), "Closing file ")
            this%ncfile_id = -1
        endif

    end subroutine

    
end submodule
