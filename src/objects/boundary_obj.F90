!>------------------------------------------------------------
!!  Handles reading boundary conditions from the forcing file(s)
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!------------------------------------------------------------
submodule(boundary_interface) boundary_implementation
    use array_utilities,        only : interpolate_in_z
    use io_routines,            only : io_getdims, io_read, io_maxDims, io_variable_is_present, io_write
    use time_io,                only : read_times, find_timestep_in_filelist
    use string,                 only : str, as_string
    use mod_atm_utilities,      only : rh_to_mr, relative_humidity, compute_3d_p, compute_3d_z, exner_function
    use geo,                    only : standardize_coordinates
    use vertical_interpolation, only : vLUT, vinterp
    use timer_interface,    only : timer_t
    use debug_module,           only : check_ncdf
    use mod_wrf_constants,      only : gravity
    use variable_interface,     only : variable_t
    use icar_constants,         only : STD_OUT_PE, kMAX_FILE_LENGTH, kVARS
    implicit none
contains

    !>------------------------------------------------------------
    !! Set default component values
    !! Reads initial conditions from the forcing file into image 1
    !!
    !! Distributes initial conditions to all other images
    !!
    !!------------------------------------------------------------
    module subroutine init_boundary(this, options, domain_lat, domain_lon, parent_options)
        class(boundary_t),    intent(inout) :: this
        type(options_t),      intent(inout) :: options
        real, dimension(:,:), intent(in)    :: domain_lat
        real, dimension(:,:), intent(in)    :: domain_lon
        type(options_t), optional, intent(in)    :: parent_options

        type(Time_type) :: strt_time
        character(len=kMAX_NAME_LENGTH), allocatable :: vars_to_read(:)
        integer,                         allocatable :: var_dimensions(:), var_indx(:)


        ! the parameters option type can't contain allocatable arrays because it is a coarray
        ! so we need to allocate the vars_to_read and var_dimensions outside of the options type
        call setup_variable_lists(options%forcing%vars_to_read, options%forcing%dim_list, options, vars_to_read, var_dimensions, var_indx)

        strt_time = options%general%start_time
        if (options%restart%restart) strt_time = options%restart%restart_time

        ! Read through forcing variable names stored in "options"
        ! needs to read each one to find the grid information for it
        ! then create grid and initialize a variable...
        ! also need to explicitly save lat and lon data
        ! if (STD_OUT_PE) then
        if (present(parent_options)) then
            call this%init_local_asnest(vars_to_read, var_dimensions, var_indx,   &
                                    domain_lat,        &
                                    domain_lon,       &
                                    parent_options)
        else
            call this%init_local(options,                           &
                                    options%forcing%boundary_files, &
                                    vars_to_read, var_dimensions, var_indx,  &
                                    strt_time,                      &
                                    options%forcing%latvar,         &
                                    options%forcing%lonvar,         &
                                    options%forcing%zvar,           &
                                    options%forcing%time_var,       &
                                    options%forcing%pvar,           &
                                    domain_lat,        &
                                    domain_lon)
        endif
        ! endif
        ! call this%distribute_initial_conditions()


        call setup_boundary_geo(this, options%domain%longitude_system)

    end subroutine init_boundary


    !>------------------------------------------------------------
    !! Set default component values
    !! Reads initial conditions from the forcing file
    !!
    !!------------------------------------------------------------
    module subroutine init_local(this, options, file_list, var_list, dim_list, var_indx, start_time, &
                                 lat_var, lon_var, z_var, time_var, p_var, domain_lat, domain_lon)
        class(boundary_t),               intent(inout)  :: this
        type(options_t),                 intent(inout)  :: options
        character(len=kMAX_FILE_LENGTH), intent(in)     :: file_list(:)
        character(len=kMAX_NAME_LENGTH), intent(in)     :: var_list(:)
        integer,                         intent(in)     :: dim_list(:), var_indx(:)
        type(Time_type),                 intent(in)     :: start_time
        character(len=kMAX_NAME_LENGTH), intent(in)     :: lat_var
        character(len=kMAX_NAME_LENGTH), intent(in)     :: lon_var
        character(len=kMAX_NAME_LENGTH), intent(in)     :: z_var
        character(len=kMAX_NAME_LENGTH), intent(in)     :: time_var
        character(len=kMAX_NAME_LENGTH), intent(in)     :: p_var
        real, dimension(:,:),            intent(in)     :: domain_lat
        real, dimension(:,:),            intent(in)     :: domain_lon

        real, allocatable :: temp_z(:,:,:), temp_z_trans(:,:,:), temp_lat(:,:), temp_lon(:,:), lat_1d(:), lon_1d(:)
        integer, allocatable :: lat_dims(:), lon_dims(:)
        real :: neg_z
        integer :: i, nx, ny, nz, PE_RANK_GLOBAL, x_len, y_len

        ! figure out while file and timestep contains the requested start_time
        call set_firstfile_firststep(this, start_time, file_list, time_var)

        !See if lat and lon are 1D or 2D
        call io_getdims(this%firstfile, lat_var, lat_dims)
        call io_getdims(this%firstfile, lon_var, lon_dims)

        if (size(lat_dims) == 1) then
            call io_read(this%firstfile, lat_var, lat_1d)

            !This will always be the case for 1D or 2D lon
            x_len = lon_dims(1)

            allocate(temp_lat(1:x_len,1:lat_dims(1)))
            do i=1,x_len
                temp_lat(i,:) = lat_1d
            end do
        elseif (size(lat_dims) == 2) then
            call io_read(this%firstfile, lat_var, temp_lat)
        else
            write(*,*) 'ERROR: lat dimension on forcing data is not 1D or 2D'
            stop
        endif

        if (size(lon_dims) == 1) then
            call io_read(this%firstfile, lon_var, lon_1d)

            if (size(lat_dims) == 1) then
                y_len = lat_dims(1)
            else
                y_len = lat_dims(2)
            endif

            allocate(temp_lon(1:lon_dims(1),1:y_len))
            do i=1,y_len
                temp_lon(:,i) = lon_1d
            end do
        elseif (size(lon_dims) == 2) then
            call io_read(this%firstfile, lon_var, temp_lon)
        else
            write(*,*) 'ERROR: lon dimension on forcing data is not 1D or 2D'
            stop
        endif

        call MPI_Comm_Rank(MPI_COMM_WORLD,PE_RANK_GLOBAL)

        if (minval(domain_lat) < minval(temp_lat) .or. maxval(domain_lat) > maxval(temp_lat)) then
            write(*,*) 'ERROR: First domain not contained within forcing data'
            write(*,*) 'Lat min/max of domain on process ',PE_RANK_GLOBAL+1,': ',minval(domain_lat),' ',maxval(domain_lat)
            write(*,*) 'Lat min/max of forcing data:         ',minval(temp_lat),' ',maxval(temp_lat)
            stop
        endif
        if (minval(domain_lon) < minval(temp_lon) .or. maxval(domain_lon) > maxval(temp_lon)) then
            write(*,*) 'ERROR: First domain not contained within forcing data'
            write(*,*) 'Lon min/max of domain on process ',PE_RANK_GLOBAL+1,': ',minval(domain_lon),' ',maxval(domain_lon)
            write(*,*) 'Lon min/max of forcing data:         ',minval(temp_lon),' ',maxval(temp_lon)
            stop
        endif

        !Here we should begin the sub-setting:
        !domain object should be fed in as an argument
        !The lat/lon bounds of the domain object are used to find the appropriate indexes of the forcing data
        !These bounds are then extended by 1 in each direction to accomodate bilinear interpolation
        
        call set_boundary_image(this, temp_lat, temp_lon, domain_lat, domain_lon)

        !After finding forcing image indices, subset the lat/lon variables
        allocate(this%lat((this%ite-this%its+1),(this%jte-this%jts+1)))
        allocate(this%lon((this%ite-this%its+1),(this%jte-this%jts+1)))

        this%lat = temp_lat(this%its:this%ite,this%jts:this%jte)
        this%lon = temp_lon(this%its:this%ite,this%jts:this%jte)

        ! read in the height coordinate of the input data
        if (.not. options%forcing%compute_z) then
            ! call io_read(file_list(this%curfile), z_var,   temp_z,   this%curstep)
            call io_read(this%firstfile, z_var,   temp_z,   this%firststep)
            nx = size(temp_z,1)
            ny = size(temp_z,2)
            nz = size(temp_z,3)
            if (allocated(this%z)) deallocate(this%z)
            allocate(this%z((this%ite-this%its+1),nz,(this%jte-this%jts+1)))
            allocate(temp_z_trans(1:nx,1:nz,1:ny))
            
            temp_z_trans(1:nx,1:nz,1:ny) = reshape(temp_z, shape=[nx,nz,ny], order=[1,3,2])
            this%z = temp_z_trans(this%its:this%ite,1:nz,this%jts:this%jte)
            this%z_is_set = .True.

            if (options%forcing%z_is_geopotential) then
                this%z = this%z / gravity
                !neg_z = minval(temp_z)/ gravity
                ! if (neg_z < 0.0) this%z = this%z - neg_z
            endif

            if (options%forcing%z_is_on_interface) then
                call interpolate_in_z(this%z)
            endif

        else
            call io_read(this%firstfile, p_var,   temp_z,   this%firststep)
            nx = size(temp_z,1)
            ny = size(temp_z,2)
            nz = size(temp_z,3)

            if (allocated(this%z)) deallocate(this%z)
            allocate(this%z((this%ite-this%its+1),nz,(this%jte-this%jts+1)))

        endif
        
        this%kts = 1
        this%kte = nz

        if (this%ite < this%its) write(*,*) 'image: ',PE_RANK_GLOBAL+1,'  its: ',this%its,'  ite: ',this%ite,'  jts: ',this%jts,'  jte: ',this%jte
        if (this%kte < this%kts) write(*,*) 'image: ',PE_RANK_GLOBAL+1,'  its: ',this%its,'  ite: ',this%ite,'  jts: ',this%jts,'  jte: ',this%jte
        if (this%jte < this%jts) write(*,*) 'image: ',PE_RANK_GLOBAL+1,'  its: ',this%its,'  ite: ',this%ite,'  jts: ',this%jts,'  jte: ',this%jte

        ! call assert(size(var_list) == size(dim_list), "list of variable dimensions must match list of variables")
        do i=1, size(var_list)
            call add_var_to_dict(this, var_list(i), dim_list(i), var_indx(i), [(this%ite-this%its+1), nz, (this%jte-this%jts+1)])
        end do


    end subroutine

    !>------------------------------------------------------------
    !! Set default component values
    !! Reads initial conditions from the forcing file
    !!
    !!------------------------------------------------------------
    module subroutine init_local_asnest(this, var_list, dim_list, var_indx, domain_lat, domain_lon, parent_options)
        class(boundary_t),               intent(inout)  :: this
        character(len=kMAX_NAME_LENGTH), intent(in)     :: var_list(:)
        integer,                         intent(in)     :: dim_list(:), var_indx(:)
        real, dimension(:,:),            intent(in)     :: domain_lat
        real, dimension(:,:),            intent(in)     :: domain_lon
        type(options_t),                 intent(in)     :: parent_options
        
        integer :: i
        real, allocatable, dimension(:,:) :: parent_nest_lat, parent_nest_lon

        !Here we should begin the sub-setting:
        !domain object should be fed in as an argument
        !The lat/lon bounds of the domain object are used to find the appropriate indexes of the forcing data
        !These bounds are then extended by 1 in each direction to accomodate bilinear interpolation

        !  read in latitude and longitude coordinate data
        call io_read(parent_options%domain%init_conditions_file, parent_options%domain%lat_hi, parent_nest_lat)
        call io_read(parent_options%domain%init_conditions_file, parent_options%domain%lon_hi, parent_nest_lon)

        if (minval(domain_lat) < minval(parent_nest_lat) .or. maxval(domain_lat) > maxval(parent_nest_lat)) then
            write(*,*) 'ERROR: Nested domain not contained within parent domain: ',trim(parent_options%domain%init_conditions_file)
            write(*,*) 'Lat min/max of nested domain: ',minval(domain_lat),' ',maxval(domain_lat)
            write(*,*) 'Lat min/max of parent domain:               ',minval(parent_nest_lat),' ',maxval(parent_nest_lat)
            stop
        endif
        if (minval(domain_lon) < minval(parent_nest_lon) .or. maxval(domain_lon) > maxval(parent_nest_lon)) then
            write(*,*) 'ERROR: Nested domain not contained within parent domain: ',trim(parent_options%domain%init_conditions_file)
            write(*,*) 'Lon min/max of nested domain: ',minval(domain_lon),' ',maxval(domain_lon)
            write(*,*) 'Lon min/max of parent domain:               ',minval(parent_nest_lon),' ',maxval(parent_nest_lon)
            stop
        endif

        call set_boundary_image(this, parent_nest_lat, parent_nest_lon, domain_lat, domain_lon)

        !After finding forcing image indices, subset the lat/lon variables
        allocate(this%lat((this%ite-this%its+1),(this%jte-this%jts+1)))
        allocate(this%lon((this%ite-this%its+1),(this%jte-this%jts+1)))

        this%lat = parent_nest_lat(this%its:this%ite,this%jts:this%jte)
        this%lon = parent_nest_lon(this%its:this%ite,this%jts:this%jte)

        ! Get the height coordinate of the parent nest
        this%kts = 1
        this%kte = parent_options%domain%nz

        if (allocated(this%z)) deallocate(this%z)
        allocate(this%z((this%ite-this%its+1),this%kte,(this%jte-this%jts+1)))
        this%z = 0.0 !parent_nest_z(this%its:this%ite,1:this%kte,this%jts:this%jte)
        ! Set this flag here so that we will set z later in domain%setup_geo_interpolation to be the z data read from forcing data
        this%z_is_set = .False.

        ! if (this%ite < this%its) write(*,*) 'image: ',PE_RANK_GLOBAL+1,'  its: ',this%its,'  ite: ',this%ite,'  jts: ',this%jts,'  jte: ',this%jte
        ! if (this%kte < this%kts) write(*,*) 'image: ',PE_RANK_GLOBAL+1,'  its: ',this%its,'  ite: ',this%ite,'  jts: ',this%jts,'  jte: ',this%jte
        ! if (this%jte < this%jts) write(*,*) 'image: ',PE_RANK_GLOBAL+1,'  its: ',this%its,'  ite: ',this%ite,'  jts: ',this%jts,'  jte: ',this%jte

        ! call assert(size(var_list) == size(dim_list), "list of variable dimensions must match list of variables")
        do i=1, size(var_list)
            call add_var_to_dict(this, var_list(i), dim_list(i), var_indx(i), [(this%ite-this%its+1), this%kte, (this%jte-this%jts+1)])
        end do

    end subroutine

    !>------------------------------------------------------------
    !! Set the boundary data structure to the correct time step / file in the list of files
    !!
    !! Reads the time_var from each file successively until it finds a timestep that matches time
    !!------------------------------------------------------------
    subroutine set_firstfile_firststep(this, time, file_list, time_var)
        implicit none
        class(boundary_t),  intent(inout) :: this
        type(Time_type),    intent(in) :: time
        character(len=*),   intent(in) :: file_list(:)
        character(len=*),   intent(in) :: time_var

        character(len=kMAX_FILE_LENGTH) :: filename
        integer          :: error, n

        this%firststep = find_timestep_in_filelist(file_list, time_var, time, this%firstfile, forward=.False., error=error)
        
        if (error==1) then
            if (STD_OUT_PE) write(*,*) 'ERROR: Could not find the first time step in forcing files'
            if (STD_OUT_PE) write(*,*) 'ERROR: Time: ',as_string(time)
            if (STD_OUT_PE) write(*,*) 'ERROR: step returned was: ',this%firststep
            stop "Ran out of files to process while searching for matching time variable!"
        endif

    end subroutine

    !>------------------------------------------------------------
    !! Determine image location for boundary object within forcing data
    !!
    !!------------------------------------------------------------

    subroutine set_boundary_image(this, temp_lat, temp_lon, domain_lat, domain_lon)
        implicit none
        type(boundary_t), intent(inout)   :: this
        real, dimension(:,:), intent(in)  :: temp_lat, temp_lon, domain_lat, domain_lon
        
        real, allocatable, dimension(:,:) ::  LL_d, UR_d
        real :: LLlat, LLlon, URlat, URlon
        real, dimension(4) :: lat_corners, lon_corners
        integer, dimension(2) :: temp_inds
        integer :: nx, ny, d_ims, d_ime, d_jms, d_jme
        
        d_ims = lbound(domain_lat,1)
        d_ime = ubound(domain_lat,1)
        d_jms = lbound(domain_lat,2)
        d_jme = ubound(domain_lat,2)

        nx = size(temp_lat,1)
        ny = size(temp_lat,2)
        
        allocate(LL_d(nx,ny))
        allocate(UR_d(nx,ny))

        ! get lower left and upper right lat/lon pairs from domain
        lat_corners(1) = domain_lat(d_ims,d_jms)
        lat_corners(2) = domain_lat(d_ime,d_jms)
        lat_corners(3) = domain_lat(d_ims,d_jme)
        lat_corners(4) = domain_lat(d_ime,d_jme)
        
        lon_corners(1) = domain_lon(d_ims,d_jms)
        lon_corners(2) = domain_lon(d_ime,d_jms)
        lon_corners(3) = domain_lon(d_ims,d_jme)
        lon_corners(4) = domain_lon(d_ime,d_jme)
        
        LLlat = minval(lat_corners)
        LLlon = minval(lon_corners)
        URlat = maxval(lat_corners)
        URlon = maxval(lon_corners)

        ! calculate distance from LL/UR lat/lon for boundary lat/lons
        LL_d = ((temp_lat-LLlat)**2+(temp_lon-LLlon)**2)
        UR_d = ((temp_lat-URlat)**2+(temp_lon-URlon)**2)

        ! find minimum distances to determine boundary image indices
        temp_inds = minloc(LL_d)
        this%its = temp_inds(1); this%jts = temp_inds(2)
        temp_inds = minloc(UR_d)
        this%ite = temp_inds(1); this%jte = temp_inds(2)

        ! increase boundary image indices by 5 as buffer to allow for interpolation
        this%its = max(this%its - 8,1)
        this%ite = min(this%ite + 8,nx)
        this%jts = max(this%jts - 8,1)
        this%jte = min(this%jte + 8,ny)

        ! if (this%ite < this%its .or. this%jte < this%jts) write(*,*) 'image: ',PE_RANK_GLOBAL+1,'  its: ',this%its,'  ite: ',this%ite,'  jts: ',this%jts,'  jte: ',this%jte
        ! if (this%ite < this%its .or. this%jte < this%jts) write(*,*) 'image: ',PE_RANK_GLOBAL+1,'  d_ims: ',d_ims,'  d_ime: ',d_ime,'  d_jms: ',d_jms,'  d_jme: ',d_jme
        ! if (this%ite < this%its .or. this%jte < this%jts) write(*,*) 'image: ',PE_RANK_GLOBAL+1,'  LLlat: ',LLlat,'  LLlon: ',LLlon,'  URlat: ',URlat,'  URlon: ',URlon
        ! if (this%ite < this%its .or. this%jte < this%jts) write(*,*) 'image: ',PE_RANK_GLOBAL+1,'  min_loc: ',minloc(LL_d),'  max_loc: ',minloc(UR_d)
        ! if (this%ite < this%its .or. this%jte < this%jts) write(*,*) 'image: ',PE_RANK_GLOBAL+1,'  min_lat: ',minval(domain_lat),'  max_lat: ',maxval(domain_lat)
        ! if (this%ite < this%its .or. this%jte < this%jts) write(*,*) 'image: ',PE_RANK_GLOBAL+1,'  lat_corners: ',lat_corners,'  lon_corners: ',lon_corners
        if (this%ite < this%its .or. this%jte < this%jts) call io_write('domain_lat.nc',"domain_lat",domain_lat)
        if (this%ite < this%its .or. this%jte < this%jts) call io_write('domain_lon.nc',"domain_lon",domain_lon)
        if (this%ite < this%its .or. this%jte < this%jts) call io_write('boundary_lat.nc',"lat",temp_lat)
        if (this%ite < this%its .or. this%jte < this%jts) call io_write('boundary_lon.nc',"lon",temp_lon)

    end subroutine set_boundary_image


    !>------------------------------------------------------------
    !! Setup the main geo structure in the boundary structure
    !!
    !!------------------------------------------------------------
    subroutine setup_boundary_geo(this, longitude_system)
        implicit none
        type(boundary_t), intent(inout) :: this
        integer,          intent(in)    :: longitude_system


        if (allocated(this%geo%lat)) deallocate(this%geo%lat)
        allocate( this%geo%lat, source=this%lat)

        if (allocated(this%geo%lon)) deallocate(this%geo%lon)
        allocate( this%geo%lon, source=this%lon)

        call standardize_coordinates(this%geo, longitude_system)

        this%geo_u   = this%geo
        this%geo_v   = this%geo
        this%geo_agl = this%geo


        if ( allocated(this%z) )  then
            ! geo%z will be interpolated from this%z to the high-res grids for vinterp in domain... not a great separation
            ! here we save the original z dataset so that it can be used to interpolate varying z through time.
            if (allocated(this%original_geo%z)) deallocate(this%original_geo%z)
            allocate( this%original_geo%z, source=this%z)
        endif

    end subroutine



    !>------------------------------------------------------------
    !! Reads and adds a variable into the variable dictionary
    !!
    !! Given a filename, varname, number of dimensions (2,3) and timestep
    !! Read the timestep of varname from filename and stores the result in a variable structure
    !!
    !! Variable is then added to a master variable dictionary
    !!
    !!------------------------------------------------------------
    subroutine add_var_to_dict(this, var_name, ndims, indx, dims)
        implicit none
        type(boundary_t), intent(inout) :: this
        character(len=*), intent(in)    :: var_name
        integer,          intent(in)    :: ndims, indx
        integer,          intent(in)    :: dims(3)

        real, allocatable :: temp_2d_data(:,:)
        real, allocatable :: temp_3d_data(:,:,:)
        type(variable_t)  :: new_variable

        if (ndims==2) then
            call new_variable%initialize( indx, [dims(1),dims(3)] )
            call this%variables%add_var(indx, new_variable)

        elseif (ndims==3) then
            call new_variable%initialize( indx, dims )
            call this%variables%add_var(indx, new_variable)

        ! these variables are computed (e.g. pressure from height or height from pressure)
        elseif (ndims==-3) then
            call new_variable%initialize( indx, dims )
            new_variable%computed = .True.

            call this%variables%add_var(indx, new_variable)
        endif

    end subroutine

    module subroutine interpolate_original_levels(this, options)
        implicit none
        class(boundary_t),   intent(inout)   :: this
        type(options_t),     intent(in)      :: options

        type(variable_t)        :: input_z, var
        type(interpolable_type) :: input_geo
        real, allocatable :: temp_3d(:,:,:)
        real :: neg_z
        integer :: id


        associate(list => this%variables)

        input_z = list%get_var(kVARS%z)

        if (options%forcing%z_is_geopotential) then
            input_z%data_3d = input_z%data_3d / gravity
            !neg_z = minval(input_z%data_3d)
            ! if (neg_z < 0.0) input_z%data_3d = input_z%data_3d - neg_z
            call list%add_var(kVARS%z, input_z)
        endif

        ! if (options%forcing%z_is_on_interface) then
        !     call interpolate_in_z(input_z%data_3d)
        ! endif

        allocate(input_geo%z, source=input_z%data_3d)

        ! create a vertical interpolation look up table for the current time step
        call vLUT(this%original_geo, input_geo)


        ! loop through the list of variables that were read in and might need to be interpolated in 3D
        call list%reset_iterator()
        do while (list%has_more_elements())
            ! get the next variable in the structure
            var = list%next(id)

            if (var%three_d) then
                ! need to vinterp this dataset to the original vertical levels (if necessary)

                temp_3d = var%data_3d
                call vinterp(var%data_3d, temp_3d, input_geo%vert_lut)

            endif
            call list%add_var(id, var)

        end do
        end associate

    end subroutine


    module subroutine update_computed_vars(this, options)
        implicit none
        class(boundary_t),   intent(inout)   :: this
        type(options_t),     intent(in)      :: options

        integer           :: err
        type(variable_t)  :: var, pvar, tvar

        integer :: nx,ny,nz, var_id
        real, allocatable :: temp_z(:,:,:)
        real :: neg_z

        associate(list => this%variables)

        if (options%forcing%qv_is_relative_humidity) then
            call compute_mixing_ratio_from_rh(list, options)
        endif

        if (options%forcing%qv_is_spec_humidity) then
            call compute_mixing_ratio_from_sh(list, options)
        endif

        if (options%forcing%t_offset /= 0) then
            tvar = list%get_var(kVARS%potential_temperature)
            tvar%data_3d = tvar%data_3d + options%forcing%t_offset
            call list%add_var(kVARS%potential_temperature, tvar)
        endif

        ! loop through the list of variables that need to be read in
        call list%reset_iterator()
        do while (list%has_more_elements())

            ! get the next variable in the structure
            var = list%next(var_id)

            if (var%computed) then

                if (var_id == kVARS%z) then
                    call compute_z_update(this, list, options)
                endif

                if (var_id == kVARS%pressure) then
                    call compute_p_update(this, list, options)
                endif

            endif
        end do

        if (.not.options%forcing%t_is_potential) then
            tvar = list%get_var(kVARS%potential_temperature)
            pvar = list%get_var(kVARS%pressure)    
            tvar%data_3d = tvar%data_3d / exner_function(pvar%data_3d)
            call list%add_var(kVARS%potential_temperature, tvar)
        endif

        if (options%forcing%limit_rh) call limit_rh(list, options)
        !limit sea surface temperature to be >= 273.15 (i.e. cannot freeze)
        call limit_2d_var(list, kVARS%sst, min_val=273.15)

        end associate

        ! if the vertical levels of the forcing data change over time, they need to be interpolated to the original levels here.
        if (options%forcing%time_varying_z) then
            call this%interpolate_original_levels(options)
        endif


    end subroutine update_computed_vars


    subroutine compute_mixing_ratio_from_rh(list, options)
        implicit none
        type(var_dict_t),   intent(inout)   :: list
        type(options_t),    intent(in)      :: options

        integer           :: err
        type(variable_t)  :: pvar, tvar, qvar

        tvar = list%get_var(kVARS%potential_temperature)
        pvar = list%get_var(kVARS%pressure)
        qvar = list%get_var(kVARS%water_vapor)

        if (maxval(qvar%data_3d) > 2) then
            qvar%data_3d = qvar%data_3d/100.0
        endif

        if (options%forcing%t_is_potential) then
            qvar%data_3d = rh_to_mr(qvar%data_3d, tvar%data_3d * exner_function(pvar%data_3d), pvar%data_3d)
        else
            qvar%data_3d = rh_to_mr(qvar%data_3d, tvar%data_3d, pvar%data_3d)
        endif

        call list%add_var(kVARS%water_vapor, qvar)

    end subroutine compute_mixing_ratio_from_rh


    subroutine limit_rh(list, options)
        implicit none
        type(var_dict_t),   intent(inout)   :: list
        type(options_t),    intent(in)      :: options

        integer           :: err
        type(variable_t)  :: pvar, tvar, qvar
        real :: rh, t
        integer :: i,j,k

        tvar = list%get_var(kVARS%potential_temperature)
        pvar = list%get_var(kVARS%pressure)
        qvar = list%get_var(kVARS%water_vapor)

        do j = lbound(tvar%data_3d, 3), ubound(tvar%data_3d, 3)
            do k = lbound(tvar%data_3d, 2), ubound(tvar%data_3d, 2)
                do i = lbound(tvar%data_3d, 1), ubound(tvar%data_3d, 1)

                    t = tvar%data_3d(i,k,j) * exner_function(pvar%data_3d(i,k,j))

                    rh = relative_humidity(t, qvar%data_3d(i,k,j), pvar%data_3d(i,k,j))

                    if (rh > 1.0) then
                        qvar%data_3d = rh_to_mr(1.0, t, pvar%data_3d(i,k,j))
                    endif

                enddo
            enddo
        enddo

        call list%add_var(kVARS%water_vapor, qvar)

    end subroutine limit_rh

    ! set minimum and/or maximum values on any 2D forcing variable (e.g. sst_var, rain_var)
    subroutine limit_2d_var(list, var_id, min_val, max_val)
        implicit none
        type(var_dict_t),   intent(inout)   :: list
        integer,            intent(in)      :: var_id
        real, optional,     intent(in)      :: min_val
        real, optional,     intent(in)      :: max_val

        type(variable_t)  :: var
        integer :: err

        var = list%get_var(var_id,err=err)
        
        if (err > 0) return

        if (present(min_val)) then
            where(var%data_2d < min_val) var%data_2d = min_val
        endif

        if (present(max_val)) then
            where(var%data_2d > max_val) var%data_2d = max_val
        endif
        call list%add_var(var_id, var)

    end subroutine limit_2d_var

    subroutine compute_mixing_ratio_from_sh(list, options)
        implicit none
        type(var_dict_t),   intent(inout)   :: list
        type(options_t),    intent(in)      :: options

        integer           :: err
        type(variable_t)  :: qvar

        qvar = list%get_var(kVARS%water_vapor)

        qvar%data_3d = qvar%data_3d / (1 - qvar%data_3d)

        call list%add_var(kVARS%water_vapor, qvar)

    end subroutine compute_mixing_ratio_from_sh


    subroutine compute_z_update(this, list, options)
        implicit none
        class(boundary_t),  intent(inout)   :: this
        type(var_dict_t),   intent(inout)   :: list
        type(options_t),    intent(in)      :: options

        integer           :: err
        type(variable_t)  :: var, pvar, zvar, tvar, qvar

        real, allocatable :: t(:,:,:)

        qvar = list%get_var(kVARS%water_vapor)
        tvar = list%get_var(kVARS%potential_temperature)
        zvar = list%get_var(kVARS%terrain)
        var = list%get_var(kVARS%pressure, err)

        pvar = list%get_var(kVARS%sea_surface_pressure, err)

        allocate(t, mold=tvar%data_3d)
        if (options%forcing%t_is_potential) then
            ! stop "Need real air temperature to compute height"
            t = exner_function(var%data_3d) * tvar%data_3d
        else
            t = tvar%data_3d
        endif

        if (err == 0) then
            call compute_3d_z(var%data_3d, pvar%data_2d, this%z, t, qvar%data_3d)

        else
            pvar = list%get_var(kVARS%surface_pressure, err)
            if (err == 0) then
                call compute_3d_z(var%data_3d, pvar%data_2d, this%z, t, qvar%data_3d, zvar%data_2d)
            else
                write(*,*) "ERROR reading surface pressure or sea level pressure, variables not found"
                error stop
            endif
        endif
        zvar = list%get_var(kVARS%z)
        zvar%data_3d = this%z
        call list%add_var(kVARS%z, zvar)

    end subroutine compute_z_update

    subroutine compute_p_update(this, list, options)
        implicit none
        class(boundary_t),  intent(inout)   :: this
        type(var_dict_t),   intent(inout)   :: list
        type(options_t),    intent(in)      :: options

        integer           :: err, hgterr
        type(variable_t)  :: pvar, zvar, tvar, qvar, psvar

        if (options%forcing%t_is_potential) stop "Need real air temperature to compute pressure"

        qvar = list%get_var(kVARS%water_vapor)
        tvar = list%get_var(kVARS%potential_temperature)
        pvar = list%get_var(kVARS%pressure)

        psvar = list%get_var(kVARS%surface_pressure, err)
        zvar = list%get_var(kVARS%terrain, hgterr)

        if ((err + hgterr) == 0) then
            call compute_3d_p(pvar%data_3d, psvar%data_2d, this%z, tvar%data_3d, qvar%data_3d, zvar%data_2d)
        else
            psvar = list%get_var(kVARS%sea_surface_pressure, err)

            if (err == 0) then
                call compute_3d_p(pvar%data_3d, psvar%data_2d, this%z, tvar%data_3d, qvar%data_3d)
            else
                write(*,*) "ERROR reading surface pressure or sea level pressure, variables not found"
                write(*,*) "ERROR or forcing height not given if sea level pressure is given"
                error stop
            endif
        endif
        call list%add_var(kVARS%pressure, pvar)


    end subroutine compute_p_update



    !>------------------------------------------------------------
    !! Setup the vars_to_read and var_dimensions arrays given a master set of variables
    !!
    !! Count the number of variables specified, then allocate and store those variables in a list just their size.
    !! The master list will have all variables, but not all will be set
    !!------------------------------------------------------------
    subroutine setup_variable_lists(master_var_list, master_dim_list, opt, vars_to_read, var_dimensions, var_indx)
        implicit none
        character(len=kMAX_NAME_LENGTH), intent(in)                 :: master_var_list(:)
        integer,                         intent(in)                 :: master_dim_list(:)
        type(options_t),                 intent(in)                 :: opt
        character(len=kMAX_NAME_LENGTH), intent(inout), allocatable :: vars_to_read(:)
        integer,                         intent(inout), allocatable :: var_dimensions(:), var_indx(:)

        integer :: n_valid_vars
        integer :: i, curvar, err

        n_valid_vars = 0
        do i=1, size(master_var_list)
            if (trim(master_var_list(i)) /= '') then
                n_valid_vars = n_valid_vars + 1
            endif
        enddo

        allocate(vars_to_read(  n_valid_vars), stat=err)
        if (err /= 0) stop "vars_to_read: Allocation request denied"

        allocate(var_dimensions(  n_valid_vars), stat=err)
        if (err /= 0) stop "var_dimensions: Allocation request denied"

        allocate(var_indx(  n_valid_vars), stat=err)
        if (err /= 0) stop "var_dimensions: Allocation request denied"

        curvar = 1
        do i=1, size(master_var_list)
            if (trim(master_var_list(i)) /= '') then
                vars_to_read(curvar) = master_var_list(i)
                var_dimensions(curvar) = master_dim_list(i)
                ! if (STD_OUT_PE) print *, "in variable list: ", vars_to_read(curvar)
                curvar = curvar + 1
            endif
        enddo

        do i=1,size(vars_to_read)
            if (vars_to_read(i) == opt%forcing%latvar) var_indx(i) = kVARS%latitude
            if (vars_to_read(i) == opt%forcing%lonvar) var_indx(i) = kVARS%longitude
            if (vars_to_read(i) == opt%forcing%uvar) var_indx(i) = kVARS%u
            if (vars_to_read(i) == opt%forcing%ulat) var_indx(i) = kVARS%u_latitude
            if (vars_to_read(i) == opt%forcing%ulon) var_indx(i) = kVARS%u_longitude
            if (vars_to_read(i) == opt%forcing%vvar) var_indx(i) = kVARS%v
            if (vars_to_read(i) == opt%forcing%vlat) var_indx(i) = kVARS%v_latitude
            if (vars_to_read(i) == opt%forcing%vlon) var_indx(i) = kVARS%v_longitude
            if (vars_to_read(i) == opt%forcing%wvar) var_indx(i) = kVARS%w_real
            if (vars_to_read(i) == opt%forcing%pvar) var_indx(i) = kVARS%pressure
            if (vars_to_read(i) == opt%forcing%tvar) var_indx(i) = kVARS%potential_temperature
            if (vars_to_read(i) == opt%forcing%qvvar) var_indx(i) = kVARS%water_vapor
            if (vars_to_read(i) == opt%forcing%qcvar) var_indx(i) = kVARS%cloud_water_mass
            if (vars_to_read(i) == opt%forcing%qivar) var_indx(i) = kVARS%ice_mass
            if (vars_to_read(i) == opt%forcing%qrvar) var_indx(i) = kVARS%rain_mass
            if (vars_to_read(i) == opt%forcing%qsvar) var_indx(i) = kVARS%snow_mass
            if (vars_to_read(i) == opt%forcing%qgvar) var_indx(i) = kVARS%graupel_mass
            if (vars_to_read(i) == opt%forcing%i2mvar) var_indx(i) = kVARS%ice2_mass
            if (vars_to_read(i) == opt%forcing%i3mvar) var_indx(i) = kVARS%ice3_mass
            if (vars_to_read(i) == opt%forcing%qncvar) var_indx(i) = kVARS%cloud_number
            if (vars_to_read(i) == opt%forcing%qnivar) var_indx(i) = kVARS%ice_number
            if (vars_to_read(i) == opt%forcing%qnrvar) var_indx(i) = kVARS%rain_number
            if (vars_to_read(i) == opt%forcing%qnsvar) var_indx(i) = kVARS%snow_number
            if (vars_to_read(i) == opt%forcing%qngvar) var_indx(i) = kVARS%graupel_number
            if (vars_to_read(i) == opt%forcing%i2nvar) var_indx(i) = kVARS%ice2_number
            if (vars_to_read(i) == opt%forcing%i3nvar) var_indx(i) = kVARS%ice3_number
            if (vars_to_read(i) == opt%forcing%i1avar) var_indx(i) = kVARS%ice1_a
            if (vars_to_read(i) == opt%forcing%i1cvar) var_indx(i) = kVARS%ice1_c
            if (vars_to_read(i) == opt%forcing%i2avar) var_indx(i) = kVARS%ice2_a
            if (vars_to_read(i) == opt%forcing%i2cvar) var_indx(i) = kVARS%ice2_c
            if (vars_to_read(i) == opt%forcing%i3avar) var_indx(i) = kVARS%ice3_a
            if (vars_to_read(i) == opt%forcing%i3cvar) var_indx(i) = kVARS%ice3_c
            if (vars_to_read(i) == opt%forcing%hgtvar) var_indx(i) = kVARS%terrain
            if (vars_to_read(i) == opt%forcing%pslvar) var_indx(i) = kVARS%sea_surface_pressure
            if (vars_to_read(i) == opt%forcing%psvar) var_indx(i) = kVARS%surface_pressure
            if (vars_to_read(i) == opt%forcing%sst_var) var_indx(i) = kVARS%sst
            if (vars_to_read(i) == opt%forcing%pblhvar) var_indx(i) = kVARS%hpbl
            if (vars_to_read(i) == opt%forcing%shvar) var_indx(i) = kVARS%sensible_heat
            if (vars_to_read(i) == opt%forcing%lhvar) var_indx(i) = kVARS%latent_heat
            if (vars_to_read(i) == opt%forcing%zvar) var_indx(i) = kVARS%z
            if (vars_to_read(i) == opt%forcing%swdown_var) var_indx(i) = kVARS%shortwave
            if (vars_to_read(i) == opt%forcing%lwdown_var) var_indx(i) = kVARS%longwave
            if (vars_to_read(i) == opt%forcing%time_var) var_indx(i) = 1000

        enddo
    end subroutine

end submodule
