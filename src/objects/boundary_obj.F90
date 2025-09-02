!>------------------------------------------------------------
!!  Handles reading boundary conditions from the forcing file(s)
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!------------------------------------------------------------
submodule(boundary_interface) boundary_implementation
    use array_utilities,        only : interpolate_in_z, array_offset_x, array_offset_y
    use io_routines,            only : io_getdims, io_read, io_maxDims, io_variable_is_present, io_write
    use time_io,                only : read_times, find_timestep_in_filelist
    use string,                 only : str
    use mod_atm_utilities,      only : rh_to_mr, relative_humidity, compute_3d_p, compute_3d_z, exner_function
    use geo,                    only : standardize_coordinates, geo_interp, geo_lut
    use vertical_interpolation, only : vLUT, vinterp
    use timer_interface,    only : timer_t
    use debug_module,           only : check_ncdf
    use mod_wrf_constants,      only : gravity
    use icar_constants,         only : STD_OUT_PE, kMAX_FILE_LENGTH
    implicit none
contains

    !>------------------------------------------------------------
    !! Set default component values
    !! Reads initial conditions from the forcing file into image 1
    !!
    !! Distributes initial conditions to all other images
    !!
    !!------------------------------------------------------------
    module subroutine init(this, options, domain_lat, domain_lon, parent_options)
        class(boundary_t),    intent(inout) :: this
        type(options_t),      intent(inout) :: options
        real, dimension(:,:), intent(in)    :: domain_lat
        real, dimension(:,:), intent(in)    :: domain_lon
        type(options_t), optional, intent(in)    :: parent_options

        type(Time_type) :: strt_time
        character(len=kMAX_NAME_LENGTH), allocatable :: vars_to_read(:)
        type(dim_arrays_type),           allocatable :: var_dimensions(:)


        ! the parameters option type can't contain allocatable arrays because it is a coarray
        ! so we need to allocate the vars_to_read and var_dimensions outside of the options type
        call setup_variable_lists(options%forcing%vars_to_read, options%forcing%dim_list, vars_to_read, var_dimensions)

        strt_time = options%general%start_time
        if (options%restart%restart) strt_time = options%restart%restart_time

        ! Read through forcing variable names stored in "options"
        ! needs to read each one to find the grid information for it
        ! then create grid and initialize a variable...
        ! also need to explicitly save lat and lon data
        ! if (STD_OUT_PE) then
        if (present(parent_options)) then
            call this%init_local_asnest(vars_to_read, var_dimensions,   &
                                    domain_lat,        &
                                    domain_lon,       &
                                    parent_options)
        else
            call this%init_local(options%forcing,                           &
                                    vars_to_read, var_dimensions,   &
                                    strt_time,                      &
                                    domain_lat,        &
                                    domain_lon)
        endif
        ! endif
        ! call this%distribute_initial_conditions()


        call setup_boundary_geo(this, options%domain%longitude_system, options)

    end subroutine

    subroutine read_latlon(file,latvar,lonvar,lat_out,lon_out)
        implicit none
        character(len=*), intent(in) :: file, latvar, lonvar
        real, allocatable, intent(out) :: lat_out(:,:), lon_out(:,:)
        real, allocatable :: lat_1d(:), lon_1d(:)
        integer, allocatable :: lat_dims(:), lon_dims(:)

        integer :: i, x_len, y_len

        !See if lat and lon are 1D or 2D
        call io_getdims(file, latvar, lat_dims)
        call io_getdims(file, lonvar, lon_dims)

        if (size(lat_dims) == 1) then
            call io_read(file, latvar, lat_1d)

            !This will always be the case for 1D or 2D lon
            x_len = lon_dims(1)

            allocate(lat_out(1:x_len,1:lat_dims(1)))
            do i=1,x_len
                lat_out(i,:) = lat_1d
            end do
        elseif (size(lat_dims) == 2) then
            call io_read(file, latvar, lat_out)
        else
            write(*,*) 'ERROR: lat dimension on forcing data is not 1D or 2D'
            stop
        endif

        if (size(lon_dims) == 1) then
            call io_read(file, lonvar, lon_1d)

            if (size(lat_dims) == 1) then
                y_len = lat_dims(1)
            else
                y_len = lat_dims(2)
            endif

            allocate(lon_out(1:lon_dims(1),1:y_len))
            do i=1,y_len
                lon_out(:,i) = lon_1d
            end do
        elseif (size(lon_dims) == 2) then
            call io_read(file, lonvar, lon_out)
        else
            write(*,*) 'ERROR: lon dimension on forcing data is not 1D or 2D'
            stop
        endif

    end subroutine read_latlon

    !>------------------------------------------------------------
    !! Set default component values
    !! Reads initial conditions from the forcing file
    !!
    !!------------------------------------------------------------
    module subroutine init_local(this, options, var_list, dim_list, start_time, domain_lat, domain_lon)
        implicit none
        class(boundary_t),               intent(inout)  :: this
        type(forcing_options_type),      intent(inout)  :: options
        character(len=kMAX_NAME_LENGTH), intent(in)     :: var_list (:)
        type(dim_arrays_type),           intent(in)     :: dim_list (:)
        type(Time_type),                 intent(in)     :: start_time
        real, dimension(:,:),            intent(in)     :: domain_lat
        real, dimension(:,:),            intent(in)     :: domain_lon

        type(variable_t)  :: zvar
        real, allocatable :: temp_3d(:,:,:), temp_z_trans(:,:,:), temp_lat(:,:), temp_lon(:,:), lat_1d(:), lon_1d(:)
        real :: neg_z
        integer :: i, nx, ny, nz, PE_RANK_GLOBAL
        logical :: z_staggered

        ! figure out while file and timestep contains the requested start_time
        call set_firstfile_firststep(this, start_time, options%boundary_files, options%time_var)

        call read_latlon(this%firstfile, options%latvar, options%lonvar, temp_lat, temp_lon)
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

        ! read in the u and v lat and lon data if specified
        if (options%ulat /= "" .and. options%ulon /= "") then

            call io_read(this%firstfile, options%uvar,   temp_3d,   this%firststep)
            if (.not.(size(temp_3d,1) == this%ide+1 .and. size(temp_3d,2) == this%jde)) then
                write(*,*) "ERROR: forcing variable uvar is not on the staggered u-grid, but lat/lon data for this grid (ulat and ulon) has been provided."
                write(*,*) "ERROR: please make the dimensions of uvar and ulat/ulon consistent"
                write(*,*) "       uvar size: ", size(temp_3d,1), " ", size(temp_3d,2)
                write(*,*) "       lat/lon size: ", this%ide, " ", this%jde
                stop
            endif

            allocate(this%ulat((this%ite-this%its+2),(this%jte-this%jts+1)))
            allocate(this%ulon((this%ite-this%its+2),(this%jte-this%jts+1)))
            call read_latlon(this%firstfile, options%ulat, options%ulon, temp_lat, temp_lon)

            this%ulat = temp_lat(this%its:this%ite+1,this%jts:this%jte)
            this%ulon = temp_lon(this%its:this%ite+1,this%jts:this%jte)
        elseif(options%ulat == "" .and. options%ulon == "") then
            !ensure that uvar has the same dimensions at lat and lon
            !read in uvar
            call io_read(this%firstfile, options%uvar,   temp_3d,   this%firststep)
            if (.not.(size(temp_3d,1) == this%ide .and. size(temp_3d,2) == this%jde)) then
                write(*,*) "ERROR: forcing variable uvar is not on the mass-grid, and no lat/lon data for this grid (ulat and ulon) has been provided"
                write(*,*) "       uvar size: ", size(temp_3d,1), " ", size(temp_3d,2)
                write(*,*) "       lat/lon size: ", this%ide, " ", this%jde
                stop
            endif
        else
            write(*,*) "ERROR: ulat and ulon must both be set or both be unset for forcing data"
            stop
        endif

        if (options%vlat /= "" .and. options%vlon /= "") then

            call io_read(this%firstfile, options%vvar,   temp_3d,   this%firststep)
            if (.not.(size(temp_3d,1) == this%ide .and. size(temp_3d,2) == this%jde+1)) then
                write(*,*) "ERROR: forcing variable vvar is not on the staggered v-grid, but lat/lon data for this grid (vlat and vlon) has been provided."
                write(*,*) "ERROR: please make the dimensions of vvar and vlat/vlon consistent"
                write(*,*) "       vvar size: ", size(temp_3d,1), " ", size(temp_3d,2)
                write(*,*) "       lat/lon size: ", this%ide, " ", this%jde
                stop
            endif

            allocate(this%vlat((this%ite-this%its+1),(this%jte-this%jts+2)))
            allocate(this%vlon((this%ite-this%its+1),(this%jte-this%jts+2)))

            call read_latlon(this%firstfile, options%vlat, options%vlon, temp_lat, temp_lon)

            this%vlat = temp_lat(this%its:this%ite,this%jts:this%jte+1)
            this%vlon = temp_lon(this%its:this%ite,this%jts:this%jte+1)
        elseif(options%vlat == "" .and. options%vlon == "") then
            !ensure that vvar has the same dimensions at lat and lon
            !read in vvar
            call io_read(this%firstfile, options%vvar,   temp_3d,   this%firststep)
            if (.not.(size(temp_3d,1) == this%ide .and. size(temp_3d,2) == this%jde)) then
                write(*,*) "ERROR: forcing variable vvar is not on the mass-grid, and no lat/lon data for this grid (vlat and vlon) has been provided"
                write(*,*) "       vvar size: ", size(temp_3d,1), " ", size(temp_3d,2)
                write(*,*) "       lat/lon size: ", this%ide, " ", this%jde
                stop
            endif
        else
            write(*,*) "ERROR: vlat and vlon must both be set or both be unset for forcing data"
            stop
        endif

        ! read in the vertical dimension coordinate of the input data
        call io_read(this%firstfile, options%tvar,   temp_3d,   this%firststep)

        nx = size(temp_3d,1)
        ny = size(temp_3d,2)
        nz = size(temp_3d,3)

        if (allocated(this%z)) deallocate(this%z)
        allocate(this%z((this%ite-this%its+1),nz,(this%jte-this%jts+1)))

        this%kts = 1
        this%kte = nz

        if (this%ite < this%its) write(*,*) 'image: ',PE_RANK_GLOBAL+1,'  its: ',this%its,'  ite: ',this%ite,'  jts: ',this%jts,'  jte: ',this%jte
        if (this%kte < this%kts) write(*,*) 'image: ',PE_RANK_GLOBAL+1,'  its: ',this%its,'  ite: ',this%ite,'  jts: ',this%jts,'  jte: ',this%jte
        if (this%jte < this%jts) write(*,*) 'image: ',PE_RANK_GLOBAL+1,'  its: ',this%its,'  ite: ',this%ite,'  jts: ',this%jts,'  jte: ',this%jte

        ! call assert(size(var_list) == size(dim_list), "list of variable dimensions must match list of variables")
        do i=1, size(var_list)
            call add_var_to_dict(this, var_list(i), dim_list(i)%num_dims, dim_list(i)%dims)
        end do

        if (.not. options%compute_z) then
            z_staggered = .False.

            call io_read(this%firstfile, options%zvar,   temp_3d,   this%firststep)
            nx = size(temp_3d,1)
            ny = size(temp_3d,2)

            ! handle case that height is staggered in z dimension
            if (size(temp_3d,3) == nz+1) then
                nz = size(temp_3d,3)
                z_staggered = .True.
            elseif(size(temp_3d,3) == nz) then
                nz = size(temp_3d,3)
            else
                !error
                write(*,*) "ERROR: Incompatible vertical dimension sizes on forcing variable: ",options%zvar
                write(*,*) "  Expected: ", nz, " or ", nz+1, " Found: ", size(temp_3d,3)
                stop
            endif

            allocate(temp_z_trans(1:nx,1:nz,1:ny))
            
            temp_z_trans(1:nx,1:nz,1:ny) = reshape(temp_3d, shape=[nx,nz,ny], order=[1,3,2])

            if (z_staggered) call interpolate_in_z(temp_z_trans)

            zvar = this%variables%get_var(options%zvar)
            zvar%data_3d = temp_z_trans(this%its:this%ite,this%kts:this%kte,this%jts:this%jte)
            call this%variables%add_var(options%zvar, zvar)
        endif

    end subroutine

    !>------------------------------------------------------------
    !! Set default component values
    !! Reads initial conditions from the forcing file
    !!
    !!------------------------------------------------------------
    module subroutine init_local_asnest(this, var_list, dim_list, domain_lat, domain_lon, parent_options)
        class(boundary_t),               intent(inout)  :: this
        character(len=kMAX_NAME_LENGTH), intent(in)     :: var_list (:)
        type(dim_arrays_type),           intent(in)     :: dim_list (:)
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
        
        call read_latlon(parent_options%domain%init_conditions_file, parent_options%domain%lat_hi, parent_options%domain%lon_hi, parent_nest_lat, parent_nest_lon)

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

        allocate(this%ulat((this%ite-this%its+2),(this%jte-this%jts+1)))
        allocate(this%ulon((this%ite-this%its+2),(this%jte-this%jts+1)))

        !see if ulat and ulon were given to parent domain
        if (parent_options%domain%ulon_hi /= "" .and. parent_options%domain%ulat_hi /= "") then
            call read_latlon(parent_options%domain%init_conditions_file, parent_options%domain%ulat_hi, parent_options%domain%ulon_hi, parent_nest_lat, parent_nest_lon)

            this%ulat = parent_nest_lat(this%its:this%ite+1,this%jts:this%jte)
            this%ulon = parent_nest_lon(this%its:this%ite+1,this%jts:this%jte)
        else
            call array_offset_x(this%lon, this%ulon)
            call array_offset_x(this%lat, this%ulat)
        endif

        allocate(this%vlat((this%ite-this%its+1),(this%jte-this%jts+2)))
        allocate(this%vlon((this%ite-this%its+1),(this%jte-this%jts+2)))

        !see if vlat and vlon were given to parent domain
        if (parent_options%domain%vlon_hi /= "" .and. parent_options%domain%vlat_hi /= "") then
            call read_latlon(parent_options%domain%init_conditions_file, parent_options%domain%vlat_hi, parent_options%domain%vlon_hi, parent_nest_lat, parent_nest_lon)

            this%vlat = parent_nest_lat(this%its:this%ite,this%jts:this%jte+1)
            this%vlon = parent_nest_lon(this%its:this%ite,this%jts:this%jte+1)
        else
            call array_offset_y(this%lon, this%vlon)
            call array_offset_y(this%lat, this%vlat)
        endif

        ! Get the height coordinate of the parent nest
        this%kts = 1
        this%kte = parent_options%domain%nz

        if (allocated(this%z)) deallocate(this%z)
        allocate(this%z((this%ite-this%its+1),this%kte,(this%jte-this%jts+1)))
        this%z = 0.0 !parent_nest_z(this%its:this%ite,1:this%kte,this%jts:this%jte)

        ! if (this%ite < this%its) write(*,*) 'image: ',PE_RANK_GLOBAL+1,'  its: ',this%its,'  ite: ',this%ite,'  jts: ',this%jts,'  jte: ',this%jte
        ! if (this%kte < this%kts) write(*,*) 'image: ',PE_RANK_GLOBAL+1,'  its: ',this%its,'  ite: ',this%ite,'  jts: ',this%jts,'  jte: ',this%jte
        ! if (this%jte < this%jts) write(*,*) 'image: ',PE_RANK_GLOBAL+1,'  its: ',this%its,'  ite: ',this%ite,'  jts: ',this%jts,'  jte: ',this%jte

        ! call assert(size(var_list) == size(dim_list), "list of variable dimensions must match list of variables")
        do i=1, size(var_list)
            call add_var_to_dict(this, var_list(i), dim_list(i)%num_dims, dim_list(i)%dims)
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
            if (STD_OUT_PE) write(*,*) 'ERROR: Time: ',time%as_string()
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
        this%ids = lbound(temp_lat,1)
        this%ide = ubound(temp_lat,1)
        this%jds = lbound(temp_lon,2)
        this%jde = ubound(temp_lon,2)

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
    subroutine setup_boundary_geo(this, longitude_system, options)
        implicit none
        type(boundary_t), intent(inout) :: this
        integer,          intent(in)    :: longitude_system
        type(options_t),  intent(in)    :: options

        if (allocated(this%geo%lat)) deallocate(this%geo%lat)
        allocate( this%geo%lat, source=this%lat)

        if (allocated(this%geo%lon)) deallocate(this%geo%lon)
        allocate( this%geo%lon, source=this%lon)

        call standardize_coordinates(this%geo, longitude_system)

        if (allocated(this%ulat) .and. allocated(this%ulon)) then
            if (allocated(this%geo_u%lat)) deallocate(this%geo_u%lat)
            allocate( this%geo_u%lat, source=this%ulat)

            if (allocated(this%geo_u%lon)) deallocate(this%geo_u%lon)
            allocate( this%geo_u%lon, source=this%ulon)

            call standardize_coordinates(this%geo_u, longitude_system)
        else
            this%geo_u   = this%geo
        endif

        if (allocated(this%vlat) .and. allocated(this%vlon)) then
            if (allocated(this%geo_v%lat)) deallocate(this%geo_v%lat)
            allocate( this%geo_v%lat, source=this%vlat)

            if (allocated(this%geo_v%lon)) deallocate(this%geo_v%lon)
            allocate( this%geo_v%lon, source=this%vlon)

            call standardize_coordinates(this%geo_v, longitude_system)
        else
            this%geo_v   = this%geo
        endif

        this%geo_agl = this%geo

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
    subroutine add_var_to_dict(this, var_name, ndims, dims)
        implicit none
        type(boundary_t), intent(inout) :: this
        character(len=*), intent(in)    :: var_name
        integer,          intent(in)    :: ndims
        integer,          intent(in)    :: dims(3)

        real, allocatable :: temp_2d_data(:,:)
        real, allocatable :: temp_3d_data(:,:,:)
        type(variable_t)  :: new_variable
        logical :: computed_flag
        integer :: nx, ny, nz

        ! these variables are computed (e.g. pressure from height or height from pressure)
        computed_flag = (len_trim(var_name) > 9 .and. var_name(max(1,len_trim(var_name)-8):len_trim(var_name)) == "_computed")

        nx = (this%ite - this%its + 1)
        nz = this%kte
        ny = (this%jte - this%jts + 1)

        !handle staggered variables
        if (.not.(computed_flag)) then
            if (dims(1) == this%ide+1) nx = nx+1
            if (dims(2) == this%jde+1) ny = ny+1
        endif

        if (ndims==2) then
            call new_variable%initialize( [nx,ny] )
            call this%variables%add_var(var_name, new_variable)

        elseif (ndims==3) then
            call new_variable%initialize( [nx,nz,ny] )
            call this%variables%add_var(var_name, new_variable)
        elseif (computed_flag) then
            call new_variable%initialize( [nx,nz,ny] )
            new_variable%computed = .True.

            call this%variables%add_var(var_name, new_variable)
        endif

    end subroutine

    module subroutine setup_z(this, options)
        implicit none
        class(boundary_t),   intent(inout)   :: this
        type(options_t),     intent(in)      :: options

        type(variable_t)        :: input_z, phbase

        associate(list => this%variables)
        input_z = list%get_var(options%forcing%zvar)

        if (options%forcing%z_is_geopotential) then
            ! see if the user has provided a base geopotential height field
            if (options%forcing%phbvar/="") then
                phbase = list%get_var(options%forcing%phbvar)
                input_z%data_3d = input_z%data_3d + phbase%data_3d
            endif
            input_z%data_3d = input_z%data_3d / gravity
            !neg_z = minval(input_z%data_3d)
            ! if (neg_z < 0.0) input_z%data_3d = input_z%data_3d - neg_z
            call list%add_var(options%forcing%zvar, input_z)
        endif

        if (allocated(this%z)) deallocate(this%z)
        allocate(this%z,source=input_z%data_3d)
        end associate

        if (.not.(allocated(this%original_geo%z)) )  then

            ! geo%z will be interpolated from this%z to the high-res grids for vinterp in domain... not a great separation
            ! here we save the original z dataset so that it can be used to interpolate varying z through time.
            allocate( this%original_geo%z, source=this%z)

            this%mass_to_u%lat = this%geo%lat
            this%mass_to_u%lon = this%geo%lon
            this%mass_to_v%lat = this%geo%lat
            this%mass_to_v%lon = this%geo%lon

            call geo_LUT(this%geo_u, this%mass_to_u)
            call geo_LUT(this%geo_v, this%mass_to_v)

            if (allocated(this%original_geo_u%z)) deallocate(this%original_geo_u%z)
            if (allocated(this%original_geo_v%z)) deallocate(this%original_geo_v%z)
            allocate(this%original_geo_u%z(lbound(this%geo_u%lon,1):ubound(this%geo_u%lon,1), this%kts:this%kte, lbound(this%geo_u%lon,2):ubound(this%geo_u%lon,2)))            
            allocate(this%original_geo_v%z(lbound(this%geo_v%lat,1):ubound(this%geo_v%lat,1), this%kts:this%kte, lbound(this%geo_v%lat,2):ubound(this%geo_v%lat,2)))     

            call geo_interp(this%original_geo_u%z,   this%z, this%mass_to_u%geolut)
            call geo_interp(this%original_geo_v%z,   this%z, this%mass_to_v%geolut)
        endif

    end subroutine setup_z

    module subroutine interpolate_original_levels(this, options)
        implicit none
        class(boundary_t),   intent(inout)   :: this
        type(options_t),     intent(in)      :: options

        type(variable_t)        :: input_z, var
        type(interpolable_type) :: input_geo, input_geo_u, input_geo_v
        real, allocatable :: temp_3d(:,:,:)
        real :: neg_z
        character(len=kMAX_NAME_LENGTH) :: name

        call this%setup_z(options)

        allocate(input_geo%z, source=this%z)
        allocate(input_geo_u%z(lbound(this%geo_u%lon,1):(ubound(this%geo_u%lon,1)),lbound(this%z,2):ubound(this%z,2),lbound(this%geo_u%lon,2):ubound(this%geo_u%lon,2)))
        allocate(input_geo_v%z(lbound(this%geo_v%lat,1):ubound(this%geo_v%lat,1),lbound(this%z,2):ubound(this%z,2),lbound(this%geo_v%lat,2):(ubound(this%geo_v%lat,2))))

        call geo_interp(input_geo_u%z, input_geo%z, this%mass_to_u%geolut)
        call geo_interp(input_geo_v%z, input_geo%z, this%mass_to_v%geolut)

        ! create a vertical interpolation look up table for the current time step
        call vLUT(this%original_geo, input_geo)
        call vLUT(this%original_geo_u, input_geo_u)
        call vLUT(this%original_geo_v, input_geo_v)

        associate(list => this%variables)

        ! loop through the list of variables that were read in and might need to be interpolated in 3D
        call list%reset_iterator()
        do while (list%has_more_elements())
            ! get the next variable in the structure
            var = list%next(name)

            if (var%three_d) then
                ! need to vinterp this dataset to the original vertical levels (if necessary)

                temp_3d = var%data_3d

                if (size(temp_3d,1) == (this%ite-this%its+2)) then
                    call vinterp(var%data_3d, temp_3d, input_geo_u%vert_lut)
                elseif (size(temp_3d,3) == (this%jte-this%jts+2)) then
                    call vinterp(var%data_3d, temp_3d, input_geo_v%vert_lut)
                else
                    call vinterp(var%data_3d, temp_3d, input_geo%vert_lut)
                endif

            endif
            call list%add_var(name, var)

        end do
        end associate

    end subroutine


    module subroutine update_computed_vars(this, options)
        implicit none
        class(boundary_t),   intent(inout)   :: this
        type(options_t),     intent(in)      :: options

        integer           :: err
        type(variable_t)  :: var, pvar, tvar, pbvar

        integer :: nx,ny,nz
        real :: neg_z
        character(len=kMAX_NAME_LENGTH) :: name
        logical :: z_is_computed = .False.

        associate(list => this%variables)

        !Add base pressure to pressure, if provided
        pbvar = list%get_var(options%forcing%pbvar, err)
        if (err == 0) then
            pvar = list%get_var(options%forcing%pvar)
            pvar%data_3d = pvar%data_3d + pbvar%data_3d
            call list%add_var(options%forcing%pvar, pvar)
        endif

        if (options%forcing%qv_is_relative_humidity) then
            call compute_mixing_ratio_from_rh(list, options)
        endif

        if (options%forcing%qv_is_spec_humidity) then
            call compute_mixing_ratio_from_sh(list, options)
        endif

        if (options%forcing%t_offset /= 0) then
            tvar = list%get_var(options%forcing%tvar)
            tvar%data_3d = tvar%data_3d + options%forcing%t_offset
            call list%add_var(options%forcing%tvar, tvar)
        endif

        ! loop through the list of variables that need to be read in
        call list%reset_iterator()
        do while (list%has_more_elements())

            ! get the next variable in the structure
            var = list%next(name)

            if (var%computed) then

                if (trim(name) == trim(options%forcing%zvar)) then
                    call compute_z_update(this, list, options)
                    z_is_computed = .True.
                endif

                if (trim(name) == trim(options%forcing%pvar)) then
                    call compute_p_update(this, list, options)
                endif

            endif
        end do

        if (.not.options%forcing%t_is_potential) then
            tvar = list%get_var(options%forcing%tvar)
            pvar = list%get_var(options%forcing%pvar)    
            tvar%data_3d = tvar%data_3d / exner_function(pvar%data_3d)
            call list%add_var(options%forcing%tvar, tvar)
        endif

        if (options%forcing%limit_rh) call limit_rh(list, options)
        !limit sea surface temperature to be >= 273.15 (i.e. cannot freeze)
        call limit_2d_var(list, options%forcing%sst_var, min_val=273.15)

        end associate

        ! if the vertical levels of the forcing data change over time, they need to be interpolated to the original levels here.
        if (options%forcing%time_varying_z .and. .not.(z_is_computed)) then
            call this%interpolate_original_levels(options)
        endif

    end subroutine update_computed_vars


    subroutine compute_mixing_ratio_from_rh(list, options)
        implicit none
        type(var_dict_t),   intent(inout)   :: list
        type(options_t),    intent(in)      :: options

        integer           :: err
        type(variable_t)  :: pvar, tvar, qvar

        tvar = list%get_var(options%forcing%tvar)
        pvar = list%get_var(options%forcing%pvar)
        qvar = list%get_var(options%forcing%qvvar)

        if (pvar%computed) stop "Need pressure as input to compute mixing ratio from relative humidity"
        if (maxval(qvar%data_3d) > 2) then
            qvar%data_3d = qvar%data_3d/100.0
        endif

        if (options%forcing%t_is_potential) then
            qvar%data_3d = rh_to_mr(qvar%data_3d, tvar%data_3d * exner_function(pvar%data_3d), pvar%data_3d)
        else
            qvar%data_3d = rh_to_mr(qvar%data_3d, tvar%data_3d, pvar%data_3d)
        endif

        call list%add_var(options%forcing%qvvar, qvar)

    end subroutine compute_mixing_ratio_from_rh


    subroutine limit_rh(list, options)
        implicit none
        type(var_dict_t),   intent(inout)   :: list
        type(options_t),    intent(in)      :: options

        integer           :: err
        type(variable_t)  :: pvar, tvar, qvar
        real :: rh, t
        integer :: i,j,k

        tvar = list%get_var(options%forcing%tvar)
        pvar = list%get_var(options%forcing%pvar)
        qvar = list%get_var(options%forcing%qvvar)

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

        call list%add_var(options%forcing%qvvar, qvar)

    end subroutine limit_rh

    ! set minimum and/or maximum values on any 2D forcing variable (e.g. sst_var, rain_var)
    subroutine limit_2d_var(list, varname, min_val, max_val)
        implicit none
        type(var_dict_t),   intent(inout)   :: list
        character(len=*),   intent(in)      :: varname
        real, optional,     intent(in)      :: min_val
        real, optional,     intent(in)      :: max_val

        type(variable_t)  :: var

        if (trim(varname) == "") return

        var = list%get_var(varname)

        if (present(min_val)) then
            where(var%data_2d < min_val) var%data_2d = min_val
        endif

        if (present(max_val)) then
            where(var%data_2d > max_val) var%data_2d = max_val
        endif
        call list%add_var(varname, var)

    end subroutine limit_2d_var

    subroutine compute_mixing_ratio_from_sh(list, options)
        implicit none
        type(var_dict_t),   intent(inout)   :: list
        type(options_t),    intent(in)      :: options

        integer           :: err
        type(variable_t)  :: qvar

        qvar = list%get_var(options%forcing%qvvar)

        qvar%data_3d = qvar%data_3d / (1 - qvar%data_3d)

        call list%add_var(options%forcing%qvvar, qvar)

    end subroutine compute_mixing_ratio_from_sh


    subroutine compute_z_update(this, list, options)
        implicit none
        class(boundary_t),  intent(inout)   :: this
        type(var_dict_t),   intent(inout)   :: list
        type(options_t),    intent(in)      :: options

        integer           :: err
        type(variable_t)  :: var, pvar, zvar, tvar, qvar

        real, allocatable :: t(:,:,:)

        qvar = list%get_var(options%forcing%qvvar)
        tvar = list%get_var(options%forcing%tvar)
        zvar = list%get_var(options%forcing%hgtvar)
        var = list%get_var(options%forcing%pvar, err)

        pvar = list%get_var(options%forcing%pslvar, err)

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
            pvar = list%get_var(options%forcing%psvar, err)
            if (err == 0) then
                call compute_3d_z(var%data_3d, pvar%data_2d, this%z, t, qvar%data_3d, zvar%data_2d)
            else
                write(*,*) "ERROR reading surface pressure or sea level pressure, variables not found"
                error stop
            endif
        endif
        zvar = list%get_var(options%forcing%zvar)
        zvar%data_3d = this%z
        call list%add_var(options%forcing%zvar, zvar)

    end subroutine compute_z_update

    subroutine compute_p_update(this, list, options)
        implicit none
        class(boundary_t),  intent(inout)   :: this
        type(var_dict_t),   intent(inout)   :: list
        type(options_t),    intent(in)      :: options

        integer           :: err, hgterr
        type(variable_t)  :: pvar, zvar, tvar, qvar, psvar

        if (options%forcing%t_is_potential) stop "Need real air temperature to compute pressure"

        qvar = list%get_var(options%forcing%qvvar)
        tvar = list%get_var(options%forcing%tvar)
        pvar = list%get_var(options%forcing%pvar)

        psvar = list%get_var(options%forcing%psvar, err)
        zvar = list%get_var(options%forcing%hgtvar, hgterr)

        if ((err + hgterr) == 0) then
            call compute_3d_p(pvar%data_3d, psvar%data_2d, this%z, tvar%data_3d, qvar%data_3d, zvar%data_2d)
        else
            psvar = list%get_var(options%forcing%pslvar, err)

            if (err == 0) then
                call compute_3d_p(pvar%data_3d, psvar%data_2d, this%z, tvar%data_3d, qvar%data_3d)
            else
                write(*,*) "ERROR reading surface pressure or sea level pressure, variables not found"
                write(*,*) "ERROR or forcing height not given if sea level pressure is given"
                error stop
            endif
        endif
        call list%add_var(options%forcing%pvar, pvar)


    end subroutine compute_p_update



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
                vars_to_read(curvar) = master_var_list(i)
                var_dimensions(curvar)%num_dims = master_dim_list(i)%num_dims
                var_dimensions(curvar)%dims = master_dim_list(i)%dims
                ! if (STD_OUT_PE) print *, "in variable list: ", vars_to_read(curvar)
                curvar = curvar + 1
            endif
        enddo
    end subroutine

end submodule
