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
    use string,                 only : str
    use mod_atm_utilities,      only : rh_to_mr, relative_humidity, compute_3d_p, compute_3d_z, exner_function
    use geo,                    only : standardize_coordinates
    use vertical_interpolation, only : vLUT, vinterp
    use timer_interface,    only : timer_t
    use debug_module,           only : check_ncdf
    use mod_wrf_constants,      only : gravity
    implicit none
contains

    !>------------------------------------------------------------
    !! Set default component values
    !! Reads initial conditions from the forcing file into image 1
    !!
    !! Distributes initial conditions to all other images
    !!
    !!------------------------------------------------------------
    module subroutine init(this, options, domain_lat, domain_lon, domain_vars)
        class(boundary_t), intent(inout) :: this
        type(options_t),   intent(inout) :: options
        real, dimension(:,:), intent(in)     :: domain_lat
        real, dimension(:,:), intent(in)     :: domain_lon
        type(var_dict_t),     intent(inout)  :: domain_vars

        character(len=kMAX_NAME_LENGTH), allocatable :: vars_to_read(:)
        integer,                         allocatable :: var_dimensions(:)


        ! the parameters option type can't contain allocatable arrays because it is a coarray
        ! so we need to allocate the vars_to_read and var_dimensions outside of the options type
        call setup_variable_lists(options%forcing%vars_to_read, options%forcing%dim_list, vars_to_read, var_dimensions)

        ! Read through forcing variable names stored in "options"
        ! needs to read each one to find the grid information for it
        ! then create grid and initialize a variable...
        ! also need to explicitly save lat and lon data
        ! if (STD_OUT_PE) then
            call this%init_local(options,                           &
                                 options%forcing%boundary_files, &
                                 vars_to_read, var_dimensions,      &
                                 options%general%start_time,     &
                                 options%forcing%latvar,         &
                                 options%forcing%lonvar,         &
                                 options%forcing%zvar,           &
                                 options%forcing%time_var,       &
                                 options%forcing%pvar,           &
                                 domain_lat, domain_lon, domain_vars)

        ! endif
        ! call this%distribute_initial_conditions()


        call setup_boundary_geo(this, options%domain%longitude_system)

    end subroutine


    !>------------------------------------------------------------
    !! Set default component values
    !! Reads initial conditions from the forcing file
    !!
    !!------------------------------------------------------------
    module subroutine init_local(this, options, file_list, var_list, dim_list, start_time, &
                                 lat_var, lon_var, z_var, time_var, p_var, domain_lat, domain_lon, domain_vars)
        class(boundary_t),               intent(inout)  :: this
        type(options_t),                 intent(inout)  :: options
        character(len=kMAX_NAME_LENGTH), intent(in)     :: file_list(:)
        character(len=kMAX_NAME_LENGTH), intent(in)     :: var_list (:)
        integer,                         intent(in)     :: dim_list (:)
        type(Time_type),                 intent(in)     :: start_time
        character(len=kMAX_NAME_LENGTH), intent(in)     :: lat_var
        character(len=kMAX_NAME_LENGTH), intent(in)     :: lon_var
        character(len=kMAX_NAME_LENGTH), intent(in)     :: z_var
        character(len=kMAX_NAME_LENGTH), intent(in)     :: time_var
        character(len=kMAX_NAME_LENGTH), intent(in)     :: p_var
        real, dimension(:,:),            intent(in)     :: domain_lat
        real, dimension(:,:),            intent(in)     :: domain_lon
        type(var_dict_t),                intent(inout)  :: domain_vars

        type(variable_t)  :: test_variable
        real, allocatable :: temp_z(:,:,:), temp_z_trans(:,:,:), temp_lat(:,:), temp_lon(:,:)
        character(len=5)       :: img_str

        integer :: i, nx, ny, nz

        ! figure out while file and timestep contains the requested start_time
        call set_firstfile_firststep(this, start_time, file_list, time_var)
        ! call read_bc_time(this%current_time, file_list(this%curfile), time_var, this%curstep)

        !  read in latitude and longitude coordinate data
        call io_read(this%firstfile, lat_var, temp_lat, this%firststep)
        call io_read(this%firstfile, lon_var, temp_lon, this%firststep)

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
            call io_read(this%firstfile, z_var,   temp_z,   1)
            nx = size(temp_z,1)
            ny = size(temp_z,2)
            nz = size(temp_z,3)
            if (allocated(this%z)) deallocate(this%z)
            allocate(this%z((this%ite-this%ite+1),nz,(this%jte-this%jts+1)))
            allocate(temp_z_trans(1:nx,1:nz,1:ny))
            
            temp_z_trans(1:nx,1:nz,1:ny) = reshape(temp_z, shape=[nx,nz,ny], order=[1,3,2])
            this%z = temp_z_trans(this%its:this%ite,1:nz,this%jts:this%jte)
        else
            call io_read(this%firstfile, p_var,   temp_z,   this%firststep)
            nx = size(temp_z,1)
            ny = size(temp_z,2)
            nz = size(temp_z,3)

            if (allocated(this%z)) deallocate(this%z)
            allocate(this%z((this%ite-this%ite+1),nz,(this%jte-this%jts+1)))

        endif
        
        this%kts = 1
        this%kte = nz

        if (this%ite < this%its) write(*,*) 'image: ',PE_RANK_GLOBAL+1,'  its: ',this%its,'  ite: ',this%ite,'  jts: ',this%jts,'  jte: ',this%jte
        if (this%kte < this%kts) write(*,*) 'image: ',PE_RANK_GLOBAL+1,'  its: ',this%its,'  ite: ',this%ite,'  jts: ',this%jts,'  jte: ',this%jte
        if (this%jte < this%jts) write(*,*) 'image: ',PE_RANK_GLOBAL+1,'  its: ',this%its,'  ite: ',this%ite,'  jts: ',this%jts,'  jte: ',this%jte

        ! call assert(size(var_list) == size(dim_list), "list of variable dimensions must match list of variables")
        do i=1, size(var_list)
            call add_var_to_dict(this, var_list(i), dim_list(i), [(this%ite-this%its+1), nz, (this%jte-this%jts+1)])
        end do

        call domain_vars%reset_iterator()

        do while (domain_vars%has_more_elements())
            test_variable = domain_vars%next()
            call add_var_hi_to_dict(this%variables_hi, test_variable%forcing_var, test_variable%grid)
        enddo

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

        character(len=MAXFILELENGTH) :: filename
        integer          :: error, n

        this%firststep = find_timestep_in_filelist(file_list, time_var, time, this%firstfile, error)

        if (error==1) then
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

        if (this%ite < this%its .or. this%jte < this%jts) write(*,*) 'image: ',PE_RANK_GLOBAL+1,'  its: ',this%its,'  ite: ',this%ite,'  jts: ',this%jts,'  jte: ',this%jte
        if (this%ite < this%its .or. this%jte < this%jts) write(*,*) 'image: ',PE_RANK_GLOBAL+1,'  d_ims: ',d_ims,'  d_ime: ',d_ime,'  d_jms: ',d_jms,'  d_jme: ',d_jme
        if (this%ite < this%its .or. this%jte < this%jts) write(*,*) 'image: ',PE_RANK_GLOBAL+1,'  LLlat: ',LLlat,'  LLlon: ',LLlon,'  URlat: ',URlat,'  URlon: ',URlon
        if (this%ite < this%its .or. this%jte < this%jts) write(*,*) 'image: ',PE_RANK_GLOBAL+1,'  min_loc: ',minloc(LL_d),'  max_loc: ',minloc(UR_d)
        if (this%ite < this%its .or. this%jte < this%jts) write(*,*) 'image: ',PE_RANK_GLOBAL+1,'  min_lat: ',minval(domain_lat),'  max_lat: ',maxval(domain_lat)
        if (this%ite < this%its .or. this%jte < this%jts) write(*,*) 'image: ',PE_RANK_GLOBAL+1,'  lat_corners: ',lat_corners,'  lon_corners: ',lon_corners
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
    subroutine add_var_to_dict(this, var_name, ndims, dims)
        implicit none
        type(boundary_t), intent(inout) :: this
        character(len=*), intent(in)    :: var_name
        integer,          intent(in)    :: ndims
        integer,          intent(in)    :: dims(3)

        real, allocatable :: temp_2d_data(:,:)
        real, allocatable :: temp_3d_data(:,:,:)
        type(variable_t)  :: new_variable

        if (ndims==2) then
            call new_variable%initialize( [dims(1),dims(3)] )
            call this%variables%add_var(var_name, new_variable)

        elseif (ndims==3) then
            call new_variable%initialize( dims )
            call this%variables%add_var(var_name, new_variable)

        ! these variables are computed (e.g. pressure from height or height from pressure)
        elseif (ndims==-3) then
            call new_variable%initialize( dims )
            new_variable%computed = .True.

            call this%variables%add_var(var_name, new_variable)
        endif

    end subroutine
    
    subroutine add_var_hi_to_dict(var_dict, var_name, grid)
        implicit none
        type(var_dict_t), intent(inout) :: var_dict
        character(len=*), intent(in)    :: var_name
        type(grid_t),     intent(in)    :: grid

        type(variable_t)  :: new_variable

        call new_variable%initialize( grid, var_name)
        call var_dict%add_var(var_name, new_variable)

    end subroutine


    module subroutine interpolate_original_levels(this, options)
        implicit none
        class(boundary_t),   intent(inout)   :: this
        type(options_t),     intent(in)      :: options

        type(variable_t)        :: input_z, var
        type(interpolable_type) :: input_geo
        real, allocatable :: temp_3d(:,:,:)
        character(len=kMAX_NAME_LENGTH) :: name


        associate(list => this%variables)

        input_z = list%get_var(options%forcing%zvar)

        if (options%forcing%z_is_geopotential) then
            input_z%data_3d = input_z%data_3d / gravity
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
            var = list%next(name)

            if (var%three_d) then
                ! need to vinterp this dataset to the original vertical levels (if necessary)

                temp_3d = var%data_3d
                call vinterp(var%data_3d, temp_3d, input_geo%vert_lut)

            endif

        end do
        end associate

    end subroutine


    module subroutine update_computed_vars(this, options, update)
        implicit none
        class(boundary_t),   intent(inout)   :: this
        type(options_t),     intent(in)      :: options
        logical,             intent(in),    optional :: update

        integer           :: err
        type(variable_t)  :: var, pvar, zvar, tvar, qvar
        logical :: update_internal

        integer :: nx,ny,nz
        real, allocatable :: temp_z(:,:,:)
        character(len=kMAX_NAME_LENGTH) :: name

        update_internal = .False.
        if (present(update)) update_internal = update

        associate(list => this%variables)

        if (options%forcing%qv_is_relative_humidity) then
            call compute_mixing_ratio_from_rh(list, options)
        endif

        if (options%forcing%qv_is_spec_humidity) then
            call compute_mixing_ratio_from_sh(list, options)
        endif

        ! because z is not updated over time, we don't want to reapply this every time, only in the initialization
        if (.not. update_internal) then
            if (options%forcing%z_is_geopotential) then
                this%z = this%z / gravity
            endif

            if (options%forcing%z_is_on_interface) then
                call interpolate_in_z(this%z)
            endif
        endif

        if (options%forcing%t_offset /= 0) then
            tvar = list%get_var(options%forcing%tvar)
            tvar%data_3d = tvar%data_3d + options%forcing%t_offset
        endif

        ! loop through the list of variables that need to be read in
        call list%reset_iterator()
        do while (list%has_more_elements())

            ! get the next variable in the structure
            var = list%next(name)

            if (var%computed) then

                if (trim(name) == trim(options%forcing%zvar)) then
                    call compute_z_update(this, list, options)
                endif

                if (trim(name) == trim(options%forcing%pvar)) then
                    call compute_p_update(this, list, options, var)
                endif

            endif
        end do

        if (.not.options%forcing%t_is_potential) then

            tvar = list%get_var(options%forcing%tvar)
            pvar = list%get_var(options%forcing%pvar)

            tvar%data_3d = tvar%data_3d / exner_function(pvar%data_3d)
        endif

        if (options%forcing%limit_rh) call limit_rh(list, options)
        !limit sea surface temperature to be >= 273.15 (i.e. cannot freeze)
        call limit_2d_var(list, options%forcing%sst_var, min_val=273.15)

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

        tvar = list%get_var(options%forcing%tvar)
        pvar = list%get_var(options%forcing%pvar)
        qvar = list%get_var(options%forcing%qvvar)

        if (maxval(qvar%data_3d) > 2) then
            qvar%data_3d = qvar%data_3d/100.0
        endif

        if (options%forcing%t_is_potential) then
            qvar%data_3d = rh_to_mr(qvar%data_3d, tvar%data_3d * exner_function(pvar%data_3d), pvar%data_3d)
        else
            qvar%data_3d = rh_to_mr(qvar%data_3d, tvar%data_3d, pvar%data_3d)
        endif


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

    end subroutine limit_2d_var

    subroutine compute_mixing_ratio_from_sh(list, options)
        implicit none
        type(var_dict_t),   intent(inout)   :: list
        type(options_t),    intent(in)      :: options

        integer           :: err
        type(variable_t)  :: qvar

        qvar = list%get_var(options%forcing%qvvar)

        qvar%data_3d = qvar%data_3d / (1 - qvar%data_3d)

    end subroutine compute_mixing_ratio_from_sh


    subroutine compute_z_update(this, list, options)
        implicit none
        class(boundary_t),  intent(inout)   :: this
        type(var_dict_t),   intent(inout)   :: list
        type(options_t),    intent(in)      :: options

        integer           :: err
        type(variable_t)  :: var, pvar, zvar, tvar, qvar
        real, allocatable, target :: real_t(:,:,:)

        real, pointer :: t(:,:,:)

        qvar = list%get_var(options%forcing%qvvar)
        tvar = list%get_var(options%forcing%tvar)
        zvar = list%get_var(options%forcing%hgtvar)
        var = list%get_var(options%forcing%pvar, err)

        pvar = list%get_var(options%forcing%pslvar, err)

        if (options%forcing%t_is_potential) then
            ! stop "Need real air temperature to compute height"
            allocate(real_t, mold=tvar%data_3d)
            real_t = exner_function(var%data_3d) * tvar%data_3d
            t => real_t
        else
            t => tvar%data_3d
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

    end subroutine compute_z_update

    subroutine compute_p_update(this, list, options, pressure_var)
        implicit none
        class(boundary_t),  intent(inout)   :: this
        type(var_dict_t),   intent(inout)   :: list
        type(options_t),    intent(in)      :: options
        type(variable_t),   intent(inout)   :: pressure_var

        integer           :: err
        type(variable_t)  :: pvar, zvar, tvar, qvar

        if (options%forcing%t_is_potential) stop "Need real air temperature to compute pressure"

        qvar = list%get_var(options%forcing%qvvar)
        tvar = list%get_var(options%forcing%tvar)
        zvar = list%get_var(options%forcing%hgtvar)

        pvar = list%get_var(options%forcing%pslvar, err)

        if (err == 0) then
            call compute_3d_p(pressure_var%data_3d, pvar%data_2d, this%z, tvar%data_3d, qvar%data_3d, zvar%data_2d)

        else
            pvar = list%get_var(options%forcing%psvar, err)

            if (err == 0) then
                call compute_3d_p(pressure_var%data_3d, pvar%data_2d, this%z, tvar%data_3d, qvar%data_3d)
            else
                write(*,*) "ERROR reading surface pressure or sea level pressure, variables not found"
                error stop
            endif
        endif


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
        integer,                         intent(in)                 :: master_dim_list(:)
        character(len=kMAX_NAME_LENGTH), intent(inout), allocatable :: vars_to_read(:)
        integer,                         intent(inout), allocatable :: var_dimensions(:)

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
                var_dimensions(curvar) = master_dim_list(i)
                ! if (STD_OUT_PE) print *, "in variable list: ", vars_to_read(curvar)
                curvar = curvar + 1
            endif
        enddo
    end subroutine



    module subroutine update_delta_fields(this, dt)
        implicit none
        class(boundary_t),    intent(inout) :: this
        type(time_delta_t), intent(in)    :: dt

        ! temporary to hold the variable to be interpolated to
        type(variable_t) :: var_to_update
        
        ! make sure the dictionary is reset to point to the first variable
        call this%variables_hi%reset_iterator()

        ! Now iterate through the dictionary as long as there are more elements present
        do while (this%variables_hi%has_more_elements())
            ! get the next variable
            var_to_update = this%variables_hi%next()

            if (var_to_update%two_d) then
                var_to_update%dqdt_2d = (var_to_update%dqdt_2d - var_to_update%data_2d) / dt%seconds()
            else if (var_to_update%three_d) then
                var_to_update%dqdt_3d = (var_to_update%dqdt_3d - var_to_update%data_3d) / dt%seconds()
            endif

        enddo

        ! w has to be handled separately because it is the only variable that can be updated using the delta fields but is not
        ! actually read from disk. Note that if we move to balancing winds every timestep, then it doesn't matter.
        !var_to_update = this%w%meta_data
        !var_to_update%dqdt_3d = (var_to_update%dqdt_3d - var_to_update%data_3d) / dt%seconds()
    end subroutine update_delta_fields

end submodule
