module output_metadata

    use icar_constants
    use meta_data_interface,    only : attribute_t, meta_data_t
    use string,       only : to_lower
    use iso_fortran_env, only : output_unit
    use options_interface, only : options_t
    implicit none

    ! type(variable_t), allocatable, target :: var_meta(:)

    !>------------------------------------------------------------
    !! Generic interface to the netcdf read routines
    !!------------------------------------------------------------
    ! interface get_metadata
    !     module procedure get_metadata_nod
    ! end interface


contains

    !>------------------------------------------------------------
    !! Get generic metadata for a no-data object
    !!------------------------------------------------------------
    ! function get_metadata_nod(var_idx) result(meta_data)
    !     implicit none
    !     integer, intent(in) :: var_idx
    !     type(variable_t), pointer :: meta_data

    !     else if (var_idx > kMAX_STORAGE_VARS) then
    !         write(*,*) "Invalid variable metadata requested, metadata index was: ", var_idx
    !         stop
    !     endif

    !     ! initialize the module level var_meta data
    !     if (.not.allocated(var_meta)) then    
    !         call init_var_meta()
    !     endif

    !     ! get the var_idx index into var_meta to return
    !     meta_data => var_meta(var_idx)

    ! end function get_metadata_nod

    ! !>------------------------------------------------------------
    ! !! Get generic metadata for a variable
    ! !!
    ! !! Sets the internal data pointer to point to the input data provided
    ! !!------------------------------------------------------------
    ! function get_metadata_var(var_idx, input_var) result(meta_data)
    !     implicit none
    !     integer, intent(in)          :: var_idx
    !     type(variable_t),intent(in)  :: input_var

    !     type(variable_t) :: meta_data       ! function result
    !     !integer          :: local_shape(2)  ! store the shape of the input data array

    !     else if (var_idx>kMAX_STORAGE_VARS) then
    !         stop "Invalid variable metadata requested"
    !     endif

    !     if (.not.allocated(var_meta)) call init_var_meta()

    !     meta_data = var_meta(var_idx)

    !     meta_data%two_d     = input_var_meta%two_d
    !     meta_data%three_d   = input_var_meta%three_d
    !     meta_data%grid      = input_var_meta%grid

    !     if(meta_data%two_d) then
    !         if (allocated(input_var_meta%dim_len)) allocate(meta_data%dim_len,source=input_var_meta%dim_len)
    !         if (allocated(input_var_meta%global_dim_len)) allocate(meta_data%global_dim_len,source=input_var_meta%global_dim_len)

    !         meta_data%data_2d => input_var_meta%data_2d
    !     else
    !         if (allocated(input_var_meta%dim_len)) then
    !             allocate(meta_data%dim_len(3))
    !             meta_data%dim_len(1) = input_var_meta%dim_len(1)
    !             meta_data%dim_len(2) = input_var_meta%dim_len(3)
    !             meta_data%dim_len(3) = input_var_meta%dim_len(2)
    !         endif
    !         if (allocated(input_var_meta%global_dim_len)) then
    !             allocate(meta_data%global_dim_len(3))
    !             meta_data%global_dim_len(1) = input_var_meta%global_dim_len(1)
    !             meta_data%global_dim_len(2) = input_var_meta%global_dim_len(3)
    !             meta_data%global_dim_len(3) = input_var_meta%global_dim_len(2)
    !         endif
    !         meta_data%data_3d => input_var_meta%data_3d
    !     endif

    ! end function get_metadata_var

    ! !>------------------------------------------------------------
    ! !! Get generic metadata for a two-dimensional variable
    ! !!
    ! !! Sets the internal data pointer to point to the input data provided
    ! !!------------------------------------------------------------
    ! function get_metadata_2d(var_idx, input_data) result(meta_data)
    !     implicit none
    !     integer, intent(in)          :: var_idx
    !     real,    intent(in), pointer :: input_data(:,:)

    !     type(variable_t) :: meta_data       ! function result
    !     integer          :: local_shape(2)  ! store the shape of the input data array

    !     else if (var_idx>kMAX_STORAGE_VARS) then
    !         stop "Invalid variable metadata requested"
    !     endif

    !     if (.not.allocated(var_meta)) call init_var_meta()

    !     meta_data = var_meta(var_idx)
    !     meta_data%dtype=kREAL

    !     if (associated(input_data)) then
    !         meta_data%data_2d   => input_data
    !         meta_data%two_d     = .True.
    !         meta_data%three_d   = .False.
    !         local_shape(1) = size(input_data, 1)
    !         local_shape(2) = size(input_data, 2)
    !         ! for some reason if shape(input_data) is passed as source, then the dim_len bounds are (0:1) instead of 1:2
    !         allocate(meta_data%dim_len, source=local_shape)
    !     endif

    ! end function get_metadata_2d

    ! function get_metadata_2dd(var_idx, input_data) result(meta_data)
    !     implicit none
    !     integer, intent(in)          :: var_idx
    !     double precision,    intent(in), pointer :: input_data(:,:)

    !     type(variable_t) :: meta_data       ! function result
    !     integer          :: local_shape(2)  ! store the shape of the input data array

    !     else if (var_idx>kMAX_STORAGE_VARS) then
    !         stop "Invalid variable metadata requested"
    !     endif

    !     if (.not.allocated(var_meta)) call init_var_meta()

    !     meta_data = var_meta(var_idx)
    !     meta_data%dtype=kDOUBLE

    !     if (associated(input_data)) then
    !         meta_data%data_2dd  => input_data
    !         meta_data%two_d     = .True.
    !         meta_data%three_d   = .False.
    !         local_shape(1) = size(input_data, 1)
    !         local_shape(2) = size(input_data, 2)
    !         ! for some reason if shape(input_data) is passed as source, then the dim_len bounds are (0:1) instead of 1:2
    !         allocate(meta_data%dim_len, source=local_shape)
    !     endif

    ! end function get_metadata_2dd

    ! !>------------------------------------------------------------
    ! !! Get generic metadata for a three-dimensional variable
    ! !!
    ! !! Sets the internal data pointer to point to the input data provided
    ! !!------------------------------------------------------------
    ! function get_metadata_3d(var_idx, input_data) result(meta_data)
    !     implicit none
    !     integer, intent(in)          :: var_idx
    !     real,    intent(in), pointer :: input_data(:,:,:)

    !     type(variable_t) :: meta_data       ! function result
    !     integer          :: local_shape(3)  ! store the shape of the input data array

    !     else if (var_idx>kMAX_STORAGE_VARS) then
    !         stop "Invalid variable metadata requested"
    !     endif

    !     ! initialize the module level constant data structure
    !     if (.not.allocated(var_meta)) call init_var_meta()

    !     meta_data = var_meta(var_idx)
    !     meta_data%dtype=kREAL

    !     if (associated(input_data)) then
    !         meta_data%data_3d   => input_data
    !         meta_data%two_d     = .False.
    !         meta_data%three_d   = .True.
    !         local_shape(1) = size(input_data, 1)
    !         local_shape(2) = size(input_data, 3)
    !         local_shape(3) = size(input_data, 2)
    !         ! for some reason if shape(input_data) is passed as source, then the dim_len bounds are (0:1) instead of 1:2
    !         allocate(meta_data%dim_len, source=local_shape)
    !     endif

    ! end function get_metadata_3d


    !>------------------------------------------------------------
    !! Get metadata variable name associated with a given index
    !!
    !! Sets the internal data pointer to point to the input data provided
    !!------------------------------------------------------------
    function get_varname(var_idx) result(name)
        implicit none
        integer, intent(in)             :: var_idx
        character(len=kMAX_NAME_LENGTH) :: name       ! function result
        type(meta_data_t)               :: var_meta

        if (var_idx>kMAX_STORAGE_VARS) then
            stop "Invalid variable metadata requested"
        endif
        ! initialize the module level constant data structure
        !if (.not.allocated(var_meta)) call init_var_meta()

        var_meta = get_varmeta(var_idx)
        name = var_meta%name

    end function get_varname

    !>------------------------------------------------------------
    !! Get metadata variable name associated with a given index
    !!
    !! Sets the internal data pointer to point to the input data provided
    !!------------------------------------------------------------
    function get_varindx(var_name) result(indx)
        implicit none
        character(len=*), intent(in) :: var_name       
        integer                                     :: indx  ! function result
        type(meta_data_t)               :: var_meta

        ! initialize the module level constant data structure
        !if (.not.allocated(var_meta)) call init_var_meta()

        do indx = 1, kMAX_STORAGE_VARS
            var_meta = get_varmeta(indx)
            if (var_meta%name == trim(var_name)) return
        enddo

        if (indx>kMAX_STORAGE_VARS) then
            if (STD_OUT_PE) write(*,*) "Invalid variable metadata requested, no metadata entry for var_name: ", trim(var_name)
        endif

    end function get_varindx

    function get_var_matchinfo(var_meta,request_string) result(match_info)
        implicit none
        type(meta_data_t), intent(in) :: var_meta
        character(len=*), intent(in) :: request_string

        logical     :: match_info
        integer :: i

        ! default is to assume no match
        match_info = .False.

        ! initialize the module level constant data structure
        !if (.not.allocated(var_meta)) call init_var_meta()
        
        ! check if any part of request_string is in the description or name of the variable
        if (index(to_lower(var_meta%name), trim(to_lower(request_string))) > 0) then
            match_info = .True.
        else
            ! loop through each element of the attribute array
            do i=1,var_meta%n_attrs
                ! ignore the "units" and "coordinates" attributes
                if (index(var_meta%attributes(i)%name, "units") > 0) cycle
                if (index(var_meta%attributes(i)%name, "coordinates") > 0) cycle

                if (index(to_lower(var_meta%attributes(i)%value), trim(to_lower(request_string))) > 0) then
                    match_info = .True.
                    exit
                endif
            enddo
        endif
    end function get_var_matchinfo

    subroutine list_output_vars(keywords)
        implicit none
        character(len=*), intent(in) :: keywords(:)

        type(meta_data_t) :: var_meta
        integer :: i, j, k, match_count
        character(len=kMAX_NAME_LENGTH) :: dimensions

        ! initialize the module level constant data structure
        !if (.not.allocated(var_meta)) call init_var_meta()

        if (STD_OUT_PE) write(*,*) "-------------------------------------------------"
        if (STD_OUT_PE) write(*,"(A)", ADVANCE='NO') "Output variables matching keywords: "
        do i = 1,size(keywords)
            if (STD_OUT_PE) write(*,"(A,A)", ADVANCE='NO') trim(keywords(i)), " "
        enddo
        if (STD_OUT_PE) write(*,"(/ A)") "-------------------------------------------------"
        if (STD_OUT_PE) write(*,*)


        do i=1,kMAX_STORAGE_VARS
            match_count = 0
            var_meta = get_varmeta(i)
            !check if name is not an empty string
            if (len_trim(var_meta%name) == 0) cycle
            do k = 1,size(keywords)
                if (get_var_matchinfo(var_meta, keywords(k))) then
                    match_count = match_count + 1
                endif
            enddo
            if (match_count == size(keywords)) then
                dimensions = ""
                do j = 1,size(var_meta%dimensions)
                    dimensions = trim(dimensions)//trim(var_meta%dimensions(j))//", "
                enddo
                ! remove the trailing space
                dimensions = trim(dimensions)
                !remove trailing comma
                dimensions = dimensions(1:len_trim(dimensions)-1)

                if (STD_OUT_PE) write(*,*) "Name: ",trim(var_meta%name)
                if (STD_OUT_PE) write(*,*) 
                if (STD_OUT_PE) write(*,"(A25,A)") "   namelist name: ",trim(var_meta%name)
                if (STD_OUT_PE) write(*,"(A25,A1,A,A1)") "   dimensions: ","[",trim(dimensions),"]"
                do j=1,var_meta%n_attrs
                    if (STD_OUT_PE) write(*,"(A25,A)") (trim(var_meta%attributes(j)%name)//": "),trim(var_meta%attributes(j)%value)
                enddo
                if (STD_OUT_PE) write(*,*) "-------------------------------------------------"
            endif
        enddo

    end subroutine list_output_vars

    !>------------------------------------------------------------
    !! Initialize the module level master data structure
    !!
    !!------------------------------------------------------------
    function get_varmeta(var_idx, opt, forcing_var, force_boundaries) result(var_meta)
        implicit none
        integer, intent(in) :: var_idx
        type(options_t), intent(in), optional :: opt
        logical, intent(out), optional :: forcing_var, force_boundaries

        type(meta_data_t) :: var_meta
        integer :: i, j

        if (kVARS%last_var/=kMAX_STORAGE_VARS) then
            stop "ERROR: variable indicies not correctly initialized"
        endif

        if (present(forcing_var)) then
            if (present(opt)) then
                forcing_var = .False.
            else
                stop 'Error in function get_varmeta: when forcing_var is passed in, the argument opt must also be passed'
            endif
        endif
        if (present(force_boundaries)) force_boundaries = .True.

        ! allocate the actual data structure to be used
        allocate(var_meta%dim_len,source=(/ 0, 0, 0/))

        !>------------------------------------------------------------
        !!  U  East West Winds
        !!------------------------------------------------------------
        if (var_idx==kVARS%u) then
            var_meta%name        = "u"
            var_meta%dimensions  = three_d_u_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "grid_eastward_wind"),              &
                               attribute_t("long_name",     "Grid relative eastward wind"),     &
                               attribute_t("units",         "m s-1"),                           &
                               attribute_t("coordinates",   "u_lat u_lon")]
            if (present(force_boundaries)) force_boundaries = .False.
            if (present(opt) .and. present(forcing_var)) forcing_var = (opt%forcing%uvar /= "")
        !>------------------------------------------------------------
        !!  V  North South Winds
        !!------------------------------------------------------------
        else if (var_idx==kVARS%v) then
            var_meta%name        = "v"
            var_meta%dimensions  = three_d_v_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "grid_northward_wind"),             &
                               attribute_t("long_name",     "Grid relative northward wind"),    &
                               attribute_t("units",         "m s-1"),                           &
                               attribute_t("coordinates",   "v_lat v_lon")]
            if (present(force_boundaries)) force_boundaries = .False.
            if (present(opt) .and. present(forcing_var)) forcing_var = (opt%forcing%vvar /= "")

        !>------------------------------------------------------------
        !!  W  Vertical Winds
        !!------------------------------------------------------------
        else if (var_idx==kVARS%w) then
            var_meta%name        = "w_grid"
            var_meta%dimensions  = three_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "grid_upward_air_velocity"),    &
                               attribute_t("long_name",     "Vertical wind"),                   &
                               attribute_t("description",   "Vertical wind relative to the grid"),&
                               attribute_t("units",         "m s-1"),                           &
                               attribute_t("coordinates",   "lat lon")]
        

        else if (var_idx==kVARS%w_real) then
            var_meta%name        = "w"
            var_meta%dimensions  = three_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "upward_air_velocity"),             &
                               attribute_t("long_name",     "Vertical wind"),                   &
                               attribute_t("description",   "Vertical wind including u/v"),     &
                               attribute_t("units",         "m s-1"),                           &
                               attribute_t("coordinates",   "lat lon")]
            if (present(force_boundaries)) force_boundaries = .False.
            if (present(opt) .and. present(forcing_var)) forcing_var = (opt%forcing%wvar /= "")

        !>------------------------------------------------------------
        !!  Brunt Vaisala frequency (squared)
        !!------------------------------------------------------------
        else if (var_idx==kVARS%nsquared) then
            var_meta%name        = "nsquared"
            var_meta%dimensions  = three_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "square_of_brunt_vaisala_frequency_in_air"),&
                               attribute_t("long_name",     "Burnt Vaisala frequency squared"), &
                               attribute_t("units",         "s-2"),                             &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Air Pressure
        !!------------------------------------------------------------
        else if (var_idx==kVARS%pressure) then
            var_meta%name        = "pressure"
            var_meta%dimensions  = three_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "air_pressure"),                    &
                               attribute_t("long_name",     "Pressure"),                        &
                               attribute_t("units",         "Pa"),                              &
                               attribute_t("coordinates",   "lat lon")]
            if (present(force_boundaries)) force_boundaries = .False.
            if (present(opt) .and. present(forcing_var)) forcing_var = (opt%forcing%pvar /= "")

        !>------------------------------------------------------------
        !!  Air Pressure on interfaces between mass levels
        !!------------------------------------------------------------
        else if (var_idx==kVARS%pressure_interface) then
            var_meta%name        = "pressure_i"
            var_meta%dimensions  = three_d_t_interface_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "air_pressure"),                    &
                               attribute_t("long_name",     "Pressure"),                        &
                               attribute_t("units",         "Pa"),                              &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Air Pressure on interfaces between mass levels
        !!------------------------------------------------------------
        else if (var_idx==kVARS%surface_pressure) then
            var_meta%name        = "psfc"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "surface_air_pressure"),            &
                               attribute_t("long_name",     "Surface Pressure"),                &
                               attribute_t("units",         "Pa"),                              &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Potential Air Temperature
        !!------------------------------------------------------------
        else if (var_idx==kVARS%potential_temperature) then
            var_meta%name        = "potential_temperature"
            var_meta%dimensions  = three_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "air_potential_temperature"),       &
                               attribute_t("long_name",     "Potential Temperature"),           &
                               attribute_t("units",         "K"),                               &
                               attribute_t("coordinates",   "lat lon")]
            if (present(opt) .and. present(forcing_var)) forcing_var = (opt%forcing%tvar /= "")

        !>------------------------------------------------------------
        !!  Real Air Temperature
        !!------------------------------------------------------------
        else if (var_idx==kVARS%temperature) then
            var_meta%name        = "temperature"
            var_meta%dimensions  = three_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "air_temperature"),                 &
                               attribute_t("long_name",     "Temperature"),                     &
                               attribute_t("units",         "K"),                               &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Water Vapor Mixing Ratio
        !!------------------------------------------------------------
        else if (var_idx==kVARS%water_vapor) then
            var_meta%name        = "qv"
            var_meta%dimensions  = three_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "mass_fraction_of_water_vapor_in_air"), &
                               attribute_t("long_name",     "Water Vapor Mixing Ratio"),            &
                               attribute_t("units",         "kg kg-1"),                             &
                               attribute_t("coordinates",   "lat lon")]
            if (present(opt) .and. present(forcing_var)) forcing_var = (opt%forcing%qvvar /= "")

        !>------------------------------------------------------------
        !!  Cloud water (liquid) mixing ratio
        !!------------------------------------------------------------
        else if (var_idx==kVARS%cloud_water_mass) then
            var_meta%name        = "qc"
            var_meta%dimensions  = three_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "cloud_liquid_water_mixing_ratio"),     &
                               attribute_t("units",         "kg kg-1"),                             &
                               attribute_t("coordinates",   "lat lon")]
            if (present(opt) .and. present(forcing_var)) forcing_var = (opt%forcing%qcvar /= "")

        !>------------------------------------------------------------
        !!  Cloud water (liquid) number concentration
        !!------------------------------------------------------------
        else if (var_idx==kVARS%cloud_number) then
            var_meta%name        = "nc"
            var_meta%dimensions  = three_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "number_concentration_of_cloud_droplets_in_air"), &
                               attribute_t("units",         "cm-3"),                                              &
                               attribute_t("coordinates",   "lat lon")]
            if (present(opt) .and. present(forcing_var)) forcing_var = (opt%forcing%qncvar /= "")
        !>------------------------------------------------------------
        !!  Cloud ice mixing ratio
        !!------------------------------------------------------------
        else if (var_idx==kVARS%ice_mass) then
            var_meta%name        = "qi"
            var_meta%dimensions  = three_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "cloud_ice_mixing_ratio"),              &
                               attribute_t("units",         "kg kg-1"),                             &
                               attribute_t("coordinates",   "lat lon")]
            if (present(opt) .and. present(forcing_var)) forcing_var = (opt%forcing%qivar /= "")
        !>------------------------------------------------------------
        !!  Cloud ice number concentration
        !!------------------------------------------------------------
        else if (var_idx==kVARS%ice_number) then
            var_meta%name        = "ni"
            var_meta%dimensions  = three_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "number_concentration_of_ice_crystals_in_air"), &
                               attribute_t("units",         "cm-3"),                                            &
                               attribute_t("coordinates",   "lat lon")]
            if (present(opt) .and. present(forcing_var)) forcing_var = (opt%forcing%qnivar /= "")        
        !>------------------------------------------------------------
        !!  Rain water mixing ratio
        !!------------------------------------------------------------
        else if (var_idx==kVARS%rain_mass) then
            var_meta%name        = "qr"
            var_meta%dimensions  = three_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "mass_fraction_of_rain_in_air"),        &
                               attribute_t("units",         "kg kg-1"),                             &
                               attribute_t("coordinates",   "lat lon")]
            if (present(opt) .and. present(forcing_var)) forcing_var = (opt%forcing%qrvar /= "")
        !>------------------------------------------------------------
        !!  Rain water number concentration
        !!------------------------------------------------------------
        else if (var_idx==kVARS%rain_number) then
            var_meta%name        = "nr"
            var_meta%dimensions  = three_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "number_concentration_of_rain_particles_in_air"), &
                               attribute_t("units",         "cm-3"),                                              &
                               attribute_t("coordinates",   "lat lon")]
            if (present(opt) .and. present(forcing_var)) forcing_var = (opt%forcing%qnrvar /= "")
        !>------------------------------------------------------------
        !!  Snow in air mixing ratio
        !!------------------------------------------------------------
        else if (var_idx==kVARS%snow_mass) then
            var_meta%name        = "qs"
            var_meta%dimensions  = three_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "mass_fraction_of_snow_in_air"),        &
                               attribute_t("units",         "kg kg-1"),                             &
                               attribute_t("coordinates",   "lat lon")]
            if (present(opt) .and. present(forcing_var)) forcing_var = (opt%forcing%qsvar /= "")
        !>------------------------------------------------------------
        !!  Snow in air number concentration
        !!------------------------------------------------------------
        else if (var_idx==kVARS%snow_number) then
            var_meta%name        = "ns"
            var_meta%dimensions  = three_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "number_concentration_of_snow_particles_in_air"), &
                               attribute_t("units",         "cm-3"),                                              &
                               attribute_t("coordinates",   "lat lon")]
            if (present(opt) .and. present(forcing_var)) forcing_var = (opt%forcing%qnsvar /= "")
        !>------------------------------------------------------------
        !!  Graupel mixing ratio
        !!------------------------------------------------------------
        else if (var_idx==kVARS%graupel_mass) then
            var_meta%name        = "qg"
            var_meta%dimensions  = three_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "mass_fraction_of_graupel_in_air"),     &
                               attribute_t("units",         "kg kg-1"),                             &
                               attribute_t("coordinates",   "lat lon")]
            if (present(opt) .and. present(forcing_var)) forcing_var = (opt%forcing%qgvar /= "")
        !>------------------------------------------------------------
        !!  Graupel number concentration
        !!------------------------------------------------------------
        else if (var_idx==kVARS%graupel_number) then
            var_meta%name        = "ng"
            var_meta%dimensions  = three_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "number_concentration_of_graupel_particles_in_air"), &
                               attribute_t("units",         "cm-3"),                                                 &
                               attribute_t("coordinates",   "lat lon")]
            if (present(opt) .and. present(forcing_var)) forcing_var = (opt%forcing%qngvar /= "")
        !>------------------------------------------------------------
        !!  planar-nucleated a^2c mixing ratio
        !!------------------------------------------------------------
        else if (var_idx==kVARS%ice1_a) then
            var_meta%name        = "ice1_a"
            var_meta%dimensions  = three_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "planar-nucleated a^2c mixing ratio"), &
                               attribute_t("units",         "m^3 kg^-1"),                                                 &
                               attribute_t("coordinates",   "lat lon")]
            if (present(opt) .and. present(forcing_var)) forcing_var = (opt%forcing%i1avar /= "")
        !>------------------------------------------------------------
        !!  planar-nucleated c^2a mixing ratio
        !!------------------------------------------------------------
        else if (var_idx==kVARS%ice1_c) then
            var_meta%name        = "ice1_c"
            var_meta%dimensions  = three_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "planar-nucleated c^2a mixing ratio"), &
                               attribute_t("units",         "m^3 kg^-1"),                                                 &
                               attribute_t("coordinates",   "lat lon")]
            if (present(opt) .and. present(forcing_var)) forcing_var = (opt%forcing%i1cvar /= "")
        !>------------------------------------------------------------
        !!  columnar-nucleated mixing ratio
        !!------------------------------------------------------------
        else if (var_idx==kVARS%ice2_mass) then
            var_meta%name        = "ice2_mass"
            var_meta%dimensions  = three_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "columnar-nucleated mixing ratio"), &
                               attribute_t("units",         "kg kg^-1"),                                                 &
                               attribute_t("coordinates",   "lat lon")]
            if (present(opt) .and. present(forcing_var)) forcing_var = (opt%forcing%i2mvar /= "")
        !>------------------------------------------------------------
        !!  columnar-nucleated number mixing ratio
        !!------------------------------------------------------------
        else if (var_idx==kVARS%ice2_number) then
            var_meta%name        = "ice2_number"
            var_meta%dimensions  = three_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "columnar-nucleated number mixing ratio"), &
                               attribute_t("units",         "kg^-1"),                                                 &
                               attribute_t("coordinates",   "lat lon")]
            if (present(opt) .and. present(forcing_var)) forcing_var = (opt%forcing%i2nvar /= "")
        !>------------------------------------------------------------
        !!  columnar-nucleated a^2c mixing ratio
        !!------------------------------------------------------------
        else if (var_idx==kVARS%ice2_a) then
            var_meta%name        = "ice2_a"
            var_meta%dimensions  = three_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "columnar-nucleated a^2c mixing ratio"), &
                               attribute_t("units",         "m^3 kg^-1"),                                                 &
                               attribute_t("coordinates",   "lat lon")]
            if (present(opt) .and. present(forcing_var)) forcing_var = (opt%forcing%i2avar /= "")
        !>------------------------------------------------------------
        !!  columnar-nucleated c^2a mixing ratio
        !!------------------------------------------------------------
        else if (var_idx==kVARS%ice2_c) then
            var_meta%name        = "ice2_c"
            var_meta%dimensions  = three_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "columnar-nucleated c^2a mixing ratio"), &
                               attribute_t("units",         "m^3 kg^-1"),                                                 &
                               attribute_t("coordinates",   "lat lon")]
            if (present(opt) .and. present(forcing_var)) forcing_var = (opt%forcing%i2cvar /= "")
        !>------------------------------------------------------------
        !!  aggregate mixing ratio
        !!------------------------------------------------------------
        else if (var_idx==kVARS%ice3_mass) then
            var_meta%name        = "ice3_mass"
            var_meta%dimensions  = three_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "aggregate mixing ratio"), &
                               attribute_t("units",         "kg kg^-1"),                                                 &
                               attribute_t("coordinates",   "lat lon")]
            if (present(opt) .and. present(forcing_var)) forcing_var = (opt%forcing%i3mvar /= "")
        !>------------------------------------------------------------
        !!  aggregate number mixing ratio
        !!------------------------------------------------------------
        else if (var_idx==kVARS%ice3_number) then
            var_meta%name        = "ice3_number"
            var_meta%dimensions  = three_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "aggregate number mixing ratio"), &
                               attribute_t("units",         "kg^-1"),                                                 &
                               attribute_t("coordinates",   "lat lon")]
            if (present(opt) .and. present(forcing_var)) forcing_var = (opt%forcing%i3nvar /= "")
        !>------------------------------------------------------------
        !!  aggregate a^2c mixing ratio
        !!------------------------------------------------------------
        else if (var_idx==kVARS%ice3_a) then
            var_meta%name        = "ice3_a"
            var_meta%dimensions  = three_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "aggregate a^2c mixing ratio"), &
                               attribute_t("units",         "m^3 kg^-1"),                                                 &
                               attribute_t("coordinates",   "lat lon")]
            if (present(opt) .and. present(forcing_var)) forcing_var = (opt%forcing%i3avar /= "")
        !>------------------------------------------------------------
        !!  aggregate c^2a mixing ratio
        !!------------------------------------------------------------
        else if (var_idx==kVARS%ice3_c) then
            var_meta%name        = "ice3_c"
            var_meta%dimensions  = three_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "aggregate c^2a mixing ratio"), &
                               attribute_t("units",         "m^3 kg^-1"),                                                 &
                               attribute_t("coordinates",   "lat lon")]
            if (present(opt) .and. present(forcing_var)) forcing_var = (opt%forcing%i3cvar /= "")

        !>------------------------------------------------------------
        !!  Precipitation rate at the surface (requires tracking past precipitation amounts)
        !!------------------------------------------------------------
        ! else if (var_idx==kVARS%precip_rate)
        !     var_meta%name        = "precip_rate"
        !     var_meta%dimensions  = two_d_t_dimensions
        !     var_meta%unlimited_dim=.True.
        !     var_meta%attributes  = [attribute_t("standard_name",   "precipitation_flux"),      &
        !                        attribute_t("units",           "kg m-2 s-1"),              &
        !                        attribute_t("coordinates",     "lat lon")]
        ! 
        !>------------------------------------------------------------
        !!  Accumulated precipitation at the surface
        !!------------------------------------------------------------
        else if (var_idx==kVARS%precipitation) then
            var_meta%name        = "precipitation"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "precipitation_amount"),                &
                               attribute_t("units",         "kg m-2"),                              &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Accumulated convective precipitation at the surface
        !!------------------------------------------------------------
        else if (var_idx==kVARS%convective_precipitation) then
            var_meta%name        = "cu_precipitation"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "convective_precipitation_amount"),     &
                               attribute_t("units",         "kg m-2"),                              &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Accumulated snowfall at the surface
        !!------------------------------------------------------------
        else if (var_idx==kVARS%snowfall) then
            var_meta%name        = "snowfall"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "snowfall_amount"),                     &
                               attribute_t("units",         "kg m-2"),                              &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Accumulated Graupel at the surface
        !!------------------------------------------------------------
        else if (var_idx==kVARS%graupel) then
            var_meta%name        = "graupel"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "graupel_amount"),                      &
                               attribute_t("units",         "kg m-2"),                              &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Exner function
        !!------------------------------------------------------------
        else if (var_idx==kVARS%exner) then
            var_meta%name        = "exner"
            var_meta%dimensions  = three_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "exner_function_result"),           &
                               attribute_t("units",         "K K-1"),                               &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Air Density
        !!------------------------------------------------------------
        else if (var_idx==kVARS%density) then
            var_meta%name        = "density"
            var_meta%dimensions  = three_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "air_density"),                         &
                               attribute_t("units",         "kg m-3"),                              &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Vertical coordinate
        !!------------------------------------------------------------
        else if (var_idx==kVARS%z) then
            var_meta%name        = "z"
            var_meta%dimensions  = three_d_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "height_above_reference_ellipsoid"),    &
                               attribute_t("units",         "m"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Vertical coordinate on the interface between mass levels
        !!------------------------------------------------------------
        else if (var_idx==kVARS%z_interface) then
            var_meta%name        = "z_i"
            var_meta%dimensions  = three_d_interface_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "height_above_reference_ellipsoid"),    &
                               attribute_t("units",         "m"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        !>------------------------------------------------------------
        !!  Vertical coordinate on the interface between mass levels on global grid
        !!------------------------------------------------------------
        else if (var_idx==kVARS%global_z_interface) then
            var_meta%name        = "global_z_i"
            var_meta%dimensions  = three_d_neighbor_interface_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "height_above_reference_ellipsoid"),    &
                                attribute_t("units",         "m"),                                   &
                                attribute_t("coordinates",   "lat lon")]
                    
        !>------------------------------------------------------------
        !!  Vertical layer thickness
        !!------------------------------------------------------------
        else if (var_idx==kVARS%dz) then
            var_meta%name        = "dz"
            var_meta%dimensions  = three_d_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "layer_thickness between mass points"),                 &
                               attribute_t("units",         "m"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Thickness of layers between interfaces
        !!------------------------------------------------------------
        else if (var_idx==kVARS%dz_interface) then
            var_meta%name        = "dz_i"
            var_meta%dimensions  = three_d_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "layer_thickness between k half-levels"),                 &
                               attribute_t("units",         "m"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        !>------------------------------------------------------------
        !!  Thickness of layers between interfaces for global grid
        !!------------------------------------------------------------
        else if (var_idx==kVARS%global_dz_interface) then
            var_meta%name        = "global_dz_i"
            var_meta%dimensions  = three_d_neighbor_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "layer_thickness between k half-levels"),                 &
                                attribute_t("units",         "m"),                                   &
                                attribute_t("coordinates",   "lat lon")]
                    
        !>------------------------------------------------------------
        !!  Change in height in x direction
        !!------------------------------------------------------------
        else if (var_idx==kVARS%dzdx) then
            var_meta%name        = "dzdx"
            var_meta%dimensions  = three_d_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "dzdx of domain mesh"),                 &
                               attribute_t("units",         "-"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Change in height in y direction
        !!------------------------------------------------------------
        else if (var_idx==kVARS%dzdy) then
            var_meta%name        = "dzdy"
            var_meta%dimensions  = three_d_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "dzdy of domain mesh"),                 &
                               attribute_t("units",         "-"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        !>------------------------------------------------------------
        !!  Change in height in x direction on the staggered u grid
        !!------------------------------------------------------------
        else if (var_idx==kVARS%dzdx_u) then
            var_meta%name        = "dzdx_u"
            var_meta%dimensions  = three_d_u_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "dzdx of domain mesh on staggered u grid"),                 &
                                attribute_t("units",         "-"),                                   &
                                attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Change in height in y direction on the staggered v grid
        !!------------------------------------------------------------
        else if (var_idx==kVARS%dzdy_v) then
            var_meta%name        = "dzdy_v"
            var_meta%dimensions  = three_d_v_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "dzdy of domain mesh on staggered v grid"),                 &
                                attribute_t("units",         "-"),                                   &
                                attribute_t("coordinates",   "lat lon")]
        !>------------------------------------------------------------
        !!  High-frequency terrain component
        !!------------------------------------------------------------
        else if (var_idx==kVARS%h1) then
            var_meta%name        = "h1"
            var_meta%dimensions  = two_d_neighbor_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "High-frequency terrain component"),                 &
                                attribute_t("units",         "m"),                                   &
                                attribute_t("coordinates",   "lat lon")]
        !>------------------------------------------------------------
        !!  Low-frequency terrain component
        !!------------------------------------------------------------
        else if (var_idx==kVARS%h2) then
            var_meta%name        = "h2"
            var_meta%dimensions  = two_d_neighbor_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "Low-frequency terrain component"),                 &
                                attribute_t("units",         "m"),                                   &
                                attribute_t("coordinates",   "lat lon")]
        !>------------------------------------------------------------
        !!  High-frequency terrain component on the staggered u grid
        !!------------------------------------------------------------
        else if (var_idx==kVARS%h1_u) then
            var_meta%name        = "h1_u"
            var_meta%dimensions  = two_d_u_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "High-frequency terrain component on staggered u grid"),                 &
                                attribute_t("units",         "m"),                                   &
                                attribute_t("coordinates",   "lat lon")]
        !>------------------------------------------------------------
        !!  Low-frequency terrain component on the staggered u grid
        !!------------------------------------------------------------
        else if (var_idx==kVARS%h2_u) then
            var_meta%name        = "h2_u"
            var_meta%dimensions  = two_d_u_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "Low-frequency terrain component on staggered u grid"),                 &
                                attribute_t("units",         "m"),                                   &
                                attribute_t("coordinates",   "lat lon")]
        !>------------------------------------------------------------
        !!  High-frequency terrain component on staggered v grid
        !!------------------------------------------------------------
        else if (var_idx==kVARS%h1_v) then
            var_meta%name        = "h1_v"
            var_meta%dimensions  = two_d_v_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "High-frequency terrain component on staggered v grid"),                 &
                                attribute_t("units",         "m"),                                   &
                                attribute_t("coordinates",   "lat lon")]
        !>------------------------------------------------------------
        !!  Low-frequency terrain component on staggered v grid
        !!------------------------------------------------------------
        else if (var_idx==kVARS%h2_v) then
            var_meta%name        = "h2_v"
            var_meta%dimensions  = two_d_v_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "Low-frequency terrain component on staggered v grid"),                 &
                                attribute_t("units",         "m"),                                   &
                                attribute_t("coordinates",   "lat lon")]
        !>------------------------------------------------------------
        !!  The Jacobian of the z-coordinate transform
        !!------------------------------------------------------------
        else if (var_idx==kVARS%jacobian) then
            var_meta%name        = "jacobian"
            var_meta%dimensions  = three_d_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "Jacobian of the z-coordinate transform"),                 &
                                attribute_t("units",         "-"),                                   &
                                attribute_t("coordinates",   "lat lon")]
        !>------------------------------------------------------------
        !!  The Jacobian of the z-coordinate transform on staggered u grid
        !!------------------------------------------------------------
        else if (var_idx==kVARS%jacobian_u) then
            var_meta%name        = "jacobian_u"
            var_meta%dimensions  = three_d_u_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "Jacobian of the z-coordinate transform on staggered u grid"),                 &
                                attribute_t("units",         "-"),                                   &
                                attribute_t("coordinates",   "lat lon")]
        !>------------------------------------------------------------
        !!  The Jacobian of the z-coordinate transform on staggered v grid
        !!------------------------------------------------------------
        else if (var_idx==kVARS%jacobian_v) then
            var_meta%name        = "jacobian_v"
            var_meta%dimensions  = three_d_v_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "Jacobian of the z-coordinate transform on staggered v grid"),                 &
                                attribute_t("units",         "-"),                                   &
                                attribute_t("coordinates",   "lat lon")]
        !>------------------------------------------------------------
        !!  The Jacobian of the z-coordinate transform on k half-levels
        !!------------------------------------------------------------
        else if (var_idx==kVARS%jacobian_w) then
            var_meta%name        = "jacobian_w"
            var_meta%dimensions  = three_d_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "Jacobian of the z-coordinate transform on k half-levels"),                 &
                                attribute_t("units",         "-"),                                   &
                                attribute_t("coordinates",   "lat lon")]
        !>------------------------------------------------------------
        !!  sin of grid rotation from true north
        !!------------------------------------------------------------
        else if (var_idx==kVARS%sintheta) then
            var_meta%name        = "sintheta"
            var_meta%dimensions  = two_d_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "sin of grid rotation from true north"),                 &
                                attribute_t("units",         "-"),                                   &
                                attribute_t("coordinates",   "lat lon")]
        !>------------------------------------------------------------
        !!  cos of grid rotation from true north
        !!------------------------------------------------------------
        else if (var_idx==kVARS%costheta) then
            var_meta%name        = "costheta"
            var_meta%dimensions  = two_d_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "cos of grid rotation from true north"),                 &
                                attribute_t("units",         "-"),                                   &
                                attribute_t("coordinates",   "lat lon")]
        !>------------------------------------------------------------
        !!  Relaxation filter used for attenuating forcing of boundary forcing variables
        !!------------------------------------------------------------
        else if (var_idx==kVARS%relax_filter_2d) then
            var_meta%name        = "relax_filter_2d"
            var_meta%dimensions  = two_d_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "2D relaxation filter used for attenuating forcing of boundary forcing variables"),                 &
                                attribute_t("units",         "-"),                                   &
                                attribute_t("coordinates",   "lat lon")]
        !>------------------------------------------------------------
        !!  Relaxation filter used for attenuating forcing of boundary forcing variables
        !!------------------------------------------------------------
        else if (var_idx==kVARS%relax_filter_3d) then
            var_meta%name        = "relax_filter_3d"
            var_meta%dimensions  = three_d_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "3D relaxation filter used for attenuating forcing of boundary forcing variables"),                 &
                                attribute_t("units",         "-"),                                   &
                                attribute_t("coordinates",   "lat lon")]
        !>------------------------------------------------------------
        !!  dz used for calculating advection
        !!------------------------------------------------------------
        else if (var_idx==kVARS%advection_dz) then
            var_meta%name        = "advection_dz"
            var_meta%dimensions  = three_d_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "dz used for calculating advection"),                 &
                                attribute_t("units",         "m"),                                   &
                                attribute_t("coordinates",   "lat lon")]                    
        !>------------------------------------------------------------
        !!  Global Terrain for the domain
        !!------------------------------------------------------------
        else if (var_idx==kVARS%global_terrain) then
            var_meta%name        = "global_terrain"
            var_meta%dimensions  = two_d_global_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "Domain terrain"),                 &
                                attribute_t("units",         "m"),                                   &
                                attribute_t("coordinates",   "lat lon")]
        !>------------------------------------------------------------
        !!  Neighbor Terrain for the domain
        !!------------------------------------------------------------
        else if (var_idx==kVARS%neighbor_terrain) then
            var_meta%name        = "neighbor_terrain"
            var_meta%dimensions  = two_d_neighbor_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "Domain terrain"),                 &
                                attribute_t("units",         "m"),                                   &
                                attribute_t("coordinates",   "lat lon")]
                    
        !>------------------------------------------------------------
        !!  Terrain slope
        !!------------------------------------------------------------
        else if (var_idx==kVARS%slope) then
            var_meta%name        = "slope"
            var_meta%dimensions  = two_d_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "slope of the terrain"),                 &
                                attribute_t("units",         "rad"),                                   &
                                attribute_t("coordinates",   "lat lon")]
        !>------------------------------------------------------------
        !!  Terrain slope in degrees
        !!------------------------------------------------------------
        else if (var_idx==kVARS%slope_angle) then
            var_meta%name        = "slope_angle"
            var_meta%dimensions  = two_d_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "slope of the terrain"),                 &
                                attribute_t("units",         "degrees"),                                   &
                                attribute_t("coordinates",   "lat lon")]
        !>------------------------------------------------------------
        !!  Terrain aspect
        !!------------------------------------------------------------
        else if (var_idx==kVARS%aspect_angle) then
            var_meta%name        = "aspect_angle"
            var_meta%dimensions  = two_d_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "Aspect of the terrain"),                 &
                                attribute_t("units",         "degrees"),                                   &
                                attribute_t("coordinates",   "lat lon")]
        !>------------------------------------------------------------
        !!  Sky view fraction
        !!------------------------------------------------------------
        else if (var_idx==kVARS%svf) then
            var_meta%name        = "svf"
            var_meta%dimensions  = two_d_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "sky view fraction"),                 &
                                attribute_t("units",         "-"),                                   &
                                attribute_t("coordinates",   "lat lon")]
        !>------------------------------------------------------------
        !!  Horizon line matrix
        !!------------------------------------------------------------
        else if (var_idx==kVARS%hlm) then
            var_meta%name        = "hlm"
            var_meta%dimensions  = three_d_hlm_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "Horizon line matrix"),                 &
                                attribute_t("units",         "degrees"),                                   &
                                attribute_t("coordinates",   "lat lon")]
        !>------------------------------------------------------------
        !!  Snow holding depth of terrain
        !!------------------------------------------------------------
        else if (var_idx==kVARS%shd) then
            var_meta%name        = "shd"
            var_meta%dimensions  = two_d_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "Snow holding depth of terrain"),                 &
                                attribute_t("units",         "m"),                                   &
                                attribute_t("coordinates",   "lat lon")]
                                       
        !>------------------------------------------------------------
        !!  Cloud cover fraction
        !!------------------------------------------------------------
        else if (var_idx==kVARS%cloud_fraction) then
            var_meta%name        = "clt"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "cloud_area_fraction"),                 &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        

        !>------------------------------------------------------------
        !!  Effective cloud droplet radius
        !!------------------------------------------------------------
        else if (var_idx==kVARS%re_cloud) then
            var_meta%name        = "re_cloud"
            var_meta%dimensions  = three_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "effective_radius_of_cloud_liquid_water_particles"), &
                               attribute_t("units",         "m"),                                                &
                               attribute_t("coordinates",   "lat lon")]
        

        !>------------------------------------------------------------
        !!  Effective cloud ice radius
        !!------------------------------------------------------------
        else if (var_idx==kVARS%re_ice) then
            var_meta%name        = "re_ice"
            var_meta%dimensions  = three_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "effective_radius_of_stratiform_cloud_ice_particles"), &
                               attribute_t("units",         "m"),                                                  &
                               attribute_t("coordinates",   "lat lon")]
        

        !>------------------------------------------------------------
        !!  Effective snow radius
        !!------------------------------------------------------------
        else if (var_idx==kVARS%re_snow) then
            var_meta%name        = "re_snow"
            var_meta%dimensions  = three_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "effective_radius_of_stratiform_snow_particles"), &
                               attribute_t("units",         "m"),                                                 &
                               attribute_t("coordinates",   "lat lon")]
        

        !>------------------------------------------------------------
        !!  Mass-weighted effective density of ice1 category
        !!------------------------------------------------------------
        else if (var_idx==kVARS%ice1_rho) then
            var_meta%name        = "ice1_rho"
            var_meta%dimensions  = three_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "mass-weighted_effective_density_of_ice1"), &
                               attribute_t("units",         "(kg m^-3)"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        

        !>------------------------------------------------------------
        !!  Number-weighted aspect ratio of ice1 category
        !!------------------------------------------------------------
        else if (var_idx==kVARS%ice1_phi) then
            var_meta%name        = "ice1_AR"
            var_meta%dimensions  = three_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "number-weighted_aspect_ratio_of_ice1_category"), &
                               attribute_t("units",         "-"),                                                 &
                               attribute_t("coordinates",   "lat lon")]
        

        !>------------------------------------------------------------
        !!  Mass-weighted fall speed of ice1 category
        !!------------------------------------------------------------
        else if (var_idx==kVARS%ice1_vmi) then
            var_meta%name        = "ice1_FS"
            var_meta%dimensions  = three_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "planar-nucleated mass-weighted fall speeds (m s^-1)"), &
                               attribute_t("units",         "-"),                                                 &
                               attribute_t("coordinates",   "lat lon")]
        


        !>------------------------------------------------------------
        !!  Mass-weighted effective density of ice2 category
        !!------------------------------------------------------------
        else if (var_idx==kVARS%ice2_rho) then
            var_meta%name        = "ice2_rho"
            var_meta%dimensions  = three_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "mass-weighted_effective_density_of_ice2"), &
                               attribute_t("units",         "(kg m^-3)"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        

        !>------------------------------------------------------------
        !!  Number-weighted aspect ratio of ice2 category
        !!------------------------------------------------------------
        else if (var_idx==kVARS%ice2_phi) then
            var_meta%name        = "ice2_AR"
            var_meta%dimensions  = three_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "number-weighted_aspect_ratio_of_ice2_category"), &
                               attribute_t("units",         "-"),                                                 &
                               attribute_t("coordinates",   "lat lon")]
        

        !>------------------------------------------------------------
        !!  Mass-weighted fall speed of ice2 category
        !!------------------------------------------------------------
        else if (var_idx==kVARS%ice2_vmi) then
            var_meta%name        = "ice2_FS"
            var_meta%dimensions  = three_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "columnar-nucleated mass-weighted fall speeds (m s^-1)"), &
                               attribute_t("units",         "-"),                                                 &
                               attribute_t("coordinates",   "lat lon")]
        

        !>------------------------------------------------------------
        !!  Mass-weighted effective density of ice3 category
        !!------------------------------------------------------------
        else if (var_idx==kVARS%ice3_rho) then
            var_meta%name        = "ice3_rho"
            var_meta%dimensions  = three_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "mass-weighted_effective_density_of_ice3"), &
                               attribute_t("units",         "(kg m^-3)"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        

        !>------------------------------------------------------------
        !!  Number-weighted aspect ratio of ice3 category
        !!------------------------------------------------------------
        else if (var_idx==kVARS%ice3_phi) then
            var_meta%name        = "ice3_AR"
            var_meta%dimensions  = three_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "number-weighted_aspect_ratio_of_ice3_category"), &
                               attribute_t("units",         "-"),                                                 &
                               attribute_t("coordinates",   "lat lon")]
        

        !>------------------------------------------------------------
        !!  Mass-weighted fall speed of ice3 category
        !!------------------------------------------------------------
        else if (var_idx==kVARS%ice3_vmi) then
            var_meta%name        = "ice3_FS"
            var_meta%dimensions  = three_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "aggregates mass-weighted fall speeds (m s^-1)"), &
                               attribute_t("units",         "-"),                                                 &
                               attribute_t("coordinates",   "lat lon")]
        


        !>------------------------------------------------------------
        !!  Alpha weighting factor in variational wind solver
        !!------------------------------------------------------------
        else if (var_idx==kVARS%wind_alpha) then
            var_meta%name        = "wind_alpha"
            var_meta%dimensions  = three_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "Alpha weighting factor in variational wind solver"), &
                               attribute_t("units",         "-"),                                                 &
                               attribute_t("coordinates",   "lat lon")]
        

        !>------------------------------------------------------------
        !!  Froude number
        !!------------------------------------------------------------
        else if (var_idx==kVARS%froude) then
            var_meta%name        = "froude"
            var_meta%dimensions  = three_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "Froude number used to calculate wind_alpha"), &
                               attribute_t("units",         "-"),                                                 &
                               attribute_t("coordinates",   "lat lon")]
        

        !>------------------------------------------------------------
        !!  Bulk richardson number
        !!------------------------------------------------------------
        else if (var_idx==kVARS%blk_ri) then
            var_meta%name        = "blk_ri"
            var_meta%dimensions  = three_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "Bulk richardson number used to calculate Sx sheltering"), &
                               attribute_t("units",         "-"),                                                 &
                               attribute_t("coordinates",   "lat lon")]
        !>------------------------------------------------------------
        !!  blocking terrain from froude number calculation
        !!------------------------------------------------------------
        else if (var_idx==kVARS%froude_terrain) then
            var_meta%name        = "froude_terrain"
            var_meta%dimensions  = four_d_azim_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "Blocking terrain height used in calculation of froude number"), &
                                attribute_t("units",         "m"),                                                 &
                                attribute_t("coordinates",   "lat lon")]
        !>------------------------------------------------------------
        !!  Sx
        !!------------------------------------------------------------
        else if (var_idx==kVARS%Sx) then
            var_meta%name        = "Sx"
            var_meta%dimensions  = four_d_azim_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "Maximum upwind slope"), &
                                attribute_t("units",         "-"),                                                 &
                                attribute_t("coordinates",   "lat lon")]
        !>------------------------------------------------------------
        !!  TPI
        !!------------------------------------------------------------
        else if (var_idx==kVARS%TPI) then
            var_meta%name        = "TPI"
            var_meta%dimensions  = two_d_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "Topographic position index"), &
                                attribute_t("units",         "-"),                                                 &
                                attribute_t("coordinates",   "lat lon")]
                                                                                                                                                                                                        

        !>------------------------------------------------------------
        !!  Outgoing longwave radiation
        !!------------------------------------------------------------
        else if (var_idx==kVARS%out_longwave_rad) then
            var_meta%name        = "rlut"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "toa_outgoing_longwave_flux"), &
                               attribute_t("units",         "W m-2"),                      &
                               attribute_t("coordinates",   "lat lon")]
        

        !>------------------------------------------------------------
        !!  Longwave cloud forcing
        !!------------------------------------------------------------
        else if (var_idx==kVARS%longwave_cloud_forcing) then
            var_meta%name        = "lwcf"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "longwave_cloud_forcing"), &
                               attribute_t("units",         "W m-2"),                      &
                               attribute_t("coordinates",   "lat lon")]
        

        !>------------------------------------------------------------
        !!  Shortwave cloud forcing
        !!------------------------------------------------------------
        else if (var_idx==kVARS%shortwave_cloud_forcing) then
            var_meta%name        = "swcf"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "shortwave_cloud_forcing"), &
                               attribute_t("units",         "W m-2"),                       &
                               attribute_t("coordinates",   "lat lon")]
        

        !>------------------------------------------------------------
        !!  Cosine solar zenith angle
        !!------------------------------------------------------------
        else if (var_idx==kVARS%cosine_zenith_angle) then
            var_meta%name        = "cosz"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "cosine_zenith_angle"), &
                               attribute_t("units",         " "),                       &
                               attribute_t("coordinates",   "lat lon")]
        

        !>------------------------------------------------------------
        !!  Tendency from short wave radiation
        !!------------------------------------------------------------
        else if (var_idx==kVARS%tend_swrad) then
            var_meta%name        = "tend_swrad"
            var_meta%dimensions  = three_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "sw_rad_tend"), &
                               attribute_t("units",         " "),               &
                               attribute_t("coordinates",   "lat lon")]
        


        !>------------------------------------------------------------
        !!  Surface emissivity
        !!------------------------------------------------------------
        else if (var_idx==kVARS%land_emissivity) then
            var_meta%name        = "emiss"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "surface_longwave_emissivity"), &
                               attribute_t("units",         " "),                           &
                               attribute_t("coordinates",   "lat lon")]
        

        !>------------------------------------------------------------
        !!  Temperature on interface
        !!------------------------------------------------------------
        else if (var_idx==kVARS%temperature_interface) then
            var_meta%name        = "temperature_i"
            var_meta%dimensions  = three_d_t_interface_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "air_temperature"),                    &
                               attribute_t("long_name",     "Temperature"),                        &
                               attribute_t("units",         "K"),                                  &
                               attribute_t("coordinates",   "lat lon")]
        



        !>------------------------------------------------------------
        !!  Downward Shortwave Radiation at the Surface (positive down)
        !!------------------------------------------------------------
        else if (var_idx==kVARS%shortwave) then
            var_meta%name        = "rsds"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "surface_downwelling_shortwave_flux_in_air"), &
                               attribute_t("units",         "W m-2"),                                     &
                               attribute_t("coordinates",   "lat lon")]
            if (present(force_boundaries)) force_boundaries = .False.
            if (present(opt) .and. present(forcing_var)) forcing_var = (opt%forcing%swdown_var /= "")

        !>------------------------------------------------------------
        !!  MJ: in OSHD as 'sdri' referring to 'direct shortwave radiation, per inclined surface area' accounted for shading and slope effects. Tobias Jonas (TJ) scheme based on swr function in metDataWizard/PROCESS_COSMO_DATA_1E2E.m and also https://github.com/Tobias-Jonas-SLF/HPEval
        !!------------------------------------------------------------
        else if (var_idx==kVARS%shortwave_direct) then
            var_meta%name        = "swtb"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "surface_direct_downwelling_shortwave_flux_in_air"), &
                               attribute_t("long_name",     "direct shortwave radiation"), &
                               attribute_t("units",         "W m-2"),                                            &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  MJ: in OSHD as 'sdfd' referring to 'diffuse shortwave radiation, per horizontal surface area' accounted for partioning (based on transmisivity) and sky view fraction. Tobias Jonas (TJ) scheme based on swr function in metDataWizard/PROCESS_COSMO_DATA_1E2E.m and also https://github.com/Tobias-Jonas-SLF/HPEval
        !!------------------------------------------------------------
        else if (var_idx==kVARS%shortwave_diffuse) then
            var_meta%name        = "swtd"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "surface_diffuse_downwelling_shortwave_flux_in_air"), &
                               attribute_t("long_name",     "diffuse shortwave radiation"), &
                               attribute_t("units",         "W m-2"),                                             &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  MJ: in OSHD as 'sdrd' referring to 'direct shortwave radiation, per horizontal surface area' only accounted for shading but not the slope effects. Tobias Jonas (TJ) scheme based on swr function in metDataWizard/PROCESS_COSMO_DATA_1E2E.m and also https://github.com/Tobias-Jonas-SLF/HPEval
        !!------------------------------------------------------------
        else if (var_idx==kVARS%shortwave_direct_above) then
            var_meta%name        = "sdrd"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "direct shortwave radiation, per horizontal surface area"), &
                               attribute_t("units",         "W m-2"),                                             &
                               attribute_t("coordinates",   "lat lon")]
                                
        !>------------------------------------------------------------
        !!  Downward Longwave Radiation at the Surface (positive down)
        !!------------------------------------------------------------
        else if (var_idx==kVARS%longwave) then
            var_meta%name        = "lwtr"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "surface_net_downward_longwave_flux"), &
                               attribute_t("units",         "W m-2"),                                    &
                               attribute_t("coordinates",   "lat lon")]
            if (present(force_boundaries)) force_boundaries = .False.
            if (present(opt) .and. present(forcing_var)) forcing_var = (opt%forcing%lwdown_var /= "")

        !>------------------------------------------------------------
        !!  Total Absorbed Solar Radiation
        !!------------------------------------------------------------
        else if (var_idx==kVARS%rad_absorbed_total) then
            var_meta%name        = "rad_absorbed_total"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "total_absorbed_radiation"),             &
                               attribute_t("units",         "W m-2"),                                    &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Solar Radiation Absorbed by Vegetation
        !!------------------------------------------------------------
        else if (var_idx==kVARS%rad_absorbed_veg) then
            var_meta%name        = "rad_absorbed_veg"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "radiation_absorbed_by_vegetation"),     &
                               attribute_t("units",         "W m-2"),                                    &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Solar Radiation Absorbed by Bare Ground
        !!------------------------------------------------------------
        else if (var_idx==kVARS%rad_absorbed_bare) then
            var_meta%name        = "rad_absorbed_bare"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "radiation_absorbed_by_bare_ground"),    &
                               attribute_t("units",         "W m-2"),                                    &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Net Longwave Radiation (positive up)
        !!------------------------------------------------------------
        else if (var_idx==kVARS%rad_net_longwave) then
            var_meta%name        = "rad_net_longwave"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "net_upward_longwave_flux_in_air"),          &
                               attribute_t("units",         "W m-2"),                                    &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Upward Longwave Radiation at the Surface (positive up)
        !!------------------------------------------------------------
        else if (var_idx==kVARS%longwave_up) then
            var_meta%name        = "rlus"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "surface_upwelling_longwave_flux_in_air"),   &
                               attribute_t("units",         "W m-2"),                                    &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Ground Heat Flux (positive down)
        !!------------------------------------------------------------
        else if (var_idx==kVARS%ground_heat_flux) then
            var_meta%name        = "hfgs"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "upward_heat_flux_at_ground_level_in_soil"), &
                               attribute_t("units",         "W m-2"),                                    &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Vegetation fraction
        !!------------------------------------------------------------
        else if (var_idx==kVARS%vegetation_fraction) then
            var_meta%name        = "vegetation_fraction"
            var_meta%dimensions  = three_d_t_month_dimensions
            var_meta%dim_len(2)  = kMONTH_GRID_Z
            var_meta%attributes  = [attribute_t("standard_name", "vegetation_area_fraction"),            &
                               attribute_t("units",         "m2 m-2"),                              &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Annual Maximum Vegetation Fraction
        !!------------------------------------------------------------
        else if (var_idx==kVARS%vegetation_fraction_max) then
            var_meta%name        = "vegetation_fraction_max"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "max_vegetation_area_fraction"),    &
                               attribute_t("units",         "m2 m-2"),                              &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Noah-MP Output Vegetation Fraction
        !!------------------------------------------------------------
        else if (var_idx==kVARS%vegetation_fraction_out) then
            var_meta%name        = "vegetation_fraction_out"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "vegetation_fraction_out"),         &
                               attribute_t("units",         "m2 m-2"),                              &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Land cover type
        !!------------------------------------------------------------
        else if (var_idx==kVARS%veg_type) then
            var_meta%name        = "veg_type"
            var_meta%dimensions  = two_d_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "vegetation_type"),                 &
                               attribute_t("units",      "1"),                                      &
                               attribute_t("coordinates",   "lat lon")]
            var_meta%dtype = kINTEGER
        !>------------------------------------------------------------
        !!  Land cover type
        !!------------------------------------------------------------
        else if (var_idx==kVARS%soil_type) then
            var_meta%name        = "soil_type"
            var_meta%dimensions  = two_d_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "soil_type"),                 &
                               attribute_t("units",      "1"),                                      &
                               attribute_t("coordinates",   "lat lon")]
            var_meta%dtype = kINTEGER

        !>------------------------------------------------------------
        !!  Leaf Mass
        !!------------------------------------------------------------
        else if (var_idx==kVARS%mass_leaf) then
            var_meta%name        = "mass_leaf"
            var_meta%dimensions  = two_d_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "leaf_mass"),                      &
                               attribute_t("units",      "g m-2"),                                 &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Root Mass
        !!------------------------------------------------------------
        else if (var_idx==kVARS%mass_root) then
            var_meta%name        = "mass_root"
            var_meta%dimensions  = two_d_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "root_mass"),                      &
                               attribute_t("units",      "g m-2"),                                 &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Stem Mass
        !!------------------------------------------------------------
        else if (var_idx==kVARS%mass_stem) then
            var_meta%name        = "mass_stem"
            var_meta%dimensions  = two_d_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "stem_mass"),                      &
                               attribute_t("units",      "g m-2"),                                 &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Wood Mass
        !!------------------------------------------------------------
        else if (var_idx==kVARS%mass_wood) then
            var_meta%name        = "mass_wood"
            var_meta%dimensions  = two_d_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "wood_mass"),                      &
                               attribute_t("units",      "g m-2"),                                 &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Leaf Area Index
        !!------------------------------------------------------------
        else if (var_idx==kVARS%lai) then
            var_meta%name        = "lai"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "leaf_area_index"),                     &
                               attribute_t("units",         "m2 m-2"),                              &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Stem Area Index
        !!------------------------------------------------------------
        else if (var_idx==kVARS%sai) then
            var_meta%name        = "sai"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "stem_area_index"),                 &
                               attribute_t("units",         "m2 m-2"),                              &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Planting Date
        !!------------------------------------------------------------
        else if (var_idx==kVARS%date_planting) then
            var_meta%name        = "date_planting"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "planting_date"),                   &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Harvesting Date
        !!------------------------------------------------------------
        else if (var_idx==kVARS%date_harvest) then
            var_meta%name        = "date_harvest"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "harvest_date"),                    &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Crop Category
        !!------------------------------------------------------------
        else if (var_idx==kVARS%crop_category) then
            var_meta%name        = "crop_category"
            var_meta%dimensions  = two_d_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "crop_category"),                   &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
            var_meta%dtype = kINTEGER

        !>------------------------------------------------------------
        !!  Crop Type (Noah-MP initialization variable)
        !!------------------------------------------------------------
        else if (var_idx==kVARS%crop_type) then
            var_meta%name        = "crop_type"
            var_meta%dimensions  = three_d_crop_dimensions
            var_meta%dim_len(2)  = kCROP_GRID_Z
            var_meta%attributes  = [attribute_t("non_standard_name", "crop_type"),                       &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
            var_meta%dtype = kINTEGER

        !>------------------------------------------------------------
        !!  Growing Season Growing Degree Days
        !!------------------------------------------------------------
        else if (var_idx==kVARS%growing_season_gdd) then
            var_meta%name        = "growing_season_gdd"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "growing_season_gdd"),              &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Irrigation Event Number, Sprinkler
        !!------------------------------------------------------------
        else if (var_idx==kVARS%irr_eventno_sprinkler) then
            var_meta%name        = "irr_eventno_sprinkler"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "irr_eventno_sprinkler"),           &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
            var_meta%dtype = kINTEGER

        !>------------------------------------------------------------
        !!  Irrigation Fraction
        !!------------------------------------------------------------
        else if (var_idx==kVARS%irr_frac_total) then
            var_meta%name        = "irr_frac_total"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "irr_frac_total"),                  &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Irrigation Fraction, Sprinkler
        !!------------------------------------------------------------
        else if (var_idx==kVARS%irr_frac_sprinkler) then
            var_meta%name        = "irr_frac_sprinkler"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "irr_frac_sprinkler"),              &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Irrigation Fraction, Micro
        !!------------------------------------------------------------
        else if (var_idx==kVARS%irr_frac_micro) then
            var_meta%name        = "irr_frac_micro"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "irr_frac_micro"),                  &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Irrigation Fraction, Flood
        !!------------------------------------------------------------
        else if (var_idx==kVARS%irr_frac_flood) then
            var_meta%name        = "irr_frac_flood"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "irr_frac_flood"),                  &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Irrigation Event Number, Micro
        !!------------------------------------------------------------
        else if (var_idx==kVARS%irr_eventno_micro) then
            var_meta%name        = "irr_eventno_micro"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "irr_eventno_micro"),               &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
            var_meta%dtype = kINTEGER

        !>------------------------------------------------------------
        !!  Irrigation Event Number, Flood
        !!------------------------------------------------------------
        else if (var_idx==kVARS%irr_eventno_flood) then
            var_meta%name        = "irr_eventno_flood"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "irr_eventno_flood"),               &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
            var_meta%dtype = kINTEGER

        !>------------------------------------------------------------
        !!  Irrigation Water Amount to be Applied, Sprinkler
        !!------------------------------------------------------------
        else if (var_idx==kVARS%irr_alloc_sprinkler) then
            var_meta%name        = "irr_alloc_sprinkler"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "irr_alloc_sprinkler"),             &
                               attribute_t("units",         "m"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Irrigation Water Amount to be Applied, Micro
        !!------------------------------------------------------------
        else if (var_idx==kVARS%irr_alloc_micro) then
            var_meta%name        = "irr_alloc_micro"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "irr_alloc_micro"),                 &
                               attribute_t("units",         "m"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Irrigation Water Amount to be Applied, Flood
        !!------------------------------------------------------------
        else if (var_idx==kVARS%irr_alloc_flood) then
            var_meta%name        = "irr_alloc_flood"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "irr_alloc_flood"),                 &
                               attribute_t("units",         "m"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Loss of Sprinkler Irrigation Water to Evaporation
        !!------------------------------------------------------------
        else if (var_idx==kVARS%irr_evap_loss_sprinkler) then
            var_meta%name        = "irr_evap_loss_sprinkler"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "irr_evap_loss_sprinkler"),         &
                               attribute_t("units",         "m timestep-1"),                        &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Amount of Irrigation, Sprinkler
        !!------------------------------------------------------------
        else if (var_idx==kVARS%irr_amt_sprinkler) then
            var_meta%name        = "irr_amt_sprinkler"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "irr_amt_sprinkler"),               &
                               attribute_t("units",         "mm"),                                  &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Amount of Irrigation, Micro
        !!------------------------------------------------------------
        else if (var_idx==kVARS%irr_amt_micro) then
            var_meta%name        = "irr_amt_micro"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "irr_amt_micro"),                   &
                               attribute_t("units",         "mm"),                                  &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Amount of Irrigation, Flood
        !!------------------------------------------------------------
        else if (var_idx==kVARS%irr_amt_flood) then
            var_meta%name        = "irr_amt_flood"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "irr_amt_flood"),                   &
                               attribute_t("units",         "mm"),                                  &
                               attribute_t("coordinates",   "lat lon")]
        !>------------------------------------------------------------
        !!  latent heating from sprinkler evaporation
        !!------------------------------------------------------------
        else if (var_idx==kVARS%evap_heat_sprinkler) then
            var_meta%name        = "evap_heat_sprinkler"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "latent heating from sprinkler evaporation"),                   &
                                attribute_t("units",         "W/m^2"),                                  &
                                attribute_t("coordinates",   "lat lon")]
                    
        
        !>------------------------------------------------------------
        !!  Mass of Agricultural Grain Produced
        !!------------------------------------------------------------
        else if (var_idx==kVARS%mass_ag_grain) then
            var_meta%name        = "mass_ag_grain"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "mass_agricultural_grain"),         &
                               attribute_t("units",         "g m-2"),                               &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Growing Degree Days (based on 10C)
        !!------------------------------------------------------------
        else if (var_idx==kVARS%growing_degree_days) then
            var_meta%name        = "growing_degree_days"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "growing_degree_days"),             &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Plant Growth Stage
        !!------------------------------------------------------------
        else if (var_idx==kVARS%plant_growth_stage) then
            var_meta%name        = "plant_growth_stage"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "plant_growth_stage"),              &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
            var_meta%dtype = kINTEGER

        !>------------------------------------------------------------
        !!  Net Ecosystem Exchange
        !!------------------------------------------------------------
        else if (var_idx==kVARS%net_ecosystem_exchange) then
            var_meta%name        = "net_ecosystem_exchange"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "net_ecosystem_exchange_expressed_as_carbon_dioxide"), &
                               attribute_t("units",         "g m-2 s-1"),                                              &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Gross Primary Productivity
        !!------------------------------------------------------------
        else if (var_idx==kVARS%gross_primary_prod) then
            var_meta%name        = "gross_primary_prod"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "gross_primary_productivity_of_biomass_expressed_as_carbon"), &
                               attribute_t("units",         "g m-2 s-1"),                                                 &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Net Primary Productivity
        !!------------------------------------------------------------
        else if (var_idx==kVARS%net_primary_prod) then
            var_meta%name        = "net_primary_prod"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "net_primary_productivity_of_biomass_expressed_as_carbon"), &
                               attribute_t("units",         "g m-2 s-1"),                                               &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Absorbed Photosynthetically Active Radiation
        !!------------------------------------------------------------
        else if (var_idx==kVARS%apar) then
            var_meta%name        = "apar"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "absorbed_photosynthetically_active_radiation"), &
                               attribute_t("units",         "W m-2"),                                            &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Total Photosynthesis
        !!------------------------------------------------------------
        else if (var_idx==kVARS%photosynthesis_total) then
            var_meta%name        = "photosynthesis_total"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "total_photosynthesis_expressed_as_carbon_dioxide"), &
                               attribute_t("units",         "mmol m-2 s-1"),                                         &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Total Leaf Stomatal Resistance
        !!------------------------------------------------------------
        else if (var_idx==kVARS%stomatal_resist_total) then
            var_meta%name        = "stomatal_resist_total"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "total_leaf_stomatal_resistance"),                   &
                               attribute_t("units",         "s m-1"),                                                &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Sunlit Leaf Stomatal Resistance
        !!------------------------------------------------------------
        else if (var_idx==kVARS%stomatal_resist_sun) then
            var_meta%name        = "stomatal_resist_sun"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "sunlif_leaf_stomatal_resistance"),                  &
                               attribute_t("units",         "s m-1"),                                                &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Shaded Leaf Stomatal Resistance
        !!------------------------------------------------------------
        else if (var_idx==kVARS%stomatal_resist_shade) then
            var_meta%name        = "stomatal_resist_shade"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "shaded_leaf_stomatal_resistance"),                 &
                               attribute_t("units",         "s m-1"),                                               &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  GECROS (Genotype-by-Envrionment interaction on CROp growth Simulator [Yin and van Laar, 2005]) crop model state
        !!------------------------------------------------------------
        else if (var_idx==kVARS%gecros_state) then
            var_meta%name        = "gecros_state"
            var_meta%dimensions  = three_d_t_gecros_dimensions
            var_meta%dim_len(2)  = kGECROS_GRID_Z
            var_meta%attributes  = [attribute_t("non_standard_name", "gecros_state"),                                    &
                               attribute_t("units",         "N/A"),                                               &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Canopy Total Water Content
        !!------------------------------------------------------------
        else if (var_idx==kVARS%canopy_water) then
            var_meta%name        = "canopy_water"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "canopy_water_amount"),                 &
                               attribute_t("units",         "kg m-2"),                              &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Canopy Frozen Water Content
        !!------------------------------------------------------------
        else if (var_idx==kVARS%canopy_water_ice) then
            var_meta%name        = "canopy_ice"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "canopy_snow_amount"),                  &
                               attribute_t("units",         "kg m-2"),                              &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Canopy Liquid Water Content (in snow)
        !!------------------------------------------------------------
        else if (var_idx==kVARS%canopy_water_liquid) then
            var_meta%name        = "canopy_liquid"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "canopy_liquid_water_amount"),      &
                               attribute_t("units",         "kg m-2"),                              &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Canopy Air Vapor Pressure
        !!------------------------------------------------------------
        else if (var_idx==kVARS%canopy_vapor_pressure) then
            var_meta%name        = "canopy_vapor_pressure"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "canopy_air_vapor_pressure"),       &
                               attribute_t("units",         "Pa"),                                  &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Canopy Air Temperature
        !!------------------------------------------------------------
        else if (var_idx==kVARS%canopy_temperature) then
            var_meta%name        = "canopy_temperature"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "canopy_air_temperature"),          &
                               attribute_t("units",         "K"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Canopy Wetted/Snowed Fraction
        !!------------------------------------------------------------
        else if (var_idx==kVARS%canopy_fwet) then
            var_meta%name        = "canopy_fwet"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "canopy_wetted_fraction"),          &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Vegetation Leaf Temperature
        !!------------------------------------------------------------
        else if (var_idx==kVARS%veg_leaf_temperature) then
            var_meta%name        = "veg_leaf_temperature"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "veg_leaf_temperature"),            &
                               attribute_t("units",         "K"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Ground Surface Temperature
        !!------------------------------------------------------------
        else if (var_idx==kVARS%ground_surf_temperature) then
            var_meta%name        = "ground_surf_temperature"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "ground_surface_temperature"),      &
                               attribute_t("units",         "K"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Between Gap Fraction
        !!------------------------------------------------------------
        else if (var_idx==kVARS%frac_between_gap) then
            var_meta%name        = "frac_between_gap"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "between_gap_fraction"),            &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Within Gap Fraction
        !!------------------------------------------------------------
        else if (var_idx==kVARS%frac_within_gap) then
            var_meta%name        = "frac_within_gap"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "within_gap_fraction"),             &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Under-Canopy Ground Temperature
        !!------------------------------------------------------------
        else if (var_idx==kVARS%ground_temperature_canopy) then
            var_meta%name        = "ground_temperature_canopy"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "under_canopy_ground_temperature"), &
                               attribute_t("units",         "K"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Bare Ground Temperature
        !!------------------------------------------------------------
        else if (var_idx==kVARS%ground_temperature_bare) then
            var_meta%name        = "ground_temperature_bare"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "bare_ground_temperature"),         &
                               attribute_t("units",         "K"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Snowfall on the Ground
        !!------------------------------------------------------------
        else if (var_idx==kVARS%snowfall_ground) then
            var_meta%name        = "snowfall_ground"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "ground_snow_rate"),                &
                               attribute_t("units",         "mm s-1"),                              &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Rainfall on the Ground
        !!------------------------------------------------------------
        else if (var_idx==kVARS%rainfall_ground) then
            var_meta%name        = "rainfall_ground"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "ground_rain_rate"),                &
                               attribute_t("units",         "mm s-1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        !>------------------------------------------------------------
        !!  Last snowfall saved for LSM
        !!------------------------------------------------------------
        else if (var_idx==kVARS%lsm_last_snow) then
            var_meta%name        = "lsm_last_snow"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "snowfall at last lsm call"),                &
                                attribute_t("units",         "mm"),                              &
                                attribute_t("coordinates",   "lat lon")]
                            
        !>------------------------------------------------------------
        !!  Last precipitation saved for LSM
        !!------------------------------------------------------------
        else if (var_idx==kVARS%lsm_last_precip) then
            var_meta%name        = "lsm_last_precip"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "precipitation at last lsm call"),                &
                                attribute_t("units",         "mm"),                                   &
                                attribute_t("coordinates",   "lat lon")]
                    
        !>------------------------------------------------------------
        !!  Snow water equivalent on the surface
        !!------------------------------------------------------------
        else if (var_idx==kVARS%snow_water_equivalent) then
            var_meta%name        = "swet"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "snow water equivalent"),                 &
                               attribute_t("units",         "kg m-2"),                              &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Snow water equivalent from previous timestep
        !!------------------------------------------------------------
        else if (var_idx==kVARS%snow_water_eq_prev) then
            var_meta%name        = "swe_0"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "surface_snow_amount_prev"),        &
                               attribute_t("units",         "kg m-2"),                              &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Grid cell albedo
        !!------------------------------------------------------------
        else if (var_idx==kVARS%albedo) then
            var_meta%name        = "albedo"
            var_meta%dimensions  = three_d_t_month_dimensions
            var_meta%dim_len(2)  = kMONTH_GRID_Z
            var_meta%attributes  = [attribute_t("non_standard_name", "albedo"),            &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        

        !>------------------------------------------------------------
        !!  Snow albedo from previous timestep
        !!------------------------------------------------------------
        else if (var_idx==kVARS%snow_albedo_prev) then
            var_meta%name        = "snow_albedo_0"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "snowpack_albedo_prev"),            &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Snow Temperature
        !!------------------------------------------------------------
        else if (var_idx==kVARS%snow_temperature) then
            var_meta%name        = "snow_temperature"
            var_meta%dimensions  = three_d_t_snow_dimensions
            var_meta%dim_len(2)  = kSNOW_GRID_Z
            var_meta%attributes  = [attribute_t("standard_name", "temperature_in_surface_snow"),         &
                               attribute_t("units",         "K"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Snow Layer Depth
        !!------------------------------------------------------------
        else if (var_idx==kVARS%snow_layer_depth) then
            var_meta%name        = "snow_layer_depth"
            var_meta%dimensions  = three_d_t_snowsoil_dimensions
            var_meta%dim_len(2)  = kSNOWSOIL_GRID_Z
            var_meta%attributes  = [attribute_t("non_standard_name", "snow_layer_depth"),                &
                               attribute_t("units",         "m"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Snow Layer Ice
        !!------------------------------------------------------------
        else if (var_idx==kVARS%snow_layer_ice) then
            var_meta%name        = "snow_layer_ice"
            var_meta%dimensions  = three_d_t_snow_dimensions
            var_meta%dim_len(2)  = kSNOW_GRID_Z
            var_meta%attributes  = [attribute_t("non_standard_name", "snow_layer_ice_content"),          &
                               attribute_t("units",         "mm"),                                  &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Snow Layer Liquid Water
        !!------------------------------------------------------------
        else if (var_idx==kVARS%snow_layer_liquid_water) then
            var_meta%name        = "snow_layer_liquid_water"
            var_meta%dimensions  = three_d_t_snow_dimensions
            var_meta%dim_len(2)  = kSNOW_GRID_Z
            var_meta%attributes  = [attribute_t("non_standard_name", "snow_layer_liquid_water_content"), &
                               attribute_t("units",         "mm"),                                  &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Snow Age Factor
        !!------------------------------------------------------------
        else if (var_idx==kVARS%snow_age_factor) then
            var_meta%name        = "tau_ss"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "snow_age_factor"),                 &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Snow height
        !!------------------------------------------------------------
        else if (var_idx==kVARS%snow_height) then
            var_meta%name        = "snow_height"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "surface_snow_height"),                 &
                               attribute_t("units",         "m"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Number of Snowpack Layers
        !!------------------------------------------------------------
        else if (var_idx==kVARS%snow_nlayers) then
            var_meta%name        = "snow_nlayers"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%dtype      = kINTEGER
            var_meta%attributes  = [attribute_t("non_standard_name", "snow_nlayers"),                    &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Soil water content
        !!------------------------------------------------------------
        else if (var_idx==kVARS%soil_water_content) then
            var_meta%name        = "soil_water_content"
            var_meta%dimensions  = three_d_t_soil_dimensions
            var_meta%dim_len(2)  = kSOIL_GRID_Z
            var_meta%attributes  = [attribute_t("standard_name", "moisture_content_of_soil_layer"),      &
                               attribute_t("units",         "m3 m-3"),                              &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Soil water content, of liquid water
        !!------------------------------------------------------------
        else if (var_idx==kVARS%soil_water_content_liq) then
            var_meta%name        = "soil_water_content_liq"
            var_meta%dimensions  = three_d_t_soil_dimensions
            var_meta%dim_len(2)  = kSOIL_GRID_Z
            var_meta%attributes  = [attribute_t("standard_name", "liquid_moisture_content_of_soil_layer"),      &
                               attribute_t("units",         "m3 m-3"),                              &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Equilibrium Volumetric Soil Moisture
        !!------------------------------------------------------------
        else if (var_idx==kVARS%eq_soil_moisture) then
            var_meta%name        = "eq_soil_moisture"
            var_meta%dimensions  = three_d_t_soil_dimensions
            var_meta%dim_len(2)  = kSOIL_GRID_Z
            var_meta%attributes  = [attribute_t("non_standard_name", "equilibrium_volumetric_soil_moisture"), &
                               attribute_t("units",         "m3 m-3"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Soil Moisture Content in the Layer Draining to Water Table when Deep
        !!------------------------------------------------------------
        else if (var_idx==kVARS%smc_watertable_deep) then
            var_meta%name        = "smc_watertable_deep"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "soil_moisture_content_in_layer_to_water_table_when_deep"), &
                               attribute_t("units",         "m3 m-3"),                                                      & !units not defined in noahmpdrv (guess) then
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Groundwater Recharge
        !!------------------------------------------------------------
        else if (var_idx==kVARS%recharge) then
            var_meta%name        = "recharge"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "groundwater_recharge"),            &
                               attribute_t("units",         "mm"),                                  & !units not defined in noahmpdrv (guess) then
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Groundwater Recharge when Water Table is Deep
        !!------------------------------------------------------------
        else if (var_idx==kVARS%recharge_deep) then
            var_meta%name        = "recharge_deep"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "groundwater_recharge_deep"),           &
                               attribute_t("units",         "mm"),                                  & !units not defined in noahmpdrv (guess) then
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Canopy Evaporation Rate
        !!------------------------------------------------------------
        else if (var_idx==kVARS%evap_canopy) then
            var_meta%name        = "evap_canopy"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "water_evaporation_flux_from_canopy"),  &
                               attribute_t("units",         "mm s-1"),                              &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Soil Surface Evaporation Rate
        !!------------------------------------------------------------
        else if (var_idx==kVARS%evap_soil_surface) then
            var_meta%name        = "evap_soil_surface"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "water_evaporation_flux_from_soil"),    &
                               attribute_t("units",         "mm s-1"),                              &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Transpiration Rate
        !!------------------------------------------------------------
        else if (var_idx==kVARS%transpiration_rate) then
            var_meta%name        = "transpiration_rate"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "transpiration_rate"),              &
                               attribute_t("units",         "mm s-1"),                              &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Monin-Obukhov Length
        !!------------------------------------------------------------
        else if (var_idx==kVARS%mol) then
            var_meta%name        = "mol"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "Monin-Obukhov Length"),                    &
                               attribute_t("units",         "m"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        !>------------------------------------------------------------
        !!  Shear velocity
        !!------------------------------------------------------------
        else if (var_idx==kVARS%ustar) then
            var_meta%name        = "ustar"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "Shear velocity"),                    &
                                attribute_t("units",         "m/s"),                                   &
                                attribute_t("coordinates",   "lat lon")]
                    
        !>------------------------------------------------------------
        !!  similarity stability function for momentum
        !!------------------------------------------------------------
        else if (var_idx==kVARS%psim) then
            var_meta%name        = "psih"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "similarity stability function for momentum"),                    &
                                attribute_t("units",         "-"),                                   &
                                attribute_t("coordinates",   "lat lon")]
                            
        !>------------------------------------------------------------
        !!  similarity stability function for heat
        !!------------------------------------------------------------
        else if (var_idx==kVARS%psih) then
            var_meta%name        = "psih"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "similarity stability function for heat"),                    &
                                attribute_t("units",         "-"),                                   &
                                attribute_t("coordinates",   "lat lon")]
                                                
        !>------------------------------------------------------------
        !!  integrated stability function for momentum
        !!------------------------------------------------------------
        else if (var_idx==kVARS%fm) then
            var_meta%name        = "fm"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "integrated stability function for momentum"),                    &
                                attribute_t("units",         "-"),                                   &
                                attribute_t("coordinates",   "lat lon")]
                                                                    
        !>------------------------------------------------------------
        !!  integrated stability function for heat
        !!------------------------------------------------------------
        else if (var_idx==kVARS%fh) then
            var_meta%name        = "fh"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "integrated stability function for heat"),                    &
                                attribute_t("units",         "-"),                                   &
                                attribute_t("coordinates",   "lat lon")]
        
                                                                                                                                                                                                        
        !>------------------------------------------------------------
        !!  Sensible Heat Exchange Coefficient, Vegetated
        !!------------------------------------------------------------
        else if (var_idx==kVARS%ch_veg) then
            var_meta%name        = "ch_veg"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "ch_vegetated"),                    &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Sensible Heat Exchange Coefficient at 2m, Vegetated
        !!------------------------------------------------------------
        else if (var_idx==kVARS%ch_veg_2m) then
            var_meta%name        = "ch_veg_2m"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "ch_vegetated_2m"),                 &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Sensible Heat Exchange Coefficient, Bare
        !!------------------------------------------------------------
        else if (var_idx==kVARS%ch_bare) then
            var_meta%name        = "ch_bare"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "ch_bare_ground"),                  &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Sensible Heat Exchange Coefficient at 2m, Bare
        !!------------------------------------------------------------
        else if (var_idx==kVARS%ch_bare_2m) then
            var_meta%name        = "ch_bare_2m"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "ch_bare_2m"),                      &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Sensible Heat Exchange Coefficient, Under Canopy
        !!------------------------------------------------------------
        else if (var_idx==kVARS%ch_under_canopy) then
            var_meta%name        = "ch_under_canopy"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "ch_under_canopy"),                 &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Sensible Heat Exchange Coefficient, Leaf
        !!------------------------------------------------------------
        else if (var_idx==kVARS%ch_leaf) then
            var_meta%name        = "ch_leaf"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "ch_leaf"),                         &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Sensible Heat Flux, Vegetated Ground (+ to atmosphere)
        !!------------------------------------------------------------
        else if (var_idx==kVARS%sensible_heat_veg) then
            var_meta%name        = "sensible_heat_veg"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "sensible_heat_veg"),               &
                               attribute_t("units",         "W m-2"),                               &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Sensible Heat Flux, Bare Ground (+ to atmosphere)
        !!------------------------------------------------------------
        else if (var_idx==kVARS%sensible_heat_bare) then
            var_meta%name        = "sensible_heat_bare"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "sensible_heat_bare"),              &
                               attribute_t("units",         "W m-2"),                               &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Sensible Heat Flux, Canopy (+ to atmosphere)
        !!------------------------------------------------------------
        else if (var_idx==kVARS%sensible_heat_canopy) then
            var_meta%name        = "sensible_heat_canopy"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "sensible_heat_canopy"),            &
                               attribute_t("units",         "W m-2"),                               &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Evaporation Heat Flux, Vegetated Ground (+ to atmosphere)
        !!------------------------------------------------------------
        else if (var_idx==kVARS%evap_heat_veg) then
            var_meta%name        = "evap_heat_veg"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "evap_heat_veg"),                   &
                               attribute_t("units",         "W m-2"),                               &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Evaporation Heat Flux, Bare Ground (+ to atmosphere)
        !!------------------------------------------------------------
        else if (var_idx==kVARS%evap_heat_bare) then
            var_meta%name        = "evap_heat_bare"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "evap_heat_bare"),                  &
                               attribute_t("units",         "W m-2"),                               &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !! Evaporation Heat Flux, Canopy (+ to atmosphere)
        !!------------------------------------------------------------
        else if (var_idx==kVARS%evap_heat_canopy) then
            var_meta%name        = "evap_heat_canopy"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "evap_heat_canopy"),                &
                               attribute_t("units",         "W m-2"),                               &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Transpiration Heat Flux (+ to atmosphere)
        !!------------------------------------------------------------
        else if (var_idx==kVARS%transpiration_heat) then
            var_meta%name        = "transpiration_heat"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "transpiration_heat"),              &
                               attribute_t("units",         "W m-2"),                               &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Ground Heat Flux, Vegetated Ground (+ to soil)
        !!------------------------------------------------------------
        else if (var_idx==kVARS%ground_heat_veg) then
            var_meta%name        = "ground_heat_veg"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "ground_heat_veg"),                 &
                               attribute_t("units",         "W m-2"),                               &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Ground Heat Flux, Bare Ground (+ to soil)
        !!------------------------------------------------------------
        else if (var_idx==kVARS%ground_heat_bare) then
            var_meta%name        = "ground_heat_bare"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "ground_heat_bare"),                &
                               attribute_t("units",         "W m-2"),                               &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Net Longwave Radiation Flux, Vegetated Ground (+ to atmosphere)
        !!------------------------------------------------------------
        else if (var_idx==kVARS%net_longwave_veg) then
            var_meta%name        = "net_longwave_veg"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "net_longwave_veg"),                &
                               attribute_t("units",         "W m-2"),                               &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Net Longwave Radiation Flux, Bare Ground (+ to atmosphere)
        !!------------------------------------------------------------
        else if (var_idx==kVARS%net_longwave_bare) then
            var_meta%name        = "net_longwave_bare"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "net_longwave_bare"),               &
                               attribute_t("units",         "W m-2"),                               &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Net Longwave Radiation Flux, Canopy (+ to atmosphere)
        !!------------------------------------------------------------
        else if (var_idx==kVARS%net_longwave_canopy) then
            var_meta%name        = "net_longwave_canopy"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "net_longwave_canopy"),             &
                               attribute_t("units",         "W m-2"),                               &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Surface Runoff Rate
        !!------------------------------------------------------------
        else if (var_idx==kVARS%runoff_surface) then
            var_meta%name        = "runoff_surface"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "surface_runoff_flux"),                 &
                               attribute_t("units",         "mm s-1"),                              &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Subsurface Runoff Rate
        !!------------------------------------------------------------
        else if (var_idx==kVARS%runoff_subsurface) then
            var_meta%name        = "runoff_subsurface"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "subsurface_runoff_flux"),          &
                               attribute_t("units",         "mm s-1"),                              &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Total Column Soil water content
        !!------------------------------------------------------------
        else if (var_idx==kVARS%soil_totalmoisture) then
            var_meta%name        = "soil_column_total_water"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "soil_moisture_content"),               &
                               attribute_t("units",         "kg m-2"),                              &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Soil Temperature
        !!------------------------------------------------------------
        else if (var_idx==kVARS%soil_temperature) then
            var_meta%name        = "soil_temperature"
            var_meta%dimensions  = three_d_t_soil_dimensions
            var_meta%dim_len(2)  = kSOIL_GRID_Z
            var_meta%attributes  = [attribute_t("standard_name", "soil_temperature"),                    &
                               attribute_t("units",         "K"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Deep Soil Temperature (time constant)
        !!------------------------------------------------------------
        else if (var_idx==kVARS%soil_deep_temperature) then
            var_meta%name        = "soil_deep_temperature"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "deep_soil_temperature"),           &
                               attribute_t("units",         "K"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Stable Carbon Mass in Deep Soil
        !!------------------------------------------------------------
        else if (var_idx==kVARS%soil_carbon_stable) then
            var_meta%name        = "soil_carbon_stable"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "slow_soil_pool_mass_content_of_carbon"), &
                               attribute_t("units",         "g m-2"),                                 &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Short-lived Carbon Mass in Shallow Soil
        !!------------------------------------------------------------
        else if (var_idx==kVARS%soil_carbon_fast) then
            var_meta%name        = "soil_carbon_fast"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "fast_soil_pool_mass_content_of_carbon"), &
                               attribute_t("units",         "g m-2"),                                 &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Soil Class, Layer 1
        !!------------------------------------------------------------
        else if (var_idx==kVARS%soil_texture_1) then
            var_meta%name        = "soil_class_1"
            var_meta%dimensions  = two_d_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "soil_class_layer1"),                 &
                               attribute_t("units",         "1"),                                     &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Soil Class, Layer 2
        !!------------------------------------------------------------
        else if (var_idx==kVARS%soil_texture_2) then
            var_meta%name        = "soil_class_2"
            var_meta%dimensions  = two_d_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "soil_class_layer2"),                 &
                               attribute_t("units",         "1"),                                     &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Soil Class, Layer 3
        !!------------------------------------------------------------
        else if (var_idx==kVARS%soil_texture_3) then
            var_meta%name        = "soil_class_3"
            var_meta%dimensions  = two_d_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "soil_class_layer3"),                 &
                               attribute_t("units",         "1"),                                     &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Soil Class, Layer 4
        !!------------------------------------------------------------
        else if (var_idx==kVARS%soil_texture_4) then
            var_meta%name        = "soil_class_4"
            var_meta%dimensions  = two_d_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "soil_class_layer4"),                 &
                               attribute_t("units",         "1"),                                     &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Soil Sand and Clay Composition by Layer
        !!------------------------------------------------------------
        else if (var_idx==kVARS%soil_sand_and_clay) then
            var_meta%name        = "soil_sand_and_clay_composition"
            var_meta%dimensions  = three_d_soilcomp_dimensions
            var_meta%dim_len(2)  = kSOILCOMP_GRID_Z
            var_meta%attributes  = [attribute_t("non_standard_name", "soil_sand_and_clay_composition"),    &
                               attribute_t("units",         "1"),                                     &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Water Table Depth
        !!------------------------------------------------------------
        else if (var_idx==kVARS%water_table_depth) then
            var_meta%name        = "water_table_depth"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "water_table_depth"),                    &
                               attribute_t("units",         "m"),                                    &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Water in Aquifer
        !!------------------------------------------------------------
        else if (var_idx==kVARS%water_aquifer) then
            var_meta%name        = "water_aquifer"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "water_in_aquifer"),                 &
                               attribute_t("units",         "mm"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Groundwater Storage
        !!------------------------------------------------------------
        else if (var_idx==kVARS%storage_gw) then
            var_meta%name        = "storage_gw"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "groundwater_storage"),              &
                               attribute_t("units",         "mm"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Lake Storage
        !!------------------------------------------------------------
        else if (var_idx==kVARS%storage_lake) then
            var_meta%name        = "storage_lake"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "lake_storage"),                     &
                               attribute_t("units",         "mm"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Surface roughness length z0
        !!------------------------------------------------------------
        else if (var_idx==kVARS%roughness_z0) then
            var_meta%name        = "surface_roughness"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "surface_roughness_length"),            &
                               attribute_t("long_name",     "Surface roughness length"),            &
                               attribute_t("units",         "m"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Surface Radiative Temperature
        !!------------------------------------------------------------
        else if (var_idx==kVARS%surface_rad_temperature) then
            var_meta%name        = "surface_rad_temperature"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "surface_radiative_temperature"),       &
                               attribute_t("units",         "K"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  2 meter air temperture
        !!------------------------------------------------------------
        else if (var_idx==kVARS%temperature_2m) then
            var_meta%name        = "taix"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "air_temperature"),                     &
                               attribute_t("long_name",     "Bulk air temperature at 2m"),          &
                               attribute_t("units",         "K"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  2 meter Air Temperature over Vegetation
        !!------------------------------------------------------------
        else if (var_idx==kVARS%temperature_2m_veg) then
            var_meta%name        = "temperature_2m_veg"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "air_temperature"),                     &
                               attribute_t("units",         "K"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  2 meter Air Temperature over Bare Ground
        !!------------------------------------------------------------
        else if (var_idx==kVARS%temperature_2m_bare) then
            var_meta%name        = "temperature_2m_bare"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "air_temperature"),                     &
                               attribute_t("units",         "K"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  2 meter Mixing Ratio over Vegetation
        !!------------------------------------------------------------
        else if (var_idx==kVARS%mixing_ratio_2m_veg) then
            var_meta%name        = "mixing_ratio_2m_veg"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "mixing_ratio"),                    &
                               attribute_t("units",         "kg kg-1"),                             &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  2 meter Mixing Ratio over Bare Ground
        !!------------------------------------------------------------
        else if (var_idx==kVARS%mixing_ratio_2m_bare) then
            var_meta%name        = "mixing_ratio_2m_bare"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "mixing_ratio"),                    &
                               attribute_t("units",         "kg kg-1"),                             &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  2 meter specific humidity
        !!------------------------------------------------------------
        else if (var_idx==kVARS%humidity_2m) then
            var_meta%name        = "hus2m"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "specific_humidity"),                   &
                               attribute_t("units",         "kg kg-2"),                             &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  10 meter height V component of wind field
        !!------------------------------------------------------------
        else if (var_idx==kVARS%v_10m) then
            var_meta%name        = "v10m"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "northward_10m_wind_speed"),            &
                               attribute_t("units",         "m s-1"),                               &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  10 meter height U component of the wind field
        !!------------------------------------------------------------
        else if (var_idx==kVARS%u_10m) then
            var_meta%name        = "u10m"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "eastward_10m_wind_speed"),             &
                               attribute_t("units",         "m s-1"),                               &
                               attribute_t("coordinates",   "lat lon")]
        !>------------------------------------------------------------
        !!  V component of wind field averaged to the mass grid
        !!------------------------------------------------------------
        else if (var_idx==kVARS%v_mass) then
            var_meta%name        = "v_mass"
            var_meta%dimensions  = three_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "northward wind, averaged to the mass grid"),            &
                                attribute_t("units",         "m s-1"),                               &
                                attribute_t("coordinates",   "lat lon")]
                            
        !>------------------------------------------------------------
        !!  U component of the wind field averaged to the mass grid
        !!------------------------------------------------------------
        else if (var_idx==kVARS%u_mass) then
            var_meta%name        = "u_mass"
            var_meta%dimensions  = three_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "eastward_ wind, averaged to the mass grid"),             &
                                attribute_t("units",         "m s-1"),                               &
                                attribute_t("coordinates",   "lat lon")]
                    
        !>------------------------------------------------------------
        !!  10 meter height wind speed magnitude, sqrt(u_10m**2+v_10m**2)
        !!------------------------------------------------------------
        else if (var_idx==kVARS%windspd_10m) then
            var_meta%name        = "wnsx"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "wind speed magnitude"),             &
                               attribute_t("units",         "m s-1"),                               &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Momentum Drag Coefficient
        !!------------------------------------------------------------
        else if (var_idx==kVARS%coeff_momentum_drag) then
            var_meta%name        = "coeff_momentum_drag"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "surface_drag_coefficient_for_momentum_in_air"), &
                               attribute_t("units",         "1"),                                            &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Sensible Heat Exchange Coefficient
        !!------------------------------------------------------------
        else if (var_idx==kVARS%chs) then
            var_meta%name        = "coeff_heat_exchange"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "sensible_heat_exchange_coefficient"), &
                               attribute_t("units",         "1"),                                      &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Sensible Heat Exchange Coefficient, 2m
        !!------------------------------------------------------------
        else if (var_idx==kVARS%chs2) then
            var_meta%name        = "coeff_heat_exchange_2m"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "sensible_heat_exchange_coefficient_at_2m"), &
                               attribute_t("units",         "1"),                                      &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Latent Heat Exchange Coefficient, 2m
        !!------------------------------------------------------------
        else if (var_idx==kVARS%cqs2) then
            var_meta%name        = "coeff_moisture_exchange"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "latent_heat_exchange_coefficient"), &
                               attribute_t("units",         "1"),                                      &
                               attribute_t("coordinates",   "lat lon")]
        

        !>------------------------------------------------------------
        !!  Sensible Heat Exchange Coefficient 3d
        !!------------------------------------------------------------
        else if (var_idx==kVARS%coeff_heat_exchange_3d) then
            var_meta%name        = "coeff_heat_exchange_3d"
            var_meta%dimensions  = three_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "sensible_heat_exchange_coefficient_3d"), &
                               attribute_t("units",         "1"),                                      &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Momentum Exchange Coefficient 3d
        !!------------------------------------------------------------
        else if (var_idx==kVARS%coeff_momentum_exchange_3d) then
            var_meta%name        = "coeff_momentum_exchange_3d"
            var_meta%dimensions  = three_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "momentum_exchange_coefficient_3d"), &
                               attribute_t("units",         "1"),                                      &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Sensible Heat Exchange Coefficient
        !!------------------------------------------------------------
        else if (var_idx==kVARS%br) then
            var_meta%name        = "sfc_Ri"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "bulk_richardson_num_from_SFC_scheme"), &
                               attribute_t("units",         "1"),                                      &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Sensible Heat Exchange Coefficient
        !!------------------------------------------------------------
        else if (var_idx==kVARS%QFX) then
            var_meta%name        = "moisture_flux"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "moisture_flux_from_surface"), &
                               attribute_t("units",         "kg/m²/s"),                                      &
                               attribute_t("coordinates",   "lat lon")]
        

        !>------------------------------------------------------------
        !!  PBL height
        !!------------------------------------------------------------
        else if (var_idx==kVARS%hpbl) then
            var_meta%name        = "hpbl"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "height_of_planetary_boundary_layer"), &
                               attribute_t("units",         "m"),                                      &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  PBL layer index
        !!------------------------------------------------------------
        else if (var_idx==kVARS%kpbl) then
            var_meta%name        = "kpbl"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "index_of_planetary_boundary_layer_height"), &
                               attribute_t("units",         "-"),                                      &
                               attribute_t("coordinates",   "lat lon")]
            var_meta%dtype = kINTEGER
        
        !>------------------------------------------------------------
        !!  Land surface radiative skin temperature
        !!------------------------------------------------------------
        else if (var_idx==kVARS%skin_temperature) then
            var_meta%name        = "tsfe"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "surface_temperature"),                 &
                               attribute_t("units",         "K"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Land surface radiative skin temperature
        !!------------------------------------------------------------
        else if (var_idx==kVARS%sst) then
            var_meta%name        = "sst"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "Sea Surface Temperature"),                 &
                               attribute_t("units",         "K"),                                   &
                               attribute_t("coordinates",   "lat lon")]
            if (present(force_boundaries)) force_boundaries = .False.
            if (present(opt) .and. present(forcing_var)) forcing_var = (opt%forcing%sst_var /= "")

        !>------------------------------------------------------------
        !!  Sensible heat flux from the surface (positive up)
        !!------------------------------------------------------------
        else if (var_idx==kVARS%sensible_heat) then
            var_meta%name        = "hfss"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "surface_upward_sensible_heat_flux"),   &
                               attribute_t("units",         "W m-2"),                               &
                               attribute_t("coordinates",   "lat lon")]
            if (present(force_boundaries)) force_boundaries = .False.
            if (present(opt) .and. present(forcing_var)) forcing_var = (opt%forcing%shvar /= "")

        !>------------------------------------------------------------
        !!  Latent heat flux from the surface (positive up)
        !!------------------------------------------------------------
        else if (var_idx==kVARS%latent_heat) then
            var_meta%name        = "hfls"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "surface_upward_latent_heat_flux"),     &
                               attribute_t("units",         "W m-2"),                               &
                               attribute_t("coordinates",   "lat lon")]
            if (present(force_boundaries)) force_boundaries = .False.
            if (present(opt) .and. present(forcing_var)) forcing_var = (opt%forcing%lhvar /= "")

        !>------------------------------------------------------------
        !!  Lake temperature 3d
        !!------------------------------------------------------------
        else if (var_idx==kVARS%t_lake3d) then
            var_meta%name        = "t_lake3d"
            var_meta%dimensions  = three_d_t_lake_dimensions
            var_meta%dim_len(2)  = kLAKE_Z
            var_meta%attributes  = [attribute_t("standard_name", "lake_water_temperature"),     &
                               attribute_t("units",         "K"),                               &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Lake lake_icefraction_3d
        !!------------------------------------------------------------
        else if (var_idx==kVARS%lake_icefrac3d) then
            var_meta%name        = "lake_icefrac3d"
            var_meta%dimensions  = three_d_t_lake_dimensions
            var_meta%dim_len(2)  = kLAKE_Z
            var_meta%attributes  = [attribute_t("standard_name", "lake_icefraction_3d"),     &
                               attribute_t("units",         "-"),                               &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Lake z_lake3d
        !!------------------------------------------------------------
        else if (var_idx==kVARS%z_lake3d) then
            var_meta%name        = "z_lake3d"
            var_meta%dimensions  = three_d_t_lake_dimensions
            var_meta%dim_len(2)  = kLAKE_Z
            var_meta%attributes  = [attribute_t("standard_name", "lake_layer_depth"),     &
                               attribute_t("units",         "m"),                               &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Lake dz_lake3d
        !!------------------------------------------------------------
        else if (var_idx==kVARS%dz_lake3d) then
            var_meta%name        = "dz_lake3d"
            var_meta%dimensions  = three_d_t_lake_dimensions
            var_meta%dim_len(2)  = kLAKE_Z
            var_meta%attributes  = [attribute_t("standard_name", "lake_layer_thickness"),     &
                               attribute_t("units",         "m"),                               &
                               attribute_t("coordinates",   "lat lon")]
        !>------------------------------------------------------------
        !!  lake_depth
        !!------------------------------------------------------------
        else if (var_idx==kVARS%lake_depth) then
            var_meta%name        = "lake_depth"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "lake_depth"),     &
                                attribute_t("units",         "m"),                               &
                                attribute_t("coordinates",   "lat lon")]
                    
        !>------------------------------------------------------------
        !!  lake snl2d
        !!------------------------------------------------------------
        else if (var_idx==kVARS%snl2d) then
            var_meta%name        = "snl2d"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "lake_snow_layer_2d"),           &
                               attribute_t("units",         "-"),                               &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  lake_t_grnd2d
        !!------------------------------------------------------------
        else if (var_idx==kVARS%t_grnd2d) then
            var_meta%name        = "t_grnd2d"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "t_grnd2d"),           &
                               attribute_t("units",         "K"),                               &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Lake t_soisno3d
        !!------------------------------------------------------------
        else if (var_idx==kVARS%t_soisno3d) then
            var_meta%name        = "t_soisno3d"
            var_meta%dimensions  = three_d_t_lake_soisno_dimensions
            var_meta%dim_len(2)  = kLAKE_SOISNO_Z
            var_meta%attributes  = [attribute_t("standard_name", "temperature_soil_snow_below_or_above_lake"),     &
                               attribute_t("units",         "K"),                               &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Lake h2osoi_ice3d
        !!------------------------------------------------------------
        else if (var_idx==kVARS%h2osoi_ice3d) then
            var_meta%name        = "h2osoi_ice3d"
            var_meta%dimensions  = three_d_t_lake_soisno_dimensions
            var_meta%dim_len(2)  = kLAKE_SOISNO_Z
            var_meta%attributes  = [attribute_t("standard_name", "h2osoi_ice3d"),     &
                               attribute_t("units",         ""),                               &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Lake soil/snowliquid water (kg/m2)
        !!------------------------------------------------------------
        else if (var_idx==kVARS%h2osoi_liq3d) then
            var_meta%name        = "h2osoi_liq3d"
            var_meta%dimensions  = three_d_t_lake_soisno_dimensions
            var_meta%dim_len(2)  = kLAKE_SOISNO_Z
            var_meta%attributes  = [attribute_t("standard_name", "lake_soil_or_snow_liquid water_content"),     &
                               attribute_t("units",         "kg/m2"),                               &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Lake h2osoi_vol3d volumetric soil water (0<=h2osoi_vol<=watsat)[m3/m3]
        !!------------------------------------------------------------
        else if (var_idx==kVARS%h2osoi_vol3d) then
            var_meta%name        = "h2osoi_vol3d"
            var_meta%dimensions  = three_d_t_lake_soisno_dimensions
            var_meta%dim_len(2)  = kLAKE_SOISNO_Z
            var_meta%attributes  = [attribute_t("standard_name", "volumetric_soil_water"),     &
                               attribute_t("units",         "m3/m3"),                               &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Lake z3d
        !!------------------------------------------------------------
        else if (var_idx==kVARS%z3d) then
            var_meta%name        = "z3d"
            var_meta%dimensions  = three_d_t_lake_soisno_dimensions
            var_meta%dim_len(2)  = kLAKE_SOISNO_Z
            var_meta%attributes  = [attribute_t("standard_name", "layer_depth_for_lake_snow&soil"),     &
                               attribute_t("units",         "m"),                               &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Lake layer_thickness_for_lake_snow&soil
        !!------------------------------------------------------------
        else if (var_idx==kVARS%dz3d) then
            var_meta%name        = "dz3d"
            var_meta%dimensions  = three_d_t_lake_soisno_dimensions
            var_meta%dim_len(2)  = kLAKE_SOISNO_Z
            var_meta%attributes  = [attribute_t("standard_name", "layer_thickness_for_lake_snow&soil"),     &
                               attribute_t("units",         "m"),                               &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Lake z3d
        !!------------------------------------------------------------
        else if (var_idx==kVARS%zi3d) then
            var_meta%name        = "zi3d"
            var_meta%dimensions  = three_d_t_lake_soisno_1_dimensions
            var_meta%dim_len(2)  = kLAKE_SOISNO_1_Z
            var_meta%attributes  = [attribute_t("standard_name", "interface_layer_depth_for_lake_snow&soil"),     &
                               attribute_t("units",         "m"),                               &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Lake watsat3d: volumetric soil water at saturation (porosity)
        !!------------------------------------------------------------
        else if (var_idx==kVARS%watsat3d) then
            var_meta%name        = "watsat3d"
            var_meta%dimensions  = three_d_t_lake_soi_dimensions
            var_meta%dim_len(2)  = kLAKE_SOI_Z
            var_meta%attributes  = [attribute_t("standard_name", "volumetric soil water at saturation (porosity)"),     &
                               attribute_t("units",         ""),                               &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Lake csol3d: heat capacity, soil solids (J/m**3/Kelvin)
        !!------------------------------------------------------------
        else if (var_idx==kVARS%csol3d) then
            var_meta%name        = "csol3d"
            var_meta%dimensions  = three_d_t_lake_soi_dimensions
            var_meta%dim_len(2)  = kLAKE_SOI_Z
            var_meta%attributes  = [attribute_t("standard_name", "heat capacity, soil solids "),     &
                               attribute_t("units",         "(J/m**3/Kelvin)"),                               &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Lake: thermal conductivity, soil minerals  [W/m-K]
        !!------------------------------------------------------------
        else if (var_idx==kVARS%tkmg3d) then
            var_meta%name        = "tkmg3d"
            var_meta%dimensions  = three_d_t_lake_soi_dimensions
            var_meta%dim_len(2)  = kLAKE_SOI_Z
            var_meta%attributes  = [attribute_t("standard_name", "thermal conductivity, soil minerals  [W/m-K]"),     &
                               attribute_t("units",         ""),                               &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Lake lakemask
        !!------------------------------------------------------------
        else if (var_idx==kVARS%lakemask) then
            var_meta%name        = "lakemask"
            var_meta%dimensions  = two_d_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "lakemask"),     &
                               attribute_t("units",         ""),                               &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Grid cell fractional sea ice lakemask
        !!------------------------------------------------------------
        else if (var_idx==kVARS%xice) then
            var_meta%name        = "xice"
            var_meta%dimensions  = two_d_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "xice"),     &
                               attribute_t("units",         ""),                               &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Lake lakedepth2d
        !!------------------------------------------------------------
        else if (var_idx==kVARS%lakedepth2d) then
            var_meta%name        = "lakedepth2d"
            var_meta%dimensions  = two_d_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "lake_depth"),     &
                               attribute_t("units",         "m"),                               &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  lake savedtke12d
        !!------------------------------------------------------------
        else if (var_idx==kVARS%savedtke12d) then
            var_meta%name        = "savedtke12d"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "savedtke12d"),           &
                               attribute_t("units",         "-?"),                               &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Lake: thermal conductivity, saturated soil [W/m-K]
        !!------------------------------------------------------------
        else if (var_idx==kVARS%tksatu3d) then
            var_meta%name        = "tksatu3d"
            var_meta%dimensions  = three_d_t_lake_soi_dimensions
            var_meta%dim_len(2)  = kLAKE_SOI_Z
            var_meta%attributes  = [attribute_t("standard_name", "thermal conductivity, saturated soil [W/m-K]"),     &
                               attribute_t("units",         ""),                               &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Lake tkdry3d: thermal conductivity, dry soil (W/m/Kelvin)
        !!------------------------------------------------------------
        else if (var_idx==kVARS%tkdry3d) then
            var_meta%name        = "tkdry3d"
            var_meta%dimensions  = three_d_t_lake_soi_dimensions
            var_meta%dim_len(2)  = kLAKE_SOI_Z
            var_meta%attributes  = [attribute_t("standard_name", "thermal conductivity, dry soil (W/m/Kelvin)"),     &
                               attribute_t("units",         "?"),                               &
                               attribute_t("coordinates",   "lat lon")]
        



        !>------------------------------------------------------------
        !!  Integrated Vapor Transport
        !!------------------------------------------------------------
        else if (var_idx==kVARS%ivt) then
            var_meta%name        = "ivt"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "integrated_vapor_transport"),  &
                               attribute_t("units",         "kg m-1 s-1"),                      &
                               attribute_t("coordinates",   "lat lon")]
        

        !>------------------------------------------------------------
        !!  Integrated Water Vapor
        !!------------------------------------------------------------
        else if (var_idx==kVARS%iwv) then
            var_meta%name        = "iwv"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "atmosphere_mass_content_of_water_vapor"),  &
                               attribute_t("units",         "kg m-2"),                      &
                               attribute_t("coordinates",   "lat lon")]
        

        !>------------------------------------------------------------
        !!  Integrated Water Liquid
        !!------------------------------------------------------------
        else if (var_idx==kVARS%iwl) then
            var_meta%name        = "iwl"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "atmosphere_mass_content_of_water_liquid"),  &
                               attribute_t("units",         "kg m-2"),                      &
                               attribute_t("coordinates",   "lat lon")]
        

        !>------------------------------------------------------------
        !!  Integrated Water Ice
        !!------------------------------------------------------------
        else if (var_idx==kVARS%iwi) then
            var_meta%name        = "iwi"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "atmosphere_mass_content_of_water_ice"),  &
                               attribute_t("units",         "kg m-2"),                      &
                               attribute_t("coordinates",   "lat lon")]
        

        !>------------------------------------------------------------
        !!  Binary land mask (water vs land)
        !!------------------------------------------------------------
        else if (var_idx==kVARS%land_mask) then
            var_meta%name        = "land_mask"
            var_meta%dimensions  = two_d_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "land_water_mask"),                 &
                               attribute_t("coordinates",       "lat lon")]
            var_meta%dtype = kINTEGER

        !>------------------------------------------------------------
        !!  Height of the terrain
        !!------------------------------------------------------------
        else if (var_idx==kVARS%terrain) then
            var_meta%name        = "terrain"
            var_meta%dimensions  = two_d_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "height_above_reference_ellipsoid"),    &
                               attribute_t("units",         "m"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  Latitude y coordinate
        !!------------------------------------------------------------
        else if (var_idx==kVARS%latitude) then
            var_meta%name        = "lat"
            var_meta%dimensions  = two_d_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "latitude"),                            &
                               attribute_t("units",         "degrees_north"),                       &
                               attribute_t("axis","Y")]
        
        !>------------------------------------------------------------
        !!  Longitude x coordinate
        !!------------------------------------------------------------
        else if (var_idx==kVARS%longitude) then
            var_meta%name        = "lon"
            var_meta%dimensions  = two_d_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "longitude"),                           &
                               attribute_t("units",         "degrees_east"),                        &
                               attribute_t("axis","X")]
        
        !>------------------------------------------------------------
        !!  Latitude y coordinate on the U-grid
        !!------------------------------------------------------------
        else if (var_idx==kVARS%u_latitude) then
            var_meta%name        = "u_lat"
            var_meta%dimensions  = two_d_u_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "latitude_on_u_grid"),              &
                               attribute_t("units",         "degrees_north"),                       &
                               attribute_t("axis","Y")]

        !>------------------------------------------------------------
        !!  Longitude x coordinate on the U-grid
        !!------------------------------------------------------------
        else if (var_idx==kVARS%u_longitude) then
            var_meta%name        = "u_lon"
            var_meta%dimensions  = two_d_u_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "longitude_on_u_grid"),             &
                               attribute_t("units",         "degrees_east"),                        &
                               attribute_t("axis","X")]

        !>------------------------------------------------------------
        !!  Latitude y coordinate on the V-grid
        !!------------------------------------------------------------
        else if (var_idx==kVARS%v_latitude) then
            var_meta%name        = "v_lat"
            var_meta%dimensions  = two_d_v_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "latitude_on_v_grid"),              &
                               attribute_t("units",         "degrees_north"),                       &
                               attribute_t("axis","Y")]

        !>------------------------------------------------------------
        !!  Longitude x coordinate on the V-grid
        !!------------------------------------------------------------
        else if (var_idx==kVARS%v_longitude) then
            var_meta%name        = "v_lon"
            var_meta%dimensions  = two_d_v_dimensions
            var_meta%attributes  = [attribute_t("non_standard_name", "longitude_on_v_grid"),             &
                               attribute_t("units",         "degrees_east"),                        &
                               attribute_t("axis","X")]


        !! MJ added for needed new vars for FSM
        !>------------------------------------------------------------
        !!  Snow Temperature
        !!------------------------------------------------------------
        else if (var_idx==kVARS%Tsnow) then
            var_meta%name        = "Tsnow"
            var_meta%dimensions  = three_d_t_snow_dimensions
            var_meta%dim_len(2)  = kSNOW_GRID_Z
            var_meta%attributes  = [attribute_t("standard_name", "snow_temperature"),                    &
                               attribute_t("units",         "K"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  ice mass content of snow 
        !!------------------------------------------------------------
        else if (var_idx==kVARS%Sice) then
            var_meta%name        = "Sice"
            var_meta%dimensions  = three_d_t_snow_dimensions
            var_meta%dim_len(2)  = kSNOW_GRID_Z
            var_meta%attributes  = [attribute_t("standard_name", "snow_ice_content"),                    &
                               attribute_t("units",         "kg m-2"),                                   &
                               attribute_t("coordinates",   "lat lon")]
                       		
        !>------------------------------------------------------------
        !!  water mass content of snow 
        !!------------------------------------------------------------
        else if (var_idx==kVARS%Sliq) then
            var_meta%name        = "Sliq"
            var_meta%dimensions  = three_d_t_snow_dimensions
            var_meta%dim_len(2)  = kSNOW_GRID_Z
            var_meta%attributes  = [attribute_t("standard_name", "snow_water_content"),                    &
                               attribute_t("units",         "kg m-2"),                                   &
                               attribute_t("coordinates",   "lat lon")]
                		
        !>------------------------------------------------------------
        !!  snow layer thicknes 
        !!------------------------------------------------------------
        else if (var_idx==kVARS%Ds) then
            var_meta%name        = "Ds"
            var_meta%dimensions  = three_d_t_snow_dimensions
            var_meta%dim_len(2)  = kSNOW_GRID_Z
            var_meta%attributes  = [attribute_t("standard_name", "snow_layer_thickness"),                    &
                               attribute_t("units",         "m"),                                   &
                               attribute_t("coordinates",   "lat lon")]
           
        !>------------------------------------------------------------
        !!  fsnow from FSM
        !!------------------------------------------------------------
        else if (var_idx==kVARS%fsnow) then
            var_meta%name        = "scfe"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "snow covered fraction"),   &
                               attribute_t("units",         "-"),                        &
                               attribute_t("coordinates",   "lat lon")]
        	
        !>------------------------------------------------------------
        !!  Nsnow from FSM
        !!------------------------------------------------------------
        else if (var_idx==kVARS%Nsnow) then
            var_meta%name        = "Nsnow"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "Nsnow"),   &
                               attribute_t("units",         "-"),                        &
                               attribute_t("coordinates",   "lat lon")]
        
        !>------------------------------------------------------------
        !!  runoff_tstep from FSM
        !!------------------------------------------------------------
        else if (var_idx==kVARS%runoff_tstep) then
            var_meta%name        = "rotc"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "total runoff, aggregated per output interval during t-1->t"),   &
                               attribute_t("units",         "kg m-2"),                        &
                               attribute_t("coordinates",   "lat lon")]
          
        !>------------------------------------------------------------
        !!  meltflux_out_tstep from FSM
        !!------------------------------------------------------------
        else if (var_idx==kVARS%meltflux_out_tstep) then
            var_meta%name        = "romc"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "runoff from snow, aggregated per output interval during t-1->t"),   &
                               attribute_t("units",         "kg m-2"),                        &
                               attribute_t("coordinates",   "lat lon")]
          
        !>------------------------------------------------------------
        !!  Sliq_out from FSM
        !!------------------------------------------------------------
        else if (var_idx==kVARS%Sliq_out) then
            var_meta%name        = "slqt"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "liquid water content"),   &
                               attribute_t("units",         "kg m-2"),                        &
                               attribute_t("coordinates",   "lat lon")]
          
        !>------------------------------------------------------------
        !!  saltation flux from FSM
        !!------------------------------------------------------------
        else if (var_idx==kVARS%dSWE_salt) then
            var_meta%name        = "salt"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "saltation flux"),   &
                               attribute_t("units",         "kg m-2 t-1"),                        &
                               attribute_t("coordinates",   "lat lon")]
          
        !>------------------------------------------------------------
        !!  suspension flux from FSM
        !!------------------------------------------------------------
        else if (var_idx==kVARS%dSWE_susp) then
            var_meta%name        = "susp"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "suspension flux"),   &
                               attribute_t("units",         "kg m-2 t-1"),                        &
                               attribute_t("coordinates",   "lat lon")]
          
        !>------------------------------------------------------------
        !!  blowing snow sublimation from FSM
        !!------------------------------------------------------------
        else if (var_idx==kVARS%dSWE_subl) then
            var_meta%name        = "blow_subl"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "blowing-snow sublimation"),   &
                               attribute_t("units",         "kg m-2 t-1"),                        &
                               attribute_t("coordinates",   "lat lon")]
          
        !>------------------------------------------------------------
        !!  sliding snow transport from FSM
        !!------------------------------------------------------------
        else if (var_idx==kVARS%dSWE_slide) then
            var_meta%name        = "slide"
            var_meta%dimensions  = two_d_t_dimensions
            var_meta%attributes  = [attribute_t("standard_name", "sliding snow transport"),   &
                               attribute_t("units",         "kg m-2 t-1"),                        &
                               attribute_t("coordinates",   "lat lon")]
        end if  

        ! loop through entire array setting n_dimensions and n_attrs based on the data that were supplied
        if (var_meta%name == "") return

        var_meta%n_attrs      = size(var_meta%attributes)
        var_meta%id = var_idx

        var_meta%unlimited_dim=.False.

        do j = 1,size(var_meta%dimensions)
            if (var_meta%dimensions(j) == "time") then
                var_meta%unlimited_dim=.True.
            endif
            if (var_meta%dimensions(j) == "lon_u") then
                var_meta%xstag = 1
            endif
            if (var_meta%dimensions(j) == "lat_v") then
                var_meta%ystag = 1
            endif
        enddo

        if (var_meta%unlimited_dim) then
                !If time is one of the dimensions, and we only have 3 dimensions, then this is a 2D variable
            if (size(var_meta%dimensions) == 2) then
                var_meta%one_d = .True.
            else if (size(var_meta%dimensions) == 3) then
                var_meta%two_d = .True.
            else if (size(var_meta%dimensions) == 4) then
                var_meta%three_d = .True.
            endif

        else
            if (size(var_meta%dimensions) == 1) then
                var_meta%one_d = .True.
            else if (size(var_meta%dimensions) == 2) then
                var_meta%two_d = .True.
            else if (size(var_meta%dimensions) == 3) then
                var_meta%three_d = .True.
            else if (size(var_meta%dimensions) == 4) then
                var_meta%four_d = .True.
            endif
        endif

    end function get_varmeta

end module output_metadata
