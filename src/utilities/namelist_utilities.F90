module namelist_utils
    use netcdf
    use string,           only  : str
    use ieee_arithmetic
    use icar_constants,    only : STD_OUT_PE, MAXLEVELS, kMAX_NESTS, kMAX_STRING_LENGTH, kMAX_FILE_LENGTH, kVERSION_STRING, &
                                  kREAL_NO_VAL, kINT_NO_VAL, kCHAR_NO_VAL
    use mod_wrf_constants, only : piconst
    use options_types,     only : forcing_options_type, domain_options_type
    use io_routines,       only : check_file_exists, check_variable_present, io_getdims, io_newunit

    implicit none

    character(len=kMAX_FILE_LENGTH), private :: nml_file = ""
    integer, private :: nml_unit = 0

    character(len=30), parameter :: BLNK_CHR_N = ' ' ! number of blank characters, used for formating variable output

    integer, parameter :: NML_BLNK_LEN = 50 ! column where the comment begins in the default .nml file
    integer, parameter :: VERBOSE_BLNK_LEN = 30 ! column where the description begins in the verbose output
    integer, parameter :: SCREEN_WIDTH = 100 ! width of the screen for formatting output

    interface set_nml_var_default
        module procedure set_real_nml_var_default, set_integer_nml_var_default, set_logical_nml_var_default, set_char_nml_var_default, &
            set_real_list_nml_var_default, set_integer_array_nml_var_default, set_logical_array_nml_var_default, &
            set_char_array_nml_var_default, set_char_list_array_nml_var_default, set_real_list_array_nml_var_default
    end interface

    interface set_nml_var
        module procedure set_real_nml_var, set_integer_nml_var, set_logical_nml_var, set_char_nml_var, &
            set_char_forcing_nml_var, set_char_domain_nml_var, set_real_list_nml_var
    end interface

contains

    subroutine set_namelist(namelist_file)
        implicit none
        character(len=*), intent(in) :: namelist_file

        nml_file = namelist_file

        ! Open the namelist file
        if (STD_OUT_PE) open(UNIT=io_newunit(nml_unit), FILE=nml_file, STATUS='NEW', ACTION='WRITE')

    end subroutine set_namelist

    subroutine set_real_nml_var(var, var_val, name, usr_default)
        implicit none
        real,    intent(inout) :: var
        real,    intent(in) :: var_val
        character(len=*), intent(in) :: name
        real, optional, intent(in) :: usr_default

        character(len=kMAX_STRING_LENGTH) :: default
        real :: minmax(2), default_val

        !convert the default value string to a real
        default = trim(get_nml_var_default(name))
        read(default,*) default_val

        if (present(usr_default)) then
            ! If the user set the first value, but not the rest, then they want this first value to apply to all nests
            if ( (usr_default /= kREAL_NO_VAL) .and. (var_val == kREAL_NO_VAL )) then
                var = usr_default
            elseif ( usr_default == kREAL_NO_VAL ) then
                ! If the user did not set the first value, then there is no entry at all, so set it to the default
                read(default,*) var
            else
                var = var_val
            endif
        ! If this is not a variable for which nested values can differ, then we just check if the value was not set in the namelist
        else if (var_val == kREAL_NO_VAL) then
            read(default,*) var
        else
            var = var_val
        endif

        minmax = get_nml_var_minmax(name)

        if (all(minmax==0)) then
            if (STD_OUT_PE) write(*,*) "Error: '", trim(name), "' does not have a min/max value"
            error stop
        endif
        if (.not.(var==default_val)) then
            if (var > minmax(2) .and. STD_OUT_PE) then
                write(*,*) "Error: '", trim(name), "' is greater than ", minmax(2), " : ", var
                stop
            endif
            if (var < minmax(1) .and. STD_OUT_PE) then
                write(*,*) "Error: '", trim(name), "' is less than ", minmax(1), " : ", var
                stop
            endif
        endif

        if (present(usr_default)) then
            if (usr_default /= kREAL_NO_VAL) then
                call check_numeric_var_entry_type(usr_default,var,name)
            else
                call check_numeric_var_entry_type(default_val,var,name)
            endif
        endif

        ! Check if the variable was not set. This would be an error
        if (var == kREAL_NO_VAL .and. default_val /= kREAL_NO_VAL) then
            if (STD_OUT_PE) write(*,*) "Error: '", trim(name), "' is not set to any value (i.e. is still ",trim(str(kREAL_NO_VAL)),")"
            error stop
        endif
    end subroutine set_real_nml_var

    subroutine set_real_list_nml_var(var, var_val, name, usr_default)
        implicit none
        real,    intent(inout) :: var(:)
        real,    intent(in) :: var_val(:)
        character(len=*), intent(in) :: name
        real, optional, intent(in) :: usr_default(:)

        character(len=kMAX_STRING_LENGTH) :: default
        real :: minmax(2), default_val

        default = trim(get_nml_var_default(name))
        read(default,*) default_val

        if (present(usr_default)) then
            ! If the user set the first value, but not the rest, then they want this first value to apply to all nests
            if ( ALL(usr_default /= kREAL_NO_VAL) .and. ALL(var_val == kREAL_NO_VAL )) then
                var = usr_default
            elseif ( ALL(usr_default == kREAL_NO_VAL) ) then
                ! If the user did not set the first value, then there is no entry at all, so set it to the default
                var = default_val
            else
                var = var_val
            endif
        ! If this is not a variable for which nested values can differ, then we just check if the value was not set in the namelist
        else
            var = var_val
            where (var == kREAL_NO_VAL)
                var = default_val
            end where
        endif

        minmax = get_nml_var_minmax(name)

        if (all(minmax==0)) then
            if (STD_OUT_PE) write(*,*) "Error: '", trim(name), "' does not have a min/max value"
            error stop
        endif
        if (size(var) == 0) then
            if (STD_OUT_PE) write(*,*) "Error: '", trim(name), "' is empty"
            error stop
        endif

        if (any(var < minmax(1))) then
            if (STD_OUT_PE) write(*,*) "Error: '", trim(name), "' is less than ", minmax(1), " : ", minval(var)
            error stop
        endif
        if (any(var > minmax(2))) then
            if (STD_OUT_PE) write(*,*) "Error: '", trim(name), "' is greater than ", minmax(2), " : ", maxval(var)
            error stop
        endif

        ! Check if the variable was not set. This would be an error
        if (ALL(var == kREAL_NO_VAL .and. default_val /= kREAL_NO_VAL)) then
            if (STD_OUT_PE) write(*,*) "Error: '", trim(name), "' is not set to any value (i.e. is still ",trim(str(kREAL_NO_VAL)),")"
            error stop
        endif

    end subroutine set_real_list_nml_var

    subroutine set_integer_nml_var(var, var_val, name, usr_default)
        implicit none
        integer,    intent(inout) :: var
        integer,    intent(in) :: var_val
        character(len=*), intent(in) :: name
        integer, optional, intent(in) :: usr_default

        integer :: minmax(2), default_val
        integer, allocatable :: values(:)
        character(len=kMAX_STRING_LENGTH) :: default

        default = trim(get_nml_var_default(name))
        read(default,*) default_val

        if (present(usr_default)) then
            ! If the user set the first value, but not the rest, then they want this first value to apply to all nests
            if ( (usr_default /= kINT_NO_VAL) .and. (var_val == kINT_NO_VAL )) then
                var = usr_default
            elseif ( usr_default == kINT_NO_VAL ) then
                ! If the user did not set the first value, then there is no entry at all, so set it to the default
                read(default,*) var
            else
                var = var_val
            endif
        ! If this is not a variable for which nested values can differ, then we just check if the value was not set in the namelist
        else if (var_val == kINT_NO_VAL) then
            read(default,*) var
        else
            var = var_val
        endif

        values = get_nml_var_values(name)
        minmax = get_nml_var_minmax(name)

        if (size(values) == 0) then
            if (all(minmax==0)) then
                if (STD_OUT_PE) write(*,*) "Error: '", trim(name), "' does not have a min/max value, or does not have a list of valid values"
                error stop
            endif
            if (var < minmax(1)) then
                if (STD_OUT_PE) write(*,*) "Error: '", trim(name), "' is less than ", minmax(1), " : ", var
                error stop
            endif
            if (var > minmax(2)) then
                if (STD_OUT_PE) write(*,*) "Error: '", trim(name), "' is greater than ", minmax(2), " : ", var
                error stop
            endif
        else
            values = get_nml_var_values(name)

            if (all(values /= var)) then
                if (STD_OUT_PE) write(*,*) "Error: '", trim(name), "' is not in the list of valid values: ", values
                error stop
            endif
        endif

        if (present(usr_default)) then
            if (usr_default /= kINT_NO_VAL) then
                call check_numeric_var_entry_type((usr_default*1.0),(var*1.0),name)
            else
                call check_numeric_var_entry_type((default_val*1.0),(var*1.0),name)
            endif
        endif

        ! Check if the variable was not set. This would be an error
        if (var == kINT_NO_VAL .and. default_val /= kINT_NO_VAL) then
            if (STD_OUT_PE) write(*,*) "Error: '", trim(name), "' is not set to any value (i.e. is still ",trim(str(kREAL_NO_VAL)),")"
            error stop
        endif

    end subroutine set_integer_nml_var

    subroutine set_logical_nml_var(var, var_val, name, usr_default)
        implicit none
        logical,    intent(inout) :: var
        logical,    intent(in) :: var_val
        character(len=*), intent(in) :: name
        logical, optional, intent(in) :: usr_default

        integer :: type

        var = var_val

        type = get_nml_var_type(name)

        if (type==1 .and. present(usr_default)) then
            var = usr_default
        endif

    end subroutine set_logical_nml_var

    subroutine set_char_forcing_nml_var(var, var_val, name, options, i, no_check)
        implicit none
        character(len=*),    intent(inout) :: var
        character(len=*),    intent(in) :: var_val
        character(len=*), intent(in) :: name
        type(forcing_options_type), intent(inout) :: options
        integer, optional, intent(inout) :: i
        logical, optional :: no_check

        character(len=kMAX_STRING_LENGTH) :: group, default, description, units
        real :: min, max
        integer, allocatable :: values(:), var_dims(:), t_var_dims(:)
        integer :: p, dim_indx, dim_len, type
        character(len=kMAX_STRING_LENGTH) :: dim_name
        character(len=1), allocatable :: dimensions(:), t_dimensions(:), dimensions_tmp(:), t_dimensions_tmp(:)
        logical :: no_check_flag, is_latlon

        if (present(no_check)) then
            no_check_flag = no_check
        else
            no_check_flag = .False.
        endif

        call get_nml_var_metadata(name,group,description,default,min,max,type,values,units,dimensions)

        ! see if user wants to use this variable anyways
        if (.not.(trim(var_val) == "") .and. .not.(no_check_flag) .and. .not.(trim(var_val) == kCHAR_NO_VAL)) then
            ! first check if variable is present in the domain file
            call check_variable_present(options%boundary_files(1), var_val)
            ! then get the number of dimensions for the variable
            call io_getdims(options%boundary_files(1), var_val, var_dims)
            ! test if the number of dimensions matches the number of dimensions required
            if (.not.(size(dimensions) == size(var_dims))) then

                is_latlon =  ( (index(name,'lat') > 0) .or. (index(name,'lon') > 0))

                if ( (size(dimensions)-1) == size(var_dims) .and. (dimensions(1) == 'T' .or. is_latlon)) then
                    allocate(dimensions_tmp,source=dimensions)
                    deallocate(dimensions)
                    allocate(dimensions(size(dimensions_tmp)-1))

                    if (.not.is_latlon) then
                        if (STD_OUT_PE) write(*,*) "Forcing variable ", trim(name), " possibly missing the time dimension."
                        if (STD_OUT_PE) write(*,*) "Proceeding, assuming that the forcing data does not contain a time dimension."
                        dimensions = dimensions_tmp(2:size(dimensions_tmp))
                    else
                        if (STD_OUT_PE) write(*,*) "Forcing variable ", trim(name), " is 1D."

                        if (index(name,'lat') > 0) then
                            dimensions = dimensions_tmp(1:size(dimensions_tmp)-1)
                        else
                            dimensions = dimensions_tmp(2:size(dimensions_tmp))
                        endif
                    endif
                    deallocate(dimensions_tmp)
                else
                    if (STD_OUT_PE) write(*,*) "Error: ", trim(name), " has the wrong number of dimensions, should be: ", size(dimensions)
                    error stop
                endif
            endif

            ! pointless if this is the temperature variable, which is used as the reference
            if (.not.(name=='tvar')) then
                ! now get the dimension lengths of the temperature variable
                call check_variable_present(options%boundary_files(1), options%tvar)
                !! then get the number of dimensions for the variable
                call io_getdims(options%boundary_files(1), options%tvar, t_var_dims)
        
                call get_nml_var_metadata('tvar',group,description,default,min,max,type,values,units,t_dimensions)

                ! reverse ordering of dimensions -- needed since fortran netCDF uses reverse ordering relative to python/nco
                allocate(t_dimensions_tmp(size(t_dimensions)))
                t_dimensions_tmp = t_dimensions
                do p=1,size(t_dimensions)
                    t_dimensions(p) = t_dimensions_tmp(size(t_dimensions_tmp)-p+1)
                end do

                ! reverse ordering of dimensions -- needed since fortran netCDF uses reverse ordering relative to python/nco
                allocate(dimensions_tmp(size(dimensions)))
                dimensions_tmp = dimensions
                do p=1,size(dimensions)
                    dimensions(p) = dimensions_tmp(size(dimensions_tmp)-p+1)
                end do

                !! now see if any of the shared dimensions (will be x and y) have the same length
                do p=1,size(t_var_dims)
                    dim_name = t_dimensions(p)
                    dim_indx = findloc(dimensions,dim_name,dim=1)
                    ! check if the dimension is present in the variable
                    if (dim_indx > 0) then
                        ! check if the dimension length matches the terrain height variable
                        dim_len = var_dims(dim_indx)
                        if (name=='ulat' .and. dim_name=='X') dim_len = dim_len - 1
                        if (name=='ulon' .and. dim_name=='X') dim_len = dim_len - 1
                        if (name=='vlat' .and. dim_name=='Y') dim_len = dim_len - 1
                        if (name=='vlon' .and. dim_name=='Y') dim_len = dim_len - 1
                        if (name=='vvar' .and. dim_name=='Y' .and. dim_len == t_var_dims(p)+1) dim_len = dim_len - 1
                        if (name=='uvar' .and. dim_name=='X' .and. dim_len == t_var_dims(p)+1) dim_len = dim_len - 1
                        if (dim_name=='Z' .and. dim_len == t_var_dims(p)+1) then
                            if (STD_OUT_PE) write(*,*) "--------------------------------------------------------------------------------"
                            if (STD_OUT_PE) write(*,*) "ATTENTION: vertical dimension ",trim(dim_name)," on forcing variable ", trim(name)
                            if (STD_OUT_PE) write(*,*) "ATTENTION: has dimension length:",var_dims(dim_indx)
                            if (STD_OUT_PE) write(*,*) "ATTENTION: variable: ",trim(options%tvar)," has length: ", t_var_dims(p), "for the same dimension"
                            if (STD_OUT_PE) write(*,*) "ATTENTION: We are interpreting this as variable ",trim(name), " being on"
                            if (STD_OUT_PE) write(*,*) "ATTENTION: the k-level interfaces, and will automatically linearly interpolate"
                            if (STD_OUT_PE) write(*,*) "ATTENTION: it to the mass-points"
                            if (STD_OUT_PE) write(*,*) "--------------------------------------------------------------------------------"
                            dim_len = dim_len - 1
                        endif

                        if (dim_len /= t_var_dims(p)) then
                            if (STD_OUT_PE) write(*,*) "Error: dimension ",trim(dim_name)," on forcing variable ", trim(name)
                            if (STD_OUT_PE) write(*,*) "has dimension length:",var_dims(dim_indx)
                            if (STD_OUT_PE) write(*,*) "But variable: ",trim(options%tvar)," has length: ", t_var_dims(p), "for the same dimension"
                            ! output special warning if it is a staggered variable
                            if (STD_OUT_PE .and. .not.(dim_len == var_dims(dim_indx))) write(*,*) "should be ",t_var_dims(p)+1," for staggered vars"
                            error stop
                        endif
                    endif
                end do
            endif
        ! Check if the variable was not set. In this case, we just set it to be blank ("")
        elseif (trim(var_val) == kCHAR_NO_VAL) then
            var = ""
            return
        endif

        ! if (var == kCHAR_NO_VAL) then
        !     var = ""
            ! if (STD_OUT_PE) write(*,*) "Error: '", trim(name), "' is not set to any value (i.e. is still ",(kCHAR_NO_VAL),")"
            ! error stop
        ! endif

        ! if i is present, then we want to add this variable to the vars-to-read list
        if (present(i)) then
            options%vars_to_read(i) = var_val
            ! only count spatial dimensions
            options%dim_list(i)%num_dims = count(.not.(dimensions=='T'))

            if (options%dim_list(i)%num_dims /= size(dimensions)) then
                options%dim_list(i)%dims = pack(var_dims, .not.(dimensions=='T'))
            else
                options%dim_list(i)%dims = var_dims
            endif
            i = i + 1
        endif
        var = var_val


    end subroutine set_char_forcing_nml_var

    subroutine set_char_domain_nml_var(var, var_val, name, options, usr_default)
        implicit none
        character(len=*),    intent(inout) :: var
        character(len=*),    intent(in) :: var_val
        character(len=*), intent(in) :: name
        type(domain_options_type),  intent(inout) :: options
        character(len=*), optional, intent(in) :: usr_default
        

        character(len=kMAX_STRING_LENGTH) :: group, default, description, units
        real :: min, max
        integer, allocatable :: values(:), var_dims(:), hgt_var_dims(:)
        integer :: i, dim_indx, dim_len, type
        character(len=kMAX_STRING_LENGTH) :: dim_name
        character(len=1), allocatable :: dimensions(:), hgt_dimensions(:), dimensions_tmp(:), hgt_dimensions_tmp(:)

        call get_nml_var_metadata(name,group,description,default,min,max,type,values,units,dimensions)

        if (present(usr_default)) then
            ! If the user set the first value, but not the rest, then they want this first value to apply to all nests
            if ( (trim(usr_default) /= kCHAR_NO_VAL) .and. (trim(var_val) == kCHAR_NO_VAL )) then
                var = usr_default
            elseif ( trim(usr_default) == kCHAR_NO_VAL ) then
                ! If the user did not set the first value, then there is no entry at all, so set it to the default
                var = default
            else
                var = var_val
            endif
        ! If this is not a variable for which nested values can differ, then we just check if the value was not set in the namelist
        else if (trim(var_val) == kCHAR_NO_VAL) then
            var = default
        else
            var = var_val
        endif

        ! see if user wants to use this variable anyways
        if (.not.(var == "")) then
            ! first check if variable is present in the domain file
            call check_variable_present(options%init_conditions_file, var)
            ! then get the number of dimensions for the variable
            call io_getdims(options%init_conditions_file, var, var_dims)
            ! test if the number of dimensions matches the number of dimensions required
            if (.not.(size(dimensions) == size(var_dims))) then
                if (STD_OUT_PE) write(*,*) "Error: ", trim(name), " has the wrong number of dimensions, should be: ", size(dimensions)
                error stop
            endif

            ! pointless if this is hgt_hi
            if (.not.(name=='hgt_hi')) then
                ! now get the dimension lengths of the terrain height variable
                call check_variable_present(options%init_conditions_file, options%hgt_hi)
                !! then get the number of dimensions for the variable
                call io_getdims(options%init_conditions_file, options%hgt_hi, hgt_var_dims)
        !
                call get_nml_var_metadata('hgt_hi',group,description,default,min,max,type,values,units,hgt_dimensions)

                ! reverse ordering of hgt_dimensions -- needed since fortran netCDF uses reverse ordering relative to python/nco
                allocate(hgt_dimensions_tmp(size(hgt_dimensions)))
                hgt_dimensions_tmp = hgt_dimensions
                do i=1,size(hgt_dimensions)
                    hgt_dimensions(i) = hgt_dimensions_tmp(size(hgt_dimensions_tmp)-i+1)
                end do

                ! reverse ordering of dimensions -- needed since fortran netCDF uses reverse ordering relative to python/nco
                allocate(dimensions_tmp(size(dimensions)))
                dimensions_tmp = dimensions
                do i=1,size(dimensions)
                    dimensions(i) = dimensions_tmp(size(dimensions_tmp)-i+1)
                end do

                !! now see if any of the shared dimensions (will be x and y) have the same length
                do i=1,size(hgt_var_dims)
                    dim_name = hgt_dimensions(i)
                    dim_indx = findloc(dimensions,dim_name,dim=1)
                    ! check if the dimension is present in the variable
                    if (dim_indx > 0) then
                        ! check if the dimension length matches the terrain height variable
                        dim_len = var_dims(dim_indx)
                        if (name=='ulat_hi' .and. dim_name=='X') dim_len = dim_len - 1
                        if (name=='ulon_hi' .and. dim_name=='X') dim_len = dim_len - 1
                        if (name=='vlat_hi' .and. dim_name=='Y') dim_len = dim_len - 1
                        if (name=='vlon_hi' .and. dim_name=='Y') dim_len = dim_len - 1
        !
                        if (dim_len /= hgt_var_dims(i)) then
                            if (STD_OUT_PE) write(*,*) "Error: dimension ",trim(dim_name)," on domain variable ", trim(name)
                            if (STD_OUT_PE) write(*,*) "has dimension length: ",var_dims(dim_indx)
                            if (STD_OUT_PE) write(*,*) "But variable: ",trim(options%hgt_hi)," has length: ", hgt_var_dims(i), "for the same dimension"
                            ! output special warning if it is a staggered variable
                            if (STD_OUT_PE .and. .not.(dim_len == var_dims(dim_indx))) write(*,*) " should be ",hgt_var_dims(i)+1," for staggered vars"
                            error stop
                        endif
                    endif
                end do
            endif
        endif

        if (present(usr_default)) then
            if (usr_default /= kCHAR_NO_VAL) then
                call check_var_entry_type(usr_default,var,name)
            else
                call check_var_entry_type(default,var,name)
            endif
        endif


        ! Check if the variable was not set. This would be an error
        if (var == kCHAR_NO_VAL .and. default /= kCHAR_NO_VAL) then
            if (STD_OUT_PE) write(*,*) "Error: '", trim(name), "' is not set to any value (i.e. is still ",(kCHAR_NO_VAL),")"
            error stop
        endif

    end subroutine set_char_domain_nml_var

    subroutine set_char_nml_var(var, var_val, name, usr_default)
        implicit none
        character(len=*),    intent(inout) :: var
        character(len=*),    intent(in) :: var_val
        character(len=*), intent(in) :: name
        character(len=*), optional, intent(in) :: usr_default

        character(len=kMAX_STRING_LENGTH) :: default

        default = trim(get_nml_var_default(name))

        if (present(usr_default)) then
            ! If the user set the first value, but not the rest, then they want this first value to apply to all nests
            if ( (trim(usr_default) /= trim(kCHAR_NO_VAL)) .and. (trim(var_val) == trim(kCHAR_NO_VAL) )) then
                var = trim(usr_default)
            elseif ( trim(usr_default) == trim(kCHAR_NO_VAL) ) then
                ! If the user did not set the first value, then there is no entry at all, so set it to the default
                var = trim(default)
            else
                var = trim(var_val)
            endif
        ! If this is not a variable for which nested values can differ, then we just check if the value was not set in the namelist
        else if (trim(var_val) == trim(kCHAR_NO_VAL)) then
            var = trim(default)
        else
            var = trim(var_val)
        endif

        ! Perform regex match on var_val to see if it ends with '.txt'
        ! If it does, then it is a file path, and we need to check if it exists
        ! If it does not, then it is a string, and we can set it directly

        if ((index(var, '.txt') > 0) .or. index(var, '.nc') > 0) then
            call check_file_exists(var, message=(trim(var)//"file does not exist."))
        endif

        if (present(usr_default)) then
            if (usr_default /= kCHAR_NO_VAL) then
                call check_var_entry_type(usr_default,var,name)
            else
                call check_var_entry_type(default,var,name)
            endif
        endif

        ! Check if the variable was not set. This would be an error
        if (var == kCHAR_NO_VAL .and. default /= kCHAR_NO_VAL) then
            if (STD_OUT_PE) write(*,*) "Error: '", trim(name), "' is not set to any value (i.e. is still ",(kCHAR_NO_VAL),")"
            error stop
        endif

    end subroutine set_char_nml_var

    subroutine set_real_nml_var_default(var, name, info, gen_nml)
        implicit none
        real,    intent(inout) :: var
        character(len=*), intent(in) :: name
        logical,          intent(in), optional :: info, gen_nml

        character(256) :: default
        logical :: print_info, gennml

        print_info = .False.
        if (present(info)) print_info = info

        gennml = .False.
        if (present(gen_nml)) gennml = gen_nml

        if (gennml .or. info) then
            !convert the default value string to a real
            default = trim(get_nml_var_default(name,print_info,gennml))
            read(default,*) var
        else
            var = kREAL_NO_VAL
        endif


    end subroutine set_real_nml_var_default

    subroutine set_integer_nml_var_default(var, name, info, gen_nml)
        implicit none
        integer,    intent(inout) :: var
        character(len=*), intent(in) :: name
        logical,          intent(in), optional :: info, gen_nml

        character(256) :: default
        logical :: print_info, gennml

        print_info = .False.
        if (present(info)) print_info = info

        gennml = .False.
        if (present(gen_nml)) gennml = gen_nml

        if (gennml .or. info) then
            !convert the default value string to an integer
            default = trim(get_nml_var_default(name,print_info,gennml))
            read(default,*) var
        else
            var = kINT_NO_VAL
        endif

    end subroutine set_integer_nml_var_default

    subroutine set_integer_array_nml_var_default(var, name, info, gen_nml)
        implicit none
        integer,    intent(inout) :: var(:)
        character(len=*), intent(in) :: name
        logical,          intent(in), optional :: info, gen_nml

        character(256) :: default
        integer :: scalar_default
        logical :: print_info, gennml

        print_info = .False.
        if (present(info)) print_info = info

        gennml = .False.
        if (present(gen_nml)) gennml = gen_nml

        if (gennml .or. info) then
            !convert the default value string to an integer
            default = trim(get_nml_var_default(name,print_info,gennml))
            read(default,*) scalar_default
            var = scalar_default
        else
            var = kINT_NO_VAL
        endif

    end subroutine set_integer_array_nml_var_default

    subroutine set_real_list_nml_var_default(var, name, info, gen_nml)
        implicit none
        real,  intent(inout) :: var(:)
        character(len=*), intent(in) :: name
        logical,          intent(in), optional :: info, gen_nml

        character(256) :: default
        real :: scalar_default
        logical :: print_info, gennml

        print_info = .False.
        if (present(info)) print_info = info

        gennml = .False.
        if (present(gen_nml)) gennml = gen_nml

        if (gennml .or. info) then
            !convert the default value string to a real
            default = trim(get_nml_var_default(name,print_info,gennml))
            read(default,*) scalar_default
            var = scalar_default
        else
            var = kREAL_NO_VAL
        endif


    end subroutine set_real_list_nml_var_default

    subroutine set_real_list_array_nml_var_default(var, name, info, gen_nml)
        implicit none
        real,  intent(inout) :: var(:,:)
        character(len=*), intent(in) :: name
        logical,          intent(in), optional :: info, gen_nml

        character(256) :: default
        real :: scalar_default
        logical :: print_info, gennml

        print_info = .False.
        if (present(info)) print_info = info

        gennml = .False.
        if (present(gen_nml)) gennml = gen_nml

        if (gennml .or. info) then
            !convert the default value string to a real
            default = trim(get_nml_var_default(name,print_info,gennml))
            read(default,*) scalar_default
            var = scalar_default
        else
            var = kREAL_NO_VAL
        endif


    end subroutine set_real_list_array_nml_var_default

    subroutine set_logical_nml_var_default(var, name, info, gen_nml)
        implicit none
        logical,    intent(inout) :: var
        character(len=*), intent(in) :: name
        logical,          intent(in), optional :: info, gen_nml

        logical :: print_info, gennml

        print_info = .False.
        if (present(info)) print_info = info

        gennml = .False.
        if (present(gen_nml)) gennml = gen_nml

        if (get_nml_var_default(name,print_info,gennml) == ".True.") then
            var = .True.
        else
            var = .False.
        endif

    end subroutine set_logical_nml_var_default

    subroutine set_logical_array_nml_var_default(var, name, info, gen_nml)
        implicit none
        logical,    intent(inout) :: var(:)
        character(len=*), intent(in) :: name
        logical,          intent(in), optional :: info, gen_nml

        logical :: print_info, gennml

        print_info = .False.
        if (present(info)) print_info = info

        gennml = .False.
        if (present(gen_nml)) gennml = gen_nml

        if (get_nml_var_default(name,print_info,gennml) == ".True.") then
            var = .True.
        else
            var = .False.
        endif

    end subroutine set_logical_array_nml_var_default

    subroutine set_char_nml_var_default(var, name, info, gen_nml)
        implicit none
        character(len=*), intent(inout) :: var
        character(len=*), intent(in) :: name
        logical,          intent(in), optional :: info, gen_nml

        logical :: print_info, gennml

        print_info = .False.
        if (present(info)) print_info = info

        gennml = .False.
        if (present(gen_nml)) gennml = gen_nml

        if (gennml .or. info) then
            var = get_nml_var_default(name,print_info,gennml)
        else
            var = kCHAR_NO_VAL
        endif

    end subroutine set_char_nml_var_default

    subroutine set_char_array_nml_var_default(var, name, info, gen_nml)
        implicit none
        character(len=*), intent(inout) :: var(:)
        character(len=*), intent(in) :: name
        logical,          intent(in), optional :: info, gen_nml

        logical :: print_info, gennml

        print_info = .False.
        if (present(info)) print_info = info

        gennml = .False.
        if (present(gen_nml)) gennml = gen_nml

        if (gennml .or. info) then
            var = get_nml_var_default(name,print_info,gennml)
        else
            var = kCHAR_NO_VAL
        endif

    end subroutine set_char_array_nml_var_default

    subroutine set_char_list_array_nml_var_default(var, name, info, gen_nml)
        implicit none
        character(len=*), intent(inout) :: var(:,:)
        character(len=*), intent(in) :: name
        logical,          intent(in), optional :: info, gen_nml

        logical :: print_info, gennml
        integer :: n

        print_info = .False.
        if (present(info)) print_info = info

        gennml = .False.
        if (present(gen_nml)) gennml = gen_nml

        if (gennml .or. info) then
            ! do n = 1, size(var,2)
            var(:,1) = get_nml_var_default(name,print_info,gennml)
            ! end do
        else
            var = kCHAR_NO_VAL
        endif

    end subroutine set_char_list_array_nml_var_default

    subroutine check_var_entry_type(val_1,val_n,name)
        implicit none
        character(len=*), intent(in) :: val_1, val_n
        character(len=*), intent(in) :: name

        character(len=kMAX_STRING_LENGTH) :: default
        integer :: type    

        !convert the default value string to a real
        default = trim(get_nml_var_default(name))
        type = get_nml_var_type(name)

        if (type==1) then
            if (trim(val_1) /= trim(val_n)) then
                if (STD_OUT_PE) write(*,*) "------------------------------------------------------------------------------------"
                if (STD_OUT_PE) write(*,*) "Error: namelist option '", trim(name), "' must be the same for all nests" 
                if (STD_OUT_PE) write(*,*) "Error: the first entry has the value: ", trim(val_1), " and another entry has the value: ", trim(val_n)
                if (STD_OUT_PE) write(*,*) "------------------------------------------------------------------------------------"
                stop
            endif
        elseif (type==2) then
            if (trim(val_n) == default) then
                if (STD_OUT_PE) write(*,*) "------------------------------------------------------------------------------------"
                if (STD_OUT_PE) write(*,*) "Error: namelist option '", trim(name), "' is not uniquely set for each of the nests" 
                if (STD_OUT_PE) write(*,*) "------------------------------------------------------------------------------------"
                stop
            endif
        endif

    end subroutine check_var_entry_type

    subroutine check_numeric_var_entry_type(val_1,val_n,name)
        implicit none
        real, intent(in) :: val_1, val_n
        character(len=*), intent(in) :: name

        character(len=kMAX_STRING_LENGTH) :: default
        real :: scalar_default
        integer :: type    

        !convert the default value string to a real
        default = trim(get_nml_var_default(name))
        read(default,*) scalar_default
        type = get_nml_var_type(name)

        if (type==1) then
            if (val_1 /= val_n) then
                if (STD_OUT_PE) write(*,*) "------------------------------------------------------------------------------------"
                if (STD_OUT_PE) write(*,*) "Error: namelist option '", trim(name), "' must be the same for all nests" 
                if (STD_OUT_PE) write(*,*) "Error: the first entry has the value: ", val_1, " and another entry has the value: ", val_n
                if (STD_OUT_PE) write(*,*) "------------------------------------------------------------------------------------"
                stop
            endif
        elseif (type==2) then
            if (val_n == scalar_default) then
                if (STD_OUT_PE) write(*,*) "------------------------------------------------------------------------------------"
                if (STD_OUT_PE) write(*,*) "Error: namelist option '", trim(name), "' is not uniquely set for each of the nests" 
                if (STD_OUT_PE) write(*,*) "------------------------------------------------------------------------------------"
                stop
            endif
        endif

    end subroutine check_numeric_var_entry_type


    function get_nml_var_minmax(name) result(minmax)
        implicit none
        character(len=*), intent(in) :: name

        character(len=kMAX_STRING_LENGTH) :: group, default, description, units
        character(len=1), allocatable :: dimensions(:)
        real :: min, max
        integer :: type
        integer, allocatable :: values(:)
        real :: minmax(2)

        call get_nml_var_metadata(name,group,description,default,min,max,type,values,units,dimensions)

        minmax = (/min, max/)

    end function get_nml_var_minmax

    function get_nml_var_values(name) result(values)
        implicit none
        character(len=*), intent(in) :: name

        character(len=kMAX_STRING_LENGTH) :: group, default, description, units
        character(len=1), allocatable :: dimensions(:)
        real :: min, max
        integer :: type
        integer, allocatable :: values(:)

        call get_nml_var_metadata(name,group,description,default,min,max,type,values,units,dimensions)

    end function get_nml_var_values

    function get_nml_var_default(name,info,gen_nml) result(default)
        implicit none
        character(len=*), intent(in) :: name
        logical,          intent(in), optional :: info, gen_nml

        character(len=kMAX_STRING_LENGTH) :: group, default, description, units
        character(len=1), allocatable :: dimensions(:)
        real :: min, max
        integer :: type
        integer, allocatable :: values(:)

        call get_nml_var_metadata(name,group,description,default,min,max,type,values,units,dimensions)

        if (present(info)) then
            if (info) call write_nml_var_info(name,group,description,default,min,max,values,units,dimensions)
        endif
        if (present(gen_nml)) then
            if (gen_nml) call write_nml_entry(name, default, description, dimensions, group)
        endif
        
    end function get_nml_var_default

    function get_nml_var_type(name) result(type)
        implicit none
        character(len=*), intent(in) :: name

        character(len=kMAX_STRING_LENGTH) :: group, default, description, units
        character(len=1), allocatable :: dimensions(:)
        real :: min, max
        integer :: type
        integer, allocatable :: values(:)

        call get_nml_var_metadata(name,group,description,default,min,max,type,values,units,dimensions)
        
    end function get_nml_var_type

    subroutine write_nml_entry(name, default, description, dimensions, group)
        implicit none
        character(len=*), intent(in) :: name, default, group
        character(len=*), intent(inout) :: description
        character(len=1), allocatable :: dimensions(:)

        character(len=30), save :: last_group = ""
        character(256) :: var_string
        integer :: indx, indx_old, print_len, num_blnks
        logical :: numeric_only, is_logical

        if (.not.(STD_OUT_PE)) return

        if (last_group /= group) then
            ! If this is not the first group, close the last group
            if (.not.(last_group == "")) then
                write(nml_unit,*) "/"
                write(nml_unit,*)
            endif
            call write_group_header(group)
        endif

        write(nml_unit,*)
        indx = index(name, "file")

        numeric_only = (0 == VERIFY(trim(default), "0123456789.+-eE"))
        is_logical = (trim(default) == ".True." .or. trim(default) == ".False.")
        ! If the name contains file, has dimensions, or has a default value of an empty string, then put quotes around the default value
        if ((.not.numeric_only .and. .not.is_logical) .or. trim(default)=="") then
            var_string = "    "//trim(name)//" = '"//trim(default)//"'"
        else
            var_string = "    "//trim(name)//" = "//trim(default)
        endif

        ! write the variable name and default value, not advancing to next line
        write(nml_unit,'(A)',advance="no") trim(var_string)

        ! Three cases: description is formatted with newlines, description is too long for one line, or description is short enough for one line

        ! calculate the number of blank characters to write for the first comment 
        num_blnks = NML_BLNK_LEN - len(trim(var_string))

        ! see if is formatted with newlines
        indx = index(trim(description),achar(10))

        ! If description is multi-line
        if (indx > 0) then
            do while (len(trim(description)) > 0)
                write(nml_unit, '('//str(num_blnks)//'X,A3,A)' ,advance="no") " ! ", trim(description(1:indx))

                ! if formatted, description will have some number of blank spaces until the next line. 
                ! scan forward to the next non-blank character
                indx = indx + 1
                do while (indx <= len(trim(description)) .and. trim(description(indx:indx)) == " ")
                    indx = indx + 1
                end do
                
                description = description(indx:)
                indx = index(trim(description),achar(10))
                ! if there is no newline character, set indx to the length of the description to print the rest
                if (indx == 0) then
                    indx = len(trim(description))
                endif
                num_blnks = NML_BLNK_LEN
            end do
            write(nml_unit,*)
        else
            ! else see if the description is too long for one line

            ! need +3 to account for the ' ! ' at the beginning of the comment
            print_len = min(SCREEN_WIDTH-len(trim(var_string))-num_blnks-3, len(trim(description)))

            ! If description is multi-line
            if (print_len < len(trim(description))) then
                do while (len(trim(description)) > 0)
                    
                    !if this is not the last line, find the last space to break the line
                    if (.not.print_len == len(trim(description))) then
                        print_len = index(description(1:print_len), " ", BACK=.TRUE.)
                    endif

                    write(nml_unit, '('//str(num_blnks)//'X,A3,A)') " ! ", trim(description(1:print_len))
                    description = description(print_len+1:)
                    num_blnks = NML_BLNK_LEN
                    print_len = min(SCREEN_WIDTH-num_blnks-3, len(trim(description)))
                end do
            else    
            ! If description is single line
                write(nml_unit, '('//str(num_blnks)//'X,A3,A)' ) " ! ", trim(description)
            endif
        endif

        last_group = group

    end subroutine write_nml_entry

    subroutine write_nml_file_end()

        if (.not.(STD_OUT_PE)) return

        ! Open the namelist file
        write(nml_unit,*) "/"
        close(nml_unit)

    end subroutine write_nml_file_end

    subroutine write_group_header(group)
        implicit none
        character(len=*), intent(in) :: group

        select case (group)
            case ("General")
                write(nml_unit,*) "! ---------------------------------------------------------------------------------------------------"
                write(nml_unit,*) "!   HICAR model namelist file"
                write(nml_unit,*) "!                                                                                                    "
                write(nml_unit,*) "!   When setting namelist variables, the user has 3 options:"
                write(nml_unit,*) "!     1. Leave the variable blank, and the default value, which is shown in this default namelist,"
                write(nml_unit,*) "!        will be used for all domains"
                write(nml_unit,*) "!     2. Set the variable to a single value, and that value will be used for all domains."
                write(nml_unit,*) "!     3. Set the variable to a list of values, where the index of each entry in the list corresponds"
                write(nml_unit,*) "!        to the index of the domain which the value should be used for."
                write(nml_unit,*) "!                                                                                                    "
                write(nml_unit,*) "!   If the variable itself is a list, for example dz_levels, and the user wishes to set different"
                write(nml_unit,*) "!   values for different domains, then the following syntax is used:"
                write(nml_unit,*) "!     dz_levels(:,1) = 0.0, 1000.0, 2000.0, 3000.0, 4000.0, 5000.0, 6000.0, 7000.0, 8000.0, 9000.0"
                write(nml_unit,*) "!     dz_levels(:,2) = 0.0, 500.0,  1500.0, 3000.0, 4000.0, 5000.0, 6000.0, 7000.0, 8000.0, 9000.0"
                write(nml_unit,*) "!                                                                                                    "
                write(nml_unit,*) "!    Most all variables can be set using any of the 3 options listed above, with two exceptions."
                write(nml_unit,*) "!    When a variable is marked with:"
                write(nml_unit,*) "!      (-) the variable can only be set with a single value, and is thus the same for all domains."
                write(nml_unit,*) "!      (*) the variable must be explicitly set for each domain."
                write(nml_unit,*) "! ---------------------------------------------------------------------------------------------------"
                write(nml_unit,*) "                                                                                                     "
                write(nml_unit,*) "                                                                                                     "
                write(nml_unit,*) "! ---------------------------------------------------------------------------------------------------"
                write(nml_unit,*) "!   General model and run meta-data"
                write(nml_unit,*) "! ---------------------------------------------------------------------------------------------------"
                write(nml_unit,*) "&general"
            case ("Domain")
                write(nml_unit,*) "! ---------------------------------------------------------------------------------------------------"
                write(nml_unit,*) "!   Information about the modeling domain, including domain file name,"
                write(nml_unit,*) "!   vertical grid structure, and variables to read from domain file"
                write(nml_unit,*) "! ---------------------------------------------------------------------------------------------------"
                write(nml_unit,*) "&domain"
            case ("Forcing")
                write(nml_unit,*) "! ---------------------------------------------------------------------------------------------------"
                write(nml_unit,*) "!   Information about the forcing data, including file list,"
                write(nml_unit,*) "!   conventions for the data, and names of variables to read"
                write(nml_unit,*) "!   Only used for domains whose parent_nest is 0"
                write(nml_unit,*) "! ---------------------------------------------------------------------------------------------------"
                write(nml_unit,*) "&forcing"
            case ("Output")
                write(nml_unit,*) "! ---------------------------------------------------------------------------------------------------"
                write(nml_unit,*) "!    Specify output files and variables"
                write(nml_unit,*) "! ---------------------------------------------------------------------------------------------------"
                write(nml_unit,*) "&output"
            case ("Restart")
                write(nml_unit,*) "! ---------------------------------------------------------------------------------------------------"
                write(nml_unit,*) "!   Specify if a restart run, and restart files"
                write(nml_unit,*) "! ---------------------------------------------------------------------------------------------------"
                write(nml_unit,*) "&restart"
            case ("Physics")
                write(nml_unit,*) "! ---------------------------------------------------------------------------------------------------"
                write(nml_unit,*) "!   Specify physics options to use for the model run"
                write(nml_unit,*) "! ---------------------------------------------------------------------------------------------------"
                write(nml_unit,*) "&physics"
            case ("MP_Parameters")
                write(nml_unit,*) "! ---------------------------------------------------------------------------------------------------"
                write(nml_unit,*) "!   Optionally specified Microphysics parameters (mostly for Thompson)"
                write(nml_unit,*) "! ---------------------------------------------------------------------------------------------------"
                write(nml_unit,*) "&mp_parameters"
            case ("LT_Parameters")
                write(nml_unit,*) "! ---------------------------------------------------------------------------------------------------"
                write(nml_unit,*) "!   Optionally specified Linear Theory parameters"
                write(nml_unit,*) "! ---------------------------------------------------------------------------------------------------"
                write(nml_unit,*) "&lt_parameters"
            case ("SFC_Parameters")
                write(nml_unit,*) "! ---------------------------------------------------------------------------------------------------"
                write(nml_unit,*) "!   Optionally specified surface layer scheme parameters"
                write(nml_unit,*) "! ---------------------------------------------------------------------------------------------------"
                write(nml_unit,*) "&sfc_parameters"
            case ("RAD_Parameters")
                write(nml_unit,*) "! ---------------------------------------------------------------------------------------------------"
                write(nml_unit,*) "!   Optionally specified radiation parameters"
                write(nml_unit,*) "! ---------------------------------------------------------------------------------------------------"
                write(nml_unit,*) "&rad_parameters"
            case ("PBL_Parameters")
                write(nml_unit,*) "! ---------------------------------------------------------------------------------------------------"
                write(nml_unit,*) "!   Optionally specified PBL parameters"
                write(nml_unit,*) "! ---------------------------------------------------------------------------------------------------"
                write(nml_unit,*) "&pbl_parameters"
            case ("LSM_Parameters")
                write(nml_unit,*) "! ---------------------------------------------------------------------------------------------------"
                write(nml_unit,*) "!   Optionally specified land surface model parameters (mostly for NoahMP)"
                write(nml_unit,*) "! ---------------------------------------------------------------------------------------------------"
                write(nml_unit,*) "&lsm_parameters"
            case ("SM_Parameters")
                write(nml_unit,*) "! ---------------------------------------------------------------------------------------------------"
                write(nml_unit,*) "!   Optionally specified snow model parameters (mostly for FSM2trans, at the moment)"
                write(nml_unit,*) "! ---------------------------------------------------------------------------------------------------"
                write(nml_unit,*) "&sm_parameters"
            case ("CU_Parameters")
                write(nml_unit,*) "! ---------------------------------------------------------------------------------------------------"
                write(nml_unit,*) "!   Optionally specified convection parameters"
                write(nml_unit,*) "! ---------------------------------------------------------------------------------------------------"
                write(nml_unit,*) "&cu_parameters"
            case ("Wind")
                write(nml_unit,*) "! ---------------------------------------------------------------------------------------------------"
                write(nml_unit,*) "!   Specify wind solver and options for downscaling"
                write(nml_unit,*) "! ---------------------------------------------------------------------------------------------------"
                write(nml_unit,*) "&wind"
            case ("ADV_Parameters")
                write(nml_unit,*) "! ---------------------------------------------------------------------------------------------------"
                write(nml_unit,*) "!   Optionally specified advection parameters"
                write(nml_unit,*) "! ---------------------------------------------------------------------------------------------------"
                write(nml_unit,*) "&adv_parameters"
            case ("Time_Parameters")
                write(nml_unit,*) "! ---------------------------------------------------------------------------------------------------"
                write(nml_unit,*) "!   Specify options for time-stepping"
                write(nml_unit,*) "! ---------------------------------------------------------------------------------------------------"
                write(nml_unit,*) "&time_parameters"
        end select

    end subroutine write_group_header

    subroutine write_nml_var_info(name, group, description, default, min, max, values, units, dimensions)
        implicit none
        character(len=*), intent(in) :: name
        character(len=*), intent(in) :: group, description, default, units
        character(len=*), intent(in) :: dimensions(:)
        real,    intent(in) :: min, max
        integer, allocatable, intent(in) :: values(:)

        logical :: is_minmax, is_values, is_units, is_dimensions
        integer :: i

        is_values = (.not.(size(values)==0)) 

        is_minmax = (.not.(min==0 .and. max==0))

        is_units = (.not.(units==""))

        is_dimensions = (.not.(size(dimensions)==0)) 

        if (STD_OUT_PE) write(*,*) "---------------------------------------------------------------------------------------------------"
        if (STD_OUT_PE) write(*,*)                        "Namelist Variable:           ", trim(name)
        if (STD_OUT_PE) write(*,*)                        "    Namelist Group:.........|", trim(group)
        if (STD_OUT_PE) write(*,*)                        "    Description:............|", trim(description)
        if (STD_OUT_PE .and. .not.default=="") write(*,*) "    Default Value:..........|", trim(default)
        if (STD_OUT_PE .and. is_minmax) write(*,*)        "    Minimum Allowed Value:..|", str(min)
        if (STD_OUT_PE .and. is_minmax) write(*,*)        "    Maximum Allowed Value:..|", str(max)
        if (STD_OUT_PE .and. is_units) write(*,*)         "    Units:..................|", trim(units)
        if (STD_OUT_PE .and. is_dimensions) then
            write(*,'(A$)', advance="no")      "     Dimensions:.............|["
            do i = 1, size(dimensions)
                if (i==size(dimensions)) then
                    WRITE(*, '(A1, A)', ADVANCE='NO') dimensions(i), "]"
                else
                    WRITE(*, '(A1, A)', ADVANCE='NO') dimensions(i), ', '
                endif
            end do
            write(*,*)
        endif
        if (STD_OUT_PE .and. is_values) then
            write(*,'(A$)', advance="no")                   "    Allowed Values:..........|"
            do i = 1, size(values)
                if (i==size(values)) then
                    WRITE(*, '(I1)', ADVANCE='NO') values(i)
                else
                    WRITE(*, '(I1, A)', ADVANCE='NO') values(i), ', '
                endif
            end do
            write(*,*)
        endif

    end subroutine write_nml_var_info


    subroutine get_nml_var_metadata(name, group, description, default, min, max, type, values, units, dimensions)
        implicit none
        character(len=*), intent(in) :: name
        character(len=*), intent(out) :: group, description, default, units
        character(len=*), allocatable, intent(out) :: dimensions(:)
        real,    intent(out) :: min, max
        integer, intent(out) :: type
        integer, allocatable, intent(out) :: values(:)

        group = ""
        description = ""
        default = ""
        units = ""
        type = 0

        min = 0
        max = 0
        select case (name)
            ! --------------------------------------
            ! --------------------------------------
            ! General namelist variables
            ! --------------------------------------
            ! --------------------------------------
            case ("debug")
                description = "Debugging flag (T/F)"
                default = ".False."
                group = "General"
                type = 1
            case ("interactive")
                description = "Interactive flag, prints out model progress during physics timesteps (T/F)"
                default = ".True."
                group = "General"
                type = 1
            case ("calendar")
                description = "Calendar type for start/end date, values: (gregorian, noleap, 360_day)"
                default = "gregorian"
                group = "General" 
                type = 1
            case ("start_date")
                description = "Start date for simulation, format: 'YYYY-MM-DD HH:MM:SS'"
                default = ""
                group = "General"
            case ("end_date")
                description = "End date for simulation, format: 'YYYY-MM-DD HH:MM:SS'"
                default = ""
                group = "General"
            case ("version")
                description = "Version of the model"
                default = kVERSION_STRING
                group = "General"
                type = 1
            case ("comment")
                description = "Comment for the run, to be added to the output file metadata"
                group = "General"
                type = 1
            case ("nests")
                description = "Number of nests to use in the simulation"
                default = "1"
                min = 1
                max = kMAX_NESTS
                group = "General"
                type = 1
            case ("parent_nest")
                description = "Parent nest indices for each nest, length 'nests'. "
                default = "0"
                min = 0
                max = kMAX_NESTS
                group = "General"
            case ("phys_suite")
                description = "Physics suite to use, current options: (HICAR)"
                group = "General"
                type = 1
            case ("use_mp_options")
                description = "Read the microphysics namelist section to set options relevant for the Thompson MP schemes (T/F)"
                default = ".True."
                group = "General"
                type = 1
            case ("use_lt_options")
                description = "Read the linear theory namelist section to set options relevant for the linear-theory wind solver (T/F)"
                default = ".True."
                group = "General"
                type = 1
            case ("use_sfc_options")
                description = "Read the surface namelist section to set options relevant for the surface model (T/F)"
                default = ".True."
                group = "General"
                type = 1
            case ("use_rad_options")
                description = "Read the radiation namelist section to set options relevant for the radiation model (T/F)"
                default = ".True."
                group = "General"
                type = 1
            case ("use_pbl_options")
                description = "Read the PBL namelist section to set options relevant for the PBL model (T/F)"
                default = ".True."
                group = "General"
                type = 1
            case ("use_lsm_options")
                description = "Read the LSM namelist section to set options relevant for the LSM model (T/F)"
                default = ".True."
                group = "General"
                type = 1
            case ("use_sm_options")
                description = "Read the SM namelist section to set options relevant for the snow model (T/F)"
                default = ".True."
                group = "General"
                type = 1
            case ("use_cu_options")
                description = "Read the cumulus namelist section to set options relevant for the cumulus scheme (T/F)"
                default = ".True."
                group = "General"
                type = 1
            case ("use_wind_options")
                description = "Read the wind namelist section to set options relevant for the wind solver (T/F)"
                default = ".True."
                group = "General"
                type = 1
            case ("use_adv_options")
                description = "Read the advection namelist section to set options relevant for the advection scheme (T/F)"
                default = ".True."
                group = "General"
                type = 1
            ! --------------------------------------
            ! --------------------------------------
            ! Domain namelist variables
            ! --------------------------------------
            ! --------------------------------------
            case("init_conditions_file")
                description = "Path to file containing initial conditions"
                default = ""
                group = "Domain"
                type = 2
            case("dx")
                description = "Horizontal grid spacing"
                min = 0
                max = 1e6
                units = "meters"
                default = "0.0"
                group = "Domain"
                type = 2
            case("nz")
                description = "Number of vertical levels"
                min = 0
                max = MAXLEVELS
                default = str(MAXLEVELS)
                group = "Domain"
            case("dz_levels")
                description = "Array of vertical grid spacing in meters, length 'nz'. Default of 0.0 is a placeholder to be replaced by the user."
                min = 0
                max = 1e6
                default = "0.0"
                group = "Domain"
            case("longitude_system")
                description = "Longitude system, values: (0=Maintain Longitude, 1=Prime Centered, 2=Dateline Centered, 3=Guess Lon)"
                allocate(values(4))
                values = (/0, 1, 2, 3/)
                default = "0"
                group = "Domain"
            case ("flat_z_height")
                description = "Height at which terrain following coordinates are flat."//achar(10)//BLNK_CHR_N// &
                    "Either height in meters above ground,"//achar(10)//BLNK_CHR_N// &
                    "positive integer representing the index of the height level,"//achar(10)//BLNK_CHR_N// &
                    "or negative integer representing number of levels below model top."
                min = -MAXLEVELS
                max = 1.e6
                default = "-1."
                group = "Domain"
            case ("sleve")
                description = "flag for using SLEVE vertical coordinates (Schär et al. 2002) (T/F)"
                default = ".False."
                group = "Domain"
            case ("terrain_smooth_windowsize")
                description = "Size of the smoothing window used to split large and small scale terrain variations"//achar(10)//BLNK_CHR_N// &
                    "for calculation of the SLEVE coordinate."
                min = 0
                max = 100
                default = "3"
                group = "Domain"
            case ("terrain_smooth_cycles")
                description = "Number of smoothing cycles used to split large and small scale terrain variations"//achar(10)//BLNK_CHR_N// &
                    "for calculation of the SLEVE coordinate."
                min = 0
                max = 200
                default = "5"
                group = "Domain"
            case ("decay_rate_L_topo")
                description = "Decay rate for the large scale topography in the SLEVE coordinate"
                min = 1
                max = 10
                default = "2"
                group = "Domain"
            case ("decay_rate_S_topo")
                description = "Decay rate for the small scale topography in the SLEVE coordinate"
                min = 1
                max = 10
                default = "6"
                group = "Domain"
            case ("sleve_n")
                description = "Exponent 'n' in SLEVE coordinate equation (see Schär et al. 2002)"
                min = 0
                max = 10
                default = "1.2"
                group = "Domain"
            case("auto_sleve")
                description = "Integer that determines whether to create levels automatically when using sleve=.True. (As used in ICON & WRF):" //achar(10)//BLNK_CHR_N// &
                    "Values: 0=No, 1=Third-order polynomial level distribution (ICON like), 2=Second-order polynomial level" //achar(10)//BLNK_CHR_N// &
                    "        distribution (COSMO like) 3=Eta style exponential level distribution (WRF like), 4=Arccosine level" //achar(10)//BLNK_CHR_N// &
                    "        distribution (COSMO like). auto_sleve=3 might be most robust solution, while options 1 and 4 allow" //achar(10)//BLNK_CHR_N// &
                    "        for forced lowest level height. Plotting the distributions beforehand is recommended," //achar(10)//BLNK_CHR_N// & 
                    "        geogebra plots can be found at https://www.geogebra.org/u/maxsesselmann."
                min = 0
                max = 4
                default = "0"
                group = "Domain"
            case("height_lowest_level")
                description = "Lowest level height in meters, can only be forced when using auto_sleve=1 or 4."
                min = 0.1
                max = 1e6
                default = "20.0"
                group = "Domain"
            case("model_top_height")
                description = "Model top height in meters, only used when auto_sleve = 1, 2, 3 or 4."
                min = 0.1
                max = 1e6
                default = "10000.0"
                group = "Domain"
            case("stretch_fac")
                description = "Factor that controls distribution of the vertical levels, only used when auto_sleve = 1, 2, 3 or 4." //achar(10)//BLNK_CHR_N// &
                    "For auto_sleve=1: stretch_fac needs to be between 0.5 and 1.0. stretch_fac -> 0.5 more linear," //achar(10)//BLNK_CHR_N// &
                    " stretch_fac -> 1.0 more level compression at the surface." //achar(10)//BLNK_CHR_N// &
                    "For auto_sleve=2: stretch_fac needs to be between 0.0 and 1.0. stretch_fac -> 0.0 more linear," //achar(10)//BLNK_CHR_N// &
                    " stretch_fac -> 1.0 more level compression at the surface." //achar(10)//BLNK_CHR_N// &
                    "For auto_sleve=3: stretch_fac needs to be > 0. stretch_fac -> 0.0 more linear," //achar(10)//BLNK_CHR_N// &
                    " stretch_fac -> higher values: more level compression at the surface." //achar(10)//BLNK_CHR_N// &
                    "For auto_sleve=4: stretch_fac needs to be > 0. stretch_fac -> 0.0 more compression at the surface," //achar(10)//BLNK_CHR_N// &
                    " stretch_fac -> higher values: more level compression at the model top." //achar(10)//BLNK_CHR_N// &
                    "auto_sleve=3 might be most robust solution, while options 1 and 4 allow for forced lowest level height." //achar(10)//BLNK_CHR_N// &
                    "Plotting the distributions is recommended, " //achar(10)//BLNK_CHR_N// &
                    " geogebra plots can be found at https://www.geogebra.org/u/maxsesselmann."
                min = 0.0001
                max = 10.0
                default = "0.7"
                group = "Domain"
            case ("use_agl_height")
                description = "Use height above ground level, instead of above sea level, "//achar(10)//BLNK_CHR_N// &
                    "for interpolating forcing data to domain grid in lower atmosphere (T/F)"
                default = ".True."
                group = "Domain"
            case ("agl_cap")
                description = "Height above ground level at which interpolation of forcing data switches "//achar(10)//BLNK_CHR_N// &
                    "to using height above sea level. Used if 'use_agl_height' is true."
                min = 0
                max = 1e6
                default = "800.0"
                group = "Domain"
            case ("hgt_hi")
                description = "Name of the high resolution terrain variable in domain file (REQUIRED)"
                units = "meters"
                allocate(dimensions(2))
                dimensions = ["Y", "X"]
                group = "Domain"
            case ("landvar")
                description = "Name of the land mask variable in domain file"
                allocate(dimensions(2))
                dimensions = ["Y", "X"]
                group = "Domain"
            case ("lakedepthvar")
                description = "Name of the lake depth variable in domain file"
                units = "meters"
                allocate(dimensions(2))
                dimensions = ["Y", "X"]
                group = "Domain"
            case ("lat_hi")
                description = "Name of the latitude variable in domain file (REQUIRED)"
                units = "degrees"
                allocate(dimensions(2))
                dimensions = ["Y", "X"]
                group = "Domain"
            case ("lon_hi")
                description = "Name of the longitude variable in domain file (REQUIRED)"
                units = "degrees"
                allocate(dimensions(2))
                dimensions = ["Y", "X"]
                group = "Domain"
            case ("ulat_hi")
                description = "Name of the latitude variable on the staggered U-grid in domain file"
                units = "degrees"
                allocate(dimensions(2))
                dimensions = ["Y", "X"]
                group = "Domain"
            case ("ulon_hi")
                description = "Name of the longitude variable on the staggered U-grid in domain file"
                units = "degrees"
                allocate(dimensions(2))
                dimensions = ["Y", "X"]
                group = "Domain"
            case ("vlat_hi")
                description = "Name of the latitude variable on the staggered V-grid in domain file"
                units = "degrees"
                allocate(dimensions(2))
                dimensions = ["Y", "X"]
                group = "Domain"
            case ("vlon_hi")
                description = "Name of the longitude variable on the staggered V-grid in domain file"
                units = "degrees"
                allocate(dimensions(2))
                dimensions = ["Y", "X"]
                group = "Domain"
            case ("soiltype_var")
                description = "Name of the soil type variable in domain file"
                allocate(dimensions(2))
                dimensions = ["Y", "X"]
                group = "Domain"
            case ("cropcategory_var")
                description = "Name of the crop category variable in domain file"
                allocate(dimensions(2))
                dimensions = ["Y", "X"]
                group = "Domain"
            case ("soil_t_var")
                description = "Name of the soil temperature variable in domain file"
                units = "K"
                allocate(dimensions(2))
                dimensions = ["Y", "X"]
                group = "Domain"
            case ("soil_vwc_var")
                description = "Name of the soil volumetric water content variable in domain file"
                allocate(dimensions(2))
                dimensions = ["Y", "X"]
                group = "Domain"
            case ("swe_var")
                description = "Name of the snow water equivalent variable in domain file"
                units = "mm"
                allocate(dimensions(2))
                dimensions = ["Y", "X"]
                group = "Domain"
            case ("snowh_var")
                description = "Name of the snow depth variable in domain file"
                units = "m"
                allocate(dimensions(2))
                dimensions = ["Y", "X"]
                group = "Domain"
            case ("soil_deept_var")
                description = "Name of the deep soil temperature variable in domain file"
                units = "K"
                allocate(dimensions(2))
                dimensions = ["Y", "X"]
                group = "Domain"
            case ("vegtype_var")
                description = "Name of the vegetation type variable in domain file"
                allocate(dimensions(2))
                dimensions = ["Y", "X"]
                group = "Domain"
            case ("vegfrac_var")
                description = "Name of the vegetation fraction variable in domain file"
                allocate(dimensions(2))
                dimensions = ["Y", "X"]
                group = "Domain"
            case ("vegfracmax_var")
                description = "Name of the maximum vegetation fraction variable in domain file"
                allocate(dimensions(2))
                dimensions = ["Y", "X"]
                group = "Domain"
            case ("albedo_var")
                description = "Name of the albedo variable in domain file"
                allocate(dimensions(2))
                dimensions = ["Y", "X"]
                group = "Domain"
            case ("lai_var")
                description = "Name of the leaf area index variable in domain file"
                allocate(dimensions(2))
                dimensions = ["Y", "X"]
                group = "Domain"
            case ("canwat_var")
                description = "Name of the canopy water variable in domain file"
                allocate(dimensions(2))
                dimensions = ["Y", "X"]
                units = "mm"
                group = "Domain"
            case ("sinalpha_var")
                description = "Name of the sine of the slope angle variable in domain file"
                allocate(dimensions(2))
                dimensions = ["Y", "X"]
                group = "Domain"
            case ("cosalpha_var")
                description = "Name of the cosine of the slope angle variable in domain file"
                allocate(dimensions(2))
                dimensions = ["Y", "X"]
                group = "Domain"
            case ("slope_var")
                description = "Name of the slope variable (0-1) in domain file"
                allocate(dimensions(2))
                dimensions = ["Y", "X"]
                group = "Domain"
            case ("aspect_angle_var")
                description = "Name of the aspect angle variable in domain file"
                allocate(dimensions(2))
                dimensions = ["Y", "X"]
                units = "radians"
                group = "Domain"
            case ("slope_angle_var")
                description = "Name of the slope angle variable in domain file"
                allocate(dimensions(2))
                dimensions = ["Y", "X"]
                units = "radians"
                group = "Domain"
            case ("svf_var")
                description = "Name of the sky view factor variable in domain file, used for radiation_downscaling=1"
                allocate(dimensions(2))
                dimensions = ["Y", "X"]
                group = "Domain"
            case ("hlm_var")
                description = "Name of the horizon line matrix variable in domain file, used for radiation_downscaling=1"
                allocate(dimensions(3))
                dimensions = ["a", "Y", "X"]
                group = "Domain"
            case ("shd_var")
                description = "Name of the snow holding depth variable in domain file, used for SLIDE=1 when running FSM2trans"
                allocate(dimensions(2))
                dimensions = ["Y", "X"]
                units = "m"
                group = "Domain"
            ! --------------------------------------
            ! --------------------------------------
            ! Forcing namelist variables
            ! --------------------------------------
            ! --------------------------------------
            case ("forcing_file_list")
                description = "Path to file containing list of forcing files. Can be generated using helpers/filelist_script.sh"
                default = ""
                group = "Forcing"
                type = 1
            case ("t_offset")
                description = "Offset to apply to temperature forcing data"
                min = -500
                max = 500
                default = "0.0"
                group = "Forcing"
                type = 1
            case ("inputinterval")
                description = "Interval between forcing data inputs in seconds"
                min = 1
                max = 86400
                default = "3600"
                group = "Forcing"
                type = 1
            case ("limit_rh")
                description = "Limit forcing relative humidity to 100% (T/F)"
                default = ".False."
                group = "Forcing"
                type = 1
            case ("t_is_potential")
                description = "Forcing temperature variable is potential temperature (T/F)"
                default = ".True."
                group = "Forcing"
                type = 1
            case ("qv_is_spec_humidity")
                description = "Forcing QV variable is specific humidity (T/F)"
                default = ".False."
                group = "Forcing"
                type = 1
            case ("relax_filters")
                description = "Use smoothly variying relaxation condition at domain boundary to nudge domain (T/F)"
                default = ".False."
                group = "Forcing"
                type = 1
            case ("qv_is_relative_humidity")
                description = "Forcing QV variable is relative humidity (T/F)"
                default = ".False."
                group = "Forcing"
                type = 1
            case ("z_is_geopotential")
                description = "Forcing Z variable is geopotential height (T/F)"
                default = ".False."
                group = "Forcing"
                type = 1
            case ("time_varying_z")
                description = "Forcing Z variable is time varying (T/F)"
                default = ".False."
                group = "Forcing"
                type = 1
            case ("hgtvar")
                description = "Name of the terrain height variable in forcing file (REQUIRED)"
                units = "meters"
                allocate(dimensions(2))
                dimensions = ["Y", "X"]
                group = "Forcing"
                type = 1
            case ("latvar")
                description = "Name of the latitude variable in forcing file (REQUIRED)"
                units = "degrees"
                allocate(dimensions(2))
                dimensions = ["Y", "X"]
                group = "Forcing"
                type = 1
            case ("lonvar")
                description = "Name of the longitude variable in forcing file (REQUIRED)"
                units = "degrees"
                allocate(dimensions(2))
                dimensions = ["Y", "X"]
                group = "Forcing"
                type = 1
            case ("time_var")
                description = "Name of the time variable in forcing file (REQUIRED)"
                allocate(dimensions(1))
                dimensions = ["T"]
                group = "Forcing"
                type = 1
            case ("uvar")
                description = "Name of the U wind variable in forcing file (REQUIRED)"
                units = "m/s"
                allocate(dimensions(4))
                dimensions = ["T", "Z", "Y", "X"]
                group = "Forcing"
                type = 1
            case ("ulat")
                description = "Name of the latitude variable on the staggered U-grid in forcing file"
                units = "degrees"
                allocate(dimensions(2))
                dimensions = ["Y", "X"]
                group = "Forcing"
                type = 1
            case ("ulon")
                description = "Name of the longitude variable on the staggered U-grid in forcing file"
                units = "degrees"
                allocate(dimensions(2))
                dimensions = ["Y", "X"]
                group = "Forcing"
                type = 1
            case ("vvar")
                description = "Name of the V wind variable in forcing file (REQUIRED)"
                units = "m/s"
                allocate(dimensions(4))
                dimensions = ["T", "Z", "Y", "X"]
                group = "Forcing"
                type = 1
            case ("vlat")
                description = "Name of the latitude variable on the staggered V-grid in forcing file"
                units = "degrees"
                allocate(dimensions(2))
                dimensions = ["Y", "X"]
                group = "Forcing"
                type = 1
            case ("vlon")
                description = "Name of the longitude variable on the staggered V-grid in forcing file"
                units = "degrees"
                allocate(dimensions(2))
                dimensions = ["Y", "X"]
                group = "Forcing"
                type = 1
            case ("wvar")
                description = "Name of the W wind variable in forcing file"
                units = "m/s"
                allocate(dimensions(4))
                dimensions = ["T", "Z", "Y", "X"]
                group = "Forcing"
                type = 1
            case ("pslvar")
                description = "Name of the sea level pressure variable in forcing file"
                units = "Pa"
                allocate(dimensions(3))
                dimensions = ["T", "Y", "X"]
                group = "Forcing"
                type = 1
            case ("psvar")
                description = "Name of the surface pressure variable in forcing file"
                units = "Pa"
                allocate(dimensions(3))
                dimensions = ["T", "Y", "X"]
                group = "Forcing"
                type = 1
            case ("pvar")
                description = "Name of the pressure variable in forcing file"
                units = "Pa"
                allocate(dimensions(4))
                dimensions = ["T", "Z", "Y", "X"]
                group = "Forcing"
                type = 1
            case ("pbvar")
                description = "Name of the base pressure variable in forcing file"
                units = "Pa"
                allocate(dimensions(4))
                dimensions = ["T", "Z", "Y", "X"]
                group = "Forcing"
                type = 1
            case ("phbvar")
                description = "Name of the base geopotential variable in forcing file"
                units = "m^2/s^2"
                allocate(dimensions(4))
                dimensions = ["T", "Z", "Y", "X"]
                group = "Forcing"
                type = 1
            case ("tvar")
                description = "Name of the temperature variable in forcing file (REQUIRED)"
                units = "K"
                allocate(dimensions(4))
                dimensions = ["T", "Z", "Y", "X"]
                group = "Forcing"
                type = 1
            case ("qvvar")
                description = "Name of the water vapor mixing ration variable in forcing file (REQUIRED)"
                units = "kg/kg"
                allocate(dimensions(4))
                dimensions = ["T", "Z", "Y", "X"]
                group = "Forcing"
                type = 1
            case ("qcvar")
                description = "Name of the cloud water mixing ratio variable in forcing file"
                units = "kg/kg"
                allocate(dimensions(4))
                dimensions = ["T", "Z", "Y", "X"]
                group = "Forcing"
                type = 1
            case ("qrvar")
                description = "Name of the rain water mixing ratio variable in forcing file"
                units = "kg/kg"
                allocate(dimensions(4))
                dimensions = ["T", "Z", "Y", "X"]
                group = "Forcing"
                type = 1
            case ("qivar")
                description = "Name of the ice water mixing ratio variable in forcing file"
                units = "kg/kg"
                allocate(dimensions(4))
                dimensions = ["T", "Z", "Y", "X"]
                group = "Forcing"
                type = 1
            case ("qgvar")
                description = "Name of the graupel mixing ratio variable in forcing file"
                units = "kg/kg"
                allocate(dimensions(4))
                dimensions = ["T", "Z", "Y", "X"]
                group = "Forcing"
                type = 1
            case ("qsvar")
                description = "Name of the snow mixing ratio variable in forcing file"
                units = "kg/kg"
                allocate(dimensions(4))
                dimensions = ["T", "Z", "Y", "X"]
                group = "Forcing"
                type = 1
            case ("i2mvar")
                description = "Name of the ice2 mixing ration variable in forcing file"
                units = "kg/kg"
                allocate(dimensions(4))
                dimensions = ["T", "Z", "Y", "X"]
                group = "Forcing"
                type = 1
            case ("i3mvar")
                description = "Name of the ice3 mixing ration variable in forcing file"
                units = "kg/kg"
                allocate(dimensions(4))
                dimensions = ["T", "Z", "Y", "X"]
                group = "Forcing"
                type = 1
            case ("qncvar")
                description = "Name of the cloud number concentration variable in forcing file"
                units = "(kg^-1)"
                allocate(dimensions(4))
                dimensions = ["T", "Z", "Y", "X"]
                group = "Forcing"
                type = 1
            case ("qnrvar")
                description = "Name of the rain number concentration variable in forcing file"
                units = "(kg^-1)"
                allocate(dimensions(4))
                dimensions = ["T", "Z", "Y", "X"]
                group = "Forcing"
                type = 1
            case ("qnivar")
                description = "Name of the ice number concentration variable in forcing file"
                units = "(kg^-1)"
                allocate(dimensions(4))
                dimensions = ["T", "Z", "Y", "X"]
                group = "Forcing"
                type = 1
            case ("qngvar")
                description = "Name of the graupel number concentration variable in forcing file"
                units = "(kg^-1)"
                allocate(dimensions(4))
                dimensions = ["T", "Z", "Y", "X"]
                group = "Forcing"
                type = 1
            case ("qnsvar")
                description = "Name of the snow number concentration variable in forcing file"
                units = "(kg^-1)"
                allocate(dimensions(4))
                dimensions = ["T", "Z", "Y", "X"]
                group = "Forcing"
                type = 1
            case ("i2nvar")
                description = "Name of the ice2 number concentration variable in forcing file"
                units = "(kg^-1)"
                allocate(dimensions(4))
                dimensions = ["T", "Z", "Y", "X"]
                group = "Forcing"
                type = 1
            case ("i3nvar")
                description = "Name of the ice3 number concentration variable in forcing file"
                units = "(kg^-1)"
                allocate(dimensions(4))
                dimensions = ["T", "Z", "Y", "X"]
                group = "Forcing"
                type = 1
            case ("i1avar")
                description = "Name of the ice1 volume mixing ratio variable in forcing file"
                units = "(m^3/kg)"
                allocate(dimensions(4))
                dimensions = ["T", "Z", "Y", "X"]
                group = "Forcing"
                type = 1
            case ("i2avar")
                description = "Name of the ice2 volume mixing ratio variable in forcing file"
                units = "(m^3/kg)"
                allocate(dimensions(4))
                dimensions = ["T", "Z", "Y", "X"]
                group = "Forcing"
                type = 1
            case ("i3avar")
                description = "Name of the ice3 volume mixing ratio variable in forcing file"
                units = "(m^3/kg)"
                allocate(dimensions(4))
                dimensions = ["T", "Z", "Y", "X"]
                group = "Forcing"
                type = 1
            case ("i1cvar")
                description = "Name of the ice1 volume x aspect ratio mixing ratio variable in forcing file"
                units = "m^3/kg"
                allocate(dimensions(4))
                dimensions = ["T", "Z", "Y", "X"]
                group = "Forcing"
                type = 1
            case ("i2cvar")
                description = "Name of the ice2 volume x aspect ratio mixing ratio variable in forcing file"
                units = "m^3/kg"
                allocate(dimensions(4))
                dimensions = ["T", "Z", "Y", "X"]
                group = "Forcing"
                type = 1
            case ("i3cvar")
                description = "Name of the ice3 volume x aspect ratio mixing ratio variable in forcing file"
                units = "m^3/kg"
                allocate(dimensions(4))
                dimensions = ["T", "Z", "Y", "X"]
                group = "Forcing"
                type = 1
            case ("zvar")
                description = "Name of the vertical height variable in forcing file"
                units = "m"
                allocate(dimensions(4))
                dimensions = ["T", "Z", "Y", "X"]
                group = "Forcing"
                type = 1
            case ("shvar")
                description = "Name of the sensible heat flux variable in forcing file"
                units = "W/m^2"
                allocate(dimensions(3))
                dimensions = ["T", "Y", "X"]
                group = "Forcing"
                type = 1
            case ("lhvar")
                description = "Name of the latent heat flux variable in forcing file"
                units = "W/m^2"
                allocate(dimensions(3))
                dimensions = ["T", "Y", "X"]
                group = "Forcing"
                type = 1
            case ("swdown_var")
                description = "Name of the shortwave downwelling radiation variable in forcing file"
                units = "W/m^2"
                allocate(dimensions(3))
                dimensions = ["T", "Y", "X"]
                group = "Forcing"
                type = 1
            case ("lwdown_var")
                description = "Name of the longwave downwelling radiation variable in forcing file"
                units = "W/m^2"
                allocate(dimensions(3))
                dimensions = ["T", "Y", "X"]
                group = "Forcing"
                type = 1
            case ("sst_var")
                description = "Name of the sea surface temperature variable in forcing file"
                units = "K"
                allocate(dimensions(3))
                dimensions = ["T", "Y", "X"]
                group = "Forcing"
                type = 1
            case ("pblhvar")
                description = "Name of the planetary boundary layer height variable in forcing file"
                units = "m"
                allocate(dimensions(3))
                dimensions = ["T", "Y", "X"]
                group = "Forcing"
                type = 1
            ! --------------------------------------
            ! --------------------------------------
            ! Restart namelist variables
            ! --------------------------------------
            ! --------------------------------------
            case ("restart_run")
                description = "Flag to control if this is a restart run or not (T/F)"
                default = '.False.'
                group = "Restart"
                type = 1
            case ("restart_folder")
                description = "Path to folder where restart files will be read from and output to"
                default = "../restart/"
                group = "Restart"
                type = 1
            case ("restartinterval")
                description = "Frequency of writing restart files, expressed in number of outputintervals for the first domain."//achar(10)//BLNK_CHR_N// &
                              "Must be set for all domains in a series of nests, such that all nests output restart files at the same time."//achar(10)//BLNK_CHR_N// &
                              "If it is not set for subequent nests, HICAR will try to calculate a restartinterval for domains"//achar(10)//BLNK_CHR_N// &
                              "that results in the same restart write time as for the first domain."
                min = 1
                max = 1000
                units ='number of output steps'
                default = "24"
                group = "Restart"
            case ("restart_date")
                description = "Date to restart from, format: 'YYYY-MM-DD HH:MM:SS'"
                default = ""
                group = "Restart"
                type = 1
            case ("restart_step")
                description = "Step in restart_in_file to restart from"
                default = "1"
                group = "Restart"
                type = 1
            ! --------------------------------------
            ! --------------------------------------
            ! Output namelist variables
            ! --------------------------------------
            ! --------------------------------------
            case ("output_folder")
                description = "Path to output folder"
                default = "../output/"
                group = "Output"
                type = 1
            case ("outputinterval")
                description = "Frequency of writing output files in seconds"
                min = 1
                max = 1000000
                units = 'seconds'
                default = "3600"
                group = "Output"
            case ("frames_per_outfile")
                description = "Number of output steps to write to each output file"
                min = 1
                max = 100
                default = "24"
                group = "Output"
            case ("output_vars")
                description = "List of variables to output. Call HICAR as: './HICAR --out-vars keyword' to get information"//achar(10)//BLNK_CHR_N// &
                                "about the available output variables matching the keyword."
                default = ""
                group = "Output"
            ! --------------------------------------
            ! --------------------------------------
            ! Physics parameterisations variables
            ! --------------------------------------
            ! --------------------------------------
            case ("pbl")
                description = "Planetary boundary layer scheme to use: "//achar(10)//BLNK_CHR_N// &
                                                                       "0 = no PBL,"//achar(10)//BLNK_CHR_N// &
                                                                       "1 = YSU PBL"
                allocate(values(2))
                values = [0, 1]
                default = "0"
                group = "Physics"
            case ("lsm")
                description = "Land surface model to use:"//achar(10)//BLNK_CHR_N// &
                                                        "0 = no LSM,"//achar(10)//BLNK_CHR_N// &
                                                        "1 = Fluxes from forcing data,"//achar(10)//BLNK_CHR_N// &
                                                        "2 = Noah LSM"//achar(10)//BLNK_CHR_N// &
                                                        "3 = Noah MP"

                allocate(values(4))
                values = [0, 1, 2, 3]
                default = "0"
                group = "Physics"
                type = 1
            case ("rad")
                description = "Radiation scheme to use: "//achar(10)//BLNK_CHR_N// &
                                                       "0 = no RAD,"//achar(10)//BLNK_CHR_N// &
                                                       "1 = Surface fluxes from forcing data"//achar(10)//BLNK_CHR_N// &
                                                       "2 = cloud fraction based radiation + radiative cooling"//achar(10)//BLNK_CHR_N// &
                                                       "3 = RRTMG"

                allocate(values(4))
                values = [0, 1, 2, 3]
                default = "0"
                group = "Physics"
            case ("conv")
                description = "Cumulus scheme to use: "//achar(10)//BLNK_CHR_N// &
                                                     "0 = no CONV,"//achar(10)//BLNK_CHR_N// &
                                                     "1 = Tiedke scheme"//achar(10)//BLNK_CHR_N// &
                                                     "2 = NSAS scheme"//achar(10)//BLNK_CHR_N// &
                                                     "3 = BMJ scheme"
                allocate(values(4))
                values = [0, 1, 2, 3]
                default = "0"
                group = "Physics"
            case ("mp")
                description = "Microphysics scheme to use: "//achar(10)//BLNK_CHR_N// &
                                                          "0 = no MP,"//achar(10)//BLNK_CHR_N// &
                                                          "1 = Thompson et al (2008),"//achar(10)//BLNK_CHR_N// &
                                                          "2 = 'Linear' microphysics"//achar(10)//BLNK_CHR_N// &
                                                          "3 = Morrison"//achar(10)//BLNK_CHR_N// &
                                                          "4 = WSM6"//achar(10)//BLNK_CHR_N// &
                                                          "5 = Thompson Aerosol"//achar(10)//BLNK_CHR_N// &
                                                          "6 = WSM3"//achar(10)//BLNK_CHR_N// &
                                                          "7 = ISHMAEL"
                allocate(values(8))
                values = [0, 1, 2, 3, 4, 5, 6, 7]
                default = "0"
                group = "Physics"
                type = 1
            case ("water")
                description = "Water model to use:  "//achar(10)//BLNK_CHR_N// &
                                                   "0 = no open water fluxes,"//achar(10)//BLNK_CHR_N// &
                                                   "1 = Simple fluxes (needs SST in forcing data)"//achar(10)//BLNK_CHR_N// &
                                                   "2 = WRF's lake model (needs lake depth in hi-res data))"
                allocate(values(3))
                values = [0, 1, 2]
                default = "0"
                group = "Physics"
            case ("wind")
                description = "Wind solver to use: "//achar(10)//BLNK_CHR_N// &
                                                   "0 = no LT,"//achar(10)//BLNK_CHR_N// &
                                                   "1 = linear theory wind perturbations"//achar(10)//BLNK_CHR_N// &
                                                   "2 = Adjustment to horizontal winds to reduce divergence, based on technique from O'brien et al., 1970"//achar(10)//BLNK_CHR_N// &
                                                   "3 = Mass-conserving wind solver based on variational calculus technique, requires PETSc"//achar(10)//BLNK_CHR_N// &
                                                   "4 = Combination of options 1 & 2"//achar(10)//BLNK_CHR_N// &
                                                   "5 = Combination of options 1 & 3"
                allocate(values(6))
                values = [0, 1, 2, 3, 4, 5]
                default = "0"
                group = "Physics"
            case ("radiation_downscaling")
                description = "0 = no downcaling"//achar(10)//BLNK_CHR_N//      &
                              "1 = terrain shading effect is considered in the radiation calculation"
                default = "0"
                allocate(values(2))
                values = [0, 1]
                group = "Physics"
            case ("adv")
                description = "Advection scheme to use:  "//achar(10)//BLNK_CHR_N// &
                                                        "0 = no ADV,"//achar(10)//BLNK_CHR_N// &
                                                        "1 = standard advection scheme"//achar(10)//BLNK_CHR_N// &
                                                        "2 = MPDATA"
                allocate(values(3))
                values = [0, 1, 2]
                default = "0"
                group = "Physics"
            case ("sfc")
                description = "Surface model to use: "//achar(10)//BLNK_CHR_N// &
                                                    "0 = no surface layer"//achar(10)//BLNK_CHR_N// &
                                                    "1 = Revised MM5 Monin-Obukhov scheme"
                allocate(values(2))
                values = [0, 1]
                default = "0"
                group = "Physics"
            case ("sm")
                description = "Snow model to use: "//achar(10)//BLNK_CHR_N// &
                                                 "0 = no snow model"//achar(10)//BLNK_CHR_N// &
                                                 "1 = FSM2trans snow model (must be compiled, see docs/compiling.md)"
                allocate(values(2))
                values = [0, 1]
                default = "0"
                group = "Physics"
            ! --------------------------------------
            ! --------------------------------------
            ! MP parameters namelist variables
            ! --------------------------------------
            ! --------------------------------------
            case ("Nt_c")
                description = "Prescribed number of cloud droplets"
                min = 0
                max = 100.e7
                units = ""
                default = "100.e6"
                group = "MP_Parameters"
            case ("TNO")
                description = "Constants in Cooper curve relation for cloud ice number"
                min = 0
                max = 100.
                units = ""
                default = "5.0"
                group = "MP_Parameters"
            case ("am_s")
                description = "Mass power law relation parameter for snow"
                min = 0
                max = 0.1
                units = ""
                default = "0.069"
                group = "MP_Parameters"
            case ("rho_g")
                description = "Density of graupel"
                min = 0
                max = 1000
                units = ""
                default = "500.0"
                group = "MP_Parameters"
            case ("av_s")
                description = "Fallspeed power law parameter for snow"
                min = 0
                max = 100.
                units = ""
                default = "40."
                group = "MP_Parameters"
            case ("bv_s")
                description = "Fallspeed power law parameter for snow"
                min = 0
                max = 1.
                units = ""
                default = "0.55"
                group = "MP_Parameters"
            case ("fv_s")
                description = "Fallspeed power law parameter for snow"
                min = 0
                max = 500.
                units = ""
                default = "100."
                group = "MP_Parameters"
            case ("av_g")
                description = "Fallspeed power law parameter for graupel"
                min = 0
                max = 1000.
                units = ""
                default = "442."
                group = "MP_Parameters"
            case ("bv_g")
                description = "Fallspeed power law parameter for graupel"
                min = 0
                max = 1.
                units = ""
                default = "0.89"
                group = "MP_Parameters"
            case ("av_i")
                description = "Fallspeed power law parameter for ice"
                min = 0
                max = 5000
                units = ""
                default = "1847.5"
                group = "MP_Parameters"
            case ("Ef_si")
                description = "Collection efficiency snow/ice"
                min = 0
                max = 1
                units = "(-)"
                default = "0.05"
                group = "MP_Parameters"
            case ("Ef_rs")
                description = "Collection efficiency rain/snow"
                min = 0
                max = 1.
                units = "(-)"
                default = "0.95"
                group = "MP_Parameters"
            case ("Ef_rg")
                description = "Collection efficiency rain/graupel"
                min = 0
                max = 1.
                units = "(-)"
                default = "0.75"
                group = "MP_Parameters"
            case ("Ef_ri")
                description = "Collection efficiency rain/ice"
                min = 0
                max = 1.
                units = "(-)"
                default = "0.95"
                group = "MP_Parameters"
            case ("C_cubes")
                description = "Capacitance of sphere and plates/aggregates D**3"
                min = 0
                max = 1.
                units = ""
                default = "0.5"
                group = "MP_Parameters"
            case ("C_sqrd")
                description = "Capacitance of sphere and plates/aggregates D**2"
                min = 0
                max = 1.
                units = ""
                default = "0.3"
                group = "MP_Parameters"
            case ("mu_r")
                description = "Generalized gamma distribution parameter for rain"
                min = 0
                max = 100
                units = ""
                default = "0.0"
                group = "MP_Parameters"
            case ("t_adjust")
                description = "trude, add tadjust, so to chane temperature for where Bigg freezing starts."//achar(10)//BLNK_CHR_N// &
                    "Follow approach in WRFV3.6 with IN"
                min = -100
                max = 100
                units = ""
                default = "0.0"
                group = "MP_Parameters"
            case ("Ef_rw_l")
                description = "True sets ef_rw = 1"
                default = ".False."
                group = "MP_Parameters"
            case ("Ef_sw_l")
                description = "True sets ef_sw = 1"
                default = ".False."
                group = "MP_Parameters"
            case ("update_interval_mp")
                description = "Time interval for updating the microphysics. If = 0, update every time step. = 0 is recommended."
                min = 0
                max = 7200
                units = "seconds"
                default = "0.0"
                group = "MP_Parameters"
            case ("top_mp_level")
                description = "Highest model level, measured from the ground up, up to which microphysics should be run."//achar(10)//BLNK_CHR_N// &
                    "If <= 0, run microphysics on all levels."
                min = 0
                max = MAXLEVELS
                default = "0"
                group = "MP_Parameters"
            ! --------------------------------------
            ! --------------------------------------
            ! Linear theory parameters namelist variables
            ! --------------------------------------
            ! --------------------------------------
            case ("variable_N")
                description = "Compute the Brunt Vaisala Frequency (N^2) every time step (T/F) (DEPRECATED, always variable now)"
                default = ".True."
                group = "LT_Parameters"
            case ("smooth_nsq")
                description = "Smooth the Brunt-Vaisala frequency over vert_smooth vertical levels (T/F)"
                default = ".True."
                group = "LT_Parameters"
            case ("vert_smooth")
                description = "Number of vertical levels over which to smooth the Brunt-Vaisala frequency"
                min = 0
                max = MAXLEVELS
                default = "10"
                group = "LT_Parameters"
            case ("buffer")
                description = "Number of vertical levels to buffer around the domain MUST be >=1"
                min = 0
                max = MAXLEVELS
                default = "50"
                group = "LT_Parameters"
            case ("stability_window_size")
                description = "window to average nsq over"
                min = 0
                max = 100
                default = "10"
                group = "LT_Parameters"
            case ("max_stability")
                description = "max limit on the calculated Brunt Vaisala Frequency"
                min = 0
                max = 1
                units = "s^-2"
                default = "6e-4" 
                group = "LT_Parameters"
            case ("min_stability")
                description = "min limit on the calculated Brunt Vaisala Frequency"
                min = 0
                max = 1
                units = "s^-2"
                default = "1e-7"
                group = "LT_Parameters"
            case ("N_squared")
                description = "static Brunt Vaisala Frequency (N^2) to use (DEPRECATED)"
                min = -10
                max = 10
                units = "s^-2"
                default = "3e-5"
                group = "LT_Parameters"
            case ("linear_contribution")
                description = "fractional contribution of linear perturbation to wind field (e.g. u_hat multiplied by this)"
                min = 0
                max = 1
                default = "1.0"
                group = "LT_Parameters"
            case ("remove_lowres_linear")
                description = "attempt to remove the linear mountain wave from the forcing low res model (T/F)"
                default = ".False."
                group = "LT_Parameters"
            case ("rm_N_squared")
                description = "static Brunt Vaisala Frequency (N^2) to use in removing linear wind field"
                min = 0
                max = 1
                default = "3e-5"
                group = "LT_Parameters"
            case ("rm_linear_contribution")
                description = "fractional contribution of linear perturbation to wind field to remove from the low-res field"
                min = 0
                max = 1
                default = "1.0"
                group = "LT_Parameters"
            case ("linear_update_fraction")
                description = "fraction of linear perturbation to add each time step"
                min = 0
                max = 1
                default = "0.2"
                group = "LT_Parameters"
            case ("dirmax")
                description = "maximum direction of the wind perturbation look up table"
                min = 0.
                max = 2.*piconst
                units = "radians"
                default = "6.283"
                group = "LT_Parameters"
            case ("dirmin")
                description = "minimum direction of the wind perturbation look up table"
                min = 0
                max = 2*piconst
                units = "radians"
                default = "0.0"
                group = "LT_Parameters"
            case ("spdmax")
                description = "maximum speed of the wind perturbation look up table"
                min = 0
                max = 100
                units = "m/s"
                default = "30.0"
                group = "LT_Parameters"
            case ("spdmin")
                description = "minimum speed of the wind perturbation look up table"
                min = 0
                max = 100
                units = "m/s"
                default = "0.0"
                group = "LT_Parameters"
            case ("nsqmax")
                description = "maximum Brunt Vaisala Frequency (N^2) of the wind perturbation look up table"
                min = -100
                max = 100
                units = "s^-2"
                default = "-3.2218487496"
                group = "LT_Parameters"
            case ("nsqmin")
                description = "minimum Brunt Vaisala Frequency (N^2) of the wind perturbation look up table"
                min = -100
                max = 100
                units = "s^-2"
                default = "-7"
                group = "LT_Parameters"
            case ("n_dir_values")
                description = "number of direction values in the wind perturbation look up table"
                min = 0
                max = 100
                default = "24"
                group = "LT_Parameters"
            case ("n_spd_values")
                description = "number of speed values in the wind perturbation look up table"
                min = 0
                max = 100
                default = "6"
                group = "LT_Parameters"
            case ("n_nsq_values")
                description = "number of Brunt Vaisala Frequency (N^2) values in the wind perturbation look up table"
                min = 0
                max = 50
                default = "5"
                group = "LT_Parameters"
            case ("minimum_layer_size")
                description = "minimum layer size for the linear wind perturbation"
                min = 0
                max = 1000
                units = "m"
                default = "100"
                group = "LT_Parameters"
            case ("read_LUT")
                description = "read in a look up table for the linear wind perturbation (T/F)"
                default = ".False."
                group = "LT_Parameters"
            case ("write_LUT")
                description = "write out a look up table for the linear wind perturbation (T/F)"
                default = ".True."
                group = "LT_Parameters"
            case ("u_LUT_Filename")
                description = "filename for the u component of the wind perturbation look up table"
                default = ""
                group = "LT_Parameters"
            case ("v_LUT_Filename")
                description = "filename for the v component of the wind perturbation look up table"
                default = ""
                group = "LT_Parameters"
            case ("LUT_Filename")
                description = "filename for the wind perturbation look up table"
                default = ""
                group = "LT_Parameters"
            case ("overwrite_lt_lut")
                description = "overwrite the linear theory look up table (T/F)"
                default = ".True."
                group = "LT_Parameters"
            ! --------------------------------------
            ! --------------------------------------
            ! Advection parameters namelist variables
            ! --------------------------------------
            ! --------------------------------------
            case ("boundary_buffer")
                description = "apply some smoothing to the x and y boundaries in MPDATA (T/F)"
                default = ".False."
                group = "ADV_Parameters"
            case ("advect_density")
                description = "use air density in the advection equations (T/F)"
                default = ".True."
                group = "ADV_Parameters"
            case ("MPDATA_FCT")
                description = "use the flux corrected transport scheme in MPDATA (T/F)"
                default = ".True."
                group = "ADV_Parameters"
            case ("mpdata_order")
                description = "order of the MPDATA scheme"
                min = 1
                max = 6
                default = "2"
                group = "ADV_Parameters"
            case ("flux_corr")
                description = "flux correction scheme to use for standard advection. Recommended use with RK3=.True."//achar(10)//BLNK_CHR_N// &
                    "0=None, 1=WRF positive definite monotonic flux correction"
                allocate(values(2))
                values = [0, 1]
                default = "0"
                group = "ADV_Parameters"
            case ("h_order")
                description = "order of the horizontal advection scheme"
                allocate(values(3))
                values = [1,3,5]
                default = "1"
                group = "ADV_Parameters"
            case ("v_order")
                description = "order of the vertical advection scheme"
                allocate(values(3))
                values = [1,3,5]
                default = "1"
                group = "ADV_Parameters"
            ! --------------------------------------
            ! --------------------------------------
            ! PBL parameters namelist variables
            ! --------------------------------------
            ! --------------------------------------
            case ("ysu_topdown_pblmix")
                description = "use the YSU PBL scheme with radiative, top-down mixing scheme. 0=off, 1=on"
                allocate(values(2))
                values = [0, 1]
                default = "0"
                group = "PBL_Parameters"
            ! --------------------------------------
            ! --------------------------------------
            ! SFC parameters namelist variables
            ! --------------------------------------
            ! --------------------------------------
            case ("isfflx")
                description = "Use surface fluxes calculated by LSM scheme."//achar(10)//BLNK_CHR_N// &
                    "If lsm is turned on, this will be changed to 1. 0=off, 1=on"
                allocate(values(2))
                values = [0, 1]
                default = "1"
                group = "SFC_Parameters"
            case ("scm_force_flux")
                description = "If lsm is turned on, this will be changed to 1. 0=on, 1=off"
                allocate(values(2))
                values = [0, 1]
                default = "0"
                group = "SFC_Parameters"
            case ("iz0tlnd")
                description = "thermal roughness length for sfclay. (0 = old, 1 = veg dependent Chen-Zhang Czil)"
                allocate(values(2))
                values = [0, 1]
                default = "1"
                group = "SFC_Parameters"
            case ("isftcflx")
                description = "not sure what this does, controls roughness length over water I think?"
                allocate(values(2))
                values = [0, 1]
                default = "0"
                group = "SFC_Parameters"
            case ("sbrlim")
                description = "parameter to limit the bulk richardson number under stable conditions"
                min = 0.
                max = 1000.
                default = "250."
                group = "SFC_Parameters"
            ! --------------------------------------
            ! --------------------------------------
            ! CU parameters namelist variables
            ! --------------------------------------
            ! --------------------------------------
            case ("tendency_fraction")
                description = "fraction of the tendency to use in the cumulus scheme"
                min = 0
                max = 1
                default = "1.0"
                group = "CU_Parameters"
            case ("tend_qv_fraction")
                description = "fraction of the tendency to use in the cumulus scheme for water vapor."//achar(10)//BLNK_CHR_N// &
                    "If < 0, use the value of 'tendency_fraction'"
                min = 0
                max = 1
                default = "-1.0"
                group = "CU_Parameters"
            case ("tend_qc_fraction")
                description = "fraction of the tendency to use in the cumulus scheme for cloud water."//achar(10)//BLNK_CHR_N// &
                    "If < 0, use the value of 'tendency_fraction'"
                min = 0
                max = 1
                default = "-1.0"
                group = "CU_Parameters"
            case ("tend_th_fraction")
                description = "fraction of the tendency to use in the cumulus scheme for potential temperature."//achar(10)//BLNK_CHR_N// &
                    "If < 0, use the value of 'tendency_fraction'"
                min = 0
                max = 1
                default = "-1.0"
                group = "CU_Parameters"
            case ("tend_qi_fraction")
                description = "fraction of the tendency to use in the cumulus scheme for ice."//achar(10)//BLNK_CHR_N// &
                    "If < 0, use the value of 'tendency_fraction'"
                min = 0
                max = 1
                default = "-1.0"
                group = "CU_Parameters"
            case ("stochastic_cu")
                description = "disturbes the W field (randomly; higher value=more disturbance). Triggers convection."//achar(10)//BLNK_CHR_N// &
                    "If set to -9999, no convective perturbation is applied."
                min = -9999
                max = 1000
                default = "-9999"
                group = "CU_Parameters"
            ! --------------------------------------
            ! --------------------------------------
            ! LSM parameters namelist variables
            ! --------------------------------------
            ! --------------------------------------
            case ("LU_Categories")
                description = "Land use category clasification to use in the LSM scheme."//achar(10)//BLNK_CHR_N// &
                    "Values: MODIFIED_IGBP_MODIS_NOAH, USGS, USGS-RUC, MODI-RUC, NLCD40"
                default = "MODIFIED_IGBP_MODIS_NOAH"
                group = "LSM_Parameters"
            case ("update_interval_lsm")
                description = "Time interval for updating the LSM. If = 0, update every time step."
                min = 0
                max = 7200
                units = "seconds"
                default = "300.0"
                group = "LSM_Parameters"
            case ("monthly_vegfrac")
                description = "Use monthly vegetation fraction data (T/F)"
                default = ".False."
                group = "LSM_Parameters"
            case ("monthly_albedo")
                description = "Use monthly albedo data (T/F)"
                default = ".False."
                group = "LSM_Parameters"
            case ("num_soil_layers")
                description = "Number of soil layers in the LSM"
                default = "4"
                min = 1
                max = 20
                group = "LSM_Parameters"
            case ("urban_category")
                description = "Land use category corresponding to urban classification."//achar(10)//BLNK_CHR_N// &
                    "Setting to -1 uses default value for the given land use classification given in 'LU_Categories'. "
                min = -1
                max = 100
                default = "-1"
                group = "LSM_Parameters"
            case ("ice_category")
                description = "Land use category corresponding to ice classification."//achar(10)//BLNK_CHR_N// &
                    "Setting to -1 uses default value for the given land use classification given in 'LU_Categories'. "
                min = -1
                max = 100
                default = "-1"
                group = "LSM_Parameters"
            case ("water_category")
                description = "Land use category corresponding to open water classification."//achar(10)//BLNK_CHR_N// &
                    "Setting to -1 uses default value for the given land use classification given in 'LU_Categories'. "
                min = -1
                max = 100
                default = "-1"
                group = "LSM_Parameters"
            case ("lake_category")
                description = "Land use category corresponding to lake classification."//achar(10)//BLNK_CHR_N// &
                    "Setting to -1 uses default value for the given land use classification given in 'LU_Categories'. "
                min = -1
                max = 100
                default = "-1"
                group = "LSM_Parameters"
            case ("sf_urban_phys")
                description = "Urban physics parameterization to use for Noah LSM (not enabled in code)."
                allocate(values(3))
                values = [0, 1, 2]
                default = "1"
                group = "LSM_Parameters"
            case ("nmp_dveg")
                description = "Dynamic vegetation type for Noah MP of vegetation types in the LSM."//achar(10)//BLNK_CHR_N// &
                                                     "1 = Off (LAI from table; FVEG = shdfac)"//achar(10)//BLNK_CHR_N// &
                                                     "2 = On  (LAI predicted;  FVEG calculated)"//achar(10)//BLNK_CHR_N// &
                                                     "3 = Off (LAI from table; FVEG calculated)"//achar(10)//BLNK_CHR_N// &
                                                     "4 = Off (LAI from table; FVEG = maximum veg. fraction)"//achar(10)//BLNK_CHR_N// &
                                                     "5 = On  (LAI predicted;  FVEG = maximum veg. fraction)"//achar(10)//BLNK_CHR_N// &
                                                     "6 = On  (use FVEG = SHDFAC from input)"//achar(10)//BLNK_CHR_N// &
                                                     "7 = Off (use input LAI; use FVEG = SHDFAC from input)"//achar(10)//BLNK_CHR_N// &
                                                     "8 = Off (use input LAI; calculate FVEG)"//achar(10)//BLNK_CHR_N// &
                                                     "9 = Off (use input LAI; use maximum vegetation fraction)"
                allocate(values(9))
                values = [1, 2, 3, 4, 5, 6, 7, 8, 9]
                default = "3"
                group = "LSM_Parameters"
            case ("nmp_opt_crs")
                description = "Noah-MP Stomatal Resistance option:"//achar(10)//BLNK_CHR_N// &
                                                    "1 = Ball-Berry"//achar(10)//BLNK_CHR_N// &
                                                    "2 = Jarvis"
                allocate(values(2))
                values = [1, 2]
                default = "1"
                group = "LSM_Parameters"
            case ("nmp_opt_sfc")
                description = "Noah-MP surface layer drag coefficient calculation"//achar(10)//BLNK_CHR_N// &
                                                     "1 = Monin-Obukhov"//achar(10)//BLNK_CHR_N// &
                                                     "2 = original Noah (Chen97)"//achar(10)//BLNK_CHR_N// &
                                                     "3 = YSU consistent"
                allocate(values(3))
                values = [-1, 1, 2, 3]
                default = "1"
                group = "LSM_Parameters"
            case ("nmp_opt_btr")
                description = "Noah-MP Soil Moisture Factor for Stomatal Resistance"//achar(10)//BLNK_CHR_N// &
                                                     "1 = Noah (soil moisture)"//achar(10)//BLNK_CHR_N// &
                                                     "2 = CLM (matric potential)"//achar(10)//BLNK_CHR_N// &
                                                     "3 = SSiB (matric potential)"
                allocate(values(3))
                values = [1, 2, 3]
                default = "2"
                group = "LSM_Parameters"
            case ("nmp_opt_run")
                description = "Noah-MP Runoff and Groundwater option"//achar(10)//BLNK_CHR_N// &
                                                     "1 = TOPMODEL with groundwater"//achar(10)//BLNK_CHR_N// &
                                                     "2 = TOPMODEL with equilibrium water table"//achar(10)//BLNK_CHR_N// &
                                                     "3 = original surface and subsurface runoff (free drainage)"//achar(10)//BLNK_CHR_N// &
                                                     "4 = BATS surface and subsurface runoff (free drainage)"//achar(10)//BLNK_CHR_N// &
                                                     "5 = Miguez-Macho&Fan groundwater scheme (Miguez-Macho et al. 2007 JGR; Fan et al. 2007 JGR)"
                allocate(values(5))
                values = [1, 2, 3, 4, 5]
                default = "1"
                group = "LSM_Parameters"
            case ("nmp_opt_frz")
                description = "Noah-MP Supercooled Liquid Water option"//achar(10)//BLNK_CHR_N// &
                                                     "1 = No iteration"//achar(10)//BLNK_CHR_N// &
                                                     "2 = Koren's iteration"
                allocate(values(2))
                values = [1, 2]
                default = "1"
                group = "LSM_Parameters"
            case ("nmp_opt_inf")
                description = "Noah-MP Soil Permeability option"//achar(10)//BLNK_CHR_N// &
                                                     "1 = Linear effects, more permeable"//achar(10)//BLNK_CHR_N// &
                                                     "2 = Non-linear effects, less permeable"
                allocate(values(2))
                values = [1, 2]
                default = "1"
                group = "LSM_Parameters"
            case ("nmp_opt_rad")
                description = "Noah-MP Radiative Transfer option"//achar(10)//BLNK_CHR_N// &
                                                     "1 = Modified two-stream (known to cause problems when vegetation fraction is small)"//achar(10)//BLNK_CHR_N// &
                                                     "2 = Two-stream applied to grid-cell"//achar(10)//BLNK_CHR_N// &
                                                     "3 = Two-stream applied to vegetated fraction"
                allocate(values(3))
                values = [1, 2, 3]
                default = "3"
                group = "LSM_Parameters"
            case ("nmp_opt_alb")
                description = "Noah-MP Ground Surface Albedo option"//achar(10)//BLNK_CHR_N// &
                                                     "1 = BATS"//achar(10)//BLNK_CHR_N// &
                                                     "2 = CLASS"
                allocate(values(2))
                values = [1, 2]
                default = "2"
                group = "LSM_Parameters"
            case ("nmp_opt_snf")
                description = "Noah-MP Precipitation Partitioning between snow and rain"//achar(10)//BLNK_CHR_N// &
                                                     "1 = Jordan (1991)"//achar(10)//BLNK_CHR_N// &
                                                     "2 = BATS:  Snow when SFCTMP < TFRZ+2.2"//achar(10)//BLNK_CHR_N// &
                                                     "3 = Snow when SFCTMP < TFRZ"//achar(10)//BLNK_CHR_N// &
                                                     "4 = Use WRF precipitation partitioning"
                allocate(values(4))
                values = [1, 2, 3, 4]
                default = "4"
                group = "LSM_Parameters"
            case ("nmp_opt_tbot")
                description = "Noah-MP Soil Temperature Lower Boundary Condition"//achar(10)//BLNK_CHR_N// &
                                                     "1 = Zero heat flux"//achar(10)//BLNK_CHR_N// &
                                                     "2 = TBOT at 8 m from input file"
                allocate(values(2))
                values = [1, 2]
                default = "2"
                group = "LSM_Parameters"
            case ("nmp_opt_stc")
                description = "Noah-MP Snow/Soil temperature time scheme"//achar(10)//BLNK_CHR_N// &
                                                     "1 = semi-implicit"//achar(10)//BLNK_CHR_N// &
                                                     "2 = full-implicit"//achar(10)//BLNK_CHR_N// &
                                                     "3 = semi-implicit where Ts uses snow cover fraction"
                allocate(values(3))
                values = [1, 2, 3]
                default = "1"
                group = "LSM_Parameters"
            case ("nmp_opt_gla")
                description = "Noah-MP glacier treatment option"//achar(10)//BLNK_CHR_N// &
                                                     "1 = includes phase change"//achar(10)//BLNK_CHR_N// &
                                                     "2 = slab ice (Noah)"
                allocate(values(2))
                values = [1, 2]
                default = "1"
                group = "LSM_Parameters"
            case ("nmp_opt_rsf")
                description = "Noah-MP surface evaporation resistance option"//achar(10)//BLNK_CHR_N// &
                                                     "1 = Sakaguchi and Zeng, 2009"//achar(10)//BLNK_CHR_N// &
                                                     "2 = Sellers (1992)"//achar(10)//BLNK_CHR_N// &
                                                     "3 = adjusted Sellers to decrease RSURF for wet soil"//achar(10)//BLNK_CHR_N// &
                                                     "4 = option 1 for non-snow; rsurf = rsurf_snow for snow (set in MPTABLE)"
                allocate(values(4))
                values = [1, 2, 3, 4]
                default = "1"
                group = "LSM_Parameters"
            case ("nmp_opt_soil")
                description = "Noah-MP options for defining soil properties"//achar(10)//BLNK_CHR_N// &
                                                     "1 = use input dominant soil texture"//achar(10)//BLNK_CHR_N// &
                                                     "2 = use input soil texture that varies with depth"
                allocate(values(2))
                values = [1, 2]
                default = "1"
                group = "LSM_Parameters"
            case ("nmp_opt_pedo")
                description = "Noah-MP options for pedotransfer functions (used when OPT_SOIL = 3; not implemented in code)"
                allocate(values(1))
                values = [1]
                default = "1"
                group = "LSM_Parameters"
            case ("nmp_opt_crop")
                description = "options for crop model"//achar(10)//BLNK_CHR_N// &
                                                     "0 = No crop model, will run default dynamic vegetation"//achar(10)//BLNK_CHR_N// &
                                                     "1 = Liu, et al. 2016"//achar(10)//BLNK_CHR_N// &
                                                     "2 = Gecros (Genotype-by-Environment interaction on CROp growth Simulator) Yin and van Laar, 2005"
                allocate(values(3))
                values = [0, 1, 2]
                default = "0"
                group = "LSM_Parameters"
            case ("nmp_opt_irr")
                description = "options for irrigation scheme"//achar(10)//BLNK_CHR_N// &
	                                                 "0 = No irrigation"//achar(10)//BLNK_CHR_N// &
                                                     "1 = Irrigation ON"//achar(10)//BLNK_CHR_N// &
                                                     "2 = irrigation trigger based on crop season Planting and harvesting dates"//achar(10)//BLNK_CHR_N// &
                                                     "3 = irrigation trigger based on LAI threshold"
                allocate(values(4))
                values = [0, 1, 2, 3]
                default = "0"
                group = "LSM_Parameters"
            case ("nmp_opt_irrm")
                description = "options for irrigation method (only if opt_irr > 0)"//achar(10)//BLNK_CHR_N// &
                                                     "0 = method based on geo_em fractions (all three methods are ON)"//achar(10)//BLNK_CHR_N// &
                                                     "1 = sprinkler method"//achar(10)//BLNK_CHR_N// &
                                                     "2 = micro/drip irrigation"//achar(10)//BLNK_CHR_N// &
                                                     "3 = surface flooding"
                allocate(values(4))
                values = [0, 1, 2, 3]
                default = "0"
                group = "LSM_Parameters"
            case ("nmp_opt_tdrn")
                description = "Noah-MP tile drainage option (currently only tested and works with opt_run=3)"//achar(10)//BLNK_CHR_N// &
                                                     "0 = No tile drainage"//achar(10)//BLNK_CHR_N// &
                                                     "1 = Simple drainage"//achar(10)//BLNK_CHR_N// &
                                                     "2 = Hooghoudt's equation based tile drainage"
                allocate(values(3))
                values = [0, 1, 2]
                default = "0"
                group = "LSM_Parameters"
            case ("nmp_soiltstep")
                description = "Noah-MP soil process timestep (seconds) for solving soil water and temperature"//achar(10)//BLNK_CHR_N// &
                                                     "0 = default, the same as main NoahMP model timestep"
                min = 0
                max = 10000
                units = "seconds"
                default = "0"
                group = "LSM_Parameters"
            case ("noahmp_output")
                description = "Noah-MP output levels"//achar(10)//BLNK_CHR_N// &
                                                     "1 = standard output"//achar(10)//BLNK_CHR_N// &
                                                     "3 = standard output with additional water and energy budget term output"
                allocate(values(2))
                values = [1, 3]
                default = "1"
                group = "LSM_Parameters"
            case ("max_swe")
                description = "Maximum snow water equivalent (SWE) for snow cover fraction"
                min = 0
                max = 1e11
                units = "mm"
                default = "1e10"
                group = "LSM_Parameters"
            case ("snow_den_const")
                description = "Snow density constant for back-calculating SWE or snow height"//achar(10)//BLNK_CHR_N// &
                    "If only one given as domain input variable."
                min = 0
                max = 1e11
                units = "kg/m^3"
                default = "100."
                group = "LSM_Parameters"
            ! --------------------------------------
            ! --------------------------------------
            ! Snow Model parameters namelist variables
            ! --------------------------------------
            ! --------------------------------------
            case ("fsm_nsnow_max")
                description = "Maximum number of snow layers to allow for a FSM simulation"
                min = 4
                max = 20
                default = "6"
                group = "SM_Parameters"
            case ("fsm_ds_min")
                description = "Minimum possible snow layer thickness to allow for a FSM simulation"
                min = 0
                max = 1
                units = "m"
                default = "0.02"
                group = "SM_Parameters"
            case ("fsm_ds_surflay")
                description = "Maximum thickness of surface fine snow layering in an FSM simulation"
                min = 0
                max = 1
                units = "m"
                default = "0.5"
                group = "SM_Parameters"
            case ("fsm_albedo")
                description = "Albedo scheme to use in an FSM simulation"//achar(10)//BLNK_CHR_N// &
                    "(0 = diagnostic, 1 = prognostic, 2 = prognostic, tuned for Switzerland)"
                allocate(values(3))
                values = [0, 1, 2]
                default = "1"
                group = "SM_Parameters"
            case ("fsm_canmod")
                description = "Canopy scheme to use in an FSM simulation"!//achar(10)//BLNK_CHR_N// &
                    !"(1 = random, 2 = maximum-random, 3 = maximum, 4 = exponential, 5 = exponential-random)"
                allocate(values(2))
                values = [0, 1]
                default = "0"
                group = "SM_Parameters"
            case ("fsm_condct")
                description = "Conductivity scheme to use in an FSM simulation"!//achar(10)//BLNK_CHR_N// &
                    !"(1 = random, 2 = maximum-random, 3 = maximum, 4 = exponential, 5 = exponential-random)"
                allocate(values(2))
                values = [0, 1]
                default = "1"
                group = "SM_Parameters"
            case ("fsm_densty")
                description = "Snow compaction scheme to use in an FSM simulation"//achar(10)//BLNK_CHR_N// &
                   "0 = Fixed snow density,"//achar(10)//BLNK_CHR_N// &
                   "1 = Snow compaction with age,"//achar(10)//BLNK_CHR_N// &
                   "2 = Snow compaction by overburden,"//achar(10)//BLNK_CHR_N// &
                   "3 = Snow compaction by overburden, dependent on liquid water content (Crocus B92))"
                allocate(values(4))
                values = [0, 1, 2, 3]
                default = "3"
                group = "SM_Parameters"
            case ("fsm_exchng")
                description = "Surface exchange scheme to use in an FSM simulation"//achar(10)//BLNK_CHR_N// &
                   "(0 = No stability adjustment, 1 = Stability adjustment, 2 = ???)"
                allocate(values(3))
                values = [0, 1, 2]
                default = "1"
                group = "SM_Parameters"
            case ("fsm_hydrol")
                description = "Snow hydraulic scheme to use in an FSM simulation"//achar(10)//BLNK_CHR_N// &
                   "(0 = free-draining snow, 1 = bucket storage, 2 = density-dependent bucket storage"
                allocate(values(3))
                values = [0, 1, 2]
                default = "2"
                group = "SM_Parameters"
            case ("fsm_snfrac")
                description = "Snow fractional area scheme to use in an FSM simulation"//achar(10)//BLNK_CHR_N// &
                    "(0 = ????, 1 = HelbigHS, 2 = HelbigHS0, 3 = Point model, 4 = Original FSM)"
                allocate(values(5))
                values = [0, 1, 2, 3, 4]
                default = "4"
                group = "SM_Parameters"
            case ("fsm_radsbg")
                description = "Radiation subgrid parameterization to use in an FSM simulation"//achar(10)//BLNK_CHR_N// &
                    "(0 = no subgrid parameterization, 1 = subgrid parameterization)"
                allocate(values(2))
                values = [0, 1]
                default = "0"
                group = "SM_Parameters"
            case ("fsm_zoffst")
                description = "Height of input provided to an FSM simulation"//achar(10)//BLNK_CHR_N// &
                    "(0 = height above terrain, 1 = height above terrain + canopy)"
                allocate(values(2))
                values = [0, 1]
                default = "0"
                group = "SM_Parameters"
            case ("fsm_oshdtn")
                description = "OSHD specific tuning options for the FSM2 model"//achar(10)//BLNK_CHR_N// &
                    "of fresh snow albedo, snow roughness lengths and fresh snow density (0 = OFF, 1 = ON)"
                allocate(values(2))
                values = [0, 1]
                default = "0"
                group = "SM_Parameters"
            case ("fsm_alradt")
                description = "Activate tuning of the albedo decay"//achar(10)//BLNK_CHR_N// &
                    "as a function of incoming direct shortwave radiation (0 = OFF, 1 = ON)"
                allocate(values(2))
                values = [0, 1]
                default = "0"
                group = "SM_Parameters"
            case ("fsm_sntran")
                description = "Flag to turn on SnowTran-3D in an FSM simulation"//achar(10)//BLNK_CHR_N// &
                    "(0 = OFF, 1 = ON)"
                allocate(values(2))
                values = [0, 1]
                default = "0"
                group = "SM_Parameters"
            case ("fsm_snslid")
                description = "Flag to turn on SnowSlide in an FSM simulation"//achar(10)//BLNK_CHR_N// &
                    "(0 = OFF, 1 = ON)"
                allocate(values(2))
                values = [0, 1]
                default = "0"
                group = "SM_Parameters"
            case ("fsm_snolay")
                description = "Snow layering scheme to use in an FSM simulation"//achar(10)//BLNK_CHR_N// &
                    "(0 = Original Layering routine, 1 = Density-dependant layering)"
                allocate(values(2))
                values = [0, 1]
                default = "1"
                group = "SM_Parameters"
            case ("fsm_checks")
                description = "Check FSM state variables at each FSM time step"//achar(10)//BLNK_CHR_N// &
                    "(0 = no checks, 1 = some checks, 2 = lotta checks)"
                allocate(values(3))
                values = [0, 1, 2]
                default = "0"
                group = "SM_Parameters"
            case("fsm_hn_on")
                description = "Flag to turn on the HN model in an FSM simulation"//achar(10)//BLNK_CHR_N// &
                    "(.False. = OFF, .True. = ON)"
                default = ".False."
                group = "SM_Parameters"
            case("fsm_for_hn")
                description = "Flag to write 18h states for the hn model"//achar(10)//BLNK_CHR_N// &
                    "(.False. = OFF, .True. = ON)"
                default = ".False."
                group = "SM_Parameters"
            case("fsm_z0pert")
                description = "Flag to turn on the z0 perturbations in an FSM simulation"//achar(10)//BLNK_CHR_N// &
                    "(.False. = OFF, .True. = ON)"
                default = ".False."
                group = "SM_Parameters"
            case("fsm_wcpert")
                description = "Flag to turn on the liquid water capacity perturbations in an FSM simulation"//achar(10)//BLNK_CHR_N// &
                    "(.False. = OFF, .True. = ON)"
                default = ".False."
                group = "SM_Parameters"
            case("fsm_fspert")
                description = "Flag to turn on the fresh snow density perturbations in an FSM simulation"//achar(10)//BLNK_CHR_N// &
                    "(.False. = OFF, .True. = ON)"
                default = ".False."
                group = "SM_Parameters"
            case("fsm_alpert")
                description = "Flag to turn on the albedo perturbations in an FSM simulation"//achar(10)//BLNK_CHR_N// &
                    "(.False. = OFF, .True. = ON)"
                default = ".False."
                group = "SM_Parameters"
            case("fsm_slpert")
                description = "Flag to turn on the settling perturbations in an FSM simulation"//achar(10)//BLNK_CHR_N// &
                    "(.False. = OFF, .True. = ON)"
                default = ".False."
                group = "SM_Parameters"
            ! --------------------------------------
            ! --------------------------------------
            ! Radiation parameters namelist variables
            ! --------------------------------------
            ! --------------------------------------
            case ("update_interval_rrtmg")
                description = "Time interval for updating the radiation. If = 0, update every time step."
                min = 0
                max = 86400
                units = "seconds"
                default = "600.0"
                group = "RAD_Parameters"
            case ("icloud")
                description = "Cloud fraction scheme for radiation"
                allocate(values(2))
                values = [0, 3]
                default = "3"
                group = "RAD_Parameters"
            case ("cldovrlp")
                description = "cloud overlap flag for radiation"//achar(10)//BLNK_CHR_N// &
                    "(1 = random, 2 = maximum-random, 3 = maximum, 4 = exponential, 5 = exponential-random)"
                allocate(values(5))
                values = [1, 2, 3, 4, 5]
                default = "2"
                group = "RAD_Parameters"
            case ("read_ghg")
                description = "read in greenhouse gas data for radiation (T/F)"
                default = ".False."
                group = "RAD_Parameters"
            case ("tzone")
                description = "time zone offset for radiation. Radiation code uses UTC, so if forcing time is not in UTC, give this offset here."
                min = -12
                max = 12
                default = "0"
                group = "RAD_Parameters"
            ! --------------------------------------
            ! --------------------------------------
            ! Wind parameters namelist variables
            ! --------------------------------------
            ! --------------------------------------
            case ("Sx")
                description = "Use Sx-terrain sheltering to modify wind field (T/F)"
                default = ".False."
                group = "Wind"
            case ("thermal")
                description = "Use thermal wind parameterization from Reynolds et al., 2024a to modify wind field (T/F)"
                default = ".False."
                group = "Wind"
            case ("smooth_wind_distance")
                description = "Distance over which to smooth the wind field. If -9999, use default of dx*2"
                min = 0
                max = 10000
                units = "m"
                default = "-9999"
                group = "Wind"
            case ("wind_iterations")
                description = "Number of iterations to use for the O'Brien iterative wind solver (wind=2)"
                min = 0
                max = 1000
                default = "800"
                group = "Wind"
            case ("Sx_dmax")
                description = "Maximum lateral distance over which to calculate the Sx parameter"
                min = 0
                max = 2000
                units = "m"
                default = "300.0"
                group = "Wind"
            case ("Sx_scale_ang")
                description = "Angle over which to scale the Sx parameter. See Reynolds et al., 2023 for details"
                min = 0
                max = 70
                units = "degrees"
                default = "30.0"
                group = "Wind"
            case ("alpha_const")
                description = "Option for setting the alpha parameter in the wind=3 euqtions to a constant"//achar(10)//BLNK_CHR_N// &
                              "(between 0.2 and 2). Default of -1.0 allows for dynamic alpha. For more information, see Reynolds et al., 2023."
                min = 0.2
                max = 2.0
                default = "-1.0"
                group = "Wind"
            case ("TPI_dmax")
                description = "Maximum lateral distance over which to calculate the TPI parameter"
                min = 0
                max = 10000
                units = "m"
                default = "1000.0"
                group = "Wind"
            case ("TPI_scale")
                description = "Scale factor for the TPI parameter. See Reynolds et al., 2023 for details"
                min = 0
                max = 1000
                default = "400.0"
                group = "Wind"
            case ("wind_only")
                description = "Flag for if this is a 'wind-only' run, i.e. only the wind solver should be run"
                default = ".False."
                group = "Wind"
                type = 1
            case ("update_frequency")
                description = "Number of times to update the wind field per input time step."
                min = 0
                max = 50
                default = "1"
                group = "Wind"
            ! --------------------------------------
            ! --------------------------------------
            ! Time parameters namelist variables
            ! --------------------------------------
            ! --------------------------------------
            case ("cfl_reduction_factor")
                description = "Multiplication factor for the CFL time step."//achar(10)//BLNK_CHR_N// &
                    "If < 1, the time step will be reduced by this factor."//achar(10)//BLNK_CHR_N//&
                    "If > 1, the time step will be increased by this factor, only to be used with RK3=.True."
                min = 0.1
                max = 2.0
                default = "0.9"
                group = "Time_Parameters"
            case ("RK3")
                description = "Use the third-order Runge-Kutta time stepping scheme for the advection scheme (T/F)"
                default = ".False."
                group = "Time_Parameters"

            case default
                if (STD_OUT_PE) write(*,*) "---------------------------------------------------------------------------"
                if (STD_OUT_PE) write(*,*) "ERROR: '", trim(name), "' is not a valid namelist variable"
                if (STD_OUT_PE) write(*,*) "---------------------------------------------------------------------------"
                STOP
        end select

        if (.not.(allocated(values))) allocate(values(0))
        if (.not.(allocated(dimensions))) allocate(dimensions(0))

        if (type==1) then
            description = "(-) "//description
        elseif (type==2) then
            description = "(*) "//description
        endif

    end subroutine get_nml_var_metadata

end module namelist_utils
