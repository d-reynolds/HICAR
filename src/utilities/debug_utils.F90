module debug_module
    use netcdf
    use domain_interface, only  : domain_t
    use string,           only  : str
    use ieee_arithmetic
    use icar_constants,    only : STD_OUT_PE, kVARS, kUNSET_REAL, kMAX_NAME_LENGTH
    use iso_fortran_env, only : output_unit
    use io_routines,            only : io_write
    use mpi_utils_module, only : get_mpi_global_rank
    use output_metadata, only : get_varname
    use variable_interface, only : variable_t
    implicit none

contains

    subroutine domain_check(domain, error_msg, fix)
        implicit none
        type(domain_t),     intent(inout)   :: domain
        character(len=*),   intent(in)      :: error_msg
        logical,            intent(in), optional :: fix
        logical :: fix_data

        fix_data = .False.
        if (present(fix)) fix_data = fix

        if (domain%var_indx(kVARS%potential_temperature)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%potential_temperature)%v), msg=error_msg, fix=fix_data)
        if (domain%var_indx(kVARS%water_vapor)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%water_vapor)%v), msg=error_msg, fix=fix_data)
        if (domain%var_indx(kVARS%cloud_water_mass)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%cloud_water_mass)%v), msg=error_msg, fix=fix_data)
        if (domain%var_indx(kVARS%ice_mass)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%ice_mass)%v), msg=error_msg, fix=fix_data)
        if (domain%var_indx(kVARS%ice_number)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%ice_number)%v), msg=error_msg, fix=fix_data)
        if (domain%var_indx(kVARS%snow_mass)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%snow_mass)%v), msg=error_msg, fix=fix_data)
        if (domain%var_indx(kVARS%snow_number)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%snow_number)%v), msg=error_msg, fix=fix_data)
        if (domain%var_indx(kVARS%rain_mass)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%rain_mass)%v), msg=error_msg, fix=fix_data)
        if (domain%var_indx(kVARS%rain_number)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%rain_number)%v), msg=error_msg, fix=fix_data)
        if (domain%var_indx(kVARS%graupel_mass)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%graupel_mass)%v), msg=error_msg, fix=fix_data)
        if (domain%var_indx(kVARS%graupel_number)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%graupel_number)%v), msg=error_msg, fix=fix_data)
        if (domain%var_indx(kVARS%ice1_a)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%ice1_a)%v), msg=error_msg, fix=fix_data)
        if (domain%var_indx(kVARS%ice1_c)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%ice1_c)%v), msg=error_msg, fix=fix_data)
        if (domain%var_indx(kVARS%ice2_mass)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%ice2_mass)%v), msg=error_msg, fix=fix_data)
        if (domain%var_indx(kVARS%ice2_number)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%ice2_number)%v), msg=error_msg, fix=fix_data)
        if (domain%var_indx(kVARS%ice2_a)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%ice2_a)%v), msg=error_msg, fix=fix_data)
        if (domain%var_indx(kVARS%ice2_c)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%ice2_c)%v), msg=error_msg, fix=fix_data)
        if (domain%var_indx(kVARS%ice3_mass)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%ice3_mass)%v), msg=error_msg, fix=fix_data)
        if (domain%var_indx(kVARS%ice3_number)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%ice3_number)%v), msg=error_msg, fix=fix_data)
        if (domain%var_indx(kVARS%ice3_a)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%ice3_a)%v), msg=error_msg, fix=fix_data)
        if (domain%var_indx(kVARS%ice3_c)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%ice3_c)%v), msg=error_msg, fix=fix_data)
        if (domain%var_indx(kVARS%sensible_heat)%v > 0) call check_var(domain%vars_2d(domain%var_indx(kVARS%sensible_heat)%v), msg=error_msg) ! check for NaN's only.
        if (domain%var_indx(kVARS%latent_heat)%v > 0) call check_var(domain%vars_2d(domain%var_indx(kVARS%latent_heat)%v), msg=error_msg)
        if (domain%var_indx(kVARS%skin_temperature)%v > 0) call check_var(domain%vars_2d(domain%var_indx(kVARS%skin_temperature)%v), msg=error_msg)
        if (domain%var_indx(kVARS%roughness_z0)%v > 0) call check_var(domain%vars_2d(domain%var_indx(kVARS%roughness_z0)%v), msg=error_msg)
        if (domain%var_indx(kVARS%surface_pressure)%v > 0) call check_var(domain%vars_2d(domain%var_indx(kVARS%surface_pressure)%v), msg=error_msg)
        if (domain%var_indx(kVARS%exner)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%exner)%v), msg=error_msg)
        if (domain%var_indx(kVARS%pressure_interface)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%pressure_interface)%v), msg=error_msg)
        if (domain%var_indx(kVARS%pressure)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%pressure)%v), msg=error_msg)
        if (domain%var_indx(kVARS%density)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%density)%v), msg=error_msg, fix=fix_data)
        call domain_check_winds(domain, error_msg)

    end subroutine domain_check

    subroutine domain_check_winds(domain, error_msg)
        implicit none
        type(domain_t), intent(inout) :: domain
        character(len=*), intent(in) :: error_msg

        associate(u => domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d, &
                  v => domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d, &
                  w => domain%vars_3d(domain%var_indx(kVARS%w)%v)%data_3d)
        !$acc update host(u, v, w)
        end associate

        if (domain%var_indx(kVARS%u)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%u)%v), msg=error_msg)
        if (domain%var_indx(kVARS%u)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%u)%v), msg=error_msg)
        if (domain%var_indx(kVARS%v)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%v)%v), msg=error_msg)
        if (domain%var_indx(kVARS%v)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%v)%v), msg=error_msg)
        if (domain%var_indx(kVARS%w)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%w)%v), msg=error_msg)
        if (domain%var_indx(kVARS%w)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%w)%v), msg=error_msg)


    end subroutine domain_check_winds

    subroutine check_var(var, msg, fix)
        implicit none
        type(variable_t),                intent(inout)      :: var
        character(len=*),   intent(in)                      :: msg
        logical,            intent(in),    optional         :: fix
        integer :: n, i,j,k
        real :: vmax, vmin, greater_than, less_than
        real, allocatable :: var_3d(:,:,:)
        real, allocatable :: var_2d(:,:)
        logical :: err_flag
        integer :: PE_rank_global
        character(len=256) :: file_name
        character(len=kMAX_NAME_LENGTH) :: name

        err_flag = .False.

        !get name for given id
        name = trim(get_varname(var%id))

        if (var%three_d) then
            var_3d = var%data_3d
            vmax = maxval(var_3d)
            vmin = minval(var_3d)
            n = COUNT(ieee_is_nan(var_3d))
        else if (var%two_d) then
            var_2d = var%data_2d
            vmax = maxval(var_2d)
            vmin = minval(var_2d)
            n = COUNT(ieee_is_nan(var_2d))
        else
            write(*,*) "check_var only works for 2D or 3D variables."
            error stop
        endif

        less_than = var%minval
        greater_than = var%maxval

        if (less_than == kUNSET_REAL) then
            write(*,*) "WARNING: Variable ", trim(name), " minval not set in metadata. Please provide a value in the metadata definition."
            less_than = -1.0e30
        endif
        if (greater_than == kUNSET_REAL) then
            write(*,*) "WARNING: Variable ", trim(name), " maxval not set in metadata. Please provide a value in the metadata definition."
            greater_than = 1.0e30
        endif

        !get rank of this MPI process on the global communicator
        PE_rank_global = get_mpi_global_rank()

        if (n > 0) then
            err_flag = .True.
            ! ALLOCATE(IsNanIdx(n))
            ! IsNanIdx = PACK( (/(i,i=1,SIZE(var))/), MASK=IsNan(var) )  ! if someone can get this to work it would be nice to have.
            write(*,*) trim(msg)
            write(*,*) trim(name)//" has", n," NaN(s) "
        endif

        if (vmax > greater_than) then
            err_flag = .True.

            write(*,*) trim(msg)
            write(*,*) trim(name)//" is greater than "//trim(str(greater_than))//" : "//trim(str(vmax))

            !print index of first occurance and break out of loop
            if (var%three_d) then
                outer: do j=lbound(var_3d,3),ubound(var_3d,3)
                    do k=lbound(var_3d,2),ubound(var_3d,2)
                        do i=lbound(var_3d,1),ubound(var_3d,1)
                            if (var_3d(i,k,j) > greater_than) then
                                print*, "First Error was in grid cell:", i,k,j, var_3d(i,k,j), " on PE: ", PE_rank_global
                                exit outer
                            endif
                        enddo
                    enddo
                enddo outer
            endif
        endif

        if (vmin < less_than) then
            err_flag = .True.
            write(*,*) trim(msg)
            write(*,*) trim(name)//" is less than "//trim(str(less_than))//" : "//trim(str(vmin))

            !print index of first occurance and break out of loop
            if (var%three_d) then
                outer2: do j=lbound(var_3d,3),ubound(var_3d,3)
                    do k=lbound(var_3d,2),ubound(var_3d,2)
                        do i=lbound(var_3d,1),ubound(var_3d,1)
                            if (var_3d(i,k,j) < less_than) then
                                print*, "First Error was in grid cell:", i,k,j, var_3d(i,k,j), " on PE: ", PE_rank_global
                                exit outer2
                            endif
                        enddo
                    enddo
                enddo outer2
            endif
        endif

        if (err_flag) then

            !write out contents of var to a .nc file for debugging. 
            ! File name is the the variable name with a _GLOBAL_PE.nc suffix
            ! Save the file to the debug/ output directory, creating it if necessary.

            ! Create the debug/ output directory if it doesn't exist
            call execute_command_line('mkdir -p debug', wait=.true.)
            file_name = 'debug/'//trim(name)//"_"//trim(str(PE_rank_global))//".nc"
            write(*,*) "Writing debug output to ", trim(file_name)
            call flush(output_unit)
            if (var%three_d) call io_write(trim(file_name), trim(name),var_3d)
            if (var%two_d)   call io_write(trim(file_name), trim(name),var_2d)

            if (present(fix)) then
                if (fix) then
                    write(*,*) "Fixing..."
                    if (var%three_d) then
                        where(var%data_3d > greater_than) var%data_3d = greater_than
                        where(var%data_3d < less_than)    var%data_3d = less_than
                    elseif (var%two_d) then
                        where(var%data_2d > greater_than) var%data_2d = greater_than
                        where(var%data_2d < less_than)    var%data_2d = less_than
                    endif
                endif
            else
                error stop
            endif
        endif

    end subroutine check_var

    !>------------------------------------------------------------
    !! Simple error handling for common netcdf file errors
    !!
    !! If status does not equal nf90_noerr, then print an error message and STOP
    !! the entire program.
    !!
    !! @param   status  integer return code from nc_* routines
    !! @param   extra   OPTIONAL string with extra context to print in case of an error
    !!
    !!------------------------------------------------------------
    subroutine check_ncdf(status,extra)
        implicit none
        integer, intent ( in) :: status
        character(len=*), optional, intent(in) :: extra

        ! check for errors
        if(status /= nf90_noerr) then
            ! print a useful message
            write(*,*) trim(nf90_strerror(status))
            flush(output_unit)
            if(present(extra)) then
                ! print any optionally provided context
                write(*,*) trim(extra)
                flush(output_unit)
            endif
            ! STOP the program execution
            stop "Stopped"
        end if
    end subroutine check_ncdf

end module debug_module
