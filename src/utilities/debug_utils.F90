module debug_module
    use netcdf
    use domain_interface, only  : domain_t
    use string,           only  : str
    use ieee_arithmetic
    use icar_constants,    only : STD_OUT_PE, kVARS
    use iso_fortran_env, only : output_unit
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

        if (domain%var_indx(kVARS%potential_temperature)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%potential_temperature)%v)%data_3d	, name="th",      msg=error_msg, less_than    =100.0, fix=fix_data)
        if (domain%var_indx(kVARS%potential_temperature)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%potential_temperature)%v)%data_3d	, name="th",      msg=error_msg, greater_than =600.0, fix=fix_data)
        if (domain%var_indx(kVARS%water_vapor)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%water_vapor)%v)%data_3d	,           name="qv",      msg=error_msg, less_than    =-1e-10,fix=fix_data)
        if (domain%var_indx(kVARS%cloud_water_mass)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%cloud_water_mass)%v)%data_3d	,      name="cloud",   msg=error_msg, less_than    =-1e-10,fix=fix_data)
        if (domain%var_indx(kVARS%ice_mass)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%ice_mass)%v)%data_3d	,        name="ice",     msg=error_msg, less_than    =-1e-10,fix=fix_data)
        if (domain%var_indx(kVARS%ice_number)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%ice_number)%v)%data_3d	,      name="nice",    msg=error_msg, less_than    =-1e-1, fix=fix_data)
        if (domain%var_indx(kVARS%snow_mass)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%snow_mass)%v)%data_3d	,             name="qsnow",   msg=error_msg, less_than    =-1e-10,fix=fix_data)
        if (domain%var_indx(kVARS%snow_number)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%snow_number)%v)%data_3d	,           name="nsnow",   msg=error_msg, less_than    =-1e-1, fix=fix_data)
        if (domain%var_indx(kVARS%rain_mass)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%rain_mass)%v)%data_3d	,             name="qrain",   msg=error_msg, less_than    =-1e-10,fix=fix_data)
        if (domain%var_indx(kVARS%rain_number)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%rain_number)%v)%data_3d	,           name="nrain",   msg=error_msg, less_than    =-1e-1, fix=fix_data)
        if (domain%var_indx(kVARS%graupel_mass)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%graupel_mass)%v)%data_3d	,          name="qgrau",   msg=error_msg, less_than    =-1e-10,fix=fix_data)
        if (domain%var_indx(kVARS%graupel_number)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%graupel_number)%v)%data_3d	,        name="ngrau",   msg=error_msg, less_than    =-1e-1, fix=fix_data)
        if (domain%var_indx(kVARS%ice1_a)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%ice1_a)%v)%data_3d	,                name="ice1a",   msg=error_msg, less_than    =-1e-10, fix=fix_data)
        if (domain%var_indx(kVARS%ice1_c)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%ice1_c)%v)%data_3d	,                name="ice1c",   msg=error_msg, less_than    =-1e-10, fix=fix_data)
        if (domain%var_indx(kVARS%ice2_mass)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%ice2_mass)%v)%data_3d	,             name="ice2mass",msg=error_msg, less_than    =-1e-10, fix=fix_data)
        if (domain%var_indx(kVARS%ice2_number)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%ice2_number)%v)%data_3d	,           name="ice2num", msg=error_msg, less_than    =-1e-1, fix=fix_data)
        if (domain%var_indx(kVARS%ice2_a)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%ice2_a)%v)%data_3d	,                name="ice2a",   msg=error_msg, less_than    =-1e-10, fix=fix_data)
        if (domain%var_indx(kVARS%ice2_c)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%ice2_c)%v)%data_3d	,                name="ice2c",   msg=error_msg, less_than    =-1e-10, fix=fix_data)
        if (domain%var_indx(kVARS%ice3_mass)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%ice3_mass)%v)%data_3d	,             name="ice3mass",msg=error_msg, less_than    =-1e-10, fix=fix_data)
        if (domain%var_indx(kVARS%ice3_number)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%ice3_number)%v)%data_3d	,           name="ice3num", msg=error_msg, less_than    =-1e-1, fix=fix_data)
        if (domain%var_indx(kVARS%ice3_a)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%ice3_a)%v)%data_3d	,                name="ice3a",   msg=error_msg, less_than    =-1e-10, fix=fix_data)
        if (domain%var_indx(kVARS%ice3_c)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%ice3_c)%v)%data_3d	,                name="ice3c",   msg=error_msg, less_than    =-1e-10, fix=fix_data)
        if (domain%var_indx(kVARS%u)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d	,                     name="u",       msg=error_msg, less_than    =-1e5,  fix=fix_data)
        if (domain%var_indx(kVARS%u)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d	,                     name="u",       msg=error_msg, greater_than =1e5,   fix=fix_data)
        if (domain%var_indx(kVARS%v)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d	,                     name="v",       msg=error_msg, less_than    =-1e5,  fix=fix_data)
        if (domain%var_indx(kVARS%v)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d	,                     name="v",       msg=error_msg, greater_than =1e5,   fix=fix_data)
        if (domain%var_indx(kVARS%w)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%w)%v)%data_3d	,                     name="w",       msg=error_msg, less_than    =-1e5,  fix=fix_data)
        if (domain%var_indx(kVARS%w)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%w)%v)%data_3d	,                     name="w",       msg=error_msg, greater_than =1e5,   fix=fix_data)
        if (domain%var_indx(kVARS%sensible_heat)%v > 0) call check_var2d(domain%vars_2d(domain%var_indx(kVARS%sensible_heat)%v)%data_2d	,       name="hfx",     msg=error_msg) ! check for NaN's only.
        if (domain%var_indx(kVARS%latent_heat)%v > 0) call check_var2d(domain%vars_2d(domain%var_indx(kVARS%latent_heat)%v)%data_2d	,         name="lfx",     msg=error_msg)
        if (domain%var_indx(kVARS%skin_temperature)%v > 0) call check_var2d(domain%vars_2d(domain%var_indx(kVARS%skin_temperature)%v)%data_2d	,    name="tskin",   msg=error_msg)
        if (domain%var_indx(kVARS%roughness_z0)%v > 0) call check_var2d(domain%vars_2d(domain%var_indx(kVARS%roughness_z0)%v)%data_2d	,        name="z0",      msg=error_msg)
        if (domain%var_indx(kVARS%surface_pressure)%v > 0) call check_var2d(domain%vars_2d(domain%var_indx(kVARS%surface_pressure)%v)%data_2d	,    name="psfc",    msg=error_msg)
        ! call check_var2d(domain%ustar,                       name="ustar", msg=error_msg)
        if (domain%var_indx(kVARS%exner)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%exner)%v)%data_3d	,   less_than    =-1e5,              name="pii",     msg=error_msg)
        if (domain%var_indx(kVARS%exner)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%exner)%v)%data_3d	,   greater_than    =5.0,              name="pii",     msg=error_msg)
        if (domain%var_indx(kVARS%pressure_interface)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%pressure_interface)%v)%data_3d	,    name="pi",      msg=error_msg)
        if (domain%var_indx(kVARS%pressure)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%pressure)%v)%data_3d	,              name="p",       msg=error_msg)
        if (domain%var_indx(kVARS%density)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%density)%v)%data_3d	,                     name="density",       msg=error_msg, less_than    =0.0,  fix=fix_data)
        if (domain%var_indx(kVARS%density)%v > 0) call check_var(domain%vars_3d(domain%var_indx(kVARS%density)%v)%data_3d	,                     name="density",       msg=error_msg, greater_than =2.0,   fix=fix_data)

    end subroutine domain_check


    subroutine check_var(var, name, msg, greater_than, less_than, fix)
        implicit none
        real,               intent(inout), allocatable      :: var(:,:,:)
        character(len=*),   intent(in)                      :: name, msg
        real,               intent(in),    optional         :: greater_than, less_than
        logical,            intent(in),    optional         :: fix
        integer :: n
        real :: vmax, vmin
        logical :: printed

        printed = .False.

        if (.not.allocated(var)) then
            return
        endif

        if (any(ieee_is_nan(var))) then
            n = COUNT(ieee_is_nan(var))
            ! ALLOCATE(IsNanIdx(n))
            ! IsNanIdx = PACK( (/(i,i=1,SIZE(var))/), MASK=IsNan(var) )  ! if someone can get this to work it would be nice to have.
            write(*,*) trim(msg)
            write(*,*) trim(name)//" has", n," NaN(s) "
        endif

        if (present(greater_than)) then
            vmax = maxval(var)
            if (vmax > greater_than) then
                write(*,*) trim(msg)
                write(*,*) trim(name)//" is greater than "//trim(str(greater_than))//" : "//trim(str(vmax))

                if (present(fix)) then
                    if (fix) then
                        write(*,*) "Fixing..."
                        where(var > greater_than) var = greater_than
                    endif
                endif
            endif
        endif

        if (present(less_than)) then
            vmin = minval(var)
            if (vmin < less_than) then
                write(*,*) trim(msg)
                write(*,*) trim(name)//" is less than "//trim(str(less_than))//" : "//trim(str(vmin))

                if ((vmin - less_than) < -1e-10) then
                ! we only want to hard stop if there is a significant difference.
                ! Numerical precision can mean that advecting hydrometeors ends up with -1e-30 type values which we can ignore
                block
                    integer :: i,j,k,nx,ny,nz

                    nx = size(var,1)
                    nz = size(var,2)
                    ny = size(var,3)
                    do j=lbound(var,3),ubound(var,3)
                        do k=lbound(var,2),ubound(var,2)
                            do i=lbound(var,1),ubound(var,1)
                                if (var(i,k,j) < less_than) then

                                    if (.not.printed) then
                                        print*, "First Error was in grid cell:", i,k,j, var(i,k,j)
                                        printed = .True.
                                    endif
                                    if (.not.present(fix)) then
                                        error stop
                                    endif
                                endif
                            enddo
                        enddo
                    enddo
                end block
                endif

                if (present(fix)) then
                    if (fix) then
                        write(*,*) "Fixing..."
                        where(var < less_than) var = less_than
                    endif
                endif
            endif
        endif

    end subroutine check_var

    subroutine check_var2d(var, name, msg, greater_than, less_than, fix)
        implicit none
        real,               intent(inout), allocatable      :: var(:,:)
        character(len=*),   intent(in)                      :: name, msg
        real,               intent(in),    optional         :: greater_than, less_than
        logical,            intent(in),    optional         :: fix
        integer :: n, i, j
        real :: vmax, vmin
        logical :: printed

        printed = .False.

        if (.not.allocated(var)) then
            return
        endif

        if (any(ieee_is_nan(var))) then
            n = COUNT(ieee_is_nan(var))
            ! ALLOCATE(IsNanIdx(n))
            ! IsNanIdx = PACK( (/(i,i=1,SIZE(var))/), MASK=IsNan(var) )  ! if someone can get this to work it would be nice to have.
            write(*,*) trim(msg)
            write(*,*) trim(name)//" has", n," NaN(s) "
            if (n < 9) then
                do j = lbound(var,2), ubound(var,2)
                    do i = lbound(var,1), ubound(var,1)
                        if (ieee_is_nan(var(i,j))) print*, "NaN in ",i,j
                    enddo
                enddo
            else
                print*, "Too many NaNs"
                error stop
            endif

        endif
    end subroutine check_var2d

    ! subroutine domain_fix(domain)
    !     implicit none
    !     type(domain_t),  intent(inout)   :: domain
    !
    !     call fix_var(domain%th,     less_than    =100.0)
    !     call fix_var(domain%th,     greater_than =600.0)
    !     call fix_var(domain%qv,     less_than    =0.0)
    !     call fix_var(domain%cloud,  less_than    =0.0)
    !     call fix_var(domain%ice,    less_than    =0.0)
    !     call fix_var(domain%nice,   less_than    =0.0)
    !     call fix_var(domain%qsnow,  less_than    =0.0)
    !     call fix_var(domain%nsnow,  less_than    =0.0)
    !     call fix_var(domain%qrain,  less_than    =0.0)
    !     call fix_var(domain%nrain,  less_than    =0.0)
    !     call fix_var(domain%qgrau,  less_than    =0.0)
    !     call fix_var(domain%ngraupel,less_than   =0.0)
    !
    ! end subroutine domain_fix

    subroutine fix_var(var, greater_than, less_than)
        implicit none
        real,               intent(inout), dimension(:,:,:), allocatable :: var
        real,               intent(in),    optional         :: greater_than, less_than

        if (.not.allocated(var)) then
            return
        endif

        if (present(greater_than)) then
            where(var > greater_than) var = greater_than
        endif

        if (present(less_than)) then
            where(var < less_than) var = less_than
        endif

    end subroutine fix_var

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
