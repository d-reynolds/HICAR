module time_io

    use time_object,        only : Time_type
    use time_delta_object,  only : time_delta_t
    use string,             only : get_integer, as_string
    use io_routines,        only : io_read, io_read_attribute, io_dimension_is_present, check
    use iso_fortran_env,    only: real64, real64
    use icar_constants,     only: kMAX_STRING_LENGTH, kMAX_NAME_LENGTH, STD_OUT_PE, kOUTPUT_FMT
    use netcdf

    implicit none

contains

    function find_timestep_in_file(filename, time_var, time, time_at_step, precision, forward, error) result(step)
        implicit none
        character(len=*),  intent(in) :: filename
        character(len=*),  intent(in) :: time_var
        type(Time_type),   intent(in) :: time
        type(Time_type),   intent(inout), optional :: time_at_step
        type(time_delta_t),intent(in),    optional :: precision
        logical,           intent(in),    optional :: forward
        integer,           intent(inout), optional :: error
        integer :: step

        type(Time_type), allocatable :: times_in_file(:)
        type(time_delta_t) :: max_dt
        integer :: i,n
        logical :: found

        call max_dt%set(seconds=1.0)
        if (present(precision)) max_dt = precision
        if (present(error)) error=0

        ! read the times for all timesteps in the specified file
        call read_times(filename, time_var, times_in_file)

        step = -1
        found= .False.
        n    = size(times_in_file)
        ! loop through times looking for a time that matches the input time to within
        ! a specified maximum delta t
        do i = 1, n
            if (.not.found) then
                if (times_in_file(i)%equals(time, precision=max_dt)) then
                    step = i
                    found=.True.
                elseif (i > 1) then
                    if ((times_in_file(i) > time) .and. (times_in_file(i-1) < time) .and. present(forward)) then
                        if (forward) then
                            step = i
                        else
                            step = i-1
                        endif
                        found=.True.
                    endif
                endif
            endif
        enddo

        if (step < 1) then
            if (present(error)) then
                error = 1
            else
                 write(*,*) "ERROR: Unable to find requested date in file."
                 write(*,*) "Filename: ",trim(filename)
                 write(*,*) "  time  : ",trim(as_string(time))
                write(*,*) "First time in file : ", trim(as_string(times_in_file(1)))
                write(*,*) " Last time in file : ", trim(as_string(times_in_file(n)))
                stop "Unable to find date in file"
            endif
        endif

        if (present(time_at_step)) time_at_step = times_in_file(step)

        deallocate(times_in_file)

    end function find_timestep_in_file


    function time_gain_from_units(units) result(gain)
        implicit none
        character(len=*), intent(in) :: units
        real(real64) :: gain

        if ((units(1:4)=="days").or.(units(1:4)=="Days")) then
            gain = 1.0D0
        else if ((units(1:4)=="hour").or.(units(1:4)=="Hour")) then
            gain = 24.0D0
        else if ((units(1:3)=="min").or.(units(1:3)=="Min")) then
            gain = 1440.0D0
        else if ((units(1:3)=="sec").or.(units(1:3)=="Sec")) then
            gain = 86400.0D0
        else if ((units(1:11)=="nanoseconds").or.(units(1:11)=="Nanoseconds")) then
            gain = 86400000000000.0D0
        else
            if (STD_OUT_PE) write(*,*) 'unknown units for input time: ', trim(units)
            stop "Error: unknown units"
        endif

    end function time_gain_from_units

    function year_from_units(units) result(year)
        implicit none
        character(len=*), intent(in) :: units
        integer :: year

        integer :: since_loc, year_loc

        since_loc = index(units,"since")

        year_loc = index(units(since_loc:)," ")
        year_loc = year_loc+since_loc

        year = get_integer(units(year_loc:year_loc+3))

    end function year_from_units

    function month_from_units(units) result(month)
        implicit none
        character(len=*), intent(in) :: units
        integer :: month

        integer :: since_loc, month_loc

        since_loc = index(units,"since")

        month_loc = index(units(since_loc:)," ") + 5
        month_loc = month_loc + since_loc

        month = get_integer(units(month_loc:month_loc+1))

    end function month_from_units

    function day_from_units(units) result(day)
        implicit none
        character(len=*), intent(in) :: units
        integer :: day

        integer :: since_loc, day_loc

        since_loc = index(units,"since")

        day_loc = index(units(since_loc:)," ") + 8
        day_loc = day_loc + since_loc

        day = get_integer(units(day_loc:day_loc+1))

    end function day_from_units

    function hour_from_units(units, error) result(hour)
        implicit none
        character(len=*), intent(in) :: units
        integer, intent(out), optional :: error
        integer :: hour

        integer :: since_loc, hour_loc
        if (present(error)) error = 0

        since_loc = index(units,"since")

        hour_loc = index(units(since_loc:)," ") + 11
        hour_loc = hour_loc + since_loc

        ! default return value if hours can't be read from the units attribute (e.g. they aren't present)
        hour = 0

        if( hour_loc+1 <= len(units) ) then
           if (trim(units(hour_loc:hour_loc+1)) /= "") then
               hour = get_integer(units(hour_loc:hour_loc+1))
           endif
        else
           if (present(error)) error = 1
        endif

    end function hour_from_units

    function find_timestep_in_filelist(filelist, time_var, time, filename, forward, error) result(step)
        implicit none
        character(len=*),  intent(in) :: filelist(:)
        character(len=*),  intent(in) :: time_var
        type(Time_type),   intent(in) :: time
        character(len=*),  intent(inout), optional :: filename
        logical,           intent(in), optional :: forward
        integer,           intent(inout), optional :: error
        integer :: step, err

        type(Time_type), allocatable :: times_in_file_up(:), times_in_file_down(:)
        type(time_delta_t) :: max_dt
        type(Time_type) :: tmp_time
        integer :: i, n_t_down, n_t_up, nup, ndown, nup_o, ndown_o, ierr
        logical :: found, fileExists, limit_loop

        call max_dt%set(seconds=1.0)
        if (present(error)) error=0
        
        step = -1
        found = .False.
        
        nup = size(filelist)
        ndown=1
        ierr = 1
        
        !Find last existing file in file list (necesarray if doing runs concurrently)
        do while (ierr==1)
            !See if nup exists
            inquire(file=filelist(nup), exist=fileExists)

            !If nup does exist, set ierr to 0
            if (fileExists) then
                ierr = 0
            !else if nup does not exist, decrement nup
            else
                nup = nup-1
            endif
            !If we reach beginning and nothing has existed, stop and error
            if (nup < 1) stop "No files in file list exist"
        enddo

        ndown_o = ndown
        nup_o = nup

        call read_times(filelist(ndown), time_var, times_in_file_down)
        call read_times(filelist(nup), time_var, times_in_file_up)

        !Quick check that time is bounded by min and max existing files in list
        if ( (times_in_file_down(1) > time .and. .not.(times_in_file_down(1)%equals(time, precision=max_dt))) .or. &
             (times_in_file_up(size(times_in_file_up)) < time .and. .not.(times_in_file_up(size(times_in_file_up))%equals(time, precision=max_dt))) ) then
            if (present(error)) then
                error = 1
                found = .True.
            else
                write(*,*) "ERROR: Requested date lays outside of filelist"
                write(*,*) "First filename:           ",trim(filelist(1))
                write(*,*) "Last (existing) filename: ",trim(filelist(nup))
                write(*,*) "  last attempted step  : ",step
                write(*,*) "  time  : ",trim(as_string(time))
                stop "Unable to find date in file"
            endif
        endif
        limit_loop = .False.

        do while(.not.(found))
            ! read the times for all timesteps in the specified file

            if (present(forward)) then
                step = find_timestep_in_file(filelist(ndown), time_var, time, precision=max_dt, forward=forward, error=err)
            else
                step = find_timestep_in_file(filelist(ndown), time_var, time, precision=max_dt, error=err)
            endif

            if (step > 0) then
                found = .True.
                filename = filelist(ndown)
            else
                if (present(forward)) then
                    step = find_timestep_in_file(filelist(nup), time_var, time, precision=max_dt, forward=forward, error=err)
                else
                    step = find_timestep_in_file(filelist(nup), time_var, time, precision=max_dt, error=err)
                endif

                if (step > 0) then
                    found = .True.
                    filename = filelist(nup)
                endif    
            endif

            if (.not.found) then
                !Pincer-move to tighten search
                !Time must be bounded, save these positions
                ndown_o = ndown
                nup_o = nup
                !Advance ndown by half of distance between nup and ndown
                ndown = ndown + nint((nup-ndown)/2.0)

                ! read the times for all timesteps in the specified file
                if (.not.(nup==nup_o)) call read_times(filelist(nup), time_var, times_in_file_up)
                if (.not.(ndown==ndown_o)) call read_times(filelist(ndown), time_var, times_in_file_down)

                ! See if we overshot (since we advance up)
                if(times_in_file_down(1) > time .and. times_in_file_up(1) > time) then
                    !set upper bound to former lower bound
                    nup = ndown
                    times_in_file_up = times_in_file_down
                    nup_o = nup
                    !Set lower bound to last bounding lower bound
                    ndown = ndown_o
                    ndown_o = 0
                endif

                !Handle the limiting case
                if (ndown==nup) then
                    ndown=max(nup-1,1)
                    if (limit_loop) then
                        if (present(forward)) then
                            write(*,*) "ERROR: Unable to find requested date in filelist."
                            write(*,*) "We have collapsed onto two files, but no time was found."
                            write(*,*) "This, despite a non-exact time search. Gnarly bug."
                            write(*,*) "First filename: ",trim(filelist(1))
                            write(*,*) "  time  : ",trim(as_string(time))
                            write(*,*) "  last attempted step  : ",step
                            stop "Unable to find date in file"
                        else
                            write(*,*) "WARNING: Unable to find requested date in filelist."
                            write(*,*) "We have collapsed onto two files, but no time was found."
                            write(*,*) "An exact time search was requested. Perhaps you want a non-exact search?"
                            write(*,*) "First filename: ",trim(filelist(1))
                            write(*,*) "  time  : ",trim(as_string(time))
                            write(*,*) "  last attempted step  : ",step
                            stop "Unable to find date in file"
                        endif
                    endif    
                    limit_loop = .True.
                endif    
                if (ndown==(nup-1)) then

                    call read_times(filelist(nup), time_var, times_in_file_up)
                    call read_times(filelist(ndown), time_var, times_in_file_down)

                    tmp_time = times_in_file_down(size(times_in_file_down))
                    if (tmp_time%equals(time, precision=max_dt)) then
                        step = size(times_in_file_down)
                        filename = filelist(ndown)
                        found = .True.
                        cycle
                    endif

                    tmp_time = times_in_file_up(1)
                    if (tmp_time%equals(time, precision=max_dt)) then
                        step = 1
                        filename = filelist(nup)
                        found = .True.
                        cycle
                    endif

                    if (times_in_file_down(size(times_in_file_down)) < time .and. times_in_file_up(1) > time) then
                        !Time is between the two
                        if (present(forward)) then
                            if (forward) then
                                step = 1
                                filename = filelist(nup)
                                found = .True.
                            else
                                step = size(times_in_file_down)
                                filename = filelist(ndown)
                                found = .True.
                            endif
                        else
                            write(*,*) "WARNING: Unable to find requested date in filelist."
                            write(*,*) "Time appears to lay between two files."
                            write(*,*) "An exact time search was requested. Perhaps you want a non-exact search?"
                            write(*,*) "'Lower' filename: ",trim(filelist(ndown))
                            write(*,*) "'Upper' filename: ",trim(filelist(nup))
                            write(*,*) "Last time in 'Upper' file: ",trim(as_string(times_in_file_up(size(times_in_file_up))))
                            write(*,*) "First time in 'Lower' file: ",trim(as_string(times_in_file_down(1)))
                            write(*,*) "  searched for time  : ",trim(as_string(time))
                            write(*,*) "  last attempted step  : ",step
                            stop "Unable to find date in file"
                        endif
                    endif
                endif
            endif
        enddo
        
        deallocate(times_in_file_down)
        deallocate(times_in_file_up)

    end function find_timestep_in_filelist


    subroutine read_times(filename, varname, times, timezone_offset, curstep)
        implicit none
        character(len=*),   intent(in) :: filename, varname
        type(Time_type),    intent(inout), allocatable, dimension(:) :: times
        real(real64),      intent(in), optional :: timezone_offset
        integer,            intent(in), optional :: curstep

        real(real64),  allocatable, dimension(:) :: temp_times_64
        real(real64), allocatable, dimension(:) :: temp_times_128
        integer :: time_idx, error
        integer :: start_year, start_month, start_day, start_hour
        character(len=kMAX_STRING_LENGTH) :: calendar, units
        real(real64) :: calendar_gain

        ! first read the time variable (presumebly a 1D real(real64) array)
        if (present(curstep)) then
            call io_read(trim(filename), trim(varname), temp_times_64, curstep=curstep)
        else
            call io_read(trim(filename), trim(varname), temp_times_64)
        endif

        ! attempt to read the calendar attribute from the time variable
        call io_read_attribute(trim(filename),"calendar", calendar, var_name=trim(varname), error=error)
        ! if time attribute it not present, set calendar to one specified in the config file
        if (error/=0) then
            if (STD_OUT_PE) write(*,*) "WARNING: assuming standard/gregorian calendar for file "//trim(filename)
            calendar = "standard"
        endif

        ! attempt to read the units for this time variable
        call io_read_attribute(trim(filename), "units", units, var_name=trim(varname), error=error)

        ! if units attribute was present, then read information from it.
        if (error==0) then
            start_year    = year_from_units(units)
            start_month   = month_from_units(units)
            start_day     = day_from_units(units)
            start_hour    = hour_from_units(units)
            ! based off of the string "Days since" (or "seconds" or...)
            calendar_gain = time_gain_from_units(units)
        else

            stop "Time variable does not have units attribute"
        endif

        ! converts the input units to "days since ..."
        ! in case it is in units of e.g. "hours since" or "seconds since"
        allocate(temp_times_128(size(temp_times_64)))
        temp_times_128 = temp_times_64 / calendar_gain

        if (present(timezone_offset)) then
            temp_times_128 = temp_times_128 + timezone_offset / 24.0
        endif

        if (allocated(times)) deallocate(times)
        allocate(times(size(temp_times_128)))

        do time_idx = 1, size(temp_times_128,1)
            call times(time_idx)%init(calendar, start_year, start_month, start_day, start_hour)
            call times(time_idx)%set(days=temp_times_128(time_idx))
        end do

        deallocate(temp_times_64, temp_times_128)

    end subroutine read_times

    function get_output_time(time, units, round_seconds) result(output_time)
        implicit none
        type(Time_type),  intent(in) :: time
        character(len=*), intent(in), optional :: units
        logical,          intent(in), optional :: round_seconds

        type(Time_type) :: output_time
        type(time_delta_t) :: half_minute

        integer :: year, month, day, hour, minute, seconds
        integer :: year0, month0, day0, hour0, minute0, seconds0
        character(len=kMAX_NAME_LENGTH) :: use_units

        if (present(units)) then
            use_units = units
        else
            write(use_units,kOUTPUT_FMT) time%year_zero,time%month_zero,time%day_zero,time%hour_zero
        endif

        call time%date(year, month, day, hour, minute, seconds)
        year0 = year_from_units(use_units)
        month0 = month_from_units(use_units)
        day0 = day_from_units(use_units)
        hour0 = hour_from_units(use_units)
        minute0 = 0 ! minute_from_units(use_units)
        seconds0 = 0 ! seconds_from_units(use_units)

        call output_time%init(time%get_calendar(), year0, month0, day0, hour0)

        if (present(round_seconds)) then
            if (round_seconds) then
                if (seconds > 30) then
                    call output_time%set(year, month, day, hour, minute, seconds)
                    call half_minute%set(seconds=30)
                    call output_time%set(output_time%mjd() + half_minute%days())

                    ! get a new date after adding 30 seconds
                    call output_time%date(year, month, day, hour, minute, seconds)
                    ! use that date after setting seconds to 0 this rounds the old date up by up to 30s
                    call output_time%set(year, month, day, hour, minute, 0)
                else
                    call output_time%set(year, month, day, hour, minute, 0)
                endif
            else
                call output_time%set(year, month, day, hour, minute, seconds)
            endif
        else
            call output_time%set(year, month, day, hour, minute, seconds)
        endif

    end function get_output_time


    function var_has_time_dim(filename, varname, time_var) result(has_time)
        implicit none
        character(len=*), intent(in) :: filename, varname, time_var
        logical :: has_time

        integer :: ncid, varid, ndims, dimids(10), i
        character(len=kMAX_NAME_LENGTH) :: time_dim_name

        call check( nf90_open(trim(filename), nf90_nowrite, ncid), "Opening file "//trim(filename))
        call check( nf90_inq_varid(ncid, trim(varname), varid), "Getting varid for "//trim(varname))

        call check( nf90_open(filename, IOR(NF90_NOWRITE,NF90_NETCDF4), ncid), " Opening file "//trim(filename))
        !get the dimension name for the single dimension on the variable this%time_var
        call check( nf90_inq_varid(ncid, trim(time_var), varid), "Getting varid for "//trim(time_var))
        call check( nf90_inquire_variable(ncid, varid, ndims=ndims, dimids=dimids), " Getting dim length for "//trim(time_var))
        if (ndims==1) then
            call check(nf90_inquire_dimension(ncid, dimids(1), name=time_dim_name), " Getting time dim name")
        elseif (ndims==0) then
            ! scalar variable, use name of time variable as dimension name
            time_dim_name = trim(time_var)
        else
            write(*,*) "ERROR: time variable "//trim(time_var)//" has more than one dimension"
            stop "Error in var_has_time_dim"
        endif

        call check( nf90_close(ncid), "Closing file "//trim(filename))

        has_time = io_dimension_is_present(trim(filename), varname, trim(time_dim_name))

    end function var_has_time_dim
end module time_io
