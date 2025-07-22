!>------------------------------------------------------------
!!  Various functions to convert a number to a string and a string
!!  to a number
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!------------------------------------------------------------
module string

    use iso_fortran_env, only : real32, real64
    use icar_constants, only : kMAX_STRING_LENGTH
    use time_object, only : Time_type
    use timer_interface, only : timer_t
    use time_delta_object, only : time_delta_t
    implicit none
    !>-----------------------------
    !!  Generic interface to various types of string conversion functions
    !!
    !!-----------------------------
    interface str
        module procedure str_d
        module procedure str_r
        module procedure str_i
    end interface

    interface as_string
        module procedure :: as_string_dt
        module procedure :: as_string_time
        module procedure :: as_string_timer
    end interface

contains
    !>------------------------------
    !! Convert a string to a double precision real number
    !!
    !!------------------------------
    elemental function get_double(str_in)
        implicit none
        character(len=*), intent(in) :: str_in      ! input string
        double precision :: get_double              ! double precision return value

        read(str_in,*) get_double

    end function get_double

    !>------------------------------
    !! Convert a string to a single precision real number
    !!
    !!------------------------------
    elemental function get_real(str_in)
        implicit none
        character(len=*), intent(in) :: str_in      ! input string
        real :: get_real                            ! return real value

        read(str_in,*) get_real

    end function get_real

    !>------------------------------
    !! Convert a string to a single precision integer
    !!
    !!------------------------------
    elemental function get_integer(str_in)
        implicit none
        character(len=*), intent(in) :: str_in      ! input string
        integer :: get_integer                      ! return integer value

        read(str_in,*) get_integer

    end function get_integer


    !>------------------------------
    !! Convert a double precision real number to a string
    !!
    !!  Optionally specify a format string
    !!
    !!------------------------------
    elemental function str_d(value,fmt) result(output_string)
        implicit none
        double precision, intent(in) :: value                       ! double precision value to be converted
        character(len=*), intent(in), optional :: fmt               ! optional format string for conversion
        character(len=kMAX_STRING_LENGTH) :: output_string ! return value

        if (present(fmt)) then
            write(output_string,fmt) value
        else
            write(output_string,*) value
        endif

        output_string=trim(adjustl(output_string))
    end function str_d

    !>------------------------------
    !! Convert a single precision real number to a string
    !!
    !!  Optionally specify a format string
    !!
    !!------------------------------
    elemental function str_r(value,fmt) result(output_string)
        implicit none
        real, intent(in) :: value                                   ! single precision value to be converted
        character(len=*), intent(in), optional :: fmt               ! optional format string for conversion
        character(len=kMAX_STRING_LENGTH) :: output_string ! return value

        if (present(fmt)) then
            write(output_string,fmt) value
        else
            write(output_string,*) value
        endif

        output_string=trim(adjustl(output_string))
    end function str_r

    !>------------------------------
    !! Convert a single precision integer to a string
    !!
    !!  Optionally specify a format string
    !!
    !!------------------------------
    elemental function str_i(value,fmt) result(output_string)
        implicit none
        integer, intent(in) :: value                                ! integer value to be converted
        character(len=*), intent(in), optional :: fmt               ! optional format string for conversion
        character(len=kMAX_STRING_LENGTH) :: output_string ! return value

        if (present(fmt)) then
            write(output_string,fmt) value
        else
            write(output_string,*) value
        endif

        output_string=trim(adjustl(output_string))
    end function str_i

    !>------------------------------
    !! Split a string into tokens according to a delimiter
    !!
    !!------------------------------
    pure function split_str(tmp_str, delimiter) result(tokens)
        implicit none
        character(len=*), intent(in)  :: tmp_str
        character(len=1), intent(in) :: delimiter
        character(len=kMAX_STRING_LENGTH), allocatable :: tokens(:)

        integer :: i, n, start, end
        character(len=kMAX_STRING_LENGTH) :: str_in

        str_in = trim(tmp_str)
        str_in = trim(str_in)

        n = 1
        do i = 1, len(str_in)
            if (str_in(i:i) == delimiter) then
                n = n + 1
            endif
        end do

        allocate(tokens(n))

        n = 0
        start = 1
        do i = 1, len(str_in)
            if (str_in(i:i) == delimiter) then
                end = i - 1
                n = n + 1
                tokens(n) = str_in(start:end)
                start = i + 1
            endif
        end do

        n = n + 1
        tokens(n) = str_in(start:)
    end function split_str

    !>------------------------------
    !! Convert a string to lower case
    !!
    !!------------------------------
    pure function to_lower(str_in) result(str_out)
        implicit none
        character(len=*), intent(in) :: str_in
        character(len=len(str_in)) :: str_out

        integer :: i

        str_out = trim(str_in)
        do i = 1, len(trim(str_in))
            if (ichar(str_in(i:i)) >= iachar('A') .and. ichar(str_in(i:i)) <= iachar('Z')) then
                str_out(i:i) = char(ichar(str_in(i:i)) + iachar('a') - iachar('A'))
            endif
        end do
        str_out = trim(str_out)

    end function to_lower

        !>------------------------------
    !! Convert a double precision real number to a string
    !!
    !!  Optionally specify a format string
    !!
    !!------------------------------
    function as_string_dt(time) result(pretty_string)
        implicit none
        type(time_delta_t), intent(in) :: time                       ! double precision value to be converted
        character(len=kMAX_STRING_LENGTH)  :: pretty_string

        if (abs(time%seconds()) < 1) then
            write(pretty_string,"(F8.4,A)") time%seconds(), " seconds"
        elseif (abs(time%seconds()) <= 60) then
            write(pretty_string,"(F6.2,A)") time%seconds(), " seconds"
        elseif (abs(time%minutes()) <= 60) then
            write(pretty_string,"(F6.2,A)") time%minutes(), " minutes"
        elseif (abs(time%hours()) <= 24) then
            write(pretty_string,"(F6.2,A)") time%hours(), " hours"
        elseif (abs(time%days()) <= 365) then
            write(pretty_string,"(F10.2,A)") time%days(), " days"
        else
            write(pretty_string,"(F10.2,A,F6.1,A)") time%days(), " days (roughly ",time%days()/365.25, " years)"
        endif
    end function as_string_dt

    elemental function as_string_time(time, input_format) result(pretty_string)
        implicit none
        type(Time_type), intent(in) :: time
        character(len=*), intent(in), optional :: input_format
        character(len=kMAX_STRING_LENGTH) :: pretty_string
        character(len=kMAX_STRING_LENGTH) :: format
        integer :: i

        associate(year  => time%year,   &
                  month => time%month,  &
                  day   => time%day,    &
                  hour  => time%hour,   &
                  minute=> time%minute, &
                  second=> time%second  &
                  )

        if (present(input_format)) then
            format = input_format
        else
            ! this is the default format string to generate "YYYY/MM/DD hh:mm:ss"
            format = '(I4,"/",I0.2,"/",I0.2," ",I0.2,":",I0.2,":",I0.2)'
        endif

        ! this and the format statement above are the important bits
        write(pretty_string, format) year, month, day, hour, minute, second

        end associate

    end function as_string_time

    function as_string_timer(time, format) result(pretty_string)
        type(timer_t),    intent(in)        :: time
        character(len=*), intent(in), optional :: format

        character(len=25) :: pretty_string ! return value

        real :: temporary_time


        temporary_time = time%get_time()

        ! if the user specified a format string, use that when creating the output
        if (present(format)) then
            write(pretty_string,format) temporary_time
        else
            write(pretty_string,*) temporary_time
        endif

    end function as_string_timer

end module string
