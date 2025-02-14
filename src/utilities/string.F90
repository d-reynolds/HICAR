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

    integer,parameter::kMAX_STRING_LENGTH=100
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

end module string
