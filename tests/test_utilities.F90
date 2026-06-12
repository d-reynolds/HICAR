!>------------------------------------------------------------
!! Test suite for small utility modules, salvaged from the ICAR-era
!! standalone tests:
!!   - string round trips        (test_string.F90)
!!   - linear_space              (test_array_utilities.F90)
!!   - var_dict_t add/get        (test_variable_dict.F90)
!!   - fftshift / ifftshift      (test_fftshift.F90 — which only
!!     printed; the assertions here are new: fftshift must equal
!!     cshift(x, n/2) and ifftshift must invert it. fftshifter is
!!     still live via linear_theory_winds.)
!!------------------------------------------------------------
module test_utilities

    use iso_fortran_env,         only : real32, real64
    use icar_constants,          only : kVARS
    use string,                  only : get_double, get_real, get_integer, str
    use array_utilities,         only : linear_space
    use variable_dict_interface, only : var_dict_t
    use variable_interface,      only : variable_t
    use grid_interface,          only : grid_t
    use fftshifter,              only : fftshift, ifftshift
    use testdrive,               only : new_unittest, unittest_type, error_type, check

    implicit none
    private

    public :: collect_utilities_suite

contains

    subroutine collect_utilities_suite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            new_unittest("string_roundtrip", test_string_roundtrip), &
            new_unittest("linear_space",     test_linear_space), &
            new_unittest("variable_dict",    test_variable_dict), &
            new_unittest("fftshift_1d",      test_fftshift_1d), &
            new_unittest("fftshift_2d",      test_fftshift_2d) &
          ]

    end subroutine collect_utilities_suite


    subroutine test_string_roundtrip(error)
        type(error_type), allocatable, intent(out) :: error

        real(real64), parameter :: real64_datum = 0.3_real64
        real(real32), parameter :: real32_datum = 0.3_real32
        integer,      parameter :: integer_datum = 1234567
        character(len=32) :: real64_string, real32_string, integer_string

        write(real64_string,*) real64_datum
        write(real32_string,*) real32_datum
        write(integer_string,*) integer_datum

        call check(error, get_double(real64_string) == real64_datum,   "get_double round trip")
        if (allocated(error)) return
        call check(error, get_real(real32_string) == real32_datum,     "get_real round trip")
        if (allocated(error)) return
        call check(error, get_integer(integer_string) == integer_datum,"get_integer round trip")
        if (allocated(error)) return
        call check(error, adjustl(real64_string) == str(real64_datum), "str(real64) matches list-directed write")
        if (allocated(error)) return
        call check(error, adjustl(real32_string) == str(real32_datum), "str(real32) matches list-directed write")
        if (allocated(error)) return
        call check(error, adjustl(integer_string) == str(integer_datum),"str(integer) matches list-directed write")
    end subroutine test_string_roundtrip


    subroutine test_linear_space(error)
        type(error_type), allocatable, intent(out) :: error

        real, allocatable :: a(:)
        real, parameter :: vmin = 1e-4, vmax = 3.0, tol = 1e-6
        integer, parameter :: n = 100
        real :: dv
        integer :: i

        call linear_space(a, vmin, vmax, n)

        call check(error, allocated(a), "linear_space must allocate the array")
        if (allocated(error)) return
        call check(error, size(a) == n, "linear_space array has the requested size")
        if (allocated(error)) return
        call check(error, abs(a(1) - vmin) < tol, "linear_space starts at vmin")
        if (allocated(error)) return
        call check(error, abs(a(n) - vmax) < tol, "linear_space ends at vmax")
        if (allocated(error)) return

        dv = (vmax - vmin) / (n - 1)
        do i = 2, n
            if (abs((a(i) - a(i-1)) - dv) > tol) then
                call check(error, .false., "linear_space step size is not uniform")
                return
            endif
        enddo
        call check(error, .true., "unreachable")
    end subroutine test_linear_space


    subroutine test_variable_dict(error)
        type(error_type), allocatable, intent(out) :: error

        type(var_dict_t) :: var_collection
        type(variable_t) :: var1, var2, var3, output_var
        type(grid_t)     :: grid

        ! initialize requires a kVARS metadata index (test_driver calls
        ! initialize_var_constants at startup)
        call grid%set_grid_dimensions(nx=10, ny=10, nz=5)
        call var1%initialize(kVARS%potential_temperature, grid)

        call grid%set_grid_dimensions(nx=2, ny=15, nz=5)
        call var2%initialize(kVARS%pressure, grid)

        call var3%initialize(kVARS%skin_temperature, [2,5])

        ! the dictionary is keyed by kVARS id (integer) in current HICAR
        call var_collection%add_var(kVARS%potential_temperature, var1)
        call var_collection%add_var(kVARS%pressure,              var2)
        call var_collection%add_var(kVARS%skin_temperature,      var3)

        ! retrieving the middle entry must return that variable's metadata
        output_var = var_collection%get_var(kVARS%pressure)

        call check(error, output_var%three_d .eqv. var2%three_d, &
                   "var_dict returned wrong dimensionality flag")
        if (allocated(error)) return
        call check(error, size(output_var%dim_len) == size(var2%dim_len), &
                   "var_dict returned wrong number of dimensions")
        if (allocated(error)) return
        call check(error, all(output_var%dim_len == var2%dim_len), &
                   "var_dict returned wrong dim_len")
    end subroutine test_variable_dict


    subroutine test_fftshift_1d(error)
        type(error_type), allocatable, intent(out) :: error

        real    :: r_odd(5), r_even(6), r_orig_odd(5), r_orig_even(6)
        complex :: c_odd(5), c_orig_odd(5)
        integer :: i

        r_orig_odd  = [(real(i), i = 1, 5)]
        r_orig_even = [(real(i), i = 1, 6)]
        c_orig_odd  = cmplx(r_orig_odd, -r_orig_odd)

        ! fftshift is a cyclic shift by n/2 (matlab/numpy convention)
        r_odd = r_orig_odd
        call fftshift(r_odd)
        call check(error, all(r_odd == cshift(r_orig_odd, size(r_odd)/2)), &
                   "1d real fftshift (odd n) is not cshift(x, n/2)")
        if (allocated(error)) return
        call ifftshift(r_odd)
        call check(error, all(r_odd == r_orig_odd), &
                   "1d real ifftshift does not invert fftshift (odd n)")
        if (allocated(error)) return

        r_even = r_orig_even
        call fftshift(r_even)
        call check(error, all(r_even == cshift(r_orig_even, size(r_even)/2)), &
                   "1d real fftshift (even n) is not cshift(x, n/2)")
        if (allocated(error)) return
        call ifftshift(r_even)
        call check(error, all(r_even == r_orig_even), &
                   "1d real ifftshift does not invert fftshift (even n)")
        if (allocated(error)) return

        c_odd = c_orig_odd
        call fftshift(c_odd)
        call check(error, all(c_odd == cshift(c_orig_odd, size(c_odd)/2)), &
                   "1d complex fftshift is not cshift(x, n/2)")
        if (allocated(error)) return
        call ifftshift(c_odd)
        call check(error, all(c_odd == c_orig_odd), &
                   "1d complex ifftshift does not invert fftshift")
    end subroutine test_fftshift_1d


    subroutine test_fftshift_2d(error)
        type(error_type), allocatable, intent(out) :: error

        real    :: r2(5,6), r2_orig(5,6), r2_expect(5,6)
        complex :: c2(5,6), c2_orig(5,6)
        integer :: i, j

        do j = 1, 6
            do i = 1, 5
                r2_orig(i,j) = i + j*50.0
            enddo
        enddo
        c2_orig = cmplx(r2_orig, -r2_orig)

        ! 2d fftshift = cyclic shift by n/2 along each dimension
        r2_expect = cshift(cshift(r2_orig, size(r2_orig,1)/2, dim=1), &
                           size(r2_orig,2)/2, dim=2)

        r2 = r2_orig
        call fftshift(r2)
        call check(error, all(r2 == r2_expect), &
                   "2d real fftshift is not cshift by n/2 in both dims")
        if (allocated(error)) return
        call ifftshift(r2)
        call check(error, all(r2 == r2_orig), &
                   "2d real ifftshift does not invert fftshift")
        if (allocated(error)) return

        c2 = c2_orig
        call fftshift(c2)
        call ifftshift(c2)
        call check(error, all(c2 == c2_orig), &
                   "2d complex ifftshift does not invert fftshift")
    end subroutine test_fftshift_2d

end module test_utilities
