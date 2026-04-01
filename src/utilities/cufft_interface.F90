!>------------------------------------------------------------
!!  Fortran interface module for cuFFT C library functions.
!!  Uses ISO_C_BINDING to call cuFFT routines directly.
!!
!!  Only used when _OPENACC is defined (GPU builds).
!!
!!------------------------------------------------------------
module cufft_interface
    use iso_c_binding
    implicit none

    ! cuFFT transform types
    integer(c_int), parameter :: CUFFT_Z2Z = 105  ! Double complex-to-complex

    ! cuFFT transform directions
    integer(c_int), parameter :: CUFFT_FORWARD = -1
    integer(c_int), parameter :: CUFFT_INVERSE =  1

    interface

        integer(c_int) function cufftPlan2d(plan, nx, ny, fft_type) bind(C, name='cufftPlan2d')
            import :: c_int
            integer(c_int), intent(out) :: plan
            integer(c_int), value :: nx, ny, fft_type
        end function

        integer(c_int) function cufftExecZ2Z(plan, idata, odata, direction) bind(C, name='cufftExecZ2Z')
            import :: c_int, c_ptr
            integer(c_int), value :: plan, direction
            type(c_ptr), value :: idata, odata
        end function

        integer(c_int) function cufftSetStream(plan, stream) bind(C, name='cufftSetStream')
            import :: c_int, c_ptr
            integer(c_int), value :: plan
            type(c_ptr), value :: stream
        end function

        integer(c_int) function cufftDestroy(plan) bind(C, name='cufftDestroy')
            import :: c_int
            integer(c_int), value :: plan
        end function

    end interface

end module cufft_interface
