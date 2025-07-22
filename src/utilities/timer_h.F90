!> -----------------------
!! Provides a timer object that can be used to monitor how long sections of code are taking
!!
!! Similar functionality can be gained from a profiling system, but sometimes profiling is too much detail
!!
!! Note that the current implementation utilizes the fortran CPU_time function.  This function only counts time actually
!! spent running on the CPU, not wallclock.  In addition, the time for all threads in the process are included.
!!
!! ------------------------
module timer_interface

    use mpi_f08

    implicit none
    private
    public :: timer_t

    type timer_t
        private
        real    :: start_time = 0
        real    :: end_time   = 0
        real    :: total_time = 0
        integer :: counter    = 0
        integer :: count_max  = 86000
        integer :: count_rate = 1
        logical :: use_cpu_time = .False.
        logical :: is_running = .False.
      contains
        procedure :: start
        procedure :: stop
        ! procedure :: as_string
        procedure :: get_time
        procedure :: reset
        procedure :: mean => timer_mean
        procedure :: min => timer_min
        procedure :: max => timer_max
    end type
interface

    !> -----------------------------------
    !! Start the timer
    !!
    !! Sets the internal start_time and marks the timer as running
    !!
    !! ------------------------------------
    module subroutine start(this, use_cpu_time)
        implicit none
        class(timer_t), intent(inout) :: this
        logical,        intent(in),   optional :: use_cpu_time
    end subroutine start

    !> -----------------------------------
    !! Stop the timer
    !!
    !! Sets the internal stop_time and marks the timer as not running
    !!
    !! Also computes the total_time the timer has been running
    !!
    !! ------------------------------------
    module subroutine stop(this)
        implicit none
        class(timer_t), intent(inout) :: this
    end subroutine stop

    !> -----------------------------------
    !! Reset the timer
    !!
    !! Resets timer internal variables as if it was never running
    !!
    !! ------------------------------------
    module subroutine reset(this)
        implicit none
        class(timer_t), intent(inout) :: this
    end subroutine reset

    !> -----------------------------------
    !! Return the time as a real
    !!
    !! If the timer is running, it includes the current time in the total reported
    !!
    !! ------------------------------------
    module function get_time(this) result(time)
        implicit none
        class(timer_t),    intent(in)        :: this
        real :: time
    end function get_time


    !> -----------------------------------
    !! Return the time as a string
    !!
    !! If the timer is running, it includes the current time in the total reported
    !!
    !! ------------------------------------
    ! module function as_string(this, format) result(time)
    !     implicit none
    !     class(timer_t),   intent(inout)        :: this
    !     character(len=*), intent(in), optional :: format
    !     character(len=25) :: time ! return value
    ! end function as_string

    !> -----------------------------------
    !! Return the mean time as a real
    !!
    !! If the timer is running, it includes the current time in the total reported
    !!
    !! ------------------------------------
    module function timer_mean(this,comms) result(mean_t)
        implicit none
        class(timer_t), intent(inout) :: this
        type(MPI_Comm), intent(in) :: comms
        real :: mean_t
    end function timer_mean
    !> -----------------------------------
    !! Return the min time as a real
    !!
    !! If the timer is running, it includes the current time in the total reported
    !!
    !! ------------------------------------
    module function timer_min(this,comms) result(min_t)
        implicit none
        class(timer_t), intent(inout) :: this
        type(MPI_Comm), intent(in) :: comms
        real :: min_t
    end function timer_min
    !> -----------------------------------
    !! Return the max time as a real
    !!
    !! If the timer is running, it includes the current time in the total reported
    !!
    !! ------------------------------------
    module function timer_max(this,comms) result(max_t)
        implicit none
        class(timer_t), intent(inout) :: this
        type(MPI_Comm), intent(in) :: comms
        real :: max_t
    end function timer_max
end interface
end module timer_interface
