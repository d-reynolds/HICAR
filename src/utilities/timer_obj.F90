submodule(timer_interface) timer_implementation

    implicit none

contains

    !>-----------------------------------
    !! Start the timer
    !!
    !! Sets the internal start_time and marks the timer as running
    !!
    !------------------------------------
    module subroutine start(this, use_cpu_time)
        class(timer_t), intent(inout) :: this
        logical,        intent(in),   optional :: use_cpu_time

        logical :: cpu_only
        if (present(use_cpu_time)) then
          cpu_only = use_cpu_time
        else
          cpu_only = .false.
        end if

        this%is_running = .True.
        this%use_cpu_time = cpu_only
        if (cpu_only) then
            call cpu_time(this%start_time)
        else
            call system_clock(this%counter, this%count_rate, this%count_max)
            this%start_time = 0
        endif
    end subroutine start

    !>-----------------------------------
    !! Stop the timer
    !!
    !! Sets the internal stop_time and marks the timer as not running
    !!
    !! Also computes the total_time the timer has been running
    !!
    !------------------------------------
    module subroutine stop(this)
        class(timer_t), intent(inout) :: this
        integer :: count_end

        ! do this before we even test anything else because we want the timer to "stop" as soon as possible
        if (this%use_cpu_time) then
            call cpu_time(this%end_time)
        else
            call system_clock(count_end)
            if (count_end<this%counter) then
                this%end_time = (count_end + (this%count_max - this%counter)) / real(this%count_rate)
            else
                this%end_time = (count_end - this%counter) / real(this%count_rate)
            endif
        endif

        if (this%is_running) then
            this%is_running = .False.
            this%total_time = this%total_time + (this%end_time - this%start_time)
        endif
    end subroutine stop

    !>-----------------------------------
    !! Reset the timer
    !!
    !! Resets timer internal variables as if it was never running
    !!
    !------------------------------------
    module subroutine reset(this)
        class(timer_t), intent(inout) :: this

        this%total_time = 0
        this%start_time = 0
        this%end_time   = 0
        this%is_running = .False.
    end subroutine reset

    !>-----------------------------------
    !! Return the time as a real
    !!
    !! If the timer is running, it includes the current time in the total reported
    !!
    !------------------------------------
    module function get_time(this) result(time)
        class(timer_t),    intent(inout)        :: this

        real :: time ! return value

        real :: current_time
        integer :: count_end

        if (this%is_running) then
            if (this%use_cpu_time) then
                call cpu_time(current_time)
            else
                call system_clock(count_end)
                if (count_end<this%counter) then
                    current_time = (count_end + (this%count_max - this%counter)) / real(this%count_rate)
                else
                    current_time = (count_end - this%counter) / real(this%count_rate)
                endif
            endif
            time = this%total_time + (current_time - this%start_time)
        else
            time = this%total_time
        endif

    end function get_time


    !>-----------------------------------
    !! Return the time as a string
    !!
    !! If the timer is running, it includes the current time in the total reported
    !!
    !------------------------------------
    module function as_string(this, format) result(time)
        class(timer_t),    intent(inout)        :: this
        character(len=*), intent(in), optional :: format

        character(len=25) :: time ! return value

        real :: temporary_time

        ! if (this%is_running) then
        !     if (this%use_cpu_time) then
        !         call cpu_time(current_time)
        !     else
        !         call system_clock(count_end)
        !         if (count_end<this%counter) then
        !             current_time = (count_end + (this%count_max - this%counter)) / real(this%count_rate)
        !         else
        !             current_time = (count_end - this%counter) / real(this%count_rate)
        !         endif
        !     endif
        !     temporary_time = this%total_time + (current_time - this%start_time)
        ! else
        !     temporary_time = this%total_time
        ! endif

        temporary_time = this%get_time()

        ! if the user specified a format string, use that when creating the output
        if (present(format)) then
            write(time,format) temporary_time
        else
            write(time,*) temporary_time
        endif

    end function as_string

    module function timer_mean(this,comms) result(mean_t)
        implicit none
        class(timer_t), intent(inout) :: this
        type(MPI_Comm), intent(in) :: comms

        real :: mean_t, t_sum
        integer :: ierr, NUM_COMPUTE
            
        call this%stop()

        t_sum = this%get_time()
        call MPI_Allreduce(MPI_IN_PLACE,t_sum,1,MPI_REAL,MPI_SUM,comms,ierr)
        call MPI_Comm_Size(comms,NUM_COMPUTE)
        mean_t = t_sum/NUM_COMPUTE
    
    end function

    module function timer_max(this,comms) result(max_t)
        implicit none
        class(timer_t), intent(inout) :: this
        type(MPI_Comm), intent(in) :: comms

        real :: max_t
        integer :: ierr
            
        call this%stop()

        max_t = this%get_time()
        call MPI_Allreduce(MPI_IN_PLACE,max_t,1,MPI_REAL,MPI_MAX,comms,ierr)
    
    end function

    module function timer_min(this,comms) result(min_t)
        implicit none
        class(timer_t), intent(inout) :: this
        type(MPI_Comm), intent(in) :: comms

        real :: min_t
        integer :: ierr
            
        call this%stop()

        min_t = this%get_time()
        call MPI_Allreduce(MPI_IN_PLACE,min_t,1,MPI_REAL,MPI_MIN,comms,ierr)
    
    end function


end submodule timer_implementation
