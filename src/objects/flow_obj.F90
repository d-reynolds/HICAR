submodule(flow_object_interface) meta_data_implementation
  implicit none

contains

    ! add a given attribute name/value pair to the metadata object
    module subroutine init_flow_obj(this, options, nest_indx)
        class(flow_obj_t), intent(inout) :: this
        type(options_t), intent(in) :: options
        integer, intent(in) :: nest_indx


        this%sim_time = options%general%start_time
        this%end_time = options%general%end_time

        this%input_dt = options%forcing%input_dt
        this%output_dt = options%output%output_dt

        this%next_output = options%general%start_time
        this%next_input = options%general%start_time

        this%started = .false.
        this%ended = .false.
        this%last_loop = .false.
        this%nest_indx = nest_indx

        call this%small_time_delta%set(1)

    end subroutine init_flow_obj

    ! increment the output time
    module subroutine increment_output_time(this)
        implicit none
        class(flow_obj_t), intent(inout) :: this

        ! check if the next output time is greater than the end time
        if (.not.(this%time_for_output())) then
            if (STD_OUT_PE) write(*,*) "For nest: ", this%nest_indx, " we were asked to increment the output time."
            if (STD_OUT_PE) write(*,*) "Currently, output_time is: ",this%next_output%as_string(), "and sim_time is: ", this%sim_time%as_string()
            if (STD_OUT_PE) write(*,*) "At this step, they should be equal."
            stop "CONTROL FLOW ERROR, EXITING"
        end if

        this%next_output = this%next_output + this%output_dt

        ! Because the output time also changes the end condition, we need to check if the end condition is now true
        call this%check_ended()

    end subroutine increment_output_time

    ! increment the input time
    module subroutine increment_input_time(this)
        implicit none
        class(flow_obj_t), intent(inout) :: this

        ! check if the next input time is greater than the end time
        if (.not.(this%time_for_input())) then
            if (STD_OUT_PE) write(*,*) "For nest: ", this%nest_indx, " we were asked to increment the input time."
            if (STD_OUT_PE) write(*,*) "Currently, next_input is: ",this%next_input%as_string(), "and sim_time is: ", this%sim_time%as_string()
            if (STD_OUT_PE) write(*,*) "At this step, they should be equal."
            stop "CONTROL FLOW ERROR, EXITING"
        end if

        this%next_input = this%next_input + this%input_dt

        ! We will always read in before outputting -- the model does not generate data
        ! therefore we only see if we should start the model here, and not in output.
        ! If we were not started before, then this is the first input. We do not want
        ! to log an increment to next_input, since we need to interpolate temporally,
        ! and therefore need the next input time already.
        if (this%started .eqv. .False.) then
            this%started = .true.
            this%next_input = this%next_input - this%input_dt
        endif

    end subroutine increment_input_time

    module subroutine increment_sim_time(this, dt)
        implicit none
        class(flow_obj_t), intent(inout) :: this
        type(time_delta_t), intent(in) :: dt

        if (this%started .eqv. .False.) then
            write(*,*) "For nest: ", this%nest_indx, " we were asked to increment the sim time."
            write(*,*) "But this nest has not started yet."
            stop "CONTROL FLOW ERROR, EXITING"
        endif

        this%sim_time = this%sim_time + dt

        call this%check_ended()

    end subroutine increment_sim_time

    module subroutine set_sim_time(this, time)
        implicit none
        class(flow_obj_t), intent(inout) :: this
        type(Time_type), intent(in) :: time

        if (this%started .eqv. .False.) then
            write(*,*) "For nest: ", this%nest_indx, " we were asked to increment the sim time."
            write(*,*) "But this nest has not started yet."
            stop "CONTROL FLOW ERROR, EXITING"
        endif

        this%sim_time = time

        call this%check_ended()

    end subroutine set_sim_time

    module subroutine check_ended(this)
        implicit none
        class(flow_obj_t), intent(inout) :: this

        if (this%sim_time + this%small_time_delta > this%end_time) then
            this%sim_time = this%end_time
            if (this%next_output > this%end_time) then
                this%ended = .true.
            else
                this%last_loop = .true.
            end if
        end if

    end subroutine check_ended

    function dead_or_asleep(this) result(doa)
        implicit none
        class(flow_obj_t), intent(in) :: this
        logical :: doa

        doa = .false.

        if (this%ended .eqv. .true.) then
            doa = .true.
        else if (this%started .eqv. .false.) then
            doa = .true.
        end if

    end function dead_or_asleep

    ! check if the input time is equal to the next input time
    module function time_for_input(this)
        implicit none
        class(flow_obj_t), intent(in) :: this
        logical :: time_for_input

        time_for_input = .false.

        if (this%ended) then
            if (STD_OUT_PE) write(*,*) "For nest: ", this%nest_indx, " we were asked to check the input time, but this nest has ended."
            time_for_input = .false.
            return
        end if

        time_for_input = ((this%sim_time  + this%small_time_delta >= this%next_input)) ! .and. (this%sim_time + this%small_time_delta <= this%next_input + this%input_dt))

    end function time_for_input

    ! check if the output time is equal to the next output time
    module function time_for_output(this)
        implicit none
        class(flow_obj_t), intent(in) :: this
        logical :: time_for_output

        time_for_output = .false.

        if (this%ended) then
            if (STD_OUT_PE) write(*,*) "For nest: ", this%nest_indx, " we were asked to check the output time, but this nest has ended."
            time_for_output = .false.
            return
        end if

        if (.not.(this%started)) then
            if (STD_OUT_PE) write(*,*) "For nest: ", this%nest_indx, " we were asked to check the output time, but this nest has not started."
            time_for_output = .false.
            return
        end if

        if (this%last_loop) then
            time_for_output = .true.
            return
        end if

        time_for_output = ((this%sim_time + this%small_time_delta >= this%next_output) .and. (this%sim_time + this%small_time_delta <= this%next_output + this%input_dt))

    end function time_for_output

    ! check if the next output time is greater than the end time
    module function next_flow_event(this)
        implicit none
        class(flow_obj_t), intent(in) :: this
        type(Time_type) :: next_flow_event

        if (this%ended) then
            if (STD_OUT_PE) write(*,*) "For nest: ", this%nest_indx, " we were asked to check the next flow event, but this nest has ended."
            stop "CONTROL FLOW ERROR, EXITING"
        end if

        if (.not.(this%started)) then
            if (STD_OUT_PE) write(*,*) "For nest: ", this%nest_indx, " we were asked to check the next flow event, but this nest has not started."
            stop "CONTROL FLOW ERROR, EXITING"
        end if

        if (this%next_input < this%next_output) then
            next_flow_event = this%next_input
        else
            next_flow_event = this%next_output
        end if
        if (next_flow_event > this%end_time) then
            next_flow_event = this%end_time
        end if

        if (next_flow_event < this%sim_time) then
            if (STD_OUT_PE) write(*,*) "For nest: ", this%nest_indx, " we were asked to check the next flow event, but the next flow event is less than the sim time."
            stop "CONTROL FLOW ERROR, EXITING"
        end if
        if (next_flow_event > this%end_time) then
            if (STD_OUT_PE) write(*,*) "For nest: ", this%nest_indx, " we were asked to check the next flow event, but the next flow event is greater than the end time."
            stop "CONTROL FLOW ERROR, EXITING"
        end if
    end function next_flow_event

end submodule