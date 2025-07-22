submodule(flow_object_interface) flow_obj_implementation
    use icar_constants,     only : STD_OUT_PE
    use iso_fortran_env,    only : output_unit
    use string,             only : as_string
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
        this%nest_indx = nest_indx

        call this%small_time_delta%set(1)

        if (options%restart%restart) then
            this%sim_time = options%restart%restart_time

            ! need to correct input and output times
            ! it should be the next time after the restart time 
            ! that is an integer multiple of the input_dt
            ! and options%general%start_time            
            this%next_input = options%general%start_time
            do while(this%next_input < options%restart%restart_time)
                call this%next_input%set(this%next_input%mjd() + this%input_dt%days())
            end do

            !By definition, a restart time has to lay on an output time, so this is valid
            call this%next_output%set(this%next_output%mjd() + this%output_dt%days())

        endif

    end subroutine init_flow_obj

    ! increment the output time
    module subroutine increment_output_time(this)
        implicit none
        class(flow_obj_t), intent(inout) :: this

        ! check if the next output time is greater than the end time
        if (.not.(this%time_for_output())) then
            if (STD_OUT_PE) write(*,*) "For nest: ", this%nest_indx, " we were asked to increment the output time."
            if (STD_OUT_PE) write(*,*) "Currently, output_time is: ",trim(as_string(this%next_output)), " and sim_time is: ", trim(as_string(this%sim_time))
            if (STD_OUT_PE) write(*,*) "At this step, they should be equal."
            stop "CONTROL FLOW ERROR, EXITING 1"
        end if

        call this%next_output%set(this%next_output%mjd() + this%output_dt%days())

        ! Because the output time also changes the end condition, we need to check if the end condition is now true
        call this%check_ended()

    end subroutine increment_output_time

    ! increment the input time
    module subroutine increment_input_time(this)
        implicit none
        class(flow_obj_t), intent(inout) :: this

        type(Time_type) :: time_tmp
        ! check if the next input time is greater than the end time
        ! if (.not.(this%time_for_input())) then
        !     if (STD_OUT_PE) write(*,*) "For nest: ", this%nest_indx, " we were asked to increment the input time."
        !     if (STD_OUT_PE) write(*,*) "Currently, next_input is: ",trim(this%next_input%as_string()), " and sim_time is: ", trim(this%sim_time%as_string())
        !     if (STD_OUT_PE) write(*,*) "At this step, they should be equal."
        !     stop "CONTROL FLOW ERROR, EXITING 2"
        ! end if

        call time_tmp%set(this%sim_time%mjd() + this%small_time_delta%days())
        if (this%next_input <= time_tmp) then
            call this%next_input%set(this%sim_time%mjd() + this%input_dt%days())
        endif

    end subroutine increment_input_time

    module subroutine increment_sim_time(this, dt)
        implicit none
        class(flow_obj_t), intent(inout) :: this
        type(time_delta_t), intent(in) :: dt

        if (this%started .eqv. .False.) then
            write(*,*) "For nest: ", this%nest_indx, " we were asked to increment the sim time."
            write(*,*) "But this nest has not started yet."
            stop "CONTROL FLOW ERROR, EXITING 3"
        endif

        call this%sim_time%set(this%sim_time%mjd() + dt%days())

        call this%check_ended()

    end subroutine increment_sim_time

    module subroutine set_sim_time(this, time)
        implicit none
        class(flow_obj_t), intent(inout) :: this
        type(Time_type), intent(in) :: time

        if (this%started .eqv. .False.) then
            write(*,*) "For nest: ", this%nest_indx, " we were asked to increment the sim time."
            write(*,*) "But this nest has not started yet."
            stop "CONTROL FLOW ERROR, EXITING 4"
        endif

        this%sim_time = time

        call this%check_ended()

    end subroutine set_sim_time

    module subroutine check_ended(this)
        implicit none
        class(flow_obj_t), intent(inout) :: this
        type(Time_type) :: time_tmp

        if (this%ended) return
        
        call time_tmp%set(this%sim_time%mjd() + this%small_time_delta%days())

        if (time_tmp > this%end_time) then
            this%sim_time = this%end_time
            if (this%next_output > this%end_time) then
                this%ended = .true.
            end if
        end if

    end subroutine check_ended


    module function dead_or_asleep(this) result(doa)
        implicit none
        class(flow_obj_t), intent(inout) :: this
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
            time_for_input = .false.
            return
        end if

        time_for_input = this%sim_time == this%next_input

    end function time_for_input

    ! check if the output time is equal to the next output time
    module function time_for_output(this)
        implicit none
        class(flow_obj_t), intent(in) :: this
        logical :: time_for_output

        time_for_output = .false.

        if (this%ended) then
            time_for_output = .false.
            return
        end if

        time_for_output = this%sim_time == this%next_output

    end function time_for_output

    ! check if the next output time is greater than the end time
    module function next_flow_event(this)
        implicit none
        class(flow_obj_t), intent(in) :: this
        type(Time_type) :: next_flow_event

        if (this%ended) then
            if (STD_OUT_PE) write(*,*) "For nest: ", this%nest_indx, " we were asked to check the next flow event, but this nest has ended."
            stop "CONTROL FLOW ERROR, EXITING 5"
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
            stop "CONTROL FLOW ERROR, EXITING 7"
        end if
        if (next_flow_event > this%end_time) then
            if (STD_OUT_PE) write(*,*) "For nest: ", this%nest_indx, " we were asked to check the next flow event, but the next flow event is greater than the end time."
            stop "CONTROL FLOW ERROR, EXITING 8"
        end if
    end function next_flow_event

end submodule