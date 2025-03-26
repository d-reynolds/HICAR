module flow_object_interface

    use time_object,        only : Time_type
    use time_delta_object,  only : time_delta_t
    use options_interface,     only : options_t
    use icar_constants,     only : STD_OUT_PE
    use iso_fortran_env,    only : output_unit
implicit none

private
public :: flow_obj_t


type flow_obj_t
!   private
    type(Time_type), public :: sim_time, end_time
    type(Time_type), public :: next_output, next_input, output_start, input_start
    type(time_delta_t) :: output_dt, input_dt

    logical :: started = .false.
    logical :: ended = .false.
    integer :: nest_indx = 0

    type(time_delta_t), public :: small_time_delta

    contains
        procedure, public :: init_flow_obj
        procedure, public :: reset_flow_obj_times
        procedure, public :: increment_input_time
        procedure, public :: increment_output_time
        procedure, public :: increment_sim_time
        procedure, public :: set_sim_time
        procedure, public :: time_for_input
        procedure, public :: time_for_output
        procedure, public :: next_flow_event
        procedure, public :: dead_or_asleep
        procedure, public :: check_ended
        procedure, public :: check_started
    end type

    interface

        module subroutine init_flow_obj(this, options, nest_indx)
            implicit none
            class(flow_obj_t), intent(inout) :: this
            type(options_t), intent(in) :: options
            integer, intent(in) :: nest_indx
        end subroutine

        module subroutine reset_flow_obj_times(this, input, output)
            implicit none
            class(flow_obj_t), intent(inout) :: this
            logical, optional, intent(in) :: input, output
        end subroutine

        module subroutine increment_output_time(this)
            implicit none
            class(flow_obj_t), intent(inout) :: this
        end subroutine

        module subroutine increment_input_time(this)
            implicit none
            class(flow_obj_t),   intent(inout)  :: this
        end subroutine

        module subroutine increment_sim_time(this, dt)
            implicit none
            class(flow_obj_t), intent(inout) :: this
            type(time_delta_t), intent(in) :: dt
        end subroutine

        module subroutine set_sim_time(this, time)
            implicit none
            class(flow_obj_t), intent(inout) :: this
            type(Time_type), intent(in) :: time
        end subroutine

        module function time_for_input(this)
            implicit none
            class(flow_obj_t), intent(in) :: this
            logical :: time_for_input
        end function

        module function time_for_output(this)
            implicit none
            class(flow_obj_t), intent(in) :: this
            logical :: time_for_output
        end function

        module function next_flow_event(this)
            implicit none
            class(flow_obj_t), intent(in) :: this
            type(Time_type) :: next_flow_event
        end function

        module function dead_or_asleep(this) result(doa)
            implicit none
            class(flow_obj_t), intent(inout) :: this
            logical :: doa
        end function

        module subroutine check_ended(this)
            implicit none
            class(flow_obj_t), intent(inout) :: this
        end subroutine

        module subroutine check_started(this)
            implicit none
            class(flow_obj_t), intent(inout) :: this
        end subroutine

    end interface

end module