!> ----------------------------------------------------------------------------
!!  Interface for the component event loop.
!!
!!  Interface-only parent module: the implementation lives in the
!!  flow_events_implementation submodule (flow_events.F90). Keeping the
!!  heavy USEs (time_step, initialization, nest_manager, the io objects)
!!  out of this module keeps its .mod file small, which matters because
!!  nvfortran's module import cost grows quadratically with .mod size and
!!  driver.F90 has to import this module.
!!
!! ----------------------------------------------------------------------------
module flow_events
    use options_interface,     only : options_t
    use boundary_interface,    only : boundary_t
    use flow_object_interface, only : flow_obj_t, comp_arr_t
    use ioclient_interface,    only : ioclient_t

    implicit none

    private
    public:: component_init, component_loop, component_program_end

    interface

        module subroutine component_init(component, options, boundary, ioclient, nest_index)
            class(flow_obj_t), intent(inout) :: component
            type(options_t), intent(inout) :: options(:)
            type(boundary_t), intent(inout) :: boundary
            type(ioclient_t), intent(inout) :: ioclient
            integer, intent(in) :: nest_index
        end subroutine component_init

        module subroutine component_loop(components, options, boundary, ioclient)
            type(comp_arr_t), intent(inout) :: components(:)
            type(options_t), intent(inout) :: options(:)
            type(boundary_t), intent(inout):: boundary(:)
            type(ioclient_t), intent(inout):: ioclient(:)
        end subroutine component_loop

        module subroutine component_program_end(component, options)
            type(comp_arr_t), intent(inout) :: component(:)
            type(options_t), intent(in) :: options(:)
        end subroutine component_program_end

    end interface

end module flow_events
