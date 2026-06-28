!> ----------------------------------------------------------------------------
!!
!!  Interface for nest context management.
!!
!!  Interface-only parent module: the implementation lives in the
!!  nest_manager_implementation submodule (nest_manager.F90). Keeping the
!!  heavy physics-driver USEs out of this module keeps its .mod file small,
!!  which matters because nvfortran's module import cost grows quadratically
!!  with .mod size and every module above this one re-embeds what it exports.
!!
!!  @author
!!  Dylan Reynolds (d.reynolds.nw@gmail.com)
!!
!! ----------------------------------------------------------------------------
module nest_manager
    use options_interface,     only : options_t
    use domain_interface,      only : domain_t
    use flow_object_interface, only : flow_obj_t, comp_arr_t

    implicit none

    private
    public:: end_nest_context, start_nest_context, wake_nest, switch_nest_context, any_nests_not_done, nest_next_up, should_update_nests, can_update_child_nest

    interface

        module subroutine switch_nest_context(domain, options)
            type(domain_t), intent(inout) :: domain
            type(options_t), intent(inout) :: options
        end subroutine switch_nest_context

        module subroutine start_nest_context(domain, options)
            type(domain_t), intent(inout) :: domain
            type(options_t), intent(inout) :: options
        end subroutine start_nest_context

        module subroutine end_nest_context()
        end subroutine end_nest_context

        module subroutine wake_nest(domain)
            type(domain_t), intent(inout) :: domain
        end subroutine wake_nest

        module function any_nests_not_done(flow_objs) result(any_not_done)
            type(comp_arr_t), intent(in) :: flow_objs(:)
            logical :: any_not_done
        end function any_nests_not_done

        module function nest_next_up(flow_objs, options) result(next)
            type(comp_arr_t), intent(in) :: flow_objs(:)
            type(options_t), intent(in) :: options
            logical :: next
        end function nest_next_up

        module function should_update_nests(flow_objs, options) result(can_update)
            type(comp_arr_t), intent(in) :: flow_objs(:)
            type(options_t), intent(in) :: options
            logical :: can_update
        end function should_update_nests

        module function can_update_child_nest(flow_obj, child_flow_obj)
            class(flow_obj_t), intent(in) :: flow_obj
            class(flow_obj_t), intent(in) :: child_flow_obj
            logical :: can_update_child_nest
        end function can_update_child_nest

    end interface

end module nest_manager
