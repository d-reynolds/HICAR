!> ----------------------------------------------------------------------------
!!  
!!  Manages changing of nest contexts. This is mostly handeled by changing 
!!  module level indices and arrays for physics drivers when entering a new
!!  nest context, and finalizing some physics routines when leaving a nest
!!  context.
!! 
!!  @author
!!  Dylan Reynolds (d.reynolds.nw@gmail.com)
!!
!! ----------------------------------------------------------------------------
module nest_manager
    use options_interface,  only : options_t
    use domain_interface,   only : domain_t
    use flow_object_interface, only : flow_obj_t, comp_arr_t
    use initialization,             only : init_model_state
    use microphysics,               only : mp_init
    use advection,                  only : adv_init
    use radiation,                  only : radiation_init
    use convection,                 only : init_convection
    use planetary_boundary_layer,   only : pbl_init
    use land_surface,               only : lsm_init
    use surface_layer,              only : sfc_init
    use wind,                       only : init_winds
    use icar_constants
    use iso_fortran_env


    ! use io_routines,                only : io_read, &
    !                                        io_write3d,io_write3di, io_write
    ! use geo,                        only : geo_LUT, geo_interp, geo_interp2d, standardize_coordinates
    ! use vertical_interpolation,     only : vLUT, vinterp
    ! use wind,                       only : init_winds
    ! use initialize_options,         only : init_options
    ! use string,                     only : str


    implicit none
    private
    integer :: current_nest = -1

    public:: end_nest_context, start_nest_context, wake_nest, switch_nest_context, any_nests_not_done, nest_next_up, should_update_nests, can_update_child_nest


contains

    subroutine switch_nest_context(domain, options)
        implicit none
        type(domain_t), intent(inout) :: domain
        type(options_t), intent(inout) :: options

        if (current_nest == -1) then
            call start_nest_context(domain, options)
        elseif (current_nest == domain%nest_indx) then
            return
        endif

        if (current_nest /= domain%nest_indx) then
            call end_nest_context()
            call start_nest_context(domain, options)
        end if

    end subroutine switch_nest_context

    !> ----------------------------------------------------------------------------
    !!  @subroutine start_nest_context(domain, options)
    !!
    !!  @brief
    !!  Initializes the physics drivers for a new nest context.
    !!
    !!  @param domain
    !!  The domain_t object for the nest context.
    !!  @param options  
    !!  The options_t object for the nest context.
    !! ----------------------------------------------------------------------------
    subroutine start_nest_context(domain, options)
        implicit none
        type(domain_t), intent(inout) :: domain
        type(options_t), intent(inout) :: options

        if (current_nest == domain%nest_indx) then
            if (STD_OUT_PE) write(*,*) "WARNING: Nest context already active for nest ", current_nest
            return
        endif

        call init_winds(domain,options, context_chng=.True.)

        ! initialize microphysics code (e.g. compute look up tables in Thompson et al)
        call mp_init(domain, options, context_chng=.True.) !this could easily be moved to init_model...
        call init_convection(domain,options, context_chng=.True.)
        call pbl_init(domain,options, context_chng=.True.)
        call radiation_init(domain,options, context_chng=.True.)
        call lsm_init(domain,options, context_chng=.True.)
        call sfc_init(domain,options, context_chng=.True.)
        call adv_init(domain,options, context_chng=.True.)

        current_nest = domain%nest_indx
    end subroutine start_nest_context

    subroutine end_nest_context()
        implicit none

        current_nest = -1
        
    end subroutine  end_nest_context

    subroutine wake_nest(domain)
        implicit none
        type(domain_t), intent(inout) :: domain


        ! If we have an active context and it is not us, end the context
        if ( (current_nest > 0) .and. (current_nest /= domain%nest_indx) ) then
            call end_nest_context()
        endif


        current_nest = domain%nest_indx

    end subroutine wake_nest

    function any_nests_not_done(flow_objs) result(any_not_done)
        implicit none
        type(comp_arr_t), intent(in) :: flow_objs(:)
        logical :: any_not_done
        integer :: i

        any_not_done = .false.
        do i = 1, size(flow_objs)
            if (allocated(flow_objs(i)%comp)) then
                if (flow_objs(i)%comp%ended .eqv. .false.) then
                    any_not_done = .true.
                    return
                end if
            endif
        end do

    end function any_nests_not_done

    function nest_next_up(flow_objs, options) result(next)
        implicit none
        type(comp_arr_t), intent(in) :: flow_objs(:)
        type(options_t), intent(in) :: options
        logical :: next

        next = .false.

        ! If we are a child nest, not at the end of our run time...
        if (flow_objs(options%nest_indx)%comp%started .eqv. .False.) then
            ! ... and if we will be running on the next iteration...
            !safety check before indexing into flow_obj
            if (options%general%parent_nest == 0 .or. options%restart%restart) then
                next = .true.
            else if (flow_objs(options%nest_indx)%comp%sim_time - flow_objs(options%nest_indx)%comp%small_time_delta <= flow_objs(options%general%parent_nest)%comp%sim_time) then
                next = .true.
            endif
        end if

    end function nest_next_up

    function should_update_nests(flow_objs, options) result(can_update)
        implicit none
        type(comp_arr_t), intent(in) :: flow_objs(:)
        type(options_t), intent(in) :: options
        logical :: can_update
        integer :: n, num_children

        can_update = .false.

        num_children = size(options%general%child_nests)
        ! if we even have nests
        if (num_children > 0) then
            ! if we are at the end of an input step, or we have just ended, or we have not started (this would mean we are in wake_component)
            if (flow_objs(options%nest_indx)%comp%time_for_input() .or. flow_objs(options%nest_indx)%comp%sim_time%equals(flow_objs(options%nest_indx)%comp%end_time) .or. (flow_objs(options%nest_indx)%comp%started .eqv. .False.)) then
                ! loop over children
                do n = 1, num_children
                    ! If we are at or ahead of our child's time, and the child has not ended, then our child needs to be updated
                    if ( flow_objs(options%nest_indx)%comp%sim_time >= flow_objs(options%general%child_nests(n))%comp%sim_time - flow_objs(options%nest_indx)%comp%small_time_delta) then
                        if (flow_objs(options%general%child_nests(n))%comp%ended .eqv. .False.) then
                            can_update = .true.
                            return
                        end if
                    end if
                end do
            end if
        endif

    end function should_update_nests

    function ahead_or_at_child(flow_obj, child_flow_obj) result(ahead)
        implicit none
        class(flow_obj_t), intent(in) :: flow_obj
        class(flow_obj_t), intent(in) :: child_flow_obj
        logical :: ahead

        ahead = .false.

        if (flow_obj%sim_time >= child_flow_obj%sim_time - flow_obj%small_time_delta) then
            ahead = .true.
        end if

    end function ahead_or_at_child

    function can_update_child_nest(flow_obj, child_flow_obj)
        implicit none
        class(flow_obj_t), intent(in) :: flow_obj
        class(flow_obj_t), intent(in) :: child_flow_obj
        logical :: can_update_child_nest

        can_update_child_nest = (flow_obj%sim_time >= child_flow_obj%sim_time - flow_obj%small_time_delta .and. .not.(child_flow_obj%ended))

    end function can_update_child_nest

end module nest_manager