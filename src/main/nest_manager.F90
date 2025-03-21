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
    use flow_object_interface, only : flow_obj_t
    use boundary_interface, only : boundary_t
    use ioclient_interface, only : ioclient_t
    use time_object,     only : Time_type
    use time_delta_object,     only : time_delta_t
    use initialization,             only : init_model_state
    use microphysics,               only : mp_init
    use advection,                  only : adv_init
    use radiation,                  only : radiation_init
    use convection,                 only : init_convection
    use planetary_boundary_layer,   only : pbl_init
    use land_surface,               only : lsm_init
    use surface_layer,              only : sfc_init
    use wind,                       only : init_winds
    use wind_iterative,             only : finalize_iter_winds
    use wind_iterative_old,         only : finalize_iter_winds_old
    use icar_constants
    use iso_fortran_env
    use mpi_f08


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

    public:: end_nest_context, start_nest_context, wake_nest
    public:: switch_nest_context
    public:: initialize_at_start
    public:: all_nests_not_done
    public:: nest_next_up
    public:: should_update_nests
    public:: can_update_child_nest


contains

    subroutine switch_nest_context(domain, options, nest_index)
        implicit none
        type(domain_t), intent(inout) :: domain(:)
        type(options_t), allocatable, intent(inout) :: options(:)
        integer, intent(in) :: nest_index

        if (current_nest == -1) then
            call start_nest_context(domain(nest_index), options(nest_index))
        endif

        if (current_nest /= nest_index) then
            call end_nest_context(domain(current_nest), options(current_nest))
            call start_nest_context(domain(nest_index), options(nest_index))
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

        call init_winds(domain,options)

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

    subroutine end_nest_context(domain, options)
        implicit none
        type(domain_t), intent(in) :: domain
        type(options_t), intent(in) :: options

        if (domain%nest_indx /= current_nest) then
            if (STD_OUT_PE) write(*,*) 'WARNING: Attempting to end nest context for nest ID: ', domain%nest_indx, ' but current nest ID is: ', current_nest
            return
        end if

        if (options%physics%windtype == kITERATIVE_WINDS .or. &
            options%physics%windtype == kLINEAR_ITERATIVE_WINDS) then
                call finalize_iter_winds_old()
        end if

        current_nest = -1
        
    end subroutine  end_nest_context

    subroutine wake_nest(options, domain, boundary, ioclient, nest_index)
        implicit none
        type(options_t), allocatable, intent(inout) :: options(:)
        type(domain_t), intent(inout) :: domain(:)
        type(boundary_t), intent(inout) :: boundary
        type(ioclient_t), intent(inout) :: ioclient
        integer, intent(in) :: nest_index


        ! If we have an active context and it is not us, end the context
        if ( (current_nest > 0) .and. (current_nest /= nest_index) ) then
            call end_nest_context(domain(current_nest), options(current_nest))
        endif

        call init_model_state(options(nest_index), domain(nest_index), boundary, ioclient)

        current_nest = nest_index

    end subroutine wake_nest

    !> ----------------------------------------------------------------------------
    !!  @function initialize_at_start(options)
    !!
    !!  @brief
    !!  Determines if we should intialize the nest context at the start of the simulation.
    !!
    !!  @param options
    !!  The options_t object for the nest context.
    !!  @return
    !!  Returns .true. if the nest context is initialized at the start of the simulation.
    !! ----------------------------------------------------------------------------
    function initialize_at_start(options, nest_indx) result(init)
        implicit none
        type(options_t), allocatable, intent(in) :: options(:)
        integer, intent(in) :: nest_indx
        logical :: init, start_time_match

        ! Initialize the nest context
        init = .false.
        start_time_match = .false.

        if (options(nest_indx)%general%parent_nest > 0) then
            start_time_match = (options(nest_indx)%general%start_time == options(options(nest_indx)%general%parent_nest)%general%start_time)
        endif

        if (options(nest_indx)%general%parent_nest == 0 .or. options(nest_indx)%restart%restart .or. start_time_match) then
            init = .true.
        end if

    end function initialize_at_start

    function all_nests_not_done(flow_objs) result(all_done)
        implicit none
        class(flow_obj_t), intent(in) :: flow_objs(:)
        logical :: all_done
        integer :: i

        all_done = ANY(flow_objs%ended .eqv. .False.)

    end function all_nests_not_done

    function nest_next_up(flow_objs, options, nest_indx) result(next)
        implicit none
        class(flow_obj_t), intent(in) :: flow_objs(:)
        type(options_t), intent(in) :: options
        integer, intent(in) :: nest_indx
        logical :: next

        next = .false.

        ! If we are a child nest, not at the end of our run time...
        if (flow_objs(nest_indx)%started .eqv. .False.) then
            ! ... and if we will be running on the next iteration...
            !safety check before indexing into flow_obj
            if (options%general%parent_nest > 0) then
                if (flow_objs(nest_indx)%sim_time - flow_objs(nest_indx)%small_time_delta <= flow_objs(options%general%parent_nest)%sim_time) then
                    next = .true.
                end if
            endif
        end if

    end function nest_next_up

    function should_update_nests(flow_objs, options, nest_indx) result(can_update)
        implicit none
        class(flow_obj_t), intent(in) :: flow_objs(:)
        type(options_t), intent(in) :: options
        integer, intent(in) :: nest_indx
        logical :: can_update

        can_update = .false.

        ! if we even have nests
        if (size(options%general%child_nests) > 0) then
            ! if we are at the end of our time step, or we have just ended
            if (flow_objs(nest_indx)%time_for_input() .or. flow_objs(nest_indx)%last_loop) then
                ! if all of our children are not done
                if (ANY(flow_objs(options%general%child_nests)%ended .eqv. .False.)) then
                    can_update = .true.
                end if
            end if
        endif

    end function should_update_nests

    function can_update_child_nest(flow_obj, child_flow_obj)
        implicit none
        class(flow_obj_t), intent(in) :: flow_obj
        class(flow_obj_t), intent(in) :: child_flow_obj
        logical :: can_update_child_nest

        can_update_child_nest = (flow_obj%sim_time >= child_flow_obj%sim_time - flow_obj%small_time_delta .and. .not.(child_flow_obj%ended))

    end function can_update_child_nest

end module nest_manager