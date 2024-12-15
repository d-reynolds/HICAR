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
    !use data_structures
    use options_interface,  only : options_t
    use domain_interface,   only : domain_t
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
    public:: end_nest_context, start_nest_context

contains

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

    end subroutine start_nest_context

    subroutine end_nest_context(domain, options)
        implicit none
        type(domain_t), intent(in) :: domain
        type(options_t), intent(in) :: options

        if (options%physics%windtype == kITERATIVE_WINDS .or. &
            options%physics%windtype == kLINEAR_ITERATIVE_WINDS) then
                call finalize_iter_winds()
        end if
        
    end subroutine  end_nest_context

end module nest_manager