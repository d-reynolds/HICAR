!> ----------------------------------------------------------------------------
!!  Model Initialization includes allocating memory for boundary and domain
!!      data structures.  It reads all of the options from the namelist
!!      file (or files).  It also reads in Lat/Lon and Terrain data.  This module
!!      also sets up geographic (and vertical) look uptables for the forcing data
!!      Finally, there is a driver routine to initialize all model physics packages
!!
!!   The module has been updated to allow arbitrary named variables
!!       this allows the use of e.g. ERAi, but still is not as flexible as it could be
!!
!!   The use of various python wrapper scripts in helpers/ makes it easy to add new
!!       datasets, and make them conform to the expectations of the current system.
!!      For now there are no plans to near term plans to substantially modify this.
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!! ----------------------------------------------------------------------------
module initialization
    !use data_structures
    use options_interface,  only : options_t
    use domain_interface,   only : domain_t
    use boundary_interface, only : boundary_t
    use microphysics,               only : mp_init
    use advection,                  only : adv_init
    use radiation,                  only : radiation_init
    use convection,                 only : init_convection
    use planetary_boundary_layer,   only : pbl_init
    use land_surface,               only : lsm_init
    use surface_layer,              only : sfc_init
    use io_routines,          only : io_read, io_write
    use mod_atm_utilities,          only : init_atm_utilities
    use wind,                       only : update_winds, init_winds
    use wind_iterative,             only : init_iter_winds
    use linear_theory_winds,        only : setup_linwinds

    use icar_constants!,             only : kITERATIVE_WINDS, kWIND_LINEAR
    use ioclient_interface,         only : ioclient_t
    use time_step,                  only : init_time_step
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
    public::init_model, init_physics, init_model_state

contains
    subroutine init_model(options,domain,boundary)
        implicit none
        type(options_t), intent(inout) :: options
        type(domain_t),  intent(inout) :: domain
        type(boundary_t),intent(inout) :: boundary ! forcing file for init conditions

        integer :: omp_get_max_threads, num_threads

#if defined(_OPENMP)
        num_threads = omp_get_max_threads()
#else
        num_threads = 1
#endif
        if (this_image()==1) call welcome_message()
        if (this_image()==1) flush(output_unit)

        if (this_image()==1) then
            write(*,*) "  Number of coarray image:",num_images()
            write(*,*) "  Max number of OpenMP Threads:",num_threads
        endif

        ! read in options file
        if (this_image()==1) write(*,*) "Initializing Options"
        if (this_image()==1) flush(output_unit)
        call options%init()
        
        if (this_image()==1) write(*,*) "Initializing Domain"
        if (this_image()==1) flush(output_unit)
        call domain%init(options)

        if (this_image()==1) write(*,*) "Initializing boundary condition data structure"
        if (this_image()==1) flush(output_unit)
        call boundary%init(options,domain%latitude%data_2d,domain%longitude%data_2d,domain%variables_to_force)

        if (this_image()==1) write(*,*) "Initializing time step helpers"
        if (this_image()==1) flush(output_unit)
        call init_time_step()

        if (this_image()==1) write(*,*) "Initializing atmospheric utilities"
        if (this_image()==1) flush(output_unit)
        ! initialize the atmospheric helper utilities
        call init_atm_utilities(options)

        if (options%physics%windtype==kITERATIVE_WINDS .or. options%physics%windtype==kLINEAR_ITERATIVE_WINDS) then
            !call init_iter_winds_old(domain)
            call init_iter_winds(domain)
        endif

        if (options%physics%windtype==kWIND_LINEAR .or. &
                 options%physics%windtype==kLINEAR_OBRIEN_WINDS .or. &
                 options%physics%windtype==kLINEAR_ITERATIVE_WINDS) then
            call setup_linwinds(domain, options, .False., options%parameters%advect_density)
        endif

        if (this_image()==1) write(*,'(/ A)') "Finished basic initialization"
        if (this_image()==1) write(*,'(A /)') "---------------------------------------"

    end subroutine init_model

    subroutine init_physics(options, domain, forcing)
        implicit none
        type(options_t), intent(inout) :: options
        type(domain_t),  intent(inout) :: domain
        type(boundary_t),intent(in)    :: forcing

        if (this_image()==1) write(*,*) "Init initial winds"
        if (this_image()==1) flush(output_unit)
        call init_winds(domain,options)

        if (this_image()==1) write(*,*) "Updating initial winds"
        if (this_image()==1) flush(output_unit)
        call update_winds(domain, forcing, options)

        ! initialize microphysics code (e.g. compute look up tables in Thompson et al)
        call mp_init(options) !this could easily be moved to init_model...
        if (this_image()==1) flush(output_unit)
        call init_convection(domain,options)
        if (this_image()==1) flush(output_unit)

        call pbl_init(domain,options)
        if (this_image()==1) flush(output_unit)

        call radiation_init(domain,options)
        if (this_image()==1) flush(output_unit)

        call lsm_init(domain,options)
        if (this_image()==1) flush(output_unit)

        call sfc_init(domain,options)
        if (this_image()==1) flush(output_unit)

        call adv_init(domain,options)
        if (this_image()==1) flush(output_unit)

    end subroutine init_physics

    subroutine init_model_state(options,domain,boundary)
        implicit none
        type(options_t), intent(inout) :: options
        type(domain_t),  intent(inout) :: domain
        type(boundary_t),intent(inout) :: boundary ! forcing file for init conditions

        if (this_image()==1) write(*,*) "Reading Initial conditions from boundary dataset"
        call domain%get_initial_conditions(boundary, options)

    end subroutine init_model_state

    subroutine welcome_message()
        implicit none

        write(*,*) ""
        write(*,*) "============================================================"
        write(*,*) "|                                                          |"
        write(*,*) "|  The Intermediate Complexity Atmospheric Research Model  |"
        write(*,*) "|                          (ICAR)                          |"
        write(*,*) "|                                                          |"
        write(*,*) "|   Developed at NCAR:                                     |"
        write(*,*) "|     The National Center for Atmospheric Research         |"
        write(*,*) "|     NCAR is sponsored by the National Science Foundation |"
        write(*,*) "|                                                          |"
        write(*,*) "|   Version: ",kVERSION_STRING,"                                         |"
        write(*,*) "|                                                          |"
        write(*,*) "============================================================"
        write(*,*) ""

    end subroutine welcome_message

end module
