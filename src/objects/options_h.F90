module options_interface

    use icar_constants,             only : kMAX_STRING_LENGTH, kMAX_STORAGE_VARS
    use options_types,              only : general_options_type, output_options_type, domain_options_type, &
                                           forcing_options_type, restart_options_type,                     &
                                           physics_type, mp_options_type, lt_options_type, sfc_options_type, &
                                           adv_options_type, lsm_options_type, pbl_options_type, &
                                           cu_options_type, rad_options_type, wind_type, time_options_type

    implicit none

    private
    public :: options_t

    type :: options_t
        character(len=kMAX_STRING_LENGTH) :: comment

        ! master list of variables for different processes... not sure if this is the best place to store this information

        ! these are the variables that the advection code should process
        integer :: vars_to_advect(   kMAX_STORAGE_VARS ) = 0
        ! these are the variables that should be exchanged only
        integer :: vars_to_exch(   kMAX_STORAGE_VARS ) = 0
        ! these are the variables that need to be allocated for the model to run given the physics options requested
        integer :: vars_to_allocate( kMAX_STORAGE_VARS ) = 0
        ! these are the variables that need to be written and read from disk for a model restart run
        integer :: vars_for_restart( kMAX_STORAGE_VARS ) = 0


        type(general_options_type)      :: general

        type(domain_options_type)       :: domain

        type(forcing_options_type)      :: forcing

        type(restart_options_type)      :: restart

        type(output_options_type)       :: output

        ! defines which physics package to be used.
        type(physics_type)              :: physics

        type(wind_type)                 :: wind
        
        type(time_options_type)         :: time

        ! physics parameterization options
        type(mp_options_type)           :: mp

        type(lt_options_type)           :: lt

        type(adv_options_type)          :: adv

        type(lsm_options_type)          :: lsm

        type(cu_options_type)           :: cu

        type(rad_options_type)          :: rad
        
        type(pbl_options_type)          :: pbl

        type(sfc_options_type)          :: sfc

    contains

        procedure, public  :: init
        procedure, public  :: alloc_vars
        procedure, public  :: restart_vars
        procedure, public  :: advect_vars
        procedure, public  :: exch_vars
    end type


interface

    module subroutine init(this)
        implicit none
        class(options_t),   intent(inout)  :: this
    end subroutine

    module subroutine alloc_vars(this, input_vars, var_idx, error)
        implicit none
        class(options_t), intent(inout):: this
        integer, optional, intent(in)  :: input_vars(:)
        integer, optional, intent(in)  :: var_idx
        integer, optional, intent(out) :: error
    end subroutine

    module subroutine restart_vars(this, input_vars, var_idx, error)
        implicit none
        class(options_t), intent(inout):: this
        integer, optional, intent(in)  :: input_vars(:)
        integer, optional, intent(in)  :: var_idx
        integer, optional, intent(out) :: error
    end subroutine

    module subroutine advect_vars(this, input_vars, var_idx, error)
        implicit none
        class(options_t), intent(inout):: this
        integer, optional, intent(in)  :: input_vars(:)
        integer, optional, intent(in)  :: var_idx
        integer, optional, intent(out) :: error
    end subroutine

    module subroutine exch_vars(this, input_vars, var_idx, error)
        implicit none
        class(options_t), intent(inout):: this
        integer, optional, intent(in)  :: input_vars(:)
        integer, optional, intent(in)  :: var_idx
        integer, optional, intent(out) :: error
    end subroutine

end interface

end module
