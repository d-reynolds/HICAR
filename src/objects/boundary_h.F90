module boundary_interface

    use mpi_f08
    use netcdf
    use options_interface,        only : options_t
    use options_types,            only : dim_arrays_type, forcing_options_type
    use variable_dict_interface,  only : var_dict_t
    use variable_interface,       only : variable_t
    use time_object,              only : Time_type
    use time_delta_object,        only : time_delta_t
    use data_structures,          only : interpolable_type, dim_arrays_type
    use icar_constants,           only : kMAX_NAME_LENGTH, kMAX_STRING_LENGTH, kMAX_FILE_LENGTH
    use grid_interface,           only : grid_t
    use flow_object_interface,     only : flow_obj_t
    implicit none

    private
    public :: boundary_t

    ! ------------------------------------------------
    ! boundary conditions type, must be linearizable so we can remove low res linear wind field
    ! ------------------------------------------------
    type :: boundary_t

        !   manage file pointer and position in file for boundary conditions
        character(len=kMAX_NAME_LENGTH) :: firstfile
        integer :: firststep
        integer :: its
        integer :: ite
        integer :: kts
        integer :: kte
        integer :: jts
        integer :: jte
        integer :: ids
        integer :: ide
        integer :: jds
        integer :: jde


        type(Time_type)                   :: current_time   ! the date/time of the forcing data in memory
        type(time_delta_t)                :: forcing_dt     ! the time step in between two forcing steps
        character(len=kMAX_STRING_LENGTH) :: time_var       ! the name of the input time variable [optional]

        type(var_dict_t)                  :: variables      ! a dictionary with all forcing data

        ! boundary data coordinate system
        real, dimension(:,:),   allocatable :: lat, lon, ulat, ulon, vlat, vlon
        real, dimension(:,:,:), allocatable :: z            ! 3D z elevation data on the forcing grid
        
        logical :: z_is_set = .false.

        type(interpolable_type) :: geo
        type(interpolable_type) :: geo_agl
        type(interpolable_type) :: geo_u
        type(interpolable_type) :: geo_v
        type(interpolable_type) :: original_geo
        type(interpolable_type) :: original_geo_u
        type(interpolable_type) :: original_geo_v
        type(interpolable_type) :: mass_to_u, mass_to_v

    contains

        procedure :: init

        ! procedure :: find_start_time
        procedure :: init_local
        procedure :: init_local_asnest
        procedure :: update_computed_vars
        procedure :: interpolate_original_levels
        procedure :: setup_z
    end type

    interface

    ! Set default component values
    module subroutine init(this, options, domain_lat, domain_lon, parent_options)
        class(boundary_t),    intent(inout) :: this
        type(options_t),      intent(inout) :: options
        real, dimension(:,:), intent(in)    :: domain_lat
        real, dimension(:,:), intent(in)    :: domain_lon
        type(options_t), optional, intent(in)    :: parent_options
    end subroutine

    module subroutine init_local(this, options, var_list, dim_list, start_time, domain_lat, domain_lon)
        implicit none
        class(boundary_t),               intent(inout)  :: this
        type(forcing_options_type),      intent(inout)  :: options
        character(len=kMAX_NAME_LENGTH), intent(in)     :: var_list (:)
        type(dim_arrays_type),           intent(in)     :: dim_list (:)
        type(Time_type),                 intent(in)     :: start_time
        real, dimension(:,:),            intent(in)     :: domain_lat
        real, dimension(:,:),            intent(in)     :: domain_lon
    end subroutine

    module subroutine init_local_asnest(this, var_list, dim_list, domain_lat, domain_lon, parent_options)
        class(boundary_t),               intent(inout)  :: this
        character(len=kMAX_NAME_LENGTH), intent(in)     :: var_list (:)
        type(dim_arrays_type),           intent(in)     :: dim_list (:)
        real, dimension(:,:),            intent(in)     :: domain_lat
        real, dimension(:,:),            intent(in)     :: domain_lon
        type(options_t),                 intent(in)     :: parent_options
    end subroutine
    
    module subroutine setup_z(this, options)
        implicit none
        class(boundary_t),   intent(inout)   :: this
        type(options_t),     intent(in)      :: options
    end subroutine

    module subroutine update_computed_vars(this, options)
        implicit none
        class(boundary_t),   intent(inout)   :: this
        type(options_t),     intent(in)      :: options
    end subroutine
    
    module subroutine interpolate_original_levels(this, options)
        implicit none
        class(boundary_t),   intent(inout)   :: this
        type(options_t),     intent(in)      :: options
    end subroutine
    
  end interface

end module
