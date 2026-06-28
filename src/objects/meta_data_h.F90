module meta_data_interface
    use icar_constants, only: kMAX_NAME_LENGTH, kMAX_ATTR_LENGTH, kMAX_DIM_LENGTH, kREAL, kUNSET_REAL

    implicit none

    private
    public :: meta_data_t
    public :: attribute_t
    public :: root_var_t

    type attribute_t
        character(len=kMAX_NAME_LENGTH) :: name
        character(len=kMAX_ATTR_LENGTH) :: value
    end type

    type root_var_t
        integer :: id = -1
        integer :: dtype = kREAL
        logical :: one_d = .False.
        logical :: three_d = .False.
        logical :: two_d = .False.
        logical :: four_d = .False.

        integer :: xstag = 0
        integer :: ystag = 0
        real    :: maxval = kUNSET_REAL   
        real    :: minval = kUNSET_REAL
        integer, allocatable :: dim_len(:)
        integer, allocatable :: global_dim_len(:)

        ! note these are used for netcdf output
        integer, allocatable    :: dim_ids(:)
        integer                 :: file_var_id = -1

    end type root_var_t

    type, extends(root_var_t) :: meta_data_t
    !   private
        character(len=kMAX_NAME_LENGTH) :: name = ""
        character(len=kMAX_DIM_LENGTH), allocatable :: dimensions(:)

        integer :: n_attrs = 0
        logical :: unlimited_dim = .False.
        logical :: output = .True.   ! Whether this variable can appear in output files
        ! Per-nest static variable: geometry / grid metric / static surface descriptor.
        ! The child reads/derives these from its own input file in read_core_variables /
        ! setup_sleve, so they must NOT be overwritten by the parent->child nest init
        ! transfer at child wake-up. Checked in unpack_init_vars_{2d,3d}.
        logical :: static_data = .False.

        type(attribute_t), allocatable :: attributes(:)
    contains

        procedure, public  :: add_attribute
        ! procedure, public  :: get_attribute
    end type

    interface

        module subroutine add_attribute(this, input_name, input_value)
            implicit none
            class(meta_data_t), intent(inout) :: this
            character(len=*),   intent(in)    :: input_name
            character(len=*),   intent(in)    :: input_value
        end subroutine

        ! module subroutine get_attribute(this, name, value)
        !     implicit none
        !     class(output_t),   intent(inout)  :: this
        !     character(len=kMAX_ATTR_LENGTH), intent(in) :: name
        !     character(len=kMAX_ATTR_LENGTH), intent(value) :: value
        ! end subroutine

    end interface
end module
