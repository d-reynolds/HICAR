module variable_interface
    use icar_constants,          only : kMAX_DIM_LENGTH, kMAX_STRING_LENGTH, kMAX_NAME_LENGTH, kINTEGER, kREAL, kDOUBLE, STD_OUT_PE, kVARS
    use grid_interface,          only : grid_t
    use meta_data_interface,     only : meta_data_t, root_var_t
    implicit none

    ! defines a variable type that can store data and attributes
    ! have to think about how to handle multiple variable types (int, 2d, etc)
    ! could add multiple "local" variables or create multiple variable types...
    type, extends(root_var_t) :: variable_t
        real, allocatable :: data_4d(:,:,:,:)! => null()
        real, allocatable :: data_3d(:,:,:)! => null()
        real, allocatable :: data_2d(:,:) !  => null()
        real, allocatable :: data_1d(:) !  => null()
        !real(kind=real64), pointer :: data_2dd(:,:)! => null()
        integer, allocatable ::           data_2di(:,:)! => null()
        real, allocatable :: dqdt_3d(:,:,:) !=> null()   ! Note these have to be pointers so they get referenced when variable_t is passed around(?)
        real, allocatable :: dqdt_2d(:,:)  ! => null()   ! Note these have to be pointers so they get referenced when variable_t is passed around(?)

        logical                         :: force_boundaries = .True.
        logical                         :: computed = .False.
        logical                         :: forcing_var = .False.

        type(grid_t)                                :: grid

    contains
        procedure, public  :: init_grid
        procedure, public  :: init_dims
        procedure, private :: set_from_metadata
        generic,   public  :: initialize => init_grid
        generic,   public  :: initialize => init_dims

    ! inherited from meta_data
    !     procedure, public : add_attribute

    end type

    interface assignment(=)
        module procedure assign_variable
    end interface

    interface

        module subroutine init_grid(this, var_idx, grid, forcing_var)
            implicit none
            class(variable_t),  intent(inout) :: this
            integer,            intent(in)    :: var_idx
            type(grid_t),       intent(in)    :: grid
            logical,            intent(in), optional :: forcing_var

        end subroutine

        module subroutine init_dims(this, var_idx, dims, forcing_var)
            implicit none
            class(variable_t),  intent(inout) :: this
            integer,            intent(in)    :: var_idx
            integer,            intent(in)    :: dims(:)
            logical,            intent(in), optional :: forcing_var
        end subroutine

        module subroutine assign_variable(dest, src)
            implicit none
            class(variable_t), intent(out) :: dest
            class(variable_t), intent(in)  :: src
        end subroutine

        module subroutine set_from_metadata(this, var_id)
            implicit none
            class(variable_t), intent(inout) :: this
            integer, intent(in) :: var_id
        end subroutine
    
    end interface

end module
