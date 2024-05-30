module halo_interface
        
    use grid_interface,           only : grid_t
    use variable_interface,       only : variable_t
    use variable_dict_interface,  only : var_dict_t
    use timer_interface,          only : timer_t

    use mpi_f08
    implicit none

    private
    public :: halo_t

    type halo_t

        integer           :: halo_size

        type(grid_t)      :: grid

        type(MPI_Win)     :: north_in_win
        type(MPI_Win)     :: south_in_win
        type(MPI_Win)     :: east_in_win
        type(MPI_Win)     :: west_in_win

        type(MPI_win)     :: north_3d_win
        type(MPI_win)     :: south_3d_win
        type(MPI_win)     :: east_3d_win
        type(MPI_win)     :: west_3d_win

        type(MPI_win)     :: north_2d_win
        type(MPI_win)     :: south_2d_win
        type(MPI_win)     :: east_2d_win
        type(MPI_win)     :: west_2d_win

        type(MPI_Datatype) :: N_3d_win_halo_type
        type(MPI_Datatype) :: S_3d_win_halo_type
        type(MPI_Datatype) :: EW_3d_win_halo_type

        type(MPI_Datatype) :: NS_2d_win_halo_type
        type(MPI_Datatype) :: EW_2d_win_halo_type

        type(MPI_Group)    :: north_neighbor_grp, south_neighbor_grp, east_neighbor_grp, west_neighbor_grp

        type(MPI_comm) :: neighbor_comm
#ifdef CRAY_PE
    	real, allocatable :: south_batch_in_3d(:,:,:,:)[:]
    	real, allocatable :: north_batch_in_3d(:,:,:,:)[:]
    	real, allocatable :: west_batch_in_3d(:,:,:,:)[:]
    	real, allocatable :: east_batch_in_3d(:,:,:,:)[:]

 	    real, allocatable :: south_batch_in_2d(:,:,:)[:]
 	    real, allocatable :: north_batch_in_2d(:,:,:)[:]
        real, allocatable :: west_batch_in_2d(:,:,:)[:]
   	    real, allocatable :: east_batch_in_2d(:,:,:)[:]
#else
        real, contiguous, pointer :: south_batch_in_3d(:,:,:,:)
        real, contiguous, pointer :: north_batch_in_3d(:,:,:,:)
        real, contiguous, pointer :: west_batch_in_3d(:,:,:,:)
        real, contiguous, pointer :: east_batch_in_3d(:,:,:,:)

        real, pointer :: south_batch_in_2d(:,:,:)
        real, pointer :: north_batch_in_2d(:,:,:)
        real, pointer :: west_batch_in_2d(:,:,:)
        real, pointer :: east_batch_in_2d(:,:,:)
#endif

#ifdef CRAY_PE
        real, allocatable :: south_in_3d(:,:,:)[:]
        real, allocatable :: north_in_3d(:,:,:)[:]
        real, allocatable :: west_in_3d(:,:,:)[:]
        real, allocatable :: east_in_3d(:,:,:)[:]
#else
        real, pointer     :: south_in_3d(:,:,:)
        real, pointer     :: north_in_3d(:,:,:)
        real, pointer     :: west_in_3d(:,:,:)
        real, pointer     :: east_in_3d(:,:,:)
#endif
   	 real, allocatable :: north_buffer_3d(:,:,:,:)
   	 real, allocatable :: south_buffer_3d(:,:,:,:)
     real, allocatable :: east_buffer_3d(:,:,:,:)
   	 real, allocatable :: west_buffer_3d(:,:,:,:)

   	  real, allocatable :: north_buffer_2d(:,:,:)
  	  real, allocatable :: south_buffer_2d(:,:,:)
  	  real, allocatable :: east_buffer_2d(:,:,:)
   	  real, allocatable :: west_buffer_2d(:,:,:)

    	  ! Neighboring images of this image
    	  integer, allocatable :: neighbors(:)
    	  integer, allocatable :: corner_neighbors(:)

  	  integer :: north_neighbor, south_neighbor, east_neighbor, west_neighbor
    	  integer :: northwest_neighbor, southwest_neighbor, northeast_neighbor, southeast_neighbor

  	  logical :: north_boundary = .True.
  	  logical :: south_boundary = .True.
   	  logical :: east_boundary = .True.
    	  logical :: west_boundary = .True.

  	  ! store the start (s) and end (e) for the i,j,k dimensions
  	  integer ::  ids,ide, jds,jde, kds,kde, & ! for the entire model domain    (d)
                ims,ime, jms,jme, kms,kme, & ! for the memory in these arrays (m)
                its,ite, jts,jte, kts,kte ! for the data tile to process   (t)

    contains
        procedure, public :: init
        procedure, public :: exch_var
        procedure, public :: batch_exch
        procedure, public :: halo_3d_send_batch
        procedure, public :: halo_3d_retrieve_batch

        procedure :: put_north
        procedure :: put_south
        procedure :: put_west
        procedure :: put_east
        procedure :: retrieve_north_halo
        procedure :: retrieve_south_halo
        procedure :: retrieve_west_halo
        procedure :: retrieve_east_halo

        procedure :: put_northwest
        procedure :: put_northeast
        procedure :: put_southwest
        procedure :: put_southeast
        procedure :: retrieve_northwest_halo
        procedure :: retrieve_northeast_halo
        procedure :: retrieve_southwest_halo
        procedure :: retrieve_southeast_halo

        procedure :: halo_2d_send_batch
        procedure :: halo_2d_retrieve_batch
        !procedure :: halo_3d_send_batch
        !procedure :: halo_3d_retrieve_batch

    end type

interface

    module subroutine init(this, exch_vars, adv_vars, grid, comms)
        implicit none
        class(halo_t), intent(inout) :: this
        type(var_dict_t), intent(inout) :: exch_vars, adv_vars
        type(grid_t), intent(in) :: grid
        type(MPI_comm), intent(inout) :: comms
    end subroutine init

    module subroutine exch_var(this, var, do_dqdt, corners)
        implicit none
        class(halo_t),     intent(inout) :: this
        type(variable_t),  intent(inout) :: var
        logical,          intent(in), optional :: do_dqdt, corners
    end subroutine exch_var


    module subroutine batch_exch(this, exch_vars, adv_vars, two_d, three_d, exch_var_only)
        implicit none
        class(halo_t), intent(inout) :: this
        class(var_dict_t), intent(inout) :: exch_vars, adv_vars
        logical, optional, intent(in) :: two_d,three_d,exch_var_only
    end subroutine

    module subroutine halo_3d_send_batch(this, exch_vars, adv_vars,exch_var_only)
        implicit none
        class(halo_t), intent(inout) :: this
        class(var_dict_t), intent(inout) :: exch_vars, adv_vars
        logical, optional, intent(in) :: exch_var_only
    end subroutine halo_3d_send_batch

    module subroutine halo_3d_retrieve_batch(this, exch_vars, adv_vars,exch_var_only, wait_timer)
        implicit none
        class(halo_t), intent(inout) :: this
        class(var_dict_t), intent(inout) :: exch_vars, adv_vars
        logical, optional, intent(in) :: exch_var_only
        type(timer_t), optional,     intent(inout)   :: wait_timer

    end subroutine halo_3d_retrieve_batch

    module subroutine halo_2d_send_batch(this, exch_vars, adv_vars)
        implicit none
        class(halo_t), intent(inout) :: this
        class(var_dict_t), intent(inout) :: exch_vars, adv_vars
    end subroutine halo_2d_send_batch

    module subroutine halo_2d_retrieve_batch(this, exch_vars, adv_vars)
        implicit none
        class(halo_t), intent(inout) :: this
        class(var_dict_t), intent(inout) :: exch_vars, adv_vars
    end subroutine halo_2d_retrieve_batch


    module subroutine put_north(this,var,do_dqdt)
        implicit none
        class(halo_t), intent(inout) :: this
        class(variable_t), intent(in) :: var
        logical, optional, intent(in) :: do_dqdt
    end subroutine

    module subroutine put_south(this,var,do_dqdt)
        implicit none
        class(halo_t), intent(inout) :: this
        class(variable_t), intent(in) :: var
        logical, optional, intent(in) :: do_dqdt
    end subroutine

    module subroutine retrieve_north_halo(this,var,do_dqdt)
        implicit none
        class(halo_t), intent(inout) :: this
        class(variable_t), intent(in) :: var
        logical, optional, intent(in) :: do_dqdt
    end subroutine

    module subroutine retrieve_south_halo(this,var,do_dqdt)
        implicit none
        class(halo_t), intent(inout) :: this
        class(variable_t), intent(in) :: var
        logical, optional, intent(in) :: do_dqdt
    end subroutine

    module subroutine put_east(this,var,do_dqdt)
        implicit none
        class(halo_t), intent(inout) :: this
        class(variable_t), intent(in) :: var
        logical, optional, intent(in) :: do_dqdt
    end subroutine

    module subroutine put_west(this,var,do_dqdt)
        implicit none
        class(halo_t), intent(inout) :: this
        class(variable_t), intent(in) :: var
        logical, optional, intent(in) :: do_dqdt

    end subroutine

    module subroutine retrieve_east_halo(this,var,do_dqdt)
        implicit none
        class(halo_t), intent(inout) :: this
        class(variable_t), intent(in) :: var
        logical, optional, intent(in) :: do_dqdt

    end subroutine

    module subroutine retrieve_west_halo(this,var,do_dqdt)
        implicit none
        class(halo_t), intent(inout) :: this
        class(variable_t), intent(in) :: var
        logical, optional, intent(in) :: do_dqdt

    end subroutine

  module subroutine put_northwest(this,var,do_dqdt)
      implicit none 
      class(halo_t), intent(inout) :: this
      class(variable_t), intent(in) :: var
      logical, optional, intent(in) :: do_dqdt
  end subroutine

  module subroutine put_northeast(this,var,do_dqdt)
      implicit none 
      class(halo_t), intent(inout) :: this
      class(variable_t), intent(in) :: var
      logical, optional, intent(in) :: do_dqdt
  end subroutine

  module subroutine put_southeast(this,var,do_dqdt)
      implicit none 
      class(halo_t), intent(inout) :: this
      class(variable_t), intent(in) :: var
      logical, optional, intent(in) :: do_dqdt
  end subroutine

  module subroutine put_southwest(this,var,do_dqdt)
      implicit none 
      class(halo_t), intent(inout) :: this
      class(variable_t), intent(in) :: var
      logical, optional, intent(in) :: do_dqdt
  end subroutine

  module subroutine retrieve_northwest_halo(this,var,do_dqdt)
      implicit none 
      class(halo_t), intent(inout) :: this
      class(variable_t), intent(in) :: var
      logical, optional, intent(in) :: do_dqdt
  end subroutine

  module subroutine retrieve_northeast_halo(this,var,do_dqdt)
      implicit none 
      class(halo_t), intent(inout) :: this
      class(variable_t), intent(in) :: var
      logical, optional, intent(in) :: do_dqdt
  end subroutine

  module subroutine retrieve_southeast_halo(this,var,do_dqdt)
      implicit none 
      class(halo_t), intent(inout) :: this
      class(variable_t), intent(in) :: var
      logical, optional, intent(in) :: do_dqdt
  end subroutine

  module subroutine retrieve_southwest_halo(this,var,do_dqdt)
      implicit none 
      class(halo_t), intent(inout) :: this
      class(variable_t), intent(in) :: var
      logical, optional, intent(in) :: do_dqdt
  end subroutine

end interface
end module
