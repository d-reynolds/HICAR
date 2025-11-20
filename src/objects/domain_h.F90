module domain_interface
  use, intrinsic :: iso_c_binding, only : C_DOUBLE_COMPLEX
  use mpi_f08
  use options_interface,        only : options_t
  use boundary_interface,       only : boundary_t
  use grid_interface,           only : grid_t
  use variable_interface,       only : variable_t
  use data_structures,          only : interpolable_type, tendencies_type, index_type
  use halo_interface,           only : halo_t
  use timer_interface,          only : timer_t
  use flow_object_interface,    only : flow_obj_t
  use icar_constants,               only : kMAX_STORAGE_VARS, kVARS, kMAX_NAME_LENGTH
  implicit none

  private
  public :: domain_t

  type , extends(flow_obj_t) :: domain_t
    type(grid_t)         :: grid,   grid8w,  u_grid,   v_grid
    type(grid_t)         :: column_grid
    type(grid_t)         :: global_grid_2d, global_grid, global_grid8w
    type(grid_t)         :: neighbor_grid_2d, neighbor_grid, neighbor_grid8w
    type(grid_t)         :: grid2d, u_grid2d, v_grid2d
    type(grid_t)         :: grid_monthly, grid_soil
    type(grid_t)         :: grid_snow, grid_snowsoil
    type(grid_t)         :: grid_soilcomp, grid_gecros, grid_croptype
    type(grid_t)         :: grid_hlm, grid_Sx !! MJ added
    type(grid_t)         :: grid_lake , grid_lake_soisno, grid_lake_soi, grid_lake_soisno_1
    type(halo_t)         :: halo

    ! note that not all variables are allocated at runtime, physics packages must request a variable be created
    ! though variables considered "required" are requested by the domain object itself (e.g. terrain)
    ! core model species to be advected

    ! these data are stored on the domain wide grid even if this process is only looking at a subgrid
    ! these variables are necessary with linear winds, especially with spatially variable dz, to compute the LUT

    type(tendencies_type) :: tend

    type(index_type) :: vars_to_out(kMAX_STORAGE_VARS)
    
    ! Array listing variables to advect with pointers to local data
    type(index_type), allocatable :: adv_vars(:), exch_vars(:)
    integer :: n_adv_2d, n_adv_3d, n_exch_2d, n_exch_3d
    
    type(interpolable_type) :: geo
    type(interpolable_type) :: geo_agl
    type(interpolable_type) :: geo_u
    type(interpolable_type) :: geo_v

    real :: smooth_height, dx
    integer :: nsmooth

    ! complex(C_DOUBLE_COMPLEX),  allocatable :: terrain_frequency(:,:) ! FFT(terrain)

    type(variable_t), allocatable :: forcing_hi(:)

    type(variable_t), allocatable :: vars_1d(:)
    type(variable_t), allocatable :: vars_2d(:)
    type(variable_t), allocatable :: vars_3d(:)
    type(variable_t), allocatable :: vars_4d(:)

    type(index_type) :: var_indx(kMAX_STORAGE_VARS), forcing_var_indx(kMAX_STORAGE_VARS)
    
    ! MPI communicator object for doing parallel communications among domain objects
    type(MPI_Comm), public :: compute_comms

    ! timers used to track the time spent doing various operations
    type(timer_t) :: initialization_timer, total_timer, input_timer, &
                        output_timer, physics_timer, wind_timer, mp_timer, &
                        adv_timer, rad_timer, lsm_timer, pbl_timer, exch_timer, &
                        send_timer, ret_timer, wait_timer, forcing_timer, diagnostic_timer, wind_bal_timer, &
                        flux_timer, flux_corr_timer, sum_timer, adv_wind_timer, cpu_gpu_timer

    ! contains the size of the domain (or the local tile?)
    integer :: nx, ny, nz, nx_global, ny_global
    integer :: ximg, ximages, yimg, yimages
    logical :: north_boundary = .True.
    logical :: south_boundary = .True.
    logical :: east_boundary = .True.
    logical :: west_boundary = .True.

    ! store the start (s) and end (e) for the i,j,k dimensions
    integer ::  ids,ide, jds,jde, kds,kde, & ! for the entire model domain    (d)
                ims,ime, jms,jme, kms,kme, & ! for the memory in these arrays (m)
                its,ite, jts,jte, kts,kte, & ! for the data tile to process   (t)
                ihs,ihe, jhs,jhe, khs,khe    ! for the neighborhood arrays for non-local calculations (h)

    integer :: neighborhood_max ! The maximum neighborhood radius in indices
    

  contains
    procedure :: init => init_domain
    procedure :: release
    procedure :: enforce_limits

    procedure :: batch_exch
    procedure :: halo_3d_send
    procedure :: halo_3d_retrieve
    procedure :: halo_2d_send
    procedure :: halo_2d_retrieve

    procedure :: get_initial_conditions
    procedure :: diagnostic_update
    procedure :: interpolate_forcing
    procedure :: update_delta_fields
    procedure :: apply_forcing

    procedure :: update_device
    procedure :: update_host

  end type

  integer, parameter :: space_dimension=3

  interface

    ! Set default component values
    module subroutine init_domain(this, options, nest_indx)
        implicit none
        class(domain_t), intent(inout) :: this
        type(options_t), intent(inout) :: options
        integer,         intent(in)    :: nest_indx
    end subroutine init_domain

    ! finalize domain object, freeing halo mpi windows
    module subroutine release(this)
        implicit none
        class(domain_t), intent(inout) :: this
    end subroutine

    module subroutine batch_exch(this, two_d, exch_only)
        implicit none
        class(domain_t), intent(inout) :: this
        logical, optional,   intent(in) :: two_d, exch_only
  end subroutine

    module subroutine halo_3d_send(this, exch_only)
        implicit none
        class(domain_t), intent(inout) :: this
        logical, optional,   intent(in) :: exch_only
    end subroutine
    
    module subroutine halo_3d_retrieve(this, exch_only)
      implicit none
      class(domain_t), intent(inout) :: this
      logical, optional,   intent(in) :: exch_only
  end subroutine

  module subroutine halo_2d_send(this)
    implicit none
    class(domain_t), intent(inout) :: this
  end subroutine

  module subroutine halo_2d_retrieve(this)
    implicit none
    class(domain_t), intent(inout) :: this
  end subroutine

    ! read initial atmospheric conditions from forcing data
    module subroutine get_initial_conditions(this, forcing, options)
        implicit none
        class(domain_t),  intent(inout) :: this
        type(boundary_t), intent(inout) :: forcing
        type(options_t),  intent(in)    :: options
    end subroutine

    module subroutine diagnostic_update(this, forcing_update)
      implicit none
      class(domain_t),  intent(inout)   :: this
      logical, intent(in), optional    :: forcing_update
    end subroutine

    module subroutine interpolate_forcing(this, forcing, update)
        implicit none
        class(domain_t),  intent(inout) :: this
        type(boundary_t), intent(inout) :: forcing
        logical,          intent(in),   optional :: update
    end subroutine

    ! Exchange subdomain boundary information
    !module subroutine halo_exchange(this, two_d, exch_only)
    !    implicit none
    !    class(domain_t), intent(inout) :: this
    !    logical, optional,   intent(in) :: two_d, exch_only
    !end subroutine

    ! Make sure no hydrometeors are getting below 0
    module subroutine enforce_limits(this,update_in)
        implicit none
        class(domain_t), intent(inout) :: this
        logical, optional, intent(in)  :: update_in
    end subroutine

    module subroutine update_delta_fields(this)
        implicit none
        class(domain_t),    intent(inout) :: this
    end subroutine

    module subroutine apply_forcing(this, options, dt)
        implicit none
        class(domain_t),    intent(inout) :: this
        type(options_t), intent(in)       :: options
        real, intent(in)                  :: dt
    end subroutine

    module subroutine update_device(this)
        implicit none
        class(domain_t), intent(inout) :: this
    end subroutine

    module subroutine update_host(this)
        implicit none
        class(domain_t), intent(inout) :: this
    end subroutine

  end interface

end module