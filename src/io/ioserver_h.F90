!>----------------------------------------------------------
!!  Define the interface for the ioserver
!!
!!  The I/O server serves as a parallelized interface between the main program
!!  and the I/O routines. The I/O server handles the buffering and exchange
!!  of data between the processes of the main program (child processes)
!!  and the I/O processes (parent processes)
!! 
!!  @author
!!  Dylan Reynolds (dylan.reynolds@slf.ch)
!!
!!----------------------------------------------------------
module ioserver_interface
  use mpi_f08
  use netcdf
  use icar_constants
  use reader_interface,   only : reader_t
  use output_interface,   only : output_t
  use options_interface,  only : options_t
  use flow_object_interface, only : flow_obj_t, comp_arr_t
  
  implicit none

  private
  public :: ioserver_t

  !>----------------------------------------------------------
  !! Output type definition
  !!
  !!----------------------------------------------------------

  type buffer_3d_t
     real, pointer, dimension(:,:,:,:) :: buff
  end type buffer_3d_t

  type buffer_2d_t
     real, pointer, dimension(:,:,:) :: buff
  end type buffer_2d_t

  type rst_pack_t
     real, allocatable :: buff_3d(:,:,:,:)
     real, allocatable :: buff_2d(:,:,:)
  end type rst_pack_t

  ! Ragged list of (i,j) cells assigned to one peer ioserver, used to build
  ! nest-transfer MPI datatypes. Populated by build_nest_geometry and consumed
  ! by build_indexed_type_from_cells.
  type :: nest_cell_list_t
     integer :: n = 0
     integer, allocatable :: i(:), j(:)
  end type nest_cell_list_t

  ! Per-child nest geometry: which (i,j) cells go to/from each peer ioserver.
  ! This information is invariant in the k axis and in n_vars, so it is shared
  ! across the 3D-per-step, 2D-init, and 3D-init datatype variants. Built once
  ! (lazily) per child and freed after all three datatype sets are committed.
  type :: nest_geometry_t
     logical :: built = .false.
     type(nest_cell_list_t), allocatable :: send_cells(:)
     type(nest_cell_list_t), allocatable :: buff_cells(:)
  end type nest_geometry_t

  type , extends(flow_obj_t) :: ioserver_t

        ! all components are private and should only be modified through procedures
        private

        ! Store the variables to be written
        ! Note n_variables may be smaller then size(variables) so that it doesn't
        ! have to keep reallocating variables whenever something is added or removed
        integer, public :: n_input_variables, n_output_variables, n_servers, n_children, n_child_ioservers
        type(MPI_Comm), public :: client_comms, IO_comms
        logical, public :: files_to_read
        logical         :: first_write = .true.
        logical, allocatable :: nest_types_initialized(:)

        ! time variable , publicis stored outside of the variable list... probably need to think about that some
        type(output_t) :: outputer
        type(reader_t) :: reader
        
        real, dimension(:,:,:,:), allocatable     :: gather_buffer, forcing_buffer
        ! Init-only nest transfer buffers (2D + extra 3D restart vars)
        real, dimension(:,:,:),   allocatable     :: gather_buffer_2d
        real, dimension(:,:,:,:), allocatable     :: gather_buffer_3d_init
        real, dimension(:,:,:),   allocatable     :: forcing_buffer_2d
        real, dimension(:,:,:,:), allocatable     :: forcing_buffer_3d_init

        type(buffer_3d_t), allocatable :: write_buffer_3d(:), read_buffer(:), child_gather_buffers(:)
        type(buffer_3d_t), allocatable :: child_gather_buffers_3d_init(:)
        type(buffer_2d_t), allocatable :: write_buffer_2d(:)
        type(buffer_2d_t), allocatable :: child_gather_buffers_2d(:)

        ! write_win_3d/2d retained for the (stubbed) restart-read path. All
        ! other RMA windows have been replaced with two-sided Isend/Irecv.
        type(MPI_Win) :: write_win_3d, write_win_2d

        ! These MPI datatypes describe the access patterns between the IO read/write buffers and
        ! the child read/write buffers
        type(MPI_Datatype), allocatable, dimension(:) :: rst_types_3d, rst_types_2d, child_rst_types_3d, child_rst_types_2d
        type(MPI_Datatype), allocatable, dimension(:,:) :: send_nest_types, buffer_nest_types
        type(MPI_Datatype), allocatable, dimension(:,:) :: send_nest_types_2d, buffer_nest_types_2d
        type(MPI_Datatype), allocatable, dimension(:,:) :: send_nest_types_3d_init, buffer_nest_types_3d_init
        logical, allocatable :: nest_types_2d_initialized(:), nest_types_3d_init_initialized(:)

        ! Per-child cached cell geometry shared across the three setup_nest_types_* variants
        type(nest_geometry_t), allocatable :: nest_geometry(:)
        
        ! coarray-indices of child io processes, indexed according to the COMPUTE_TEAM
        integer, allocatable :: children_ranks(:)
        
        ! coarray-index of parent io process, indexed according to the IO_TEAM
        integer :: parent_id = 0
            
        ! indices of process extent. For child processes, this is just the tile extents, 
        ! except for edges, where we go to the domain extents
        !Convention is:
        ! (i/k/j) -- x, z, or y index
        ! (s/e) -- start or end index
        ! (r/w) -- index for the read or write buffer
        ! c -- defines this index as the index of a particular child, 
        ! where the index in this array corresponds to the child image at the same index in the children array

        integer, allocatable, dimension(:) :: isrc, ierc, ksrc, kerc, jsrc, jerc, iswc, iewc, kswc, kewc, jswc, jewc, isrec, ierec, jsrec, jerec
        
        integer :: i_s_w, i_e_w, k_s_w, k_e_w, j_s_w, j_e_w, i_s_r, i_e_r, k_s_r, k_e_r, k_e_f, j_s_r, j_e_r, i_s_re, i_e_re, j_s_re, j_e_re, n_restart
        integer :: nx_w, nz_w, ny_w, n_w_3d, n_w_2d, nx_r, nz_r, ny_r, n_r, n_f, n_f_2d, n_f_3d_init, nx_re, ny_re
        integer, public :: nz_init_3d = 0
        logical :: initial_nest_done = .false.
        integer         :: ide, kde, jde
        integer :: restart_counter = 0
        integer :: output_counter = 0
        integer :: frames_per_outfile, restart_count
        
        !the indices of the output buffer corresponding to the restart vars
        integer, public, allocatable :: out_var_indices(:), rst_var_indices(:)
        !the names of the restart vars, indexed the same as the above array
        character(len=kMAX_NAME_LENGTH), public, allocatable :: rst_var_names(:)

        ! Classification counts for two-pass buffer packing
        integer :: n_out_3d = 0, n_out_2d = 0
        integer :: n_rst_only_3d = 0, n_rst_only_2d = 0
        integer, allocatable :: out_ordered_indices(:)
        integer, allocatable :: rst_only_ordered_indices(:)

        ! The filename of the netcdf file to write
        character(len=kMAX_FILE_LENGTH) :: output_file_name, restart_file_name

        ! number of dimensions in the file
        integer :: n_dims = 0

        ! list of netcdf dimension IDs
        integer :: dim_ids(kMAX_DIMENSIONS)
        ! name of the dimensions in the file
        character(len=kMAX_DIM_LENGTH) :: dimensions(kMAX_DIMENSIONS)

  contains

        procedure, public  :: write_file
        procedure, public  :: read_file
        procedure, public  :: read_restart_file
        procedure, public  :: close_files
        procedure, public  :: init => init_ioserver
        procedure, private  :: setup_nest_types
        procedure, private  :: setup_nest_types_2d
        procedure, private  :: setup_nest_types_3d_init
        procedure, private  :: build_nest_geometry
        procedure, public  :: gather_forcing
        procedure, public  :: gather_forcing_2d
        procedure, public  :: gather_forcing_3d_init
        procedure, public  :: distribute_forcing
        procedure, public  :: distribute_init_forcing
        procedure, private  :: scatter_forcing
  end type

  interface

        !>----------------------------------------------------------
        !! Initialize the object (e.g. allocate the variables array)
        !!
        !!----------------------------------------------------------
        module subroutine init_ioserver(this, options, nest_indx)
            implicit none
            class(ioserver_t),   intent(inout) :: this
            type(options_t), intent(in)    :: options(:)
            integer,             intent(in)    :: nest_indx

        end subroutine init_ioserver

        !>----------------------------------------------------------
        !! Increase the size of the variables array if necessary
        !!
        !!----------------------------------------------------------
        module subroutine write_file(this)
            implicit none
            class(ioserver_t),   intent(inout)  :: this
        end subroutine

        !>----------------------------------------------------------
        !! Set the nest datatypes used to gather forcing data for child ioserver
        !!
        !!----------------------------------------------------------
        module subroutine setup_nest_types(this, child_ioserver, child_indx, send_nest_types, buffer_nest_types)
            class(ioserver_t), intent(inout) :: this
            type(ioserver_t), intent(in)    :: child_ioserver
            integer,            intent(in)  :: child_indx
            type(MPI_Datatype), intent(out) :: send_nest_types(:), buffer_nest_types(:)
        end subroutine

        module subroutine setup_nest_types_2d(this, child_ioserver, child_indx, send_nest_types, buffer_nest_types)
            class(ioserver_t), intent(inout) :: this
            type(ioserver_t), intent(in)    :: child_ioserver
            integer,            intent(in)  :: child_indx
            type(MPI_Datatype), intent(out) :: send_nest_types(:), buffer_nest_types(:)
        end subroutine

        module subroutine setup_nest_types_3d_init(this, child_ioserver, child_indx, send_nest_types, buffer_nest_types)
            class(ioserver_t), intent(inout) :: this
            type(ioserver_t), intent(in)    :: child_ioserver
            integer,            intent(in)  :: child_indx
            type(MPI_Datatype), intent(out) :: send_nest_types(:), buffer_nest_types(:)
        end subroutine

        ! One-time init transfer for a parent nest: resolve the per-family
        ! nz_init_3d (max across this parent and its children), grow this
        ! parent's gather_buffer_3d_init if needed, run gather_forcing_2d +
        ! gather_forcing_3d_init, then push the init 2D + 3D restart vars
        ! to every child via MPI_Alltoallw + MPI_Send. Called once per
        ! parent, guarded by initial_nest_done.
        module subroutine distribute_init_forcing(this, components, child_nests)
            class(ioserver_t), intent(inout) :: this
            type(comp_arr_t),  intent(inout) :: components(:)
            integer,           intent(in)    :: child_nests(:)
        end subroutine

        module subroutine build_nest_geometry(this, child_ioserver, geom)
            class(ioserver_t),     intent(inout) :: this
            type(ioserver_t),      intent(in)    :: child_ioserver
            type(nest_geometry_t), intent(inout) :: geom
        end subroutine

        !>----------------------------------------------------------
        !! Set the domain data structure to be used when writing
        !!
        !!----------------------------------------------------------
        module subroutine read_file(this)
            implicit none
            class(ioserver_t), intent(inout)   :: this
        end subroutine
        
        ! Same as above, but for restart file
        module subroutine read_restart_file(this, options)
            class(ioserver_t),   intent(inout) :: this
            type(options_t),     intent(in)    :: options
        end subroutine

        !>----------------------------------------------------------
        !! Gather the forcing data from the parent domain ioservers, and pass to child processes ioservers 
        !!
        !!----------------------------------------------------------
        module subroutine gather_forcing(this)
            class(ioserver_t), intent(inout) :: this
        end subroutine

        module subroutine gather_forcing_2d(this)
            class(ioserver_t), intent(inout) :: this
        end subroutine

        module subroutine gather_forcing_3d_init(this)
            class(ioserver_t), intent(inout) :: this
        end subroutine

        module subroutine distribute_forcing(this, child_ioserver, child_indx)
            class(ioserver_t), intent(inout) :: this
            class(ioserver_t), intent(inout) :: child_ioserver
            integer, intent(in) :: child_indx
        end subroutine

        !>----------------------------------------------------------
        !! Scatter the forcing data to the child processes ioclients
        !!
        !!----------------------------------------------------------
        module subroutine scatter_forcing(this)
            class(ioserver_t), intent(inout) :: this
        end subroutine

        !>----------------------------------------------------------
        !! Close the files
        !!
        !!----------------------------------------------------------
        module subroutine close_files(this)
            class(ioserver_t), intent(inout) :: this
        end subroutine
      
      
  end interface
end module
