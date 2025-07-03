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
  use variable_interface, only : variable_t
  use reader_interface,   only : reader_t
  use output_interface,   only : output_t
  use options_interface,  only : options_t
  use time_object,        only : Time_type
  use time_delta_object,  only : time_delta_t
  use boundary_interface, only : boundary_t
  use flow_object_interface, only : flow_obj_t
  
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
        
        real, dimension(:,:,:,:), pointer :: gather_buffer

        type(buffer_3d_t), allocatable :: write_buffer_3d(:), read_buffer(:), child_gather_buffers(:)
        type(buffer_2d_t), allocatable :: write_buffer_2d(:)

        type(MPI_Win) :: write_win_3d, write_win_2d, read_win, nest_win
        type(MPI_Group) :: children_group

        ! These MPI datatypes describe the access patterns between the IO read/write buffers and
        ! the child read/write buffers
        type(MPI_Datatype), allocatable, dimension(:) :: get_types_3d, get_types_2d, rst_types_3d, rst_types_2d, put_types, child_get_types_3d, child_get_types_2d, &
                                                        child_put_types, force_types, child_force_types, child_rst_types_3d, child_rst_types_2d
        type(MPI_Datatype), allocatable, dimension(:,:) :: send_nest_types, buffer_nest_types
        
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
        integer :: nx_w, nz_w, ny_w, n_w_3d, n_w_2d, nx_r, nz_r, ny_r, n_r, n_f, nx_re, ny_re
        integer         :: ide, kde, jde
        integer :: restart_counter = 0
        integer :: output_counter = 0
        integer :: frames_per_outfile, restart_count
        
        !the indices of the output buffer corresponding to the restart vars
        integer, public, allocatable :: out_var_indices(:), rst_var_indices(:)
        !the names of the restart vars, indexed the same as the above array
        character(len=kMAX_NAME_LENGTH), public, allocatable :: rst_var_names(:)

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
        procedure, public  :: gather_forcing
        procedure, public  :: distribute_forcing
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
        module subroutine setup_nest_types(this, child_ioserver, send_nest_types, buffer_nest_types)
            class(ioserver_t), intent(inout) :: this
            type(ioserver_t), intent(in)    :: child_ioserver
            type(MPI_Datatype), intent(out) :: send_nest_types(:), buffer_nest_types(:)
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

        module subroutine distribute_forcing(this, child_ioserver, child_indx)
            class(ioserver_t), intent(inout) :: this
            class(ioserver_t), intent(inout) :: child_ioserver
            integer, intent(in) :: child_indx
        end subroutine

        !>----------------------------------------------------------
        !! Scatter the forcing data to the child processes ioclients
        !!
        !!----------------------------------------------------------
        module subroutine scatter_forcing(this, forcing_buffer)
            class(ioserver_t), intent(inout) :: this
            real, dimension(:,:,:,:), intent(in) :: forcing_buffer
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
