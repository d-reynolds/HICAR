!>----------------------------------------------------------
!!  Define the interface for the output object
!!
!!  Output objects store all of the data and references to data necessary to write
!!  an output file.  This includes primarily internal netcdf related IDs.
!!  Output objects also store an array of variables to output.
!!  These variables maintain pointers to the data to be output as well as
!!  Metadata (e.g. dimension names, units, other attributes)
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!----------------------------------------------------------
module output_interface
  use mpi_f08
  use netcdf
  use icar_constants
  use variable_interface,       only : variable_t
  use meta_data_interface,      only : meta_data_t
  use time_object,              only : Time_type
  use options_interface,        only : options_t

  implicit none

  private
  public :: output_t

  !>----------------------------------------------------------
  !! Output type definition
  !!
  !!----------------------------------------------------------
  type, extends(meta_data_t) :: output_t

      ! all components are private and should only be modified through procedures
      private

      ! Store the variables to be written
      ! Note n_variables may be smaller then size(variables) so that it doesn't
      ! have to keep reallocating variables whenever something is added or removed
      integer, public :: n_vars = 0
      type(meta_data_t), public, allocatable :: var_meta(:)
      type(variable_t), public, allocatable :: variables(:)
      ! time variable , publicis stored outside of the variable list... probably need to think about that some
      type(meta_data_t) :: time

      ! store status of the object
      logical :: is_initialized = .false.
      logical :: creating = .false.
      logical :: block_checked = .false.
      logical :: is_blocked = .false.
      logical :: blocked_LL = .false.
      logical :: blocked_UR = .false.

      ! The filename of the netcdf file to write
      character(len=kMAX_FILE_LENGTH), public :: output_fn, base_out_file_name, base_rst_file_name, restart_fn

      character(len=49)   :: file_date_format = '(I4,"-",I0.2,"-",I0.2,"_",I0.2,"-",I0.2,"-",I0.2)'

      ! the netcdf ID for an open file
      integer :: out_ncfile_id = -1
      integer :: rst_ncfile_id = -1
      
      integer :: active_nc_id = -1 !Used during save events to control which file is written to
      
      ! number of dimensions in the file
      integer :: n_dims = 0
      
      integer :: start_3d(3), cnt_3d(3), cnt_2d(2)
      
      !same as above, but for block, if present
      integer :: start_3d_b(3), cnt_3d_b(3), cnt_2d_b(2)
      integer :: start_3d_b2(3), cnt_3d_b2(3), cnt_2d_b2(2)

      integer :: its, ite, kts, kte, jts, jte

      integer :: output_counter, output_count

      integer :: file_dim_len(3)

      character(len=kMAX_NAME_LENGTH) :: time_units

  contains

      procedure, public  :: set_attrs
      procedure, public  :: save_out_file
      procedure, public  :: save_rst_file
      procedure, public  :: close_output_files

      procedure, public  :: init => init_output
      procedure, public  :: init_restart
      procedure, private :: increase_var_capacity
      procedure, private :: add_to_output
      procedure, private :: add_variables
  end type

  interface

      !>----------------------------------------------------------
      !! Initialize the object (e.g. allocate the variables array)
      !!
      !!----------------------------------------------------------
      module subroutine init_output(this, options, its, ite, kts, kte, jts, jte, ide, kde, jde)
        implicit none
        class(output_t),  intent(inout)  :: this
        type(options_t),  intent(in)     :: options
        integer,          intent(in)     :: its, ite, kts, kte, jts, jte, ide, kde, jde

      end subroutine init_output

      !>----------------------------------------------------------
      !! Initialize the output counter and file after a restart
      !!
      !!----------------------------------------------------------
      module subroutine init_restart(this, options, par_comms, out_var_indices)
        implicit none
        class(output_t),  intent(inout)  :: this
        type(options_t),  intent(in)     :: options
        type(MPI_Comm),   intent(in)     :: par_comms
        integer,          intent(in)     :: out_var_indices(:)
      end subroutine

      !>----------------------------------------------------------
      !! Increase the size of the variables array if necessary
      !!
      !!----------------------------------------------------------
      module subroutine increase_var_capacity(this)
          implicit none
          class(output_t),   intent(inout)  :: this
      end subroutine

      !>----------------------------------------------------------
      !! Set the domain data structure to be used when writing
      !!
      !!----------------------------------------------------------
      module subroutine set_attrs(this, options)
          implicit none
          class(output_t),  intent(inout)  :: this
          type(options_t),  intent(in)     :: options
      end subroutine

      !>----------------------------------------------------------
      !! Add a variable to the list of output variables
      !!
      !!----------------------------------------------------------
      module subroutine add_to_output(this, in_variable)
          implicit none
          class(output_t),   intent(inout)  :: this
          type(meta_data_t), intent(in)     :: in_variable
      end subroutine
      
      !>----------------------------------------------------------
      !! Essentially copy domain vars_to_out to a local array
      !!
      !!----------------------------------------------------------
      module subroutine add_variables(this, vars_to_out)
          class(output_t),  intent(inout)  :: this
          integer, dimension(:), intent(in):: vars_to_out
        end subroutine
      !>----------------------------------------------------------
      !! Save a new timestep (time) to the output file 
      !!
      !!----------------------------------------------------------
      module subroutine save_out_file(this, time, par_comms, out_var_indices)
          implicit none
          class(output_t),  intent(inout) :: this
          type(Time_type),  intent(in)    :: time
          type(MPI_Comm),   intent(in)    :: par_comms
          integer,          intent(in)    :: out_var_indices(:)
      end subroutine
      
      !>----------------------------------------------------------
      !! Save a new timestep (time) to the restart file 
      !!
      !!----------------------------------------------------------
      module subroutine save_rst_file(this, time, par_comms, rst_var_indices)
        implicit none
        class(output_t),  intent(inout) :: this
        type(Time_type),  intent(in)    :: time
        type(MPI_Comm),   intent(in)    :: par_comms
        integer,          intent(in)    :: rst_var_indices(:)
    end subroutine

      module subroutine close_output_files(this)
        implicit none
        class(output_t),   intent(inout)  :: this
      end subroutine

  end interface
end module
