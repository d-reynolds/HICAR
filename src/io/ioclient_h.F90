!>----------------------------------------------------------
!!  Define the interface for the ioclient
!!
!!  The I/O client serves as a parallelized interface between the main program
!!  and the I/O routines. The I/O client handles the buffering and exchange
!!  of data between the processes of the main program (child processes)
!!  and the I/O processes (parent processes)
!! 
!!  @author
!!  Dylan Reynolds (dylan.reynolds@slf.ch)
!!
!!----------------------------------------------------------
module ioclient_interface
  use mpi_f08
  use icar_constants
  use variable_interface, only : variable_t
  use boundary_interface, only : boundary_t
  use domain_interface,   only : domain_t
  use options_interface,  only : options_t

!  use time_object,        only : Time_type, THREESIXTY, GREGORIAN, NOCALENDAR, NOLEAP

  implicit none

  private
  public :: ioclient_t

  !>----------------------------------------------------------
  !! Output type definition
  !!
  !!----------------------------------------------------------
  type :: ioclient_t

      ! all components are private and should only be modified through procedures
      private

      ! Store the variables to be written
      ! Note n_variables may be smaller then size(variables) so that it doesn't
      ! have to keep reallocating variables whenever something is added or removed
      integer, public :: n_input_variables, n_output_variables
      
      type(MPI_Comm), public :: parent_comms

      type(variable_t), public, allocatable :: variables(:)
      ! time variable , publicis stored outside of the variable list... probably need to think about that some
                  
      integer, public :: server
      
      integer, public ::  i_s_w, i_e_w, k_s_w, k_e_w, j_s_w, j_e_w, n_w, i_s_r, i_e_r, k_s_r, k_e_r, j_s_r, j_e_r, n_r, i_s_re, i_e_re, j_s_re, j_e_re, ide, kde, jde
      integer :: restart_counter = 0
      integer :: output_counter = 0
      integer :: frames_per_outfile, restart_count

      logical :: written = .False.
      logical :: nest_updated = .False.

      real, dimension(:,:,:,:), pointer :: read_buffer, write_buffer_3d, forcing_buffer
      real, dimension(:,:,:),   pointer :: write_buffer_2d
      type(MPI_Win) :: write_win_3d, write_win_2d, read_win, forcing_win
      type(MPI_Group) :: parent_group
      character(len=kMAX_NAME_LENGTH) :: vars_for_nest(kMAX_STORAGE_VARS)

  contains

      procedure, public  :: push
      procedure, public  :: receive
      procedure, public  :: receive_rst
      procedure, public  :: update_nest
      procedure, public  :: init => init_ioclient
  end type

  interface

      !>----------------------------------------------------------
      !! Initialize the object (e.g. allocate the variables array)
      !!
      !!----------------------------------------------------------
    module subroutine init_ioclient(this, domain, forcing, options, n_indx)
        implicit none
        class(ioclient_t),  intent(inout)  :: this
        type(domain_t),     intent(inout)  :: domain
        type(boundary_t),   intent(in)     :: forcing
        type(options_t),    intent(in)     :: options(:)
        integer,            intent(in)     :: n_indx
    end subroutine init_ioclient

      !>----------------------------------------------------------
      !! Push output data to IO buffer
      !!
      !!----------------------------------------------------------
      module subroutine push(this, domain)
          implicit none
          class(ioclient_t),   intent(inout) :: this
          type(domain_t),   intent(inout)    :: domain
      end subroutine

      !>----------------------------------------------------------
      !! Receive input data
      !!
      !!----------------------------------------------------------
      module subroutine receive(this, forcing, domain)
          implicit none
          class(ioclient_t), intent(inout) :: this
          type(boundary_t), intent(inout)  :: forcing
          type(domain_t),   intent(inout)  :: domain
      end subroutine

      !>----------------------------------------------------------
      !! Receive restart data
      !!
      !!----------------------------------------------------------
      module subroutine receive_rst(this, domain, options)
          implicit none
          class(ioclient_t), intent(inout) :: this
          type(domain_t),   intent(inout)  :: domain
          type(options_t),  intent(in)     :: options

      end subroutine

        !>----------------------------------------------------------
        !! Update the nest
        !!
        !!----------------------------------------------------------
        module subroutine update_nest(this, domain)
            implicit none
            class(ioclient_t), intent(inout) :: this
            type(domain_t),    intent(in) :: domain
        end subroutine

  end interface
end module
