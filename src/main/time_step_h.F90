!> ----------------------------------------------------------------------------
!!  Interface for the main time stepping module.
!!
!!  Interface-only parent module: the implementation lives in the
!!  time_step_implementation submodule (time_step.F90). Keeping the heavy
!!  physics-driver USEs out of this module keeps its .mod file small, which
!!  matters because nvfortran's module import cost grows quadratically with
!!  .mod size and every module above this one re-embeds what it exports.
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!! ----------------------------------------------------------------------------
module time_step
    use options_interface,  only : options_t
    use domain_interface,   only : domain_t
    use time_object,        only : Time_type

    implicit none

    private
    public :: step, compute_dt

    interface

        !>------------------------------------------------------------
        !!  Calculate the maximum stable time step given some CFL criteria
        !!
        !! @param dx  [ scalar ]        horizontal grid cell width  [m]
        !! @param u   [nx+1 x nz x ny]  east west wind speeds       [m/s]
        !! @param v   [nx x nz x ny+1]  North South wind speed      [m/s]
        !! @param w   [nx x nz x ny]    vertical wind speed         [m/s]
        !! @param CFL [ scalar ]        CFL limit to use (e.g. 1.0)
        !! @return dt [ scalar ]        Maximum stable time step    [s]
        !!
        !!------------------------------------------------------------
        module function compute_dt(dx, u, v, w, rho, dz, ims, ime, kms, kme, jms, jme, its, ite, jts, jte, CFL, max_mapfac, err_msg) result(dt)
            real,       intent(in)                   :: dx
            real,       intent(in), dimension(ims:ime+1,kms:kme,jms:jme) :: u
            real,       intent(in), dimension(ims:ime,kms:kme,jms:jme+1) :: v
            real,       intent(in), dimension(ims:ime,kms:kme,jms:jme)   :: w, rho, dz
            integer,    intent(in)                   :: ims, ime, kms, kme, jms, jme, its, ite, jts, jte
            real,       intent(in)                   :: CFL, max_mapfac
            character(len=*), intent(in), optional  :: err_msg
            real :: dt
        end function compute_dt

        !>------------------------------------------------------------
        !!  Step forward one IO time step.
        !!
        !!  Calculated the internal model time step to satisfy the CFL criteria,
        !!  then updates all forcing update increments for that dt and loops through
        !!  time calling physics modules.
        !!  Also checks to see if it is time to write a model output file.
        !!
        !! @param domain    domain data structure containing model state
        !! @param end_time  time to advance the model to
        !! @param options   model options structure
        !!
        !!------------------------------------------------------------
        module subroutine step(domain, end_time, options)
            type(domain_t),     intent(inout)   :: domain
            type(Time_type),    intent(in)      :: end_time
            type(options_t),    intent(in)      :: options
        end subroutine step

    end interface

end module time_step
