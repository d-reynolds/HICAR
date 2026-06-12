!> ----------------------------------------------------------------------------
!!  Interface for model initialization.
!!
!!  Interface-only parent module: the implementation lives in the
!!  initialization_implementation submodule (init.F90). Keeping the heavy
!!  physics-driver USEs out of this module keeps its .mod file small, which
!!  matters because nvfortran's module import cost grows quadratically with
!!  .mod size and every module above this one re-embeds what it exports.
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!! ----------------------------------------------------------------------------
module initialization
    use options_interface,     only : options_t
    use domain_interface,      only : domain_t
    use boundary_interface,    only : boundary_t
    use flow_object_interface, only : comp_arr_t
    use ioclient_interface,    only : ioclient_t

    implicit none

    private
    public :: split_processes, init_options, init_model, init_physics, init_model_state

    interface

        module subroutine split_processes(components, ioclient, n_nests, options)
            type(comp_arr_t), intent(inout) :: components(:)
            type(ioclient_t), intent(inout) :: ioclient(:)
            integer, intent(in) :: n_nests
            type(options_t), intent(in) :: options
        end subroutine split_processes

        module subroutine init_options(options, namelist_file, info_only, gen_nml, only_namelist_check)
            type(options_t), allocatable, intent(out) :: options(:)
            character(len=*), intent(in) :: namelist_file
            logical, intent(in), optional :: info_only, gen_nml, only_namelist_check
        end subroutine init_options

        module subroutine init_model(options,domain,boundary,ioclient,nest_indx)
            type(options_t), intent(inout) :: options(:)
            type(domain_t),  intent(inout) :: domain
            type(boundary_t),intent(inout) :: boundary ! forcing file for init conditions
            type(ioclient_t),intent(inout) :: ioclient
            integer, intent(in) :: nest_indx
        end subroutine init_model

        module subroutine init_model_state(options, domain, boundary, ioclient)
            type(options_t), intent(inout) :: options
            type(domain_t),  intent(inout) :: domain
            type(boundary_t),intent(inout) :: boundary ! forcing file for init conditions
            type(ioclient_t),intent(inout) :: ioclient
        end subroutine init_model_state

        module subroutine init_physics(options, domain)
            type(options_t), intent(inout) :: options
            type(domain_t),  intent(inout) :: domain
        end subroutine init_physics

    end interface

end module initialization
