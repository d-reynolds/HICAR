!>------------------------------------------------
!! Contains type definitions for a variety of model data strucutres
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!>------------------------------------------------
module data_structures
    use, intrinsic :: iso_c_binding ! needed for fftw compatible complex types
    use icar_constants
    implicit none


! ------------------------------------------------
!   various data structures for use in geographic interpolation routines
! ------------------------------------------------

    type index_type
        integer :: v = -1
        integer :: id = -1
    end type index_type
    ! contains the location of a specific grid point
    type position
        integer::x,y
    end type position
    ! contains location of surrounding 4 grid cells
    type fourpos
        integer::x(4),y(4)
    end type fourpos

    ! a geographic look up table for spatial interpolation, from x,y with weight w
    type geo_look_up_table
        ! x,y index positions, [n by m by 4] where there are 4 surrounding low-res points
        ! for every high resolution point grid point to interpolate to
        integer,allocatable, dimension(:,:,:)   :: x, y
        ! weights to use for each of the 4 surrounding gridpoints.  Sum(over axis 3) must be 1.0
        real,   allocatable, dimension(:,:,:)   :: w
    end type geo_look_up_table

    ! ------------------------------------------------
    ! A look up table for vertical interpolation. from z with weight w
    ! ------------------------------------------------
    type vert_look_up_table
        ! z index positions for all x,y,z points (x 2 for above and below z levels)
        integer,allocatable, dimension(:,:,:,:) :: z

        ! weights to use for each of the two surrounding points.  Sum (over axis 1) must be 1.0
        real,   allocatable, dimension(:,:,:,:) :: w
    end type vert_look_up_table

    ! ------------------------------------------------
    ! generic interpolable type so geo interpolation routines will work on winds, domain, or boundary conditions.
    ! ------------------------------------------------
    type interpolable_type
        ! all interpolables must have position (lat, lon, z)
        real, allocatable, dimension(:,:) :: lat,lon
        real, allocatable, dimension(:,:,:) :: z

        ! these are the look up tables that describe how to interpolate vertically (vert_lut) and horizontally (geolut)
        type(vert_look_up_table)::vert_lut
        type(geo_look_up_table)::geolut

        ! used to keep track of whether or not a particular error has been printed yet for this structure
        logical :: dx_errors_printed=.False.
        logical :: dy_errors_printed=.False.
    end type interpolable_type

    ! ------------------------------------------------
    ! Data type to hold all of the array temporaries required by the lineary theory calculations
    ! e.g. k and l wave number arrays
    ! ------------------------------------------------
    type linear_theory_type
        real,                       allocatable, dimension(:,:) :: sig, k, l, kl
        complex(C_DOUBLE_COMPLEX),  allocatable, dimension(:,:) :: denom, msq, mimag, m, ineta

        complex(C_DOUBLE_COMPLEX),  pointer,     dimension(:,:) :: uhat, vhat
        complex(C_DOUBLE_COMPLEX),  pointer,     dimension(:,:) :: u_perturb, v_perturb
        complex(C_DOUBLE_COMPLEX),  pointer,     dimension(:,:) :: u_accumulator, v_accumulator
        type(C_PTR) :: uh_aligned_data, up_aligned_data, ua_aligned_data
        type(C_PTR) :: vh_aligned_data, vp_aligned_data, va_aligned_data

        type(C_PTR) :: uplan, vplan
    end type linear_theory_type

    ! ------------------------------------------------
    ! Tendency terms output by various physics subroutines
    ! ------------------------------------------------
    type tendencies_type
        ! 3D atmospheric field tendencies
        ! These are used by various physics parameterizations 
        real,   allocatable, dimension(:,:,:) :: th,qv,qc,qi,u,v,qr,qs

        ! advection and pbl tendencies that need to be saved for the cumulus/pbl scheme
        real, allocatable, dimension(:,:,:) :: qv_adv,qv_pbl, th_pbl, qi_pbl, qc_pbl
        real, allocatable, dimension(:,:,:) :: th_lwrad, th_swrad
    end type tendencies_type

end module data_structures
