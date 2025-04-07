! test_advect.F90
! This program initializes a grid and sets up MPI communication for a simple advection test.
! It tests all possible configurations of advection, of varying numerical orders, and use of the monotonic flux limiter

module test_advect
          
    use variable_interface,      only: variable_t
    use mpi_f08
    use icar_constants
    use testdrive, only : new_unittest, unittest_type, error_type, check, test_failed
    use domain_interface, only: domain_t
    use options_interface, only: options_t
    use advection, only: adv_init, advect, adv_var_request
    use timer_interface, only: timer_t
    use io_routines,     only : check_file_exists
    use time_step, only: compute_dt
    use grid_interface, only: grid_t
    implicit none
    private
    
    ! Hard coded time step factor to be used with RK3 time stepping for 3rd and 5th order advection
    real :: time_step_factor(2) = (/ 1.4, 1.6/)
    integer :: h_order(3) = (/ 1, 3, 5/)
    integer :: v_order(3) = (/ 1, 3, 5/)
    integer :: i, j, k, l

    public :: collect_advect_suite
    
    contains
    


    !> Collect all exported unit tests
    subroutine collect_advect_suite(testsuite)
        !> Collection of tests
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
        
        testsuite = [ &
            new_unittest("adv_h1_v1", test_adv_h1_v1), &
            ! new_unittest("adv_h3_v1", test_adv_h3_v1), &
            ! new_unittest("adv_h5_v1", test_adv_h5_v1), &
            new_unittest("adv_h3_v3", test_adv_h3_v3), &
            ! new_unittest("adv_h5_v3", test_adv_h5_v3), &
            new_unittest("adv_h5_v5", test_adv_h5_v5)  &
            ]
    
    end subroutine collect_advect_suite
    
    subroutine test_adv_h1_v1(error)
        type(error_type), allocatable, intent(out) :: error

        call advect_standard(h_order=1,v_order=1,error=error)

    end subroutine test_adv_h1_v1
    
    subroutine test_adv_h3_v1(error)
        type(error_type), allocatable, intent(out) :: error

        call advect_standard(h_order=3,v_order=1,error=error)

    end subroutine test_adv_h3_v1

    subroutine test_adv_h5_v1(error)
        type(error_type), allocatable, intent(out) :: error

        call advect_standard(h_order=5,v_order=1,error=error)

    end subroutine test_adv_h5_v1

    subroutine test_adv_h3_v3(error)
        type(error_type), allocatable, intent(out) :: error

        call advect_standard(h_order=3,v_order=3,error=error)

    end subroutine test_adv_h3_v3

    subroutine test_adv_h5_v3(error)
        type(error_type), allocatable, intent(out) :: error

        call advect_standard(h_order=5,v_order=3,error=error)

    end subroutine test_adv_h5_v3

    subroutine test_adv_h3_v5(error)
        type(error_type), allocatable, intent(out) :: error

        call advect_standard(h_order=3,v_order=5,error=error)

    end subroutine test_adv_h3_v5

    subroutine test_adv_h5_v5(error)
        type(error_type), allocatable, intent(out) :: error

        call advect_standard(h_order=5,v_order=5,error=error)

    end subroutine test_adv_h5_v5

    subroutine advect_standard(h_order,v_order,error)
        integer, intent(in) :: h_order,v_order
        type(error_type), allocatable, intent(out) :: error

        type(domain_t) :: domain
        type(options_t) :: options
        type(variable_t) :: var
        type(grid_t)     :: var_grid
        type(timer_t) :: flux_time, flux_up_time, flux_corr_time, sum_time, adv_wind_time
        real :: max_u, max_v, max_w, max_wind
        real :: dt, t_step_f, initial_val, spatial_constraint
        integer :: fluxcorr
        integer :: i,j,k,l

        initial_val = 273.0
        STD_OUT_PE = .False.

        call options%init()
        options%domain%init_conditions_file = '../tests/Test_Cases/domains/flat_plane_250m.nc'
        options%domain%hgt_hi = 'topo'
        options%domain%lat_hi = 'lat'
        options%domain%lon_hi = 'lon'
        ! check that the init_conditions_file exists
        if (trim(options%domain%init_conditions_file) /= '') then
            call check_file_exists(trim(options%domain%init_conditions_file), message='The test domain file does not exist. Ensure that the HICAR_test repo was installed to build/test when building HICAR.')
        endif

        options%domain%dx = 250.0
        options%domain%nz = 20
        options%domain%dz_levels = 50.0
        options%domain%sleve = .True.

        options%domain%lat_hi = 'lat'
        options%domain%lon_hi = 'lon'
        options%domain%hgt_hi = 'topo'
        options%physics%advection = kADV_STD

        ! setup t_step and fluxcorr based on the horizontal and vertical orders
        if (h_order == 1) then
            options%time%cfl_reduction_factor = 1.0
            options%adv%flux_corr = 0
            options%time%RK3 = .False.
        else if (h_order == 3) then
            options%time%cfl_reduction_factor = time_step_factor(1)
            options%adv%flux_corr = 1
            options%time%RK3 = .true.
        else if (h_order == 5) then
            options%time%cfl_reduction_factor = time_step_factor(2)
            options%adv%flux_corr = 1
            options%time%RK3 = .true.
        end if

        options%adv%h_order = h_order
        options%adv%v_order = v_order
        options%adv%advect_density = .False.

        ! Initialize domain and advection modules
        domain%compute_comms = MPI_COMM_WORLD

        call adv_var_request(options)
        call domain%init(options,1)

        call adv_init(domain,options)

        ! Set values for u,v,w and the tracer to advect
        domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d =  10.0
        domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d =  -10.0
        domain%vars_3d(domain%var_indx(kVARS%w)%v)%data_3d =  0.0

        do i = 1, size(domain%adv_vars)
            var_grid = domain%vars_3d(domain%adv_vars(i)%v)%grid

            domain%vars_3d(domain%adv_vars(i)%v)%data_3d = -99999.0
            domain%vars_3d(domain%adv_vars(i)%v)%data_3d(var_grid%its:var_grid%ite,:,var_grid%jts:var_grid%jte) = initial_val

            !copy over initial value like "forcing" at the domain edges
            if (domain%north_boundary) domain%vars_3d(domain%adv_vars(i)%v)%data_3d(:,:,var_grid%jte+1:var_grid%jme) = initial_val
            if (domain%south_boundary) domain%vars_3d(domain%adv_vars(i)%v)%data_3d(:,:,var_grid%jms:var_grid%jts-1) = initial_val
            if (domain%east_boundary) domain%vars_3d(domain%adv_vars(i)%v)%data_3d(var_grid%ite+1:var_grid%ime,:,:) = initial_val
            if (domain%west_boundary) domain%vars_3d(domain%adv_vars(i)%v)%data_3d(var_grid%ims:var_grid%its-1,:,:) = initial_val
        end do

        ! exchange to get each others values in the halo regions
        call domain%batch_exch()

        dt = compute_dt(domain%dx, domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d, domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d, &
                        domain%vars_3d(domain%var_indx(kVARS%w)%v)%data_3d, domain%vars_3d(domain%var_indx(kVARS%density)%v)%data_3d, options%domain%dz_levels, &
                        domain%ims, domain%ime, domain%kms, domain%kme, domain%jms, domain%jme, &
                        domain%its, domain%ite, domain%jts, domain%jte, &
                        options%time%cfl_reduction_factor, &
                        use_density=.false.)

        STD_OUT_PE = .True.

        ! loop 10 times, calling advect and halo_exchange
        do l = 1, 100
            call advect(domain, options, dt, flux_time, flux_up_time, flux_corr_time, sum_time, adv_wind_time)
            call domain%batch_exch()
        end do

        ! test that all of the tile indices of temperature for the domain object are the same as the intial value
        do i = 1, size(domain%adv_vars)
            var = domain%vars_3d(domain%adv_vars(i)%v)
            if (.not.(ALL(var%data_3d(var%grid%its:var%grid%ite,:,var%grid%jts:var%grid%jte) == initial_val))) then
                write(*,*) "Initial value (expected) is: ", initial_val
                write(*,*) "bounds of var are: ", var%grid%its, " ", var%grid%ite, " ", var%grid%jts, " ", var%grid%jte
                write(*,*) "Max val of advected variable in tile: ", maxval(var%data_3d(var%grid%its:var%grid%ite,:,var%grid%jts:var%grid%jte))
                write(*,*) "Min val of advected variable in tile: ", minval(var%data_3d(var%grid%its:var%grid%ite,:,var%grid%jts:var%grid%jte))    
                call test_failed(error,"advect_standard","Tile data is not equal to initial value")
                return
            end if
        end do

    end subroutine advect_standard


end module test_advect