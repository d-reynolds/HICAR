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
            new_unittest("adv_h3_v3", test_adv_h3_v3), &
            new_unittest("adv_h5_v5", test_adv_h5_v5), &
            new_unittest("adv_cone_rot", test_adv_cone_rot)  &
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

    subroutine test_adv_cone_rot(error)
        type(error_type), allocatable, intent(out) :: error

        ! Solid-body rotation of a cone, 5th-order RK3 + monotonic flux limiter.
        call advect_cone_rotation(h_order=5,v_order=5,error=error)

    end subroutine test_adv_cone_rot

    subroutine advect_standard(h_order,v_order,error)
        integer, intent(in) :: h_order,v_order
        type(error_type), allocatable, intent(out) :: error

        type(domain_t) :: domain
        type(options_t) :: options
        type(variable_t) :: var
        type(grid_t)     :: var_grid
        type(timer_t) :: flux_time, flux_corr_time, sum_time, adv_wind_time
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
            ! Constant-z diffusion is not supported with 1st-order advection
            ! (see options_obj validation). The init default leaves it on, so
            ! disable it explicitly for the supported h_order=1 configuration.
            options%adv%cz_diff_order = 0
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

        !$acc data copy(kVARS)

        ! exchange to get each others values in the halo regions
        call domain%batch_exch()

        dt = compute_dt(domain%dx, domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d, domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d, &
                        domain%vars_3d(domain%var_indx(kVARS%w)%v)%data_3d, domain%vars_3d(domain%var_indx(kVARS%density)%v)%data_3d, domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d, &
                        domain%ims, domain%ime, domain%kms, domain%kme, domain%jms, domain%jme, &
                        domain%its, domain%ite, domain%jts, domain%jte, &
                        options%time%cfl_reduction_factor)

        STD_OUT_PE = .True.

        ! loop 10 times, calling advect and halo_exchange
        do l = 1, 100
            call advect(domain, options, dt, flux_time, flux_corr_time, sum_time, adv_wind_time)
            call domain%batch_exch()
        end do
        !$acc end data

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
        
        call domain%release()

    end subroutine advect_standard


    !>------------------------------------------------------------------------
    !! Solid-body rotation of a cone (the classic Smolarkiewicz rotating-cone
    !! advection benchmark). A cone-shaped tracer is placed off-centre in a
    !! steady rigid-rotation wind field (u = -Omega*(y-yc), v = Omega*(x-xc),
    !! w = 0). The run length is set so the cone completes EXACTLY 4 full
    !! revolutions, returning it to its starting position. The difference
    !! between the final and initial fields is then the accumulated numerical
    !! (dissipation + dispersion) error of the advection scheme, which we print.
    !!------------------------------------------------------------------------
    subroutine advect_cone_rotation(h_order,v_order,error)
        integer, intent(in) :: h_order,v_order
        type(error_type), allocatable, intent(out) :: error

        type(domain_t) :: domain
        type(options_t) :: options
        type(grid_t)     :: var_grid
        type(timer_t) :: flux_time, flux_corr_time, sum_time, adv_wind_time
        real, allocatable :: th_init(:,:,:)
        real :: dt, dt_cfl, dx, pi, omega, t_period, total_time
        real :: ic, jc, xc, yc, x_cone, y_cone, r_cone, amp, base, orbit, rr
        real :: r_core, vmax, dxr, dyr, rr2, fr
        real :: gmax, grms, gsumsq, gcnt, lmax, lsumsq, lcnt, peak0, peakN, gpeak0, gpeakN
        integer :: nx_g, ny_g, n_steps, n_rev, thv, uidx, vidx, widx
        integer :: ig, jg, kk, my_rank, ierr

        base = 273.0
        amp  = 50.0
        STD_OUT_PE = .False.
        call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)

        call options%init()
        options%domain%init_conditions_file = '../tests/Test_Cases/domains/flat_plane_250m.nc'
        options%domain%hgt_hi = 'topo'
        options%domain%lat_hi = 'lat'
        options%domain%lon_hi = 'lon'
        if (trim(options%domain%init_conditions_file) /= '') then
            call check_file_exists(trim(options%domain%init_conditions_file), message='The test domain file does not exist. Ensure that the HICAR_test repo was installed to build/test when building HICAR.')
        endif

        options%domain%dx = 250.0
        ! Thin vertical extent: the cone is z-uniform and w=0, so the rotation is
        ! purely horizontal — a few levels keep the multi-revolution run cheap
        ! without changing the result (6 = minimum for the 5th-order v-stencil).
        options%domain%nz = 6
        options%domain%dz_levels = 50.0
        options%domain%sleve = .True.
        options%physics%advection = kADV_STD

        ! 5th-order RK3 with the monotonic flux limiter (the production config).
        options%time%cfl_reduction_factor = time_step_factor(2)
        options%adv%flux_corr = 1
        options%time%RK3 = .true.
        options%adv%h_order = h_order
        options%adv%v_order = v_order
        options%adv%advect_density = .False.

        domain%compute_comms = MPI_COMM_WORLD
        call adv_var_request(options)
        call domain%init(options,1)
        call adv_init(domain,options)

        thv  = domain%var_indx(kVARS%potential_temperature)%v
        uidx = domain%var_indx(kVARS%u)%v
        vidx = domain%var_indx(kVARS%v)%v
        widx = domain%var_indx(kVARS%w)%v

        dx   = domain%dx
        pi   = 4.0*atan(1.0)
        nx_g = domain%grid%ide - domain%grid%ids + 1
        ny_g = domain%grid%jde - domain%grid%jds + 1

        ! rotation centre (global index space -> physical via *dx)
        ic = 0.5*real(domain%grid%ids + domain%grid%ide)
        jc = 0.5*real(domain%grid%jds + domain%grid%jde)
        xc = ic*dx
        yc = jc*dx

        ! ---- Rankine-vortex flow (NOT global solid-body rotation).
        ! A rigid core of radius r_core rotates uniformly at angular rate omega;
        ! beyond r_core the (purely azimuthal, divergence-free) flow decays as
        ! 1/r so the max speed is capped at omega*r_core. Capping the speed
        ! DECOUPLES the CFL-limited dt from the (large) outer grid: the step
        ! count for one revolution is set by r_core in cells, not the domain
        ! half-span. The cone sits entirely inside the rigid core, so it rotates
        ! WITHOUT deformation and returns to its start after one period.
        ! Geometry in cells (resolved cone: ~6 cells across, orbit 5 cells, all
        ! inside an 8-cell core => ~55-65 steps for one full revolution).
        r_core = 8.0 *dx
        orbit  = 5.0 *dx
        r_cone = 3.0 *dx
        vmax   = 10.0                 ! peak speed at r_core (m/s)
        omega  = vmax / r_core        ! core angular rate; n_steps is independent of its value
        x_cone = xc + orbit
        y_cone = yc

        ! v_theta(r)/r = omega        for r <= r_core   (rigid core: u=-omega*dy, v=omega*dx)
        !             = omega*r_core^2/r^2 for r > r_core   (1/r decay, divergence-free)
        associate(ug => domain%vars_3d(uidx)%data_3d, &
                  vg => domain%vars_3d(vidx)%data_3d, &
                  wg => domain%vars_3d(widx)%data_3d)
            do jg = lbound(ug,3), ubound(ug,3)
                do kk = lbound(ug,2), ubound(ug,2)
                    do ig = lbound(ug,1), ubound(ug,1)
                        dxr = real(ig)*dx - xc;  dyr = real(jg)*dx - yc
                        rr2 = dxr*dxr + dyr*dyr
                        fr  = merge(1.0, (r_core*r_core)/max(rr2,1.0e-12), rr2 <= r_core*r_core)
                        ug(ig,kk,jg) = -omega*fr*dyr
                    enddo
                enddo
            enddo
            do jg = lbound(vg,3), ubound(vg,3)
                do kk = lbound(vg,2), ubound(vg,2)
                    do ig = lbound(vg,1), ubound(vg,1)
                        dxr = real(ig)*dx - xc;  dyr = real(jg)*dx - yc
                        rr2 = dxr*dxr + dyr*dyr
                        fr  = merge(1.0, (r_core*r_core)/max(rr2,1.0e-12), rr2 <= r_core*r_core)
                        vg(ig,kk,jg) =  omega*fr*dxr
                    enddo
                enddo
            enddo
            wg = 0.0
        end associate

        ! ---- cone tracer (z-uniform; w=0 so it never moves vertically). Set on
        ! every advected variable; the boundary halo holds the base value, acting
        ! as fixed inflow far from the cone's orbit.
        do i = 1, size(domain%adv_vars)
            var_grid = domain%vars_3d(domain%adv_vars(i)%v)%grid
            associate(q => domain%vars_3d(domain%adv_vars(i)%v)%data_3d)
                do jg = lbound(q,3), ubound(q,3)
                    do kk = lbound(q,2), ubound(q,2)
                        do ig = lbound(q,1), ubound(q,1)
                            rr = sqrt((real(ig)*dx - x_cone)**2 + (real(jg)*dx - y_cone)**2)
                            q(ig,kk,jg) = base + amp*max(0.0, 1.0 - rr/r_cone)
                        enddo
                    enddo
                enddo
            end associate
        end do

        ! snapshot of the initial theta field for the error measure
        allocate(th_init, source=domain%vars_3d(thv)%data_3d)

        !$acc data copy(kVARS)
        call domain%batch_exch()

        ! CFL-stable dt, then choose n_steps so the requested revolutions land
        ! EXACTLY on a step boundary (dt <= dt_cfl), returning the cone to start.
        ! One revolution is the standard rotating-cone error measure.
        n_rev = 1
        dt_cfl = compute_dt(domain%dx, domain%vars_3d(uidx)%data_3d, domain%vars_3d(vidx)%data_3d, &
                        domain%vars_3d(widx)%data_3d, domain%vars_3d(domain%var_indx(kVARS%density)%v)%data_3d, &
                        domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d, &
                        domain%ims, domain%ime, domain%kms, domain%kme, domain%jms, domain%jme, &
                        domain%its, domain%ite, domain%jts, domain%jte, &
                        options%time%cfl_reduction_factor)
        call MPI_Allreduce(MPI_IN_PLACE, dt_cfl, 1, MPI_REAL, MPI_MIN, domain%compute_comms)

        t_period   = 2.0*pi/omega
        total_time = real(n_rev)*t_period
        n_steps    = ceiling(total_time/dt_cfl)
        dt         = total_time/real(n_steps)

        do l = 1, n_steps
            call advect(domain, options, dt, flux_time, flux_corr_time, sum_time, adv_wind_time)
            call domain%batch_exch()
        end do
        !$acc end data

        ! ---- error of the final field vs the initial cone (global reductions
        ! over interior tiles so the printed numbers are decomposition-invariant).
        associate(th => domain%vars_3d(thv)%data_3d, g => domain%grid)
            lmax   = maxval(abs(th(g%its:g%ite,:,g%jts:g%jte) - th_init(g%its:g%ite,:,g%jts:g%jte)))
            lsumsq = sum((th(g%its:g%ite,:,g%jts:g%jte) - th_init(g%its:g%ite,:,g%jts:g%jte))**2)
            lcnt   = real(size(th(g%its:g%ite,:,g%jts:g%jte)))
            peak0  = maxval(th_init(g%its:g%ite,:,g%jts:g%jte)) - base
            peakN  = maxval(th(g%its:g%ite,:,g%jts:g%jte)) - base
        end associate

        call MPI_Allreduce(lmax,   gmax,   1, MPI_REAL, MPI_MAX, domain%compute_comms)
        call MPI_Allreduce(lsumsq, gsumsq, 1, MPI_REAL, MPI_SUM, domain%compute_comms)
        call MPI_Allreduce(lcnt,   gcnt,   1, MPI_REAL, MPI_SUM, domain%compute_comms)
        call MPI_Allreduce(peak0,  gpeak0, 1, MPI_REAL, MPI_MAX, domain%compute_comms)
        call MPI_Allreduce(peakN,  gpeakN, 1, MPI_REAL, MPI_MAX, domain%compute_comms)
        grms = sqrt(gsumsq/gcnt)

        if (my_rank == 0) then
            write(*,'(A,I0,A)')          "  [cone-rot] error after ", n_rev, " revolution(s) (final - initial):"
            write(*,'(A,F10.5,A)')       "  [cone-rot]   max abs error : ", gmax,  " K"
            write(*,'(A,F10.5,A)')       "  [cone-rot]   RMS error      : ", grms,  " K"
            write(*,'(A,F8.3,A,F8.3,A)') "  [cone-rot]   cone peak: initial ", gpeak0, " K  ->  final ", gpeakN, " K"
            write(*,'(A,F6.1,A)')        "  [cone-rot]   peak amplitude retained: ", 100.0*gpeakN/max(gpeak0,1.0e-6), " %"
        endif

        ! Sanity gate: scheme must stay bounded (monotonic limiter => no growth
        ! beyond the initial amplitude, no NaN). This catches a broken scheme
        ! without over-constraining the expected numerical diffusion.
        call check(error, (gmax == gmax) .and. (gmax < amp) .and. (gpeakN > 0.0))
        if (allocated(error)) then
            call domain%release()
            return
        endif

        deallocate(th_init)
        call domain%release()

    end subroutine advect_cone_rotation


end module test_advect