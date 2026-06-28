! test_wind_iterative.F90
! Unit tests for the hand-rolled BiCGStab iterative wind solver (wind_iterative.F90).
!
! The motivating regression: the solver's Krylov halo-exchange buffers
! (east_send/east_recv/...) were allocated with non-1-based bounds (k_s-1:k_e+1,
! ...) and passed whole to the mpi_f08 MPI_Isend/Irecv assumed-rank async buffer
! dummy. gfortran rebinds such an actual's descriptor to 1-based across the call,
! so the pack/unpack then ran off the ends of the buffers. The bug only fires when
! the domain is decomposed (>1 compute rank) so a halo exchange actually runs, and
! it is silent except under -fcheck=bounds (debug). This suite therefore EXISTS to
! be run multi-rank (the CI full-test runs `mpiexec -np 4 HICAR-tester`): with >=2
! ranks at least one of the four face exchanges executes and the bug crashes here.
!
! The test also checks the solve is functional (finite winds, reduced divergence),
! so it doubles as a guard against silent halo-correctness regressions on GPU/release
! where there is no bounds check.

module test_wind_iterative

    use mpi_f08
    use icar_constants
    use testdrive,          only : new_unittest, unittest_type, error_type, test_failed
    use domain_interface,   only : domain_t
    use options_interface,  only : options_t
    use wind,               only : wind_var_request, init_winds, calc_divergence
    use wind_iterative,     only : calc_iter_winds, finalize_iter_winds
    use advection,          only : adv_var_request
    use io_routines,        only : check_file_exists
    implicit none
    private

    real, parameter :: PI = 3.14159265358979

    public :: collect_wind_iterative_suite

contains

    !> Collect all exported unit tests
    subroutine collect_wind_iterative_suite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            new_unittest("iter_wind_solve_decomp", test_iter_wind_solve) &
            ]
    end subroutine collect_wind_iterative_suite


    !> Build the iterative wind solver on the flat test domain, drive one solve of a
    !> seeded divergent wind field, and verify it runs across the MPI decomposition
    !> (exercising exchange_krylov_halos) and reduces the divergence to a finite field.
    subroutine test_iter_wind_solve(error)
        type(error_type), allocatable, intent(out) :: error

        type(domain_t)  :: domain
        type(options_t) :: options
        real, allocatable :: div(:,:,:)
        real    :: div0_max, div1_max, div0_max_g, div1_max_g
        integer :: ierr
        integer :: ims, ime, jms, jme, kms, kme, its, ite, jts, jte
        logical :: ok
        character(len=256) :: msg

        STD_OUT_PE = .False.

        ! --- options: iterative wind solver with a CONSTANT alpha so the solve does
        !     not depend on the Froude / calc_alpha machinery, and with the optional
        !     wind add-ons (Sx, thermal, linear) disabled. ----------------------------
        call options%init()
        options%domain%init_conditions_file = '../tests/Test_Cases/domains/flat_plane_250m.nc'
        options%domain%hgt_hi = 'topo'
        options%domain%lat_hi = 'lat'
        options%domain%lon_hi = 'lon'
        if (trim(options%domain%init_conditions_file) /= '') then
            call check_file_exists(trim(options%domain%init_conditions_file), &
                message='The test domain file does not exist. Ensure the HICAR Test-Data repo was '// &
                        'installed and ../tests/Test_Cases/domains/flat_plane_250m.nc exists.')
        endif
        options%domain%dx              = 250.0
        options%domain%nz              = 20
        options%domain%sleve           = .True.
        options%domain%use_map_factors = .False.
        options%physics%advection      = kADV_STD
        options%physics%windtype       = kITERATIVE_WINDS
        options%wind%alpha_const       = 1.0       ! constant alpha -> skip Froude / calc_alpha
        options%wind%Sx                = .False.
        options%wind%thermal           = .False.
        options%wind%linear_theory     = .False.

        ! --- build the domain over MPI_COMM_WORLD (this is what splits it into tiles,
        !     so the solver's halo exchange actually runs) -----------------------------
        domain%compute_comms = MPI_COMM_WORLD
        call adv_var_request(options)              ! core dynamics vars (u/v/w/density/jaco/dz)
        call wind_var_request(options)             ! wind_alpha / w_real
        call domain%init(options, 1)
        call init_winds(domain, options)           ! -> init_iter_winds: solver state + neighbours

        ims = domain%ims; ime = domain%ime; jms = domain%jms; jme = domain%jme
        kms = domain%kms; kme = domain%kme
        its = domain%its; ite = domain%ite; jts = domain%jts; jte = domain%jte

        allocate(div(ims:ime, kms:kme, jms:jme)); div = 0.0

        ! The iterative solver operates on the forcing-tendency winds (dqdt_3d), which
        ! are auto-allocated only for variables flagged as forcing inputs (i.e. when
        ! the forcing %uvar/%vvar/%wvar names are set). This test drives the solver
        ! directly with no forcing read, so allocate u/v/w dqdt_3d to mirror data_3d.
        call ensure_dqdt(domain%var_indx(kVARS%u)%v)
        call ensure_dqdt(domain%var_indx(kVARS%v)%v)
        call ensure_dqdt(domain%var_indx(kVARS%w)%v)

        ! --- seed a smooth, tile-continuous, DIVERGENT forcing wind in dqdt_3d. HICAR
        !     memory indices are global, so a function of (i,j) is automatically
        !     continuous across tile boundaries and the halo cells are well-defined. --
        call seed_divergent_winds()
        domain%vars_3d(domain%var_indx(kVARS%density)%v)%data_3d    = 1.0
        domain%vars_3d(domain%var_indx(kVARS%wind_alpha)%v)%data_3d = 1.0

        !$acc data copy(kVARS) copy(div)
        ! divergence of the seeded field
        call calc_divergence(div, domain, advect_density=.False., horz_only=.False., use_dqdt=.True.)
        !$acc update host(div)
        div0_max = maxval(abs(div(its:ite, kms:kme, jts:jte)))

        ! THE solve: drives bicgstab_solve -> exchange_krylov_halos on every iteration
        call calc_iter_winds(domain, &
            domain%vars_3d(domain%var_indx(kVARS%wind_alpha)%v)%data_3d, div, .False.)

        ! divergence of the corrected field
        call calc_divergence(div, domain, advect_density=.False., horz_only=.False., use_dqdt=.True.)
        !$acc update host(div)
        div1_max = maxval(abs(div(its:ite, kms:kme, jts:jte)))
        !$acc end data

        ! reduce to a global max so the assertion is decomposition-independent
        call MPI_Allreduce(div0_max, div0_max_g, 1, MPI_REAL, MPI_MAX, domain%compute_comms, ierr)
        call MPI_Allreduce(div1_max, div1_max_g, 1, MPI_REAL, MPI_MAX, domain%compute_comms, ierr)

        ! --- evaluate BEFORE tearing down (so cleanup always runs) -------------------
        ok = .True.; msg = ''
        if (div1_max_g /= div1_max_g) then                       ! NaN
            ok = .False.
            msg = 'corrected wind divergence is NaN'
        else if (.not. (div1_max_g < div0_max_g)) then
            ! The single uncalibrated analytic-operator solve should still reduce the
            ! peak divergence. If this proves marginal in CI, relax to a factor
            ! (e.g. div1 < 0.9*div0) rather than removing it -- a broken halo crashes
            ! under -fcheck before reaching here anyway.
            ok = .False.
            write(msg,'(A,ES12.4,A,ES12.4)') &
                'iterative solver did not reduce divergence: initial=', div0_max_g, &
                ' final=', div1_max_g
        endif

        ! reset the solver's module-level state so it does not leak into other suites
        call finalize_iter_winds()
        call domain%release()

        if (.not. ok) call test_failed(error, "test_iter_wind_solve", trim(msg))

    contains

        !> Allocate a variable's dqdt_3d to exactly mirror its data_3d bounds (incl.
        !> the staggered +1 face) and zero it, if not already allocated.
        subroutine ensure_dqdt(vidx)
            integer, intent(in) :: vidx
            integer :: l1,u1,l2,u2,l3,u3
            if (.not. allocated(domain%vars_3d(vidx)%dqdt_3d)) then
                l1 = lbound(domain%vars_3d(vidx)%data_3d,1); u1 = ubound(domain%vars_3d(vidx)%data_3d,1)
                l2 = lbound(domain%vars_3d(vidx)%data_3d,2); u2 = ubound(domain%vars_3d(vidx)%data_3d,2)
                l3 = lbound(domain%vars_3d(vidx)%data_3d,3); u3 = ubound(domain%vars_3d(vidx)%data_3d,3)
                allocate(domain%vars_3d(vidx)%dqdt_3d(l1:u1, l2:u2, l3:u3))
            endif
            domain%vars_3d(vidx)%dqdt_3d = 0.0
        end subroutine ensure_dqdt

        !> Fill u/v/w dqdt with a smooth, globally-continuous field whose horizontal
        !> divergence is non-zero. Loop over each array's OWN bounds (the staggered
        !> u/v carry an extra face) so the seeding itself is bounds-clean under -fcheck.
        subroutine seed_divergent_winds()
            integer :: ii, jj, kk
            real    :: lx, ly

            lx = real(max(domain%grid%ide - domain%grid%ids, 1))
            ly = real(max(domain%grid%jde - domain%grid%jds, 1))

            associate(ud => domain%vars_3d(domain%var_indx(kVARS%u)%v)%dqdt_3d)
                do jj = lbound(ud,3), ubound(ud,3)
                    do kk = lbound(ud,2), ubound(ud,2)
                        do ii = lbound(ud,1), ubound(ud,1)
                            ud(ii,kk,jj) = 8.0 * cos(2.0*PI*real(ii - domain%grid%ids)/lx)
                        enddo
                    enddo
                enddo
            end associate

            associate(vd => domain%vars_3d(domain%var_indx(kVARS%v)%v)%dqdt_3d)
                do jj = lbound(vd,3), ubound(vd,3)
                    do kk = lbound(vd,2), ubound(vd,2)
                        do ii = lbound(vd,1), ubound(vd,1)
                            vd(ii,kk,jj) = 6.0 * sin(2.0*PI*real(jj - domain%grid%jds)/ly)
                        enddo
                    enddo
                enddo
            end associate

            domain%vars_3d(domain%var_indx(kVARS%w)%v)%dqdt_3d = 0.0
        end subroutine seed_divergent_winds

    end subroutine test_iter_wind_solve

end module test_wind_iterative
