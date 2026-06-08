! test_snow_drift.F90
! Decomposition-reproducibility test for the FINE-MESH advection used by
! src/physics/snow_drift.F90 (blowing-snow suspension transport).
!
! The fine-mesh horizontal advection (flux_2d_fm -> WRF_flux_corr_fm ->
! sum_kernel_2d_fm, with the exch_fine_mesh_3d halo exchange) is the only
! MPI-coupled part of snow_drift; the vertical Thomas solver is column-local.
! This test replicates that exact call sequence on a NON-CONSTANT field and
! prints qs_fm at fixed GLOBAL interior points ("FMPROBE i j k val"). Running
! the HICAR-tester at two process counts and diffing those lines (see
! tests/test_fm_decomposition.sh) checks bit-for-bit decomposition
! reproducibility — a constant field would mask halo/FCT bugs, so a varying
! field is used deliberately.
module test_snow_drift

    use variable_interface,  only: variable_t
    use mpi_f08
    use icar_constants
    use testdrive,           only: new_unittest, unittest_type, error_type, check, test_failed, get_argument
    use string,              only: to_lower
    use domain_interface,    only: domain_t
    use options_interface,   only: options_t
    use io_routines,         only: check_file_exists
    use time_step,           only: compute_dt
    use grid_interface,      only: grid_t
    use data_structures,     only: index_type
    use snow_drift,          only: snow_drift_var_request
    use advection,           only: adv_init, adv_var_request
    use adv_std,             only: adv_std_compute_wind_2d_fm, flux_2d_fm, sum_kernel_2d_fm, &
                                   adv_std_clean_wind_arrays_fm, flux_x_fm, flux_y_fm
    use adv_fluxcorr,        only: init_fluxcorr_fm, set_sign_arrays_fm, WRF_flux_corr_fm
    implicit none
    private

    public :: collect_snow_drift_suite

contains

    subroutine collect_snow_drift_suite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
        testsuite = [ &
            new_unittest("fm_advect_decomp", test_fm_advect_decomp) &
            ]
    end subroutine collect_snow_drift_suite

    subroutine test_fm_advect_decomp(error)
        type(error_type), allocatable, intent(out) :: error

        type(domain_t)  :: domain
        type(options_t) :: options
        type(index_type), allocatable :: exch_vars(:)
        real, allocatable :: u_fm(:,:,:), v_fm(:,:,:), rho_fm(:,:,:)
        real, allocatable :: jaco_fm(:,:,:), jaco_u_fm(:,:,:), jaco_v_fm(:,:,:)
        real, allocatable :: U_m_fm(:,:,:), V_m_fm(:,:,:), denom_fm(:,:,:)
        real, allocatable :: qs_fm_old(:,:,:)
        integer :: N, qsv, ii, ig, jg, kk, iter, max_iters, fc, my_rank, ierr
        integer :: ims,ime,jms,jme,its,ite,jts,jte
        real :: dx, dt, t_fac, gmax, gmin
        logical :: probe_output
        character(len=:), allocatable :: arg1

        N = 6                 ! fine-mesh levels (small for speed; >=3 exercises the k-loop)
        kFM_GRID_Z = N        ! qs_fm/ns_fm are allocated on grid_fm with this many levels
        call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)

        ! The per-cell FMPROBE dump is only consumed by the decomposition gate;
        ! for a normal unit-test run it is just noise (thousands of lines). Emit it
        ! only in verbose mode — i.e. when "-v" is the first argument to the test
        ! driver (the same flag test_driver.F90 keys verbose output off of). Detect
        ! it directly from the command line rather than via the global STD_OUT_PE,
        ! which other suites mutate. The pass/fail check below runs either way.
        call get_argument(1, arg1)
        probe_output = allocated(arg1)
        if (probe_output) probe_output = (to_lower(arg1) == "-v")

        call options%init()
        options%domain%init_conditions_file = '../tests/Test_Cases/domains/flat_plane_250m.nc'
        options%domain%hgt_hi='topo'; options%domain%lat_hi='lat'; options%domain%lon_hi='lon'
        options%domain%dx = 250.0
        options%domain%nz = 10
        options%domain%dz_levels = 50.0
        options%domain%sleve = .True.
        options%physics%advection = kADV_STD
        options%sm%suspension_fine_mesh_levels = N
        ! 3rd-order + monotonic flux limiter + RK3: the production fine-mesh config
        options%adv%h_order = 3
        options%adv%v_order = 3
        options%adv%flux_corr = kFLUXCOR_MONO
        options%time%RK3 = .true.
        options%adv%advect_density = .False.

        domain%compute_comms = MPI_COMM_WORLD
        call snow_drift_var_request(options)
        call adv_var_request(options)
        call domain%init(options,1)
        call adv_init(domain,options)

        ims=domain%ims; ime=domain%ime; jms=domain%jms; jme=domain%jme
        its=domain%its; ite=domain%ite; jts=domain%jts; jte=domain%jte
        dx = domain%dx
        dt = 15.0
        qsv = domain%var_indx(kVARS%qs_fm)%v

        ! halo-exchange descriptor for qs_fm (mirrors snow_drift's exch_vars_fm)
        allocate(exch_vars(1))
        exch_vars(1)%v  = qsv
        exch_vars(1)%id = domain%vars_3d(qsv)%id

        ! fine-mesh scratch (snow_drift keeps these module-level; the test owns them)
        allocate(u_fm(ims:ime,N,jms:jme), v_fm(ims:ime,N,jms:jme), rho_fm(ims:ime,N,jms:jme))
        allocate(jaco_fm(ims:ime,N,jms:jme), jaco_u_fm(ims:ime,N,jms:jme), jaco_v_fm(ims:ime,N,jms:jme))
        allocate(denom_fm(ims:ime,N,jms:jme), qs_fm_old(ims:ime,N,jms:jme))
        ! U_m_fm / V_m_fm are now allocated inside adv_std_compute_wind_2d_fm on its
        ! first call (mirroring the regular scheme); pass them in unallocated.
        rho_fm = 1.0; jaco_fm = 1.0; jaco_u_fm = 1.0; jaco_v_fm = 1.0
        u_fm = 8.0; v_fm = -6.0      ! constant => divergence-free on the flat fine mesh

        ! SHARP blob (steep gradients + near-zero background) so the FCT limiter
        ! actively clips. Set over full memory (incl halos) so domain-edge cells
        ! act as fixed inflow; set on EVERY advected fine-mesh level.
        associate(q => domain%vars_3d(qsv)%data_3d)
            do jg = lbound(q,3), ubound(q,3)
                do kk = 1, N
                    do ig = lbound(q,1), ubound(q,1)
                        q(ig,kk,jg) = 1.0e-6 + 2.0e-2 &
                            * exp(-((real(ig)-72.0)**2 + (real(jg)-64.0)**2)/(2.0*36.0)) &
                            * (1.0 + 0.1*real(kk))
                    enddo
                enddo
            enddo
        end associate
        if (my_rank == 0 .and. probe_output) &
            write(*,'(A,4(I0,1X))') "FMRANK0 its ite jts jte = ", its, ite, jts, jte

        !$acc data copy(kVARS)
        ! Advance the fine-mesh advection several steps (each step == one
        ! snow_drift_integrate horizontal pass). Constant winds => the blob
        ! translates diagonally across interior tile boundaries.
        do ii = 1, 8
            ! Exchange the halo BEFORE snapshotting qs_fm_old, exactly as
            ! snow_drift_integrate does (exch_fine_mesh_3d precedes the
            ! qs_fm_old = qs_fm save). The FCT limiter reads qs_fm_old's halo
            ! (compute_upwind_fluxes_fm reaches q(its-3) at a tile boundary), so a
            ! stale qs_fm_old halo would make the limiter decomposition-dependent.
            call domain%halo%halo_3d_send_batch(exch_vars, domain%vars_3d)
            call domain%halo%halo_3d_retrieve_batch(exch_vars, domain%vars_3d)
            associate(q => domain%vars_3d(qsv)%data_3d)
                qs_fm_old(:,:,:) = q(:,:,:)
            end associate
            call adv_std_compute_wind_2d_fm(u_fm, v_fm, rho_fm, jaco_fm, jaco_u_fm, jaco_v_fm, &
                dx, dt, U_m_fm, V_m_fm, denom_fm, 1, N)
            if (options%adv%flux_corr == kFLUXCOR_MONO) then
                call init_fluxcorr_fm(1, N)
                call set_sign_arrays_fm(U_m_fm, V_m_fm, 1, N)
            endif

            max_iters = 3
            do iter = 1, max_iters
                select case(iter)
                case (1); t_fac = 1.0/3.0; fc = 0
                case (2); t_fac = 0.5;     fc = 0
                case (3); t_fac = 1.0;     fc = options%adv%flux_corr
                end select
                if (iter > 1) then
                    call domain%halo%halo_3d_send_batch(exch_vars, domain%vars_3d)
                    call domain%halo%halo_3d_retrieve_batch(exch_vars, domain%vars_3d)
                endif
                call flux_2d_fm(domain%vars_3d(qsv)%data_3d, U_m_fm, V_m_fm, t_fac, 1, N)
                if (fc == kFLUXCOR_MONO) then
                    call WRF_flux_corr_fm(qs_fm_old, U_m_fm, V_m_fm, flux_x_fm, flux_y_fm, denom_fm, 1, N)
                endif
                call sum_kernel_2d_fm(qs_fm_old, domain%vars_3d(qsv)%data_3d, denom_fm, 1, N)
            enddo
        end do
        ! Fine-mesh wind/flux arrays are allocated once (first compute_wind call)
        ! and reused across every step; free them once here at the end.
        call adv_std_clean_wind_arrays_fm()
        !$acc end data

        ! probe qs_fm at fixed GLOBAL interior points (whichever rank owns them).
        ! Only emitted for the decomposition gate (HICAR_FM_PROBE=1); a normal
        ! unit-test run skips the dump and just evaluates the pass/fail check.
        associate(q => domain%vars_3d(qsv)%data_3d)
            kk = 2
            if (probe_output) then
                do jg = 5, 90
                    do ig = 5, 95
                        if (ig>=its .and. ig<=ite .and. jg>=jts .and. jg<=jte) then
                            write(*,'(A,I3,1X,I3,1X,I2,1X,ES24.16)') "FMPROBE ", ig, jg, kk, q(ig,kk,jg)
                        endif
                    enddo
                enddo
            endif
            gmax = maxval(q(its:ite,1:N,jts:jte))
            gmin = minval(q(its:ite,1:N,jts:jte))
        end associate

        ! testdrive pass/fail: scheme must stay finite and (monotonic FCT) bounded
        call check(error, (gmax == gmax) .and. (gmax < 1.0) .and. (gmin > -1.0e-9))

        call domain%release()
    end subroutine test_fm_advect_decomp

end module test_snow_drift
