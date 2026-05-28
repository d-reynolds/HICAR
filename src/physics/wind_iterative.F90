!>------------------------------------------------------------
!! Native HICAR iterative wind solver.
!!
!! Right-preconditioned BiCGStab + Block-Jacobi (diagonal scaling)
!! over a 15-point stencil with terrain-following boundary
!! conditions. Pure OpenACC + MPI — vendor-portable, no external
!! solver library.
!!
!! Phase 1 (this file): classic BiCGStab, MPI_Allreduce on dot
!! products, MPI Sendrecv halos, pure-stencil SpMV (no CSR).
!! Multi-nest cache and pipelined/NCCL paths land in later phases.
!!
!! Matrix structure mirrors wind_iterative_amgx.F90 exactly so the
!! two solvers are interchangeable for validation.
!!
!!  @author
!!  Dylan Reynolds (dylan.reynolds@slf.ch)
!!------------------------------------------------------------

module wind_iterative
    use iso_c_binding
    use domain_interface,  only : domain_t
    use icar_constants,    only : STD_OUT_PE, kVARS
    use options_interface, only : options_t
    use iso_fortran_env
    use mpi_f08
    use openacc
    use string,        only : str
    use debug_module,  only : domain_check_winds
#ifdef USE_NCCL
    use nccl_interface, only : nccl_comm_init, nccl_comm_destroy, &
                               nccl_group_start, nccl_group_end, &
                               nccl_send_double, nccl_recv_double, &
                               nccl_allreduce_double_sum
#endif

    implicit none
    private
    public :: init_iter_winds, calc_iter_winds, finalize_iter_winds

    logical :: initialized_iter_winds = .false.
    logical :: structure_uploaded = .false.

    real, parameter :: deg2rad = 0.017453293
    real, parameter :: rad2deg = 57.2957779371

    integer :: wind_solver_max_iters = 1500
    real(c_double), parameter :: bicg_tol_abs = 1.0e-10_c_double
    real(c_double), parameter :: bicg_tol_rel = 1.0e-5_c_double
    real(c_double), parameter :: breakdown_eps = 1.0e-30_c_double

    ! Number of Richardson sweeps per Block-Jacobi preconditioner apply.
    ! AMGX's PBICGSTAB+BLOCK_JACOBI config uses prec:max_iters=2 — matching that
    ! gives ~2x stronger preconditioning per outer iter at ~2x per-iter cost,
    ! but the *outer* iter count typically drops 1.5-2x for Poisson-like problems.
    ! Net: similar total time, but fewer outer iters means fewer halo exchanges
    ! and allreduces — strict win at multi-rank scale.
    !
    ! Adaptive retry: if a solve diverges/stagnates, calc_iter_winds bumps
    ! this up to MAX_PREC_SWEEPS, retrying after each bump. On successful
    ! convergence we reset to BASE_PREC_SWEEPS for the next call. Mirrors the
    ! AMGX module's prec_max_iters retry. Per-nest state is cached.
    integer, parameter :: BASE_PREC_SWEEPS = 2
    integer, parameter :: MAX_PREC_SWEEPS  = 4
    integer :: precond_n_sweeps = BASE_PREC_SWEEPS

    ! 15-point stencil coefficients (same names as AMGX module)
    real, allocatable, dimension(:,:,:) :: A_coef, B_coef, C_coef, D_coef, E_coef, F_coef, G_coef, &
                                           H_coef, I_coef, J_coef, K_coef, L_coef, M_coef, N_coef, O_coef
    real    :: dx
    real, allocatable, dimension(:,:,:) :: div, dz_if, jaco, dzdx, dzdy, sigma, alpha
    real, allocatable, dimension(:,:)   :: dzdx_surf, dzdy_surf
    real, allocatable, dimension(:,:,:) :: jaco_w, dzdx_u_stag, jaco_u_stag, dzdy_v_stag, jaco_v_stag
    real, allocatable, dimension(:)     :: adv_dz_col

    ! Domain decomposition / boundary-cell scalars (PETSc/AMGX convention)
    integer :: hs, i_s, i_e, k_s, k_e, j_s, j_e
    integer :: ims, ime, jms, jme, ids, ide, jds, jde
    integer :: xs, ys, zs, xm, ym, zm, mx, my, mz
    integer :: n_rows, n_rows_global

    ! Krylov / RHS / preconditioner state — 3D with halos so SpMV reads neighbours directly.
    ! Allocation is (i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1) which covers:
    !   - the owned cell range xs..xs+xm-1 (incl. BC layer at boundary ranks),
    !   - one lateral halo cell in i and j on non-boundary edges.
    !
    ! Classic right-preconditioned BiCGStab vectors. Earlier Phase-2 versions tracked
    ! u = M^{-1} r and q = M^{-1} v via recurrences to skip fresh preconditioner applies,
    ! but the recurrences accumulated roundoff and could stall convergence at high iter
    ! counts. Reverted to fresh M^{-1} on p and s. Merged 5-value 2nd allreduce is kept
    ! so we still pay only 2 allreduces per iter.
    real(c_double), allocatable, dimension(:,:,:) :: x_sol      ! solution
    real(c_double), allocatable, dimension(:,:,:) :: r_vec      ! residual
    real(c_double), allocatable, dimension(:,:,:) :: r_hat      ! shadow residual r̂_0
    real(c_double), allocatable, dimension(:,:,:) :: p_vec, p_hat
    real(c_double), allocatable, dimension(:,:,:) :: v_vec
    real(c_double), allocatable, dimension(:,:,:) :: s_vec, s_hat
    real(c_double), allocatable, dimension(:,:,:) :: t_vec
    real(c_double), allocatable, dimension(:,:,:) :: rhs        ! right-hand side b
    real(c_double), allocatable, dimension(:,:,:) :: D_inv      ! Block-Jacobi: 1/diag(A) per row type
    real(c_double), allocatable, dimension(:,:,:) :: prec_res   ! scratch for multi-sweep Richardson residual

    ! Persistent halo face buffers — allocated once at first solve, reused every iter.
    ! 8 buffers (4 directions x send/recv) so all 4 sends and all 4 recvs can be in flight
    ! simultaneously via Isend/Irecv + Waitall. Cuts halo latency from 4*alpha (sequential
    ! Sendrecvs) to ~1*alpha — major win at high MPI rank counts.
    real(c_double), allocatable, target, dimension(:,:) :: east_send,  east_recv     ! shape (k_s-1:k_e+1, j_s-1:j_e+1)
    real(c_double), allocatable, target, dimension(:,:) :: west_send,  west_recv     ! shape (k_s-1:k_e+1, j_s-1:j_e+1)
    real(c_double), allocatable, target, dimension(:,:) :: north_send, north_recv    ! shape (i_s-1:i_e+1, k_s-1:k_e+1)
    real(c_double), allocatable, target, dimension(:,:) :: south_send, south_recv    ! shape (i_s-1:i_e+1, k_s-1:k_e+1)

    type(MPI_Comm) :: solver_comm
    integer :: solver_rank = -1
    integer :: east_neighbor = -1, west_neighbor = -1, north_neighbor = -1, south_neighbor = -1

    ! Multi-nest cache: each nest's full state lives in its own slot. On nest
    ! context switch (nest_manager.F90:switch_nest_context), init_iter_winds
    ! saves the active state and restores the target. O(1) switching once each
    ! nest has been visited once. Mirrors the AMGX module's cache (which is
    ! itself based on the long-standing physics module pattern in HICAR).
    integer, parameter :: MAX_NESTS = 4
    integer :: active_nest_indx = -1

    type :: hicar_cache_t
        logical :: valid = .false.
        ! Grid scalars
        integer :: i_s, i_e, k_s, k_e, j_s, j_e
        integer :: ims, ime, jms, jme, ids, ide, jds, jde
        integer :: xs, ys, zs, xm, ym, zm, mx, my, mz
        integer :: hs, n_rows, n_rows_global
        real    :: dx
        integer :: solver_rank, east_neighbor, west_neighbor, north_neighbor, south_neighbor
        type(MPI_Comm) :: solver_comm
        ! Flags
        logical :: structure_uploaded
        integer :: wind_solver_max_iters
        integer :: precond_n_sweeps = BASE_PREC_SWEEPS
        ! Geometry (single precision, owned-range allocations)
        real, allocatable, dimension(:,:,:) :: A_coef, B_coef, C_coef, D_coef, E_coef
        real, allocatable, dimension(:,:,:) :: F_coef, G_coef, H_coef, I_coef, J_coef
        real, allocatable, dimension(:,:,:) :: K_coef, L_coef, M_coef, N_coef, O_coef
        real, allocatable, dimension(:,:,:) :: div, dz_if, jaco, dzdx, dzdy, sigma, alpha
        real, allocatable, dimension(:,:)   :: dzdx_surf, dzdy_surf
        real, allocatable, dimension(:,:,:) :: jaco_w, dzdx_u_stag, jaco_u_stag, dzdy_v_stag, jaco_v_stag
        real, allocatable, dimension(:)     :: adv_dz_col
        ! Krylov / RHS / precond (double precision, halo'd allocations)
        real(c_double), allocatable, dimension(:,:,:) :: x_sol, r_vec, r_hat
        real(c_double), allocatable, dimension(:,:,:) :: p_vec, p_hat, v_vec
        real(c_double), allocatable, dimension(:,:,:) :: s_vec, s_hat, t_vec
        real(c_double), allocatable, dimension(:,:,:) :: rhs, D_inv, prec_res
        ! Halo face buffers (target attribute only needed on the live module-level
        ! variables for c_loc — not on cache slots, which are pure host storage)
        real(c_double), allocatable, dimension(:,:) :: east_send,  east_recv
        real(c_double), allocatable, dimension(:,:) :: west_send,  west_recv
        real(c_double), allocatable, dimension(:,:) :: north_send, north_recv
        real(c_double), allocatable, dimension(:,:) :: south_send, south_recv
#ifdef USE_NCCL
        real(c_double), allocatable, dimension(:) :: sigma_dev, red5_dev, rho0_dev
#endif
    end type hicar_cache_t

    type(hicar_cache_t) :: domain_cache(MAX_NESTS)

#ifdef USE_NCCL
    ! NCCL communicator + stream. Stream is aliased to OpenACC's sync queue so
    ! all NCCL ops are stream-ordered with the OpenACC pack/unpack/reduction kernels
    ! that produce inputs and consume outputs — no explicit cross-stream sync needed.
    type(c_ptr) :: nccl_comm   = c_null_ptr
    type(c_ptr) :: nccl_stream = c_null_ptr
    logical     :: nccl_initialized = .false.
    ! Device-resident allreduce buffers (one for sigma, one for the 5-value reduction).
    ! Allocated once at first solve so NCCL operates entirely on device pointers.
    real(c_double), allocatable, target, dimension(:) :: sigma_dev   ! length 1
    real(c_double), allocatable, target, dimension(:) :: red5_dev    ! length 5
    real(c_double), allocatable, target, dimension(:) :: rho0_dev    ! length 2 (rho_0 + ||r0||²)
#endif

    ! Per-solve timing accumulators (set to 0 at start of each bicgstab_solve, summary
    ! printed at end). MPI_Wtime is portable across CPU and GPU paths. With synchronous
    ! OpenACC kernels, MPI_Wtime around a kernel call measures the kernel's wall time.
    real(c_double) :: t_spmv_acc      = 0.0_c_double
    real(c_double) :: t_halo_acc      = 0.0_c_double
    real(c_double) :: t_precond_acc   = 0.0_c_double
    real(c_double) :: t_allreduce_acc = 0.0_c_double
    real(c_double) :: t_vecops_acc    = 0.0_c_double
    real(c_double) :: t_total_acc     = 0.0_c_double

contains

    !>------------------------------------------------------------
    !! Per-nest initialisation — called once per nest at startup.
    !! Phase 1: no multi-nest cache, so this just allocates state
    !! the first time and otherwise returns. Calling for a second
    !! nest will currently overwrite state for the first.
    !!------------------------------------------------------------
    subroutine init_iter_winds(domain, options)
        implicit none
        type(domain_t),  intent(in) :: domain
        type(options_t), intent(in) :: options
        integer :: ierr
        integer :: nprocs, device_num
        integer :: target_nest
        integer :: my_rank_in_comm
#ifdef USE_NCCL
        integer(c_int) :: nccl_rc
#endif

        ! NCCL communicator is shared across nests (created once on first call).
        ! It's tied to compute_comms which we assume is the same across nests.
#ifdef USE_NCCL
        if (.not. nccl_initialized) then
            call MPI_Comm_size(domain%compute_comms, nprocs,         ierr)
            call MPI_Comm_rank(domain%compute_comms, my_rank_in_comm, ierr)
            device_num = acc_get_device_num(acc_device_nvidia)
            nccl_rc = nccl_comm_init(nccl_comm, nprocs, my_rank_in_comm, &
                                      domain%compute_comms%MPI_VAL, device_num)
            if (nccl_rc /= 0 .and. STD_OUT_PE) then
                print*, "WARNING: nccl_comm_init failed in wind_iterative, rc=", nccl_rc
            endif
            nccl_stream = transfer(acc_get_cuda_stream(acc_async_sync), nccl_stream)
            nccl_initialized = .true.
        endif
#endif

        target_nest = domain%nest_indx

        ! Already on this nest — nothing to do
        if (target_nest == active_nest_indx) return

        ! Save current state to its cache slot (if any)
        if (active_nest_indx > 0 .and. active_nest_indx <= MAX_NESTS) then
            call save_to_cache(active_nest_indx)
        endif

        ! Restore from cache if we've seen this nest before — O(1) context switch
        if (target_nest > 0 .and. target_nest <= MAX_NESTS .and. &
            domain_cache(target_nest)%valid) then
            call restore_from_cache(target_nest)
            active_nest_indx = target_nest
            return
        endif

        ! --- Fresh init for a nest seen for the first time ---
        solver_comm = domain%compute_comms

        ! Cache MPI rank + neighbour ranks (exchange_krylov_halos hits these ~200x/solve)
        call MPI_Comm_rank(solver_comm, solver_rank, ierr)
        east_neighbor  = solver_rank + 1
        west_neighbor  = solver_rank - 1
        north_neighbor = solver_rank + domain%grid%ximages
        south_neighbor = solver_rank - domain%grid%ximages

        wind_solver_max_iters = options%wind%wind_solver_iterations

        call init_module_vars(domain)

        n_rows_global = mx * my * mz
        n_rows        = xm * ym * zm

        initialized_iter_winds  = .true.
        structure_uploaded = .false.
        active_nest_indx   = target_nest

        if (STD_OUT_PE) print*, "HICAR native wind solver initialised for nest ", target_nest
    end subroutine init_iter_winds


    !>------------------------------------------------------------
    !! Main entry: run BiCGStab solve and apply the wind correction.
    !! Signature matches calc_iter_winds_amgx so wind.F90 dispatch
    !! is a one-line swap.
    !!------------------------------------------------------------
    subroutine calc_iter_winds(domain, alpha_in, div_in, adv_den)
        implicit none
        type(domain_t), intent(inout) :: domain
        real, dimension(ims:ime, domain%kms:domain%kme, jms:jme), intent(in) :: alpha_in, div_in
        logical, intent(in) :: adv_den

        integer :: i, j, k
        real    :: alpha_min, alpha_max
        logical :: varying_alpha
        integer :: status, n_iters
        integer :: nan_count
        real(c_double) :: res0, res_final

        ! Copy alpha and div into module-resident arrays on GPU
        !$acc parallel loop gang vector collapse(3) present(alpha, div, alpha_in, div_in)
        do j = j_s, j_e
            do k = k_s, k_e
                do i = i_s, i_e
                    div(i,k,j)   = div_in(i,k,j)
                    alpha(i,k,j) = alpha_in(i,k,j)
                enddo
            enddo
        enddo

        ! Debug: if the input divergence contains any NaN, dump the min/max of the
        ! state fields that feed the divergence so the source can be traced.
        nan_count = 0
        !$acc parallel loop gang vector collapse(3) reduction(+:nan_count) present(div)
        do j = j_s, j_e
            do k = k_s, k_e
                do i = i_s, i_e
                    if (div(i,k,j) /= div(i,k,j)) nan_count = nan_count + 1
                enddo
            enddo
        enddo
        if (nan_count > 0) then
            associate(density     => domain%vars_3d(domain%var_indx(kVARS%density)%v)%data_3d,     &
                      pressure    => domain%vars_3d(domain%var_indx(kVARS%pressure)%v)%data_3d,    &
                      temperature => domain%vars_3d(domain%var_indx(kVARS%temperature)%v)%data_3d, &
                      qv          => domain%vars_3d(domain%var_indx(kVARS%water_vapor)%v)%data_3d,  &
                      u           => domain%vars_3d(domain%var_indx(kVARS%u)%v)%dqdt_3d,            &
                      v           => domain%vars_3d(domain%var_indx(kVARS%v)%v)%dqdt_3d,            &
                      w_grid      => domain%vars_3d(domain%var_indx(kVARS%w)%v)%dqdt_3d)
                !$acc update host(density, pressure, temperature, qv, u, v, w_grid)
                write(*,*) '--------------- HICAR wind solver: NaN in input divergence (rank=', &
                           solver_rank, ', count=', nan_count, ') ---------------'
                write(*,*) "  density:     min=", minval(density),     " max=", maxval(density)
                write(*,*) "  pressure:    min=", minval(pressure),    " max=", maxval(pressure)
                write(*,*) "  temperature: min=", minval(temperature), " max=", maxval(temperature)
                write(*,*) "  qv:          min=", minval(qv),          " max=", maxval(qv)
                write(*,*) "  u:           min=", minval(u),           " max=", maxval(u)
                write(*,*) "  v:           min=", minval(v),           " max=", maxval(v)
                write(*,*) "  w:           min=", minval(w_grid),      " max=", maxval(w_grid)
            end associate
        endif

        if (.not. structure_uploaded) then
            ! First call: build coefficient arrays (CPU host needs alpha)
            !$acc update host(alpha)

            call initialize_coefs(domain)

            ! Allocate Krylov / RHS / D_inv state on host then push to device.
            ! Range matches AMGX's lambda_3d: (i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1).
            allocate(x_sol  (i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1))
            allocate(r_vec  (i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1))
            allocate(r_hat  (i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1))
            allocate(p_vec  (i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1))
            allocate(p_hat  (i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1))
            allocate(v_vec  (i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1))
            allocate(s_vec  (i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1))
            allocate(s_hat  (i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1))
            allocate(t_vec  (i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1))
            allocate(rhs     (i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1))
            allocate(D_inv   (i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1))
            allocate(prec_res(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1))
#ifdef USE_NCCL
            ! Device-resident scalar buffers for NCCL allreduces.
            allocate(sigma_dev(1)); sigma_dev = 0.0_c_double
            allocate(red5_dev (5)); red5_dev  = 0.0_c_double
            allocate(rho0_dev (2)); rho0_dev  = 0.0_c_double
            !$acc enter data copyin(sigma_dev, red5_dev, rho0_dev)
#endif

            ! Persistent halo face buffers — 8 separate buffers so all 4 directions can be
            ! in flight at once (Isend/Irecv + Waitall pattern).
            allocate(east_send (k_s-1:k_e+1, j_s-1:j_e+1))
            allocate(east_recv (k_s-1:k_e+1, j_s-1:j_e+1))
            allocate(west_send (k_s-1:k_e+1, j_s-1:j_e+1))
            allocate(west_recv (k_s-1:k_e+1, j_s-1:j_e+1))
            allocate(north_send(i_s-1:i_e+1, k_s-1:k_e+1))
            allocate(north_recv(i_s-1:i_e+1, k_s-1:k_e+1))
            allocate(south_send(i_s-1:i_e+1, k_s-1:k_e+1))
            allocate(south_recv(i_s-1:i_e+1, k_s-1:k_e+1))

            ! Zero everything to make halo cells well-defined
            x_sol = 0.0_c_double; r_vec = 0.0_c_double; r_hat = 0.0_c_double
            p_vec = 0.0_c_double; p_hat = 0.0_c_double; v_vec = 0.0_c_double
            s_vec = 0.0_c_double; s_hat = 0.0_c_double; t_vec = 0.0_c_double
            rhs   = 0.0_c_double; D_inv = 1.0_c_double; prec_res = 0.0_c_double
            east_send  = 0.0_c_double; east_recv  = 0.0_c_double
            west_send  = 0.0_c_double; west_recv  = 0.0_c_double
            north_send = 0.0_c_double; north_recv = 0.0_c_double
            south_send = 0.0_c_double; south_recv = 0.0_c_double

            !$acc enter data copyin(x_sol, r_vec, r_hat, p_vec, p_hat, v_vec, &
            !$acc                   s_vec, s_hat, t_vec, rhs, D_inv, prec_res, &
            !$acc                   east_send, east_recv, west_send, west_recv, &
            !$acc                   north_send, north_recv, south_send, south_recv)
            !$acc enter data copyin(A_coef, B_coef, C_coef, D_coef, E_coef, F_coef, G_coef, &
            !$acc                   H_coef, I_coef, J_coef, K_coef, L_coef, M_coef, N_coef, O_coef)

            structure_uploaded = .true.
        else
            ! Subsequent call: detect varying alpha on GPU; refresh coefs only when needed
            alpha_min =  HUGE(1.0)
            alpha_max = -HUGE(1.0)
            !$acc parallel loop gang vector collapse(3) reduction(min:alpha_min) reduction(max:alpha_max) present(alpha)
            do j = j_s, j_e
                do k = k_s, k_e
                    do i = i_s, i_e
                        alpha_min = min(alpha_min, alpha(i,k,j))
                        alpha_max = max(alpha_max, alpha(i,k,j))
                    enddo
                enddo
            enddo
            varying_alpha = (alpha_max > alpha_min)
            if (varying_alpha) call update_coefs_gpu()
        endif

        ! Build / refresh the diagonal preconditioner from current coefficients
        call extract_diagonal()

        ! Build RHS on GPU (3D layout: rhs(i,k,j) = -2*div for interior, 0 at BCs)
        call compute_rhs_3d()

        ! Reset initial guess to avoid BiCGStab breakdown from stale state
        call vec_zero(x_sol)

        !$acc wait

        ! Solve A x = b at the current preconditioner sweep count
        call bicgstab_solve(domain, wind_solver_max_iters, status, n_iters, res0, res_final)

        if (STD_OUT_PE) then
            write(*,*) ' HICAR BiCGStab status=', status, ' iterations=', n_iters, &
                       ' precond_n_sweeps=', precond_n_sweeps
            write(*,*) '   Residual at iter 0:    ', res0
            write(*,*) '   Residual at final iter:', res_final
            if (t_total_acc > 0.0_c_double) then
                write(*,'(A)')          '   --- per-solve timing breakdown (rank 0, seconds) ---'
                write(*,'(A,F10.4,A)')  '     total           : ', t_total_acc, ''
                write(*,'(A,F10.4,A,F6.1,A)') '     SpMV            : ', t_spmv_acc, &
                    '   (', 100.0_c_double * t_spmv_acc / t_total_acc, ' %)'
                write(*,'(A,F10.4,A,F6.1,A)') '     halo exchange   : ', t_halo_acc, &
                    '   (', 100.0_c_double * t_halo_acc / t_total_acc, ' %)'
                write(*,'(A,F10.4,A,F6.1,A)') '     preconditioner  : ', t_precond_acc, &
                    '   (', 100.0_c_double * t_precond_acc / t_total_acc, ' %)'
                write(*,'(A,F10.4,A,F6.1,A)') '     allreduce       : ', t_allreduce_acc, &
                    '   (', 100.0_c_double * t_allreduce_acc / t_total_acc, ' %)'
                write(*,'(A,F10.4,A,F6.1,A)') '     vector ops      : ', t_vecops_acc, &
                    '   (', 100.0_c_double * t_vecops_acc / t_total_acc, ' %)'
                write(*,'(A,F10.4,A,F6.1,A)') '     other           : ', &
                    t_total_acc - (t_spmv_acc + t_halo_acc + t_precond_acc + t_allreduce_acc + t_vecops_acc), &
                    '   (', 100.0_c_double * (t_total_acc - (t_spmv_acc + t_halo_acc + t_precond_acc + t_allreduce_acc + t_vecops_acc)) / t_total_acc, ' %)'
                if (n_iters > 0) then
                    write(*,'(A,F10.6,A)')  '     per-iter total  : ', t_total_acc / real(n_iters, c_double), ' s'
                endif
            endif
        endif

        ! Adaptive preconditioner retry: if the solver diverged/stagnated
        ! (status non-zero AND residual didn't shrink by 100x), bump the
        ! Block-Jacobi sweep count and retry. Each retry restarts from
        ! x_sol = 0 (bicgstab_solve assumes x0=0 — see r_vec = rhs init).
        ! Mirrors wind_iterative_amgx.F90's prec_max_iters retry.
        if (status /= 0 .and. n_iters > 0 .and. res_final > 0.01_c_double * res0) then
            do while (status /= 0 .and. res_final > 0.01_c_double * res0 &
                      .and. precond_n_sweeps < MAX_PREC_SWEEPS)
                precond_n_sweeps = precond_n_sweeps + 1
                if (STD_OUT_PE) write(*,*) ' Convergence unsatisfactory, retrying solve with precond_n_sweeps=', &
                                           precond_n_sweeps

                ! Restart from x = 0 so bicgstab_solve's r0 = b - A*x0 = b assumption holds.
                call vec_zero(x_sol)
                !$acc wait

                call bicgstab_solve(domain, wind_solver_max_iters, status, n_iters, res0, res_final)

                if (STD_OUT_PE) then
                    write(*,*) '  Retry status=', status, ' iterations=', n_iters
                    write(*,*) '  Residual at iter 0:     ', res0
                    write(*,*) '  Residual at final iter: ', res_final
                endif
            end do

            ! Retry succeeded — reset sweep count to base for the next call
            if (status == 0 .or. res_final <= 0.01_c_double * res0) then
                if (STD_OUT_PE) write(*,*) ' Retry converged, resetting precond_n_sweeps=', BASE_PREC_SWEEPS
            endif
            precond_n_sweeps = BASE_PREC_SWEEPS

            ! Still diverged after max retries and residual is worse than start — fatal
            if (status /= 0 .and. res_final > res0) then
                if (STD_OUT_PE) then
                    associate(density => domain%vars_3d(domain%var_indx(kVARS%density)%v)%data_3d, &
                              u       => domain%vars_3d(domain%var_indx(kVARS%u)%v)%dqdt_3d, &
                              v       => domain%vars_3d(domain%var_indx(kVARS%v)%v)%dqdt_3d, &
                              w_grid  => domain%vars_3d(domain%var_indx(kVARS%w)%v)%dqdt_3d, &
                              w_real  => domain%vars_3d(domain%var_indx(kVARS%w_real)%v)%data_3d, &
                              alpha_d => domain%vars_3d(domain%var_indx(kVARS%wind_alpha)%v)%data_3d)
                        write(*,*) '---------------HICAR BiCGStab DIVERGED (precond_n_sweeps=', &
                                   precond_n_sweeps, ')-----------'
                        !$acc update host(density, u, v, w_grid, w_real, alpha_d, div)
                        write(*,*) "max abs val of density: ", maxval(abs(density))
                        write(*,*) "max abs val of u:       ", maxval(abs(u))
                        write(*,*) "max abs val of v:       ", maxval(abs(v))
                        write(*,*) "max abs val of w_grid:  ", maxval(abs(w_grid))
                        write(*,*) "max abs val of w_real:  ", maxval(abs(w_real))
                        write(*,*) "max abs val of alpha:   ", maxval(abs(alpha_d))
                        write(*,*) "max abs val of div:     ", maxval(abs(div))
                        call domain_check_winds(domain, "Solver Diverged: ", dqdt=.True.)
                    end associate
                endif
                stop
            endif
        endif

        call calc_updated_winds(domain, adv_den)
    end subroutine calc_iter_winds


    !>------------------------------------------------------------
    !! Right-preconditioned BiCGStab with merged reductions.
    !!
    !! Solves A x = b with Block-Jacobi (diagonal) preconditioner M.
    !!
    !! Per iteration: 2 SpMVs, 2 halo exchanges, 2 MPI_Iallreduce
    !! calls, 2 preconditioner applies (fresh M^{-1} on p and s — no
    !! auxiliary recurrences, so no roundoff accumulation that could
    !! stall convergence at high iter count).
    !!
    !! Convergence test mirrors AMGX's COMBINED_REL_INI_ABS:
    !!   converged when ||r|| <= max(tol_abs, tol_rel * ||r0||).
    !! ||r_new||^2 is computed locally from the packed reduction:
    !!   ||r_new||^2 = <s,s> - 2*omega*<t,s> + omega^2*<t,t>
    !! so no separate convergence allreduce is needed.
    !!
    !! 2 allreduces/iter:
    !!   1) sigma = <r̂, v>                                  (1 double)
    !!   2) <t,s>, <t,t>, <r̂,s>, <r̂,t>, <s,s>               (5 doubles)
    !!------------------------------------------------------------
    subroutine bicgstab_solve(domain, max_iters, status_out, n_iters_out, res0_out, res_final_out)
        implicit none
        type(domain_t), intent(in)  :: domain
        integer,        intent(in)  :: max_iters
        integer,        intent(out) :: status_out      ! 0 = converged, 1 = iter cap, 2 = breakdown
        integer,        intent(out) :: n_iters_out
        real(c_double), intent(out) :: res0_out, res_final_out

        integer :: it, ierr
        real(c_double) :: rho, rho_new, alpha_s, beta_s, omega
        real(c_double) :: sigma_local, sigma_global
        real(c_double) :: red5_local(5), red5_global(5)
        real(c_double) :: ts, tt, rs, rt, ss
        real(c_double) :: rnorm_global, rnorm_squared, target_norm
        real(c_double) :: rho0_pack(2)
        type(MPI_Request) :: req_sigma, req_red5
        real(c_double) :: t0_solve, t0_region

        status_out    = 1
        n_iters_out   = 0
        res0_out      = 0.0_c_double
        res_final_out = 0.0_c_double

        ! Reset timing accumulators for this solve
        t_spmv_acc      = 0.0_c_double
        t_halo_acc      = 0.0_c_double
        t_precond_acc   = 0.0_c_double
        t_allreduce_acc = 0.0_c_double
        t_vecops_acc    = 0.0_c_double
        t0_solve        = MPI_Wtime()

        ! --- Initial setup ---
        t0_region = MPI_Wtime()
        ! r0 = b - A*x0.  x0 is zero (set by caller), so r0 = b.
        call vec_copy(r_vec, rhs)
        ! r̂0 = r0 (shadow residual; held fixed across the solve)
        call vec_copy(r_hat, r_vec)
        ! p_0 = r_0
        call vec_copy(p_vec, r_vec)

        ! Initial reductions: rho_0 = <r̂, r_0>, ||r_0||²
        call vec_dot_local(r_hat, r_vec, rho0_pack(1))
        call vec_norm2_local(r_vec, rho0_pack(2))
        t_vecops_acc = t_vecops_acc + (MPI_Wtime() - t0_region)

        t0_region = MPI_Wtime()
#ifdef USE_NCCL
        rho0_dev(1) = rho0_pack(1); rho0_dev(2) = rho0_pack(2)
        !$acc update device(rho0_dev)
        !$acc host_data use_device(rho0_dev)
        ierr = nccl_allreduce_double_sum(c_loc(rho0_dev), c_loc(rho0_dev), 2, nccl_comm, nccl_stream)
        !$acc end host_data
        !$acc update host(rho0_dev)
        rho0_pack(1) = rho0_dev(1); rho0_pack(2) = rho0_dev(2)
#else
        call MPI_Allreduce(MPI_IN_PLACE, rho0_pack, 2, MPI_DOUBLE_PRECISION, &
                           MPI_SUM, solver_comm, ierr)
#endif
        t_allreduce_acc = t_allreduce_acc + (MPI_Wtime() - t0_region)

        rho           = rho0_pack(1)
        rnorm_global  = sqrt(rho0_pack(2))
        res0_out      = rnorm_global
        res_final_out = rnorm_global
        target_norm   = max(bicg_tol_abs, bicg_tol_rel * rnorm_global)

        if (rnorm_global <= bicg_tol_abs) then
            status_out  = 0
            n_iters_out = 0
            t_total_acc = MPI_Wtime() - t0_solve
            return
        endif

        omega   = 1.0_c_double
        alpha_s = 1.0_c_double

        do it = 1, max_iters

            ! ----- First half-step: p_hat = M^{-1} p, v = A p_hat -----
            t0_region = MPI_Wtime()
            call apply_precond(p_vec, p_hat)
            t_precond_acc = t_precond_acc + (MPI_Wtime() - t0_region)

            t0_region = MPI_Wtime()
            call exchange_krylov_halos(p_hat, domain)
            t_halo_acc = t_halo_acc + (MPI_Wtime() - t0_region)

            t0_region = MPI_Wtime()
            call spmv(p_hat, v_vec)
            t_spmv_acc = t_spmv_acc + (MPI_Wtime() - t0_region)

            ! ============== ALLREDUCE 1: sigma = <r̂, v> ==================
            t0_region = MPI_Wtime()
            call vec_dot_local(r_hat, v_vec, sigma_local)
            t_vecops_acc = t_vecops_acc + (MPI_Wtime() - t0_region)

            t0_region = MPI_Wtime()
#ifdef USE_NCCL
            sigma_dev(1) = sigma_local
            !$acc update device(sigma_dev)
            !$acc host_data use_device(sigma_dev)
            ierr = nccl_allreduce_double_sum(c_loc(sigma_dev), c_loc(sigma_dev), 1, nccl_comm, nccl_stream)
            !$acc end host_data
            !$acc update host(sigma_dev)
            sigma_global = sigma_dev(1)
#else
            call MPI_Iallreduce(sigma_local, sigma_global, 1, MPI_DOUBLE_PRECISION, &
                                MPI_SUM, solver_comm, req_sigma, ierr)
            ! No useful overlap candidate — wait
            call MPI_Wait(req_sigma, MPI_STATUS_IGNORE, ierr)
#endif
            t_allreduce_acc = t_allreduce_acc + (MPI_Wtime() - t0_region)

            if (abs(sigma_global) < breakdown_eps) then
                if (STD_OUT_PE) write(*,*) ' BiCGStab breakdown: <r̂, v> ~ 0 at iter ', it
                status_out  = 2
                n_iters_out = it
                exit
            endif
            alpha_s = rho / sigma_global

            ! s = r - alpha v
            t0_region = MPI_Wtime()
            call vec_axpby_into(s_vec, 1.0_c_double, r_vec, -alpha_s, v_vec)
            t_vecops_acc = t_vecops_acc + (MPI_Wtime() - t0_region)

            ! ----- Second half-step: s_hat = M^{-1} s, t = A s_hat -----
            t0_region = MPI_Wtime()
            call apply_precond(s_vec, s_hat)
            t_precond_acc = t_precond_acc + (MPI_Wtime() - t0_region)

            t0_region = MPI_Wtime()
            call exchange_krylov_halos(s_hat, domain)
            t_halo_acc = t_halo_acc + (MPI_Wtime() - t0_region)

            t0_region = MPI_Wtime()
            call spmv(s_hat, t_vec)
            t_spmv_acc = t_spmv_acc + (MPI_Wtime() - t0_region)

            ! ============== ALLREDUCE 2: 5-value packed =================
            ! <t,s>, <t,t>, <r̂,s>, <r̂,t>, <s,s>
            ! From these we get omega, rho_new (via <r̂,r_new> = <r̂,s> - omega <r̂,t>),
            ! and ||r_new||² (via ss - 2 omega ts + omega² tt) — convergence rides for free.
            t0_region = MPI_Wtime()
            call vec_dots_fused(t_vec, s_vec, r_hat, &
                                red5_local(1), red5_local(2), red5_local(3), &
                                red5_local(4), red5_local(5))
            t_vecops_acc = t_vecops_acc + (MPI_Wtime() - t0_region)

            t0_region = MPI_Wtime()
#ifdef USE_NCCL
            red5_dev(1:5) = red5_local(1:5)
            !$acc update device(red5_dev)
            !$acc host_data use_device(red5_dev)
            ierr = nccl_allreduce_double_sum(c_loc(red5_dev), c_loc(red5_dev), 5, nccl_comm, nccl_stream)
            !$acc end host_data
            !$acc update host(red5_dev)
            red5_global(1:5) = red5_dev(1:5)
#else
            call MPI_Iallreduce(red5_local, red5_global, 5, MPI_DOUBLE_PRECISION, &
                                MPI_SUM, solver_comm, req_red5, ierr)
            ! No useful overlap candidate — wait
            call MPI_Wait(req_red5, MPI_STATUS_IGNORE, ierr)
#endif
            t_allreduce_acc = t_allreduce_acc + (MPI_Wtime() - t0_region)

            ts = red5_global(1)
            tt = red5_global(2)
            rs = red5_global(3)
            rt = red5_global(4)
            ss = red5_global(5)

            if (tt < breakdown_eps) then
                if (STD_OUT_PE) write(*,*) ' BiCGStab breakdown: <t,t> ~ 0 at iter ', it
                status_out  = 2
                n_iters_out = it
                exit
            endif
            omega   = ts / tt
            rho_new = rs - omega * rt

            ! ||r_new||² = <s - omega t, s - omega t> = ss - 2*omega*ts + omega²*tt
            rnorm_squared = ss - 2.0_c_double*omega*ts + omega*omega*tt
            if (rnorm_squared < 0.0_c_double) rnorm_squared = 0.0_c_double  ! roundoff guard
            rnorm_global  = sqrt(rnorm_squared)
            res_final_out = rnorm_global
            n_iters_out   = it

            ! x = x + alpha*p_hat + omega*s_hat   (single fused kernel)
            t0_region = MPI_Wtime()
            call vec_axpy2(x_sol, alpha_s, p_hat, omega, s_hat)
            t_vecops_acc = t_vecops_acc + (MPI_Wtime() - t0_region)

            if (rnorm_global <= target_norm) then
                status_out = 0
                exit
            endif

            if (abs(omega) < breakdown_eps) then
                if (STD_OUT_PE) write(*,*) ' BiCGStab breakdown: omega ~ 0 at iter ', it
                status_out  = 2
                exit
            endif

            ! r_new = s - omega*t
            t0_region = MPI_Wtime()
            call vec_axpby_into(r_vec, 1.0_c_double, s_vec, -omega, t_vec)
            t_vecops_acc = t_vecops_acc + (MPI_Wtime() - t0_region)

            if (abs(rho) < breakdown_eps) then
                if (STD_OUT_PE) write(*,*) ' BiCGStab breakdown: rho ~ 0 at iter ', it
                status_out  = 2
                exit
            endif
            beta_s = (rho_new / rho) * (alpha_s / omega)

            ! p_new = r + beta * (p - omega * v)
            t0_region = MPI_Wtime()
            call vec_p_update(p_vec, r_vec, beta_s, omega, v_vec)
            t_vecops_acc = t_vecops_acc + (MPI_Wtime() - t0_region)

            rho = rho_new
        enddo

        t_total_acc = MPI_Wtime() - t0_solve

    end subroutine bicgstab_solve


    !>------------------------------------------------------------
    !! 15-point stencil SpMV: y = A x.
    !!
    !! Split into 4 kernels by row category to remove branches from
    !! the dominant interior loop:
    !!   1) Lateral identity rows (i=0, i=mx-1, j=0, j=my-1)
    !!   2) Top BC k=mz-1 (interior i,j only)
    !!   3) Bottom BC k=0 (interior i,j only)
    !!   4) Interior 15-pt stencil — branchless, vectorizable
    !!
    !! Interior optimisations:
    !!   - D/E/F/G coefs are all 1/dx² (constants) — use a scalar,
    !!     skip the 4 array reads.
    !!   - A coef is derivable as -4/dx² - B - C — compute on the
    !!     fly, skip the A_coef array read.
    !! Net: 5 fewer single-precision reads per interior cell out of
    !! ~140 bytes/cell traffic (~14% reduction), plus the branchless
    !! inner loop auto-vectorises cleanly.
    !!------------------------------------------------------------
    subroutine spmv(x, y)
        implicit none
        real(c_double), dimension(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1), intent(in)    :: x
        real(c_double), dimension(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1), intent(inout) :: y
        integer :: i, j, k
        integer :: i_int_lo, i_int_hi, j_int_lo, j_int_hi
        real(c_double) :: denom, dzdx_s, dzdy_s, dz_kp1, two_dx
        real(c_double) :: inv_dx2, neg_4_inv_dx2
        real(c_double) :: b_val, c_val, a_val

        two_dx        = 2.0_c_double * real(dx, c_double)
        inv_dx2       = 1.0_c_double / (real(dx, c_double) * real(dx, c_double))
        neg_4_inv_dx2 = -4.0_c_double * inv_dx2

        ! Interior (i,j) ranges (avoiding lateral BC cells)
        i_int_lo = max(xs,         1)
        i_int_hi = min(xs + xm - 1, mx - 2)
        j_int_lo = max(ys,         1)
        j_int_hi = min(ys + ym - 1, my - 2)

        ! ============== 1. Lateral identity rows ==============
        ! West face (i=0): only owned at west-boundary ranks
        if (xs == 0) then
            !$acc parallel loop gang vector collapse(2) present(x, y)
            do j = ys, ys + ym - 1
                do k = zs, zs + zm - 1
                    y(0, k, j) = x(0, k, j)
                enddo
            enddo
        endif
        ! East face (i=mx-1)
        if (xs + xm - 1 == mx - 1) then
            !$acc parallel loop gang vector collapse(2) present(x, y)
            do j = ys, ys + ym - 1
                do k = zs, zs + zm - 1
                    y(mx-1, k, j) = x(mx-1, k, j)
                enddo
            enddo
        endif
        ! South face (j=0)
        if (ys == 0) then
            !$acc parallel loop gang vector collapse(2) present(x, y)
            do k = zs, zs + zm - 1
                do i = xs, xs + xm - 1
                    y(i, k, 0) = x(i, k, 0)
                enddo
            enddo
        endif
        ! North face (j=my-1)
        if (ys + ym - 1 == my - 1) then
            !$acc parallel loop gang vector collapse(2) present(x, y)
            do k = zs, zs + zm - 1
                do i = xs, xs + xm - 1
                    y(i, k, my-1) = x(i, k, my-1)
                enddo
            enddo
        endif

        ! ============== 2. Top BC k=mz-1 ==============
        ! 2-pt stencil: y = (x(i,k,j) - x(i,k-1,j)) / dz_if(i,k,j)
        ! Restrict to interior (i,j) so we don't double-write lateral cells.
        if (i_int_lo <= i_int_hi .and. j_int_lo <= j_int_hi) then
            !$acc parallel loop gang vector collapse(2) present(x, y, dz_if)
            do j = j_int_lo, j_int_hi
                do i = i_int_lo, i_int_hi
                    y(i, mz-1, j) = (x(i, mz-1, j) - x(i, mz-2, j)) / real(dz_if(i, mz-1, j), c_double)
                enddo
            enddo
        endif

        ! ============== 3. Bottom BC k=0 ==============
        ! 10-pt terrain-following stencil
        if (i_int_lo <= i_int_hi .and. j_int_lo <= j_int_hi) then
            !$acc parallel loop gang vector collapse(2) &
            !$acc present(x, y, dzdx_surf, dzdy_surf, alpha, jaco, dz_if) &
            !$acc private(denom, dzdx_s, dzdy_s, dz_kp1)
            do j = j_int_lo, j_int_hi
                do i = i_int_lo, i_int_hi
                    dzdx_s = real(dzdx_surf(i,j), c_double)
                    dzdy_s = real(dzdy_surf(i,j), c_double)
                    dz_kp1 = real(dz_if(i,1,j), c_double)
                    denom  = 2.0_c_double * (dzdx_s**2 + dzdy_s**2 + real(alpha(i,k_s,j),c_double)**2) &
                             / real(jaco(i,k_s,j), c_double)
                    y(i,0,j) =                                                                       &
                          ( dzdy_s/(denom*two_dx)) * x(i,   0, j-1)                                  &
                        + ( dzdy_s/(denom*two_dx)) * x(i,   1, j-1)                                  &
                        + ( dzdx_s/(denom*two_dx)) * x(i-1, 0, j  )                                  &
                        + (-1.0_c_double / dz_kp1) * x(i,   0, j  )                                  &
                        + (-dzdx_s/(denom*two_dx)) * x(i+1, 0, j  )                                  &
                        + ( dzdx_s/(denom*two_dx)) * x(i-1, 1, j  )                                  &
                        + ( 1.0_c_double / dz_kp1) * x(i,   1, j  )                                  &
                        + (-dzdx_s/(denom*two_dx)) * x(i+1, 1, j  )                                  &
                        + (-dzdy_s/(denom*two_dx)) * x(i,   0, j+1)                                  &
                        + (-dzdy_s/(denom*two_dx)) * x(i,   1, j+1)
                enddo
            enddo
        endif

        ! ============== 4. Interior 15-pt stencil (branchless, vectorisable) ==============
        ! Skip A/D/E/F/G array reads:
        !   D = E = F = G = 1/dx²  (use scalar)
        !   A = -4/dx² - B - C     (compute on the fly)
        if (i_int_lo <= i_int_hi .and. j_int_lo <= j_int_hi) then
            !$acc parallel loop gang vector collapse(3) &
            !$acc present(x, y, B_coef, C_coef, H_coef, I_coef, J_coef, K_coef, &
            !$acc         L_coef, M_coef, N_coef, O_coef) &
            !$acc private(b_val, c_val, a_val)
            do j = j_int_lo, j_int_hi
                do k = 1, mz - 2
                    do i = i_int_lo, i_int_hi
                        b_val = real(B_coef(i,k,j), c_double)
                        c_val = real(C_coef(i,k,j), c_double)
                        a_val = neg_4_inv_dx2 - b_val - c_val

                        y(i,k,j) =                                                       &
                              real(O_coef(i,k,j), c_double) * x(i,   k-1, j-1)           &
                            + inv_dx2                       * x(i,   k,   j-1)           &
                            + real(M_coef(i,k,j), c_double) * x(i,   k+1, j-1)           &
                            + real(K_coef(i,k,j), c_double) * x(i-1, k-1, j  )           &
                            + c_val                         * x(i,   k-1, j  )           &
                            + real(J_coef(i,k,j), c_double) * x(i+1, k-1, j  )           &
                            + inv_dx2                       * x(i-1, k,   j  )           &
                            + a_val                         * x(i,   k,   j  )           &
                            + inv_dx2                       * x(i+1, k,   j  )           &
                            + real(I_coef(i,k,j), c_double) * x(i-1, k+1, j  )           &
                            + b_val                         * x(i,   k+1, j  )           &
                            + real(H_coef(i,k,j), c_double) * x(i+1, k+1, j  )           &
                            + real(N_coef(i,k,j), c_double) * x(i,   k-1, j+1)           &
                            + inv_dx2                       * x(i,   k,   j+1)           &
                            + real(L_coef(i,k,j), c_double) * x(i,   k+1, j+1)
                    enddo
                enddo
            enddo
        endif

    end subroutine spmv


    !>------------------------------------------------------------
    !! Block-Jacobi diagonal extraction.
    !!
    !! D(row) per row category:
    !!   lateral identity row -> 1
    !!   top BC k=mz-1        -> 1/dz_if(i,k,j)
    !!   bottom BC k=0        -> -1/dz_if(i,k+1,j)
    !!   interior             -> A_coef(i,k,j)
    !!
    !! Stored as the inverse so the apply step is a single multiply.
    !!------------------------------------------------------------
    subroutine extract_diagonal()
        implicit none
        integer :: i, j, k
        real(c_double) :: d

        !$acc parallel loop gang vector collapse(3) &
        !$acc present(D_inv, A_coef, dz_if) private(d)
        do j = ys, ys + ym - 1
            do k = zs, zs + zm - 1
                do i = xs, xs + xm - 1
                    if (i <= 0 .or. j <= 0 .or. i >= mx-1 .or. j >= my-1) then
                        d = 1.0_c_double
                    else if (k >= mz-1) then
                        d = 1.0_c_double / real(dz_if(i,k,j), c_double)
                    else if (k <= 0) then
                        d = -1.0_c_double / real(dz_if(i,k+1,j), c_double)
                    else
                        d = real(A_coef(i,k,j), c_double)
                    endif
                    ! Guard against zero diagonal (shouldn't happen for this matrix)
                    if (abs(d) < breakdown_eps) d = sign(breakdown_eps, d)
                    D_inv(i,k,j) = 1.0_c_double / d
                enddo
            enddo
        enddo
    end subroutine extract_diagonal


    !>------------------------------------------------------------
    !! Block-Jacobi preconditioner apply with multi-sweep Richardson.
    !!
    !! Applies M^{-1} to in_vec, with M^{-1} approximated by `precond_n_sweeps`
    !! Richardson sweeps using the diagonal D as the smoother:
    !!   y_0 = 0
    !!   y_{k+1} = y_k + D^{-1} (in - A y_k)
    !!
    !! Sweep 1 reduces to y = D^{-1} in (single diagonal scale).
    !! Subsequent sweeps add 1 SpMV + 1 fused update each.
    !!
    !! Block-Jacobi convention: no halo exchange between inner sweeps.
    !! Halo cells of out_vec stay at zero (initialised at first solve, never
    !! written by any solver routine) → inner SpMV uses zero halo, equivalent
    !! to homogeneous Dirichlet at the rank boundary. This is the standard
    !! approximation used by AMGX BLOCK_JACOBI and Hypre's PCJacobi.
    !!------------------------------------------------------------
    subroutine apply_precond(in_vec, out_vec)
        implicit none
        real(c_double), dimension(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1), intent(in)    :: in_vec
        real(c_double), dimension(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1), intent(inout) :: out_vec
        integer :: i, j, k, sweep

        ! Zero entire out_vec INCLUDING halo cells before sweeps. The inner SpMV
        ! in sweep 2+ reads out_vec at halo positions; for block-Jacobi convention
        ! those must be zero (homogeneous Dirichlet at rank boundary). Without
        ! this, halo cells leak previous outer-iter's exchange_krylov_halos data
        ! into the smoother and convergence drifts on multi-rank runs (bug
        ! manifested as slow residual growth at ~0.25%/iter, 2-GPU 2026-05-10).
        ! Halo cells of out_vec stay at zero throughout the precond apply since
        ! sweeps only write to owned cells.
        !$acc parallel loop gang vector collapse(3) present(out_vec)
        do j = j_s-1, j_e+1
            do k = k_s-1, k_e+1
                do i = i_s-1, i_e+1
                    out_vec(i,k,j) = 0.0_c_double
                enddo
            enddo
        enddo

        ! Sweep 1: y_1 = D^{-1} in   (starting from y_0 = 0)
        !$acc parallel loop gang vector collapse(3) present(in_vec, out_vec, D_inv)
        do j = ys, ys + ym - 1
            do k = zs, zs + zm - 1
                do i = xs, xs + xm - 1
                    out_vec(i,k,j) = D_inv(i,k,j) * in_vec(i,k,j)
                enddo
            enddo
        enddo

        ! Sweeps 2..n: y += D^{-1} (in - A y).  Inner SpMV uses out_vec's
        ! zero halo (no MPI exchange — Jacobi convention).
        do sweep = 2, precond_n_sweeps
            call spmv(out_vec, prec_res)
            !$acc parallel loop gang vector collapse(3) &
            !$acc present(in_vec, out_vec, D_inv, prec_res)
            do j = ys, ys + ym - 1
                do k = zs, zs + zm - 1
                    do i = xs, xs + xm - 1
                        out_vec(i,k,j) = out_vec(i,k,j) + D_inv(i,k,j) * (in_vec(i,k,j) - prec_res(i,k,j))
                    enddo
                enddo
            enddo
        enddo
    end subroutine apply_precond


    !>------------------------------------------------------------
    !! Vector ops over owned cell range (xs..xs+xm-1 etc.).
    !! Local reductions only — global allreduce is the caller's job.
    !!------------------------------------------------------------

    subroutine vec_zero(v)
        real(c_double), dimension(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1), intent(inout) :: v
        integer :: i, j, k
        !$acc parallel loop gang vector collapse(3) present(v)
        do j = ys, ys + ym - 1
            do k = zs, zs + zm - 1
                do i = xs, xs + xm - 1
                    v(i,k,j) = 0.0_c_double
                enddo
            enddo
        enddo
    end subroutine vec_zero

    subroutine vec_copy(dst, src)
        real(c_double), dimension(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1), intent(in)    :: src
        real(c_double), dimension(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1), intent(inout) :: dst
        integer :: i, j, k
        !$acc parallel loop gang vector collapse(3) present(dst, src)
        do j = ys, ys + ym - 1
            do k = zs, zs + zm - 1
                do i = xs, xs + xm - 1
                    dst(i,k,j) = src(i,k,j)
                enddo
            enddo
        enddo
    end subroutine vec_copy

    !> y = a*x + y
    subroutine vec_axpy(y, a, x)
        real(c_double), dimension(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1), intent(inout) :: y
        real(c_double), dimension(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1), intent(in)    :: x
        real(c_double), intent(in) :: a
        integer :: i, j, k
        !$acc parallel loop gang vector collapse(3) present(y, x)
        do j = ys, ys + ym - 1
            do k = zs, zs + zm - 1
                do i = xs, xs + xm - 1
                    y(i,k,j) = a * x(i,k,j) + y(i,k,j)
                enddo
            enddo
        enddo
    end subroutine vec_axpy

    !> y = y + a*x + b*z  (fused 3-term update — used for the x_sol += alpha*p_hat + omega*s_hat
    !> step in BiCGStab, replaces two sequential vec_axpy calls and saves one full pass through y).
    subroutine vec_axpy2(y, a, x, b, z)
        real(c_double), dimension(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1), intent(inout) :: y
        real(c_double), dimension(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1), intent(in)    :: x, z
        real(c_double), intent(in) :: a, b
        integer :: i, j, k
        !$acc parallel loop gang vector collapse(3) present(y, x, z)
        do j = ys, ys + ym - 1
            do k = zs, zs + zm - 1
                do i = xs, xs + xm - 1
                    y(i,k,j) = y(i,k,j) + a * x(i,k,j) + b * z(i,k,j)
                enddo
            enddo
        enddo
    end subroutine vec_axpy2

    !> p = r + beta * (p - omega * v)  (classic BiCGStab search-direction update; kept for reference)
    subroutine vec_p_update(p, r, beta_s, omega_s, v)
        real(c_double), dimension(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1), intent(inout) :: p
        real(c_double), dimension(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1), intent(in)    :: r, v
        real(c_double), intent(in) :: beta_s, omega_s
        integer :: i, j, k
        !$acc parallel loop gang vector collapse(3) present(p, r, v)
        do j = ys, ys + ym - 1
            do k = zs, zs + zm - 1
                do i = xs, xs + xm - 1
                    p(i,k,j) = r(i,k,j) + beta_s * (p(i,k,j) - omega_s * v(i,k,j))
                enddo
            enddo
        enddo
    end subroutine vec_p_update

    !> p_hat = u + beta * (p_hat - omega * q)  (pipelined-BiCGStab fused recurrence — replaces
    !> two AXPYs and avoids fresh M^{-1} apply on p)
    subroutine vec_p_hat_update(p_hat_v, u, beta_s, omega_s, q)
        real(c_double), dimension(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1), intent(inout) :: p_hat_v
        real(c_double), dimension(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1), intent(in)    :: u, q
        real(c_double), intent(in) :: beta_s, omega_s
        integer :: i, j, k
        !$acc parallel loop gang vector collapse(3) present(p_hat_v, u, q)
        do j = ys, ys + ym - 1
            do k = zs, zs + zm - 1
                do i = xs, xs + xm - 1
                    p_hat_v(i,k,j) = u(i,k,j) + beta_s * (p_hat_v(i,k,j) - omega_s * q(i,k,j))
                enddo
            enddo
        enddo
    end subroutine vec_p_hat_update

    !> y = a*x + b*y  (axpby — needed for the s_hat = u - alpha*q recurrence)
    subroutine vec_axpby_into(out_v, a, x, b, y)
        real(c_double), dimension(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1), intent(inout) :: out_v
        real(c_double), dimension(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1), intent(in)    :: x, y
        real(c_double), intent(in) :: a, b
        integer :: i, j, k
        !$acc parallel loop gang vector collapse(3) present(out_v, x, y)
        do j = ys, ys + ym - 1
            do k = zs, zs + zm - 1
                do i = xs, xs + xm - 1
                    out_v(i,k,j) = a * x(i,k,j) + b * y(i,k,j)
                enddo
            enddo
        enddo
    end subroutine vec_axpby_into

    !> result = sum_i x_i * y_i  (LOCAL — caller MPI_Allreduces)
    subroutine vec_dot_local(x, y, result_local)
        real(c_double), dimension(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1), intent(in) :: x, y
        real(c_double), intent(out) :: result_local
        integer :: i, j, k
        real(c_double) :: s
        s = 0.0_c_double
        !$acc parallel loop gang vector collapse(3) reduction(+:s) present(x, y)
        do j = ys, ys + ym - 1
            do k = zs, zs + zm - 1
                do i = xs, xs + xm - 1
                    s = s + x(i,k,j) * y(i,k,j)
                enddo
            enddo
        enddo
        result_local = s
    end subroutine vec_dot_local

    !> Single-pass fused inner products used in BiCGStab's omega allreduce:
    !>   ts = <t, s>, tt = <t, t>, rs = <r_hat, s>, rt = <r_hat, t>, ss = <s, s>
    !> Reads each of t, s, r_hat once instead of 5× sequentially — saves both kernel
    !> launch overhead and ~3× the bandwidth of separate dot products.
    subroutine vec_dots_fused(t_v, s_v, r_h, ts_l, tt_l, rs_l, rt_l, ss_l)
        real(c_double), dimension(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1), intent(in) :: t_v, s_v, r_h
        real(c_double), intent(out) :: ts_l, tt_l, rs_l, rt_l, ss_l
        integer :: i, j, k
        real(c_double) :: ts_acc, tt_acc, rs_acc, rt_acc, ss_acc
        real(c_double) :: tv, sv, rv

        ts_acc = 0.0_c_double; tt_acc = 0.0_c_double
        rs_acc = 0.0_c_double; rt_acc = 0.0_c_double
        ss_acc = 0.0_c_double

        !$acc parallel loop gang vector collapse(3) &
        !$acc reduction(+:ts_acc,tt_acc,rs_acc,rt_acc,ss_acc) &
        !$acc private(tv, sv, rv) present(t_v, s_v, r_h)
        do j = ys, ys + ym - 1
            do k = zs, zs + zm - 1
                do i = xs, xs + xm - 1
                    tv = t_v(i,k,j)
                    sv = s_v(i,k,j)
                    rv = r_h(i,k,j)
                    ts_acc = ts_acc + tv * sv
                    tt_acc = tt_acc + tv * tv
                    rs_acc = rs_acc + rv * sv
                    rt_acc = rt_acc + rv * tv
                    ss_acc = ss_acc + sv * sv
                enddo
            enddo
        enddo

        ts_l = ts_acc
        tt_l = tt_acc
        rs_l = rs_acc
        rt_l = rt_acc
        ss_l = ss_acc
    end subroutine vec_dots_fused

    !> result = sum_i x_i^2  (LOCAL — caller MPI_Allreduces, then sqrt)
    subroutine vec_norm2_local(x, result_local)
        real(c_double), dimension(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1), intent(in) :: x
        real(c_double), intent(out) :: result_local
        integer :: i, j, k
        real(c_double) :: s
        s = 0.0_c_double
        !$acc parallel loop gang vector collapse(3) reduction(+:s) present(x)
        do j = ys, ys + ym - 1
            do k = zs, zs + zm - 1
                do i = xs, xs + xm - 1
                    s = s + x(i,k,j) * x(i,k,j)
                enddo
            enddo
        enddo
        result_local = s
    end subroutine vec_norm2_local


    !>------------------------------------------------------------
    !! RHS vector: rhs(i,k,j) = -2*div(i,k,j) for interior; 0 at BCs.
    !! Halo cells are left at their previous value (won't be read).
    !!------------------------------------------------------------
    subroutine compute_rhs_3d()
        implicit none
        integer :: i, j, k

        !$acc parallel loop gang vector collapse(3) present(rhs, div)
        do j = ys, ys + ym - 1
            do k = zs, zs + zm - 1
                do i = xs, xs + xm - 1
                    if (i <= 0 .or. j <= 0 .or. i >= mx-1 .or. j >= my-1 .or. &
                        k <= 0 .or. k >= mz-1) then
                        rhs(i,k,j) = 0.0_c_double
                    else
                        rhs(i,k,j) = real(-2.0 * div(i,k,j), c_double)
                    endif
                enddo
            enddo
        enddo
    end subroutine compute_rhs_3d


    !>------------------------------------------------------------
    !! Halo exchange for a Krylov vector.
    !!
    !! Pack all 4 faces on GPU into persistent buffers, post Irecv
    !! on all non-boundary directions, post Isend on the same, wait
    !! for completion of all 8 (or fewer) requests, unpack into
    !! halo cells on GPU.
    !!
    !! The Isend/Irecv/Waitall pattern lets all 4 directions of comm
    !! happen in parallel — total halo time ~ max(4 directions) instead
    !! of sum. Combined with face-only PCIe and persistent buffers,
    !! halo cost is now latency-bound rather than allocator-bound.
    !!
    !! Tag convention:
    !!   tag_ew (200): data going east (R sends east face, R+1 recvs as west halo)
    !!   tag_we (201): data going west (R sends west face, R-1 recvs as east halo)
    !!   tag_ns (202): data going north (j+ direction)
    !!   tag_sn (203): data going south (j- direction)
    !!------------------------------------------------------------
    subroutine exchange_krylov_halos(v, domain)
        implicit none
        real(c_double), dimension(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1), intent(inout) :: v
        type(domain_t), intent(in) :: domain
        integer :: ierr, i, j, k, nreq
        type(MPI_Request) :: reqs(8)
        integer :: nz_w, nx_w, ny_w
        integer, parameter :: tag_ew = 200, tag_we = 201, tag_ns = 202, tag_sn = 203

        nz_w = (k_e+1) - (k_s-1) + 1
        nx_w = (i_e+1) - (i_s-1) + 1
        ny_w = (j_e+1) - (j_s-1) + 1

        ! ===== Pack all four faces on GPU =====
        if (.not. domain%east_boundary) then
            !$acc parallel loop collapse(2) present(v, east_send)
            do j = j_s-1, j_e+1
                do k = k_s-1, k_e+1
                    east_send(k, j) = v(i_e, k, j)
                enddo
            enddo
        endif
        if (.not. domain%west_boundary) then
            !$acc parallel loop collapse(2) present(v, west_send)
            do j = j_s-1, j_e+1
                do k = k_s-1, k_e+1
                    west_send(k, j) = v(i_s, k, j)
                enddo
            enddo
        endif
        if (.not. domain%north_boundary) then
            !$acc parallel loop collapse(2) present(v, north_send)
            do k = k_s-1, k_e+1
                do i = i_s-1, i_e+1
                    north_send(i, k) = v(i, k, j_e)
                enddo
            enddo
        endif
        if (.not. domain%south_boundary) then
            !$acc parallel loop collapse(2) present(v, south_send)
            do k = k_s-1, k_e+1
                do i = i_s-1, i_e+1
                    south_send(i, k) = v(i, k, j_s)
                enddo
            enddo
        endif

#ifdef USE_NCCL
        ! ===== NCCL path: device-to-device send/recv (no PCIe) =====
        ! Pack already happened on GPU into east_send/west_send/etc. above.
        ! Wrap all 4 directions of send+recv in nccl_group_start/end so NCCL
        ! fuses them into a single stream-ordered op.
        ! host_data use_device gives us device pointers for the buffers.
        !$acc host_data use_device(east_send, east_recv, west_send, west_recv, &
        !$acc                      north_send, north_recv, south_send, south_recv)
        call nccl_group_start()
        if (.not. domain%east_boundary) then
            ierr = nccl_recv_double(c_loc(east_recv),  nz_w*ny_w, east_neighbor,  nccl_comm, nccl_stream)
            ierr = nccl_send_double(c_loc(east_send),  nz_w*ny_w, east_neighbor,  nccl_comm, nccl_stream)
        endif
        if (.not. domain%west_boundary) then
            ierr = nccl_recv_double(c_loc(west_recv),  nz_w*ny_w, west_neighbor,  nccl_comm, nccl_stream)
            ierr = nccl_send_double(c_loc(west_send),  nz_w*ny_w, west_neighbor,  nccl_comm, nccl_stream)
        endif
        if (.not. domain%north_boundary) then
            ierr = nccl_recv_double(c_loc(north_recv), nx_w*nz_w, north_neighbor, nccl_comm, nccl_stream)
            ierr = nccl_send_double(c_loc(north_send), nx_w*nz_w, north_neighbor, nccl_comm, nccl_stream)
        endif
        if (.not. domain%south_boundary) then
            ierr = nccl_recv_double(c_loc(south_recv), nx_w*nz_w, south_neighbor, nccl_comm, nccl_stream)
            ierr = nccl_send_double(c_loc(south_send), nx_w*nz_w, south_neighbor, nccl_comm, nccl_stream)
        endif
        call nccl_group_end()
        !$acc end host_data
        ! NCCL ops are stream-ordered — subsequent OpenACC kernels on the same
        ! sync queue will see the recv buffers populated.

#else
        ! ===== MPI path: face data via host (face-only PCIe) =====
        ! H2D for send buffers (face-only — small)
        if (.not. domain%east_boundary)  then; !$acc update host(east_send)
        endif
        if (.not. domain%west_boundary)  then; !$acc update host(west_send)
        endif
        if (.not. domain%north_boundary) then; !$acc update host(north_send)
        endif
        if (.not. domain%south_boundary) then; !$acc update host(south_send)
        endif

        ! Post all Irecvs first (to maximise overlap with Isend)
        nreq = 0
        if (.not. domain%east_boundary) then
            nreq = nreq + 1
            call MPI_Irecv(east_recv,  nz_w*ny_w, MPI_DOUBLE_PRECISION, east_neighbor,  tag_we, &
                           solver_comm, reqs(nreq), ierr)
        endif
        if (.not. domain%west_boundary) then
            nreq = nreq + 1
            call MPI_Irecv(west_recv,  nz_w*ny_w, MPI_DOUBLE_PRECISION, west_neighbor,  tag_ew, &
                           solver_comm, reqs(nreq), ierr)
        endif
        if (.not. domain%north_boundary) then
            nreq = nreq + 1
            call MPI_Irecv(north_recv, nx_w*nz_w, MPI_DOUBLE_PRECISION, north_neighbor, tag_sn, &
                           solver_comm, reqs(nreq), ierr)
        endif
        if (.not. domain%south_boundary) then
            nreq = nreq + 1
            call MPI_Irecv(south_recv, nx_w*nz_w, MPI_DOUBLE_PRECISION, south_neighbor, tag_ns, &
                           solver_comm, reqs(nreq), ierr)
        endif

        ! Post all Isends
        if (.not. domain%east_boundary) then
            nreq = nreq + 1
            call MPI_Isend(east_send,  nz_w*ny_w, MPI_DOUBLE_PRECISION, east_neighbor,  tag_ew, &
                           solver_comm, reqs(nreq), ierr)
        endif
        if (.not. domain%west_boundary) then
            nreq = nreq + 1
            call MPI_Isend(west_send,  nz_w*ny_w, MPI_DOUBLE_PRECISION, west_neighbor,  tag_we, &
                           solver_comm, reqs(nreq), ierr)
        endif
        if (.not. domain%north_boundary) then
            nreq = nreq + 1
            call MPI_Isend(north_send, nx_w*nz_w, MPI_DOUBLE_PRECISION, north_neighbor, tag_ns, &
                           solver_comm, reqs(nreq), ierr)
        endif
        if (.not. domain%south_boundary) then
            nreq = nreq + 1
            call MPI_Isend(south_send, nx_w*nz_w, MPI_DOUBLE_PRECISION, south_neighbor, tag_sn, &
                           solver_comm, reqs(nreq), ierr)
        endif

        ! Wait for all
        if (nreq > 0) call MPI_Waitall(nreq, reqs(1:nreq), MPI_STATUSES_IGNORE, ierr)

        ! D2H for recv buffers (face-only)
        if (.not. domain%east_boundary)  then; !$acc update device(east_recv)
        endif
        if (.not. domain%west_boundary)  then; !$acc update device(west_recv)
        endif
        if (.not. domain%north_boundary) then; !$acc update device(north_recv)
        endif
        if (.not. domain%south_boundary) then; !$acc update device(south_recv)
        endif
#endif

        ! ===== Unpack on GPU =====
        if (.not. domain%east_boundary) then
            !$acc parallel loop collapse(2) present(v, east_recv)
            do j = j_s-1, j_e+1
                do k = k_s-1, k_e+1
                    v(i_e+1, k, j) = east_recv(k, j)
                enddo
            enddo
        endif
        if (.not. domain%west_boundary) then
            !$acc parallel loop collapse(2) present(v, west_recv)
            do j = j_s-1, j_e+1
                do k = k_s-1, k_e+1
                    v(i_s-1, k, j) = west_recv(k, j)
                enddo
            enddo
        endif
        if (.not. domain%north_boundary) then
            !$acc parallel loop collapse(2) present(v, north_recv)
            do k = k_s-1, k_e+1
                do i = i_s-1, i_e+1
                    v(i, k, j_e+1) = north_recv(i, k)
                enddo
            enddo
        endif
        if (.not. domain%south_boundary) then
            !$acc parallel loop collapse(2) present(v, south_recv)
            do k = k_s-1, k_e+1
                do i = i_s-1, i_e+1
                    v(i, k, j_s-1) = south_recv(i, k)
                enddo
            enddo
        endif
    end subroutine exchange_krylov_halos


    !>------------------------------------------------------------
    !! Apply the BiCGStab solution `x_sol` to the wind field.
    !!
    !! Direct port of calc_updated_winds from wind_iterative_amgx.F90:1117-1423,
    !! but consumes x_sol in 3D form (no flat→3D reshape needed).
    !!------------------------------------------------------------
    subroutine calc_updated_winds(domain, adv_den)
        implicit none
        type(domain_t), intent(inout) :: domain
        logical,        intent(in)    :: adv_den

        real, allocatable, dimension(:,:,:) :: u_dlambdz, v_dlambdz, dlambdz, u_temp, v_temp
        real, allocatable, dimension(:,:,:) :: lambda_3d
        real, allocatable, dimension(:,:,:) :: rho, rho_u, rho_v, rho_w
        integer :: i, j, k, i_start, i_end, j_start, j_end

        i_start = i_s
        i_end   = i_e + 1
        j_start = j_s
        j_end   = j_e + 1

        allocate(u_temp   (i_start:i_end,    k_s-1:k_e+1, j_s:j_e))
        allocate(v_temp   (i_s:i_e,          k_s-1:k_e+1, j_start:j_end))
        allocate(lambda_3d(i_s-1:i_e+1,      k_s-1:k_e+1, j_s-1:j_e+1))

        allocate(u_dlambdz(i_start:i_end,    k_s:k_e, j_s:j_e))
        allocate(v_dlambdz(i_s:i_e,          k_s:k_e, j_start:j_end))
        allocate(dlambdz  (i_s-1:i_e+1,      k_s:k_e, j_s-1:j_e+1))

        allocate(rho   (domain%ims:domain%ime, k_s:k_e, domain%jms:domain%jme))
        allocate(rho_u (i_start:i_end, k_s:k_e, j_s:j_e))
        allocate(rho_v (i_s:i_e,       k_s:k_e, j_start:j_end))
        allocate(rho_w (i_s:i_e,       k_s:k_e, j_s:j_e))

        !$acc enter data create(lambda_3d)

        ! Cast x_sol (double) into the single-precision lambda_3d on GPU
        !$acc parallel loop gang vector collapse(3) present(x_sol, lambda_3d)
        do j = j_s-1, j_e+1
            do k = k_s-1, k_e+1
                do i = i_s-1, i_e+1
                    lambda_3d(i,k,j) = real(x_sol(i,k,j))
                enddo
            enddo
        enddo

        ! Halo exchange for lambda_3d (single precision — reuse AMGX path's pattern)
        !$acc update host(lambda_3d)
        call exchange_lambda_halos_real(lambda_3d, domain)
        !$acc update device(lambda_3d)

        associate(density => domain%vars_3d(domain%var_indx(kVARS%density)%v)%data_3d, &
                  u       => domain%vars_3d(domain%var_indx(kVARS%u)%v)%dqdt_3d, &
                  v       => domain%vars_3d(domain%var_indx(kVARS%v)%v)%dqdt_3d, &
                  w       => domain%vars_3d(domain%var_indx(kVARS%w_real)%v)%data_3d, &
                  dz      => domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d, &
                  alpha_d => domain%vars_3d(domain%var_indx(kVARS%wind_alpha)%v)%data_3d, &
                  jaco_u_domain => domain%vars_3d(domain%var_indx(kVARS%jacobian_u)%v)%data_3d, &
                  jaco_v_domain => domain%vars_3d(domain%var_indx(kVARS%jacobian_v)%v)%data_3d, &
                  jaco_domain   => domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d, &
                  dzdx_u => domain%vars_3d(domain%var_indx(kVARS%dzdx_u)%v)%data_3d, &
                  dzdy_v => domain%vars_3d(domain%var_indx(kVARS%dzdy_v)%v)%data_3d)

        !$acc enter data create(u_temp, v_temp, u_dlambdz, v_dlambdz, dlambdz, rho, rho_u, rho_v, rho_w)
        !$acc data present(density, u, v, jaco_u_domain, jaco_v_domain, jaco_domain, dzdx_u, dzdy_v, &
        !$acc              u_temp, v_temp, u_dlambdz, v_dlambdz, dlambdz, rho, rho_u, rho_v, rho_w, &
        !$acc              alpha, dz_if, lambda_3d)

        !$acc kernels
        rho   = 1.0
        rho_w = 1.0
        rho_u = 1.0
        rho_v = 1.0
        !$acc end kernels

        if (adv_den) then
            !$acc parallel loop gang vector collapse(3)
            do j = jms, jme
                do k = k_s, k_e
                    do i = ims, ime
                        rho(i,k,j) = density(i,k,j)
                    enddo
                enddo
            enddo
        endif

        !$acc parallel
        if (i_s == ids .and. i_e == ide) then
            !$acc loop gang vector collapse(3)
            do j = j_s, j_e
                do k = k_s, k_e
                    do i = i_start+1, i_end-1
                        rho_u(i,k,j) = 0.5 * (rho(i,k,j) + rho(i-1,k,j))
                    enddo
                enddo
            enddo
            !$acc loop gang vector collapse(2)
            do j = j_s, j_e
                do k = k_s, k_e
                    rho_u(i_end,k,j)   = rho(i_end-1,k,j)
                    rho_u(i_start,k,j) = rho(i_start,k,j)
                enddo
            enddo
        else if (i_s == ids) then
            !$acc loop gang vector collapse(3)
            do j = j_s, j_e
                do k = k_s, k_e
                    do i = i_start+1, i_end
                        rho_u(i,k,j) = 0.5 * (rho(i,k,j) + rho(i-1,k,j))
                    enddo
                enddo
            enddo
            !$acc loop gang vector collapse(2)
            do j = j_s, j_e
                do k = k_s, k_e
                    rho_u(i_start,k,j) = rho(i_start,k,j)
                enddo
            enddo
        else if (i_e == ide) then
            !$acc loop gang vector collapse(3)
            do j = j_s, j_e
                do k = k_s, k_e
                    do i = i_start, i_end-1
                        rho_u(i,k,j) = 0.5 * (rho(i,k,j) + rho(i-1,k,j))
                    enddo
                enddo
            enddo
            !$acc loop gang vector collapse(2)
            do j = j_s, j_e
                do k = k_s, k_e
                    rho_u(i_end,k,j) = rho(i_end-1,k,j)
                enddo
            enddo
        else
            !$acc loop gang vector collapse(3)
            do j = j_s, j_e
                do k = k_s, k_e
                    do i = i_start, i_end
                        rho_u(i,k,j) = 0.5 * (rho(i,k,j) + rho(i-1,k,j))
                    enddo
                enddo
            enddo
        endif
        !$acc end parallel

        !$acc parallel
        if (j_s == jds .and. j_e == jde) then
            !$acc loop gang vector collapse(3)
            do j = j_start+1, j_end-1
                do k = k_s, k_e
                    do i = i_s, i_e
                        rho_v(i,k,j) = 0.5 * (rho(i,k,j) + rho(i,k,j-1))
                    enddo
                enddo
            enddo
            !$acc loop gang vector collapse(2)
            do k = k_s, k_e
                do i = i_s, i_e
                    rho_v(i,k,j_start) = rho(i,k,j_start)
                    rho_v(i,k,j_end)   = rho(i,k,j_end-1)
                enddo
            enddo
        else if (j_s == jds) then
            !$acc loop gang vector collapse(3)
            do j = j_start+1, j_end
                do k = k_s, k_e
                    do i = i_s, i_e
                        rho_v(i,k,j) = 0.5 * (rho(i,k,j) + rho(i,k,j-1))
                    enddo
                enddo
            enddo
            !$acc loop gang vector collapse(2)
            do k = k_s, k_e
                do i = i_s, i_e
                    rho_v(i,k,j_start) = rho(i,k,j_start)
                enddo
            enddo
        else if (j_e == jde) then
            !$acc loop gang vector collapse(3)
            do j = j_start, j_end-1
                do k = k_s, k_e
                    do i = i_s, i_e
                        rho_v(i,k,j) = 0.5 * (rho(i,k,j) + rho(i,k,j-1))
                    enddo
                enddo
            enddo
            !$acc loop gang vector collapse(2)
            do k = k_s, k_e
                do i = i_s, i_e
                    rho_v(i,k,j_end) = rho(i,k,j_end-1)
                enddo
            enddo
        else
            !$acc loop gang vector collapse(3)
            do j = j_start, j_end
                do k = k_s, k_e
                    do i = i_s, i_e
                        rho_v(i,k,j) = 0.5 * (rho(i,k,j) + rho(i,k,j-1))
                    enddo
                enddo
            enddo
        endif
        !$acc end parallel

        !$acc parallel loop gang vector tile(32,2,1)
        do j = j_s, j_e
            do k = k_s, k_e-1
                do i = i_s, i_e
                    rho_w(i,k,j) = ( rho(i,k,j)*dz(i,k+1,j) + rho(i,k+1,j)*dz(i,k,j) ) / &
                                   (dz(i,k,j) + dz(i,k+1,j))
                enddo
            enddo
        enddo

        !$acc parallel loop gang vector collapse(2)
        do j = j_s, j_e
            do i = i_s, i_e
                rho_w(i,k_e,j) = rho(i,k_e,j)
            enddo
        enddo

        !$acc parallel loop gang vector collapse(3)
        do j = j_s, j_e
            do k = k_s-1, k_e+1
                do i = i_start, i_end
                    u_temp(i,k,j) = (lambda_3d(i,k,j) + lambda_3d(i-1,k,j)) / 2
                enddo
            enddo
        enddo

        !$acc parallel loop gang vector collapse(3)
        do j = j_start, j_end
            do k = k_s-1, k_e+1
                do i = i_s, i_e
                    v_temp(i,k,j) = (lambda_3d(i,k,j) + lambda_3d(i,k,j-1)) / 2
                enddo
            enddo
        enddo

        !$acc parallel loop gang vector collapse(3)
        do j = j_s, j_e
            do k = k_s, k_e
                do i = i_start, i_end
                    u_dlambdz(i,k,j) = u_temp(i,k+1,j) - u_temp(i,k-1,j)
                    u_dlambdz(i,k,j) = u_dlambdz(i,k,j) / (dz_if(i_s,k+1,j_s) + dz_if(i_s,k,j_s))
                enddo
            enddo
        enddo

        !$acc parallel loop gang vector collapse(3)
        do j = j_start, j_end
            do k = k_s, k_e
                do i = i_s, i_e
                    v_dlambdz(i,k,j) = v_temp(i,k+1,j) - v_temp(i,k-1,j)
                    v_dlambdz(i,k,j) = v_dlambdz(i,k,j) / (dz_if(i_s,k+1,j_s) + dz_if(i_s,k,j_s))
                enddo
            enddo
        enddo

        !$acc parallel loop gang vector collapse(3)
        do j = j_s-1, j_e+1
            do k = k_s, k_e
                do i = i_s-1, i_e+1
                    dlambdz(i,k,j) = lambda_3d(i,k+1,j) - lambda_3d(i,k-1,j)
                    dlambdz(i,k,j) = dlambdz(i,k,j) / (dz_if(i_s,k+1,j_s) + dz_if(i_s,k,j_s))
                enddo
            enddo
        enddo

        !$acc parallel loop gang vector collapse(3)
        do j = j_s, j_e
            do k = k_s, k_e
                do i = i_start, i_end
                    u(i,k,j) = u(i,k,j) + 0.5 * ( (lambda_3d(i,k,j) - lambda_3d(i-1,k,j)) / dx - &
                                                   dzdx_u(i,k,j) * (u_dlambdz(i,k,j)) / jaco_u_domain(i,k,j) ) &
                                                / (rho_u(i,k,j))
                enddo
            enddo
        enddo

        !$acc parallel loop gang vector collapse(3)
        do j = j_start, j_end
            do k = k_s, k_e
                do i = i_s, i_e
                    v(i,k,j) = v(i,k,j) + 0.5 * ( (lambda_3d(i,k,j) - lambda_3d(i,k,j-1)) / dx - &
                                                   dzdy_v(i,k,j) * (v_dlambdz(i,k,j)) / jaco_v_domain(i,k,j) ) &
                                                / (rho_v(i,k,j))
                enddo
            enddo
        enddo

        !$acc parallel loop gang vector collapse(3)
        do j = j_s, j_e
            do k = k_s, k_e
                do i = i_s, i_e
                    w(i,k,j) = w(i,k,j) + 0.5 * (alpha(i,k,j)**2) * dlambdz(i,k,j) / jaco_domain(i,k,j) / rho_w(i,k,j)
                enddo
            enddo
        enddo

        !$acc end data
        !$acc exit data delete(u_temp, v_temp, u_dlambdz, v_dlambdz, dlambdz, rho, rho_u, rho_v, rho_w, lambda_3d)
        end associate

        deallocate(u_temp, v_temp, lambda_3d, u_dlambdz, v_dlambdz, dlambdz)
        deallocate(rho, rho_u, rho_v, rho_w)
    end subroutine calc_updated_winds


    !>------------------------------------------------------------
    !! Single-precision halo exchange for lambda_3d (used by
    !! calc_updated_winds). Direct mirror of
    !! exchange_lambda_halos in wind_iterative_amgx.F90:2214-2310.
    !!------------------------------------------------------------
    subroutine exchange_lambda_halos_real(lambda_3d, domain)
        implicit none
        real, dimension(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1), intent(inout) :: lambda_3d
        type(domain_t), intent(in) :: domain
        integer :: ierr, my_rank
        type(MPI_Status) :: mpi_stat
        integer :: nz, nx, ny
        integer :: east_rank, west_rank, north_rank, south_rank
        integer, parameter :: tag_ew = 100, tag_we = 101, tag_ns = 102, tag_sn = 103
        real, allocatable :: send_buff(:,:), rec_buff(:,:)

        nz = (k_e+1) - (k_s-1) + 1
        nx = (i_e+1) - (i_s-1) + 1
        ny = (j_e+1) - (j_s-1) + 1

        call MPI_Comm_rank(domain%compute_comms, my_rank, ierr)
        east_rank  = my_rank + 1
        west_rank  = my_rank - 1
        north_rank = my_rank + domain%grid%ximages
        south_rank = my_rank - domain%grid%ximages

        if (.not. domain%east_boundary) then
            allocate(send_buff(nz, ny), rec_buff(nz, ny))
            send_buff = lambda_3d(i_e, k_s-1:k_e+1, j_s-1:j_e+1)
            call MPI_Sendrecv(send_buff, nz*ny, MPI_REAL, east_rank, tag_ew, &
                              rec_buff,  nz*ny, MPI_REAL, east_rank, tag_we, &
                              domain%compute_comms, mpi_stat, ierr)
            lambda_3d(i_e+1, k_s-1:k_e+1, j_s-1:j_e+1) = rec_buff
            deallocate(send_buff, rec_buff)
        endif
        if (.not. domain%west_boundary) then
            allocate(send_buff(nz, ny), rec_buff(nz, ny))
            send_buff = lambda_3d(i_s, k_s-1:k_e+1, j_s-1:j_e+1)
            call MPI_Sendrecv(send_buff, nz*ny, MPI_REAL, west_rank, tag_we, &
                              rec_buff,  nz*ny, MPI_REAL, west_rank, tag_ew, &
                              domain%compute_comms, mpi_stat, ierr)
            lambda_3d(i_s-1, k_s-1:k_e+1, j_s-1:j_e+1) = rec_buff
            deallocate(send_buff, rec_buff)
        endif
        if (.not. domain%north_boundary) then
            allocate(send_buff(nx, nz), rec_buff(nx, nz))
            send_buff = lambda_3d(i_s-1:i_e+1, k_s-1:k_e+1, j_e)
            call MPI_Sendrecv(send_buff, nz*nx, MPI_REAL, north_rank, tag_ns, &
                              rec_buff,  nz*nx, MPI_REAL, north_rank, tag_sn, &
                              domain%compute_comms, mpi_stat, ierr)
            lambda_3d(i_s-1:i_e+1, k_s-1:k_e+1, j_e+1) = rec_buff
            deallocate(send_buff, rec_buff)
        endif
        if (.not. domain%south_boundary) then
            allocate(send_buff(nx, nz), rec_buff(nx, nz))
            send_buff = lambda_3d(i_s-1:i_e+1, k_s-1:k_e+1, j_s)
            call MPI_Sendrecv(send_buff, nz*nx, MPI_REAL, south_rank, tag_sn, &
                              rec_buff,  nz*nx, MPI_REAL, south_rank, tag_ns, &
                              domain%compute_comms, mpi_stat, ierr)
            lambda_3d(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1) = rec_buff
            deallocate(send_buff, rec_buff)
        endif
    end subroutine exchange_lambda_halos_real


    !>------------------------------------------------------------
    !! Cleanup: device delete + host deallocate.
    !!------------------------------------------------------------
    subroutine finalize_iter_winds()
        implicit none
        integer :: n

        if (.not. initialized_iter_winds) return

        ! Finalize all cached nests by restoring each one and cleaning it up.
        ! Mirrors finalize_amgx pattern.
        do n = 1, MAX_NESTS
            if (n /= active_nest_indx .and. domain_cache(n)%valid) then
                call restore_from_cache(n)
                call finalize_active_state()
            endif
        enddo

        ! Finalize the currently-active state
        call finalize_active_state()

#ifdef USE_NCCL
        if (nccl_initialized) then
            call nccl_comm_destroy(nccl_comm)
            nccl_comm = c_null_ptr
            ! nccl_stream is aliased to OpenACC's sync queue — do NOT destroy
            nccl_stream = c_null_ptr
            nccl_initialized = .false.
        endif
#endif

        initialized_iter_winds  = .false.
        structure_uploaded = .false.
        active_nest_indx   = -1
    end subroutine finalize_iter_winds


    !>------------------------------------------------------------
    !! Clean up the currently-active state (device data + host alloc).
    !! Used by finalize_iter_winds both for the live state and for each
    !! cached nest after restoring it.
    !!------------------------------------------------------------
    subroutine finalize_active_state()
        implicit none

        if (structure_uploaded) then
            !$acc exit data delete(x_sol, r_vec, r_hat, p_vec, p_hat, v_vec, &
            !$acc                  s_vec, s_hat, t_vec, rhs, D_inv, prec_res, &
            !$acc                  east_send, east_recv, west_send, west_recv, &
            !$acc                  north_send, north_recv, south_send, south_recv)
            !$acc exit data delete(A_coef, B_coef, C_coef, D_coef, E_coef, F_coef, G_coef, &
            !$acc                  H_coef, I_coef, J_coef, K_coef, L_coef, M_coef, N_coef, O_coef)
        endif

        if (allocated(x_sol))       deallocate(x_sol)
        if (allocated(r_vec))       deallocate(r_vec)
        if (allocated(r_hat))       deallocate(r_hat)
        if (allocated(p_vec))       deallocate(p_vec)
        if (allocated(p_hat))       deallocate(p_hat)
        if (allocated(v_vec))       deallocate(v_vec)
        if (allocated(s_vec))       deallocate(s_vec)
        if (allocated(s_hat))       deallocate(s_hat)
        if (allocated(t_vec))       deallocate(t_vec)
        if (allocated(rhs))         deallocate(rhs)
        if (allocated(D_inv))       deallocate(D_inv)
        if (allocated(prec_res))    deallocate(prec_res)
#ifdef USE_NCCL
        if (allocated(sigma_dev) .or. allocated(red5_dev) .or. allocated(rho0_dev)) then
            !$acc exit data delete(sigma_dev, red5_dev, rho0_dev)
        endif
        if (allocated(sigma_dev)) deallocate(sigma_dev)
        if (allocated(red5_dev))  deallocate(red5_dev)
        if (allocated(rho0_dev))  deallocate(rho0_dev)

        if (nccl_initialized) then
            call nccl_comm_destroy(nccl_comm)
            nccl_comm = c_null_ptr
            ! nccl_stream is aliased to OpenACC's sync queue — do NOT destroy
            nccl_stream = c_null_ptr
            nccl_initialized = .false.
        endif
#endif
        if (allocated(east_send))   deallocate(east_send)
        if (allocated(east_recv))   deallocate(east_recv)
        if (allocated(west_send))   deallocate(west_send)
        if (allocated(west_recv))   deallocate(west_recv)
        if (allocated(north_send))  deallocate(north_send)
        if (allocated(north_recv))  deallocate(north_recv)
        if (allocated(south_send))  deallocate(south_send)
        if (allocated(south_recv))  deallocate(south_recv)

        if (allocated(A_coef)) deallocate(A_coef)
        if (allocated(B_coef)) deallocate(B_coef)
        if (allocated(C_coef)) deallocate(C_coef)
        if (allocated(D_coef)) deallocate(D_coef)
        if (allocated(E_coef)) deallocate(E_coef)
        if (allocated(F_coef)) deallocate(F_coef)
        if (allocated(G_coef)) deallocate(G_coef)
        if (allocated(H_coef)) deallocate(H_coef)
        if (allocated(I_coef)) deallocate(I_coef)
        if (allocated(J_coef)) deallocate(J_coef)
        if (allocated(K_coef)) deallocate(K_coef)
        if (allocated(L_coef)) deallocate(L_coef)
        if (allocated(M_coef)) deallocate(M_coef)
        if (allocated(N_coef)) deallocate(N_coef)
        if (allocated(O_coef)) deallocate(O_coef)

        if (allocated(dz_if)) then
            !$acc exit data delete(dz_if, div, alpha, dzdx, dzdy, jaco, dzdx_surf, dzdy_surf, sigma)
            !$acc exit data delete(jaco_w, dzdx_u_stag, jaco_u_stag, dzdy_v_stag, jaco_v_stag, adv_dz_col)
            deallocate(dz_if)
        endif
        if (allocated(jaco))         deallocate(jaco)
        if (allocated(dzdx))         deallocate(dzdx)
        if (allocated(dzdy))         deallocate(dzdy)
        if (allocated(dzdx_surf))    deallocate(dzdx_surf)
        if (allocated(dzdy_surf))    deallocate(dzdy_surf)
        if (allocated(sigma))        deallocate(sigma)
        if (allocated(alpha))        deallocate(alpha)
        if (allocated(div))          deallocate(div)
        if (allocated(jaco_w))       deallocate(jaco_w)
        if (allocated(dzdx_u_stag))  deallocate(dzdx_u_stag)
        if (allocated(jaco_u_stag))  deallocate(jaco_u_stag)
        if (allocated(dzdy_v_stag))  deallocate(dzdy_v_stag)
        if (allocated(jaco_v_stag))  deallocate(jaco_v_stag)
        if (allocated(adv_dz_col))   deallocate(adv_dz_col)

        structure_uploaded = .false.
    end subroutine finalize_active_state


    !>------------------------------------------------------------
    !! Physics setup — copied verbatim from wind_iterative_amgx.F90
    !! (same source-of-truth for matrix coefficients).
    !!------------------------------------------------------------
    subroutine init_module_vars(domain)
        implicit none
        type(domain_t), intent(in) :: domain
        integer :: i, j, k, i_s_bnd, i_e_bnd, j_s_bnd, j_e_bnd

        i_s = domain%its;  i_e = domain%ite
        k_s = domain%kts;  k_e = domain%kte
        j_s = domain%jts;  j_e = domain%jte
        ims = domain%ims;  ime = domain%ime
        jms = domain%jms;  jme = domain%jme
        ids = domain%grid%ids;  ide = domain%grid%ide
        jds = domain%grid%jds;  jde = domain%grid%jde

        if (ims == ids) i_s = ids
        if (ime == ide) i_e = ide
        if (jms == jds) j_s = jds
        if (jme == jde) j_e = jde

        xs = i_s - merge(1, 0, i_s == ids)
        ys = j_s - merge(1, 0, j_s == jds)
        zs = 0

        mx = ide + 2
        my = jde + 2
        mz = domain%kde + 2

        xm = (i_e + merge(1, 0, i_e == ide)) - xs + 1
        ym = (j_e + merge(1, 0, j_e == jde)) - ys + 1
        zm = mz - zs

        i_s_bnd = i_s; i_e_bnd = i_e; j_s_bnd = j_s; j_e_bnd = j_e
        if (i_s == ids) i_s_bnd = i_s + 1
        if (i_e == ide) i_e_bnd = i_e - 1
        if (j_s == jds) j_s_bnd = j_s + 1
        if (j_e == jde) j_e_bnd = j_e - 1

        hs = domain%grid%halo_size

        if (.not. allocated(dzdx)) then
            allocate(dzdx     (i_s:i_e, k_s:k_e,   j_s:j_e))
            allocate(dzdy     (i_s:i_e, k_s:k_e,   j_s:j_e))
            allocate(jaco     (i_s:i_e, k_s:k_e,   j_s:j_e))
            allocate(dzdx_surf(i_s:i_e,            j_s:j_e))
            allocate(dzdy_surf(i_s:i_e,            j_s:j_e))
            allocate(sigma    (i_s:i_e, k_s:k_e,   j_s:j_e))
            allocate(dz_if    (i_s:i_e, k_s:k_e+1, j_s:j_e))
            allocate(alpha    (i_s:i_e, k_s:k_e,   j_s:j_e))
            allocate(div      (i_s:i_e, k_s:k_e,   j_s:j_e))
            allocate(jaco_w     (i_s:i_e,   k_s:k_e, j_s:j_e))
            allocate(dzdx_u_stag(i_s:i_e+1, k_s:k_e, j_s:j_e))
            allocate(jaco_u_stag(i_s:i_e+1, k_s:k_e, j_s:j_e))
            allocate(dzdy_v_stag(i_s:i_e,   k_s:k_e, j_s:j_e+1))
            allocate(jaco_v_stag(i_s:i_e,   k_s:k_e, j_s:j_e+1))
            allocate(adv_dz_col (k_s:k_e))

            !$acc enter data create(dz_if, div, alpha, dzdx, dzdy, jaco, dzdx_surf, dzdy_surf, sigma)
            !$acc enter data create(jaco_w, dzdx_u_stag, jaco_u_stag, dzdy_v_stag, jaco_v_stag, adv_dz_col)

            dx = domain%dx

            associate(advection_dz_var    => domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d, &
                      neighbor_terrain_var=> domain%vars_2d(domain%var_indx(kVARS%neighbor_terrain)%v)%data_2d, &
                      adv_dz_dom         => domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d, &
                      dzdx_domain         => domain%vars_3d(domain%var_indx(kVARS%dzdx)%v)%data_3d, &
                      dzdy_domain         => domain%vars_3d(domain%var_indx(kVARS%dzdy)%v)%data_3d, &
                      dzdx_u_domain      => domain%vars_3d(domain%var_indx(kVARS%dzdx_u)%v)%data_3d, &
                      dzdy_v_domain      => domain%vars_3d(domain%var_indx(kVARS%dzdy_v)%v)%data_3d, &
                      jaco_u_domain      => domain%vars_3d(domain%var_indx(kVARS%jacobian_u)%v)%data_3d, &
                      jaco_v_domain      => domain%vars_3d(domain%var_indx(kVARS%jacobian_v)%v)%data_3d, &
                      jaco_w_domain      => domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d, &
                      jaco_domain         => domain%vars_3d(domain%var_indx(kVARS%jacobian)%v)%data_3d)

            !$acc data present(advection_dz_var, neighbor_terrain_var, adv_dz_dom, dzdx_domain, dzdy_domain, dzdx_u_domain, dzdy_v_domain, &
            !$acc                  jaco_u_domain, jaco_v_domain, jaco_w_domain, jaco_domain, dz_if, div, alpha, dzdx, dzdy, jaco, dzdx_surf, dzdy_surf, sigma, &
            !$acc                  jaco_w, dzdx_u_stag, jaco_u_stag, dzdy_v_stag, jaco_v_stag, adv_dz_col)

            !$acc parallel loop gang vector collapse(3)
            do j = j_s, j_e
                do k = k_s+1, k_e
                    do i = i_s, i_e
                        dz_if(i,k,j) = (advection_dz_var(i,k,j) + advection_dz_var(i,k-1,j)) / 2
                    enddo
                enddo
            enddo
            !$acc parallel loop gang vector collapse(2)
            do j = j_s, j_e
                do i = i_s, i_e
                    dz_if(i,k_s,j)   = advection_dz_var(i,k_s,j)
                    dz_if(i,k_e+1,j) = advection_dz_var(i,k_e,j)
                    dzdx_surf(i,j) = dzdx_domain(i,k_s,j)
                    dzdy_surf(i,j) = dzdy_domain(i,k_s,j)
                enddo
            enddo
            !$acc parallel loop gang vector collapse(3)
            do j = j_s, j_e
                do k = k_s, k_e
                    do i = i_s, i_e
                        dzdx(i,k,j)  = dzdx_domain(i,k,j)
                        dzdy(i,k,j)  = dzdy_domain(i,k,j)
                        jaco(i,k,j)  = jaco_domain(i,k,j)
                        jaco_w(i,k,j)  = jaco_w_domain(i,k,j)
                        sigma(i,k,j) = dz_if(i,k,j) / dz_if(i,k+1,j)
                    enddo
                enddo
            enddo
            !$acc parallel loop gang vector collapse(3)
            do j = j_s, j_e
                do k = k_s, k_e
                    do i = i_s, i_e+1
                        dzdx_u_stag(i,k,j)  = dzdx_u_domain(i,k,j)
                        jaco_u_stag(i,k,j)  = jaco_u_domain(i,k,j)
                    enddo
                enddo
            enddo
            !$acc parallel loop gang vector collapse(3)
            do j = j_s, j_e+1
                do k = k_s, k_e
                    do i = i_s, i_e
                        dzdy_v_stag(i,k,j)  = dzdy_v_domain(i,k,j)
                        jaco_v_stag(i,k,j)  = jaco_v_domain(i,k,j)
                    enddo
                enddo
            enddo

            !$acc parallel loop gang vector collapse(2)
            do j = j_s, j_e
                do i = i_s_bnd, i_e_bnd
                    dzdx_surf(i,j) = (neighbor_terrain_var(i+1,j) - neighbor_terrain_var(i-1,j)) / (2*dx)
                enddo
            enddo
            !$acc parallel loop gang vector collapse(2)
            do j = j_s_bnd, j_e_bnd
                do i = i_s, i_e
                    dzdy_surf(i,j) = (neighbor_terrain_var(i,j+1) - neighbor_terrain_var(i,j-1)) / (2*dx)
                enddo
            enddo
            !$acc parallel loop
            do k = k_s, k_e
                adv_dz_col(k) = adv_dz_dom(i_s, k, j_s)
            enddo
            !$acc end data
            end associate

        endif
        !$acc update host(dz_if, div, alpha, dzdx, dzdy, jaco, dzdx_surf, dzdy_surf, sigma)
        !$acc update host(jaco_w, dzdx_u_stag, jaco_u_stag, dzdy_v_stag, jaco_v_stag, adv_dz_col)

    end subroutine init_module_vars


    subroutine initialize_coefs(domain)
        implicit none
        type(domain_t), intent(in) :: domain
        integer :: k
        real :: mixed_denom

        allocate(A_coef(i_s:i_e, k_s:k_e, j_s:j_e))
        allocate(B_coef(i_s:i_e, k_s:k_e, j_s:j_e))
        allocate(C_coef(i_s:i_e, k_s:k_e, j_s:j_e))
        allocate(D_coef(i_s:i_e, k_s:k_e, j_s:j_e))
        allocate(E_coef(i_s:i_e, k_s:k_e, j_s:j_e))
        allocate(F_coef(i_s:i_e, k_s:k_e, j_s:j_e))
        allocate(G_coef(i_s:i_e, k_s:k_e, j_s:j_e))
        allocate(H_coef(i_s:i_e, k_s:k_e, j_s:j_e))
        allocate(I_coef(i_s:i_e, k_s:k_e, j_s:j_e))
        allocate(J_coef(i_s:i_e, k_s:k_e, j_s:j_e))
        allocate(K_coef(i_s:i_e, k_s:k_e, j_s:j_e))
        allocate(L_coef(i_s:i_e, k_s:k_e, j_s:j_e))
        allocate(M_coef(i_s:i_e, k_s:k_e, j_s:j_e))
        allocate(N_coef(i_s:i_e, k_s:k_e, j_s:j_e))
        allocate(O_coef(i_s:i_e, k_s:k_e, j_s:j_e))

        A_coef = 1; B_coef = 0; C_coef = 0; D_coef = 0; E_coef = 0
        F_coef = 0; G_coef = 0; H_coef = 0; I_coef = 0; J_coef = 0
        K_coef = 0; L_coef = 0; M_coef = 0; N_coef = 0; O_coef = 0

        D_coef = 1.0/(domain%dx**2)
        E_coef = 1.0/(domain%dx**2)
        F_coef = 1.0/(domain%dx**2)
        G_coef = 1.0/(domain%dx**2)

        do k = k_s, k_e
            mixed_denom = 2*domain%dx*(dz_if(i_s,k+1,j_s) + dz_if(i_s,k,j_s))
            H_coef(:,k,:) = -(dzdx(:,k,:)/jaco(:,k,:) + domain%vars_3d(domain%var_indx(kVARS%dzdx_u)%v)%data_3d(i_s+1:i_e+1,k,j_s:j_e)/domain%vars_3d(domain%var_indx(kVARS%jacobian_u)%v)%data_3d(i_s+1:i_e+1,k,j_s:j_e))/mixed_denom
            I_coef(:,k,:) =  (dzdx(:,k,:)/jaco(:,k,:) + domain%vars_3d(domain%var_indx(kVARS%dzdx_u)%v)%data_3d(i_s:i_e,k,j_s:j_e)/domain%vars_3d(domain%var_indx(kVARS%jacobian_u)%v)%data_3d(i_s:i_e,k,j_s:j_e))/mixed_denom
            L_coef(:,k,:) = -(dzdy(:,k,:)/jaco(:,k,:) + domain%vars_3d(domain%var_indx(kVARS%dzdy_v)%v)%data_3d(i_s:i_e,k,j_s+1:j_e+1)/domain%vars_3d(domain%var_indx(kVARS%jacobian_v)%v)%data_3d(i_s:i_e,k,j_s+1:j_e+1))/mixed_denom
            M_coef(:,k,:) =  (dzdy(:,k,:)/jaco(:,k,:) + domain%vars_3d(domain%var_indx(kVARS%dzdy_v)%v)%data_3d(i_s:i_e,k,j_s:j_e)/domain%vars_3d(domain%var_indx(kVARS%jacobian_v)%v)%data_3d(i_s:i_e,k,j_s:j_e))/mixed_denom
            J_coef(:,k,:) =  (dzdx(:,k,:)/jaco(:,k,:) + domain%vars_3d(domain%var_indx(kVARS%dzdx_u)%v)%data_3d(i_s+1:i_e+1,k,j_s:j_e)/domain%vars_3d(domain%var_indx(kVARS%jacobian_u)%v)%data_3d(i_s+1:i_e+1,k,j_s:j_e))/mixed_denom
            K_coef(:,k,:) = -(dzdx(:,k,:)/jaco(:,k,:) + domain%vars_3d(domain%var_indx(kVARS%dzdx_u)%v)%data_3d(i_s:i_e,k,j_s:j_e)/domain%vars_3d(domain%var_indx(kVARS%jacobian_u)%v)%data_3d(i_s:i_e,k,j_s:j_e))/mixed_denom
            N_coef(:,k,:) =  (dzdy(:,k,:)/jaco(:,k,:) + domain%vars_3d(domain%var_indx(kVARS%dzdy_v)%v)%data_3d(i_s:i_e,k,j_s+1:j_e+1)/domain%vars_3d(domain%var_indx(kVARS%jacobian_v)%v)%data_3d(i_s:i_e,k,j_s+1:j_e+1))/mixed_denom
            O_coef(:,k,:) = -(dzdy(:,k,:)/jaco(:,k,:) + domain%vars_3d(domain%var_indx(kVARS%dzdy_v)%v)%data_3d(i_s:i_e,k,j_s:j_e)/domain%vars_3d(domain%var_indx(kVARS%jacobian_v)%v)%data_3d(i_s:i_e,k,j_s:j_e))/mixed_denom
        enddo

        ! B / C / A — host computation (used on first call before update_coefs_gpu)
        call update_coefs_host(domain)
    end subroutine initialize_coefs


    subroutine update_coefs_host(domain)
        implicit none
        type(domain_t), intent(in) :: domain
        real, allocatable, dimension(:,:,:) :: mixed_denom, X_coef
        real, allocatable, dimension(:,:)   :: M_up, M_dwn
        integer :: k

        allocate(mixed_denom(i_s:i_e, k_s:k_e, j_s:j_e))
        allocate(X_coef     (i_s:i_e, k_s:k_e, j_s:j_e))
        allocate(M_up       (i_s:i_e,          j_s:j_e))
        allocate(M_dwn      (i_s:i_e,          j_s:j_e))

        mixed_denom = (dz_if(:,k_s+1:k_e+1,:) + dz_if(:,k_s:k_e,:)) * 2 * domain%dx

        do k = k_s, k_e
            if (k == k_s) then
                M_up = dzdx(:,k,:)*(dzdx(:,k,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k+1,j_s) + &
                        dzdx(:,k+1,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s))/ &
                        (domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s)+domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k+1,j_s))/domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d(i_s:i_e,k,j_s:j_e)
                M_dwn = dzdx(:,k,:)*dzdx_surf
                M_up = M_up + dzdy(:,k,:)*(dzdy(:,k,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k+1,j_s) + &
                        dzdy(:,k+1,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s))/ &
                        (domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s)+domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k+1,j_s))/domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d(i_s:i_e,k,j_s:j_e)
                M_dwn = M_dwn + dzdy(:,k,:)*dzdy_surf
                M_up = M_up + (alpha(:,k,:)**2*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k+1,j_s) + &
                        alpha(:,k+1,:)**2*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s))/ &
                        (domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s)+domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k+1,j_s))/domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d(i_s:i_e,k,j_s:j_e)
                M_dwn = M_dwn + alpha(:,k,:)**2
            else if (k == k_e) then
                M_up = dzdx(:,k,:)*(dzdx(:,k,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s) + &
                        dzdx(:,k,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s))/ &
                        (domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s)+domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s))
                M_dwn = dzdx(:,k,:)*(dzdx(:,k,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k-1,j_s) + &
                        dzdx(:,k-1,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s))/ &
                        (domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s)+domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k-1,j_s))/domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d(i_s:i_e,k-1,j_s:j_e)
                M_up = M_up + dzdy(:,k,:)*(dzdy(:,k,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s) + &
                        dzdy(:,k,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s))/ &
                        (domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s)+domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s))
                M_dwn = M_dwn + dzdy(:,k,:)*(dzdy(:,k,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k-1,j_s) + &
                        dzdy(:,k-1,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s))/ &
                        (domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s)+domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k-1,j_s))/domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d(i_s:i_e,k-1,j_s:j_e)
                M_up = M_up + (alpha(:,k,:)**2*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s) + &
                        alpha(:,k,:)**2*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s))/ &
                        (domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s)+domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s))
                M_dwn = M_dwn + (alpha(:,k,:)**2*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k-1,j_s) + &
                        alpha(:,k-1,:)**2*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s))/ &
                        (domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s)+domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k-1,j_s))/domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d(i_s:i_e,k-1,j_s:j_e)
            else
                M_up = dzdx(:,k,:)*(dzdx(:,k,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k+1,j_s) + &
                        dzdx(:,k+1,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s))/ &
                        (domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s)+domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k+1,j_s))/domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d(i_s:i_e,k,j_s:j_e)
                M_dwn = dzdx(:,k,:)*(dzdx(:,k,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k-1,j_s) + &
                        dzdx(:,k-1,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s))/ &
                        (domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s)+domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k-1,j_s))/domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d(i_s:i_e,k-1,j_s:j_e)
                M_up = M_up + dzdy(:,k,:)*(dzdy(:,k,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k+1,j_s) + &
                        dzdy(:,k+1,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s))/ &
                        (domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s)+domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k+1,j_s))/domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d(i_s:i_e,k,j_s:j_e)
                M_dwn = M_dwn + dzdy(:,k,:)*(dzdy(:,k,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k-1,j_s) + &
                        dzdy(:,k-1,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s))/ &
                        (domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s)+domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k-1,j_s))/domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d(i_s:i_e,k-1,j_s:j_e)
                M_up = M_up + (alpha(:,k,:)**2*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k+1,j_s) + &
                        alpha(:,k+1,:)**2*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s))/ &
                        (domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s)+domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k+1,j_s))/domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d(i_s:i_e,k,j_s:j_e)
                M_dwn = M_dwn + (alpha(:,k,:)**2*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k-1,j_s) + &
                        alpha(:,k-1,:)**2*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s))/ &
                        (domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s)+domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k-1,j_s))/domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d(i_s:i_e,k-1,j_s:j_e)
            endif
            B_coef(:,k,:) = M_up /(jaco(:,k,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s)*dz_if(i_s,k+1,j_s))
            C_coef(:,k,:) = M_dwn/(jaco(:,k,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s)*dz_if(i_s,k,  j_s))
        enddo

        A_coef = -4/(domain%dx**2) - B_coef - C_coef

        X_coef = -((domain%vars_3d(domain%var_indx(kVARS%dzdx_u)%v)%data_3d(i_s+1:i_e+1,:,j_s:j_e)/domain%vars_3d(domain%var_indx(kVARS%jacobian_u)%v)%data_3d(i_s+1:i_e+1,:,j_s:j_e) - domain%vars_3d(domain%var_indx(kVARS%dzdx_u)%v)%data_3d(i_s:i_e,:,j_s:j_e)/domain%vars_3d(domain%var_indx(kVARS%jacobian_u)%v)%data_3d(i_s:i_e,:,j_s:j_e)) + &
                   (domain%vars_3d(domain%var_indx(kVARS%dzdy_v)%v)%data_3d(i_s:i_e,:,j_s+1:j_e+1)/domain%vars_3d(domain%var_indx(kVARS%jacobian_v)%v)%data_3d(i_s:i_e,:,j_s+1:j_e+1) - domain%vars_3d(domain%var_indx(kVARS%dzdy_v)%v)%data_3d(i_s:i_e,:,j_s:j_e)/domain%vars_3d(domain%var_indx(kVARS%jacobian_v)%v)%data_3d(i_s:i_e,:,j_s:j_e)))/mixed_denom

        B_coef = B_coef + X_coef
        C_coef = C_coef - X_coef

        deallocate(mixed_denom, X_coef, M_up, M_dwn)
    end subroutine update_coefs_host


    subroutine update_coefs_gpu()
        implicit none
        integer :: i, j, k
        real :: M_up_val, M_dwn_val, X_coef_val, mixed_denom_val
        real :: adv_dz_k, adv_dz_kp1, adv_dz_km1, sum_dz_up, sum_dz_dwn

        !$acc parallel loop gang vector collapse(3) &
        !$acc present(B_coef, C_coef, A_coef, dzdx, dzdy, jaco, alpha, dz_if, &
        !$acc         dzdx_surf, dzdy_surf, jaco_w, adv_dz_col, &
        !$acc         dzdx_u_stag, jaco_u_stag, dzdy_v_stag, jaco_v_stag)
        do j = j_s, j_e
            do k = k_s, k_e
                do i = i_s, i_e
                    adv_dz_k = adv_dz_col(k)

                    if (k == k_e) then
                        M_up_val = dzdx(i,k,j)**2 + dzdy(i,k,j)**2 + alpha(i,k,j)**2
                    else
                        adv_dz_kp1 = adv_dz_col(k+1)
                        sum_dz_up  = adv_dz_k + adv_dz_kp1
                        M_up_val = dzdx(i,k,j) * (dzdx(i,k,j)*adv_dz_kp1 + dzdx(i,k+1,j)*adv_dz_k) / sum_dz_up / jaco_w(i,k,j) &
                                 + dzdy(i,k,j) * (dzdy(i,k,j)*adv_dz_kp1 + dzdy(i,k+1,j)*adv_dz_k) / sum_dz_up / jaco_w(i,k,j) &
                                 + (alpha(i,k,j)**2*adv_dz_kp1 + alpha(i,k+1,j)**2*adv_dz_k) / sum_dz_up / jaco_w(i,k,j)
                    endif

                    if (k == k_s) then
                        M_dwn_val = dzdx(i,k,j)*dzdx_surf(i,j) + dzdy(i,k,j)*dzdy_surf(i,j) + alpha(i,k,j)**2
                    else
                        adv_dz_km1 = adv_dz_col(k-1)
                        sum_dz_dwn = adv_dz_k + adv_dz_km1
                        M_dwn_val = dzdx(i,k,j) * (dzdx(i,k,j)*adv_dz_km1 + dzdx(i,k-1,j)*adv_dz_k) / sum_dz_dwn / jaco_w(i,k-1,j) &
                                  + dzdy(i,k,j) * (dzdy(i,k,j)*adv_dz_km1 + dzdy(i,k-1,j)*adv_dz_k) / sum_dz_dwn / jaco_w(i,k-1,j) &
                                  + (alpha(i,k,j)**2*adv_dz_km1 + alpha(i,k-1,j)**2*adv_dz_k) / sum_dz_dwn / jaco_w(i,k-1,j)
                    endif

                    B_coef(i,k,j) = M_up_val  / (jaco(i,k,j) * adv_dz_k * dz_if(i_s,k+1,j_s))
                    C_coef(i,k,j) = M_dwn_val / (jaco(i,k,j) * adv_dz_k * dz_if(i_s,k,  j_s))

                    mixed_denom_val = (dz_if(i,k+1,j) + dz_if(i,k,j)) * 2.0 * dx
                    X_coef_val = -((dzdx_u_stag(i+1,k,j)/jaco_u_stag(i+1,k,j) - dzdx_u_stag(i,k,j)/jaco_u_stag(i,k,j)) + &
                                   (dzdy_v_stag(i,k,j+1)/jaco_v_stag(i,k,j+1) - dzdy_v_stag(i,k,j)/jaco_v_stag(i,k,j))) / mixed_denom_val

                    B_coef(i,k,j) = B_coef(i,k,j) + X_coef_val
                    C_coef(i,k,j) = C_coef(i,k,j) - X_coef_val
                    A_coef(i,k,j) = -4.0/(dx**2) - B_coef(i,k,j) - C_coef(i,k,j)
                enddo
            enddo
        enddo
    end subroutine update_coefs_gpu


    !>------------------------------------------------------------
    !! Save current module state into a cache slot via MOVE_ALLOC.
    !! Host arrays are transferred (no copy); device data must be
    !! removed first because OpenACC tracks attachments by host
    !! address. After this routine, the cache slot owns all the
    !! state; the module's allocatables are unallocated.
    !!------------------------------------------------------------
    subroutine save_to_cache(slot)
        implicit none
        integer, intent(in) :: slot

        ! Grid scalars
        domain_cache(slot)%i_s = i_s; domain_cache(slot)%i_e = i_e
        domain_cache(slot)%k_s = k_s; domain_cache(slot)%k_e = k_e
        domain_cache(slot)%j_s = j_s; domain_cache(slot)%j_e = j_e
        domain_cache(slot)%ims = ims; domain_cache(slot)%ime = ime
        domain_cache(slot)%jms = jms; domain_cache(slot)%jme = jme
        domain_cache(slot)%ids = ids; domain_cache(slot)%ide = ide
        domain_cache(slot)%jds = jds; domain_cache(slot)%jde = jde
        domain_cache(slot)%xs  = xs;  domain_cache(slot)%ys  = ys;  domain_cache(slot)%zs  = zs
        domain_cache(slot)%xm  = xm;  domain_cache(slot)%ym  = ym;  domain_cache(slot)%zm  = zm
        domain_cache(slot)%mx  = mx;  domain_cache(slot)%my  = my;  domain_cache(slot)%mz  = mz
        domain_cache(slot)%hs  = hs;  domain_cache(slot)%dx  = dx
        domain_cache(slot)%n_rows        = n_rows
        domain_cache(slot)%n_rows_global = n_rows_global
        domain_cache(slot)%solver_rank   = solver_rank
        domain_cache(slot)%east_neighbor  = east_neighbor
        domain_cache(slot)%west_neighbor  = west_neighbor
        domain_cache(slot)%north_neighbor = north_neighbor
        domain_cache(slot)%south_neighbor = south_neighbor
        domain_cache(slot)%solver_comm    = solver_comm
        domain_cache(slot)%structure_uploaded   = structure_uploaded
        domain_cache(slot)%wind_solver_max_iters = wind_solver_max_iters
        domain_cache(slot)%precond_n_sweeps      = precond_n_sweeps

        ! Sync device → host before removing device attachments
        !$acc wait
        if (structure_uploaded) then
            !$acc update host(A_coef, B_coef, C_coef, D_coef, E_coef, F_coef, G_coef, &
            !$acc             H_coef, I_coef, J_coef, K_coef, L_coef, M_coef, N_coef, O_coef)
            !$acc update host(x_sol, r_vec, r_hat, p_vec, p_hat, v_vec, &
            !$acc             s_vec, s_hat, t_vec, rhs, D_inv, prec_res)
            !$acc update host(east_send, east_recv, west_send, west_recv, &
            !$acc             north_send, north_recv, south_send, south_recv)
            !$acc exit data delete(A_coef, B_coef, C_coef, D_coef, E_coef, F_coef, G_coef, &
            !$acc                  H_coef, I_coef, J_coef, K_coef, L_coef, M_coef, N_coef, O_coef)
            !$acc exit data delete(x_sol, r_vec, r_hat, p_vec, p_hat, v_vec, &
            !$acc                  s_vec, s_hat, t_vec, rhs, D_inv, prec_res)
            !$acc exit data delete(east_send, east_recv, west_send, west_recv, &
            !$acc                  north_send, north_recv, south_send, south_recv)
#ifdef USE_NCCL
            if (allocated(sigma_dev)) then
                !$acc update host(sigma_dev, red5_dev, rho0_dev)
                !$acc exit data delete(sigma_dev, red5_dev, rho0_dev)
            endif
#endif
        endif
        if (allocated(dz_if)) then
            !$acc exit data delete(dz_if, div, alpha, dzdx, dzdy, jaco, dzdx_surf, dzdy_surf, sigma)
            !$acc exit data delete(jaco_w, dzdx_u_stag, jaco_u_stag, dzdy_v_stag, jaco_v_stag, adv_dz_col)
        endif

        ! Geometry / coefficient arrays
        if (allocated(A_coef))       call move_alloc(A_coef,       domain_cache(slot)%A_coef)
        if (allocated(B_coef))       call move_alloc(B_coef,       domain_cache(slot)%B_coef)
        if (allocated(C_coef))       call move_alloc(C_coef,       domain_cache(slot)%C_coef)
        if (allocated(D_coef))       call move_alloc(D_coef,       domain_cache(slot)%D_coef)
        if (allocated(E_coef))       call move_alloc(E_coef,       domain_cache(slot)%E_coef)
        if (allocated(F_coef))       call move_alloc(F_coef,       domain_cache(slot)%F_coef)
        if (allocated(G_coef))       call move_alloc(G_coef,       domain_cache(slot)%G_coef)
        if (allocated(H_coef))       call move_alloc(H_coef,       domain_cache(slot)%H_coef)
        if (allocated(I_coef))       call move_alloc(I_coef,       domain_cache(slot)%I_coef)
        if (allocated(J_coef))       call move_alloc(J_coef,       domain_cache(slot)%J_coef)
        if (allocated(K_coef))       call move_alloc(K_coef,       domain_cache(slot)%K_coef)
        if (allocated(L_coef))       call move_alloc(L_coef,       domain_cache(slot)%L_coef)
        if (allocated(M_coef))       call move_alloc(M_coef,       domain_cache(slot)%M_coef)
        if (allocated(N_coef))       call move_alloc(N_coef,       domain_cache(slot)%N_coef)
        if (allocated(O_coef))       call move_alloc(O_coef,       domain_cache(slot)%O_coef)
        if (allocated(div))          call move_alloc(div,          domain_cache(slot)%div)
        if (allocated(dz_if))        call move_alloc(dz_if,        domain_cache(slot)%dz_if)
        if (allocated(jaco))         call move_alloc(jaco,         domain_cache(slot)%jaco)
        if (allocated(dzdx))         call move_alloc(dzdx,         domain_cache(slot)%dzdx)
        if (allocated(dzdy))         call move_alloc(dzdy,         domain_cache(slot)%dzdy)
        if (allocated(sigma))        call move_alloc(sigma,        domain_cache(slot)%sigma)
        if (allocated(alpha))        call move_alloc(alpha,        domain_cache(slot)%alpha)
        if (allocated(dzdx_surf))    call move_alloc(dzdx_surf,    domain_cache(slot)%dzdx_surf)
        if (allocated(dzdy_surf))    call move_alloc(dzdy_surf,    domain_cache(slot)%dzdy_surf)
        if (allocated(jaco_w))       call move_alloc(jaco_w,       domain_cache(slot)%jaco_w)
        if (allocated(dzdx_u_stag))  call move_alloc(dzdx_u_stag,  domain_cache(slot)%dzdx_u_stag)
        if (allocated(jaco_u_stag))  call move_alloc(jaco_u_stag,  domain_cache(slot)%jaco_u_stag)
        if (allocated(dzdy_v_stag))  call move_alloc(dzdy_v_stag,  domain_cache(slot)%dzdy_v_stag)
        if (allocated(jaco_v_stag))  call move_alloc(jaco_v_stag,  domain_cache(slot)%jaco_v_stag)
        if (allocated(adv_dz_col))   call move_alloc(adv_dz_col,   domain_cache(slot)%adv_dz_col)

        ! Krylov / RHS / precond
        if (allocated(x_sol))    call move_alloc(x_sol,    domain_cache(slot)%x_sol)
        if (allocated(r_vec))    call move_alloc(r_vec,    domain_cache(slot)%r_vec)
        if (allocated(r_hat))    call move_alloc(r_hat,    domain_cache(slot)%r_hat)
        if (allocated(p_vec))    call move_alloc(p_vec,    domain_cache(slot)%p_vec)
        if (allocated(p_hat))    call move_alloc(p_hat,    domain_cache(slot)%p_hat)
        if (allocated(v_vec))    call move_alloc(v_vec,    domain_cache(slot)%v_vec)
        if (allocated(s_vec))    call move_alloc(s_vec,    domain_cache(slot)%s_vec)
        if (allocated(s_hat))    call move_alloc(s_hat,    domain_cache(slot)%s_hat)
        if (allocated(t_vec))    call move_alloc(t_vec,    domain_cache(slot)%t_vec)
        if (allocated(rhs))      call move_alloc(rhs,      domain_cache(slot)%rhs)
        if (allocated(D_inv))    call move_alloc(D_inv,    domain_cache(slot)%D_inv)
        if (allocated(prec_res)) call move_alloc(prec_res, domain_cache(slot)%prec_res)

        ! Halo buffers
        if (allocated(east_send))   call move_alloc(east_send,   domain_cache(slot)%east_send)
        if (allocated(east_recv))   call move_alloc(east_recv,   domain_cache(slot)%east_recv)
        if (allocated(west_send))   call move_alloc(west_send,   domain_cache(slot)%west_send)
        if (allocated(west_recv))   call move_alloc(west_recv,   domain_cache(slot)%west_recv)
        if (allocated(north_send))  call move_alloc(north_send,  domain_cache(slot)%north_send)
        if (allocated(north_recv))  call move_alloc(north_recv,  domain_cache(slot)%north_recv)
        if (allocated(south_send))  call move_alloc(south_send,  domain_cache(slot)%south_send)
        if (allocated(south_recv))  call move_alloc(south_recv,  domain_cache(slot)%south_recv)

#ifdef USE_NCCL
        if (allocated(sigma_dev))   call move_alloc(sigma_dev,   domain_cache(slot)%sigma_dev)
        if (allocated(red5_dev))    call move_alloc(red5_dev,    domain_cache(slot)%red5_dev)
        if (allocated(rho0_dev))    call move_alloc(rho0_dev,    domain_cache(slot)%rho0_dev)
#endif

        domain_cache(slot)%valid = .true.
    end subroutine save_to_cache


    !>------------------------------------------------------------
    !! Restore module state from a cache slot via MOVE_ALLOC.
    !! After this routine, the cache slot is empty and the module
    !! owns all the state. Device data is re-registered on GPU.
    !!------------------------------------------------------------
    subroutine restore_from_cache(slot)
        implicit none
        integer, intent(in) :: slot

        ! Grid scalars
        i_s = domain_cache(slot)%i_s; i_e = domain_cache(slot)%i_e
        k_s = domain_cache(slot)%k_s; k_e = domain_cache(slot)%k_e
        j_s = domain_cache(slot)%j_s; j_e = domain_cache(slot)%j_e
        ims = domain_cache(slot)%ims; ime = domain_cache(slot)%ime
        jms = domain_cache(slot)%jms; jme = domain_cache(slot)%jme
        ids = domain_cache(slot)%ids; ide = domain_cache(slot)%ide
        jds = domain_cache(slot)%jds; jde = domain_cache(slot)%jde
        xs  = domain_cache(slot)%xs;  ys  = domain_cache(slot)%ys;  zs  = domain_cache(slot)%zs
        xm  = domain_cache(slot)%xm;  ym  = domain_cache(slot)%ym;  zm  = domain_cache(slot)%zm
        mx  = domain_cache(slot)%mx;  my  = domain_cache(slot)%my;  mz  = domain_cache(slot)%mz
        hs  = domain_cache(slot)%hs;  dx  = domain_cache(slot)%dx
        n_rows        = domain_cache(slot)%n_rows
        n_rows_global = domain_cache(slot)%n_rows_global
        solver_rank   = domain_cache(slot)%solver_rank
        east_neighbor  = domain_cache(slot)%east_neighbor
        west_neighbor  = domain_cache(slot)%west_neighbor
        north_neighbor = domain_cache(slot)%north_neighbor
        south_neighbor = domain_cache(slot)%south_neighbor
        solver_comm    = domain_cache(slot)%solver_comm
        structure_uploaded   = domain_cache(slot)%structure_uploaded
        wind_solver_max_iters = domain_cache(slot)%wind_solver_max_iters
        precond_n_sweeps      = domain_cache(slot)%precond_n_sweeps

        ! Geometry / coefficient arrays
        if (allocated(domain_cache(slot)%A_coef))       call move_alloc(domain_cache(slot)%A_coef,       A_coef)
        if (allocated(domain_cache(slot)%B_coef))       call move_alloc(domain_cache(slot)%B_coef,       B_coef)
        if (allocated(domain_cache(slot)%C_coef))       call move_alloc(domain_cache(slot)%C_coef,       C_coef)
        if (allocated(domain_cache(slot)%D_coef))       call move_alloc(domain_cache(slot)%D_coef,       D_coef)
        if (allocated(domain_cache(slot)%E_coef))       call move_alloc(domain_cache(slot)%E_coef,       E_coef)
        if (allocated(domain_cache(slot)%F_coef))       call move_alloc(domain_cache(slot)%F_coef,       F_coef)
        if (allocated(domain_cache(slot)%G_coef))       call move_alloc(domain_cache(slot)%G_coef,       G_coef)
        if (allocated(domain_cache(slot)%H_coef))       call move_alloc(domain_cache(slot)%H_coef,       H_coef)
        if (allocated(domain_cache(slot)%I_coef))       call move_alloc(domain_cache(slot)%I_coef,       I_coef)
        if (allocated(domain_cache(slot)%J_coef))       call move_alloc(domain_cache(slot)%J_coef,       J_coef)
        if (allocated(domain_cache(slot)%K_coef))       call move_alloc(domain_cache(slot)%K_coef,       K_coef)
        if (allocated(domain_cache(slot)%L_coef))       call move_alloc(domain_cache(slot)%L_coef,       L_coef)
        if (allocated(domain_cache(slot)%M_coef))       call move_alloc(domain_cache(slot)%M_coef,       M_coef)
        if (allocated(domain_cache(slot)%N_coef))       call move_alloc(domain_cache(slot)%N_coef,       N_coef)
        if (allocated(domain_cache(slot)%O_coef))       call move_alloc(domain_cache(slot)%O_coef,       O_coef)
        if (allocated(domain_cache(slot)%div))          call move_alloc(domain_cache(slot)%div,          div)
        if (allocated(domain_cache(slot)%dz_if))        call move_alloc(domain_cache(slot)%dz_if,        dz_if)
        if (allocated(domain_cache(slot)%jaco))         call move_alloc(domain_cache(slot)%jaco,         jaco)
        if (allocated(domain_cache(slot)%dzdx))         call move_alloc(domain_cache(slot)%dzdx,         dzdx)
        if (allocated(domain_cache(slot)%dzdy))         call move_alloc(domain_cache(slot)%dzdy,         dzdy)
        if (allocated(domain_cache(slot)%sigma))        call move_alloc(domain_cache(slot)%sigma,        sigma)
        if (allocated(domain_cache(slot)%alpha))        call move_alloc(domain_cache(slot)%alpha,        alpha)
        if (allocated(domain_cache(slot)%dzdx_surf))    call move_alloc(domain_cache(slot)%dzdx_surf,    dzdx_surf)
        if (allocated(domain_cache(slot)%dzdy_surf))    call move_alloc(domain_cache(slot)%dzdy_surf,    dzdy_surf)
        if (allocated(domain_cache(slot)%jaco_w))       call move_alloc(domain_cache(slot)%jaco_w,       jaco_w)
        if (allocated(domain_cache(slot)%dzdx_u_stag))  call move_alloc(domain_cache(slot)%dzdx_u_stag,  dzdx_u_stag)
        if (allocated(domain_cache(slot)%jaco_u_stag))  call move_alloc(domain_cache(slot)%jaco_u_stag,  jaco_u_stag)
        if (allocated(domain_cache(slot)%dzdy_v_stag))  call move_alloc(domain_cache(slot)%dzdy_v_stag,  dzdy_v_stag)
        if (allocated(domain_cache(slot)%jaco_v_stag))  call move_alloc(domain_cache(slot)%jaco_v_stag,  jaco_v_stag)
        if (allocated(domain_cache(slot)%adv_dz_col))   call move_alloc(domain_cache(slot)%adv_dz_col,   adv_dz_col)

        ! Krylov / RHS / precond
        if (allocated(domain_cache(slot)%x_sol))    call move_alloc(domain_cache(slot)%x_sol,    x_sol)
        if (allocated(domain_cache(slot)%r_vec))    call move_alloc(domain_cache(slot)%r_vec,    r_vec)
        if (allocated(domain_cache(slot)%r_hat))    call move_alloc(domain_cache(slot)%r_hat,    r_hat)
        if (allocated(domain_cache(slot)%p_vec))    call move_alloc(domain_cache(slot)%p_vec,    p_vec)
        if (allocated(domain_cache(slot)%p_hat))    call move_alloc(domain_cache(slot)%p_hat,    p_hat)
        if (allocated(domain_cache(slot)%v_vec))    call move_alloc(domain_cache(slot)%v_vec,    v_vec)
        if (allocated(domain_cache(slot)%s_vec))    call move_alloc(domain_cache(slot)%s_vec,    s_vec)
        if (allocated(domain_cache(slot)%s_hat))    call move_alloc(domain_cache(slot)%s_hat,    s_hat)
        if (allocated(domain_cache(slot)%t_vec))    call move_alloc(domain_cache(slot)%t_vec,    t_vec)
        if (allocated(domain_cache(slot)%rhs))      call move_alloc(domain_cache(slot)%rhs,      rhs)
        if (allocated(domain_cache(slot)%D_inv))    call move_alloc(domain_cache(slot)%D_inv,    D_inv)
        if (allocated(domain_cache(slot)%prec_res)) call move_alloc(domain_cache(slot)%prec_res, prec_res)

        ! Halo buffers
        if (allocated(domain_cache(slot)%east_send))   call move_alloc(domain_cache(slot)%east_send,   east_send)
        if (allocated(domain_cache(slot)%east_recv))   call move_alloc(domain_cache(slot)%east_recv,   east_recv)
        if (allocated(domain_cache(slot)%west_send))   call move_alloc(domain_cache(slot)%west_send,   west_send)
        if (allocated(domain_cache(slot)%west_recv))   call move_alloc(domain_cache(slot)%west_recv,   west_recv)
        if (allocated(domain_cache(slot)%north_send))  call move_alloc(domain_cache(slot)%north_send,  north_send)
        if (allocated(domain_cache(slot)%north_recv))  call move_alloc(domain_cache(slot)%north_recv,  north_recv)
        if (allocated(domain_cache(slot)%south_send))  call move_alloc(domain_cache(slot)%south_send,  south_send)
        if (allocated(domain_cache(slot)%south_recv))  call move_alloc(domain_cache(slot)%south_recv,  south_recv)

#ifdef USE_NCCL
        if (allocated(domain_cache(slot)%sigma_dev))   call move_alloc(domain_cache(slot)%sigma_dev,   sigma_dev)
        if (allocated(domain_cache(slot)%red5_dev))    call move_alloc(domain_cache(slot)%red5_dev,    red5_dev)
        if (allocated(domain_cache(slot)%rho0_dev))    call move_alloc(domain_cache(slot)%rho0_dev,    rho0_dev)
#endif

        ! Re-register on device
        if (allocated(dz_if)) then
            !$acc enter data copyin(dz_if, div, alpha, dzdx, dzdy, jaco, dzdx_surf, dzdy_surf, sigma)
            !$acc enter data copyin(jaco_w, dzdx_u_stag, jaco_u_stag, dzdy_v_stag, jaco_v_stag, adv_dz_col)
        endif
        if (structure_uploaded) then
            !$acc enter data copyin(A_coef, B_coef, C_coef, D_coef, E_coef, F_coef, G_coef, &
            !$acc                   H_coef, I_coef, J_coef, K_coef, L_coef, M_coef, N_coef, O_coef)
            !$acc enter data copyin(x_sol, r_vec, r_hat, p_vec, p_hat, v_vec, &
            !$acc                   s_vec, s_hat, t_vec, rhs, D_inv, prec_res)
            !$acc enter data copyin(east_send, east_recv, west_send, west_recv, &
            !$acc                   north_send, north_recv, south_send, south_recv)
#ifdef USE_NCCL
            if (allocated(sigma_dev)) then
                !$acc enter data copyin(sigma_dev, red5_dev, rho0_dev)
            endif
#endif
        endif

        domain_cache(slot)%valid = .false.
    end subroutine restore_from_cache


end module wind_iterative
