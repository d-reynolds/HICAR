!> ----------------------------------------------------------------------------
!!  Standard advection scheme with variable order
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!  Dylan Reynolds (dylan.reynolds@slf.ch)
!!
!! ----------------------------------------------------------------------------
module adv_std
    use icar_constants
    use options_interface, only: options_t
    use domain_interface,  only: domain_t
    use adv_fluxcorr,      only: WRF_flux_corr
    use timer_interface,   only: timer_t
    use vertical_interpolation, only: find_match, weights   ! Zaengl Variant-B constant-z LUT
    use mpi_f08,           only: MPI_Allreduce, MPI_IN_PLACE, MPI_SUM, MPI_REAL, MPI_MIN, MPI_MAX, MPI_DOUBLE_PRECISION
    implicit none
    private

    integer :: ims, ime, jms, jme, kms, kme, its, ite, jts, jte
    integer :: i_s_w, i_e_w, j_s_w, j_e_w, i_s, i_e, j_s, j_e, horder, vorder, cz_diff_order
    !type(timer_t) :: flux_time, flux_up_time, flux_corr_time, sum_time
    ! For use advecting a (convective?) wind field
    ! real,dimension(:,:,:),allocatable :: U_4cu_u, V_4cu_u, W_4cu_u
    ! real,dimension(:,:,:),allocatable :: U_4cu_v, V_4cu_v, W_4cu_v
    real, dimension(:,:,:), allocatable   :: flux_x, flux_y, flux_z

    ! adv_theta_ref(i,k,j) = theta_bar(z(i,k,j))
    ! mean is a static snapshot, so this is filled once and freed on nest
    ! switch, mirroring the diagnostic arrays.
    real, dimension(:,:,:), allocatable   :: adv_theta_ref

    ! ----- constant-z horizontal-reconstruction LUT -----
    ! Layout: (ims:ime, kms:kme, jms:jme, 2). The pair-index (1=hi, 2=lo) is the
    ! slowest-varying dim so consecutive-i threads in a warp load fully coalesced.
    ! x: self = column i, west = column i-1 (u-face between i-1 and i).
    ! y: self = column j, south = column j-1 (v-face between j-1 and j).
    integer, dimension(:,:,:,:), allocatable :: zb_xs_lev, zb_xw_lev, zb_ys_lev, zb_yn_lev
    real,    dimension(:,:,:,:), allocatable :: zb_xs_wt,  zb_xw_wt,  zb_ys_wt,  zb_yn_wt
    ! Zaengl-2012 vertical alpha ramp: purely a function of k (=1 at surface,
    ! linearly -> 0 at band top). Stored 1D; previously was two 3D copies.
    real,    dimension(:),       allocatable :: zb_alpha_k
    ! F2 (zb_diff_order=4): wider-stencil constz reconstruction for the 4-point
    ! biharmonic. x: xw2 = column i-2, xe = column i+1 (to the u-face z).
    ! y: yn2 = column j-2, yp = column j+1 (to the v-face z). Built only when
    ! zb_corr/zb_diff; consumed only when zb_diff_order==4
    integer, dimension(:,:,:,:), allocatable :: zb_xw2_lev, zb_xe_lev, zb_yn2_lev, zb_yp_lev
    real,    dimension(:,:,:,:), allocatable :: zb_xw2_wt,  zb_xe_wt,  zb_yn2_wt,  zb_yp_wt

    ! Fine-mesh flux arrays (used by flux3_fm / sum_kernel_fm)
    real, dimension(:,:,:), allocatable   :: flux_x_fm, flux_y_fm, flux_z_fm

    public :: adv_std_init, adv_std_var_request, adv_std_advect3d, adv_std_compute_wind
    public :: adv_std_compute_wind_2d_fm, flux_2d_fm, sum_kernel_2d_fm, adv_std_clean_wind_arrays_fm
    public :: adv_std_clean_wind_arrays
    public :: flux_x_fm, flux_y_fm, flux_z_fm

    public :: adv_theta_ref, adv_std_apply_ref_vert

contains

    subroutine adv_std_init(domain,options)
        implicit none
        type(domain_t), intent(in) :: domain
        type(options_t),intent(in) :: options


        ims = domain%ims
        ime = domain%ime
        jms = domain%jms
        jme = domain%jme
        kms = domain%kms
        kme = domain%kme
        its = domain%its
        ite = domain%ite
        jts = domain%jts
        jte = domain%jte
        
        !set order of advection
        horder = options%adv%h_order
        vorder = options%adv%v_order

        !order of constant-z diffusion. If > 0, turned on, otherwise does nothing
        cz_diff_order = options%adv%cz_diff_order
        
        !Define bounds of advection computation. If using monotonic flux-limiter, it is necesarry to increase
        !advection bounds by 1. The necesarry extension of the halo is handeled in domain_object
        i_s = its
        i_e = ite
        j_s = jts
        j_e = jte

        i_s_w = i_s 
        i_e_w = i_e 
        j_s_w = j_s 
        j_e_w = j_e 

        if (options%adv%flux_corr==kFLUXCOR_MONO) then
            i_s = its - 1
            i_e = ite + 1
            j_s = jts - 1
            j_e = jte + 1
            
            ! wind arrays need to be extended by 1 in each direction for monotonic flux correction
            ! this allows for 2 upwind advection steps to be computed without needing a halo exchange
            i_s_w = i_s - 1
            i_e_w = i_e + 1
            j_s_w = j_s - 1
            j_e_w = j_e + 1
        endif
        
        call adv_std_init_theta_ref(domain)
        ! Zaengl const-z diffusion: build the constant-z LUT once per nest
        if (options%adv%cz_diff_order > 0) then
            call adv_std_init_zb_lut(domain)
        endif
    end subroutine

    subroutine adv_std_var_request(options)
        implicit none
        type(options_t), intent(inout) :: options

        ! List the variables that are required to be allocated for adv4 advection
        call options%alloc_vars( &
                        [kVARS%u,    kVARS%v,   kVARS%w,     kVARS%dz_interface, kVARS%water_vapor])

        ! List the variables that are required for restarts with adv4 advection
        call options%restart_vars( &
                        [kVARS%u,    kVARS%v,   kVARS%w,     kVARS%dz_interface, kVARS%water_vapor])

    end subroutine

    subroutine flux3(q,U_m,V_m,W_m,flux_x,flux_z,flux_y,t_factor,apply_cz_diff)
        implicit none
        real, dimension(ims:ime,  kms:kme,jms:jme),   intent(in)       :: q
        real, dimension(i_s_w:i_e_w+1,kms:kme,j_s_w:j_e_w+1), intent(in)       :: U_m, V_m, W_m
        real, dimension(i_s:i_e+1,kms:kme,j_s:j_e+1),   intent(out)    :: flux_x
        real, dimension(i_s:i_e+1,  kms:kme,j_s:j_e+1), intent(out)    :: flux_y
        real, dimension(i_s:i_e+1,  kms:kme+1,j_s:j_e+1), intent(out)    :: flux_z
        real, intent(in) :: t_factor
        ! Per-variable opt-in to constant-z diffusion. When .false. the upwind
        ! kernel keeps its implicit |u|-diffusion (diff_fac=1) and the cz_diff
        ! correction block is skipped — saves the bulk of the cz_diff cost for
        ! species that don't need horizontal-surface diffusion (everything
        ! except theta, and qv if opted in). Default .true. preserves legacy
        ! behavior of any caller that doesn't pass the argument.
        logical, intent(in), optional :: apply_cz_diff
        integer :: i, j, k
        real :: tmp, coef, u, v, w, q0, q1, q2, qn1, qn2, qn3, abs_u, t_factor_compact
        real :: qin1, qin2, qi1, qk1, qj1, qkn1, qkn2, qjn1, qjn2
        integer :: top, bot, bot2, nor, sou, sou2, eas, wes, wes2
        real :: ax, ay, q0z, qn1z, qn3z
        real :: diff_fac   ! 1.0 = keep scheme's implicit diffusion; 0.0 = even-order centered (zb_diff)
        real :: deta, dcz, dlim   ! Phase D: explicit constant-z diffusion
        real :: qm2z, qp1z                ! F3: wider-stencil constz values (4th-order)
        real :: gate                      ! branchless minmod sign-agreement gate
        logical :: cz_diff_on             ! effective cz_diff dispatch for this call

        ! GPU optimization variables. qc1..qc7 are stencil-cache *scalars*
        ! (one register each) rather than a small array — NVHPC spills
        ! private q_cache(:) arrays to local memory even with literal indices.
        real :: u_val, v_val, w_val, abs_u_val
        real :: qc1, qc2, qc3, qc4, qc5, qc6, qc7

        cz_diff_on = (cz_diff_order > 0)
        if (present(apply_cz_diff)) cz_diff_on = cz_diff_on .and. apply_cz_diff

        ! If cz_diff is active for this variable we skip the upwind kernel's
        ! implicit |u|-diffusion (diff_fac=0) and compute the constant-z
        ! diffusion below. Otherwise leave the scheme's original diffusion in.
        diff_fac = 1.0
        if (cz_diff_on) diff_fac = 0.0

        !$acc data present(q,U_m,V_m,W_m,flux_x,flux_y,flux_z)

        if (horder==1) then
            t_factor_compact = 0.5  * t_factor
            !$acc parallel async(1)
            ! Using tile(32,2) for better cache locality in 3D loop with moderate k dimension
            !$acc loop gang vector tile(32,4,1) private(u_val, v_val, abs_u_val, q0, qn1, qn3)
            do j = j_s, j_e
                do k = kms+1, kme
                    do i = i_s, i_e
                        ! Cache frequently accessed values to reduce memory traffic
                        q0 = q(i,k,j)
                        qn1 = q(i-1,k,j)
                        qn3 = q(i,k,j-1)

                        ! X-direction flux computation with cached values
                        u_val = U_m(i,k,j)
                        abs_u_val = ABS(u_val)
                        flux_x(i, k, j) = (u_val + abs_u_val) * qn1 + (u_val - abs_u_val) * q0
                        flux_x(i, k, j) = flux_x(i, k, j) * t_factor_compact

                        ! Y-direction flux computation with cached values
                        v_val = V_m(i,k,j)
                        abs_u_val = ABS(v_val)  ! reuse abs_u_val variable
                        flux_y(i, k, j) = (v_val + abs_u_val) * qn3 + (v_val - abs_u_val) * q0
                        flux_y(i, k, j) = flux_y(i, k, j) * t_factor_compact
                    enddo
                enddo
            enddo

            ! Using tile(32) for 2D loops - good balance for horizontal dimensions
            !$acc loop gang vector tile(32,4) private(u_val, abs_u_val, q0, qn1)
            do j = j_s, j_e
                do k = kms, kme
                    qn1 = q(i_e,k,j)
                    q0 = q(i_e+1,k,j)

                    u_val = U_m(i_e+1,k,j)
                    abs_u_val = ABS(u_val)
                    flux_x(i_e+1, k, j) = ((u_val + abs_u_val) * qn1 + (u_val - abs_u_val) * q0) * t_factor_compact
                enddo
            enddo

            !$acc loop gang vector tile(32,4) private(v_val, abs_u_val, q0, qn3)
            do k = kms, kme
                do i = i_s, i_e
                    qn3 = q(i,k,j_e)
                    q0 = q(i,k,j_e+1)

                    v_val = V_m(i,k,j_e+1)
                    abs_u_val = ABS(v_val)  ! reuse variable
                    flux_y(i, k, j_e+1) = ((v_val + abs_u_val) * qn3 + (v_val - abs_u_val) * q0) * t_factor_compact
                enddo
            enddo

            !$acc loop gang vector tile(64,2) private(u_val, v_val, abs_u_val, q0, qn1, qn3)
            do j = j_s, j_e
                do i = i_s, i_e
                    q0 = q(i,kms,j)
                    qn1 = q(i-1,kms,j)
                    qn3 = q(i,kms,j-1)

                    u_val = U_m(i,kms,j)
                    abs_u_val = ABS(u_val)
                    flux_x(i, kms, j) = ((u_val + abs_u_val) * qn1 + (u_val - abs_u_val) * q0) * t_factor_compact

                    v_val = V_m(i,kms,j)
                    abs_u_val = ABS(v_val)  ! reuse variable
                    flux_y(i, kms, j) = ((v_val + abs_u_val) * qn3 + (v_val - abs_u_val) * q0) * t_factor_compact
                enddo
            enddo
            !$acc end parallel

        else if (horder==3) then
            coef = (1./12) * t_factor

            ! Stencil scalars (qc1..qc7) kept in registers via named scalars;
            ! a private q_cache(:) array gets spilled to local memory by
            ! NVHPC even when indices are compile-time literals.
            !$acc parallel loop gang vector tile(32,4,1) async(1) &
            !$acc   private(u_val, abs_u_val, tmp, qc1, qc2, qc3, qc4, qc5, qc6, qc7)
            do j = j_s,j_e+1
            do k = kms,kme
            do i = i_s,i_e+1
                qc1 = q(i,k,j-2)  ! qn2 (j-2)
                qc2 = q(i,k,j-1)  ! qn1 (j-1)
                qc3 = q(i-2,k,j)  ! qn2 (i-2)
                qc4 = q(i-1,k,j)  ! qn1 (i-1)
                qc5 = q(i,k,j)    ! q0
                qc6 = q(i+1,k,j)  ! q1 (i+1)
                qc7 = q(i,k,j+1)  ! q1 (j+1)

                u_val = U_m(i,k,j)
                abs_u_val = ABS(u_val)
                tmp = 7.0 * (qc5 + qc4) - (qc6 + qc3)
                tmp = u_val * tmp
                tmp = tmp - diff_fac * abs_u_val * (3.0 * (qc5 - qc4) - (qc6 - qc3))
                flux_x(i,k,j) = tmp * coef

                u_val = V_m(i,k,j)
                abs_u_val = ABS(u_val)
                tmp = 7.0 * (qc5 + qc2) - (qc7 + qc1)
                tmp = u_val * tmp
                tmp = tmp - diff_fac * abs_u_val * (3.0 * (qc5 - qc2) - (qc7 - qc1))
                flux_y(i,k,j) = tmp * coef
            enddo
            enddo
            enddo

        else if (horder==5) then
            coef = (1./60)*t_factor
            !$acc parallel async(1)
            ! Using tile(32,4) for complex 5th order stencil computations
            !$acc loop gang vector tile(32,4,1) &
            !$acc   private(u_val, abs_u_val, tmp, qc1, qc2, qc3, qc4, qc5)
            do j = j_s,j_e
            do k = kms,kme
            do i = i_s,i_e+1
                qc1 = q(i-3,k,j)  ! qn3
                qc2 = q(i-2,k,j)  ! qn2
                qc3 = q(i-1,k,j)  ! qn1
                qc4 = q(i,k,j)    ! q0
                qc5 = q(i+1,k,j)  ! q1
                ! q2 = q(i+2,k,j) accessed directly when needed

                u_val = U_m(i,k,j)
                abs_u_val = ABS(u_val)

                ! 6th order centered flux: 37*(q0+qn1) - 8*(q1+qn2) + (q2+qn3)
                tmp = 37.0 * (qc4 + qc3) - 8.0 * (qc5 + qc2) + (q(i+2,k,j) + qc1)
                tmp = u_val * tmp
                ! 5th order diffusive correction (gated by diff_fac):
                ! abs(u) * (10*(q0-qn1) - 5*(q1-qn2) + (q2-qn3))
                tmp = tmp - diff_fac * abs_u_val * (10.0 * (qc4 - qc3) - 5.0 * (qc5 - qc2) + (q(i+2,k,j) - qc1))
                flux_x(i,k,j) = tmp * coef
            enddo
            enddo
            enddo
            
            !$acc loop gang vector tile(32,4,1) &
            !$acc   private(u_val, abs_u_val, tmp, qc1, qc2, qc3, qc4, qc5)
            do j = j_s,j_e+1
            do k = kms,kme
            do i = i_s,i_e
                qc1 = q(i,k,j-3)  ! qn3
                qc2 = q(i,k,j-2)  ! qn2
                qc3 = q(i,k,j-1)  ! qn1
                qc4 = q(i,k,j)    ! q0
                qc5 = q(i,k,j+1)  ! q1
                ! q2 = q(i,k,j+2) accessed directly when needed

                u_val = V_m(i,k,j)
                abs_u_val = ABS(u_val)

                tmp = 37.0 * (qc4 + qc3) - 8.0 * (qc5 + qc2) + (q(i,k,j+2) + qc1)
                tmp = u_val * tmp
                tmp = tmp - diff_fac * abs_u_val * (10.0 * (qc4 - qc3) - 5.0 * (qc5 - qc2) + (q(i,k,j+2) - qc1))
                flux_y(i,k,j) = tmp * coef
            enddo
            enddo
            enddo
            !$acc end parallel
        endif

        ! ---- CONSTANT-Z horizontal diffusion
        ! Replaces the implicit |u|-diffusion removed in earlier C (diff_fac=0) with
        ! the scheme's own diffusion, evaluated along constant-z via the zb LUT,
        ! alpha-ramp blended with the along-eta difference, then MINMOD-limited
        ! against the scheme's own (proven-monotone) eta-diffusion. Below split
        ! per-face for register pressure, branchless minmod, and the alpha=0
        ! top-slab fast path that skips all LUT loads. Whole block is gated on
        ! cz_diff_on, which combines cz_diff_order>0 with the per-variable
        ! apply_cz_diff dispatch (default off for non-theta species).
        if (cz_diff_on .and. cz_diff_order == 4) then
            ! ---- F3: faithful 4-point CONSTANT-Z biharmonic (zb_diff_order=4).
            ! The 3rd-order-upwind scheme's own d4 diffusive stencil
            !   3*(q_i-q_{i-1}) - (q_{i+1}-q_{i-2})
            ! evaluated constant-z (each of the 4 columns reconstructed to the
            ! face z via the F2 LUT), alpha-blended with along-eta, minmod-
            ! limited against the scheme's own (proven-stable) d4 eta-stencil.
            coef = (1.0/12.0) * t_factor

            ! x-face, k=kms..kme-1 (alpha > 0 region: LUT path)
            !$acc parallel loop gang vector tile(32,4,1) async(1) &
            !$acc   present(q,U_m,flux_x, &
            !$acc           zb_xs_lev,zb_xw_lev,zb_xs_wt,zb_xw_wt, &
            !$acc           zb_xw2_lev,zb_xe_lev,zb_xw2_wt,zb_xe_wt, &
            !$acc           zb_alpha_k) &
            !$acc   private(u_val,abs_u_val,ax,q0z,qn1z,qm2z,qp1z,deta,dcz,dlim,gate)
            do j = j_s, j_e
                do k = kms, kme-1
                    do i = i_s, i_e+1     ! match flux3 x-face: correct the i_e+1 face sum_kernel consumes
                        ax   = zb_alpha_k(k)
                        deta = 3.0*(q(i,k,j)-q(i-1,k,j)) - (q(i+1,k,j)-q(i-2,k,j))
                        q0z  = zb_xs_wt (i,k,j,1)*q(i,   zb_xs_lev (i,k,j,1), j) &
                             + zb_xs_wt (i,k,j,2)*q(i,   zb_xs_lev (i,k,j,2), j)
                        qn1z = zb_xw_wt (i,k,j,1)*q(i-1, zb_xw_lev (i,k,j,1), j) &
                             + zb_xw_wt (i,k,j,2)*q(i-1, zb_xw_lev (i,k,j,2), j)
                        qm2z = zb_xw2_wt(i,k,j,1)*q(i-2, zb_xw2_lev(i,k,j,1), j) &
                             + zb_xw2_wt(i,k,j,2)*q(i-2, zb_xw2_lev(i,k,j,2), j)
                        qp1z = zb_xe_wt (i,k,j,1)*q(i+1, zb_xe_lev (i,k,j,1), j) &
                             + zb_xe_wt (i,k,j,2)*q(i+1, zb_xe_lev (i,k,j,2), j)
                        dcz  = 3.0*(q0z-qn1z) - (qp1z-qm2z)
                        dlim = (1.0 - ax)*deta + ax*dcz
                        ! Branchless minmod: gate = 1 if dlim,deta same sign, else 0.
                        gate = 0.5*(1.0 + sign(1.0, dlim*deta))
                        dlim = gate * sign(min(abs(dlim), abs(deta)), deta)
                        u_val = U_m(i,k,j); abs_u_val = ABS(u_val)
                        flux_x(i,k,j) = flux_x(i,k,j) - abs_u_val * coef * dlim
                    enddo
                enddo
            enddo

            ! y-face, k=kms..kme-1 (alpha > 0 region: LUT path)
            !$acc parallel loop gang vector tile(32,4,1) async(1) &
            !$acc   present(q,V_m,flux_y, &
            !$acc           zb_ys_lev,zb_yn_lev,zb_ys_wt,zb_yn_wt, &
            !$acc           zb_yn2_lev,zb_yp_lev,zb_yn2_wt,zb_yp_wt, &
            !$acc           zb_alpha_k) &
            !$acc   private(v_val,abs_u_val,ay,q0z,qn3z,qm2z,qp1z,deta,dcz,dlim,gate)
            do j = j_s, j_e+1     ! match flux3 y-face: correct the j_e+1 face sum_kernel consumes
                do k = kms, kme-1
                    do i = i_s, i_e
                        ay   = zb_alpha_k(k)
                        deta = 3.0*(q(i,k,j)-q(i,k,j-1)) - (q(i,k,j+1)-q(i,k,j-2))
                        q0z  = zb_ys_wt (i,k,j,1)*q(i, zb_ys_lev (i,k,j,1), j  ) &
                             + zb_ys_wt (i,k,j,2)*q(i, zb_ys_lev (i,k,j,2), j  )
                        qn3z = zb_yn_wt (i,k,j,1)*q(i, zb_yn_lev (i,k,j,1), j-1) &
                             + zb_yn_wt (i,k,j,2)*q(i, zb_yn_lev (i,k,j,2), j-1)
                        qm2z = zb_yn2_wt(i,k,j,1)*q(i, zb_yn2_lev(i,k,j,1), j-2) &
                             + zb_yn2_wt(i,k,j,2)*q(i, zb_yn2_lev(i,k,j,2), j-2)
                        qp1z = zb_yp_wt (i,k,j,1)*q(i, zb_yp_lev (i,k,j,1), j+1) &
                             + zb_yp_wt (i,k,j,2)*q(i, zb_yp_lev (i,k,j,2), j+1)
                        dcz  = 3.0*(q0z-qn3z) - (qp1z-qm2z)
                        dlim = (1.0 - ay)*deta + ay*dcz
                        gate = 0.5*(1.0 + sign(1.0, dlim*deta))
                        dlim = gate * sign(min(abs(dlim), abs(deta)), deta)
                        v_val = V_m(i,k,j); abs_u_val = ABS(v_val)
                        flux_y(i,k,j) = flux_y(i,k,j) - abs_u_val * coef * dlim
                    enddo
                enddo
            enddo

            ! Top slab k=kme where alpha=0: dlim collapses to deta, so no LUT
            ! loads are needed. x and y fused (tiny per-thread work).
            !$acc parallel loop gang vector tile(32,1) async(1) &
            !$acc   present(q,U_m,V_m,flux_x,flux_y) &
            !$acc   private(u_val,v_val,abs_u_val,deta)
            do j = j_s, j_e+1     ! match flux3: correct the i_e+1/j_e+1 faces sum_kernel consumes
                do i = i_s, i_e+1
                    deta = 3.0*(q(i,kme,j)-q(i-1,kme,j)) - (q(i+1,kme,j)-q(i-2,kme,j))
                    u_val = U_m(i,kme,j); abs_u_val = ABS(u_val)
                    flux_x(i,kme,j) = flux_x(i,kme,j) - abs_u_val * coef * deta
                    deta = 3.0*(q(i,kme,j)-q(i,kme,j-1)) - (q(i,kme,j+1)-q(i,kme,j-2))
                    v_val = V_m(i,kme,j); abs_u_val = ABS(v_val)
                    flux_y(i,kme,j) = flux_y(i,kme,j) - abs_u_val * coef * deta
                enddo
            enddo
        else if (cz_diff_on .and. cz_diff_order == 2) then
            ! ---- 2nd-order Laplacian ----
            if (horder == 1) then
                coef = 0.5 * t_factor
            else if (horder == 3) then
                coef = 0.25 * t_factor
            else
                coef = (1.0/6.0) * t_factor
            endif

            ! x-face, k=kms..kme-1 (alpha > 0 region: LUT path)
            !$acc parallel loop gang vector tile(32,4,1) async(1) &
            !$acc   present(q,U_m,flux_x, &
            !$acc           zb_xs_lev,zb_xw_lev,zb_xs_wt,zb_xw_wt, &
            !$acc           zb_alpha_k) &
            !$acc   private(u_val,abs_u_val,ax,q0z,qn1z,deta,dcz,dlim,gate)
            do j = j_s, j_e
                do k = kms, kme-1
                    do i = i_s, i_e+1     ! match flux3 x-face: correct the i_e+1 face sum_kernel consumes
                        ax   = zb_alpha_k(k)
                        deta = q(i,k,j) - q(i-1,k,j)
                        q0z  = zb_xs_wt(i,k,j,1)*q(i,   zb_xs_lev(i,k,j,1), j) &
                             + zb_xs_wt(i,k,j,2)*q(i,   zb_xs_lev(i,k,j,2), j)
                        qn1z = zb_xw_wt(i,k,j,1)*q(i-1, zb_xw_lev(i,k,j,1), j) &
                             + zb_xw_wt(i,k,j,2)*q(i-1, zb_xw_lev(i,k,j,2), j)
                        dcz  = q0z - qn1z
                        dlim = (1.0 - ax)*deta + ax*dcz
                        gate = 0.5*(1.0 + sign(1.0, dlim*deta))
                        dlim = gate * sign(min(abs(dlim), abs(deta)), deta)
                        u_val = U_m(i,k,j); abs_u_val = ABS(u_val)
                        flux_x(i,k,j) = flux_x(i,k,j) - abs_u_val * coef * dlim
                    enddo
                enddo
            enddo

            ! y-face, k=kms..kme-1 (alpha > 0 region: LUT path)
            !$acc parallel loop gang vector tile(32,4,1) async(1) &
            !$acc   present(q,V_m,flux_y, &
            !$acc           zb_ys_lev,zb_yn_lev,zb_ys_wt,zb_yn_wt, &
            !$acc           zb_alpha_k) &
            !$acc   private(v_val,abs_u_val,ay,q0z,qn3z,deta,dcz,dlim,gate)
            do j = j_s, j_e+1     ! match flux3 y-face: correct the j_e+1 face sum_kernel consumes
                do k = kms, kme-1
                    do i = i_s, i_e
                        ay   = zb_alpha_k(k)
                        deta = q(i,k,j) - q(i,k,j-1)
                        q0z  = zb_ys_wt(i,k,j,1)*q(i, zb_ys_lev(i,k,j,1), j) &
                             + zb_ys_wt(i,k,j,2)*q(i, zb_ys_lev(i,k,j,2), j)
                        qn3z = zb_yn_wt(i,k,j,1)*q(i, zb_yn_lev(i,k,j,1), j-1) &
                             + zb_yn_wt(i,k,j,2)*q(i, zb_yn_lev(i,k,j,2), j-1)
                        dcz  = q0z - qn3z
                        dlim = (1.0 - ay)*deta + ay*dcz
                        gate = 0.5*(1.0 + sign(1.0, dlim*deta))
                        dlim = gate * sign(min(abs(dlim), abs(deta)), deta)
                        v_val = V_m(i,k,j); abs_u_val = ABS(v_val)
                        flux_y(i,k,j) = flux_y(i,k,j) - abs_u_val * coef * dlim
                    enddo
                enddo
            enddo

            ! Top slab k=kme where alpha=0: dlim = deta. No LUT loads needed.
            !$acc parallel loop gang vector tile(32,1) async(1) &
            !$acc   present(q,U_m,V_m,flux_x,flux_y) &
            !$acc   private(u_val,v_val,abs_u_val,deta)
            do j = j_s, j_e+1     ! match flux3: correct the i_e+1/j_e+1 faces sum_kernel consumes
                do i = i_s, i_e+1
                    deta = q(i,kme,j) - q(i-1,kme,j)
                    u_val = U_m(i,kme,j); abs_u_val = ABS(u_val)
                    flux_x(i,kme,j) = flux_x(i,kme,j) - abs_u_val * coef * deta
                    deta = q(i,kme,j) - q(i,kme,j-1)
                    v_val = V_m(i,kme,j); abs_u_val = ABS(v_val)
                    flux_y(i,kme,j) = flux_y(i,kme,j) - abs_u_val * coef * deta
                enddo
            enddo
        endif

        if (vorder==1) then
            t_factor_compact = 0.5  * t_factor
            !$acc parallel async(2)
            ! Using tile(32,4,1) for vertical flux computation
            !$acc loop gang vector tile(32,4,1) private(w_val, abs_u_val, q0, qn1)
            do j = j_s,j_e
               do k = kms+1,kme
                   do i = i_s,i_e
                       ! Cache values to reduce memory accesses
                       w_val = W_m(i,k-1,j)
                       abs_u_val = ABS(w_val)
                       q0 = q(i,k,j)
                       qn1 = q(i,k-1,j)
                       
                       flux_z(i,k,j) = ((w_val + abs_u_val) * qn1 + (w_val - abs_u_val) * q0) * t_factor_compact
                   enddo
               enddo
            enddo
            !$acc loop gang vector tile(64,2) private(w_val, q0)
            do j = j_s,j_e
                do i = i_s,i_e
                    ! flux_z(i,kms,j) = 0
                    q0 = q(i,kme,j)
                    w_val = W_m(i,kme,j)
                    flux_z(i,kme+1,j) = q0 * w_val * t_factor
                enddo
            enddo
            !$acc end parallel
        else if (vorder==3) then
            coef = (1./12)*t_factor
            t_factor_compact = 0.5 * t_factor
            !$acc parallel async(2)

            ! Interior k=kms+2..kme-1 — branchless 4-pt 3rd-order stencil
            !$acc loop gang vector tile(32,4,1) private(u, qn1, q0, q1, qn2, tmp)
            do j = j_s,j_e
            do k = kms+2,kme-1
            do i = i_s,i_e
                u   = W_m(i,k-1,j)
                qn2 = q(i,k-2,j)
                qn1 = q(i,k-1,j)
                q0  = q(i,k,j)
                q1  = q(i,k+1,j)
                tmp = 7*(q0+qn1) - (q1+qn2)
                tmp = u*tmp
                tmp = tmp - abs(u) * (3 * (q0 - qn1) - (q1 - qn2))
                flux_z(i,k,j) = tmp*coef
            enddo
            enddo
            enddo

            ! Boundary slabs k=kms+1 and k=kme — 2-pt upwind (insufficient stencil for 3rd order)
            !$acc loop gang vector tile(64,2) private(u, qn1, q0)
            do j = j_s,j_e
            do i = i_s,i_e
                u = W_m(i,kms,j)
                qn1 = q(i,kms,j)
                q0  = q(i,kms+1,j)
                flux_z(i,kms+1,j) = ((u + ABS(u)) * qn1 + (u - ABS(u)) * q0) * t_factor_compact

                u = W_m(i,kme-1,j)
                qn1 = q(i,kme-1,j)
                q0  = q(i,kme,j)
                flux_z(i,kme,j) = ((u + ABS(u)) * qn1 + (u - ABS(u)) * q0) * t_factor_compact
            enddo
            enddo

            ! Top slab k=kme+1 — outflow-only
            !$acc loop gang vector tile(64,2) private(u, qn1)
            do j = j_s,j_e
            do i = i_s,i_e
                u = W_m(i,kme,j)
                qn1 = q(i,kme,j)
                flux_z(i,kme+1,j) = qn1 * u * t_factor
            enddo
            enddo

            !$acc end parallel
        else if (vorder==5) then
            coef = (1./60)*t_factor
            !$acc parallel async(2)
            ! Using tile(32,2) for complex 5th order vertical computations
            !$acc loop gang vector tile(32,4,1) &
            !$acc   private(u_val, abs_u_val, tmp, qc1, qc2, qc3, qc4, qc5)
            do j = j_s,j_e
            do k = kms+3,kme-2
            do i = i_s,i_e
                qc1 = q(i,k-3,j)  ! qn3
                qc2 = q(i,k-2,j)  ! qn2
                qc3 = q(i,k-1,j)  ! qn1
                qc4 = q(i,k,j)    ! q0
                qc5 = q(i,k+1,j)  ! q1
                ! q2 = q(i,k+2,j) accessed directly when needed

                u_val = W_m(i,k-1,j)
                abs_u_val = ABS(u_val)

                tmp = 37.0 * (qc4 + qc3) - 8.0 * (qc5 + qc2) + (q(i,k+2,j) + qc1)
                tmp = u_val * tmp
                tmp = tmp - abs_u_val * (10.0 * (qc4 - qc3) - 5.0 * (qc5 - qc2) + (q(i,k+2,j) - qc1))
                flux_z(i,k,j) = tmp * coef
            enddo
            enddo
            enddo
            
            coef = (1./12)*t_factor
            !$acc loop gang vector tile(64,2) &
            !$acc   private(u_val, abs_u_val, tmp, q0, qn1, qc1, qc2, qc3, qc4)
            do j = j_s,j_e
                do i = i_s,i_e

                    ! flux_z(i,kms,j) = 0

                    ! Boundary k=kms+1: 1st-order upwind (insufficient stencil)
                    u_val = W_m(i,kms,j)
                    abs_u_val = ABS(u_val)
                    q0 = q(i,kms+1,j)
                    qn1 = q(i,kms,j)
                    flux_z(i,kms+1,j) = ((u_val + abs_u_val) * qn1 + (u_val - abs_u_val) * q0) * 0.5 * t_factor

                    ! Use 3rd order for kms+2 cell
                    u_val = W_m(i,kms+1,j)
                    abs_u_val = ABS(u_val)
                    qc1 = q(i,kms,j)    ! qn2
                    qc2 = q(i,kms+1,j)  ! qn1
                    qc3 = q(i,kms+2,j)  ! q0
                    qc4 = q(i,kms+3,j)  ! q1

                    tmp = 7.0 * (qc3 + qc2) - (qc4 + qc1)
                    tmp = u_val * tmp
                    tmp = tmp - abs_u_val * (3.0 * (qc3 - qc2) - (qc4 - qc1))
                    flux_z(i,kms+2,j) = tmp * coef

                    ! Use 3rd order for kme-1 cell
                    u_val = W_m(i,kme-2,j)
                    abs_u_val = ABS(u_val)
                    qc1 = q(i,kme-3,j)  ! qn2
                    qc2 = q(i,kme-2,j)  ! qn1
                    qc3 = q(i,kme-1,j)  ! q0
                    qc4 = q(i,kme,j)    ! q1

                    tmp = 7.0 * (qc3 + qc2) - (qc4 + qc1)
                    tmp = u_val * tmp
                    tmp = tmp - abs_u_val * (3.0 * (qc3 - qc2) - (qc4 - qc1))
                    flux_z(i,kme-1,j) = tmp * coef

                    u_val = W_m(i,kme-1,j)
                    abs_u_val = ABS(u_val)
                    q0 = q(i,kme,j)
                    qn1 = q(i,kme-1,j)
                    flux_z(i,kme,j) = ((u_val + abs_u_val) * qn1 + (u_val - abs_u_val) * q0) * 0.5 * t_factor

                    flux_z(i,kme+1,j) = q(i,kme,j) * W_m(i,kme,j) * t_factor
                enddo
            enddo
            !$acc end parallel
        endif

        !$acc wait(1,2)
        !$acc end data
    end subroutine flux3

    subroutine adv_std_advect3d(qfluxes,qold,U_m,V_m,W_m,denom,dz, flux_time, flux_corr_time, sum_time, t_factor_in,flux_corr_in,q_id_in,apply_cz_diff_in)
        ! !DIR$ INLINEALWAYS adv_std_advect3d
        implicit none
        real, dimension(ims:ime,  kms:kme,jms:jme),  intent(inout)   :: qfluxes
        real, dimension(ims:ime,  kms:kme,jms:jme),  intent(in)      :: qold
        real, dimension(ims:ime,  kms:kme,jms:jme),  intent(in)      :: dz, denom
        real, dimension(i_s_w:i_e_w+1,kms:kme,j_s_w:j_e_w+1), intent(in)       :: U_m, V_m, W_m
        type(timer_t), optional, intent(inout) :: flux_time, flux_corr_time, sum_time
        real, optional,                              intent(in)      :: t_factor_in
        integer, optional,                           intent(in)      :: flux_corr_in, q_id_in
        logical, optional,                           intent(in)      :: apply_cz_diff_in

        ! interal parameters
        real    :: t_factor, flux_diff_x, flux_diff_y, denom_val, dz_val
        integer :: i, k, j, flux_corr, q_id
        logical :: apply_cz_diff

        
        q_id = 1
        if (present(q_id_in)) q_id = q_id_in
        !Initialize t_factor, which is used during RK time stepping to scale the time step
        t_factor = 1.0
        if (present(t_factor_in)) t_factor = t_factor_in
        
        flux_corr = 0
        if (present(flux_corr_in)) flux_corr = flux_corr_in

        ! Per-variable cz_diff dispatch. Default .true. preserves legacy
        ! behavior of any caller that doesn't pass it (current callers should
        ! always pass it explicitly for clarity).
        apply_cz_diff = .true.
        if (present(apply_cz_diff_in)) apply_cz_diff = apply_cz_diff_in

        ! Choose optimized path based on whether flux correction is needed
        if (flux_corr == 0) then
            ! if ((horder==vorder) .or. (horder==5 .and. vorder==3)) then
            !     ! FUSED PATH: No flux correction needed
            !     ! Compute fluxes and directly update qfluxes

            !     call flux_time%start()
            !     !$acc wait(q_id) !qold is needed for following function
            !     call flux_and_advect_fused(qfluxes,qold,U_m,V_m,W_m,denom,dz,t_factor,q_id)
            !     call flux_time%stop()
            ! else
                if(present(flux_time)) call flux_time%start()
                !$acc wait(q_id) !qold is needed for following function
                call flux3(qfluxes,U_m,V_m,W_m,flux_x,flux_z,flux_y,t_factor,apply_cz_diff=apply_cz_diff)
                if(present(flux_time)) call flux_time%stop()

                if(present(sum_time)) call sum_time%start()
                call sum_kernel(flux_x, flux_y, flux_z, qold, qfluxes, denom, dz, q_id)
                if(present(sum_time)) call sum_time%stop()
            ! endif
            
        else
            ! STANDARD PATH: Flux correction enabled
            ! Must store fluxes for correction step
            if(present(flux_time)) call flux_time%start()
            call flux3(qfluxes,U_m, V_m, W_m, flux_x,flux_z,flux_y,t_factor,apply_cz_diff=apply_cz_diff)
            !$acc wait(q_id) !qold is not needed until after this point
            if(present(flux_time)) call flux_time%stop()
            
            if(present(flux_corr_time)) call flux_corr_time%start()
            ! Use async version which waits on q_id+100 then applies corrections
            call WRF_flux_corr(qold,U_m, V_m, W_m, flux_x,flux_z,flux_y,dz,denom,(q_id+100))
            if(present(flux_corr_time)) call flux_corr_time%stop()

            if(present(sum_time)) call sum_time%start()
            call sum_kernel(flux_x, flux_y, flux_z, qold, qfluxes, denom, dz, q_id)
            if(present(sum_time)) call sum_time%stop()
        endif

    end subroutine adv_std_advect3d

    !>------------------------------------------------------------
    !! Fix 3: apply theta_bar's vertical transport after the
    !! full RK3 loop AND the FCT flux corrector
    !!------------------------------------------------------------
    subroutine adv_std_apply_ref_vert(qprime, w_real, zheight, dt)
        implicit none
        real, dimension(ims:ime,kms:kme,jms:jme), intent(inout) :: qprime
        real, dimension(ims:ime,kms:kme,jms:jme), intent(in)    :: w_real, zheight
        real, intent(in) :: dt
        integer :: i, j, k, ku, kd
        real :: dthdz

        ! Branchless interior + boundary: clamp k+1 / k-1 to [kms,kme]. At
        ! k=kms (kd=kms, ku=kms+1) and k=kme (ku=kme, kd=kme-1) this collapses
        ! to the one-sided form; interior remains centered.
        !$acc parallel loop gang vector collapse(3) default(present) private(dthdz,ku,kd)
        do j = jts, jte
            do k = kms, kme
                do i = its, ite
                    ku = min(k+1, kme)
                    kd = max(k-1, kms)
                    dthdz = (adv_theta_ref(i,ku,j) - adv_theta_ref(i,kd,j)) &
                          / (zheight(i,ku,j) - zheight(i,kd,j))
                    qprime(i,k,j) = qprime(i,k,j) &
                        - w_real(i,k,j) * dthdz * dt
                enddo
            enddo
        enddo

    end subroutine adv_std_apply_ref_vert


    subroutine sum_kernel(flux_x, flux_y, flux_z, qold, qfluxes, denom, dz, q_id)
        implicit none
        real, dimension(i_s:i_e+1,kms:kme,j_s:j_e+1), intent(in) :: flux_x, flux_y
        real, dimension(i_s:i_e+1,kms:kme+1,j_s:j_e+1), intent(in) :: flux_z
        real, dimension(ims:ime,kms:kme,jms:jme), intent(in) :: qold, denom, dz
        real, dimension(ims:ime,kms:kme,jms:jme), intent(out) :: qfluxes
        integer, intent(in) :: q_id

        integer :: i, j, k
        real :: flux_diff_z

        !$acc parallel loop gang vector async(q_id) tile(64,2,1) present(qfluxes, qold, flux_x, flux_y, flux_z, denom, dz) private(flux_diff_z)
        do j = jts, jte
            do k = kms, kme
                do i = its, ite
                    flux_diff_z = flux_z(i,k+1,j) - merge(0.0, flux_z(i,k,j), k==kms)
                    qfluxes(i,k,j) = qold(i,k,j) - ((flux_x(i+1,k,j) - flux_x(i,k,j) + &
                                        flux_y(i,k,j+1) - flux_y(i,k,j) + &
                                        flux_diff_z / dz(i,k,j)) * denom(i,k,j))
                enddo
            enddo
        enddo
        !$acc wait(q_id)

    end subroutine sum_kernel
    

    ! subroutine test_divergence(dz)
    !     implicit none
    !     real, intent(in) :: dz(ims:ime,kms:kme,jms:jme)

    !     real, allocatable :: du(:,:), dv(:,:), dw(:,:)
    !     integer :: i,j,k

    !     allocate(du(i_s:i_e,j_s:j_e))
    !     allocate(dv(i_s:i_e,j_s:j_e))
    !     allocate(dw(i_s:i_e,j_s:j_e))

    !     do concurrent (j = j_s:j_e, k = kms:kme, i = i_s:i_e)

    !         du(i,j) = (U_m(i+1,k,j)-U_m(i,k,j))
    !         dv(i,j) = (V_m(i,k,j+1)-V_m(i,k,j))
    !         if (k==kms) then
    !             dw(i,j) = (W_m(i,k,j))/dz(i,k,j)
    !         else
    !             dw(i,j) = (W_m(i,k,j)-W_m(i,k-1,j))/dz(i,k,j)
    !         endif
    !         if (abs(du(i,j) + dv(i,j) + dw(i,j)) > 1e-3) then
    !             print*,  i,k,j , abs(du(i,j) + dv(i,j) + dw(i,j))
    !             print*, "Winds are not balanced on entry to advect"
    !             !error stop
    !         endif
    !     enddo

    ! end subroutine test_divergence

    subroutine adv_std_compute_wind(u,v,w,density,jaco,jaco_u,jaco_v,jaco_w,dz, dx, options, dt, U_m, V_m, W_m, denom)
        implicit none

        real, dimension(ims:ime+1,kms:kme,jms:jme), intent(in) :: u, jaco_u
        real, dimension(ims:ime,kms:kme,jms:jme+1), intent(in) :: v, jaco_v
        real, dimension(ims:ime,kms:kme,jms:jme), intent(in) :: w, jaco_w, density, dz, jaco
        type(options_t),    intent(in)  :: options
        real, allocatable, intent(inout) :: U_m(:,:,:), V_m(:,:,:), W_m(:,:,:), denom(:,:,:)
        real,intent(in)::dt, dx
        
        integer :: i, j, k
        real, dimension(ims:ime,kms:kme,jms:jme) :: rho
        logical :: advect_density

        ! Allocate the arrays on first call; subsequent calls reuse them. Per-timestep
        ! alloc/free cycle was wasteful (and consistent with the wind solver lesson:
        ! per-call alloc churn hurts at scale). Caller's intent changed from `out` to
        ! `inout` to allow reuse; first call still allocates from an unallocated argument.
        ! Pair each host allocation with a one-shot `enter data create` — repeating
        ! `enter data create` on every call inflates the dynamic refcount and turns
        ! cleanup at nest-switch into a partial-present trap.
        if (.not. allocated(U_m))   then
            allocate(U_m   (i_s_w:i_e_w+1, kms:kme,   j_s_w:j_e_w+1))
            !$acc enter data create(U_m)
        endif
        if (.not. allocated(V_m))   then
            allocate(V_m   (i_s_w:i_e_w+1, kms:kme,   j_s_w:j_e_w+1))
            !$acc enter data create(V_m)
        endif
        if (.not. allocated(W_m))   then
            allocate(W_m   (i_s_w:i_e_w+1, kms:kme,   j_s_w:j_e_w+1))
            !$acc enter data create(W_m)
        endif
        if (.not. allocated(denom)) then
            allocate(denom (ims:ime,       kms:kme,   jms:jme))
            !$acc enter data create(denom)
        endif
        if (.not. allocated(flux_x)) then
            allocate(flux_x(i_s:i_e+1, kms:kme,   j_s:j_e+1))
            !$acc enter data create(flux_x)
        endif
        if (.not. allocated(flux_y)) then
            allocate(flux_y(i_s:i_e+1, kms:kme,   j_s:j_e+1))
            !$acc enter data create(flux_y)
        endif
        if (.not. allocated(flux_z)) then
            allocate(flux_z(i_s:i_e+1, kms:kme+1, j_s:j_e+1))
            !$acc enter data create(flux_z)
        endif

        advect_density = options%adv%advect_density


        !$acc data present(u,v,w,density,jaco,jaco_u,jaco_v,jaco_w,dz, dx, U_m, V_m, W_m, denom) create(rho)
        if (advect_density) then
            !$acc parallel loop gang vector collapse(3)
            do j = jms,jme
                do k = kms,kme
                    do i = ims,ime
                        rho(i,k,j) = density(i,k,j)  
                    enddo
                enddo
            enddo
        else
            !$acc parallel loop gang vector collapse(3)
            do j = jms,jme
                do k = kms,kme
                    do i = ims,ime
                        rho(i,k,j) = 1
                    enddo
                enddo
            enddo
        endif

        !Compute the denomenator for all of the flux summation terms here once

        !$acc parallel loop gang vector collapse(3) async(1)
        do j = jms,jme
            do k = kms,kme
                do i = ims,ime
                    denom(i,k,j) = 1/(rho(i,k,j)*jaco(i,k,j))
                enddo
            enddo
        enddo

        !$acc parallel loop gang vector tile(32,4,1) async(2)
        do j = j_s_w,j_e_w+1
            do k = kms,kme
                do i = i_s_w,i_e_w+1
                    U_m(i,k,j) = u(i,k,j) * dt * (rho(i,k,j)+rho(i-1,k,j))*0.5 * &
                        jaco_u(i,k,j) / dx
                
                    V_m(i,k,j) = v(i,k,j) * dt * (rho(i,k,j)+rho(i,k,j-1))*0.5 * &
                        jaco_v(i,k,j) / dx
                enddo
            enddo
        enddo

        !$acc parallel loop gang vector tile(32,4,1) async(3)
        do j = j_s_w,j_e_w+1
            do k = kms,kme-1
                do i = i_s_w,i_e_w+1
                    W_m(i,k,j) = w(i,k,j) * dt * jaco_w(i,k,j) * &
                        ( rho(i,k,j)*dz(i,k+1,j) + &
                        rho(i,k+1,j)*dz(i,k,j) ) / &
                        (dz(i,k,j)+dz(i,k+1,j))
                enddo
            enddo
        enddo
        
        !$acc parallel loop gang vector collapse(2) async(4)
        do j = j_s_w,j_e_w+1
            do i = i_s_w,i_e_w+1
                W_m(i,kme,j) = w(i,kme,j) * dt * jaco_w(i,kme,j) * rho(i,kme,j)
            enddo
        enddo
        !$acc wait(1,2,3,4)
        !$acc end data

    end subroutine adv_std_compute_wind



    !>------------------------------------------------------------
    !! Compute 3D wind Courant numbers for fine-mesh advection.
    !! Extends adv_std_compute_wind_horiz with vertical wind + Jacobian.
    !! All bounds are explicit (no module-level state).
    !!------------------------------------------------------------
    subroutine adv_std_compute_wind_2d_fm(u_cell, v_cell, rho, &
        jaco, jaco_u, jaco_v, dx, dt, &
        U_m, V_m, denom, &
        ims_l, ime_l, ks, ke, jms_l, jme_l, its_l, ite_l, jts_l, jte_l)
        implicit none
        integer, intent(in) :: ims_l, ime_l, ks, ke, jms_l, jme_l
        integer, intent(in) :: its_l, ite_l, jts_l, jte_l
        real, dimension(ims_l:ime_l, ks:ke, jms_l:jme_l), intent(in) :: u_cell, v_cell, rho
        real, dimension(ims_l:ime_l, ks:ke, jms_l:jme_l), intent(in) :: jaco, jaco_u, jaco_v
        real, intent(in) :: dx, dt
        ! U_m, V_m, denom must be pre-allocated by the caller and already on
        ! the device (via !$acc enter data create). The caller manages their
        ! lifecycle to avoid OpenACC dummy-argument tracking issues.
        real, intent(inout) :: U_m   (its_l-2:ite_l+3, ks:ke, jts_l-2:jte_l+3)
        real, intent(inout) :: V_m   (its_l-2:ite_l+3, ks:ke, jts_l-2:jte_l+3)
        real, intent(inout) :: denom (ims_l:ime_l,     ks:ke, jms_l:jme_l)

        integer :: i, j, k
        integer :: i_s_w_l, i_e_w_l, j_s_w_l, j_e_w_l
        integer :: i_s_l, i_e_l, j_s_l, j_e_l

        ! Extended bounds for flux correction:
        ! Flux computation range: its-1:ite+1  (extend interior by 1)
        ! Wind array range: its-2:ite+2        (extend flux range by 1)
        i_s_l   = its_l - 1
        i_e_l   = ite_l + 1
        j_s_l   = jts_l - 1
        j_e_l   = jte_l + 1

        i_s_w_l = max(its_l - 2, ims_l+1)  ! Ensure we don't go out of bounds
        i_e_w_l = min(ite_l + 2, ime_l-1)
        j_s_w_l = max(jts_l - 2, jms_l+1)
        j_e_w_l = min(jte_l + 2, jme_l-1)

        ! Allocate module-level fine-mesh flux arrays
        allocate(flux_x_fm(i_s_l:i_e_l+1, ks:ke,   j_s_l:j_e_l+1))
        allocate(flux_y_fm(i_s_l:i_e_l+1, ks:ke,   j_s_l:j_e_l+1))
        !$acc enter data create(flux_x_fm, flux_y_fm)

        ! Compute 1/(rho * jaco) denominator
        !$acc parallel loop gang vector collapse(3) default(present)
        do j = jms_l, jme_l
            do k = ks, ke
                do i = ims_l, ime_l
                    denom(i,k,j) = 1.0 / (rho(i,k,j) * jaco(i,k,j))
                enddo
            enddo
        enddo

        ! U_m: face-staggered in x (mirrors adv_std_compute_wind line 1203)
        !$acc parallel loop gang vector collapse(3) default(present)
        do j = j_s_w_l, j_e_w_l+1
            do k = ks, ke
                do i = i_s_w_l, i_e_w_l+1
                    U_m(i,k,j) = u_cell(i,k,j) * dt * &
                        0.5 * (rho(i-1,k,j) + rho(i,k,j)) * jaco_u(i,k,j) / dx
                enddo
            enddo
        enddo

        ! V_m: face-staggered in y (mirrors line 1206)
        !$acc parallel loop gang vector collapse(3) default(present)
        do j = j_s_w_l, j_e_w_l+1
            do k = ks, ke
                do i = i_s_w_l, i_e_w_l+1
                    V_m(i,k,j) = v_cell(i,k,j) * dt * &
                        0.5 * (rho(i,k,j-1) + rho(i,k,j)) * jaco_v(i,k,j) / dx
                enddo
            enddo
        enddo


    end subroutine adv_std_compute_wind_2d_fm


    !>------------------------------------------------------------
    !! 3rd-order flux computation for fine-mesh advection.
    !! Stores fluxes in module-level flux_x_fm, flux_y_fm, flux_z_fm.
    !! Mirrors flux3 h_order=3 / v_order=3 with explicit bounds
    !! and fine-mesh vertical BCs (zero flux bottom, zero-gradient top).
    !!------------------------------------------------------------
    subroutine flux_2d_fm(q, U_m, V_m, t_factor, &
        ims_l, ime_l, ks, ke, jms_l, jme_l, its_l, ite_l, jts_l, jte_l)
        implicit none
        integer, intent(in) :: ims_l, ime_l, ks, ke, jms_l, jme_l
        integer, intent(in) :: its_l, ite_l, jts_l, jte_l
        real, dimension(ims_l:ime_l, ks:ke, jms_l:jme_l), intent(in) :: q
        real, dimension(its_l-2:ite_l+3, ks:ke, jts_l-2:jte_l+3), intent(in) :: U_m, V_m
        real, intent(in) :: t_factor

        integer :: i, j, k
        integer :: i_s_l, i_e_l, j_s_l, j_e_l
        real :: coef, t_factor_up
        real :: u_val, v_val, w_val, abs_u_val, tmp
        real :: q0, qn1, qn2, q1

        ! Flux computation bounds (extended for flux correction)
        i_s_l = its_l - 1
        i_e_l = ite_l + 1
        j_s_l = jts_l - 1
        j_e_l = jte_l + 1

        if (horder == 1) then
            t_factor_up = 0.5 * t_factor

            ! ==========================================
            ! Horizontal fluxes (1st order upwind)
            ! ==========================================
            !$acc parallel loop gang vector collapse(3) default(present)
            do j = j_s_l, j_e_l+1
                do k = ks, ke
                    do i = i_s_l, i_e_l+1
                        ! X-direction flux
                        u_val = U_m(i,k,j)
                        abs_u_val = ABS(u_val)
                        flux_x_fm(i,k,j) = ((u_val + abs_u_val) * q(i-1,k,j) + &
                                            (u_val - abs_u_val) * q(i,k,j)) * t_factor_up

                        ! Y-direction flux
                        v_val = V_m(i,k,j)
                        abs_u_val = ABS(v_val)
                        flux_y_fm(i,k,j) = ((v_val + abs_u_val) * q(i,k,j-1) + &
                                            (v_val - abs_u_val) * q(i,k,j)) * t_factor_up
                    enddo
                enddo
            enddo
        else if (horder==3) then
            coef = (1.0/12.0) * t_factor

            ! ==========================================
            ! Horizontal fluxes (3rd order)
            ! ==========================================
            !$acc parallel loop gang vector collapse(3) default(present)
            do j = j_s_l, j_e_l+1
                do k = ks, ke
                    do i = i_s_l, i_e_l+1
                        ! X-direction flux
                        u_val = U_m(i,k,j)
                        abs_u_val = ABS(u_val)
                        tmp = 7.0 * (q(i,k,j) + q(i-1,k,j)) - (q(i+1,k,j) + q(i-2,k,j))
                        tmp = u_val * tmp
                        tmp = tmp - abs_u_val * (3.0 * (q(i,k,j) - q(i-1,k,j)) - (q(i+1,k,j) - q(i-2,k,j)))
                        flux_x_fm(i,k,j) = tmp * coef

                        ! Y-direction flux
                        v_val = V_m(i,k,j)
                        abs_u_val = ABS(v_val)
                        tmp = 7.0 * (q(i,k,j) + q(i,k,j-1)) - (q(i,k,j+1) + q(i,k,j-2))
                        tmp = v_val * tmp
                        tmp = tmp - abs_u_val * (3.0 * (q(i,k,j) - q(i,k,j-1)) - (q(i,k,j+1) - q(i,k,j-2)))
                        flux_y_fm(i,k,j) = tmp * coef
                    enddo
                enddo
            enddo
        else if (horder==5) then
            coef = (1./60)*t_factor

            ! ==========================================
            ! Horizontal fluxes (5th order)
            ! ==========================================
            !$acc parallel loop gang vector collapse(3) default(present)
            do j = j_s_l, j_e_l+1
                do k = ks, ke
                    do i = i_s_l, i_e_l+1
                        ! X-direction flux
                        u_val = U_m(i,k,j)
                        abs_u_val = ABS(u_val)
                        tmp = 37.0 * (q(i,k,j) + q(i-1,k,j)) - 8.0 * (q(i+1,k,j) + q(i-2,k,j)) + (q(i+2,k,j) + q(i-3,k,j))
                        tmp = u_val * tmp
                        tmp = tmp - abs_u_val * (10.0 * (q(i,k,j) - q(i-1,k,j)) - 5.0 * (q(i+1,k,j) - q(i-2,k,j)) + (q(i+2,k,j) - q(i-3,k,j)))
                        flux_x_fm(i,k,j) = tmp * coef

                        ! Y-direction flux
                        v_val = V_m(i,k,j)
                        abs_u_val = ABS(v_val)
                        tmp = 37.0 * (q(i,k,j) + q(i,k,j-1)) - 8.0 * (q(i,k,j+1) + q(i,k,j-2)) + (q(i,k,j+2) + q(i,k,j-3))
                        tmp = v_val * tmp
                        tmp = tmp - abs_u_val * (10.0 * (q(i,k,j) - q(i,k,j-1)) - 5.0 * (q(i,k,j+1) - q(i,k,j-2)) + (q(i,k,j+2) - q(i,k,j-3)))
                        flux_y_fm(i,k,j) = tmp * coef
                    enddo
                enddo
            enddo
        endif
        ! ==========================================
        ! Vertical fluxes (3rd order with fine-mesh BCs)
        ! ==========================================

        ! do j = j_s_l, j_e_l
        !     do i = i_s_l, i_e_l
        !         ! k=ks: zero flux at bottom (ground)
        !         flux_z_fm(i,ks,j) = 0.0

        !         ! k=ks+1: 1st-order upwind (insufficient stencil below)
        !         w_val = W_m(i,ks,j)
        !         abs_u_val = ABS(w_val)
        !         flux_z_fm(i,ks+1,j) = ((w_val + abs_u_val) * q(i,ks,j) + &
        !             (w_val - abs_u_val) * q(i,ks+1,j)) * t_factor_up
        !     enddo
        ! enddo

        ! ! Interior vertical fluxes (3rd order): ks+2 to ke-1
        ! if (ke >= ks+3) then
        !     do j = j_s_l, j_e_l
        !         do k = ks+2, ke-1
        !             do i = i_s_l, i_e_l
        !                 w_val = W_m(i,k-1,j)
        !                 abs_u_val = ABS(w_val)
        !                 q0  = q(i,k,j)
        !                 qn1 = q(i,k-1,j)
        !                 q1  = q(i,k+1,j)
        !                 qn2 = q(i,k-2,j)
        !                 tmp = 7.0 * (q0 + qn1) - (q1 + qn2)
        !                 tmp = w_val * tmp
        !                 tmp = tmp - abs_u_val * (3.0 * (q0 - qn1) - (q1 - qn2))
        !                 flux_z_fm(i,k,j) = tmp * coef
        !             enddo
        !         enddo
        !     enddo
        ! endif

        ! do j = j_s_l, j_e_l
        !     do i = i_s_l, i_e_l
        !         ! k=ke: 1st-order upwind (insufficient stencil above)
        !         if (ke > ks+1) then
        !             w_val = W_m(i,ke-1,j)
        !             abs_u_val = ABS(w_val)
        !             flux_z_fm(i,ke,j) = ((w_val + abs_u_val) * q(i,ke-1,j) + &
        !                 (w_val - abs_u_val) * q(i,ke,j)) * t_factor_up
        !         else
        !             flux_z_fm(i,ke,j) = 0.0
        !         endif

        !         ! k=ke+1: zero-gradient outflow at top
        !         ! Use W_m(ke) (top boundary interface)
        !         flux_z_fm(i,ke+1,j) = q(i,ke,j) * W_m(i,ke,j) * t_factor
        !     enddo
        ! enddo

    end subroutine flux_2d_fm


    !>------------------------------------------------------------
    !! Apply flux divergence for fine-mesh advection.
    !! Mirrors sum_kernel with explicit bounds and fine-mesh BCs.
    !!------------------------------------------------------------
    subroutine sum_kernel_2d_fm(qold, qfluxes, denom, &
        ims_l, ime_l, ks, ke, jms_l, jme_l, its_l, ite_l, jts_l, jte_l)
        implicit none
        integer, intent(in) :: ims_l, ime_l, ks, ke, jms_l, jme_l
        integer, intent(in) :: its_l, ite_l, jts_l, jte_l
        real, dimension(ims_l:ime_l, ks:ke, jms_l:jme_l), intent(in)  :: qold, denom
        real, dimension(ims_l:ime_l, ks:ke, jms_l:jme_l), intent(out) :: qfluxes

        integer :: i, j, k

        !$acc parallel loop gang vector collapse(3) default(present)
        do j = jts_l, jte_l
            do k = ks, ke
                do i = its_l, ite_l
                    qfluxes(i,k,j) = qold(i,k,j) - ((flux_x_fm(i+1,k,j) - flux_x_fm(i,k,j) + &
                                        flux_y_fm(i,k,j+1) - flux_y_fm(i,k,j)) * denom(i,k,j))
                enddo
            enddo
        enddo

    end subroutine sum_kernel_2d_fm


    !>------------------------------------------------------------
    !! Deallocate fine-mesh wind and flux arrays.
    !! Mirrors adv_std_clean_wind_arrays for the fine mesh.
    !!------------------------------------------------------------
    subroutine adv_std_clean_wind_arrays_fm()
        implicit none

        !$acc exit data delete(flux_x_fm, flux_y_fm)
        if (allocated(flux_x_fm)) deallocate(flux_x_fm)
        if (allocated(flux_y_fm)) deallocate(flux_y_fm)

    end subroutine adv_std_clean_wind_arrays_fm

    !>------------------------------------------------------------
    !! Deallocate persistent module-level flux arrays + device copies.
    !! Required on nest context switch: bounds (i_s/i_e/...) change in
    !! adv_std_init, but the persistent arrays remain at the previous
    !! nest's shape, leading to OOB on device.
    !!------------------------------------------------------------
    subroutine adv_std_clean_wind_arrays()
        implicit none

        ! `finalize` forces the dynamic reference count to zero regardless of how
        ! many `enter data create` calls accumulated in adv_std_compute_wind across
        ! timesteps. Plain `exit data delete` only decrements by 1, leaving stale
        ! device mappings tied to a host pointer that deallocate() is about to
        ! free — those leftovers later collide with other arrays the allocator
        ! places in the same host range ("partially present" fatal in rte_lw etc).
        if (allocated(flux_x)) then
            !$acc exit data delete(flux_x) finalize
            deallocate(flux_x)
        endif
        if (allocated(flux_y)) then
            !$acc exit data delete(flux_y) finalize
            deallocate(flux_y)
        endif
        if (allocated(flux_z)) then
            !$acc exit data delete(flux_z) finalize
            deallocate(flux_z)
        endif
        if (allocated(adv_theta_ref)) then
            !$acc exit data delete(adv_theta_ref) finalize
            deallocate(adv_theta_ref)
        endif
        if (allocated(zb_xs_lev)) then
            !$acc exit data delete(zb_xs_lev, zb_xw_lev, zb_ys_lev, zb_yn_lev, &
            !$acc                   zb_xs_wt,  zb_xw_wt,  zb_ys_wt,  zb_yn_wt,  &
            !$acc                   zb_alpha_k, &
            !$acc                   zb_xw2_lev, zb_xe_lev, zb_yn2_lev, zb_yp_lev, &
            !$acc                   zb_xw2_wt,  zb_xe_wt,  zb_yn2_wt,  zb_yp_wt) finalize
            deallocate(zb_xs_lev, zb_xw_lev, zb_ys_lev, zb_yn_lev, &
                       zb_xs_wt,  zb_xw_wt,  zb_ys_wt,  zb_yn_wt,  &
                       zb_alpha_k, &
                       zb_xw2_lev, zb_xe_lev, zb_yn2_lev, zb_yp_lev, &
                       zb_xw2_wt,  zb_xe_wt,  zb_yn2_wt,  zb_yp_wt)
        endif

    end subroutine adv_std_clean_wind_arrays

    !>------------------------------------------------------------
    !! Fix 3: build the PARAMETER-FREE reference profile theta_bar(z) =
    !! the domain horizontal-MEAN potential temperature binned by PHYSICAL
    !! HEIGHT (NOT by terrain-following index k — a constant-k surface spans
    !! valley floor to mountaintop, so a per-k mean smears air across a wide
    !! height range and is exact only for a linear profile). theta_bar is a
    !! single, horizontally-uniform-in-physical-space profile (the
    !! load-bearing property: zero horizontal advection of theta_bar);
    !! evaluated per column at the true height:
    !!   adv_theta_ref(i,k,j) = theta_bar( z(i,k,j) )
    !! Built from a TRUE global mean (MPI_Allreduce over domain%compute_comms
    !! of interior height-bin sums) so it is decomposition-invariant. z is
    !! static per nest and the mean is a one-time snapshot, so this is filled
    !! once; idempotent + freed on nest switch (adv_std_clean_wind_arrays),
    !! so it re-fills with the new nest's state after a context switch.
    !!------------------------------------------------------------
    subroutine adv_std_init_theta_ref(domain)
        implicit none
        type(domain_t), intent(in) :: domain
        integer :: zv, thv, i, j, k, nz, nbins, b, nc, m
        ! th_sum/z_sum/bcnt accumulate in DOUBLE precision so the per-bin mean is
        ! decomposition-invariant: single-precision summation of the same cells in
        ! a different partition order (GPU atomics + MPI_SUM grouping) drifts by
        ! ~1e-4 K, which taints theta = theta_bar + theta' and breaks bit-for-bit
        ! MPI reproducibility. Double accumulation shrinks that drift below a
        ! single-precision ULP, so th_c/zb_c (single) round identically on any rank
        ! count. zb_c/th_c stay single (the downstream interpolation is single).
        double precision, allocatable :: th_sum(:), z_sum(:), bcnt(:)
        real, allocatable :: zb_c(:), th_c(:)
        real :: zmin, zmax, dzbin, zt, w1

        if (allocated(adv_theta_ref)) return
        zv  = domain%var_indx(kVARS%z)%v
        thv = domain%var_indx(kVARS%potential_temperature)%v

        nz    = kme - kms + 1
        nbins = nz                      ! parameter-free: tie bin count to the
                                        ! model's own vertical resolution
        allocate(adv_theta_ref(ims:ime, kms:kme, jms:jme))
        allocate(th_sum(nbins), z_sum(nbins), bcnt(nbins), zb_c(nbins), th_c(nbins))
        !$acc enter data create(adv_theta_ref)

        associate(z  => domain%vars_3d(zv )%data_3d, &
                  th => domain%vars_3d(thv)%data_3d)

        ! ---- theta_bar(z) = horizontal mean of theta binned by PHYSICAL
        ! HEIGHT (NOT by terrain-following index k). A constant-k surface
        ! spans valley-floor to mountaintop, so a per-k mean blends air
        ! across a huge height range (coordinate smear, exact only for a
        ! linear profile). Binning by geometric height gives a clean,
        ! single, horizontally-uniform-in-physical-space theta_bar(z) (the
        ! load-bearing property: zero horizontal advection of theta_bar).
        ! Global physical-height range over interior real cells (its:ite,
        ! jts:jte; no halo double-count), MPI-reduced so the bin axis is
        ! decomposition-invariant.
        zmin =  huge(1.0)
        zmax = -huge(1.0)
        !$acc parallel loop gang vector collapse(3) default(present) &
        !$acc          reduction(min:zmin) reduction(max:zmax)
        do j = jts, jte
            do k = kms, kme
                do i = its, ite
                    zmin = min(zmin, z(i,k,j))
                    zmax = max(zmax, z(i,k,j))
                enddo
            enddo
        enddo
        call MPI_Allreduce(MPI_IN_PLACE, zmin, 1, MPI_REAL, MPI_MIN, domain%compute_comms)
        call MPI_Allreduce(MPI_IN_PLACE, zmax, 1, MPI_REAL, MPI_MAX, domain%compute_comms)

        dzbin = (zmax - zmin) / real(nbins)
        if (dzbin <= 0.0) dzbin = 1.0   ! defensive (degenerate single-height)

        ! ---- histogram: accumulate theta and z into fixed geometric-height
        ! bins (device atomics; once-per-nest init, not perf-critical), then
        ! a TRUE global sum via MPI_Allreduce over the compute communicator.
        th_sum = 0.0; z_sum = 0.0; bcnt = 0.0
        !$acc enter data copyin(th_sum, z_sum, bcnt)
        !$acc parallel loop gang vector collapse(3) default(present) &
        !$acc          present(th_sum,z_sum,bcnt) private(b,zt)
        do j = jts, jte
            do k = kms, kme
                do i = its, ite
                    zt = z(i,k,j)
                    b  = int((zt - zmin)/dzbin) + 1
                    if (b < 1)     b = 1
                    if (b > nbins) b = nbins
                    !$acc atomic update
                    th_sum(b) = th_sum(b) + th(i,k,j)
                    !$acc atomic update
                    z_sum(b)  = z_sum(b)  + zt
                    !$acc atomic update
                    bcnt(b)   = bcnt(b)   + 1.0
                enddo
            enddo
        enddo
        !$acc exit data copyout(th_sum, z_sum, bcnt)

        call MPI_Allreduce(MPI_IN_PLACE, th_sum, nbins, MPI_DOUBLE_PRECISION, MPI_SUM, domain%compute_comms)
        call MPI_Allreduce(MPI_IN_PLACE, z_sum,  nbins, MPI_DOUBLE_PRECISION, MPI_SUM, domain%compute_comms)
        call MPI_Allreduce(MPI_IN_PLACE, bcnt,   nbins, MPI_DOUBLE_PRECISION, MPI_SUM, domain%compute_comms)

        ! ---- compact NON-EMPTY bins into a monotone (zb_c, th_c) table.
        ! Bin mean height z_sum/bcnt (not the geometric bin centre) so each
        ! sample sits where the data actually is; disjoint ordered bins =>
        ! zb_c strictly increasing across the nc kept points.
        nc = 0
        do b = 1, nbins
            if (bcnt(b) > 0.0d0) then
                nc = nc + 1
                zb_c(nc) = real(z_sum(b)  / bcnt(b))
                th_c(nc) = real(th_sum(b) / bcnt(b))
            endif
        enddo

        ! ---- adv_theta_ref(i,k,j) = theta_bar at physical height z(i,k,j),
        ! linearly interpolated on (zb_c, th_c). Out-of-range -> nearest
        ! endpoint (zero-gradient clamp: weight-1 copy of an existing mean
        ! value -> always well-defined, no new extrema).
        !$acc parallel loop gang vector collapse(3) default(present) &
        !$acc          copyin(zb_c,th_c) private(m,zt,w1)
        do j = jms, jme
            do k = kms, kme
                do i = ims, ime
                    zt = z(i,k,j)
                    if (nc == 1) then
                        adv_theta_ref(i,k,j) = th_c(1)
                    else if (zt <= zb_c(1)) then
                        adv_theta_ref(i,k,j) = th_c(1)
                    else if (zt >= zb_c(nc)) then
                        adv_theta_ref(i,k,j) = th_c(nc)
                    else
                        m = 1
                        do while (m < nc-1 .and. zb_c(m+1) <= zt)
                            m = m + 1
                        enddo
                        if (zb_c(m+1) > zb_c(m)) then
                            w1 = (zt - zb_c(m)) / (zb_c(m+1) - zb_c(m))
                        else
                            w1 = 0.0
                        endif
                        adv_theta_ref(i,k,j) = (1.0 - w1)*th_c(m) + w1*th_c(m+1)
                    endif
                enddo
            enddo
        enddo
        end associate

        deallocate(th_sum, z_sum, bcnt, zb_c, th_c)

    end subroutine adv_std_init_theta_ref

    !>------------------------------------------------------------
    !! Zaengl Variant-B: build the constant-z reconstruction LUT +
    !! slope-limited (alpha) blend from static geometry. Native to
    !! adv_std (same module as flux3). Idempotent; freed on nest
    !! switch (adv_std_clean_wind_arrays). RUNTIME-INERT until flux3
    !! consumes it (Step 3). Prints a coverage summary for validation.
    !!------------------------------------------------------------
    subroutine adv_std_init_zb_lut(domain)
        implicit none
        type(domain_t),  intent(in) :: domain

        integer :: zv, i, j, k, ncov_x, ntot_x, ncov_y, ntot_y
        integer :: lh, ll
        real    :: tgt, ramp, wh, wl
        logical :: bself, bnbr

        ! Routine-local helper holding zc transposed to (k, i, j) layout so
        ! that resolve_col's assumed-shape colz argument receives a stride-1
        ! section. zc(i, :, j) is strided in column-major Fortran (stride =
        ! ime-ims+1), which NVHPC's acc-routine-seq runtime rejects; zc_kij(:,
        ! i, j) is contiguous. One-shot LUT init -> allocated, transposed,
        ! used, then freed within this subroutine.
        real, allocatable :: zc_kij(:,:,:)

        if (allocated(zb_xs_lev)) return

        zv = domain%var_indx(kVARS%z)%v

        if (.not. allocated(domain%geo_u%z) .or. .not. allocated(domain%geo_v%z)) then
            if (STD_OUT_PE) write(*,*) "ERROR adv_std_init_zb_lut: geo_u%z/geo_v%z not allocated."
            error stop
        endif

        allocate(zb_xs_lev(ims:ime,kms:kme,jms:jme,2), zb_xw_lev(ims:ime,kms:kme,jms:jme,2))
        allocate(zb_ys_lev(ims:ime,kms:kme,jms:jme,2), zb_yn_lev(ims:ime,kms:kme,jms:jme,2))
        allocate(zb_xs_wt (ims:ime,kms:kme,jms:jme,2), zb_xw_wt (ims:ime,kms:kme,jms:jme,2))
        allocate(zb_ys_wt (ims:ime,kms:kme,jms:jme,2), zb_yn_wt (ims:ime,kms:kme,jms:jme,2))
        allocate(zb_alpha_k(kms:kme))
        allocate(zb_xw2_lev(ims:ime,kms:kme,jms:jme,2), zb_xe_lev(ims:ime,kms:kme,jms:jme,2))
        allocate(zb_yn2_lev(ims:ime,kms:kme,jms:jme,2), zb_yp_lev(ims:ime,kms:kme,jms:jme,2))
        allocate(zb_xw2_wt (ims:ime,kms:kme,jms:jme,2), zb_xe_wt (ims:ime,kms:kme,jms:jme,2))
        allocate(zb_yn2_wt (ims:ime,kms:kme,jms:jme,2), zb_yp_wt (ims:ime,kms:kme,jms:jme,2))

        !$acc enter data create(zb_xs_lev, zb_xw_lev, zb_ys_lev, zb_yn_lev, &
        !$acc                    zb_xs_wt,  zb_xw_wt,  zb_ys_wt,  zb_yn_wt,  &
        !$acc                    zb_alpha_k, &
        !$acc                    zb_xw2_lev, zb_xe_lev, zb_yn2_lev, zb_yp_lev, &
        !$acc                    zb_xw2_wt,  zb_xe_wt,  zb_yn2_wt,  zb_yp_wt)

        ! Safe defaults everywhere: identity bracket, pure along-eta (alpha=0).
        !$acc kernels default(present)
        zb_xs_lev = kms; zb_xw_lev = kms; zb_ys_lev = kms; zb_yn_lev = kms
        zb_xs_wt(:,:,:,1) = 1.0; zb_xs_wt(:,:,:,2) = 0.0
        zb_xw_wt(:,:,:,1) = 1.0; zb_xw_wt(:,:,:,2) = 0.0
        zb_ys_wt(:,:,:,1) = 1.0; zb_ys_wt(:,:,:,2) = 0.0
        zb_yn_wt(:,:,:,1) = 1.0; zb_yn_wt(:,:,:,2) = 0.0
        zb_alpha_k = 0.0
        ! F2 wider-stencil defaults (identity): xw2/xe default to the i-1/i
        ! columns' fallback so the biharmonic outer difference degrades safely
        ! at build-loop edges (those edge faces are not in the consumption range).
        zb_xw2_lev = kms; zb_xe_lev = kms; zb_yn2_lev = kms; zb_yp_lev = kms
        zb_xw2_wt(:,:,:,1) = 1.0; zb_xw2_wt(:,:,:,2) = 0.0
        zb_xe_wt (:,:,:,1) = 1.0; zb_xe_wt (:,:,:,2) = 0.0
        zb_yn2_wt(:,:,:,1) = 1.0; zb_yn2_wt(:,:,:,2) = 0.0
        zb_yp_wt (:,:,:,1) = 1.0; zb_yp_wt (:,:,:,2) = 0.0
        !$acc end kernels
        ncov_x = 0; ntot_x = 0; ncov_y = 0; ntot_y = 0

        ! Zaengl-2012 vertical alpha ramp: =1 at surface (k=kms),
        ! linearly -> 0 at band top (k=kme). Removes the binary
        ! 0/1 hard edge that otherwise seeds grid-scale NaN.
        !$acc parallel loop present(zb_alpha_k)
        do k = kms, kme
            ramp = real(kme - k) / real(kme - kms)
            zb_alpha_k(k) = max(0.0, min(1.0, ramp))
        end do

        ! Transpose zc into a k-leading helper so column slices passed to
        ! resolve_col are stride-1 on the device. One-shot at LUT-init.
        allocate(zc_kij(kms:kme, ims:ime, jms:jme))
        !$acc enter data create(zc_kij)

        associate(zc => domain%vars_3d(zv)%data_3d)
        !$acc parallel loop gang vector collapse(3) default(present)
        do j = jms, jme
            do k = kms, kme
                do i = ims, ime
                    zc_kij(k, i, j) = zc(i, k, j)
                enddo
            enddo
        enddo
        end associate

        associate(geo_u_z => domain%geo_u%z, geo_v_z => domain%geo_v%z)
        !$acc parallel loop gang vector collapse(3) default(present) copyin(geo_u_z, geo_v_z)
        do j = jms, jme
            do k = kms, kme
                do i = ims, ime
                    ! x-direction: u-face (i,k,j) between cells (i-1) and (i)
                    if (i > ims) then
                        ntot_x = ntot_x + 1
                        tgt = geo_u_z(i,k,j)
                        call resolve_col(tgt, zc_kij(:,i,j),   kms, kme, lh, ll, wh, wl, bself)
                        zb_xs_lev(i,k,j,1) = lh; zb_xs_lev(i,k,j,2) = ll
                        zb_xs_wt (i,k,j,1) = wh; zb_xs_wt (i,k,j,2) = wl
                        call resolve_col(tgt, zc_kij(:,i-1,j), kms, kme, lh, ll, wh, wl, bnbr)
                        zb_xw_lev(i,k,j,1) = lh; zb_xw_lev(i,k,j,2) = ll
                        zb_xw_wt (i,k,j,1) = wh; zb_xw_wt (i,k,j,2) = wl
                        if (bself .and. bnbr) ncov_x = ncov_x + 1
                        ! F2: wider-stencil columns for the 4-point biharmonic
                        ! (i-2, i+1) reconstructed to the same u-face z. Guarded
                        ! at the memory edges (edge faces not consumed anyway).
                        if (i-2 >= ims) then
                            call resolve_col(tgt, zc_kij(:,i-2,j), kms, kme, lh, ll, wh, wl, bnbr)
                            zb_xw2_lev(i,k,j,1) = lh; zb_xw2_lev(i,k,j,2) = ll
                            zb_xw2_wt (i,k,j,1) = wh; zb_xw2_wt (i,k,j,2) = wl
                        endif
                        if (i+1 <= ime) then
                            call resolve_col(tgt, zc_kij(:,i+1,j), kms, kme, lh, ll, wh, wl, bnbr)
                            zb_xe_lev(i,k,j,1) = lh; zb_xe_lev(i,k,j,2) = ll
                            zb_xe_wt (i,k,j,1) = wh; zb_xe_wt (i,k,j,2) = wl
                        endif
                    endif
                    ! y-direction: v-face (i,k,j) between cells (j-1) and (j)
                    if (j > jms) then
                        ntot_y = ntot_y + 1
                        tgt = geo_v_z(i,k,j)
                        call resolve_col(tgt, zc_kij(:,i,j),   kms, kme, lh, ll, wh, wl, bself)
                        zb_ys_lev(i,k,j,1) = lh; zb_ys_lev(i,k,j,2) = ll
                        zb_ys_wt (i,k,j,1) = wh; zb_ys_wt (i,k,j,2) = wl
                        call resolve_col(tgt, zc_kij(:,i,j-1), kms, kme, lh, ll, wh, wl, bnbr)
                        zb_yn_lev(i,k,j,1) = lh; zb_yn_lev(i,k,j,2) = ll
                        zb_yn_wt (i,k,j,1) = wh; zb_yn_wt (i,k,j,2) = wl
                        if (bself .and. bnbr) ncov_y = ncov_y + 1
                        ! F2: wider-stencil columns for the 4-point biharmonic
                        ! (j-2, j+1) reconstructed to the same v-face z.
                        if (j-2 >= jms) then
                            call resolve_col(tgt, zc_kij(:,i,j-2), kms, kme, lh, ll, wh, wl, bnbr)
                            zb_yn2_lev(i,k,j,1) = lh; zb_yn2_lev(i,k,j,2) = ll
                            zb_yn2_wt (i,k,j,1) = wh; zb_yn2_wt (i,k,j,2) = wl
                        endif
                        if (j+1 <= jme) then
                            call resolve_col(tgt, zc_kij(:,i,j+1), kms, kme, lh, ll, wh, wl, bnbr)
                            zb_yp_lev(i,k,j,1) = lh; zb_yp_lev(i,k,j,2) = ll
                            zb_yp_wt (i,k,j,1) = wh; zb_yp_wt (i,k,j,2) = wl
                        endif
                    endif
                enddo
            enddo
        enddo
        end associate

        !$acc exit data delete(zc_kij)
        deallocate(zc_kij)

    end subroutine adv_std_init_zb_lut

    ! Resolve target physical height tgt_in against a 1-based column of
    ! monotone-increasing heights colz. q~ = w_hi*q(lev_hi)+w_lo*q(lev_lo):
    !  - true bracket -> linear interpolation
    !  - tgt below the column's lowest level -> zero-gradient bottom clamp
    !    (w=1 on q(kms); no new extrema -> proven-stable Phase-3 form)
    !  - tgt above the column top -> top-value clamp
    subroutine resolve_col(tgt_in, colz, kms_loc, kme_loc, lev_hi, lev_lo, w_hi, w_lo, bracketed)
        !$acc routine seq
        real,    intent(in)  :: tgt_in
        real,    intent(in)  :: colz(:)
        integer, intent(in)  :: kms_loc, kme_loc
        integer, intent(out) :: lev_hi, lev_lo
        real,    intent(out) :: w_hi, w_lo
        logical, intent(out) :: bracketed
        integer :: fmm(2), nz
        real    :: ww(2), dzc, t

        nz  = size(colz)
        fmm = find_match(tgt_in, colz)
        if (fmm(1) > 0) then
            ww     = weights(tgt_in, colz(fmm(2)), colz(fmm(1)))
            lev_hi = kms_loc + fmm(2) - 1;  lev_lo = kms_loc + fmm(1) - 1
            w_hi   = ww(1);             w_lo   = ww(2)
            bracketed = .true.
        else if (fmm(1) == -1) then
            ! below the column's lowest level
            bracketed = .false.
            if (nz >= 2) then
                dzc = colz(2) - colz(1)
                if (dzc <= 0.0) then
                    lev_hi = kms_loc;  lev_lo = kms_loc;  w_hi = 1.0;  w_lo = 0.0
                else
                    t = (tgt_in - colz(1)) / dzc      ! < 0
                    t = max(t, -1.0)                  ! slope-limit: <= 1 cell
                    lev_hi = kms_loc + 1;  lev_lo = kms_loc
                    w_hi   = t;        w_lo   = 1.0 - t
                endif
            else
                lev_hi = kms_loc;  lev_lo = kms_loc;  w_hi = 1.0;  w_lo = 0.0
            endif
        else
            ! above the column's top level
            bracketed = .false.
            if (nz >= 2) then
                dzc = colz(nz) - colz(nz-1)
                if (dzc <= 0.0) then
                    lev_hi = kme_loc;  lev_lo = kme_loc;  w_hi = 1.0;  w_lo = 0.0
                else
                    t = (tgt_in - colz(nz-1)) / dzc   ! > 1
                    t = min(t, 2.0)                   ! slope-limit: <= 1 cell
                    lev_hi = kme_loc;  lev_lo = kme_loc - 1
                    w_hi   = t;     w_lo   = 1.0 - t
                endif
            else
                lev_hi = kme_loc;  lev_lo = kme_loc;  w_hi = 1.0;  w_lo = 0.0
            endif
        endif
    end subroutine resolve_col

end module adv_std
