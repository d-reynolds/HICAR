!>------------------------------------------------------------
!! Prognostic RANS wind solver — momentum dynamics.
!!
!! Activated by namelist `wind = "RANS"` (windtype == kRANS_WINDS).
!! This module owns the momentum side of the RANS step:
!!
!!   rans_momentum_step:
!!     (1) WS-RK3 advection of u, v (staggered control volumes) and
!!         w_real (cell-centred, reuses the standard scalar kernel),
!!         with buoyancy g*theta'/theta_bar on w_real.
!!     (2) Horizontal Smagorinsky deformation diffusion of u, v, w_real.
!!     (3) Rayleigh damping of w_real near the model lid.
!!
!! The pressure projection that closes the time step lives in
!! wind.F90::rans_project (it needs calc_divergence / balance_uvw /
!! calc_w_real, which live there — keeping it in wind.F90 avoids a
!! module dependency cycle). Call sequence per model time step is in
!! advection_driver.F90::advect; the full formulation is documented in
!! docs/rans_solver_math.md.
!!
!! Design notes (see the math doc for the reasoning):
!!   - Momentum transport uses the SAME rho*J*dt/dx contravariant
!!     Courant-flux weighting as scalar advection. Vertical flux
!!     divergence is divided by the control-volume dz exactly as in
!!     adv_std's sum_kernel — no bare w_real/dz terms anywhere.
!!   - Transport is frozen at time n (the projected, divergence-free
!!     field), the same policy scalar RK3 uses, which makes flux form
!!     and advective form equivalent and conservative.
!!   - 3rd-order WRF upwind fluxes, degrading to 1st order where the
!!     stencil would leave the domain (global boundaries, top/bottom).
!!   - No total-pressure or total-temperature arithmetic: buoyancy is
!!     computed from theta' = theta - adv_theta_ref in single
!!     precision safely.
!!
!! Requirements checked at init (wind.F90::init_winds):
!!   - advection = "std" (kADV_STD) and h_order >= 3 (halo width 2).
!!------------------------------------------------------------
module wind_rans

    use icar_constants,    only : STD_OUT_PE, kVARS
    use domain_interface,  only : domain_t
    use options_interface, only : options_t
    use adv_std,           only : adv_std_advect3d, adv_theta_ref
    use mod_wrf_constants, only : gravity
    use timer_interface,   only : timer_t

    implicit none
    private

    public :: rans_momentum_step

    ! Horizontal Smagorinsky constant and explicit-diffusion stability clip
    ! (K*dt/dx^2 <= 0.12 keeps the forward-Euler 2D Laplacian comfortably
    ! inside its 0.25 limit). WRF uses C_s ~ 0.25 for its 2D Smagorinsky.
    real, parameter :: RANS_SMAG_CS       = 0.21
    real, parameter :: RANS_SMAG_CLIP     = 0.12

    ! Rayleigh damping of w_real: active above RANS_RAYLEIGH_FRAC of the
    ! column depth (above terrain), minimum relaxation time at the lid.
    real, parameter :: RANS_RAYLEIGH_FRAC = 0.75
    real, parameter :: RANS_RAYLEIGH_TAU  = 30.0   ! [s]

    real, parameter :: pi = 3.14159265358979

    ! Domain extents, refreshed on every call (multi-nest safety; the
    ! module holds no per-nest state).
    integer :: ims, ime, jms, jme, kms, kme
    integer :: its, ite, jts, jte, kts, kte

contains

    !>--------------------------------------------------------
    !! One full momentum dynamics step (everything except the
    !! projection). Called from advection_driver::advect after
    !! scalar advection; U_m/V_m/W_m/denom are the scalar
    !! transport arrays built from the time-n winds.
    !!--------------------------------------------------------
    subroutine rans_momentum_step(domain, options, dt, U_m, V_m, W_m, denom, &
                                  flux_time, flux_corr_time, sum_time)
        implicit none
        type(domain_t),  intent(inout) :: domain
        type(options_t), intent(in)    :: options
        real,            intent(in)    :: dt
        real, allocatable, intent(in)  :: U_m(:,:,:), V_m(:,:,:), W_m(:,:,:), denom(:,:,:)
        type(timer_t),   intent(inout) :: flux_time, flux_corr_time, sum_time

        ! time-n state (RK3 base), previous-stage scratch, transport fluxes,
        ! buoyancy
        real, allocatable :: u_n(:,:,:), v_n(:,:,:), wr_n(:,:,:)
        real, allocatable :: u_p(:,:,:), v_p(:,:,:)
        real, allocatable :: fu(:,:,:), fv(:,:,:), fw(:,:,:)
        real, allocatable :: buoy(:,:,:)

        integer :: i, j, k, rk3_step
        real    :: t_fac

        ims = domain%ims; ime = domain%ime
        jms = domain%jms; jme = domain%jme
        kms = domain%kms; kme = domain%kme
        its = domain%its; ite = domain%ite
        jts = domain%jts; jte = domain%jte
        kts = domain%kts; kte = domain%kte

        allocate(u_n (ims:ime+1, kms:kme, jms:jme))
        allocate(u_p (ims:ime+1, kms:kme, jms:jme))
        allocate(v_n (ims:ime,   kms:kme, jms:jme+1))
        allocate(v_p (ims:ime,   kms:kme, jms:jme+1))
        allocate(wr_n(ims:ime,   kms:kme, jms:jme))
        allocate(fu  (ims:ime+1, kms:kme, jms:jme))
        allocate(fv  (ims:ime,   kms:kme, jms:jme+1))
        allocate(fw  (ims:ime,   kms:kme, jms:jme))
        allocate(buoy(ims:ime,   kms:kme, jms:jme))
        !$acc enter data create(u_n, u_p, v_n, v_p, wr_n, fu, fv, fw, buoy)

        call build_transport_and_sources(domain, options, dt, fu, fv, fw, buoy)

        associate(u_data  => domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d, &
                  v_data  => domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d, &
                  wr_data => domain%vars_3d(domain%var_indx(kVARS%w_real)%v)%data_3d)

        ! Save the time-n state (RK3 always restarts from it)
        !$acc kernels default(present)
        u_n  = u_data
        v_n  = v_data
        wr_n = wr_data
        !$acc end kernels

        do rk3_step = 1, 3
            select case(rk3_step)
            case(1); t_fac = 1.0/3.0
            case(2); t_fac = 0.5
            case(3); t_fac = 1.0
            end select

            ! Stage scratch: fluxes are evaluated from the previous stage's
            ! field while the update writes u_data in place.
            !$acc kernels default(present)
            u_p = u_data
            v_p = v_data
            !$acc end kernels

            call advect_u_stage(domain, options, t_fac, u_n, u_p, v_p, fu, fv, fw, u_data)
            call advect_v_stage(domain, options, t_fac, v_n, u_p, v_p, fu, fv, fw, v_data)

            ! w_real is cell-centred: the standard scalar kernel handles it
            ! (no flux corrector — velocities are not sign-preserving — and
            ! no constant-z diffusion).
            call adv_std_advect3d(wr_data, wr_n, U_m, V_m, W_m, denom, &
                    domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d, &
                    flux_time, flux_corr_time, sum_time, &
                    t_factor_in=t_fac, q_id_in=1, flux_corr_in=0, apply_cz_diff_in=.false.)

            ! Buoyancy on w_real, same stage weighting as the advective
            ! tendency (constant across stages: theta does not change
            ! inside the momentum RK3).
            !$acc parallel loop gang vector collapse(3) default(present)
            do j = jts, jte
                do k = kms, kme
                    do i = its, ite
                        wr_data(i,k,j) = wr_data(i,k,j) + t_fac * dt * buoy(i,k,j)
                    enddo
                enddo
            enddo

            call domain%halo%exch_var(domain%vars_3d(domain%var_indx(kVARS%u)%v), corners=.True.)
            call domain%halo%exch_var(domain%vars_3d(domain%var_indx(kVARS%v)%v), corners=.True.)
            call domain%halo%exch_var(domain%vars_3d(domain%var_indx(kVARS%w_real)%v), corners=.True.)
        enddo

        end associate

        call apply_smagorinsky(domain, options, dt)
        call apply_rayleigh_w(domain, dt)

        ! ------------------------------------------------------------
        ! Blending zone: taper the full dynamics increment to zero across
        ! the lateral-boundary relaxation ring (relax_filter: 1 at the
        ! boundary, decaying inward). In the ring the solution is
        ! forcing-dominated, in the interior dynamics-dominated — the
        ! standard specified/blended-zone treatment for limited-area
        ! prognostic cores. Without it, the shear line between the
        ! forcing-pinned ring and the free interior slowly amplifies w at
        ! the lid in the boundary corners (found in validation).
        ! ------------------------------------------------------------
        associate(u_data  => domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d, &
                  v_data  => domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d, &
                  wr_data => domain%vars_3d(domain%var_indx(kVARS%w_real)%v)%data_3d, &
                  rf      => domain%vars_3d(domain%var_indx(kVARS%relax_filter_3d)%v)%data_3d)
        !$acc parallel default(present)
        !$acc loop gang vector collapse(3) private(t_fac)
        do j = jts, jte
            do k = kms, kme
                do i = its, ite+1
                    t_fac = 0.5 * (rf(max(i-1,ims),k,j) + rf(min(i,ime),k,j))
                    u_data(i,k,j) = u_n(i,k,j) + (1.0 - t_fac) * (u_data(i,k,j) - u_n(i,k,j))
                enddo
            enddo
        enddo
        !$acc loop gang vector collapse(3) private(t_fac)
        do j = jts, jte+1
            do k = kms, kme
                do i = its, ite
                    t_fac = 0.5 * (rf(i,k,max(j-1,jms)) + rf(i,k,min(j,jme)))
                    v_data(i,k,j) = v_n(i,k,j) + (1.0 - t_fac) * (v_data(i,k,j) - v_n(i,k,j))
                enddo
            enddo
        enddo
        !$acc loop gang vector collapse(3)
        do j = jts, jte
            do k = kms, kme
                do i = its, ite
                    wr_data(i,k,j) = wr_n(i,k,j) + (1.0 - rf(i,k,j)) * (wr_data(i,k,j) - wr_n(i,k,j))
                enddo
            enddo
        enddo
        !$acc end parallel
        end associate

        ! Refresh halos after diffusion/damping/blending so the grid-w
        ! increment and the projection read a consistent field.
        call domain%halo%exch_var(domain%vars_3d(domain%var_indx(kVARS%u)%v), corners=.True.)
        call domain%halo%exch_var(domain%vars_3d(domain%var_indx(kVARS%v)%v), corners=.True.)
        call domain%halo%exch_var(domain%vars_3d(domain%var_indx(kVARS%w_real)%v), corners=.True.)

        ! ------------------------------------------------------------
        ! Map this step's PHYSICAL vertical-momentum increment onto the
        ! prognostic grid-w. The base grid-w carries exact discrete
        ! continuity from the previous projection closure; only the
        ! increment (advection + buoyancy + diffusion + damping of
        ! w_real, all small and smooth) is converted. Re-DERIVING grid-w
        ! from w_real wholesale (calc_idealized_wgrid) is NOT the inverse
        ! of the closure's w_real diagnostic: over steep cells the round
        ! trip re-injects O(w*slope) divergence every step, which the
        ! Poisson solve then "corrects" by persistently accelerating the
        ! horizontal wind — a secular instability found the hard way.
        ! Same dz-weighted staggering as calc_idealized_wgrid.
        ! ------------------------------------------------------------
        associate(w_grid => domain%vars_3d(domain%var_indx(kVARS%w)%v)%data_3d, &
                  wr     => domain%vars_3d(domain%var_indx(kVARS%w_real)%v)%data_3d, &
                  jaco_w => domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d, &
                  dz     => domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d)
        !$acc parallel loop gang vector collapse(3) default(present)
        do j = jms, jme
            do k = kms, kme-1
                do i = ims, ime
                    w_grid(i,k,j) = w_grid(i,k,j) + &
                        ( (wr(i,k,j)  - wr_n(i,k,j) ) * dz(i,k+1,j) + &
                          (wr(i,k+1,j)- wr_n(i,k+1,j)) * dz(i,k,j) ) / &
                        ( (dz(i,k,j) + dz(i,k+1,j)) * jaco_w(i,k,j) )
                enddo
            enddo
        enddo
        end associate
        call domain%halo%exch_var(domain%vars_3d(domain%var_indx(kVARS%w)%v), corners=.True.)

        !$acc exit data delete(u_n, u_p, v_n, v_p, wr_n, fu, fv, fw, buoy)
        deallocate(u_n, u_p, v_n, v_p, wr_n, fu, fv, fw, buoy)

    end subroutine rans_momentum_step


    !>--------------------------------------------------------
    !! Time-n contravariant Courant fluxes on the native staggered
    !! points (same weighting as adv_std_compute_wind), plus the
    !! buoyancy field.
    !!
    !!   fu = u * rho_u * J_u * dt/dx          (u-points)
    !!   fv = v * rho_v * J_v * dt/dx          (v-points)
    !!   fw = w_grid * rho_w * J_w * dt        (interfaces; /dz at use site)
    !!   buoy = g * (theta - theta_bar)/theta_bar   (mass points)
    !!
    !! rho stencils are clamped at the memory edges (only reached on
    !! global-boundary ranks, where the affected faces are BC-owned
    !! or first-order anyway).
    !!--------------------------------------------------------
    subroutine build_transport_and_sources(domain, options, dt, fu, fv, fw, buoy)
        implicit none
        type(domain_t),  intent(inout) :: domain
        type(options_t), intent(in)    :: options
        real,            intent(in)    :: dt
        real, intent(inout) :: fu  (ims:ime+1, kms:kme, jms:jme)
        real, intent(inout) :: fv  (ims:ime,   kms:kme, jms:jme+1)
        real, intent(inout) :: fw  (ims:ime,   kms:kme, jms:jme)
        real, intent(inout) :: buoy(ims:ime,   kms:kme, jms:jme)

        integer :: i, j, k
        real    :: rho_face, dtdx
        logical :: adv_den

        adv_den = options%adv%advect_density
        dtdx    = dt / domain%dx

        associate(u      => domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d, &
                  v      => domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d, &
                  w      => domain%vars_3d(domain%var_indx(kVARS%w)%v)%data_3d, &
                  rho    => domain%vars_3d(domain%var_indx(kVARS%density)%v)%data_3d, &
                  jaco_u => domain%vars_3d(domain%var_indx(kVARS%jacobian_u)%v)%data_3d, &
                  jaco_v => domain%vars_3d(domain%var_indx(kVARS%jacobian_v)%v)%data_3d, &
                  jaco_w => domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d, &
                  dz     => domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d, &
                  theta  => domain%vars_3d(domain%var_indx(kVARS%potential_temperature)%v)%data_3d)

        !$acc parallel loop gang vector collapse(3) default(present) private(rho_face)
        do j = jms, jme
            do k = kms, kme
                do i = ims, ime+1
                    if (adv_den) then
                        rho_face = 0.5 * (rho(max(i-1,ims),k,j) + rho(min(i,ime),k,j))
                    else
                        rho_face = 1.0
                    endif
                    fu(i,k,j) = u(i,k,j) * rho_face * jaco_u(i,k,j) * dtdx
                enddo
            enddo
        enddo

        !$acc parallel loop gang vector collapse(3) default(present) private(rho_face)
        do j = jms, jme+1
            do k = kms, kme
                do i = ims, ime
                    if (adv_den) then
                        rho_face = 0.5 * (rho(i,k,max(j-1,jms)) + rho(i,k,min(j,jme)))
                    else
                        rho_face = 1.0
                    endif
                    fv(i,k,j) = v(i,k,j) * rho_face * jaco_v(i,k,j) * dtdx
                enddo
            enddo
        enddo

        !$acc parallel loop gang vector collapse(3) default(present) private(rho_face)
        do j = jms, jme
            do k = kms, kme
                do i = ims, ime
                    if (.not. adv_den) then
                        rho_face = 1.0
                    else if (k < kme) then
                        rho_face = ( rho(i,k,j)*dz(i,k+1,j) + rho(i,k+1,j)*dz(i,k,j) ) / &
                                   ( dz(i,k,j) + dz(i,k+1,j) )
                    else
                        rho_face = rho(i,kme,j)
                    endif
                    fw(i,k,j) = w(i,k,j) * rho_face * jaco_w(i,k,j) * dt

                    buoy(i,k,j) = gravity * (theta(i,k,j) - adv_theta_ref(i,k,j)) &
                                          / adv_theta_ref(i,k,j)
                enddo
            enddo
        enddo

        end associate
    end subroutine build_transport_and_sources


    !>--------------------------------------------------------
    !! 3rd-order WRF upwind face flux (1st-order handled by the
    !! callers where the stencil is incomplete):
    !!   F = m*(7(qr+ql) - (qrr+qll))/12 - |m|*(3(qr-ql) - (qrr-qll))/12
    !! with ql/qr the values adjacent to the face and qll/qrr the
    !! next ones outward.
    !!--------------------------------------------------------
    pure function flux3rd(m, qll, ql, qr, qrr) result(f)
        !$acc routine seq
        real, intent(in) :: m, qll, ql, qr, qrr
        real :: f
        f = ( m * (7.0*(qr + ql) - (qrr + qll)) &
              - abs(m) * (3.0*(qr - ql) - (qrr - qll)) ) / 12.0
    end function flux3rd

    pure function flux1st(m, ql, qr) result(f)
        !$acc routine seq
        real, intent(in) :: m, ql, qr
        real :: f
        f = 0.5 * ( (m + abs(m))*ql + (m - abs(m))*qr )
    end function flux1st

    ! Centered and dissipative parts of the 3rd-order upwind face flux,
    ! exposed separately for the per-column vertical flux divergence:
    !   F(m) = m*cent - |m|*diss
    pure function cent3rd(qll, ql, qr, qrr) result(c)
        !$acc routine seq
        real, intent(in) :: qll, ql, qr, qrr
        real :: c
        c = (7.0*(qr + ql) - (qrr + qll)) / 12.0
    end function cent3rd

    pure function diss3rd(qll, ql, qr, qrr) result(dval)
        !$acc routine seq
        real, intent(in) :: qll, ql, qr, qrr
        real :: dval
        dval = (3.0*(qr - ql) - (qrr - qll)) / 12.0
    end function diss3rd


    !>--------------------------------------------------------
    !! One RK3 stage of u advection on the staggered u control
    !! volumes. Reads the previous-stage fields (u_p, v_p),
    !! restarts from the time-n state u_n, writes u_out.
    !!
    !! u-CV at u-point i: x-faces at mass centres i-1, i;
    !! y-faces at corners (i,j), (i,j+1); z-faces at u-column
    !! interfaces k-1, k.
    !!--------------------------------------------------------
    subroutine advect_u_stage(domain, options, t_fac, u_n, u_p, v_p, fu, fv, fw, u_out)
        implicit none
        type(domain_t),  intent(in) :: domain
        type(options_t), intent(in) :: options
        real, intent(in)    :: t_fac
        real, intent(in)    :: u_n(ims:ime+1, kms:kme, jms:jme)
        real, intent(in)    :: u_p(ims:ime+1, kms:kme, jms:jme)
        real, intent(in)    :: v_p(ims:ime,   kms:kme, jms:jme+1)
        real, intent(in)    :: fu (ims:ime+1, kms:kme, jms:jme)
        real, intent(in)    :: fv (ims:ime,   kms:kme, jms:jme+1)
        real, intent(in)    :: fw (ims:ime,   kms:kme, jms:jme)
        real, intent(inout) :: u_out(ims:ime+1, kms:kme, jms:jme)

        integer :: i, j, k, iu_s, iu_e, ju_s, ju_e, c1, c2
        real    :: m_e, m_w, m_n, m_s
        real    :: f_e, f_w, f_n, f_s
        real    :: c_t, d_t, c_b, d_b
        real    :: fwt1, fwt2, fwb1, fwb2, zdiv
        real    :: denom_u, rho_u, fdiv
        logical :: adv_den

        adv_den = options%adv%advect_density

        ! Global boundary u-faces are owned by the lateral BC relaxation.
        iu_s = its;   if (domain%west_boundary)  iu_s = its + 1
        iu_e = ite+1; if (domain%east_boundary)  iu_e = ite
        ju_s = jts;   if (domain%south_boundary) ju_s = jts + 1
        ju_e = jte;   if (domain%north_boundary) ju_e = jte - 1

        associate(rho    => domain%vars_3d(domain%var_indx(kVARS%density)%v)%data_3d, &
                  jaco_u => domain%vars_3d(domain%var_indx(kVARS%jacobian_u)%v)%data_3d, &
                  dz     => domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d)

        !$acc parallel loop gang vector collapse(3) default(present) &
        !$acc   private(m_e, m_w, m_n, m_s, f_e, f_w, f_n, f_s, c_t, d_t, c_b, d_b, &
        !$acc           c1, c2, fwt1, fwt2, fwb1, fwb2, zdiv, denom_u, rho_u, fdiv)
        do j = ju_s, ju_e
            do k = kms, kme
                do i = iu_s, iu_e

                    ! --- x: CV faces at mass centres i-1 (west) and i (east).
                    ! Face at centre c sits between u(c) and u(c+1).
                    m_e = 0.5 * (fu(i,k,j)   + fu(i+1,k,j))
                    m_w = 0.5 * (fu(i-1,k,j) + fu(i,k,j))
                    if (i+2 <= ime+1) then
                        f_e = flux3rd(m_e, u_p(i-1,k,j), u_p(i,k,j), u_p(i+1,k,j), u_p(i+2,k,j))
                    else
                        f_e = flux1st(m_e, u_p(i,k,j), u_p(i+1,k,j))
                    endif
                    if (i-2 >= ims) then
                        f_w = flux3rd(m_w, u_p(i-2,k,j), u_p(i-1,k,j), u_p(i,k,j), u_p(i+1,k,j))
                    else
                        f_w = flux1st(m_w, u_p(i-1,k,j), u_p(i,k,j))
                    endif

                    ! --- y: CV faces at corners (i,j) (south) and (i,j+1)
                    ! (north). Corner jc sits between u(.,jc-1) and u(.,jc).
                    m_n = 0.5 * (fv(i-1,k,j+1) + fv(i,k,j+1))
                    m_s = 0.5 * (fv(i-1,k,j)   + fv(i,k,j))
                    if (j+2 <= jme) then
                        f_n = flux3rd(m_n, u_p(i,k,j-1), u_p(i,k,j), u_p(i,k,j+1), u_p(i,k,j+2))
                    else
                        f_n = flux1st(m_n, u_p(i,k,j), u_p(i,k,j+1))
                    endif
                    if (j-2 >= jms) then
                        f_s = flux3rd(m_s, u_p(i,k,j-2), u_p(i,k,j-1), u_p(i,k,j), u_p(i,k,j+1))
                    else
                        f_s = flux1st(m_s, u_p(i,k,j-1), u_p(i,k,j))
                    endif

                    ! --- z: CV faces at u-column interfaces k (top) and k-1
                    ! (bottom). Interface kf sits between u(kf) and u(kf+1);
                    ! the surface flux (below kms) is zero, the lid face uses
                    ! the outflow-upwind form of flux3.
                    !
                    ! The reconstructed face value (f_t/f_b WITHOUT the mass
                    ! flux) is shared between the two columns the CV spans,
                    ! but each column's flux difference is divided by ITS OWN
                    ! layer thickness before averaging. This makes the CV sum
                    ! telescope to the average of the two columns' per-cell
                    ! mass divergences — which the projection closure has
                    ! zeroed — so constant flow is preserved exactly. Using a
                    ! CV-averaged dz instead leaves a residual source
                    ! ~ u*w*(d(dz)/dx)/dz that runs away over sharp peaks
                    ! where near-surface w is large.
                    ! Centered (c) and dissipative (d) parts of the shared
                    ! interface reconstructions of u at the CV's top (t) and
                    ! bottom (b) faces. Face flux for column c is
                    ! fw(c)*cent - |fw(c)|*diss.
                    if (k < kme) then
                        if (k >= kms+1 .and. k+2 <= kme) then
                            c_t = cent3rd(u_p(i,k-1,j), u_p(i,k,j), u_p(i,k+1,j), u_p(i,k+2,j))
                            d_t = diss3rd(u_p(i,k-1,j), u_p(i,k,j), u_p(i,k+1,j), u_p(i,k+2,j))
                        else
                            c_t = 0.5 * (u_p(i,k,j) + u_p(i,k+1,j))
                            d_t = 0.5 * (u_p(i,k+1,j) - u_p(i,k,j))
                        endif
                    else
                        ! lid face: outflow-upwind form (matches flux3's top slab)
                        c_t = u_p(i,kme,j)
                        d_t = 0.0
                    endif
                    if (k > kms) then
                        if (k-2 >= kms .and. k+1 <= kme) then
                            c_b = cent3rd(u_p(i,k-2,j), u_p(i,k-1,j), u_p(i,k,j), u_p(i,k+1,j))
                            d_b = diss3rd(u_p(i,k-2,j), u_p(i,k-1,j), u_p(i,k,j), u_p(i,k+1,j))
                        else
                            c_b = 0.5 * (u_p(i,k-1,j) + u_p(i,k,j))
                            d_b = 0.5 * (u_p(i,k,j) - u_p(i,k-1,j))
                        endif
                    else
                        c_b = 0.0
                        d_b = 0.0
                    endif

                    c1 = max(i-1,ims); c2 = min(i,ime)
                    fwt1 = fw(c1,k,j);  fwt2 = fw(c2,k,j)
                    if (k > kms) then
                        fwb1 = fw(c1,k-1,j); fwb2 = fw(c2,k-1,j)
                    else
                        fwb1 = 0.0; fwb2 = 0.0
                    endif
                    zdiv = 0.5 * ( (fwt1*c_t - abs(fwt1)*d_t - (fwb1*c_b - abs(fwb1)*d_b)) / dz(c1,k,j) &
                                 + (fwt2*c_t - abs(fwt2)*d_t - (fwb2*c_b - abs(fwb2)*d_b)) / dz(c2,k,j) )

                    if (adv_den) then
                        rho_u = 0.5 * (rho(c1,k,j) + rho(c2,k,j))
                    else
                        rho_u = 1.0
                    endif
                    denom_u = 1.0 / (rho_u * jaco_u(i,k,j))

                    fdiv = (f_e - f_w) + (f_n - f_s) + zdiv
                    u_out(i,k,j) = u_n(i,k,j) - t_fac * fdiv * denom_u
                enddo
            enddo
        enddo

        end associate
    end subroutine advect_u_stage


    !>--------------------------------------------------------
    !! One RK3 stage of v advection on the staggered v control
    !! volumes (mirror of advect_u_stage).
    !!--------------------------------------------------------
    subroutine advect_v_stage(domain, options, t_fac, v_n, u_p, v_p, fu, fv, fw, v_out)
        implicit none
        type(domain_t),  intent(in) :: domain
        type(options_t), intent(in) :: options
        real, intent(in)    :: t_fac
        real, intent(in)    :: v_n(ims:ime, kms:kme, jms:jme+1)
        real, intent(in)    :: u_p(ims:ime+1, kms:kme, jms:jme)
        real, intent(in)    :: v_p(ims:ime,   kms:kme, jms:jme+1)
        real, intent(in)    :: fu (ims:ime+1, kms:kme, jms:jme)
        real, intent(in)    :: fv (ims:ime,   kms:kme, jms:jme+1)
        real, intent(in)    :: fw (ims:ime,   kms:kme, jms:jme)
        real, intent(inout) :: v_out(ims:ime, kms:kme, jms:jme+1)

        integer :: i, j, k, jv_s, jv_e, iv_s, iv_e, c1, c2
        real    :: m_e, m_w, m_n, m_s
        real    :: f_e, f_w, f_n, f_s
        real    :: c_t, d_t, c_b, d_b
        real    :: fwt1, fwt2, fwb1, fwb2, zdiv
        real    :: denom_v, rho_v, fdiv
        logical :: adv_den

        adv_den = options%adv%advect_density

        jv_s = jts;   if (domain%south_boundary) jv_s = jts + 1
        jv_e = jte+1; if (domain%north_boundary) jv_e = jte
        iv_s = its;   if (domain%west_boundary)  iv_s = its + 1
        iv_e = ite;   if (domain%east_boundary)  iv_e = ite - 1

        associate(rho    => domain%vars_3d(domain%var_indx(kVARS%density)%v)%data_3d, &
                  jaco_v => domain%vars_3d(domain%var_indx(kVARS%jacobian_v)%v)%data_3d, &
                  dz     => domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d)

        !$acc parallel loop gang vector collapse(3) default(present) &
        !$acc   private(m_e, m_w, m_n, m_s, f_e, f_w, f_n, f_s, c_t, d_t, c_b, d_b, &
        !$acc           c1, c2, fwt1, fwt2, fwb1, fwb2, zdiv, denom_v, rho_v, fdiv)
        do j = jv_s, jv_e
            do k = kms, kme
                do i = iv_s, iv_e

                    ! --- y: CV faces at mass centres j-1 (south) and j
                    ! (north). Face at centre c sits between v(c) and v(c+1).
                    m_n = 0.5 * (fv(i,k,j)   + fv(i,k,j+1))
                    m_s = 0.5 * (fv(i,k,j-1) + fv(i,k,j))
                    if (j+2 <= jme+1) then
                        f_n = flux3rd(m_n, v_p(i,k,j-1), v_p(i,k,j), v_p(i,k,j+1), v_p(i,k,j+2))
                    else
                        f_n = flux1st(m_n, v_p(i,k,j), v_p(i,k,j+1))
                    endif
                    if (j-2 >= jms) then
                        f_s = flux3rd(m_s, v_p(i,k,j-2), v_p(i,k,j-1), v_p(i,k,j), v_p(i,k,j+1))
                    else
                        f_s = flux1st(m_s, v_p(i,k,j-1), v_p(i,k,j))
                    endif

                    ! --- x: CV faces at corners (i,j) (west) and (i+1,j)
                    ! (east). Corner ic sits between v(ic-1,.) and v(ic,.).
                    m_e = 0.5 * (fu(i+1,k,j-1) + fu(i+1,k,j))
                    m_w = 0.5 * (fu(i,k,j-1)   + fu(i,k,j))
                    if (i+2 <= ime) then
                        f_e = flux3rd(m_e, v_p(i-1,k,j), v_p(i,k,j), v_p(i+1,k,j), v_p(i+2,k,j))
                    else
                        f_e = flux1st(m_e, v_p(i,k,j), v_p(i+1,k,j))
                    endif
                    if (i-2 >= ims) then
                        f_w = flux3rd(m_w, v_p(i-2,k,j), v_p(i-1,k,j), v_p(i,k,j), v_p(i+1,k,j))
                    else
                        f_w = flux1st(m_w, v_p(i-1,k,j), v_p(i,k,j))
                    endif

                    ! --- z: per-column vertical flux divergence (see
                    ! advect_u_stage for the telescoping rationale).
                    if (k < kme) then
                        if (k >= kms+1 .and. k+2 <= kme) then
                            c_t = cent3rd(v_p(i,k-1,j), v_p(i,k,j), v_p(i,k+1,j), v_p(i,k+2,j))
                            d_t = diss3rd(v_p(i,k-1,j), v_p(i,k,j), v_p(i,k+1,j), v_p(i,k+2,j))
                        else
                            c_t = 0.5 * (v_p(i,k,j) + v_p(i,k+1,j))
                            d_t = 0.5 * (v_p(i,k+1,j) - v_p(i,k,j))
                        endif
                    else
                        c_t = v_p(i,kme,j)
                        d_t = 0.0
                    endif
                    if (k > kms) then
                        if (k-2 >= kms .and. k+1 <= kme) then
                            c_b = cent3rd(v_p(i,k-2,j), v_p(i,k-1,j), v_p(i,k,j), v_p(i,k+1,j))
                            d_b = diss3rd(v_p(i,k-2,j), v_p(i,k-1,j), v_p(i,k,j), v_p(i,k+1,j))
                        else
                            c_b = 0.5 * (v_p(i,k-1,j) + v_p(i,k,j))
                            d_b = 0.5 * (v_p(i,k,j) - v_p(i,k-1,j))
                        endif
                    else
                        c_b = 0.0
                        d_b = 0.0
                    endif

                    c1 = max(j-1,jms); c2 = min(j,jme)
                    fwt1 = fw(i,k,c1);  fwt2 = fw(i,k,c2)
                    if (k > kms) then
                        fwb1 = fw(i,k-1,c1); fwb2 = fw(i,k-1,c2)
                    else
                        fwb1 = 0.0; fwb2 = 0.0
                    endif
                    zdiv = 0.5 * ( (fwt1*c_t - abs(fwt1)*d_t - (fwb1*c_b - abs(fwb1)*d_b)) / dz(i,k,c1) &
                                 + (fwt2*c_t - abs(fwt2)*d_t - (fwb2*c_b - abs(fwb2)*d_b)) / dz(i,k,c2) )

                    if (adv_den) then
                        rho_v = 0.5 * (rho(i,k,c1) + rho(i,k,c2))
                    else
                        rho_v = 1.0
                    endif
                    denom_v = 1.0 / (rho_v * jaco_v(i,k,j))

                    fdiv = (f_e - f_w) + (f_n - f_s) + zdiv
                    v_out(i,k,j) = v_n(i,k,j) - t_fac * fdiv * denom_v

                enddo
            enddo
        enddo

        end associate
    end subroutine advect_v_stage


    !>--------------------------------------------------------
    !! Horizontal Smagorinsky deformation diffusion (along model
    !! levels) for u, v, w_real. Forward Euler with the eddy
    !! viscosity clipped to its explicit stability bound. The
    !! outermost global ring is skipped (BC-owned).
    !!--------------------------------------------------------
    subroutine apply_smagorinsky(domain, options, dt)
        implicit none
        type(domain_t),  intent(inout) :: domain
        type(options_t), intent(in)    :: options
        real,            intent(in)    :: dt

        real, allocatable :: Ksmag(:,:,:)
        real, allocatable :: u_o(:,:,:), v_o(:,:,:), wr_o(:,:,:)
        integer :: i, j, k, i_s, i_e, j_s, j_e
        real    :: dx, d11, d12, k_max, k_fac, kf, lap
        real    :: uc_jp, uc_jm, vc_ip, vc_im

        dx    = domain%dx
        k_max = RANS_SMAG_CLIP * dx * dx / dt
        ! K is used as K*dt/dx^2 * (4-point Laplacian sum) below
        k_fac = dt / (dx * dx)

        allocate(Ksmag(ims:ime,   kms:kme, jms:jme))
        allocate(u_o  (ims:ime+1, kms:kme, jms:jme))
        allocate(v_o  (ims:ime,   kms:kme, jms:jme+1))
        allocate(wr_o (ims:ime,   kms:kme, jms:jme))
        !$acc enter data create(Ksmag, u_o, v_o, wr_o)

        associate(u  => domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d, &
                  v  => domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d, &
                  wr => domain%vars_3d(domain%var_indx(kVARS%w_real)%v)%data_3d)

        ! Deformation magnitude at mass points (index-clamped at memory
        ! edges; only the diffusion of halo-adjacent faces sees the
        ! degraded estimate).
        !$acc parallel loop gang vector collapse(3) default(present) &
        !$acc   private(d11, d12, uc_jp, uc_jm, vc_ip, vc_im)
        do j = jms, jme
            do k = kms, kme
                do i = ims, ime
                    d11 = ( (u(i+1,k,j) - u(i,k,j)) - (v(i,k,j+1) - v(i,k,j)) ) / dx

                    uc_jp = 0.5 * (u(i,k,min(j+1,jme)) + u(i+1,k,min(j+1,jme)))
                    uc_jm = 0.5 * (u(i,k,max(j-1,jms)) + u(i+1,k,max(j-1,jms)))
                    vc_ip = 0.5 * (v(min(i+1,ime),k,j) + v(min(i+1,ime),k,j+1))
                    vc_im = 0.5 * (v(max(i-1,ims),k,j) + v(max(i-1,ims),k,j+1))
                    d12 = (uc_jp - uc_jm + vc_ip - vc_im) / (2.0 * dx)

                    Ksmag(i,k,j) = min( (RANS_SMAG_CS*dx)**2 * sqrt(d11*d11 + d12*d12), k_max )
                enddo
            enddo
        enddo

        ! Snapshot the fields: the Laplacian must read pre-update
        ! neighbours (in-place update would race on the GPU).
        !$acc kernels default(present)
        u_o  = u
        v_o  = v
        wr_o = wr
        !$acc end kernels

        ! u faces (skip global boundary ring)
        i_s = its;   if (domain%west_boundary)  i_s = its + 1
        i_e = ite+1; if (domain%east_boundary)  i_e = ite
        j_s = jts;   if (domain%south_boundary) j_s = jts + 1
        j_e = jte;   if (domain%north_boundary) j_e = jte - 1
        !$acc parallel loop gang vector collapse(3) default(present) private(kf, lap)
        do j = j_s, j_e
            do k = kms, kme
                do i = i_s, i_e
                    kf  = 0.5 * (Ksmag(max(i-1,ims),k,j) + Ksmag(min(i,ime),k,j))
                    lap = u_o(i+1,k,j) + u_o(i-1,k,j) + u_o(i,k,j+1) + u_o(i,k,j-1) - 4.0*u_o(i,k,j)
                    u(i,k,j) = u_o(i,k,j) + kf * k_fac * lap
                enddo
            enddo
        enddo

        ! v faces
        i_s = its;   if (domain%west_boundary)  i_s = its + 1
        i_e = ite;   if (domain%east_boundary)  i_e = ite - 1
        j_s = jts;   if (domain%south_boundary) j_s = jts + 1
        j_e = jte+1; if (domain%north_boundary) j_e = jte
        !$acc parallel loop gang vector collapse(3) default(present) private(kf, lap)
        do j = j_s, j_e
            do k = kms, kme
                do i = i_s, i_e
                    kf  = 0.5 * (Ksmag(i,k,max(j-1,jms)) + Ksmag(i,k,min(j,jme)))
                    lap = v_o(i+1,k,j) + v_o(i-1,k,j) + v_o(i,k,j+1) + v_o(i,k,j-1) - 4.0*v_o(i,k,j)
                    v(i,k,j) = v_o(i,k,j) + kf * k_fac * lap
                enddo
            enddo
        enddo

        ! w_real (mass points)
        i_s = its; if (domain%west_boundary)  i_s = its + 1
        i_e = ite; if (domain%east_boundary)  i_e = ite - 1
        j_s = jts; if (domain%south_boundary) j_s = jts + 1
        j_e = jte; if (domain%north_boundary) j_e = jte - 1
        !$acc parallel loop gang vector collapse(3) default(present) private(kf, lap)
        do j = j_s, j_e
            do k = kms, kme
                do i = i_s, i_e
                    kf  = Ksmag(i,k,j)
                    lap = wr_o(i+1,k,j) + wr_o(i-1,k,j) + wr_o(i,k,j+1) + wr_o(i,k,j-1) - 4.0*wr_o(i,k,j)
                    wr(i,k,j) = wr_o(i,k,j) + kf * k_fac * lap
                enddo
            enddo
        enddo

        end associate

        !$acc exit data delete(Ksmag, u_o, v_o, wr_o)
        deallocate(Ksmag, u_o, v_o, wr_o)
    end subroutine apply_smagorinsky


    !>--------------------------------------------------------
    !! Implicit Rayleigh damping of w_real over the top of the
    !! column:  w <- w / (1 + dt/tau(z)),
    !! tau(z)^-1 = tau_max^-1 * sin^2((pi/2)*(z-z_d)/(z_top-z_d))
    !! for z > z_d = frac * column depth above terrain.
    !!--------------------------------------------------------
    subroutine apply_rayleigh_w(domain, dt)
        implicit none
        type(domain_t), intent(inout) :: domain
        real,           intent(in)    :: dt

        integer :: i, j, k
        real    :: zagl, ztop, zd, s, inv_tau

        associate(wr => domain%vars_3d(domain%var_indx(kVARS%w_real)%v)%data_3d, &
                  z  => domain%vars_3d(domain%var_indx(kVARS%z)%v)%data_3d)

        !$acc parallel loop gang vector collapse(3) default(present) &
        !$acc   private(zagl, ztop, zd, s, inv_tau)
        do j = jts, jte
            do k = kms, kme
                do i = its, ite
                    ztop = z(i,kme,j) - z(i,kms,j)
                    zd   = RANS_RAYLEIGH_FRAC * ztop
                    zagl = z(i,k,j) - z(i,kms,j)
                    if (zagl > zd) then
                        s = sin( 0.5 * pi * (zagl - zd) / (ztop - zd) )
                        inv_tau = (s*s) / RANS_RAYLEIGH_TAU
                        wr(i,k,j) = wr(i,k,j) / (1.0 + dt * inv_tau)
                    endif
                enddo
            enddo
        enddo

        end associate
    end subroutine apply_rayleigh_w

end module wind_rans
