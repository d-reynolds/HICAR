!>------------------------------------------------------------
!! FLake bulk lake model — HICAR coupling
!!
!! FLake (Freshwater Lake model) is a two-layer self-similar bulk model that
!! predicts the time evolution of a lake's mixed-layer temperature, mixed-layer
!! depth, mean-water-column temperature, bottom temperature, thermocline shape
!! factor, and ice/snow/sediment cover with no vertical grid. Reference:
!!
!!   Mironov, D. V., 2008: Parameterization of lakes in numerical weather
!!     prediction. Description of a lake model. COSMO Tech. Rep. No. 11, DWD.
!!
!! This module is the HICAR-side interface around the FLake physics core
!! (vendored under LGPL into `src/physics/flake_core.F90`). Prognostic state
!! is held in 2D kVARS fields:
!!     flake_t_snow, flake_t_ice, flake_t_mnw, flake_t_wml, flake_t_bot,
!!     flake_t_b1,   flake_c_t,   flake_h_ice, flake_h_ml,  flake_h_b1
!! plus the reused kVARS%skin_temperature (T_sfc), kVARS%snow_height (h_snow),
!! kVARS%albedo (output diagnostic), and kVARS%xice (output diagnostic), with
!! the static external parameter kVARS%lake_depth (depth_w).
!!
!! Turbulent surface fluxes are computed bulk-style using the chs exchange
!! coefficient produced by sfclayrev (same path as water_simple) with the
!! freshwater 1.0 saturation factor (no 0.98 salinity scaling). Net longwave
!! at the surface is computed here as glw - emiss·σ·T_sfc^4 and combined with
!! SH+LH into the FLake Q_*_flk partition (snow vs ice vs water).
!!
!! Non-lake water cells (lakemask == 0) are still handled by water_simple in
!! lsm_driver.F90 before this routine is called — this module only touches
!! cells where lakemask == 1.
!!------------------------------------------------------------
module module_water_flake

    use data_structures
    use icar_constants,         only : kVARS, kWATER_FLAKE, kLC_WATER, STD_OUT_PE
    use options_interface,      only : options_t
    use domain_interface,       only : domain_t
    use mod_wrf_constants,      only : XLV
    use mod_atm_utilities,      only : sat_mr
    use flake_data_parameters,  only : ireals
    use flake_parameters,       only : tpl_T_f, h_Ice_min_flk, h_Snow_min_flk
    use flake_albedo_ref,       only : albedo_water_ref, albedo_whiteice_ref, &
                                       albedo_drysnow_ref
    use flake_paramoptic_ref,   only : opticpar_water_ref, opticpar_whiteice_ref, &
                                       opticpar_drysnow_ref
    use flake_core,             only : flake_state_t, flake_radflux, flake_driver
    implicit none
    private
    public :: flake_init, flake_step

    ! FLake external parameters — defaults; promote to namelist if/when needed.
    real, parameter :: FLAKE_DEPTH_BS_DEFAULT = 10.0     ! [m]  thermally active sediment layer depth
    real, parameter :: FLAKE_T_BS_DEFAULT     = 277.0    ! [K]  climatological sediment temperature (~4 °C)
    real, parameter :: FLAKE_FETCH_DEFAULT    = 1.0e3    ! [m]  typical wind fetch
    real, parameter :: FLAKE_TFRZ             = 273.15   ! [K]  freezing point of freshwater
    real, parameter :: FLAKE_MIN_DEPTH        = 1.0      ! [m]  hard floor for legitimately shallow lakes
    real, parameter :: FLAKE_DEFAULT_DEPTH    = 50.0     ! [m]  fallback when lake_depth input is 0/missing
                                                         !      (matches WRF's lakedepth_default; without it
                                                         !       lakes collapse to FLAKE_MIN_DEPTH and T_sfc
                                                         !       equilibrates to T_air within hours, killing SH)
    real, parameter :: FLAKE_CT_INIT          = 0.5      ! [-]  initial thermocline shape factor

    real, parameter :: PI_FLAKE       = 3.14159265
    real, parameter :: OMEGA_EARTH    = 7.2921e-5        ! [s-1] Earth rotation rate
    real, parameter :: SIGMA_SB       = 5.670374e-8      ! [W m-2 K-4] Stefan-Boltzmann constant

    integer :: ims, ime, jms, jme, kms, kme

contains

    !>------------------------------------------------------------
    !! Initialize FLake bulk state for every lake cell.
    !!
    !! Seeds the 10 prognostic kVARS plus the reused skin_temperature and
    !! snow_height to physically reasonable starting values based on the
    !! prescribed water-surface temperature (kVARS%sst, set from the init_sst
    !! namelist or the sst_var forcing field) and the lake depth
    !! (kVARS%lake_depth). Skipped on restart — state is read back from the
    !! restart file via the kVARS restart machinery.
    !!
    !! Multi-nest safety: fully (re)initializes on every nest-context switch.
    !! No save-state sentinels here (HICAR multi-nest leak pattern).
    !!------------------------------------------------------------
    subroutine flake_init(domain, options, restart, context_change)
        implicit none
        type(domain_t),  intent(inout) :: domain
        type(options_t), intent(in)    :: options
        logical,         intent(in)    :: restart
        logical,         intent(in)    :: context_change

        integer :: i, j, its, ite, jts, jte
        real    :: depth_w, t_init

        ! On restart, prognostic state comes from the restart file - do not seed.
        if (restart) return

        its = domain%its ; ite = domain%ite
        jts = domain%jts ; jte = domain%jte
        ims = domain%ims ; ime = domain%ime
        kms = domain%kms ; kme = domain%kme
        jms = domain%jms ; jme = domain%jme
        ! lakemask is populated by lsm_init's shared `kWATER_LAKE .OR. kWATER_FLAKE`
        ! block (lsm_driver.F90:721-750), which handles both ISLAKE>0 (vegtype==ISLAKE)
        ! and ISLAKE==-1 (vegtype==ISWATER .AND. terrain>=1m) cases. We rely on that
        ! here rather than maintaining a duplicate, less-complete derivation.

        ! Seed prognostic state per lake cell.
        associate(lakemask  => domain%vars_2d(domain%var_indx(kVARS%lakemask)%v)%data_2d,           &
                  lakedepth => domain%vars_2d(domain%var_indx(kVARS%lake_depth)%v)%data_2d,         &
                  sst       => domain%vars_2d(domain%var_indx(kVARS%sst)%v)%data_2d,                &
                  tskin     => domain%vars_2d(domain%var_indx(kVARS%skin_temperature)%v)%data_2d,   &
                  hsnow     => domain%vars_2d(domain%var_indx(kVARS%snow_height)%v)%data_2d,        &
                  t_snow    => domain%vars_2d(domain%var_indx(kVARS%flake_t_snow)%v)%data_2d,       &
                  t_ice     => domain%vars_2d(domain%var_indx(kVARS%flake_t_ice)%v)%data_2d,        &
                  t_mnw     => domain%vars_2d(domain%var_indx(kVARS%flake_t_mnw)%v)%data_2d,        &
                  t_wml     => domain%vars_2d(domain%var_indx(kVARS%flake_t_wml)%v)%data_2d,        &
                  t_bot     => domain%vars_2d(domain%var_indx(kVARS%flake_t_bot)%v)%data_2d,        &
                  t_b1      => domain%vars_2d(domain%var_indx(kVARS%flake_t_b1)%v)%data_2d,         &
                  c_t       => domain%vars_2d(domain%var_indx(kVARS%flake_c_t)%v)%data_2d,          &
                  h_ice     => domain%vars_2d(domain%var_indx(kVARS%flake_h_ice)%v)%data_2d,        &
                  h_ml      => domain%vars_2d(domain%var_indx(kVARS%flake_h_ml)%v)%data_2d,         &
                  h_b1      => domain%vars_2d(domain%var_indx(kVARS%flake_h_b1)%v)%data_2d,         &
                  xice      => domain%vars_2d(domain%var_indx(kVARS%xice)%v)%data_2d)

            !$acc parallel loop gang vector collapse(2)                                              &
            !$acc&  present(lakemask, lakedepth, sst, tskin, hsnow,                                  &
            !$acc&          t_snow, t_ice, t_mnw, t_wml, t_bot, t_b1,                                &
            !$acc&          c_t, h_ice, h_ml, h_b1, xice)
            do j = jts, jte
                do i = its, ite
                    if (lakemask(i,j) > 0.5) then
                        ! Persist the effective depth back to kVARS%lake_depth so
                        ! the output field shows what FLake actually uses and so
                        ! flake_step doesn't have to re-apply the fallback. Cells
                        ! with missing input (lake_depth <= 0) get the 50 m
                        ! default; a 1 m lake would equilibrate with the
                        ! atmosphere within hours and collapse SH to zero.
                        if (lakedepth(i,j) <= 0.0) lakedepth(i,j) = FLAKE_DEFAULT_DEPTH
                        depth_w = max(lakedepth(i,j), FLAKE_MIN_DEPTH)
                        t_init  = max(sst(i,j), FLAKE_TFRZ)

                        ! Bulk water state — well-mixed start.
                        t_wml(i,j) = t_init
                        t_mnw(i,j) = t_init
                        t_bot(i,j) = t_init
                        h_ml (i,j) = depth_w
                        c_t  (i,j) = FLAKE_CT_INIT

                        ! Ice / snow on ice — none at init. Safe seed temps for
                        ! FLake to overwrite once ice nucleates.
                        h_ice (i,j) = 0.0
                        hsnow (i,j) = 0.0
                        t_ice (i,j) = t_init
                        t_snow(i,j) = t_init

                        ! Sediment — initialize toward climatological boundary.
                        t_b1(i,j) = FLAKE_T_BS_DEFAULT
                        h_b1(i,j) = FLAKE_DEPTH_BS_DEFAULT

                        ! Atmosphere-facing diagnostics.
                        tskin(i,j) = t_init
                        xice (i,j) = 0.0
                    endif
                end do
            end do
        end associate

    end subroutine flake_init


    !>------------------------------------------------------------
    !! Advance FLake by one LSM update interval (dt).
    !!
    !! For each lake cell (lakemask == 1):
    !!   1. Assemble atmospheric forcing from the lowest model level plus
    !!      surface radiation and precipitation.
    !!   2. Compute bulk turbulent fluxes from sfclayrev's chs coefficient
    !!      and FLake's predicted skin temperature (water_simple style,
    !!      freshwater 1.0 saturation factor).
    !!   3. Partition the net surface heat flux into FLake's Q_w / Q_ice /
    !!      Q_snow channels per the upstream FLake convention.
    !!   4. Step FLake (flake_radflux + flake_driver) on the column.
    !!   5. Write atmosphere-side outputs: skin_temperature, sensible_heat,
    !!      latent_heat, qfx, ground_heat_flux, albedo, xice.
    !!------------------------------------------------------------
    subroutine flake_step(domain, options, dt, current_precipitation, windspd)
        implicit none
        type(domain_t),  intent(inout) :: domain
        type(options_t), intent(in)    :: options
        real,            intent(in)    :: dt
        real, dimension(ims:ime,jms:jme), intent(in) :: current_precipitation  ! [mm] accumulated over dt
        real, dimension(ims:ime,jms:jme), intent(in) :: windspd                ! [m/s] sqrt(u_10m^2 + v_10m^2)

        integer :: i, j, its, ite, jts, jte

        its = domain%its ; ite = domain%ite
        jts = domain%jts ; jte = domain%jte

        associate(lakemask  => domain%vars_2d(domain%var_indx(kVARS%lakemask)%v)%data_2d,             &
                  lakedepth => domain%vars_2d(domain%var_indx(kVARS%lake_depth)%v)%data_2d,           &
                  tskin     => domain%vars_2d(domain%var_indx(kVARS%skin_temperature)%v)%data_2d,     &
                  hsnow     => domain%vars_2d(domain%var_indx(kVARS%snow_height)%v)%data_2d,          &
                  swdown    => domain%vars_2d(domain%var_indx(kVARS%shortwave)%v)%data_2d,            &
                  glw       => domain%vars_2d(domain%var_indx(kVARS%longwave)%v)%data_2d,             &
                  psfc      => domain%vars_2d(domain%var_indx(kVARS%surface_pressure)%v)%data_2d,     &
                  chs       => domain%vars_2d(domain%var_indx(kVARS%chs)%v)%data_2d,                  &
                  sh_out    => domain%vars_2d(domain%var_indx(kVARS%sensible_heat)%v)%data_2d,        &
                  lh_out    => domain%vars_2d(domain%var_indx(kVARS%latent_heat)%v)%data_2d,          &
                  qfx_out   => domain%vars_2d(domain%var_indx(kVARS%qfx)%v)%data_2d,                  &
                  grdflx    => domain%vars_2d(domain%var_indx(kVARS%ground_heat_flux)%v)%data_2d,     &
                  albedo    => domain%vars_2d(domain%var_indx(kVARS%albedo)%v)%data_2d,               &
                  xice      => domain%vars_2d(domain%var_indx(kVARS%xice)%v)%data_2d,                 &
                  emiss     => domain%vars_2d(domain%var_indx(kVARS%land_emissivity)%v)%data_2d,      &
                  latitude  => domain%vars_2d(domain%var_indx(kVARS%latitude)%v)%data_2d,             &
                  T3d       => domain%vars_3d(domain%var_indx(kVARS%temperature)%v)%data_3d,          &
                  q3d       => domain%vars_3d(domain%var_indx(kVARS%water_vapor)%v)%data_3d,          &
                  t_snow    => domain%vars_2d(domain%var_indx(kVARS%flake_t_snow)%v)%data_2d,         &
                  t_ice     => domain%vars_2d(domain%var_indx(kVARS%flake_t_ice)%v)%data_2d,          &
                  t_mnw     => domain%vars_2d(domain%var_indx(kVARS%flake_t_mnw)%v)%data_2d,          &
                  t_wml     => domain%vars_2d(domain%var_indx(kVARS%flake_t_wml)%v)%data_2d,          &
                  t_bot     => domain%vars_2d(domain%var_indx(kVARS%flake_t_bot)%v)%data_2d,          &
                  t_b1      => domain%vars_2d(domain%var_indx(kVARS%flake_t_b1)%v)%data_2d,           &
                  c_t       => domain%vars_2d(domain%var_indx(kVARS%flake_c_t)%v)%data_2d,            &
                  h_ice     => domain%vars_2d(domain%var_indx(kVARS%flake_h_ice)%v)%data_2d,          &
                  h_ml      => domain%vars_2d(domain%var_indx(kVARS%flake_h_ml)%v)%data_2d,           &
                  h_b1      => domain%vars_2d(domain%var_indx(kVARS%flake_h_b1)%v)%data_2d)

            !$acc parallel loop gang vector collapse(2)                                                      &
            !$acc&  present(lakemask, lakedepth, tskin, hsnow, swdown, glw, psfc, chs,                       &
            !$acc&          sh_out, lh_out, qfx_out, grdflx, albedo, xice, emiss, latitude,                  &
            !$acc&          T3d, q3d, windspd, current_precipitation,                                        &
            !$acc&          t_snow, t_ice, t_mnw, t_wml, t_bot, t_b1, c_t, h_ice, h_ml, h_b1)
            do j = jts, jte
                do i = its, ite
                    if (lakemask(i,j) > 0.5) then
                        ! lake_depth was populated by flake_init (0/missing -> 50 m),
                        ! so this just floors at 1 m for genuinely shallow lakes.
                        call flake_cell_step(                                              &
                            ! geometry / external params
                            max(lakedepth(i,j), FLAKE_MIN_DEPTH),                          &
                            FLAKE_DEPTH_BS_DEFAULT, FLAKE_T_BS_DEFAULT,                    &
                            FLAKE_FETCH_DEFAULT, latitude(i,j),                            &
                            ! atmospheric forcing
                            T3d(i,kms,j), q3d(i,kms,j), windspd(i,j), psfc(i,j),           &
                            swdown(i,j), glw(i,j), emiss(i,j),                             &
                            current_precipitation(i,j), dt,                                &
                            ! exchange coefficient
                            chs(i,j),                                                      &
                            ! prognostic state (inout)
                            t_snow(i,j), t_ice(i,j),  t_mnw(i,j), t_wml(i,j), t_bot(i,j),  &
                            t_b1(i,j),   c_t(i,j),    hsnow(i,j), h_ice(i,j), h_ml(i,j),   &
                            h_b1(i,j),   tskin(i,j),                                       &
                            ! diagnostics out
                            sh_out(i,j), lh_out(i,j), qfx_out(i,j), grdflx(i,j),           &
                            xice(i,j),   albedo(i,j))
                    endif
                end do
            end do
        end associate

    end subroutine flake_step


    !>------------------------------------------------------------
    !! Per-cell FLake update.
    !!
    !! Builds a `flake_state_t` from the input scalars, sets the surface
    !! flux partition, runs flake_radflux + flake_driver, writes the new
    !! state back to the args. Device-callable (`acc routine seq`) so the
    !! calling kernel is a flat gang-vector loop over the tile.
    !!
    !! Sign convention: per FLake's upstream interface, Q_*_flk is positive
    !! when heat leaves the surface (cooling). SH+LH from water_simple are
    !! computed as positive when heat leaves the water → matches directly.
    !! Net LW = glw - emiss·σ·T_sfc^4 (positive INTO surface) is subtracted
    !! from the assembled Q_w, then reassigned to Q_ice or Q_snow based on
    !! the current surface state (water / ice / ice+snow).
    !!------------------------------------------------------------
    subroutine flake_cell_step(                                                  &
            depth_w, depth_bs, T_bs, fetch, lat_deg,                             &
            T_air, q_air, U_a, P_a, swdown, glw, emiss_in,                       &
            precip_acc, dt, chs_val,                                             &
            T_snow, T_ice, T_mnw, T_wML, T_bot,                                  &
            T_B1,   C_T,   h_snow, h_ice, h_ML,                                  &
            H_B1,   T_sfc,                                                       &
            SH_out, LH_out, qfx_out, grdflx_out, xice_out, albedo_out)
        !$acc routine seq
        implicit none
        real, intent(in)    :: depth_w, depth_bs, T_bs, fetch, lat_deg
        real, intent(in)    :: T_air, q_air, U_a, P_a, swdown, glw, emiss_in
        real, intent(in)    :: precip_acc, dt, chs_val
        real, intent(inout) :: T_snow, T_ice, T_mnw, T_wML, T_bot
        real, intent(inout) :: T_B1,   C_T,   h_snow, h_ice, h_ML
        real, intent(inout) :: H_B1,   T_sfc
        real, intent(out)   :: SH_out, LH_out, qfx_out, grdflx_out, xice_out, albedo_out

        type(flake_state_t)   :: s
        real(KIND=ireals)     :: qv_s, sh, lh, evap, lw_up, lw_net
        real(KIND=ireals)     :: dMsnowdt, fcor, T_sfc_p_loc, T_sfc_n_loc
        real(KIND=ireals)     :: q_surface_net, extincoef_water_typ

        ! ----- 1. Coriolis & snowfall partition -----
        fcor = 2.0_ireals * real(OMEGA_EARTH, ireals) * &
               sin(real(lat_deg, ireals) * real(PI_FLAKE, ireals) / 180.0_ireals)
        if (T_air < FLAKE_TFRZ) then
            dMsnowdt = real(precip_acc, ireals) / max(real(dt, ireals), 1.0_ireals)
        else
            dMsnowdt = 0.0_ireals
        endif

        ! ----- 2. Bulk turbulent fluxes (water_simple style, freshwater factor 1.0) -----
        qv_s = real(sat_mr(T_sfc, P_a), ireals)
        ! Floor only - do NOT clamp qv_s to q_air; the qv_s > q_air gradient
        ! drives evaporation. (Mirrors the water_simple fix.)
        qv_s = max(qv_s, 1.0e-6_ireals)
        ! Widen T_sfc and T_air to ireals BEFORE subtracting. With default-real
        ! (single) HICAR kinds the subtraction at ~280 K loses everything below
        ! ~3e-5 K, so once flake_driver brings T_sfc within that band of T_air
        ! the single-precision difference collapses to exactly 0 and SH pins
        ! to 0 even though the underlying values still differ in double.
        sh   = real(chs_val, ireals) * real(U_a, ireals) * &
               (real(T_sfc, ireals) - real(T_air, ireals))     ! W m-2 (positive = heat OUT of water)
        evap = real(chs_val, ireals) * real(U_a, ireals) * &
               (qv_s - real(q_air, ireals))                    ! kg m-2 s-1
        lh   = evap * real(XLV, ireals)

        ! ----- 3. Net longwave (downward positive into surface) -----
        lw_up  = real(emiss_in, ireals) * real(SIGMA_SB, ireals) * real(T_sfc, ireals)**4
        lw_net = real(glw, ireals) - lw_up

        ! ----- 4. Load previous-step state into s -----
        s%T_mnw_p_flk  = real(T_mnw,  ireals)
        s%T_snow_p_flk = real(T_snow, ireals)
        s%T_ice_p_flk  = real(T_ice,  ireals)
        s%T_wML_p_flk  = real(T_wML,  ireals)
        s%T_bot_p_flk  = real(T_bot,  ireals)
        s%T_B1_p_flk   = real(T_B1,   ireals)
        s%h_snow_p_flk = real(h_snow, ireals)
        s%h_ice_p_flk  = real(h_ice,  ireals)
        s%h_ML_p_flk   = real(h_ML,   ireals)
        s%H_B1_p_flk   = real(H_B1,   ireals)
        s%C_T_p_flk    = real(C_T,    ireals)
        s%dMsnowdt_flk = dMsnowdt
        s%I_atm_flk    = real(swdown, ireals)    ! incident SW (no albedo - flake_radflux applies)

        ! ----- 5. Surface heat flux partition (FLake sign: positive = heat OUT) -----
        ! Assemble the net surface heat-loss term as if the surface were water,
        ! then reassign to ice or snow channel based on actual surface state.
        q_surface_net = sh + lh - lw_net    ! lw_net was positive-into-surface; this is positive-out
        if (s%h_ice_p_flk >= h_Ice_min_flk .and. s%h_snow_p_flk >= h_Snow_min_flk) then
            s%Q_snow_flk = q_surface_net
            s%Q_ice_flk  = 0.0_ireals
            s%Q_w_flk    = 0.0_ireals      ! computed by flake_driver from ice-water balance
        elseif (s%h_ice_p_flk >= h_Ice_min_flk) then
            s%Q_snow_flk = 0.0_ireals
            s%Q_ice_flk  = q_surface_net
            s%Q_w_flk    = 0.0_ireals      ! computed by flake_driver
        else
            s%Q_snow_flk = 0.0_ireals
            s%Q_ice_flk  = 0.0_ireals
            s%Q_w_flk    = q_surface_net   ! direct air-water flux
        endif

        ! ----- 6. FLake radiation extinction through snow/ice/water -----
        call flake_radflux(s, real(depth_w, ireals),                                 &
                           real(albedo_water_ref,   ireals),                         &
                           real(albedo_whiteice_ref, ireals),                        &
                           real(albedo_drysnow_ref,  ireals),                        &
                           opticpar_water_ref, opticpar_whiteice_ref,                &
                           opticpar_drysnow_ref)

        ! ----- 7. FLake thermodynamic driver -----
        T_sfc_p_loc         = real(T_sfc, ireals)
        extincoef_water_typ = opticpar_water_ref%extincoef_optic(1)
        call flake_driver(s, real(depth_w, ireals), real(depth_bs, ireals),          &
                          real(T_bs, ireals), fcor, extincoef_water_typ,             &
                          real(dt, ireals), T_sfc_p_loc, T_sfc_n_loc)

        ! ----- 8. Read updated state back to args -----
        T_mnw  = real(s%T_mnw_n_flk)
        T_snow = real(s%T_snow_n_flk)
        T_ice  = real(s%T_ice_n_flk)
        T_wML  = real(s%T_wML_n_flk)
        T_bot  = real(s%T_bot_n_flk)
        T_B1   = real(s%T_B1_n_flk)
        h_snow = real(s%h_snow_n_flk)
        h_ice  = real(s%h_ice_n_flk)
        h_ML   = real(s%h_ML_n_flk)
        H_B1   = real(s%H_B1_n_flk)
        C_T    = real(s%C_T_n_flk)
        T_sfc  = real(T_sfc_n_loc)

        ! ----- 9. Atmosphere-side outputs -----
        SH_out     = real(sh)
        LH_out     = real(lh)
        qfx_out    = real(evap)
        grdflx_out = real(s%Q_bot_flk)       ! sediment-to-water heat flux (W m-2)
        if (s%h_ice_n_flk > 0.0_ireals) then
            xice_out = 1.0
        else
            xice_out = 0.0
        endif
        ! Diagnostic surface albedo based on final-state surface.
        if (s%h_ice_n_flk >= h_Ice_min_flk .and. s%h_snow_n_flk >= h_Snow_min_flk) then
            albedo_out = real(albedo_drysnow_ref)
        elseif (s%h_ice_n_flk >= h_Ice_min_flk) then
            albedo_out = real(albedo_whiteice_ref)
        else
            albedo_out = real(albedo_water_ref)
        endif

    end subroutine flake_cell_step

end module module_water_flake
