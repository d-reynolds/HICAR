# HICAR prognostic RANS wind solver — formulation and discretization

Status: design + implementation reference for `wind = "RANS"` (`kRANS_WINDS`).
Code: `src/physics/wind_rans.F90` (momentum dynamics),
`src/physics/wind.F90::rans_project` (projection cascade),
`src/physics/advection_driver.F90::advect` (call sequence).

This document is the contract between the math and the code. If the code and
this document disagree, one of them is wrong — fix whichever it is.

## 1. Motivation and design principles

The standard HICAR wind solvers are *diagnostic*: every forcing interval they
take the interpolated forcing winds and minimally adjust them to satisfy
`∇·(ρ u) = 0` (variational solver) — the wind field carries no memory and no
momentum physics. The RANS solver makes `u, v, w` *prognostic*: momentum is
advected, buoyancy acts, PBL turbulence drags, and an anelastic pressure
projection enforces mass continuity every time step.

Design principles, learned partly from the failure of the first attempt:

1. **Never form a derivative of a full (unperturbed) thermodynamic field.**
   Single precision cannot represent the dynamically active O(1 Pa) signal on
   a 10⁵ Pa background. Pressure exists in this solver *only* as the
   projection's Lagrange multiplier; buoyancy uses θ′ = θ − θ̄(z) with the
   existing `adv_theta_ref` reference profile.
2. **Momentum advection uses the same terrain-following contravariant mass
   fluxes as scalar advection** (`U_m, V_m, W_m` weighting: ρ·J·Δt/Δx). The
   vertical flux divergence is `Δ(flux_z)/Δz · 1/(ρJ)` exactly as in
   `sum_kernel` — there is no bare `w_real/dz` anywhere, which is what blew up
   over steep terrain in the first attempt.
3. **Close continuity with the proven diagnostic cascade**: Poisson solve →
   correct (u, v) → recompute grid-w by column integration (`calc_w`) →
   diagnose `w_real` kinematically (`calc_w_real`). The post-projection field
   is divergence-free *to machine precision per column* regardless of the
   BiCGStab tolerance, so no divergence residual is carried into the next
   step.
4. **Reuse, don't fork.** The Poisson backend is `wind_iterative` with
   α ≡ 1 — the same BiCGStab + Block-Jacobi machinery as the diagnostic
   solver, unchanged.

## 2. Governing equations

Anelastic momentum equations in advective form (equivalently flux form, since
the transport field is divergence-free — see §4.1):

    ∂u/∂t = −(u·∇)u            + f(v − v_e)   − ∂π/∂x + F_u
    ∂v/∂t = −(u·∇)v            + f(u_e − u)   − ∂π/∂y + F_v
    ∂w/∂t = −(u·∇)w + g·θ′/θ̄                 − ∂π/∂z + F_w

subject to

    ∇·(ρ u) = 0

where `w` here is the *physical* vertical velocity (`kVARS%w_real`,
cell-centred), π is the pressure perturbation potential (never stored; it is
the projection multiplier λ up to scaling), θ′ = θ − θ̄ with θ̄ the
horizontally uniform reference profile `adv_theta_ref(i,k,j) = θ̄(z(i,k,j))`
already constructed by `adv_std_init_theta_ref`, and F the turbulent stress
divergence:

* vertical: YSU PBL momentum tendencies (`pbl_driver.F90`, applied to the
  staggered u/v grids each physics step under `kRANS_WINDS`) — this is the
  "RANS closure" in the vertical, including surface drag;
* horizontal: 2-D Smagorinsky deformation diffusion (§6).

### Geostrophic forcing (Coriolis + synoptic PGF in one term)

The Coriolis term is the EULAG ambient-state form **f·k×(u_e − u)** with
f = 2Ω·sin(lat) and u_e the *balanced forcing target* — the same
grid-relative, projection-balanced, time-interpolated forcing state the
Davies relaxation nudges toward (§7; `forcing_hi` u/v, advanced each step
by `apply_forcing`). The ambient state is assumed to carry the synoptic
pressure gradient (geostrophic/thermal-wind balance), so this single term
supplies both the Coriolis force and the large-scale PGF that the
anelastic core cannot otherwise represent: π is only the projection
multiplier, so internally *generated* meso-scale pressure gradients
(thermal circulations, wave-induced) are produced by the projection, but
an externally *imposed* synoptic gradient has no carrier without this
term. Only deviations from the balanced state feel rotation; a deviation
executes a bounded inertial oscillation (period 2π/f) instead of secular
drift, and combined with the YSU drag this supports an Ekman-balanced
PBL wind profile. Two degenerate alternatives are recorded as wrong:
full Coriolis f·k×u without a PGF rotates the *mean* flow inertially
(worse than omitting the term), and a forcing-diagnosed PGF without
Coriolis accelerates the flow with nothing to balance it.
Implementation: `wind_rans.F90::apply_geostrophic_forcing`, applied per
RK3 stage from the previous stage's winds, BC-owned global faces
skipped, increment tapered across the relaxation ring by the §7
blending. Precedent: EULAG's f×(u − u_e) about geostrophically balanced
environmental states (Prusa, Smolarkiewicz & Wyszogrodzki 2008);
WRF/ICON need no such term because their compressible cores prognose
the full mass field and regenerate the synoptic PGF from the lateral
boundary conditions.

The thermodynamic equation is untouched: θ is advected by the existing
scalar RK3, including the θ′-split and the `w_real·dθ̄/dz` background term
(`adv_std_apply_ref_vert`), which together with the buoyancy term supports
gravity waves with the correct N².

## 3. Grid and notation

Arakawa-C, terrain-following:

* `u(i,k,j)` on x-faces, `ims:ime+1` — face i sits between mass cells i−1, i
* `v(i,k,j)` on y-faces, `jms:jme+1`
* `w(i,k,j)` = **grid-relative** w on z-interfaces (interface k = top of cell
  k); diagnostic only in this solver (derived from continuity)
* `w_real(i,k,j)` = physical w at **mass centres**; prognostic
* mass fluxes from `adv_std_compute_wind` (units: Courant-like, Δt and 1/Δx
  baked in):
  * `U_m(i,k,j) = u·ρ_u·J_u·Δt/Δx` at u-points
  * `V_m(i,k,j) = v·ρ_v·J_v·Δt/Δx` at v-points
  * `W_m(i,k,j) = w·ρ_w·J_w·Δt` at interface k (divided by Δz at use site)
* cell update convention (`sum_kernel`):
  `q ← q_old − [ΔF_x + ΔF_y + ΔF_z/Δz]·denom`, `denom = 1/(ρJ)`.

The momentum solver computes its own staggered-CV transport fluxes directly
from `u, v, w, ρ, J_*` over the halo (the `U_m` arrays are tile-bounded and
do not extend far enough for staggered control volumes at tile edges), using
*identical* weighting formulas, so momentum and scalar transport are
consistent to interpolation order.

### Map-scale factors

The nominal `dx` is a grid distance; on a projected grid the true ground
distance is `dx/m` with m the map-scale factor (`use_map_factors`
namelist, default on; host-model machinery, not RANS-specific). m is
computed once at init from the hi-res lat/lon using local WGS84 radii of
curvature in double precision (`domain_obj.F90::init_map_factors`) — a
spherical-earth formula carries an up-to-~0.3% latitude bias, larger than
the distortion of a well-centred projection. Direct two-point factors on
the natural stagger (m_x at u-points, m_y at v-points), 4-point averages
for the transverse factors, `m_x·m_y` as the cell-area factor at mass
points. Degenerate lat/lon (idealized grids) falls back to m ≡ 1 with a
warning; with m ≡ 1 every use site multiplies by exactly 1.0, so the off
state is bit-identical to the pre-map-factor code.

The finite-volume placement keeps `sum_kernel`/`flux3` untouched: face
fluxes are divided by their transverse factor (true face length = dx/m —
`U_m /= m_y|_u`, `V_m /= m_x|_v`), the vertical flux by the cell-area
factor (`W_m /= m_x·m_y`), and the cell denominator multiplied by it
(`denom = m_x·m_y/(ρJ)`); the vertical-flux factors cancel per column
(m is column-constant), recovering

    dq/dt = −(m_x·m_y/ρJ)[ Δ(ρJu q/m_y)/Δx + Δ(ρJv q/m_x)/Δy ]
            −(1/ρJ) Δ(ρJ_w w q)/Δz.

`calc_divergence` carries the same form (faces ÷ transverse m, cell sum ×
m_x·m_y, vertical unchanged), and the correction operator G multiplies
its whole horizontal bracket — λ-gradient *and* dzdx cross term, since
the true slope is m·(grid slope) — by m_x|_u (resp. m_y|_v). **The probed
operator A = D∘G absorbs all of this automatically on the next probe**;
no stencil work. CFL: dt is divided by max(m) (conservative one-liner in
`update_dt`).

The snow-drift fine-mesh advection (`adv_std_compute_wind_2d_fm`) carries
the same factors — it shares the parent horizontal grid, and its vertical
transport is the implicit column solver, so only the horizontal
flux/denominator placement applies.

*Not* m-corrected (O(m−1) on diagnostics, accepted): the slope terms in
`calc_w_real` (kinematic w_real diagnosis), Smagorinsky deformation, Sx/
thermal-wind adjustments, and the MPDATA advection path (`adv = "std"`
only). The vertical λ-gradient in the w correction carries no factor *by
construction*: z is true meters (map factors are horizontal-projection
metrics), and the divergence operator's vertical term is equally m-free
(column-constant m cancels between W_m and the cell denominator), so D
and G remain mutually consistent in the vertical.

Constant-flow preservation: for a uniform velocity field the staggered-CV
flux divergence telescopes to the average of the two adjacent mass cells'
mass-flux divergences, which the projection closure has zeroed — exactly in
x and y; the vertical term uses the CV-averaged layer thickness
`dz_u = ½(dz(i−1)+dz(i))`, leaving an O(∂dz/∂x) metric residual over slopes
(bounded, advective, re-absorbed by the next projection). This is the
standard level of metric approximation for staggered momentum advection in
terrain-following models.

### Control volumes

* **u-CV** centred at u-point i: x-faces at mass centres i−1, i; y-faces at
  the cell corners (i, j), (i, j+1); z-faces at u-column interfaces.
  Face mass fluxes are 2-point averages of the native fluxes, e.g.
  x-face flux at mass centre i: `½(Û_m(i) + Û_m(i+1))` where `Û_m` is the
  u-point flux recomputed locally.
* **v-CV**: symmetric.
* **w_real-CV**: the mass cell itself — `w_real` is cell-centred, so it is
  advected by the *unmodified scalar kernel* (`adv_std_advect3d` with
  `U_m, V_m, W_m, denom`), with no flux corrector and no cz-diff.

## 4. Time integration

One RANS dynamics step runs inside `advect()`, once per model time step Δt:

```
1. projection (wind.F90::rans_project)     — FIRST, before any transport
2. scalar advection (existing RK3)
3. momentum RK3 (wind_rans):  u,v,w_real advected; buoyancy on w_real;
   geostrophic forcing f·k×(u_e − u) on u,v (§2)
4. horizontal Smagorinsky + Rayleigh lid damping (wind_rans)
5. grid-w increment update (§4.2a)
```

**Ordering is load-bearing: project, THEN advect.** The PBL momentum
tendencies and the lateral boundary nudge modify u/v between advect calls,
leaving divergence in the field. Flux-form momentum advection telescopes to
zero spurious source *only* when the transport satisfies discrete
continuity; advecting with divergent transport self-amplifies momentum
wherever div ≠ 0 (this mirrors diagnostic mode, which calls `balance_uvw`
before advect every step). The small divergence the dynamics step itself
introduces is projected out at the start of the next step.

### 4.1 RK3 momentum advection

Wicker–Skamarock RK3 matching the scalar driver: stage factors
t_fac = 1/3, 1/2, 1, each stage evaluating fluxes from the *previous stage's*
field, always restarting from the time-n state. Transport mass fluxes are
frozen at time n (built from the projected, divergence-free time-n
velocities) — identical policy to scalar advection. Because the transport
field satisfies `∇·(ρu)=0` discretely (guaranteed by the previous step's
column-integration closure), flux form and advective form coincide and the
scheme is conservative with no spurious mass source.

Spatial order: 3rd-order WRF upwind fluxes in the horizontal and vertical
interior,

    F_face = m·[7(q_up + q_dn) − (q_upup + q_dndn)]/12
             − |m|·[3(q_dn − q_up) − (q_dndn − q_upup)]/12

degrading to 1st-order upwind at the first interior face adjacent to the
top/bottom boundaries (same policy as `flux3`'s vertical boundary slabs).
Horizontal stencils near lateral tile edges read across the halo (halo width
≥ 2 holds whenever scalar `h_order ≥ 3`).

Buoyancy `B = g·θ′/θ̄` (cell-centred, evaluated once per step from the
already-advected θ) is added to the w_real stage update with the same
t_fac·Δt weight. The geostrophic forcing (§2) is likewise added per stage
with t_fac·Δt weighting, evaluated from the *previous stage's* winds (the
same policy as the advective tendency); the ambient state u_e is frozen
across the step (f·Δt ~ 1e-4, so sub-step variation is negligible).

Global-boundary faces (u at i = ids, ide+1; v at j = jds, jde+1) are owned by
the lateral BC relaxation and are never updated by the dynamics. w_real and
interior faces in the relaxation zone *are* integrated; the Davies nudge
(§7) blends them toward forcing afterwards.

After each stage: halo exchange of u, v, w_real (`exch_var`, corners).

### 4.2a Grid-w increment (vertical-momentum carrier)

The prognostic divergence-bearing vertical variable is **grid-w**, not
w_real. At the end of the momentum step, the *increment* of physical w over
the step (advection + buoyancy + diffusion + damping, small and smooth) is
staggered to interfaces (dz-weighted, divided by J_w) and **added** to the
grid-w left by the previous closure.

This replaces an earlier design that re-derived grid-w wholesale from
w_real via the contravariant identity (`calc_idealized_wgrid`) at every
projection. That conversion is *not* the discrete inverse of the closure's
`calc_w_real` diagnostic: over steep cells the round trip re-injects
O(w·slope) divergence at the same cells every step, and the Poisson solve
"corrects" it by persistently accelerating the horizontal wind — a secular
instability (~1 m/s per step at a 250 m Gaudergrat ridge crest) that no
amount of advection-side dissipation removes. Carrying the
continuity-exact grid-w and converting only increments kills the feedback;
the first-step predictor divergence drops from O(1e-4) to O(1e-10).

### 4.2b The exact discrete operator (probed)

A projection is the composition of the model's divergence operator D
(`calc_divergence`) with the velocity-correction operator G
(`calc_updated_winds`, `w_to_grid` form). The matrix the solver inverts
must equal **A = 2·D∘G as discrete operators** — every interpolation and
metric placement included — or the corrections leave an
O(operator-mismatch) divergence residual in the transport that flux-form
momentum advection cannot tolerate (an "approximate projection"). The
legacy analytic stencil assembly (`initialize_coefs`) differs from D∘G
through single-column `dz_if`/`advection_dz` shortcuts, ρ placement, and
independently discretized cross terms; Richardson iteration over the
approximate operator was tested and is NOT sufficient over steep terrain.

A is therefore extracted numerically at init (once per nest): apply G
then D to 27 lattice-colored indicator fields (3×3×3 coloring — no two
stencil supports overlap) and read off all 15 coefficients per cell
(`probe_*` routines in `wind_iterative.F90`, orchestrated at the end of
the first `rans_project`). Off-stencil leakage is measured during
probing and is exactly zero. Under the probed operator (`rans_probed`):
ghost planes are identity rows; G uses one-sided vertical λ-gradients at
the bottom/top levels and a compact interface gradient for the w
correction (no ghost-λ reads; the lid interface gets no correction —
rigid lid). The probe is bit-consistent with whatever D and G implement,
so future changes to either are absorbed by re-probing.

With the exact operator, **one solve per step** projects to solver
tolerance (verified: the re-derived divergence after one solve is ~1e-5
of the predictor divergence; with the analytic operator it was O(1)).
Health metric: the per-step solve's initial residual must decay over the
first simulated minutes — a constant res₀ indicates a standing
(boundary) divergence source.

**The diagnostic solver uses the same probed operator** (the correction
operator is unified: corrections land on u, v, and grid-w with the
compact interface gradient; the legacy w_real-correction form is
removed). Its Froude-dependent α enters the composition, so the operator
is re-probed at each wind update when α varies (skipped for
`alpha_const > 0` after the first probe). The probe requires the solver
workspace, which the first solve allocates — so the first-ever update
runs solve → calibrate → polishing solve, and later updates re-probe up
front. The historical `wind_iterations` loop and `balance_uvw` final
closure remain only in diagnostic mode semantics: the loop is gone (it
compensated for operator mismatch) while `balance_uvw` stays as the
diagnostic w contract.

**Warm start + normalized tolerance**: `bicgstab_solve` computes
r₀ = b − A·x₀ properly (one extra SpMV per solve), the convergence
target is `max(tol_abs, tol_rel·‖b‖)` — normalized by the RHS, not by
r₀, so a good initial guess directly reduces iterations — and the RANS
path keeps the previous step's λ as x₀. Measured effect: per-step
iterations fall from a constant ~250–340 (cold) to ~90–120 as the flow
settles.

### 4.2 Projection

Chorin-style non-incremental projection, once per step (the multiplier is
not carried between steps; with RK3 advection the splitting error is O(Δt)
in the pressure term, acceptable at HICAR's Δt ≈ O(1 s)):

1. `u*, v*, w_grid* → dqdt_3d` (full memory range; halos are current).
2. `calc_divergence(div, use_dqdt=.True., horz_only=.False.)` →
   `div = ∇·(ρ J u*)` in the model's discrete operator.
3. `calc_iter_winds(domain, α≡1, div, adv_den)` — the existing BiCGStab
   solve of the variational system, which with α ≡ 1 collapses to the
   pressure-Poisson problem. Applies
   `u ← u + (∇λ − metric)/2ρ` to `u%dqdt, v%dqdt` and the λ_z correction to
   `w_real%data_3d`. λ relates to the perturbation pressure by
   λ = −2Δt·π (not stored).
4. Halo exchange (u, v dqdt; w_real).
5. `balance_uvw(update=.True.)` — recompute `w%dqdt` by vertical integration
   of the horizontal divergence: exact discrete continuity per column.
6. `calc_w_real(u%dqdt, v%dqdt, w%dqdt → w_real%data_3d)` — final physical w
   consistent with the divergence-free field. This *is* the prognostic w
   carried to the next step: it equals w_real* + projection correction up to
   the BiCGStab residual, so vertical momentum (incl. buoyancy work) is
   retained.
7. dqdt → data_3d for u, v, w over the full memory range; final halo
   exchange.

Boundary conditions of the Poisson solve are inherited unchanged from the
diagnostic solver (λ = 0 ring at global lateral boundaries, terrain/lid
handling in the stencil assembly). Mass imbalance imposed by the lateral
boundaries exits through the lid exactly as in diagnostic HICAR.

### 4.3 Time step

`compute_dt`'s CFL (|u|/Δx + |v|/Δx + |w|/Δz, with the user CFL factor)
covers 3rd-order RK3 advection of momentum just as it does scalars. Because
the winds now evolve *between* forcing updates, dt is recomputed **every
time step** under `kRANS_WINDS` (in `time_step.F90`), not only at wind
updates. The explicit Smagorinsky diffusion is clipped to its own stability
bound (§6) rather than constraining dt.

Recommended `cfl_reduction_factor ≤ 1.2` for RANS runs: the momentum RK3 +
3rd-order upwind combination has a 1-D stability limit of ~1.6, and the
scalar-tuned 1.6 used by existing namelists sits at that edge with no FCT
safety net on momentum.

Debug kill switches (environment variables, read once at the first step):
`HICAR_RANS_NO_ADV`, `NO_BUOY`, `NO_SMAG`, `NO_RAYLEIGH`, `NO_ADV_Z`,
`NO_ADV_XY`, `FIRST_ORDER` (=1 to activate) disable individual terms for
instability bisection without rebuilding.

## 5. Buoyancy and reference state

θ̄ = `adv_theta_ref` (domain-mean profile mapped through local z; zero
horizontal gradients in physical space by construction). Buoyancy uses dry θ
(moisture buoyancy `0.61 qv − qc` deferred; document if/when added). Because
the projection immediately absorbs the hydrostatically balanced part of any
horizontally uniform buoyancy error into λ, only *horizontal variations* of
θ′ drive circulation — which is the physically correct behavior.

### Thermodynamic pressure

The dynamics never use a thermodynamic pressure (anelastic: the dynamic
perturbation pressure is the projection multiplier), but the physics do —
exner/temperature conversion, density, microphysics saturation. Pinning the
interior pressure to the forcing (the pre-existing treatment) lets it drift
out of consistency with the prognostic θ. Instead, under `kRANS_WINDS` the
pressure is re-diagnosed hydrostatically from the model state
(`domain_obj.F90::rans_diagnose_pressure`, run at the top of every
non-forcing-mode `diagnostic_update`, i.e. at the step-start call after
`apply_forcing` has advanced the top anchor, at the initialization ingest,
and at the cheap `thermo_only` refresh before microphysics — so p, exner,
T, ρ track the post-advection θ_v the microphysics is given): column-wise
Exner integration

    dπ_E/dz = −g / (cp·θ_v),   θ_v = θ·(1 + 0.61·qv)

downward from the model top, whose pressure stays pinned to the
time-advanced forcing. The synoptic pressure signal thus enters through
the top anchor and the model's own thermal structure. This respects design
principle 1: π_E is O(1) with O(1e-3) per-layer increments — safe in
single precision — and the resulting p is only consumed pointwise, never
differenced horizontally.

## 6. Horizontal Smagorinsky diffusion

2-D deformation-based eddy viscosity at mass points, time-n velocities:

    D11 = ∂u/∂x − ∂v/∂y,  D12 = ∂u/∂y + ∂v/∂x      (along model levels)
    K   = (C_s Δx)² · √(D11² + D12²),   C_s = 0.21
    K   ← min(K, 0.12·Δx²/Δt)                        (explicit stability clip)

Applied as a forward-Euler horizontal Laplacian (along model levels) to u, v
and w_real after the RK3 loop. Along-level rather than constant-z diffusion
is a known approximation over steep slopes; θ already has the constant-z
machinery (`cz_diff`) and extending it to momentum is possible follow-up
work. This term is not optional — without explicit horizontal dissipation,
2Δx energy accumulates and any prognostic core eventually goes unstable.

Vertical mixing is owned entirely by the PBL scheme (YSU applies u/v
tendencies under `kRANS_WINDS`); the dynamics core adds none, avoiding
double-counting.

## 7. Lateral and top boundaries

* **Lateral**: u, v get `force_boundaries = .True.` under `kRANS_WINDS`
  (metadata), so `apply_forcing` runs its boundary-only Davies relaxation,
  hard-set on the outermost ring; the interior is never touched by forcing.
  The relaxation RATE is `relax_filter·Δt/τ` with τ = `kRANS_LATERAL_TAU`
  (300 s), NOT the scalars' per-hour rate: the ring must absorb gravity-wave
  energy (period 2π/N ≈ 10 min), and at the per-hour rate it cannot.
  Matching state damping for w_real: implicit lateral Rayleigh
  `w_real ← w_real / (1 + relax_filter·Δt/τ)` applied with the blending
  (below) in `rans_momentum_step`. Found in the Agnesi validation: the
  increment blending slows wave propagation through the taper, so without
  a state sink on the wave-period timescale, incident wave energy SHOALS
  there — |w| at the inflow taper grew linearly (~0.16 m/s per hour, to
  1 m/s in 4 h on the hm=100 case; ring-edge jets and a shear instability
  within 1–2 h on hm=500). With the absorber both cases reach a clean
  steady state; τ = 100 s gave results identical to 300 s (the residual
  ~0.13 m/s taper dipole is the steady incident-flux/absorber
  equilibrium, not absorber-limited).
  **`relax_filters = .True.` is required in practice** (the solver warns at
  init if not): with it off, the filter is a hard step (1 in the outer halo
  ring, 0 elsewhere) and the resulting shear line destabilizes the
  prognostic dynamics.
* **Blending zone**: at the end of the momentum step, the full dynamics
  increment of u, v, w_real is multiplied by `(1 − relax_filter)` — the
  solution is forcing-dominated in the ring, dynamics-dominated inside
  (the standard specified/blended-zone treatment).
* **Balanced forcing targets (option D, EULAG-style — the adopted
  design)**: at every wind update, `rans_balance_forcing` projects the
  interpolated forcing wind field at BOTH ends of the update interval
  (rotated to grid-relative first), so the Davies relaxation targets
  themselves satisfy the model's discrete continuity constraint —
  relaxing toward them injects no divergence, exactly as EULAG's
  "ambient state" absorbers relax toward anelastically balanced states.
  Linear time-interpolation between balanced endpoints stays balanced
  (the divergence operator is linear in velocity at fixed ρ, J). The
  matching grid-w target pair (`w_forc_now/w_forc_tend` in wind.F90)
  provides the ring w relaxation target. With balanced targets, the
  ring needs NO special treatment by the solve: no RHS masking, no
  specified-zone exclusion, no diagnostic/O'Brien ring w. Validated:
  24+ simulated minutes on Gaudergrat 250 m, fields saturating (ridge
  jet ~50–56 m/s free-slip-like, |w| ≤ 18 at mid-levels, T(lid)
  constant to ±0.1 K).
  Failed alternatives, recorded so they are not retried: relaxing
  toward RAW forcing winds creates a standing divergence source the
  solve pumps against (zone-edge runaway); zeroing ring w + masking the
  projection RHS leaves divergent transport in the ring (momentum
  source at the zone edge, eventual NaN); raw column closure dumps the
  column divergence into ring lid-w → lid heating/updraft feedback
  (+63 K, 94 m/s in 30 min); O'Brien-adjusted closure without masking
  feeds the distributed residual to the solve (column-wide pumping);
  O'Brien + masking ran 12 min but grew deep w columns at the zone's
  inner edge. The specified-zone-exclusion formulation (λ = 0 identity
  rows in the ring) was implemented and then retired unused once
  option D made it unnecessary.
  The forcing winds are rotated to grid-relative (`make_winds_grid_relative`)
  *before* balancing, inside `rans_balance_forcing`, so the relaxation
  targets — and the geostrophic ambient state u_e (§2), which reads the
  same arrays — are grid-relative balanced states.
* **w_real**: not forced after initialization (`apply_base_from_forcing` is
  *only* called at the first wind update under RANS — on later calls it
  would overwrite the prognostic w_real with forcing data).
* **Top**: Rayleigh damping of w_real over the top of the domain,

      w_real ← w_real / (1 + Δt/τ(z)),
      τ(z)⁻¹ = τ_max⁻¹ · sin²( (π/2)·(z − z_d)/(z_top − z_d) ),  z > z_d

  with z_d = 0.75·z_top and τ_max = 30 s (implicit form, unconditionally
  stable). The projection afterwards keeps the field non-divergent.

## 8. Initialization (first wind update)

1. `apply_base_from_forcing` seeds u, v dqdt with the interpolated forcing
   (and w_real with forcing w or 0).
2. Rotate dqdt to grid-relative (`make_winds_grid_relative`), halo exchange.
3. Copy dqdt → data for u, v.
4. `balance_uvw(update=.False.)`: initial grid-w from column integration —
   the initial state is discretely divergence-free in the column sense.
5. `calc_w_real` for the initial physical w.

No special spin-up is required: the first projection removes any remaining
horizontal-divergence imbalance.

## 9. What is deliberately NOT in v1

**Resolved — large-scale forcing of the interior.** Earlier versions carried
no synoptic pressure-gradient force, so the interior drifted from the
boundary inflow on the advective timescale and the contrast re-intensified
at the relaxation-zone edge (observed in the first Gaudergrat validations).
This is now addressed by the geostrophic ambient-state forcing (§2), which
was option (c) of the remedies considered — Coriolis + geostrophic forcing
folded into one term. The weak interior nudge (option (a),
`HICAR_RANS_NUDGE=1`, τ = 1800 s in `domain_obj.F90::apply_forcing`)
remains as an opt-in experiment/mitigation switch only.

| Item | Why deferred | Where it would go |
|---|---|---|
| Moist buoyancy (qv, q_hydrometeor loading) | small vs θ′ for first validation | buoyancy term, §5 |
| Incremental projection (carry λ) | O(Δt) pressure splitting acceptable; needs λ persistence across steps | `rans_project` |
| Momentum FCT / monotone limiting | velocities are not sign-preserving quantities; 3rd-order upwind dissipation suffices | momentum kernels |
| Constant-z momentum diffusion over steep slopes | reuse of `cz_diff` LUTs for staggered grids is non-trivial; error shrinks with Δx | §6 |
| Namelist exposure of C_s, τ, z_d | validate defaults first | `opt_types` + `namelist_utilities` |
| 3-D turbulence closure (D13/D23 cross stresses) | split 2-D-horizontal + 1-D-PBL closure is the standard at the target Δx > 1 km | §6 |
| Vertical turbulent mixing of w | negligible at mesoscale; WRF PBL schemes omit it too | momentum step |
| Map-scale factors | HICAR-wide approximation (~1–2% at few-100-km extents) | host model |

## 10. Validation sequence

1. **Compile + unit tests** (`make check`).
2. **Uniform flow, flat terrain, neutral profile**: the dynamics must
   preserve the initial state to roundoff (advection of a constant field,
   zero buoyancy, projection a no-op).
3. **Idealized hill (Witch of Agnesi)**: stationary mountain-wave response;
   compare against the diagnostic variational solver and linear theory;
   check w amplitude/phase and absence of 2Δx noise.
4. **Real-terrain case**: stability over ≥ 24 h, energy bounded, boundary
   inflow/outflow well-behaved; compare against diagnostic-solver run.
