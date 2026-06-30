#!/usr/bin/env python3
"""Seed a HICAR domain file with an initial SNOWPACK state.

The SNOWPACK driver (both the C++-bindings path, src/physics/sm_SNOWPACK.F90, and
the native-Fortran path, snowpack_driver.F90) reads its initial per-layer snow
state from the domain / init_conditions_file via read_snowpack_state, using the
namelist `snowpack_*_var` names. This script copies an existing domain file and
adds those variables so the model starts from a prescribed, multi-layer snowpack
instead of cold-starting. It is used by the C++-vs-Fortran SNOWPACK comparison
test (tests/scripts/snowpack/test_snowpack_compare.sh) so both builds initialize from a bit-
identical state.

Layout expected by HICAR's reader (verified against try_read_snp_3d):
  * 3-D per-layer variables are stored as (snow_layer, y, x) in CDL order, so
    io_read returns (x, y, snow_layer) — matching the model's internal
    (i=x, k=layer, j=y) after its [1,3,2] reshape.
  * snow_layer index 1 is the SURFACE element; index nlayers is the bottom
    (ground) element. HICAR's internal convention is k=1 = surface (e.g. the
    skin diagnostic is tn_w(i,1,j), and the C++-bindings wrapper reverses with
    (n_elem:1:-1) when filling the ground-up-numbered C++ SnowStation). The
    reader uses the data as-read, so the file must be written surface-first.
    [2026-06-11: previous versions wrote base-first, giving both builds an
    upside-down pack.]
  * snowpack_nlayers is a 2-D int variable (y, x).
  * depositionDate holds ABSOLUTE Modified Julian Dates (both drivers compute
    layer age as now_mjd - depositionDate; the port uses domain%sim_time%mjd()
    and the C++ wrapper passes domain MJD into Mdata.date). Writing relative
    ages makes age ~ 2.4e6 days and silently disables the loading-rate stress
    Sig0 in both builds.
  * snow_stress holds the model-convention NEGATIVE Cauchy stress
    -g*cos(slope)*(mass_above + m_layer/2), so the step-1 stress increment
    (SigC - oldStress) is ~0 and no phantom-loading CDot/Sig0 is generated.

The seeded pack is a physically self-consistent, settled profile (so neither
implementation immediately restructures it): density / grain size / temperature
ramp from a warmer, denser, coarser base to a colder, lighter, finer surface.
Seeded on every land cell (landmask == 1).
"""
import argparse
import datetime
import sys

import numpy as np
import netCDF4

# Physical constants / pack parameters
RHO_ICE = 917.0          # kg/m^3, density of pure ice (SNOWPACK Constants::density_ice)
T_FREEZE = 273.15        # K
MJD_EPOCH = datetime.datetime(1858, 11, 17)  # Modified Julian Date epoch

# Base (ground) -> surface ramp endpoints for the layered profile
RHO_BASE, RHO_SURF = 330.0, 180.0      # bulk snow density [kg/m^3]
T_BASE, T_SURF = 272.65, 265.15        # snow temperature [K] (all sub-freezing)
RG_BASE, RG_SURF = 0.60, 0.25          # grain radius [mm]
DD_BASE, DD_SURF = 0.0, 0.20           # dendricity [-] (0 = fully rounded)
SP_BASE, SP_SURF = 1.0, 0.55           # sphericity [-]
RB_FRAC = 0.20                          # bond radius as a fraction of grain radius
N3_VALUE = 3.0                          # coordination number [-]
GRAVITY = 9.80665                       # m/s^2 (SN_G / SNOWPACK Constants::g)

# Map of (kVARS internal / file variable name) -> namelist option key.
# The namelist generator (SNOWPACK_Compare.sh) points each snowpack_*_var at the
# file variable name on the left, so keep these names in sync with that script.
VAR3D_NAMES = [
    "Ds", "Vol_Frac_I", "Vol_Frac_W", "Vol_Frac_A", "Vol_Frac_S", "Vol_Frac_WP",
    "snow_temperature", "Rg", "Rb", "Dd", "Sp", "mk", "CDot", "snow_stress",
    "depositionDate", "N3",
]


def build_profiles(n_layers, total_depth, start_mjd):
    """Return per-layer 1-D profiles in FILE order (index 0 = SURFACE, index
    n-1 = base, matching HICAR's k=1 = surface convention). Internally built
    base-first and flipped on return.

    snow_stress is returned WITHOUT the cos(slope) factor (flat-terrain
    Cauchy stress); the caller multiplies by cos(slope) per cell."""
    # Fractional position from base (0) to surface (1) at layer centers.
    frac = (np.arange(n_layers) + 0.5) / n_layers

    def ramp(base, surf):
        return base + (surf - base) * frac

    ds = np.full(n_layers, total_depth / n_layers)        # uniform thickness
    rho = ramp(RHO_BASE, RHO_SURF)
    vfi = rho / RHO_ICE                                   # volumetric ice fraction
    vfw = np.zeros(n_layers)                              # dry pack
    vfs = np.zeros(n_layers)                              # no soil in snow elements
    vfa = 1.0 - vfi - vfw - vfs                           # air fills the remainder
    vfwp = np.zeros(n_layers)
    tsnow = ramp(T_BASE, T_SURF)
    rg = ramp(RG_BASE, RG_SURF)
    rb = RB_FRAC * rg
    dd = ramp(DD_BASE, DD_SURF)
    sp = ramp(SP_BASE, SP_SURF)
    mk = np.zeros(n_layers)                               # marker: ordinary snow
    cdot = np.zeros(n_layers)                             # d(stress)/dt
    n3 = np.full(n_layers, N3_VALUE)

    # Cauchy stress [Pa], model convention (negative = compression), flat
    # terrain: SigC = -g * (mass strictly above + half this layer's own mass).
    # This is exactly what both drivers recompute each step (compSnowCreep /
    # the port settling loop), so the step-1 increment (SigC - oldStress) ~ 0
    # and no phantom-loading CDot/Sig0 is generated. cos(slope) applied per
    # cell by the caller.
    layer_mass = rho * ds                                 # kg/m^2 per layer
    mass_above = np.cumsum(layer_mass[::-1])[::-1] - layer_mass
    snow_stress = -GRAVITY * (mass_above + 0.5 * layer_mass)

    # Deposition date as absolute MJD: deeper layers were deposited earlier
    # (base ~60 days before the run start, surface ~0). Both drivers compute
    # age = now_mjd - depositionDate, so this yields real ages of 2-58 days.
    deposition = start_mjd - ramp(60.0, 0.0)

    # Interface temperatures (nlayers+1 values: one below each element plus one
    # above the top). Provided explicitly so the driver's interface-fill path is
    # skipped — that path indexes snow_temperature(:,nz,:) with nz=0 on water
    # cells and goes out of bounds. interface[k] mirrors what the driver would
    # compute: base = element 1, interior = element average, top = element nz.
    tsnow_i = np.empty(n_layers + 1)
    tsnow_i[0] = tsnow[0]
    tsnow_i[1:n_layers] = 0.5 * (tsnow[:-1] + tsnow[1:])
    tsnow_i[n_layers] = tsnow[-1]

    profiles = {
        "Ds": ds, "Vol_Frac_I": vfi, "Vol_Frac_W": vfw, "Vol_Frac_A": vfa,
        "Vol_Frac_S": vfs, "Vol_Frac_WP": vfwp, "snow_temperature": tsnow,
        "Rg": rg, "Rb": rb, "Dd": dd, "Sp": sp, "mk": mk, "CDot": cdot,
        "snow_stress": snow_stress, "depositionDate": deposition, "N3": n3,
    }
    # Flip everything to FILE order: index 0 = surface (HICAR k=1 = surface).
    profiles = {k: v[::-1].copy() for k, v in profiles.items()}
    return profiles, tsnow_i[::-1].copy()


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("src", help="source domain NetCDF (e.g. Gaudergrat_250m.nc)")
    p.add_argument("dst", help="output domain NetCDF to write (seeded)")
    p.add_argument("--n-layers", type=int, default=14,
                   help="snow elements (default 14). Kept moderate (~0.07 m thick "
                        "for a 1 m pack) so the layers sit comfortably inside "
                        "SNOWPACK's rezoning thresholds and neither implementation "
                        "merges/splits them differently on the first step.")
    p.add_argument("--depth", type=float, default=1.0, help="total snow depth [m] (default 1.0)")
    p.add_argument("--layer-dim", default="snow_layer", help="name of the new snow-layer dimension")
    p.add_argument("--landmask-var", default="landmask", help="2-D land mask variable (1 = land)")
    p.add_argument("--slope-var", default="slope_rad",
                   help="2-D slope variable [radians] used for the cos(slope) factor "
                        "in the seeded Cauchy stress (default slope_rad)")
    p.add_argument("--start-date", default="2017-02-14_10:00",
                   help="simulation start (YYYY-MM-DD_HH:MM, UTC); deposition dates are "
                        "written as MJD relative to this (default 2017-02-14_10:00, the "
                        "shifted-midday SNOWPACK comparison window)")
    args = p.parse_args()

    start_dt = datetime.datetime.strptime(args.start_date, "%Y-%m-%d_%H:%M")
    start_mjd = (start_dt - MJD_EPOCH).total_seconds() / 86400.0

    profiles, tsnow_i = build_profiles(args.n_layers, args.depth, start_mjd)

    with netCDF4.Dataset(args.src) as src, netCDF4.Dataset(args.dst, "w") as dst:
        # --- copy the source verbatim (dims, attrs, vars) -------------------
        dst.setncatts({k: src.getncattr(k) for k in src.ncattrs()})
        for name, dim in src.dimensions.items():
            dst.createDimension(name, len(dim) if not dim.isunlimited() else None)
        for name, var in src.variables.items():
            out = dst.createVariable(name, var.datatype, var.dimensions,
                                     zlib=getattr(var, "_FillValue", None) is not None)
            out.setncatts({k: var.getncattr(k) for k in var.ncattrs()})
            out[:] = var[:]

        if args.landmask_var not in dst.variables:
            sys.exit(f"ERROR: land mask variable '{args.landmask_var}' not in {args.src}")
        land = np.asarray(dst.variables[args.landmask_var][:]) == 1   # (y, x)
        ny, nx = land.shape
        nlay = args.n_layers

        # --- new snow-layer dimension --------------------------------------
        if args.layer_dim in dst.dimensions:
            sys.exit(f"ERROR: dimension '{args.layer_dim}' already exists in {args.src}")
        dst.createDimension(args.layer_dim, nlay)

        # Reference an existing (y, x) variable's dimension names so we match the
        # file's spatial dim ordering exactly.
        yx_dims = dst.variables[args.landmask_var].dimensions   # ('y', 'x')

        # --- snowpack_nlayers (2-D int): nlay on land, 0 elsewhere ---------
        nl = dst.createVariable("snowpack_nlayers", "i4", yx_dims)
        nl.long_name = "SNOWPACK number of snow elements (initial state)"
        nl[:] = np.where(land, nlay, 0).astype("i4")

        # cos(slope) for the per-cell Cauchy stress
        if args.slope_var in dst.variables:
            cos_sl = np.cos(np.asarray(dst.variables[args.slope_var][:], dtype="f8"))
        else:
            print(f"WARNING: slope variable '{args.slope_var}' not found; "
                  "seeding snow_stress with cos(slope)=1")
            cos_sl = np.ones((ny, nx))

        # --- 3-D per-layer state vars: (snow_layer, y, x) ------------------
        layer_dims = (args.layer_dim,) + yx_dims
        for vname in VAR3D_NAMES:
            prof = profiles[vname].astype("f4")                 # (nlay,) surface->base
            field = np.zeros((nlay, ny, nx), dtype="f4")
            # broadcast the column profile over all land cells; 0 over water
            if vname == "snow_stress":
                field[:, land] = prof[:, None] * cos_sl[land][None, :].astype("f4")
            else:
                field[:, land] = prof[:, None]
            v = dst.createVariable(vname, "f4", layer_dims)
            v[:] = field

        # --- snow_temperature_i: (snow_layer_i = nlay+1, y, x) -------------
        # The element interface temperatures live on a grid one larger than the
        # elements. Provided so the driver skips its interface-fill loop (which
        # is out-of-bounds on nz=0 water cells).
        layer_i_dim = args.layer_dim + "_i"
        dst.createDimension(layer_i_dim, nlay + 1)
        field_i = np.zeros((nlay + 1, ny, nx), dtype="f4")
        field_i[:, land] = tsnow_i.astype("f4")[:, None]
        vi = dst.createVariable("snow_temperature_i", "f4", (layer_i_dim,) + yx_dims)
        vi[:] = field_i

        # provenance
        dst.snowpack_seed = (
            f"seeded by make_snowpack_init.py: {nlay} layers, {args.depth} m, "
            "surface-first layer order (HICAR k=1 = surface), depositionDate as MJD "
            f"(start {args.start_date}), snow_stress = -g*cos(slope)*(mass_above+m/2), "
            "on all land cells")

    print(f"Wrote {args.dst}: {nlay}-layer / {args.depth} m snowpack on "
          f"{int(land.sum())} land cells (of {land.size}).")


if __name__ == "__main__":
    main()
