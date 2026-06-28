#!/usr/bin/env python3
"""Spatial difference maps (C++ reference - Fortran port) for the SNOWPACK
comparison runs. Produces, for each key variable, a 3-panel row:
  [ C++ field | Fortran field | difference (cpp - port) ]
at the final output frame, masked to interior snow cells.

Usage:
  plot_snowpack_diff.py <cpp.nc> <fortran.nc> <out_dir>
"""
import sys, os
import numpy as np
import xarray as xr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

cpp_f, for_f, out_dir = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(out_dir, exist_ok=True)
c = xr.open_dataset(cpp_f)
f = xr.open_dataset(for_f)

# interior-snow mask from final frame
sh = c["snow_height"].isel(time=-1).values
mask = sh > 0.1

# (name, vertical-layer-index or None, label, units)
VARS = [
    ("rsds",                None, "downwelling shortwave",       "W m$^{-2}$"),
    ("albedo",              None, "surface albedo",              "-"),
    ("tsfe",                None, "skin temperature",            "K"),
    ("snow_temperature_i",  0,    "surface node temp (node 0)",  "K"),
    ("snow_temperature_i",  1,    "node 1 temp",                 "K"),
    ("snow_temperature_i",  2,    "node 2 temp",                 "K"),
    ("snow_temperature",    1,    "element 1 temp",              "K"),
    ("ustar",               None, "friction velocity",           "m s$^{-1}$"),
    ("coeff_heat_exchange", None, "heat-exchange coeff",         "-"),
    ("hfss",                None, "sensible heat flux",          "W m$^{-2}$"),
    ("hfls",                None, "latent heat flux",            "W m$^{-2}$"),
    ("snow_height",         None, "snow depth",                  "m"),
]


def get(ds, name, lev):
    a = ds[name].isel(time=-1).values
    if lev is not None:
        a = a[lev]
    return np.where(mask, a, np.nan)


def robust_lim(*arrs):
    v = np.concatenate([a[np.isfinite(a)].ravel() for a in arrs])
    return np.nanpercentile(v, 1), np.nanpercentile(v, 99)


def sym_lim(d):
    v = np.nanpercentile(np.abs(d[np.isfinite(d)]), 99)
    return -v, v


for name, lev, label, units in VARS:
    if name not in c or name not in f:
        print(f"  skip {name} (absent)")
        continue
    A = get(c, name, lev)
    B = get(f, name, lev)
    D = A - B
    tag = f"{name}" + (f"_L{lev}" if lev is not None else "")

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6), constrained_layout=True)
    vmin, vmax = robust_lim(A, B)
    for k, (arr, ttl) in enumerate([(A, "C++ reference"), (B, "Fortran port")]):
        im = ax[k].imshow(arr, origin="lower", cmap="viridis", vmin=vmin, vmax=vmax)
        ax[k].set_title(ttl)
        fig.colorbar(im, ax=ax[k], shrink=0.8, label=units)
    dmin, dmax = sym_lim(D)
    im = ax[2].imshow(D, origin="lower", cmap="RdBu_r", vmin=dmin, vmax=dmax)
    mae = np.nanmean(np.abs(D))
    ax[2].set_title(f"difference (C++ - port)\nmean|Δ|={mae:.3e} {units}")
    fig.colorbar(im, ax=ax[2], shrink=0.8, label=f"Δ {units}")
    for a in ax:
        a.set_xticks([]); a.set_yticks([])
    fig.suptitle(f"{label}  [{tag}]  — final frame (1 h)", fontsize=13)
    fp = os.path.join(out_dir, f"diff_{tag}.png")
    fig.savefig(fp, dpi=110)
    plt.close(fig)
    print(f"  wrote {fp}  (mean|Δ|={mae:.3e} {units})")

print(f"\nDone -> {out_dir}")
