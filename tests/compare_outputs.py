#!/usr/bin/env python3
"""Compare two HICAR NetCDF output files for reproducibility testing.

Usage:
    compare_outputs.py FILE1 FILE2 --tolerance TOL [--last-timestep-only] [--figures-dir DIR]
"""
import argparse
import os
import sys
import numpy as np
import xarray as xr


class bcolors:
    BLUE = '\033[0;36m'
    GREEN = '\033[0;32m'
    RED = '\033[0;31m'
    NC = '\033[0m'


def plot_diff(a1, a2, var_name, label1, label2, figures_dir):
    """Save a comparison figure for a variable that failed the tolerance check."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print(f"  {bcolors.RED}matplotlib not available, skipping figure for {var_name}{bcolors.NC}")
        return

    diff = a1 - a2

    # For 3D+ arrays, take a 2D slice: last time, mid-level
    if diff.ndim >= 4:
        # (time, level, lat, lon) — last time, mid-level
        a1 = a1[-1, a1.shape[1] // 2, :, :]
        a2 = a2[-1, a2.shape[1] // 2, :, :]
        diff = diff[-1, diff.shape[1] // 2, :, :]
        slice_label = f"t=-1, level={a1.shape[0] // 2 if a1.ndim > 2 else 'mid'}"
    elif diff.ndim == 3:
        # Could be (time, lat, lon) or (level, lat, lon)
        a1 = a1[-1, :, :]
        a2 = a2[-1, :, :]
        diff = np.sum(diff[:, :, :],axis=0)
        slice_label = "t=-1"
    elif diff.ndim == 2:
        slice_label = "2D field"
    else:
        # 1D or scalar — skip plotting
        return

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    vmin = min(np.nanmin(a1), np.nanmin(a2))
    vmax = max(np.nanmax(a1), np.nanmax(a2))

    im0 = axes[0].imshow(a1, origin='lower', vmin=vmin, vmax=vmax)
    axes[0].set_title(f"{label1}")
    plt.colorbar(im0, ax=axes[0], shrink=0.8)

    im1 = axes[1].imshow(a2, origin='lower', vmin=vmin, vmax=vmax)
    axes[1].set_title(f"{label2}")
    plt.colorbar(im1, ax=axes[1], shrink=0.8)

    im2 = axes[2].imshow(diff, origin='lower', cmap='RdBu_r')
    axes[2].set_title("Difference (file1 - file2)")
    plt.colorbar(im2, ax=axes[2], shrink=0.8)

    fig.suptitle(f"{var_name}  ({slice_label})", fontsize=13)
    fig.tight_layout()

    outpath = os.path.join(figures_dir, f"{var_name}.png")
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    print(f"    Saved figure: {outpath}")


def main():
    parser = argparse.ArgumentParser(description="Compare two HICAR output files")
    parser.add_argument("file1", help="First NetCDF file")
    parser.add_argument("file2", help="Second NetCDF file")
    parser.add_argument("--tolerance", type=float, default=0.0,
                        help="Maximum allowed absolute difference (0.0 = exact match)")
    parser.add_argument("--last-timestep-only", action="store_true",
                        help="Only compare the last timestep common to both files")
    parser.add_argument("--figures-dir", default=None,
                        help="Directory to save comparison plots for failing variables")
    args = parser.parse_args()

    try:
        ds1 = xr.open_dataset(args.file1)
        ds2 = xr.open_dataset(args.file2)
    except Exception as e:
        print(f"{bcolors.RED}Error opening files: {e}{bcolors.NC}")
        return 1

    # Find common data variables (skip coordinate-only vars)
    common_vars = sorted(set(ds1.data_vars) & set(ds2.data_vars))
    if not common_vars:
        print(f"{bcolors.RED}No common variables found between the two files{bcolors.NC}")
        return 1

    # Create figures directory if requested
    if args.figures_dir:
        os.makedirs(args.figures_dir, exist_ok=True)

    # Labels for plot titles (basename of each file's parent dir)
    label1 = os.path.basename(os.path.dirname(args.file1))
    label2 = os.path.basename(os.path.dirname(args.file2))

    print(f"Comparing {bcolors.BLUE}{args.file1}{bcolors.NC}")
    print(f"     with {bcolors.BLUE}{args.file2}{bcolors.NC}")
    if args.last_timestep_only:
        print("  (last timestep only)")
    print("-" * 60)

    error_flag = False

    for var_name in common_vars:
        v1 = ds1[var_name]
        v2 = ds2[var_name]

        if args.last_timestep_only and "time" in v1.dims:
            v1 = v1.isel(time=-1)
            v2 = v2.isel(time=-1)

        a1 = v1.values
        a2 = v2.values

        # Shapes may differ (e.g. different number of timesteps); trim to common size
        if a1.shape != a2.shape:
            common_shape = tuple(min(s1, s2) for s1, s2 in zip(a1.shape, a2.shape))
            slices = tuple(slice(0, s) for s in common_shape)
            a1 = a1[slices]
            a2 = a2[slices]

        diff = np.abs(a1 - a2)
        max_diff = float(np.nanmax(diff)) if diff.size > 0 else 0.0
        sum_diff = float(np.nansum(diff)) if diff.size > 0 else 0.0

        status = f"{bcolors.GREEN}OK{bcolors.NC}" if max_diff <= args.tolerance else f"{bcolors.RED}DIFF{bcolors.NC}"
        print(f"  {var_name:20s}: max|diff| = {max_diff:.6e} sum|diff| = {sum_diff:.6e}  [{status}]")

        if max_diff > args.tolerance:
            error_flag = True
            if args.figures_dir:
                plot_diff(a1, a2, var_name, label1, label2, args.figures_dir)

    print("-" * 60)

    ds1.close()
    ds2.close()

    if error_flag:
        print(f"{bcolors.RED}FAILED: Files differ beyond tolerance {args.tolerance}{bcolors.NC}")
        return 1
    else:
        print(f"{bcolors.GREEN}PASSED: All variables match within tolerance {args.tolerance}{bcolors.NC}")
        return 0


if __name__ == "__main__":
    sys.exit(main())
