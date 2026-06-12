#!/usr/bin/env python3
"""Compare two HICAR NetCDF output files for regression / reproducibility / CPU-GPU testing.

This is the single comparison engine used by all HICAR test lanes:

  * exact      — bit-for-bit (CPU regression vs a CPU baseline on a fixed toolchain)
  * tolerance  — per-variable relative+absolute tolerances (CPU<->GPU, gnu<->nvfortran),
                 driven by a YAML spec (see tests/tolerances.yaml)

It is deliberately strict about the things that matter for "the model won't crash and
behaves as expected":

  * NaN / Inf introduced by FILE2 (the candidate) that are not present in FILE1 (the
    baseline) are ALWAYS a failure, regardless of tolerance. The previous implementation
    used nanmax/nansum and silently ignored NaNs — a run that blew up to NaN could pass.
  * A variable present in the baseline but missing from the candidate is a failure
    (a dropped/renamed output). A variable added by the candidate is reported as a warning.
  * On failure it localizes the worst offender (variable, dimension indices, both values,
    abs and rel error) so a developer learns *what* broke and *where*.

Usage:
    compare_outputs.py FILE1 FILE2 [options]

    FILE1 is treated as the reference / baseline, FILE2 as the candidate.

Common invocations:
    # bit-for-bit (legacy default; --tolerance 0.0 == exact)
    compare_outputs.py baseline.nc candidate.nc --mode exact

    # per-variable tolerance for the GPU lane
    compare_outputs.py cpu_ref.nc gpu.nc --tolerance-spec tests/tolerances.yaml \
        --github-annotations --report-json report.json
"""
import argparse
import json
import os
import sys

import numpy as np
import xarray as xr


class bcolors:
    BLUE = '\033[0;36m'
    GREEN = '\033[0;32m'
    YELLOW = '\033[0;33m'
    RED = '\033[0;31m'
    NC = '\033[0m'


# ----------------------------------------------------------------------------
# Tolerance spec
# ----------------------------------------------------------------------------
class ToleranceSpec:
    """Resolves the (rtol, atol, frac) triple to apply to a given variable.

    Spec file (YAML) format:

        defaults:
          rtol: 1.0e-5
          atol: 1.0e-8
        variables:
          pressure: {rtol: 1.0e-6, atol: 1.0e-2}
          qv:       {rtol: 1.0e-4, atol: 1.0e-12}
          hfls:     {rtol: 1.0e-3, atol: 5.0e-1, frac: 1.0e-4}
        ignore:
          - time

    `frac` (default 0) is the maximum allowed FRACTION of finite values that may
    violate atol+rtol*|baseline| before the variable is flagged DIFF. Use it for
    fields with rare, legitimate single-cell outliers (e.g. threshold/regime
    flips such as melt onset) where a max-based tolerance would have to be so
    loose it loses all discriminating power for the bulk of the field.
    """

    def __init__(self, default_rtol=0.0, default_atol=0.0, default_frac=0.0,
                 per_var=None, ignore=None):
        self.default_rtol = float(default_rtol)
        self.default_atol = float(default_atol)
        self.default_frac = float(default_frac)
        self.per_var = per_var or {}
        self.ignore = set(ignore or [])

    @classmethod
    def from_file(cls, path):
        import yaml  # local import so exact-mode runs need no pyyaml
        with open(path) as fh:
            spec = yaml.safe_load(fh) or {}
        defaults = spec.get("defaults", {}) or {}
        per_var = {}
        for name, tol in (spec.get("variables", {}) or {}).items():
            tol = tol or {}
            per_var[name] = (
                float(tol.get("rtol", defaults.get("rtol", 0.0))),
                float(tol.get("atol", defaults.get("atol", 0.0))),
                float(tol.get("frac", defaults.get("frac", 0.0))),
            )
        return cls(
            default_rtol=defaults.get("rtol", 0.0),
            default_atol=defaults.get("atol", 0.0),
            default_frac=defaults.get("frac", 0.0),
            per_var=per_var,
            ignore=spec.get("ignore", []) or [],
        )

    def for_var(self, name):
        return self.per_var.get(name,
                                (self.default_rtol, self.default_atol, self.default_frac))


# ----------------------------------------------------------------------------
# Per-variable comparison
# ----------------------------------------------------------------------------
def _index_label(flat_index, shape, dims):
    """Turn a flat index into a human-readable 'dim=i, dim=j, ...' label."""
    if not shape:
        return "(scalar)"
    idx = np.unravel_index(flat_index, shape)
    return ", ".join(f"{d}={int(i)}" for d, i in zip(dims, idx))


def compare_variable(name, v1, v2, rtol, atol, frac=0.0):
    """Compare one variable. Returns a result dict with status in
    {'ok', 'diff', 'nan', 'shape', 'skip'} and supporting diagnostics.
    Up to frac*n_finite violating values are tolerated (rare-outlier allowance,
    see ToleranceSpec); the tolerated count is reported in the result."""
    dims = list(v1.dims)
    a1 = np.asarray(v1.values)
    a2 = np.asarray(v2.values)

    # Non-numeric data (characters, etc.) — only exact-equality is meaningful.
    if not (np.issubdtype(a1.dtype, np.number) and np.issubdtype(a2.dtype, np.number)):
        equal = a1.shape == a2.shape and bool(np.all(a1 == a2))
        return {"name": name, "status": "ok" if equal else "diff",
                "message": "non-numeric", "rtol": rtol, "atol": atol}

    a1 = a1.astype(np.float64)
    a2 = a2.astype(np.float64)

    if a1.shape != a2.shape:
        return {"name": name, "status": "shape", "rtol": rtol, "atol": atol,
                "shape1": list(a1.shape), "shape2": list(a2.shape)}

    # --- NaN / Inf: introduced by the candidate is always a failure ----------
    bad1 = ~np.isfinite(a1)
    bad2 = ~np.isfinite(a2)
    introduced = bad2 & ~bad1
    n_introduced = int(np.count_nonzero(introduced))
    if n_introduced > 0:
        first = int(np.argmax(introduced.ravel()))
        return {"name": name, "status": "nan", "rtol": rtol, "atol": atol,
                "n_introduced": n_introduced,
                "where": _index_label(first, a1.shape, dims),
                "baseline_value": float(a1.ravel()[first])}

    # Compare only positions finite in both (pre-existing fills are ignored).
    finite = np.isfinite(a1) & np.isfinite(a2)
    if not np.any(finite):
        return {"name": name, "status": "ok", "rtol": rtol, "atol": atol,
                "max_abs": 0.0, "max_rel": 0.0, "n_violations": 0,
                "message": "no finite values to compare"}

    abs_err = np.abs(a2 - a1)
    abs_err = np.where(finite, abs_err, 0.0)

    denom = np.abs(a1)
    rel_err = np.where(finite & (denom > 0), abs_err / np.where(denom > 0, denom, 1.0), 0.0)

    max_abs = float(np.max(abs_err))
    max_rel = float(np.max(rel_err))

    # isclose semantics: violation where |a2-a1| > atol + rtol*|a1|
    threshold = atol + rtol * denom
    violations = finite & (abs_err > threshold)
    n_violations = int(np.count_nonzero(violations))

    n_finite = int(np.count_nonzero(finite))
    result = {"name": name, "rtol": rtol, "atol": atol, "frac": frac,
              "max_abs": max_abs, "max_rel": max_rel,
              "n_violations": n_violations,
              "n_finite": n_finite}

    if n_violations <= int(frac * n_finite):
        result["status"] = "ok"
        return result

    # Localize the worst offender by how far it exceeds its own threshold.
    margin = np.where(violations, abs_err - threshold, -np.inf)
    worst = int(np.argmax(margin.ravel()))
    result["status"] = "diff"
    result["where"] = _index_label(worst, a1.shape, dims)
    result["baseline_value"] = float(a1.ravel()[worst])
    result["candidate_value"] = float(a2.ravel()[worst])
    result["worst_abs"] = float(abs_err.ravel()[worst])
    result["worst_rel"] = float(rel_err.ravel()[worst])
    return result


# ----------------------------------------------------------------------------
# Plotting (unchanged behaviour, used for failing variables)
# ----------------------------------------------------------------------------
def plot_diff(a1, a2, var_name, label1, label2, figures_dir):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print(f"  {bcolors.RED}matplotlib not available, skipping figure for {var_name}{bcolors.NC}")
        return

    a1 = np.asarray(a1, dtype=np.float64)
    a2 = np.asarray(a2, dtype=np.float64)
    diff = a1 - a2

    if diff.ndim >= 4:
        # Physically-3D field (time, layer, y, x): show the mean across all
        # layers at the last timestep, rather than a single mid-level slice
        # (which can be ~0 for fields whose signal lives in other layers).
        a1 = np.nanmean(a1[-1], axis=0)
        a2 = np.nanmean(a2[-1], axis=0)
        diff = np.nanmean(diff[-1], axis=0)
        slice_label = "t=-1, layer-mean"
    elif diff.ndim == 3:
        # 2D surface field (time, y, x): last timestep on all three panels.
        a1 = a1[-1, :, :]
        a2 = a2[-1, :, :]
        diff = diff[-1, :, :]
        slice_label = "t=-1"
    elif diff.ndim == 2:
        slice_label = "2D field"
    else:
        return

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    vmin = min(np.nanmin(a1), np.nanmin(a2))
    vmax = max(np.nanmax(a1), np.nanmax(a2))
    im0 = axes[0].imshow(a1, origin='lower', vmin=vmin, vmax=vmax)
    axes[0].set_title(label1); plt.colorbar(im0, ax=axes[0], shrink=0.8)
    im1 = axes[1].imshow(a2, origin='lower', vmin=vmin, vmax=vmax)
    axes[1].set_title(label2); plt.colorbar(im1, ax=axes[1], shrink=0.8)
    im2 = axes[2].imshow(diff, origin='lower', cmap='RdBu_r')
    axes[2].set_title("Difference (file1 - file2)"); plt.colorbar(im2, ax=axes[2], shrink=0.8)
    fig.suptitle(f"{var_name}  ({slice_label})", fontsize=13)
    fig.tight_layout()
    outpath = os.path.join(figures_dir, f"{var_name}.png")
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    print(f"    Saved figure: {outpath}")


# ----------------------------------------------------------------------------
def gh_annotation(level, title, message):
    """Emit a GitHub Actions workflow annotation (rendered in the PR checks UI)."""
    msg = message.replace('\n', '%0A')
    print(f"::{level} title={title}::{msg}")


def main():
    p = argparse.ArgumentParser(description="Compare two HICAR output files")
    p.add_argument("file1", help="Reference / baseline NetCDF file")
    p.add_argument("file2", help="Candidate NetCDF file")
    p.add_argument("--mode", choices=["exact", "tolerance"], default=None,
                   help="exact = bit-for-bit; tolerance = use spec/flags. "
                        "Default: exact unless a spec or nonzero --tolerance is given.")
    p.add_argument("--tolerance", type=float, default=None,
                   help="Legacy global absolute tolerance (atol). 0.0 = exact match.")
    p.add_argument("--tolerance-spec", default=None,
                   help="YAML file of per-variable rtol/atol (see tests/tolerances.yaml)")
    p.add_argument("--default-rtol", type=float, default=0.0)
    p.add_argument("--default-atol", type=float, default=0.0)
    p.add_argument("--last-timestep-only", action="store_true")
    p.add_argument("--allow-trim", action="store_true",
                   help="Trim mismatched array shapes to the common size instead of failing "
                        "(needed for continuous-vs-restart reproducibility comparisons).")
    p.add_argument("--ignore-added", action="store_true",
                   help="Do not warn about variables present only in the candidate.")
    p.add_argument("--figures-dir", default=None)
    p.add_argument("--report-json", default=None,
                   help="Write a machine-readable JSON report to this path.")
    p.add_argument("--github-annotations", action="store_true",
                   help="Emit ::error:: / ::warning:: annotations for the PR checks UI.")
    args = p.parse_args()

    # Resolve the tolerance spec.
    if args.tolerance_spec:
        spec = ToleranceSpec.from_file(args.tolerance_spec)
        mode = args.mode or "tolerance"
    else:
        # Legacy --tolerance maps onto a global atol.
        atol = args.tolerance if args.tolerance is not None else args.default_atol
        rtol = args.default_rtol
        spec = ToleranceSpec(default_rtol=rtol, default_atol=atol)
        if args.mode == "exact":
            spec = ToleranceSpec(0.0, 0.0)
        mode = args.mode or ("tolerance" if (atol > 0 or rtol > 0) else "exact")

    try:
        ds1 = xr.open_dataset(args.file1)
        ds2 = xr.open_dataset(args.file2)
    except Exception as e:
        print(f"{bcolors.RED}Error opening files: {e}{bcolors.NC}")
        return 1

    vars1 = set(ds1.data_vars)
    vars2 = set(ds2.data_vars)
    missing = sorted((vars1 - vars2) - spec.ignore)   # in baseline, gone from candidate
    added = sorted((vars2 - vars1) - spec.ignore)     # new in candidate
    common = sorted((vars1 & vars2) - spec.ignore)

    if args.figures_dir:
        os.makedirs(args.figures_dir, exist_ok=True)
    label1 = os.path.basename(os.path.dirname(args.file1)) or "file1"
    label2 = os.path.basename(os.path.dirname(args.file2)) or "file2"

    print(f"Comparing {bcolors.BLUE}{args.file1}{bcolors.NC}  (baseline)")
    print(f"     with {bcolors.BLUE}{args.file2}{bcolors.NC}  (candidate)")
    print(f"  mode = {mode}" + (f"   spec = {args.tolerance_spec}" if args.tolerance_spec else ""))
    if args.last_timestep_only:
        print("  (last timestep only)")
    print("-" * 78)

    results = []
    failed = False

    # Dropped variables are regressions.
    for name in missing:
        failed = True
        msg = f"variable '{name}' present in baseline but MISSING from candidate"
        print(f"  {name:24s}: {bcolors.RED}MISSING{bcolors.NC}")
        results.append({"name": name, "status": "missing"})
        if args.github_annotations:
            gh_annotation("error", "Dropped output variable", msg)

    for name in added:
        results.append({"name": name, "status": "added"})
        if not args.ignore_added:
            print(f"  {name:24s}: {bcolors.YELLOW}ADDED (candidate only){bcolors.NC}")
            if args.github_annotations:
                gh_annotation("warning", "New output variable",
                              f"variable '{name}' present only in candidate")

    for name in common:
        v1 = ds1[name]
        v2 = ds2[name]
        if args.last_timestep_only and "time" in v1.dims:
            v1 = v1.isel(time=-1)
            v2 = v2.isel(time=-1)
        if args.allow_trim and v1.shape != v2.shape:
            common_shape = tuple(min(s1, s2) for s1, s2 in zip(v1.shape, v2.shape))
            slc = tuple(slice(0, s) for s in common_shape)
            v1 = v1[slc]
            v2 = v2[slc]

        rtol, atol, frac = spec.for_var(name)
        r = compare_variable(name, v1, v2, rtol, atol, frac)
        results.append(r)

        st = r["status"]
        if st == "ok":
            extra = (f"max|abs|={r['max_abs']:.3e} max|rel|={r['max_rel']:.3e}"
                     if "max_abs" in r else r.get("message", ""))
            if r.get("n_violations", 0) > 0:
                extra += (f"  [{r['n_violations']}/{r['n_finite']} outliers "
                          f"tolerated, frac<={r['frac']:.1e}]")
            print(f"  {name:24s}: {bcolors.GREEN}OK{bcolors.NC}   {extra}")
        elif st == "nan":
            failed = True
            msg = (f"{r['n_introduced']} new NaN/Inf value(s); first at {r['where']} "
                   f"(baseline was {r['baseline_value']:.6e})")
            print(f"  {name:24s}: {bcolors.RED}NaN/Inf{bcolors.NC}  {msg}")
            if args.github_annotations:
                gh_annotation("error", f"NaN/Inf introduced in {name}", msg)
            if args.figures_dir:
                plot_diff(ds1[name].values, ds2[name].values, name, label1, label2, args.figures_dir)
        elif st == "shape":
            failed = True
            msg = f"shape {r['shape1']} (baseline) != {r['shape2']} (candidate)"
            print(f"  {name:24s}: {bcolors.RED}SHAPE{bcolors.NC}   {msg}")
            if args.github_annotations:
                gh_annotation("error", f"Shape mismatch in {name}", msg)
        elif st == "diff":
            failed = True
            msg = (f"{r['n_violations']}/{r['n_finite']} values exceed tol "
                   f"(rtol={rtol:.1e}, atol={atol:.1e}); worst at {r['where']}: "
                   f"baseline={r['baseline_value']:.6e} candidate={r['candidate_value']:.6e} "
                   f"(abs={r['worst_abs']:.3e}, rel={r['worst_rel']:.3e})")
            print(f"  {name:24s}: {bcolors.RED}DIFF{bcolors.NC}    {msg}")
            if args.github_annotations:
                gh_annotation("error", f"{name} differs beyond tolerance", msg)
            if args.figures_dir:
                plot_diff(ds1[name].values, ds2[name].values, name, label1, label2, args.figures_dir)

    print("-" * 78)
    ds1.close()
    ds2.close()

    n_fail = sum(1 for r in results if r["status"] in ("diff", "nan", "shape", "missing"))
    summary = {"passed": not failed, "mode": mode,
               "n_compared": len(common), "n_failed": n_fail,
               "n_missing": len(missing), "n_added": len(added),
               "results": results}
    if args.report_json:
        with open(args.report_json, "w") as fh:
            json.dump(summary, fh, indent=2)
        print(f"  JSON report: {bcolors.BLUE}{args.report_json}{bcolors.NC}")

    if failed:
        print(f"{bcolors.RED}FAILED: {n_fail} variable(s) differ beyond tolerance / "
              f"introduced NaN / changed shape / were dropped.{bcolors.NC}")
        return 1
    print(f"{bcolors.GREEN}PASSED: all {len(common)} compared variable(s) within tolerance.{bcolors.NC}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
