#!/usr/bin/env python3
"""Compare two SNOWPACK-parity comparison reports (parity_report.json, written
by tests/scripts/snowpack/test_snowpack_compare.sh via compare_outputs.py --report-json and
archived on each `snowpack-parity/*` GitHub release).

Shows how the C++-vs-Fortran residual changed between two blesses:

    # fetch two archived reports
    gh release download snowpack-parity/20260612-abc1234 -p parity_report.json -O old.json
    gh release download snowpack-parity/20260701-def5678 -p parity_report.json -O new.json
    tests/scripts/snowpack/parity_trend.py old.json new.json

Variables are ranked by the growth of max|abs| difference; growth beyond
--flag-ratio (default 3x) is marked as a potential regression worth a
tests/scripts/snowpack/snowpack_divergence_report.sh look.
"""
import argparse
import json


def load(path):
    with open(path) as fh:
        rep = json.load(fh)
    return {r["name"]: r for r in rep.get("results", []) if "max_abs" in r}


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("old", help="earlier parity_report.json")
    p.add_argument("new", help="later parity_report.json")
    p.add_argument("--flag-ratio", type=float, default=3.0,
                   help="flag variables whose max|abs| grew by more than this factor")
    p.add_argument("--top", type=int, default=25, help="rows to print")
    args = p.parse_args()

    old, new = load(args.old), load(args.new)
    common = sorted(set(old) & set(new))
    rows = []
    for nm in common:
        a, b = old[nm]["max_abs"], new[nm]["max_abs"]
        ratio = (b / a) if a > 0 else (float("inf") if b > 0 else 1.0)
        rows.append((ratio, nm, a, b,
                     old[nm].get("n_violations", 0), new[nm].get("n_violations", 0)))
    rows.sort(reverse=True)

    print(f"{len(common)} common variables "
          f"(old-only: {len(set(old)-set(new))}, new-only: {len(set(new)-set(old))})")
    print(f"{'variable':>26} {'max|abs| old':>13} {'max|abs| new':>13} "
          f"{'ratio':>8} {'viol old':>9} {'viol new':>9}")
    flagged = 0
    for ratio, nm, a, b, vo, vn in rows[:args.top]:
        mark = "  <-- REGRESSED?" if ratio > args.flag_ratio and b > 0 else ""
        if mark:
            flagged += 1
        rs = f"{ratio:8.2f}" if ratio != float("inf") else "     inf"
        print(f"{nm:>26} {a:13.3e} {b:13.3e} {rs} {vo:9d} {vn:9d}{mark}")

    if flagged:
        print(f"\n{flagged} variable(s) grew >{args.flag_ratio}x — check upstream drift:")
        print("  tests/scripts/snowpack/snowpack_divergence_report.sh <hicar_repo>")
    else:
        print(f"\nNo variable grew >{args.flag_ratio}x — parity level stable.")


if __name__ == "__main__":
    main()
