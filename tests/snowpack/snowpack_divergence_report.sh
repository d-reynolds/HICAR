#!/usr/bin/env bash
# ===========================================================================
# SNOWPACK C++ vs Fortran divergence report.
#
# Companion to tests/snowpack/test_snowpack_compare.sh. When the nightly/PR comparison
# between the C++ SNOWPACK build and the native-Fortran port starts FAILING,
# the usual cause is new physics landing on the upstream `fortran-bindings`
# branch (C++ side) that has not yet been ported to the Fortran driver. This
# script turns that into a concrete porting to-do list:
#
#   1. Resolve the most recent HICAR commit blessed with
#      `hicar-snowpack-parity-blessed=success` and extract the upstream
#      SNOWPACK SHA recorded in its status description (snowpack=<sha>) —
#      the last upstream state at which the two implementations provably
#      matched.
#   2. Resolve the CURRENT upstream SNOWPACK SHA (from the fetched checkout,
#      or the remote branch head).
#   3. Print `git log` / `git diff --stat` between the two, split into the
#      C++ core (snowpack/, meteoio/) — changes that may need PORTING — and
#      fortran/ — changes to the port itself.
#
# Usage:
#   snowpack_divergence_report.sh <hicar_repo> [options]
#
#     --blessed-snowpack SHA  skip the gh lookup and use this anchor SHA
#     --snowpack-dir DIR      fetched SNOWPACK checkout
#                             (default: <hicar_repo>/build/external/SNOWPACK)
#     --current SHA|remote    compare the anchor against this SHA, or against
#                             the remote fortran-bindings head ('remote');
#                             default: HEAD of the fetched checkout
#     --full-diff FILE        also write the complete diff of the C++ core to FILE
#
# Requires: git; gh (authenticated) unless --blessed-snowpack is given.
# ===========================================================================
set -uo pipefail

BLUE='\033[0;36m'; GREEN='\033[0;32m'; RED='\033[0;31m'; NC='\033[0m'

hicar_repo="${1:?usage: snowpack_divergence_report.sh <hicar_repo> [options]}"
hicar_repo="$(cd "$hicar_repo" && pwd)"
shift

blessed_sha=""
snowpack_dir="$hicar_repo/build/external/SNOWPACK"
current="checkout"
full_diff=""
while [ $# -gt 0 ]; do
    case "$1" in
        --blessed-snowpack) blessed_sha="$2"; shift 2;;
        --snowpack-dir)     snowpack_dir="$2"; shift 2;;
        --current)          current="$2"; shift 2;;
        --full-diff)        full_diff="$2"; shift 2;;
        *) echo "Unknown option: $1"; exit 2;;
    esac
done

[ -d "$snowpack_dir/.git" ] || [ -f "$snowpack_dir/.git" ] || {
    echo -e "${RED}No SNOWPACK git checkout at ${snowpack_dir}.${NC}"
    echo "Configure a build first, or pass --snowpack-dir."
    exit 1
}

# --- 1. the parity anchor (last blessed upstream SHA) -----------------------
hicar_blessed=""
if [ -z "$blessed_sha" ]; then
    line=$(bash "$hicar_repo/tests/resolve_blessed_commit.sh" "$hicar_repo" \
            hicar-snowpack-parity-blessed --with-description) || {
        echo -e "${RED}No hicar-snowpack-parity-blessed commit found.${NC}"
        echo "Bless one first:  tests/snowpack/test_snowpack_compare.sh <repo> --bless"
        echo "or pass --blessed-snowpack <sha> directly."
        exit 1
    }
    hicar_blessed="${line%%$'\t'*}"
    desc="${line#*$'\t'}"
    blessed_sha=$(printf '%s' "$desc" | sed -n 's/.*snowpack=\([0-9a-f]\{7,40\}\).*/\1/p')
    [ -n "$blessed_sha" ] || {
        echo -e "${RED}Blessed commit ${hicar_blessed:0:12} has no snowpack=<sha> in its status description:${NC}"
        echo "  \"$desc\""
        exit 1
    }
fi

# --- 2. the current upstream SHA ---------------------------------------------
# The fetched checkout is shallow (GIT_SHALLOW), so deepen enough to contain the
# anchor before diffing.
git -C "$snowpack_dir" fetch --quiet origin fortran-bindings 2>/dev/null || true
if ! git -C "$snowpack_dir" cat-file -e "$blessed_sha" 2>/dev/null; then
    echo "Anchor not in the shallow clone — unshallowing..."
    git -C "$snowpack_dir" fetch --quiet --unshallow origin 2>/dev/null || \
        git -C "$snowpack_dir" fetch --quiet origin 2>/dev/null || true
fi
git -C "$snowpack_dir" cat-file -e "$blessed_sha" 2>/dev/null || {
    echo -e "${RED}Cannot find the anchor commit ${blessed_sha} in ${snowpack_dir} even after fetching.${NC}"
    exit 1
}

case "$current" in
    checkout) cur_sha=$(git -C "$snowpack_dir" rev-parse HEAD);;
    remote)   cur_sha=$(git -C "$snowpack_dir" rev-parse FETCH_HEAD 2>/dev/null) || \
              cur_sha=$(git -C "$snowpack_dir" ls-remote origin fortran-bindings | cut -f1);;
    *)        cur_sha="$current";;
esac

echo "======================================================="
echo -e "  ${BLUE}SNOWPACK upstream divergence report${NC}"
echo "======================================================="
[ -n "$hicar_blessed" ] && echo -e "  blessed HICAR commit:     ${hicar_blessed}"
echo -e "  parity anchor (SNOWPACK): ${GREEN}${blessed_sha}${NC}"
echo -e "  current  (SNOWPACK):      ${BLUE}${cur_sha}${NC}"
echo

if [ "$(git -C "$snowpack_dir" rev-parse "$blessed_sha")" = "$(git -C "$snowpack_dir" rev-parse "$cur_sha")" ]; then
    echo -e "${GREEN}No upstream movement since the last parity bless.${NC}"
    echo "If the comparison test fails anyway, the divergence is on the HICAR side"
    echo "(sm_SNOWPACK.F90 / sm_driver.F90 / seed or tolerance changes), not upstream."
    exit 0
fi

echo -e "${BLUE}--- Upstream C++ core changes since parity (PORTING CANDIDATES) ---${NC}"
git -C "$snowpack_dir" log --oneline --no-decorate "$blessed_sha..$cur_sha" -- snowpack/ meteoio/ | sed 's/^/  /'
echo
git -C "$snowpack_dir" diff --stat "$blessed_sha..$cur_sha" -- snowpack/ meteoio/ | tail -20 | sed 's/^/  /'
echo
echo -e "${BLUE}--- fortran/ (port-side) changes since parity ---${NC}"
git -C "$snowpack_dir" log --oneline --no-decorate "$blessed_sha..$cur_sha" -- fortran/ | sed 's/^/  /'
echo

if [ -n "$full_diff" ]; then
    git -C "$snowpack_dir" diff "$blessed_sha..$cur_sha" -- snowpack/ meteoio/ > "$full_diff"
    echo -e "Full C++ core diff written to ${BLUE}${full_diff}${NC}"
fi

echo "Next steps: port the relevant C++ changes into fortran/snowpack_driver.F90"
echo "(and friends), re-run tests/snowpack/test_snowpack_compare.sh until it passes, then"
echo "re-bless:  tests/snowpack/test_snowpack_compare.sh <repo> --bless --reason '...'"
