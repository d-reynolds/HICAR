#!/bin/bash
# =============================================================================
# HICAR regression test — compare integration outputs against a "blessed" commit.
#
# The blessed reference is NOT a stored .nc file and NOT a file in the repo: it is
# the most recent commit carrying a `hicar-regression-blessed=success` commit
# status (resolved by tests/resolve_blessed_commit.sh). This script:
#   1. builds the blessed commit's HICAR exe (cached by hash, in a git worktree),
#   2. runs it on the requested integration case(s) to REGENERATE the reference,
#   3. diffs the *current* integration outputs (already produced by
#      test_case_runner.sh) against the regenerated reference with compare_outputs.py.
#
# Because both the current and the blessed exe are built on the same runner with
# the same toolchain, the comparison defaults to bit-for-bit (--mode exact): any
# difference is a real change in model output.
#
# Run integration FIRST (so output/<case>/ exists), then this:
#   test_regression.sh <hicar_repo> <build_dir> <cases> [options]
#     <cases>             comma-separated base cases, e.g. "Standard,Nested"
#                         (_restart variants are reproducibility cases; skipped)
#     --blessed-commit SHA  use this commit as the reference (else auto-resolve the
#                           most recent hicar-regression-blessed ancestor via gh)
#     --mode exact|tolerance   default exact
#     --tolerance-spec PATH     spec for tolerance mode (default tests/tolerances.yaml)
#     --bless              instead of comparing, post hicar-regression-blessed=success
#                          to HEAD (needs gh + statuses:write)
# =============================================================================
set -euo pipefail
BLUE='\033[0;36m'; GREEN='\033[0;32m'; RED='\033[0;31m'; YELLOW='\033[0;33m'; NC='\033[0m'

if [ $# -lt 3 ]; then
    echo "Usage: $0 <hicar_repo> <build_dir> <cases> [--blessed-commit SHA] [--mode exact|tolerance] [--tolerance-spec PATH] [--bless]"
    exit 1
fi
hicar_repo=$(cd "$1" && pwd); shift
BUILD_DIR=$(cd "$1" && pwd); shift
CASES="$1"; shift
MODE="exact"
TOL_SPEC="$hicar_repo/tests/tolerances.yaml"
BLESS=false
BLESSED_COMMIT=""
while [ $# -gt 0 ]; do
    case "$1" in
        --mode) MODE="$2"; shift 2;;
        --tolerance-spec) TOL_SPEC="$2"; shift 2;;
        --blessed-commit) BLESSED_COMMIT="$2"; shift 2;;
        --bless) BLESS=true; shift;;
        *) echo -e "${RED}Unknown arg: $1${NC}"; exit 1;;
    esac
done

COMPARE="$hicar_repo/tests/compare_outputs.py"
TC="$hicar_repo/tests/Test_Cases"

# The regression reference is the most recent commit carrying a
# `hicar-regression-blessed=success` commit status (see resolve_blessed_commit.sh
# and bless-baseline.yml) — there is no stored reference file.

# --- bless mode: post the regression-blessed status to HEAD, done -----------
if [ "$BLESS" = true ]; then
    head_sha=$(git -C "$hicar_repo" rev-parse HEAD)
    slug="${GITHUB_REPOSITORY:-$(gh repo view --json nameWithOwner -q .nameWithOwner 2>/dev/null || true)}"
    [ -n "$slug" ] || { echo -e "${RED}Cannot determine owner/repo for gh${NC}"; exit 1; }
    gh api -X POST "/repos/$slug/statuses/$head_sha" \
        -f state=success -f context=hicar-regression-blessed \
        -f description="Blessed as the regression reference" >/dev/null
    echo -e "${GREEN}Blessed ${head_sha} (hicar-regression-blessed=success)${NC}"
    exit 0
fi

# --- resolve the blessed reference commit -----------------------------------
REF_HASH="$BLESSED_COMMIT"
if [ -z "$REF_HASH" ]; then
    REF_HASH=$(bash "$hicar_repo/tests/resolve_blessed_commit.sh" "$hicar_repo" --exclude-head) || {
        echo -e "${RED}Could not resolve a blessed commit (none found, or gh not authed).${NC}"
        echo -e "${RED}Bless one first (bless-baseline.yml or test_regression.sh --bless), or pass --blessed-commit.${NC}"
        exit 1
    }
fi
git -C "$hicar_repo" cat-file -t "$REF_HASH" &>/dev/null || \
    git -C "$hicar_repo" fetch origin "$REF_HASH" 2>/dev/null || true
git -C "$hicar_repo" cat-file -t "$REF_HASH" &>/dev/null || \
    { echo -e "${RED}Blessed commit '$REF_HASH' is not a valid commit${NC}"; exit 1; }
if [ "$(git -C "$hicar_repo" rev-parse "$REF_HASH")" = "$(git -C "$hicar_repo" rev-parse HEAD)" ]; then
    echo -e "${RED}Blessed commit equals HEAD — regression would trivially pass.${NC}"; exit 1
fi

echo "======================================================="
echo -e "  HICAR regression vs blessed commit ${BLUE}${REF_HASH:0:12}${NC}"
echo -e "  Cases: ${BLUE}${CASES}${NC}   Mode: ${BLUE}${MODE}${NC}"
echo "======================================================="

# --- build (or reuse cached) blessed exe ------------------------------------
OLD_EXE="$hicar_repo/bin/HICAR_blessed"
OLD_HASH_FILE="$hicar_repo/bin/.regression_hash"
WORKTREE_DIR=""
cleanup() { [ -n "$WORKTREE_DIR" ] && [ -d "$WORKTREE_DIR" ] && \
    git -C "$hicar_repo" worktree remove "$WORKTREE_DIR" --force 2>/dev/null || true; }
trap cleanup EXIT ERR INT TERM

if [ -f "$OLD_EXE" ] && [ "$(head -n1 "$OLD_HASH_FILE" 2>/dev/null | tr -d '[:space:]')" = "$REF_HASH" ]; then
    echo -e "${GREEN}Cached blessed exe matches ${REF_HASH:0:12} — skipping rebuild${NC}"
else
    echo -e "${BLUE}Building blessed exe from ${REF_HASH:0:12}...${NC}"
    WORKTREE_DIR=$(mktemp -d "${TMPDIR:-/tmp}/hicar_regression_XXXXXX"); rmdir "$WORKTREE_DIR"
    git -C "$hicar_repo" worktree add "$WORKTREE_DIR" "$REF_HASH" --detach
    # Reproduce the current build's CMake config so the comparison is apples-to-apples.
    cmake_args=()
    if [ -f "$BUILD_DIR/CMakeCache.txt" ]; then
        for var in MODE FC OPENACC FSM SNOWPACK NETCDF_DIR FFTW_DIR PETSC_DIR NCCL; do
            val=$(grep "^${var}:" "$BUILD_DIR/CMakeCache.txt" 2>/dev/null | head -n1 | sed 's/^[^=]*=//')
            [ -n "$val" ] && cmake_args+=("-D${var}=${val}")
        done
    fi
    OLD_BUILD="$WORKTREE_DIR/build_regression"; mkdir -p "$OLD_BUILD"
    cmake -S "$WORKTREE_DIR" -B "$OLD_BUILD" "${cmake_args[@]}" 2>&1 | tail -n 5
    cmake --build "$OLD_BUILD" --target HICAR --parallel 2>&1 | tail -n 5
    cmake --install "$OLD_BUILD" 2>&1 | tail -n 3 || true
    old_built=""
    for n in HICAR_debug_gpu HICAR_debug HICAR_gpu HICAR; do
        [ -f "$WORKTREE_DIR/bin/$n" ] && { old_built="$WORKTREE_DIR/bin/$n"; break; }
    done
    [ -n "$old_built" ] || { echo -e "${RED}No exe built from blessed commit${NC}"; exit 1; }
    mkdir -p "$hicar_repo/bin"; cp "$old_built" "$OLD_EXE"; echo "$REF_HASH" > "$OLD_HASH_FILE"
    git -C "$hicar_repo" worktree remove "$WORKTREE_DIR" --force 2>/dev/null || true; WORKTREE_DIR=""
fi
echo -e "  Blessed exe: ${BLUE}${OLD_EXE}${NC}"

# --- MPI launcher / rank count ----------------------------------------------
mpiexec_path=""
IFS=':' read -ra PD <<< "$PATH"
for d in "${PD[@]}"; do
    [ -x "$d/mpiexec" ] && ! [[ "$d" =~ python|conda ]] && { mpiexec_path="$d/mpiexec"; break; }
done
[ -z "$mpiexec_path" ] && { echo -e "${RED}mpiexec not found${NC}"; exit 1; }
if [[ "$OSTYPE" == "darwin"* ]]; then total_np=$(sysctl -n hw.logicalcpu); else total_np=$(nproc --all); fi
np=$((total_np / 2)); np=$((np > 2 ? np : 2)); np=$((np < 21 ? np : 21))
export OMP_NUM_THREADS=1

cd "$TC"
[ ! -f ./input/file_list_TestCase.txt ] && \
    "$hicar_repo/helpers/filelist_script.sh" "forcing/*" input/file_list_TestCase.txt

python_exe=$(command -v python3 || command -v python)
overall_status=0
IFS=',' read -ra CASE_ARR <<< "$CASES"
for raw in "${CASE_ARR[@]}"; do
    case_name=$(echo "$raw" | xargs)
    # _restart variants are reproducibility cases (run-to-restart), not regression.
    [[ "$case_name" == *"_restart"* ]] && { echo -e "${YELLOW}Skipping restart case ${case_name} (reproducibility, not regression)${NC}"; continue; }
    gen="input/nml_gen_scripts/${case_name}.sh"
    [ -f "$gen" ] || { echo -e "${YELLOW}No gen script for ${case_name}; skipping${NC}"; continue; }

    # The integration output for this case must already exist (run integration first).
    cur_dir="output/${case_name}"
    cur_nc=$(ls "${cur_dir}"/*.nc 2>/dev/null | head -n1 || true)
    [ -n "$cur_nc" ] || { echo -e "${RED}No integration output in ${cur_dir} — run test_case_runner.sh first${NC}"; overall_status=1; continue; }
    out_name=$(basename "$cur_nc")

    # Regenerate the blessed reference for this case.
    bless_dir="output/Blessed_${case_name}"
    rm -rf "$bless_dir" "restart/Blessed_${case_name}"
    mkdir -p "$bless_dir" "restart/Blessed_${case_name}"
    ( cd input
      nml="Blessed_${case_name}.nml"
      "$OLD_EXE" --gen-nml "$nml"
      bash "nml_gen_scripts/${case_name}.sh" "$nml"
      sed -i'.bak' "s|output_folder = .*|output_folder = '../${bless_dir}/'|" "$nml"
      sed -i'.bak' "s|restart_folder = .*|restart_folder = '../restart/Blessed_${case_name}/'|" "$nml"
      rm -f "${nml}.bak"
      echo -e "Regenerating blessed ${BLUE}${case_name}${NC} with ${np} ranks..."
      "$mpiexec_path" -np "$np" "$OLD_EXE" "$nml" 1>"Blessed_${case_name}.out" 2>"Blessed_${case_name}.err" || {
          echo -e "${RED}Blessed run failed for ${case_name}; last stderr:${NC}"; tail -n 20 "Blessed_${case_name}.err"; exit 1; }
    )
    bless_nc="${bless_dir}/${out_name}"
    [ -f "$bless_nc" ] || { echo -e "${RED}Blessed output not produced for ${case_name} (${bless_nc})${NC}"; overall_status=1; continue; }

    echo "-------------------------------------------------------"
    cmp_args=(--github-annotations --report-json "$hicar_repo/tests/regression_report_${case_name}.json" \
              --figures-dir "$hicar_repo/tests/figures/regression/${case_name}" --allow-trim)
    [ "$MODE" = "tolerance" ] && cmp_args+=(--tolerance-spec "$TOL_SPEC") || cmp_args+=(--mode exact)
    if "$python_exe" "$COMPARE" "$bless_nc" "$cur_nc" "${cmp_args[@]}"; then
        echo -e "${GREEN}Regression PASSED for ${case_name}${NC}"
    else
        echo -e "${RED}Regression FAILED for ${case_name}${NC}"; overall_status=1
    fi
done

exit $overall_status
