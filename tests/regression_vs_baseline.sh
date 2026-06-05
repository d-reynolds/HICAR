#!/bin/bash
# =============================================================================
# HICAR regression against a STORED baseline (replaces the rebuild-old-binary
# approach in test_regression.sh).
#
# Runs the current HICAR executable on the Standard test case with
# output_vars='all', then either:
#   * compares the output against a stored baseline NetCDF (default), or
#   * --bless: copies the output to the baseline path (to refresh the baseline).
#
# Baselines live in the Test-Data repo (git-LFS) under baselines/<lane>/ and are
# checked out alongside the test data at tests/Test_Cases/baselines/<lane>/.
# The bless-baseline.yml workflow regenerates them and opens a PR to Test-Data,
# giving an explicit audit trail of intentional, result-changing commits.
#
# Usage:
#   regression_vs_baseline.sh <hicar_repo> <build_dir> [options]
#     --lane NAME          baseline lane subdir (default: cpu-gnu)
#     --bless              write the output as the new baseline instead of comparing
#     --mode exact|tolerance   comparison mode (default: exact)
#     --tolerance-spec PATH    YAML spec for tolerance mode (default: tests/tolerances.yaml)
#     --baseline-root DIR  override baseline location
#                          (default: <repo>/tests/Test_Cases/baselines)
# =============================================================================
set -euo pipefail

BLUE='\033[0;36m'; GREEN='\033[0;32m'; RED='\033[0;31m'; YELLOW='\033[0;33m'; NC='\033[0m'

if [ $# -lt 2 ]; then
    echo "Usage: $0 <hicar_repo> <build_dir> [--lane NAME] [--bless] [--mode exact|tolerance] [--tolerance-spec PATH] [--baseline-root DIR]"
    exit 1
fi
hicar_repo=$(cd "$1" && pwd); shift
BUILD_DIR=$(cd "$1" && pwd); shift

LANE="cpu-gnu"
BLESS=false
MODE="exact"
TOL_SPEC="$hicar_repo/tests/tolerances.yaml"
BASELINE_ROOT="$hicar_repo/tests/Test_Cases/baselines"
while [ $# -gt 0 ]; do
    case "$1" in
        --lane) LANE="$2"; shift 2;;
        --bless) BLESS=true; shift;;
        --mode) MODE="$2"; shift 2;;
        --tolerance-spec) TOL_SPEC="$2"; shift 2;;
        --baseline-root) BASELINE_ROOT="$2"; shift 2;;
        *) echo -e "${RED}Unknown arg: $1${NC}"; exit 1;;
    esac
done

COMPARE_SCRIPT="$hicar_repo/tests/compare_outputs.py"
STANDARD_GEN="$hicar_repo/tests/Test_Cases/input/nml_gen_scripts/Standard.sh"
FIGURES_DIR="$hicar_repo/tests/figures/regression"
OUTPUT_FILENAME="Gaudergrat_250m_2017-02-14_00-00-00.nc"
BASELINE_FILE="$BASELINE_ROOT/$LANE/$OUTPUT_FILENAME"

# --- Locate the current executable -------------------------------------------
new_exe=""
for name in HICAR_debug_gpu HICAR_debug HICAR HICAR_gpu; do
    if [ -f "$hicar_repo/bin/$name" ]; then new_exe="$hicar_repo/bin/$name"; break; fi
done
[ -z "$new_exe" ] && { echo -e "${RED}No HICAR executable in $hicar_repo/bin${NC}"; exit 1; }

use_gpu=false
[[ "$new_exe" == *"gpu"* ]] && use_gpu=true

echo "======================================================="
echo -e "  HICAR Regression vs stored baseline"
echo -e "  Executable: ${BLUE}${new_exe}${NC}"
echo -e "  Lane:       ${BLUE}${LANE}${NC}"
echo -e "  Mode:       ${BLUE}${MODE}${NC}   $([ "$BLESS" = true ] && echo '(BLESS)')"
echo -e "  Baseline:   ${BLUE}${BASELINE_FILE}${NC}"
echo "======================================================="

# --- Support files, forcing list, GPU radiation data -------------------------
cd "$hicar_repo/tests/Test_Cases"
[ ! -f ./input/VEGPARM.TBL ]      && cp "$hicar_repo/run/"*.TBL ./input/
[ ! -d ./input/rrtmg_support ]    && cp -r "$hicar_repo/run/rrtmg_support" ./input
[ ! -d ./input/mp_support ]       && cp -r "$hicar_repo/run/mp_support" ./input
if [ "$use_gpu" = true ] && [ ! -d ./input/rrtmgp_support ]; then
    cp -r "$hicar_repo/run/rrtmgp_support" ./input
fi
[ ! -f ./input/file_list_TestCase.txt ] && \
    "$hicar_repo/helpers/filelist_script.sh" "forcing/*" input/file_list_TestCase.txt

# --- Generate the Standard namelist (canonical) + output_vars='all' override -
out_dir="../output/Regression_${LANE}"
restart_dir="../restart/Regression_${LANE}"
rm -rf "output/Regression_${LANE}" "restart/Regression_${LANE}"
mkdir -p "output/Regression_${LANE}" "restart/Regression_${LANE}"

cd input
nml="Regression_${LANE}.nml"
"$new_exe" --gen-nml "$nml"
# Apply the canonical Standard config, then override for regression.
bash "$STANDARD_GEN" "$nml"
sed -i'.bak' "s|output_vars = .*|output_vars = 'all'|" "$nml"
sed -i'.bak' "s|output_folder = '../output/Standard/'|output_folder = '${out_dir}/'|" "$nml"
sed -i'.bak' "s|restart_folder = '../restart/Standard/'|restart_folder = '${restart_dir}/'|" "$nml"
[ "$use_gpu" = true ] && sed -i'.bak' "s/rad = 'rrtmg'/rad = 'rrtmgp'/g" "$nml"
rm -f "${nml}.bak"

# --- MPI launcher + rank count (mirrors test_case_runner.sh) ------------------
mpiexec_path=""
IFS=':' read -ra PATH_DIRS <<< "$PATH"
for dir in "${PATH_DIRS[@]}"; do
    if [ -x "$dir/mpiexec" ] && ! [[ "$dir" =~ python|conda ]]; then mpiexec_path="$dir/mpiexec"; break; fi
done
if [ "$use_gpu" = true ]; then
    num_gpus=$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | wc -l | tr -d ' ')
    [ "${num_gpus:-0}" -lt 1 ] && { echo -e "${RED}GPU mode but no GPUs detected${NC}"; exit 1; }
    np=$((num_gpus + 1))
else
    if [[ "$OSTYPE" == "darwin"* ]]; then total_np=$(sysctl -n hw.logicalcpu); else total_np=$(nproc --all); fi
    np=$((total_np / 2)); np=$((np > 2 ? np : 2)); np=$((np < 21 ? np : 21))
fi
export OMP_NUM_THREADS=1

echo -e "Running Standard (output_vars=all) with ${np} ranks..."
if [ -n "$mpiexec_path" ]; then
    "$mpiexec_path" -np "$np" "$new_exe" "$nml" 1>"Regression_${LANE}.out" 2>"Regression_${LANE}.err" || {
        echo -e "${RED}HICAR run failed. Last stderr:${NC}"; tail -n 20 "Regression_${LANE}.err"; exit 1; }
elif command -v srun &>/dev/null; then
    srun -N 1 -n "$np" "$new_exe" "$nml" 1>"Regression_${LANE}.out" 2>"Regression_${LANE}.err" || {
        echo -e "${RED}HICAR run failed. Last stderr:${NC}"; tail -n 20 "Regression_${LANE}.err"; exit 1; }
else
    echo -e "${RED}No MPI launcher (mpiexec/srun) found${NC}"; exit 1
fi
cd ..

PRODUCED="output/Regression_${LANE}/${OUTPUT_FILENAME}"
[ ! -f "$PRODUCED" ] && { echo -e "${RED}Expected output not produced: ${PRODUCED}${NC}"; exit 1; }

# --- Bless or compare --------------------------------------------------------
if [ "$BLESS" = true ]; then
    mkdir -p "$(dirname "$BASELINE_FILE")"
    cp "$PRODUCED" "$BASELINE_FILE"
    echo -e "${GREEN}Blessed baseline: ${BASELINE_FILE}${NC}"
    echo -e "${YELLOW}Commit this to the Test-Data repo (git-LFS) to persist.${NC}"
    exit 0
fi

if [ ! -f "$BASELINE_FILE" ]; then
    echo -e "${RED}Baseline not found: ${BASELINE_FILE}${NC}"
    echo -e "${RED}Generate it once with --bless (and commit to Test-Data).${NC}"
    exit 1
fi

python_exe=$(command -v python3 || command -v python)
cmp_args=(--github-annotations --report-json "$hicar_repo/tests/regression_report_${LANE}.json" --figures-dir "${FIGURES_DIR}/${LANE}")
if [ "$MODE" = "tolerance" ]; then
    cmp_args+=(--tolerance-spec "$TOL_SPEC")
else
    cmp_args+=(--mode exact)
fi

echo "-------------------------------------------------------"
"$python_exe" "$COMPARE_SCRIPT" "$BASELINE_FILE" "$PRODUCED" "${cmp_args[@]}"
