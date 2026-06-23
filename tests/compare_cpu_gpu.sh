#!/bin/bash
# =============================================================================
# GPU-vs-CPU equivalence / regression for the self-hosted GPU lane.
#
# Runs the Standard test case with TWO HICAR executables built from the SAME
# commit by the SAME compiler (NVHPC) — one CPU/host (OPENACC=OFF) and one
# GPU (OPENACC=ON) — using an IDENTICAL namelist, then compares the outputs by
# tolerance. Building both with NVHPC isolates the GPU port: the only difference
# is CPU vs GPU code generation/execution, not the compiler.
#
# Usage:
#   compare_cpu_gpu.sh <hicar_repo> <cpu_exe> <gpu_exe> [--tolerance-spec PATH]
# =============================================================================
set -euo pipefail
BLUE='\033[0;36m'; GREEN='\033[0;32m'; RED='\033[0;31m'; NC='\033[0m'

if [ $# -lt 3 ]; then
    echo "Usage: $0 <hicar_repo> <cpu_exe> <gpu_exe> [--tolerance-spec PATH]"; exit 1
fi
hicar_repo=$(cd "$1" && pwd); CPU_EXE="$2"; GPU_EXE="$3"; shift 3
TOL_SPEC="$hicar_repo/tests/tolerances.yaml"
while [ $# -gt 0 ]; do
    case "$1" in
        --tolerance-spec) TOL_SPEC="$2"; shift 2;;
        *) echo -e "${RED}Unknown arg: $1${NC}"; exit 1;;
    esac
done

[ -x "$CPU_EXE" ] || { echo -e "${RED}CPU exe not found/executable: $CPU_EXE${NC}"; exit 1; }
[ -x "$GPU_EXE" ] || { echo -e "${RED}GPU exe not found/executable: $GPU_EXE${NC}"; exit 1; }

COMPARE="$hicar_repo/tests/compare_outputs.py"
STANDARD_GEN="$hicar_repo/tests/Test_Cases/input/nml_gen_scripts/Standard.sh"
OUTPUT_FILENAME="Gaudergrat_250m_2017-02-14_00-00-00.nc"

# set_var <nml_file> <name> <value> <group>: set a namelist variable by name,
# inserting it into &<group> if absent (replaces the old `sed -i` edits). See
# helpers/example_namelists/set_nml_var.py.
set_var() { "${PYTHON:-python3}" "$hicar_repo/helpers/example_namelists/set_nml_var.py" "$1" "$2" "$3" --group "$4" --insert || exit 1; }
# Radiation scheme used for BOTH runs so physics is identical. rrtmgp is the
# portable RTE+RRTMGP scheme (runs on host with OPENACC=OFF and on the GPU). If a
# clean run needs the legacy rrtmg instead, set RAD=rrtmg — but keep it the SAME
# for both exes, or the comparison is meaningless.
RAD="${RAD:-rrtmgp}"

cd "$hicar_repo/tests/Test_Cases"
# Support files + forcing list.
[ ! -f ./input/VEGPARM.TBL ]   && cp "$hicar_repo/run/"*.TBL ./input/
[ ! -d ./input/rrtmg_support ] && cp -r "$hicar_repo/run/rrtmg_support" ./input
[ ! -d ./input/mp_support ]    && cp -r "$hicar_repo/run/mp_support" ./input
[ ! -d ./input/rrtmgp_support ] && cp -r "$hicar_repo/run/rrtmgp_support" ./input
[ ! -f ./input/file_list_TestCase.txt ] && \
    "$hicar_repo/helpers/filelist_script.sh" "forcing/*" input/file_list_TestCase.txt

# One rank per GPU; use the SAME rank count for the CPU run so the decomposition
# matches and the comparison is apples-to-apples (reproducibility is verified
# separately in the CPU full-test).
np=$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | wc -l | tr -d ' ')
[ "${np:-0}" -lt 1 ] && np=1
export OMP_NUM_THREADS=1

mpiexec_path=""
IFS=':' read -ra PD <<< "$PATH"
for d in "${PD[@]}"; do
    if [ -x "$d/mpiexec" ] && ! [[ "$d" =~ python|conda ]]; then mpiexec_path="$d/mpiexec"; break; fi
done
[ -z "$mpiexec_path" ] && { echo -e "${RED}mpiexec not found${NC}"; exit 1; }

# ---- run one exe on the Standard case, output to a named folder -------------
run_one() {  # $1=exe  $2=tag
    local exe="$1" tag="$2"
    local out_dir="../output/cmp_${tag}" rst_dir="../restart/cmp_${tag}"
    rm -rf "output/cmp_${tag}" "restart/cmp_${tag}"
    mkdir -p "output/cmp_${tag}" "restart/cmp_${tag}"
    cd input
    local nml="cmp_${tag}.nml"
    "$exe" --gen-nml "$nml"
    bash "$STANDARD_GEN" "$nml"
    set_var "$nml" output_vars    "'all'"          output
    set_var "$nml" output_folder  "'${out_dir}/'"  output
    set_var "$nml" restart_folder "'${rst_dir}/'"  restart
    set_var "$nml" rad            "'${RAD}'"       physics
    echo -e "Running ${BLUE}${tag}${NC} (${exe##*/}) with ${np} ranks, rad=${RAD}..."
    "$mpiexec_path" -np "$np" "$exe" "$nml" 1>"cmp_${tag}.out" 2>"cmp_${tag}.err" || {
        echo -e "${RED}${tag} run failed; last stderr:${NC}"; tail -n 20 "cmp_${tag}.err"; exit 1; }
    cd ..
    [ -f "output/cmp_${tag}/${OUTPUT_FILENAME}" ] || \
        { echo -e "${RED}${tag}: expected output not produced${NC}"; exit 1; }
}

run_one "$CPU_EXE" cpu
run_one "$GPU_EXE" gpu

echo "-------------------------------------------------------"
echo -e "Comparing GPU vs CPU (reference) within tolerance..."
python_exe=$(command -v python3 || command -v python)
"$python_exe" "$COMPARE" \
    "output/cmp_cpu/${OUTPUT_FILENAME}" \
    "output/cmp_gpu/${OUTPUT_FILENAME}" \
    --tolerance-spec "$TOL_SPEC" \
    --github-annotations \
    --report-json "$hicar_repo/tests/gpu_cpu_report.json" \
    --figures-dir "$hicar_repo/tests/figures/gpu_cpu"
