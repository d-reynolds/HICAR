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
#   compare_cpu_gpu.sh <hicar_repo> --build <ref_build_dir> [--tolerance-spec PATH]
#
# --build: compile BOTH exes fresh from the CURRENT source (CPU = OPENACC=OFF,
#   GPU = OPENACC=ON, same compiler) into sub-build trees under <ref_build_dir>,
#   inheriting the compiler/library config from that build's CMakeCache.txt. This
#   is what `make test_gpu` uses, so the comparison is fully self-contained.
#   Needs an NVHPC-configured reference build (FC=nvfortran) for the GPU exe.
# =============================================================================
set -euo pipefail
BLUE='\033[0;36m'; GREEN='\033[0;32m'; RED='\033[0;31m'; NC='\033[0m'

# check if this system has a GPU (nvidia-smi) and an NVHPC compiler (nvfortran) in PATH.
if ! command -v nvidia-smi >/dev/null 2>&1; then
    echo -e "${RED}nvidia-smi not found; no GPU detected${NC}"; exit 1
fi
if ! command -v nvfortran >/dev/null 2>&1; then
    echo -e "${RED}nvfortran not found; NVHPC compiler required for GPU build${NC}"; exit 1
fi

if [ $# -lt 3 ]; then
    echo "Usage: $0 <hicar_repo> <cpu_exe> <gpu_exe> [--tolerance-spec PATH]"
    echo "       $0 <hicar_repo> --build <ref_build_dir> [--tolerance-spec PATH]"
    exit 1
fi
hicar_repo=$(cd "$1" && pwd); shift
TOL_SPEC="$hicar_repo/tests/tolerances/tolerances_gpu_cpu.yaml"

# Default build trees for --build mode. The GPU workflow stages its pre-built
# CPU-reference and GPU exes here so the comparison reuses them; a local run
# builds them. Removed on exit (see the trap below) so the next run builds fresh.
BUILD_REF=""
GPU_COMPARE_CPU_BUILD="${hicar_repo}/build_gpu_compare_cpu"
GPU_COMPARE_GPU_BUILD="${hicar_repo}/build_gpu_compare_gpu"

# find_hicar_exe <build_tree>: print the HICAR exe in <build_tree> (any MODE
# variant), or return 1 if none is present.
find_hicar_exe() {
    local tree="$1" n
    for n in HICAR_debug_gpu HICAR_debug HICAR_gpu HICAR_profile HICAR; do
        [ -x "$tree/$n" ] && { echo "$tree/$n"; return 0; }
    done
    return 1
}

# build_hicar_variant <out_tree> <label> [extra -D...]: print the HICAR exe in
# <out_tree>, building it first ONLY if it is not already there — so an exe the
# GPU workflow already staged is reused (skip compiling). Config is inherited
# from the reference build's cache ($BUILD_REF/CMakeCache.txt); cmake chatter
# goes to <out_tree>.{cfg,build}.log.
build_hicar_variant() {
    local out_tree="$1" label="$2"; shift 2
    local existing
    if existing=$(find_hicar_exe "$out_tree"); then
        echo -e "${GREEN}Reusing pre-built ${label} exe: ${existing}${NC}" >&2
        echo "$existing"; return 0
    fi
    local args=() v val
    if [ -f "$BUILD_REF/CMakeCache.txt" ]; then
        for v in MODE FC NETCDF_DIR FFTW_DIR MPI_DIR PETSC_DIR FSM FSM_DIR ASSERTIONS NCCL SNOWPACK_DIR; do
            val=$(grep "^${v}:" "$BUILD_REF/CMakeCache.txt" 2>/dev/null | head -n1 | sed 's/^[^=]*=//' || true)
            [ -n "$val" ] && args+=("-D${v}=${val}")
        done
    fi
    args+=("$@")
    echo -e "${BLUE}Building HICAR (${label}) in ${out_tree} ...${NC}" >&2
    cmake -S "$hicar_repo" -B "$out_tree" "${args[@]}" >"${out_tree}.cfg.log" 2>&1 \
        || { echo -e "${RED}configure failed (${label}); tail ${out_tree}.cfg.log:${NC}" >&2; tail -n 25 "${out_tree}.cfg.log" >&2; return 1; }
    cmake --build "$out_tree" --target HICAR --parallel >"${out_tree}.build.log" 2>&1 \
        || { echo -e "${RED}build failed (${label}); tail ${out_tree}.build.log:${NC}" >&2; tail -n 25 "${out_tree}.build.log" >&2; return 1; }
    find_hicar_exe "$out_tree" || { echo -e "${RED}no HICAR exe produced in ${out_tree}${NC}" >&2; return 1; }
}

if [ "${1:-}" = "--build" ]; then
    BUILD_REF="${2:?--build requires a reference build dir}"; shift 2
else
    CPU_EXE="${1:?need <cpu_exe> or --build}"; GPU_EXE="${2:?need <gpu_exe>}"; shift 2
fi
while [ $# -gt 0 ]; do
    case "$1" in
        --tolerance-spec) TOL_SPEC="$2"; shift 2;;
        *) echo -e "${RED}Unknown arg: $1${NC}"; exit 1;;
    esac
done

if [ -n "$BUILD_REF" ]; then
    # Clean the default build trees on exit (success OR failure) so the next
    # `make test_gpu` builds fresh exes.
    trap 'rm -rf "$GPU_COMPARE_CPU_BUILD" "$GPU_COMPARE_GPU_BUILD"' EXIT
    CPU_EXE=$(build_hicar_variant "$GPU_COMPARE_CPU_BUILD" "CPU OPENACC=OFF" -DOPENACC=OFF) || exit 1
    GPU_EXE=$(build_hicar_variant "$GPU_COMPARE_GPU_BUILD" "GPU OPENACC=ON"  -DOPENACC=ON)  || exit 1
fi

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
# number of total ranks for the run is num gpu * 2, to get one IO process per GPU
np=$((np * 2))
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
# pyyaml is needed by compare_outputs.py --tolerance-spec.
python_exe=$(command -v python3 || command -v python)
if [ -z "$python_exe" ] || ! "$python_exe" -c "import xarray,numpy,netCDF4,yaml" 2>/dev/null; then
    venv="${hicar_repo}/tests/Test_Cases/venv"
    if [ ! -x "${venv}/bin/python" ]; then
        echo -e "Creating venv at ${venv}"
        "${python_exe:-python3}" -m venv "$venv"
        "${venv}/bin/pip" -q install numpy netCDF4 xarray matplotlib pyyaml
    else
        "${venv}/bin/pip" -q install pyyaml >/dev/null 2>&1 || true
    fi
    python_exe="${venv}/bin/python"
fi
"$python_exe" "$COMPARE" \
    "output/cmp_cpu/${OUTPUT_FILENAME}" \
    "output/cmp_gpu/${OUTPUT_FILENAME}" \
    --tolerance-spec "$TOL_SPEC" \
    --github-annotations \
    --report-json "$hicar_repo/tests/gpu_cpu_report.json" \
    --figures-dir "$hicar_repo/tests/figures/gpu_cpu"
