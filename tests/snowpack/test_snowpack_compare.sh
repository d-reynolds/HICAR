#!/bin/bash
# ===========================================================================
# SNOWPACK implementation comparison: C++ bindings  vs  native Fortran
# ===========================================================================
# Runs the full HICAR model twice on the Standard Gaudergrat domain with the
# SNOWPACK snow model, from a prescribed multi-layer snowpack, using two
# executables that differ ONLY in the snow driver:
#
#   reference (A): built with -DSNOWPACK_CPP=ON  (C++ Alpine3D/SNOWPACK)
#   candidate (B): default build                 (native-Fortran port)
#
# Both are CPU/gfortran builds reading the SAME seeded domain through the SAME
# namelist, so any output difference is the snow implementation. Outputs are
# compared over ALL variables against tests/tolerances/tolerances_snowpack.yaml.
#
# Usage:
#   test_snowpack_compare.sh <hicar_repo> <cpp_exe> <fortran_exe> [np]
#   test_snowpack_compare.sh <hicar_repo> --build <ref_build_dir> [np]
#
#   <hicar_repo>  repo root
#   <cpp_exe>     HICAR executable built with the C++ SNOWPACK   (reference)
#   <fortran_exe> HICAR executable built with the Fortran SNOWPACK (candidate)
#   [np]          MPI ranks for each run (default 10; both runs use the same)
#
#   --build <ref_build_dir>  compile BOTH exes fresh from the CURRENT source
#                 (C++ = -DSNOWPACK_CPP=ON, Fortran = the default port; both CPU)
#                 into sub-build trees under <ref_build_dir>, inheriting the
#                 compiler/library config from that build's CMakeCache.txt. This
#                 is what `make test_snowpack` uses, so it is self-contained.
#
# Exit 0 = within tolerance; 1 = a run failed or outputs diverged.
# ===========================================================================
set -u

BLUE='\033[0;36m'; GREEN='\033[0;32m'; RED='\033[0;31m'; NC='\033[0m'

if [ $# -lt 2 ]; then
    echo "Usage: $0 <hicar_repo> <cpp_exe> <fortran_exe> [np]"
    echo "       $0 <hicar_repo> --build <ref_build_dir> [np]"
    exit 2
fi

hicar_repo="$(cd "$1" && pwd)"

# set_var <nml_file> <name> <value> <group>: set a namelist variable by name,
# inserting it into &<group> if absent (replaces the old `sed -i` edits). See
# helpers/example_namelists/set_nml_var.py.
set_var() { "${PYTHON:-python3}" "$hicar_repo/helpers/example_namelists/set_nml_var.py" "$1" "$2" "$3" --group "$4" --insert || exit 1; }

# Default build trees for --build mode. The snowpack-compare workflow stages its
# pre-built exes here so the comparison step reuses them; a local run builds them.
# Removed on exit (see the trap below) so the next run always builds fresh.
BUILD_REF=""
SNOWPACK_CPP_BUILD="${hicar_repo}/build_snowpack_cpp"
SNOWPACK_FORTRAN_BUILD="${hicar_repo}/build_snowpack_fortran"

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
# snowpack-compare workflow already staged is reused (skip compiling). Config is
# inherited from the reference build's cache ($BUILD_REF/CMakeCache.txt); cmake
# chatter goes to <out_tree>.{cfg,build}.log.
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
            val=$(grep "^${v}:" "$BUILD_REF/CMakeCache.txt" 2>/dev/null | head -n1 | sed 's/^[^=]*=//')
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

# --build: reuse/compile both exes at the default trees; otherwise use the two
# provided exe paths.
if [ "$2" = "--build" ]; then
    BUILD_REF="${3:?--build requires a reference build dir}"
    np_in="${4:-}"
    # Always clean the default build trees on exit (success OR failure) so the
    # next `make test_snowpack` builds fresh exes.
    trap 'rm -rf "$SNOWPACK_CPP_BUILD" "$SNOWPACK_FORTRAN_BUILD"' EXIT
    cpp_exe=$(build_hicar_variant "$SNOWPACK_CPP_BUILD"     "C++ SNOWPACK"           -DSNOWPACK_CPP=ON  -DOPENACC=OFF) || exit 1
    fortran_exe=$(build_hicar_variant "$SNOWPACK_FORTRAN_BUILD" "native-Fortran SNOWPACK" -DSNOWPACK_CPP=OFF -DOPENACC=OFF) || exit 1
else
    if [ $# -lt 3 ]; then
        echo "Usage: $0 <hicar_repo> <cpp_exe> <fortran_exe> [np]"
        echo "       $0 <hicar_repo> --build <ref_build_dir> [np]"
        exit 2
    fi
    cpp_exe="$2"
    fortran_exe="$3"
    np_in="${4:-}"
fi
# Ranks per run. Default to min(10, available cores) so the run is fast on a
# beefy machine but does NOT oversubscribe a small GitHub CI runner (2-4 cores).
if [ -n "${np_in:-}" ]; then
    np="$np_in"
else
    if [[ "$OSTYPE" == "darwin"* ]]; then
        ncores=$(sysctl -n hw.logicalcpu 2>/dev/null || echo 4)
    else
        ncores=$(nproc 2>/dev/null || echo 4)
    fi
    np=$(( ncores < 10 ? ncores : 10 ))
fi

for exe in "$cpp_exe" "$fortran_exe"; do
    if [ ! -x "$exe" ]; then
        echo -e "${RED}HICAR executable not found/executable: ${exe}${NC}"
        exit 1
    fi
done
# Resolve to absolute paths (we cd into the input dir to run).
cpp_exe="$(cd "$(dirname "$cpp_exe")" && pwd)/$(basename "$cpp_exe")"
fortran_exe="$(cd "$(dirname "$fortran_exe")" && pwd)/$(basename "$fortran_exe")"

input_dir="${hicar_repo}/tests/Test_Cases/input"
domains_dir="${hicar_repo}/tests/Test_Cases/domains"
nmlgen_dir="${input_dir}/nml_gen_scripts"
compare_script="${hicar_repo}/tests/compare_outputs.py"
tol_spec="${hicar_repo}/tests/tolerances/tolerances_snowpack.yaml"
seeder="${hicar_repo}/tests/snowpack/make_snowpack_init.py"
seeded_domain="${domains_dir}/Gaudergrat_250m_snowpack.nc"

# --- MPI launcher (same heuristic as test_reproducibility.sh) ---------------
mpiexec_path=""
IFS=':' read -ra PATH_DIRS <<< "$PATH"
for dir in "${PATH_DIRS[@]}"; do
    if [ -x "$dir/mpiexec" ] && ! [[ "$dir" =~ python|conda ]]; then
        mpiexec_path="$dir/mpiexec"; break
    fi
done
[ -z "$mpiexec_path" ] && mpiexec_path="$(command -v mpiexec || command -v mpirun || true)"
if [ -z "$mpiexec_path" ]; then
    echo -e "${RED}No mpiexec/mpirun found in PATH${NC}"; exit 1
fi

# --- Python with xarray/numpy/netCDF4 (+pyyaml for the tolerance spec) -------
python_exe="$(command -v python3 || command -v python || true)"
if [ -z "$python_exe" ] || ! "$python_exe" -c "import xarray,numpy,netCDF4,yaml" 2>/dev/null; then
    venv="${hicar_repo}/tests/Test_Cases/venv"
    if [ ! -x "${venv}/bin/python" ]; then
        echo -e "Creating venv at ${BLUE}${venv}${NC}"
        "${python_exe:-python3}" -m venv "$venv"
        "${venv}/bin/pip" -q install numpy netCDF4 xarray matplotlib pyyaml
    else
        "${venv}/bin/pip" -q install pyyaml >/dev/null 2>&1 || true
    fi
    python_exe="${venv}/bin/python"
fi

export OMP_NUM_THREADS=1

echo
echo "======================================================="
echo -e "  ${BLUE}SNOWPACK C++ vs native-Fortran comparison${NC}"
echo "======================================================="
echo "  reference (C++):     $cpp_exe"
echo "  candidate (Fortran): $fortran_exe"
echo "  ranks each run:      $np"

# --- forcing list + support files (mirror test_reproducibility.sh) ----------
# Shift the shipped night-time forcing (00-03 UTC) to midday (10-13 UTC) so the
# run exercises the shortwave path, then build the file list from the shifted
# copies. Regenerated every run so it tracks the namelist's start/end window.
midday_forcing="${hicar_repo}/tests/Test_Cases/forcing_midday"
echo
echo -e "Shifting forcing +10 h (night -> midday) for the shortwave-on comparison..."
"$python_exe" "${hicar_repo}/tests/shift_forcing_timestamps.py" \
    "${hicar_repo}/tests/Test_Cases/forcing" "$midday_forcing" 10 \
    || { echo -e "${RED}forcing time-shift failed${NC}"; exit 1; }
"${hicar_repo}/helpers/filelist_script.sh" \
    "${midday_forcing}/*" "${input_dir}/file_list_TestCase.txt"
[ -f "${input_dir}/VEGPARM.TBL" ]        || cp "${hicar_repo}/run/"*.TBL "${input_dir}/"
[ -d "${input_dir}/rrtmg_support" ]      || cp -r "${hicar_repo}/run/rrtmg_support" "${input_dir}/"
[ -d "${input_dir}/mp_support" ]         || cp -r "${hicar_repo}/run/mp_support"   "${input_dir}/"

# --- seed the domain with the prescribed snowpack ---------------------------
echo
echo -e "Seeding domain with a 14-layer / 1 m snowpack on all land cells..."
"$python_exe" "$seeder" "${domains_dir}/Gaudergrat_250m.nc" "$seeded_domain" \
    --n-layers 14 --depth 1.0 || { echo -e "${RED}domain seeding failed${NC}"; exit 1; }

# The comparison namelist is generated from the committed alpine_realdata
# example by nml_gen_scripts/SNOWPACK_Compare.sh (same pattern as the other test
# cases), so there is no separate base-namelist step here.

# run_one <exe> <label>  -> sets RUN_OUT to the produced output .nc, returns status
run_one() {
    local exe="$1" label="$2"
    local out_dir="${hicar_repo}/tests/Test_Cases/output/SNOWPACK_${label}"
    local rst_dir="${hicar_repo}/tests/Test_Cases/restart/SNOWPACK_${label}"
    rm -rf "$out_dir" "$rst_dir"; mkdir -p "$out_dir" "$rst_dir"

    local nml="${input_dir}/SNOWPACK_${label}.nml"
    ( cd "$input_dir" && bash "${nmlgen_dir}/SNOWPACK_Compare.sh" "$(basename "$nml")" )
    # per-build output/restart folders
    set_var "$nml" output_folder  "'../output/SNOWPACK_${label}/'"  output
    set_var "$nml" restart_folder "'../restart/SNOWPACK_${label}/'" restart

    echo
    echo -e "Running ${BLUE}${label}${NC} ($exe) with ${np} ranks ..."
    ( cd "$input_dir" && "$mpiexec_path" -np "$np" "$exe" "$(basename "$nml")" ) \
        > "${input_dir}/SNOWPACK_${label}.out" 2> "${input_dir}/SNOWPACK_${label}.err"
    local status=$?
    # HICAR writes its fatal messages (incl. a bare `stop`, which exits 0) to
    # STDOUT, not stderr -- so always surface the tail of the .out on failure.
    if [ $status -ne 0 ]; then
        echo -e "${RED}${label} run failed (exit ${status}). Last stdout:${NC}"
        tail -n 25 "${input_dir}/SNOWPACK_${label}.out"
        echo -e "${RED}Last stderr:${NC}"
        tail -n 15 "${input_dir}/SNOWPACK_${label}.err"
        RUN_OUT=""
        return 1
    fi
    RUN_OUT="$(ls -1 "${out_dir}"/*.nc 2>/dev/null | head -1)"
    if [ -z "$RUN_OUT" ]; then
        # Exit 0 but no output -> typically a bare `stop` (exit 0) during init.
        # The reason is in stdout; show it so CI logs are self-diagnosing.
        echo -e "${RED}${label}: no output .nc produced in ${out_dir}. Last stdout:${NC}"
        tail -n 25 "${input_dir}/SNOWPACK_${label}.out"
        return 1
    fi
    echo -e "${GREEN}${label} completed${NC} -> $(basename "$RUN_OUT")"
    return 0
}

run_one "$cpp_exe"     "cpp"     || exit 1
ref_out="$RUN_OUT"
run_one "$fortran_exe" "fortran" || exit 1
cand_out="$RUN_OUT"

# --- compare ----------------------------------------------------------------
echo
echo "Comparing all output variables (tol spec: $(basename "$tol_spec")) ..."
# The JSON report records max|abs|/max|rel|/violation counts for EVERY variable
# (passing ones included) — on a main-line pass the snowpack-compare workflow's
# archive job publishes it as a GitHub release so the parity level can be tracked.
report_dir="${hicar_repo}/tests/figures/snowpack_compare"
mkdir -p "$report_dir"
"$python_exe" "$compare_script" "$ref_out" "$cand_out" \
    --mode tolerance --tolerance-spec "$tol_spec" \
    --figures-dir "$report_dir" \
    --report-json "${report_dir}/parity_report.json"
cmp_status=$?

# Spatial difference maps for the key snow/surface fields, generated on PASS as
# well as FAIL (compare_outputs only plots failures): these are the figures the
# archive job snapshots, and the visual baseline for "how did the residual change".
"$python_exe" "${hicar_repo}/tests/snowpack/plot_snowpack_diff.py" \
    "$ref_out" "$cand_out" "${report_dir}/diffmaps" >/dev/null 2>&1 \
    && echo -e "Difference maps: ${BLUE}${report_dir}/diffmaps${NC}" \
    || echo -e "${RED}(diff-map generation failed — non-fatal)${NC}"


echo
echo "======================================================="
if [ $cmp_status -eq 0 ]; then
    echo -e "  SNOWPACK C++ vs Fortran: ${GREEN}PASS${NC}"
    echo "======================================================="
    exit 0
else
    echo -e "  SNOWPACK C++ vs Fortran: ${RED}FAIL${NC}"
    echo "======================================================="
    exit 1
fi
