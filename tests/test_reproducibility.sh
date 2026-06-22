#!/bin/bash

# Reproducibility tests for HICAR:
#   decomposition: MPI Reproducibility — run with different rank counts, compare output (must be identical)
#   restart:       Restart Reproducibility — continuous run vs restart run, compare final timestep
#   all:           Run both tests
#
# GPU mode: If HICAR_debug_gpu is found, uses all available GPUs (N) and half (N/2) for decomposition.
#           Requires at least 2 GPUs for the decomposition test.
# CPU mode: Falls back to HICAR_debug with 5 vs 10 ranks.
#
# Usage: test_reproducibility.sh <hicar_repo> [decomposition|restart|all]

# Define colors for output styling
BLUE='\033[0;36m'
GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Get location of HICAR repo (macOS-compatible, no readlink -f)
hicar_repo=$(cd "$1" && pwd)
test_mode="${2:-all}"

# Validate test mode
case "$test_mode" in
    decomposition|restart|all) ;;
    *)
        echo -e "${RED}Unknown test mode: ${test_mode}${NC}"
        echo "Usage: $0 <hicar_repo> [decomposition|restart|all]"
        exit 1
        ;;
esac

# Get the location of a HICAR executable
# Prefer HICAR_debug_gpu (OpenACC release build) if available, else fall back to HICAR_debug
use_gpu=false
if [ -f "$hicar_repo/bin/HICAR_debug_gpu" ]; then
    hicar_exe="$hicar_repo/bin/HICAR_debug_gpu"
    use_gpu=true
    echo -e "${GREEN}Found HICAR_debug_gpu — running tests in GPU mode${NC}"
elif [ -f "$hicar_repo/bin/HICAR_debug" ]; then
    hicar_exe="$hicar_repo/bin/HICAR_debug"
    echo -e "${GREEN}Found HICAR_debug — running tests in CPU mode${NC}"
else
    echo -e "${RED}No HICAR executable found in bin directory.${NC}"
    echo -e "${RED}Build with MODE=debug (CPU) or MODE=release with OpenACC (GPU) and run 'make install'.${NC}"
    exit 1
fi

# Detect GPU count if running in GPU mode
num_gpus=0
if [ "$use_gpu" = true ]; then
    if command -v nvidia-smi &> /dev/null; then
        num_gpus=$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | wc -l | tr -d ' ')
    fi
    if [ "$num_gpus" -lt 1 ]; then
        echo -e "${RED}GPU mode selected but no GPUs detected via nvidia-smi.${NC}"
        exit 1
    fi
    echo -e "${GREEN}Detected ${num_gpus} GPU(s)${NC}"
fi

# Path to the compare script and figures output directory
compare_script="$hicar_repo/tests/compare_outputs.py"
figures_dir="$hicar_repo/tests/figures"

# Output file name pattern (from Standard.sh: domain=Gaudergrat_250m, start=2017-02-14)
OUTPUT_FILENAME="Gaudergrat_250m_2017-02-14_00-00-00.nc"

# Make forcing file list
if [ ! -f "${hicar_repo}/tests/Test_Cases/input/file_list_TestCase.txt" ]; then
    "$hicar_repo/helpers/filelist_script.sh" "${hicar_repo}/tests/Test_Cases/forcing/*" "${hicar_repo}/tests/Test_Cases/input/file_list_TestCase.txt"
fi

# Copy supporting files if missing
if [ ! -f "${hicar_repo}/tests/Test_Cases/input/VEGPARM.TBL" ]; then
    echo "Copying .TBL files to ${hicar_repo}/tests/Test_Cases/input"
    cp "${hicar_repo}/run/"*.TBL "${hicar_repo}/tests/Test_Cases/input/"
fi
if [ ! -d "${hicar_repo}/tests/Test_Cases/input/rrtmg_support" ]; then
    cp -r "${hicar_repo}/run/rrtmg_support" "${hicar_repo}/tests/Test_Cases/input"
fi
if [ ! -d "${hicar_repo}/tests/Test_Cases/input/mp_support" ]; then
    cp -r "${hicar_repo}/run/mp_support" "${hicar_repo}/tests/Test_Cases/input"
fi

if [ "$use_gpu" = true ]; then
    if [ ! -d "${hicar_repo}/tests/Test_Cases/input/rrtmgp_support" ]; then
        cp -r "${hicar_repo}/run/rrtmgp_support" "${hicar_repo}/tests/Test_Cases/input"
    fi
fi
# Generate default namelist
default_file="${hicar_repo}/tests/Test_Cases/input/default_hicar_options.nml"
if [ -f "$default_file" ]; then
    rm "$default_file"
fi
echo "Generating default namelist to ${hicar_repo}/tests/Test_Cases/input"
$hicar_exe --gen-nml "$default_file"

# Find mpiexec (same logic as test_case_runner.sh)
mpiexec_path=""
IFS=':' read -ra PATH_DIRS <<< "$PATH"
for dir in "${PATH_DIRS[@]}"; do
    if [ -x "$dir/mpiexec" ] && ! [[ "$dir" =~ python|conda ]]; then
        mpiexec_path="$dir/mpiexec"
        break
    fi
done

if [ -z "$mpiexec_path" ]; then
    if command -v srun &> /dev/null; then
        echo -e "${GREEN}Using srun to run HICAR${NC}"
        # Get any srun flags
        SRUN_FLAGS=""
        for arg in "$@"; do
            if [[ "$arg" == -SRUN_FLAGS=* ]]; then
                SRUN_FLAGS="${arg#-SRUN_FLAGS=}"
                SRUN_FLAGS="${SRUN_FLAGS#\'}"
                SRUN_FLAGS="${SRUN_FLAGS%\'}"
                break
            fi
        done
    else
        echo -e "${RED}mpiexec not found. Please install mpiexec.${NC}"
        exit 1
    fi
else
    echo -e "${GREEN}Using mpiexec from: $mpiexec_path${NC}"
fi

# Setup Python environment
python_exe=$(which python 2>/dev/null || which python3 2>/dev/null)
if [ -z "$python_exe" ]; then
    echo -e "${RED}Python is not installed. Please install Python.${NC}"
    exit 1
fi

if ! $python_exe -c "import xarray, numpy, netCDF4" &> /dev/null; then
    PY_ENV_PATH=${hicar_repo}/tests/Test_Cases/venv
    echo -e "Creating virtual environment at ${BLUE}${PY_ENV_PATH}${NC}"
    mkdir -p "$PY_ENV_PATH"
    $python_exe -m venv "${PY_ENV_PATH}"
    "${PY_ENV_PATH}/bin/pip" install numpy netCDF4 xarray matplotlib
    python_exe="${PY_ENV_PATH}/bin/python"
fi

export OMP_NUM_THREADS=1

# Cleanup trap for temp files
cleanup() {
    rm -f "${hicar_repo}/tests/Test_Cases/input"/*.bak
}
trap cleanup EXIT

# -------------------------------------------------------
# Helper: generate_standard_nml <output_nml> <output_dir> <restart_dir>
# -------------------------------------------------------
generate_standard_nml() {
    local out_nml="$1"
    local output_dir="$2"
    local restart_dir="$3"
    # Optional 4th arg: end_date. Defaults to the 10-min window the restart test
    # needs (its checkpoints land at 00:05 and 00:10). The decomposition test only
    # runs-and-compares (no restart), so it passes a shorter window to keep the
    # slow debug build fast.
    local end_date="${4:-2017-02-14 00:10:00}"

    cd ${hicar_repo}/tests/Test_Cases/input
    cp "$default_file" "$out_nml"

    # Apply Standard.sh settings
    ${hicar_repo}/tests/Test_Cases/input/nml_gen_scripts/Standard.sh "$(basename "$out_nml")"
    cd ..

    # Override wind, Sx, and output_vars for reproducibility testing. wind is at
    # its default in the example (so absent from the namelist) — insert it.
    "${PYTHON:-python3}" "${hicar_repo}/helpers/example_namelists/set_nml_var.py" "$out_nml" wind "'none'" --group physics --insert
    sed -i'.bak' 's/Sx = .True./Sx = .False./g' "$out_nml"
    sed -i'.bak' "s/output_vars = .*$/output_vars = 'all'/g" "$out_nml"

    if [ $use_gpu = true ]; then
        # GPU mode: use all available GPUs and half for decomposition test, so set end time to 10 min for both runs
        sed -i'.bak' "s/rad = 'rrtmg'/rad = 'rrtmgp'/g" "$out_nml"
    fi
    # Shorten run (default 10 min, or the override), output every 5 min for restart checkpoints
    sed -i'.bak' "s/end_date = '2017-02-14 00:20:00'/end_date = '${end_date}'/g" "$out_nml"
    sed -i'.bak' 's/outputinterval = 600/outputinterval = 300/g' "$out_nml"

    # Override output and restart folders
    sed -i'.bak' "s|output_folder = '../output/Standard/'|output_folder = '${output_dir}'|g" "$out_nml"
    sed -i'.bak' "s|restart_folder = '../restart/Standard/'|restart_folder = '${restart_dir}'|g" "$out_nml"
    rm -f "${out_nml}.bak"

    cd $hicar_repo

}


# -------------------------------------------------------
# Helper: generate_restart_nml <output_nml> <output_dir> <restart_dir>
# -------------------------------------------------------
generate_restart_nml() {
    local out_nml="$1"
    local output_dir="$2"
    local restart_dir="$3"

    cd ${hicar_repo}/tests/Test_Cases/input
    cp "$default_file" "$out_nml"

    # Apply Standard_restart.sh settings
    ${hicar_repo}/tests/Test_Cases/input/nml_gen_scripts/Standard_restart.sh "$(basename "$out_nml")"
    cd ..

    # Override wind, Sx, and output_vars for reproducibility testing. wind is at
    # its default in the example (so absent from the namelist) — insert it.
    "${PYTHON:-python3}" "${hicar_repo}/helpers/example_namelists/set_nml_var.py" "$out_nml" wind "'none'" --group physics --insert
    sed -i'.bak' 's/Sx = .True./Sx = .False./g' "$out_nml"
    sed -i'.bak' "s/output_vars = .*$/output_vars = 'all'/g" "$out_nml"

    if [ $use_gpu = true ]; then
        # GPU mode: use all available GPUs and half for decomposition test, so set end time to 10 min for both runs
        sed -i'.bak' "s/rad = 'rrtmg'/rad = 'rrtmgp'/g" "$out_nml"
    fi

    # Shorten run: end at 10 min, restart from 5 min checkpoint, output every 5 min
    sed -i'.bak' "s/end_date = '2017-02-14 00:20:00'/end_date = '2017-02-14 00:10:00'/g" "$out_nml"
    sed -i'.bak' "s/restart_date = '2017-02-14 00:10:00'/restart_date = '2017-02-14 00:05:00'/g" "$out_nml"
    #sed -i'.bak' 's/outputinterval = 600/outputinterval = 300/g' "$out_nml"

    # Override output and restart folders
    sed -i'.bak' "s|output_folder = '../output/Standard/'|output_folder = '${output_dir}'|g" "$out_nml"
    sed -i'.bak' "s|restart_folder = '../restart/Standard/'|restart_folder = '${restart_dir}'|g" "$out_nml"
    rm -f "${out_nml}.bak"

    cd $hicar_repo
}

# -------------------------------------------------------
# Helper: run_hicar <nml_file> <np> <label>
# -------------------------------------------------------
run_hicar() {
    local nml_file="$1"
    local run_np="$2"
    local label="$3"

    echo
    echo -e "Running: ${BLUE}${label}${NC} with ${run_np} MPI ranks"

    local base_label
    base_label=$(echo "$label" | tr ' ' '_')

    cd $hicar_repo/tests/Test_Cases/input

    if [ ! -z "$mpiexec_path" ]; then
        $mpiexec_path -np "$run_np" $hicar_exe "$nml_file" 1>"${base_label}.out" 2>"${base_label}.err" &
        hicar_pid=$!
    elif command -v srun &> /dev/null; then
        srun -N 1 -n "$run_np" $SRUN_FLAGS $hicar_exe "$nml_file" 1>"${base_label}.out" 2>"${base_label}.err" &
        hicar_pid=$!
    else
        echo -e "${RED}No MPI launcher available${NC}"
        return 1
    fi

    echo -n "Initializing..."

    # Monitor progress
    local last_line_count=0
    local wait_counter=0
    local total_last_lines=0
    while kill -0 $hicar_pid 2>/dev/null; do
        sleep 1
        if [ $wait_counter -gt 300 ]; then
            echo
            echo -e "${RED}HICAR has not produced output in 5 minutes. Possible hang.${NC}"
            kill -9 $hicar_pid 2>/dev/null
            return 1
        fi
        if [ -f "${base_label}.out" ]; then
            total_current_lines=$(wc -l < "${base_label}.out" | tr -d ' ')
            if [ "$total_current_lines" -eq "$total_last_lines" ]; then
                wait_counter=$((wait_counter + 1))
            else
                wait_counter=0
            fi
            total_last_lines=$total_current_lines
            current_lines=$(grep -c "^ *Model time" "${base_label}.out" 2>/dev/null) || true
            if [ "$current_lines" -gt "$last_line_count" ]; then
                if [ $last_line_count -eq 0 ]; then
                    echo -e "\r\033[KRunning..."
                fi
                latest_line=$(grep -a "^ *Model time" "${base_label}.out" | tail -n 1)
                end_time=$(grep -a "^ *End  time" "${base_label}.out" | tail -n 1)
                echo -e "\r\033[K$latest_line"
                echo -n "$end_time"
                last_line_count=$current_lines
            fi
        fi
    done
    echo

    wait $hicar_pid
    local status=$?

    if [ $status -ne 0 ]; then
        echo -e "${RED}HICAR exited with error code $status for ${label}${NC}"
        if [ -f "${base_label}.err" ]; then
            echo -e "${RED}Last 10 lines of stderr:${NC}"
            tail -n 10 "${base_label}.err"
        fi
    else
        echo -e "${GREEN}${label} completed successfully${NC}"
    fi
    cd $hicar_repo

    return $status
}

# -------------------------------------------------------
# Helper: run_fm_decomposition  (sets test_fm_decomp_result)
# -------------------------------------------------------
# Fine-mesh (snow_drift) advection decomposition check. Unlike the full-model
# decomposition test above, this exercises the fine-mesh advection + flux limiter
# in isolation via the HICAR-tester unit executable. Run verbose ("-v"), the
# snow_drift suite prints qs_fm at fixed GLOBAL indices as "FMPROBE i j k val"
# (a normal, non-verbose unit run stays quiet); we run it at two rank counts and
# diff those lines. In verbose mode the driver does NOT redirect non-root stdout,
# so every rank prints its own interior cells — together the full domain — and the
# join below compares every cell both decompositions own. Independent of the
# full-model run and runs even when the full-model decomposition test is skipped
# (e.g. <2 GPUs), since the tester is a CPU executable. SKIPs if HICAR-tester was
# not built (the GPU lane skips its device link).
run_fm_decomposition() {
    echo
    echo "-------------------------------------------------------"
    echo -e "  ${BLUE}Fine-Mesh Advection Decomposition Test${NC}"
    echo "-------------------------------------------------------"

    # Locate the build dir whose tests/HICAR-tester exists.
    local build_dir="${HICAR_BUILD_DIR:-}"
    if [ -z "$build_dir" ] || [ ! -x "$build_dir/tests/HICAR-tester" ]; then
        build_dir=""
        for cand in "$hicar_repo/build" "$hicar_repo/Build" "$hicar_repo"; do
            if [ -x "$cand/tests/HICAR-tester" ]; then
                build_dir="$cand"
                break
            fi
        done
    fi

    if [ -z "$build_dir" ]; then
        echo -e "${BLUE}HICAR-tester not found (build it with 'make HICAR-tester'); skipping fine-mesh decomposition test.${NC}"
        test_fm_decomp_result="SKIP"
        return 0
    fi

    # MPI launcher prefix (everything up to the rank-count flag). Reuse the
    # launcher already discovered above for the full-model runs.
    local launcher=""
    if [ -n "$mpiexec_path" ]; then
        launcher="$mpiexec_path -np"
    elif command -v srun &> /dev/null; then
        launcher="srun -N 1 $SRUN_FLAGS -n"
    else
        echo -e "${BLUE}No MPI launcher available; skipping fine-mesh decomposition test.${NC}"
        test_fm_decomp_result="SKIP"
        return 0
    fi

    # 1 vs 4 ranks. With "-v" every rank prints its own interior, so the union is
    # the full domain regardless of rank count and the comparison covers every cell.
    local np_low=2 np_high=4
    local tmp
    tmp="$(mktemp -d)"

    # Run the snow_drift suite verbose at <np> ranks; emit the FMPROBE lines as a
    # sorted "ZZZ_ZZZ_ZZ value" stream keyed by global (i,j,k) so they can be
    # joined across decompositions. "-v" enables the dump (off by default).
    run_fm_one() {
        local np="$1" out="$2"
        ( cd "$build_dir" && OMP_NUM_THREADS=1 \
            $launcher "$np" ./tests/HICAR-tester -v snow_drift ) 2>&1 \
            | grep '^FMPROBE' | awk '{printf "%03d_%03d_%02d %s\n",$2,$3,$4,$5}' | sort > "$out"
    }

    echo -e "Running HICAR-tester snow_drift at ${BLUE}${np_low}${NC} and ${BLUE}${np_high}${NC} ranks"
    run_fm_one "$np_low"  "$tmp/low.txt"
    run_fm_one "$np_high" "$tmp/high.txt"

    local n_low n_high
    n_low=$(wc -l < "$tmp/low.txt"); n_high=$(wc -l < "$tmp/high.txt")
    if [ "$n_low" -eq 0 ] || [ "$n_high" -eq 0 ]; then
        echo -e "  ${RED}FAIL: a run produced no FMPROBE output (crash?). low=${n_low} high=${n_high}${NC}"
        test_fm_decomp_result="FAIL"
        rm -rf "$tmp"
        return 0
    fi

    # Compare the cells owned by BOTH decompositions (join on the global i_j_k key).
    join "$tmp/low.txt" "$tmp/high.txt" > "$tmp/joined.txt"
    awk '$2 != $3 {print}' "$tmp/joined.txt" > "$tmp/mismatch.txt"
    local n_common n_bad
    n_common=$(wc -l < "$tmp/joined.txt"); n_bad=$(wc -l < "$tmp/mismatch.txt")

    echo "  compared ${n_common} common interior cells"
    if [ "$n_bad" -eq 0 ]; then
        echo -e "  ${GREEN}PASS: fine-mesh advection is bit-for-bit decomposition-reproducible${NC}"
        test_fm_decomp_result="PASS"
    else
        echo -e "  ${RED}FAIL: ${n_bad} cell(s) differ between ${np_low} and ${np_high} ranks${NC}"
        echo "  first mismatches (cell  ${np_low}-rank  ${np_high}-rank):"
        head -8 "$tmp/mismatch.txt" | sed 's/^/    /'
        test_fm_decomp_result="FAIL"
    fi

    rm -rf "$tmp"
}

# ===============================================================
# Test Results Tracking
# ===============================================================
test_decomp_result=""
test_fm_decomp_result=""
test_restart_result=""

echo
echo "======================================================="
echo "  HICAR Reproducibility Tests (mode: ${test_mode})"
echo "======================================================="

# ===============================================================
# Test: MPI Domain Decomposition (5 ranks vs 10 ranks)
# ===============================================================
if [ "$test_mode" == "decomposition" ] || [ "$test_mode" == "all" ]; then

    # Determine rank counts for decomposition test
    if [ "$use_gpu" = true ]; then
        if [ "$num_gpus" -lt 2 ]; then
            echo
            echo "-------------------------------------------------------"
            echo -e "  ${BLUE}Domain Decomposition Test${NC}"
            echo "-------------------------------------------------------"
            echo -e "${RED}Only ${num_gpus} GPU detected. The decomposition test requires at least 2 GPUs.${NC}"
            echo -e "${RED}Skipping decomposition test.${NC}"
            test_decomp_result="SKIP"
        else
            decomp_np_full=$((num_gpus + 1))
            decomp_np_half=$(( (num_gpus / 2) + 1 ))
        fi
    else
        # Once all of the .nml files have been created, run the HICAR executable
        # with each of them
        # detect if this is running on mac OS
        if [[ "$OSTYPE" == "darwin"* ]]; then
            decomp_np_full=$(sysctl -n hw.logicalcpu)
        else
            decomp_np_full=$(nproc --all)
        fi

        decomp_np_half=$(( decomp_np_full/2 ))
    fi

    # Fine-mesh advection decomposition check (HICAR-tester; CPU, GPU-independent).
    run_fm_decomposition

    if [ "$test_decomp_result" != "SKIP" ]; then

    echo
    echo "-------------------------------------------------------"
    echo -e "  ${BLUE}Domain Decomposition Test (${decomp_np_half} vs ${decomp_np_full} ranks)${NC}"
    echo "-------------------------------------------------------"

    # Clear previous figures for this test
    rm -rf "${hicar_repo}/tests/Test_Cases/figures/decomposition"

    # Create output/restart directories
    for suffix in Repro_MPI_${decomp_np_half} Repro_MPI_${decomp_np_full}; do
        rm -rf "${hicar_repo}/tests/Test_Cases/output/${suffix}" "${hicar_repo}/tests/Test_Cases/restart/${suffix}"
        mkdir -p "${hicar_repo}/tests/Test_Cases/output/${suffix}" "${hicar_repo}/tests/Test_Cases/restart/${suffix}"
    done

    # Generate and run with half ranks (5-min window — decomposition test does not restart)
    generate_standard_nml "${hicar_repo}/tests/Test_Cases/input/Repro_MPI_${decomp_np_half}.nml" "${hicar_repo}/tests/Test_Cases/output/Repro_MPI_${decomp_np_half}/" "${hicar_repo}/tests/Test_Cases/restart/Repro_MPI_${decomp_np_half}/" "2017-02-14 00:05:00"

    run_hicar "Repro_MPI_${decomp_np_half}.nml" "$decomp_np_half" "MPI_${decomp_np_half}ranks"
    mpi_half_status=$?


    if [ $mpi_half_status -ne 0 ]; then
        echo -e "${RED}Domain Decomposition: FAILED (${decomp_np_half}-rank run failed)${NC}"
        test_decomp_result="FAIL"
    else
        # Generate and run with full ranks (5-min window — decomposition test does not restart)
        generate_standard_nml "${hicar_repo}/tests/Test_Cases/input/Repro_MPI_${decomp_np_full}.nml" "${hicar_repo}/tests/Test_Cases/output/Repro_MPI_${decomp_np_full}/" "${hicar_repo}/tests/Test_Cases/restart/Repro_MPI_${decomp_np_full}/" "2017-02-14 00:05:00"

        run_hicar "Repro_MPI_${decomp_np_full}.nml" "$decomp_np_full" "MPI_${decomp_np_full}ranks"
        mpi_full_status=$?


        if [ $mpi_full_status -ne 0 ]; then
            echo -e "${RED}Domain Decomposition: FAILED (${decomp_np_full}-rank run failed)${NC}"
            test_decomp_result="FAIL"
        else
            # Compare outputs
            echo
            echo "Comparing MPI outputs..."
            ${python_exe} ${compare_script} \
                "${hicar_repo}/tests/Test_Cases/output/Repro_MPI_${decomp_np_half}/${OUTPUT_FILENAME}" \
                "${hicar_repo}/tests/Test_Cases/output/Repro_MPI_${decomp_np_full}/${OUTPUT_FILENAME}" \
                --tolerance 0.0 \
                --figures-dir "${hicar_repo}/tests/Test_Cases/figures/decomposition"
            if [ $? -eq 0 ]; then
                test_decomp_result="PASS"
            else
                test_decomp_result="FAIL"
            fi
        fi
    fi

    fi # end skip check

fi

# ===============================================================
# Test: Restart Reproducibility
# ===============================================================
if [ "$test_mode" == "restart" ] || [ "$test_mode" == "all" ]; then
    echo
    echo "-------------------------------------------------------"
    echo -e "  ${BLUE}Restart Reproducibility${NC}"
    echo "-------------------------------------------------------"

    # Clear previous figures for this test
    rm -rf "${hicar_repo}/tests/Test_Cases/figures/restart"

    # Determine np for restart test
    if [ "$use_gpu" = true ]; then
        restart_np=$num_gpus
    else
        # CPU mode: use half of available cores (same logic as test_case_runner.sh)
        if [[ "$OSTYPE" == "darwin"* ]]; then
            total_np=$(sysctl -n hw.logicalcpu)
        else
            total_np=$(nproc --all)
        fi
        restart_np=$((total_np / 2))
        restart_np=$((restart_np > 2 ? restart_np : 2))
        restart_np=$((restart_np < 21 ? restart_np : 21))
    fi

    echo -e "Using ${restart_np} MPI ranks for restart test"

    # Create output/restart directories
    # Both runs share the same output/restart dirs (HICAR restart appends to existing output)
    rm -rf "${hicar_repo}/tests/Test_Cases/output/Repro_Restart" "${hicar_repo}/tests/Test_Cases/restart/Repro_Restart" "${hicar_repo}/tests/Test_Cases/output/Repro_Restart_continuous_backup"
    mkdir -p "${hicar_repo}/tests/Test_Cases/output/Repro_Restart" "${hicar_repo}/tests/Test_Cases/restart/Repro_Restart"

    # Run continuous (full 10 min)
    generate_standard_nml "${hicar_repo}/tests/Test_Cases/input/Repro_Restart_cont.nml" "${hicar_repo}/tests/Test_Cases/output/Repro_Restart/" "${hicar_repo}/tests/Test_Cases/restart/Repro_Restart/"

    run_hicar "Repro_Restart_cont.nml" "$restart_np" "Restart_continuous"
    cont_status=$?


    if [ $cont_status -ne 0 ]; then
        echo -e "${RED}Restart Reproducibility: FAILED (continuous run failed)${NC}"
        test_restart_result="FAIL"
    else
        # Save backup of continuous output before restart overwrites it
        mkdir -p "${hicar_repo}/tests/Test_Cases/output/Repro_Restart_continuous_backup"
        cp "${hicar_repo}/tests/Test_Cases/output/Repro_Restart/${OUTPUT_FILENAME}" "${hicar_repo}/tests/Test_Cases/output/Repro_Restart_continuous_backup/${OUTPUT_FILENAME}"

        # Run restart (picks up from checkpoint at 5 min)
        # Uses the same output and restart dirs (HICAR appends to existing output file)
        generate_restart_nml "${hicar_repo}/tests/Test_Cases/input/Repro_Restart_rst.nml" "${hicar_repo}/tests/Test_Cases/output/Repro_Restart/" "${hicar_repo}/tests/Test_Cases/restart/Repro_Restart/"

        run_hicar "Repro_Restart_rst.nml" "$restart_np" "Restart_from_checkpoint"
        rst_status=$?


        if [ $rst_status -ne 0 ]; then
            echo -e "${RED}Restart Reproducibility: FAILED (restart run failed)${NC}"
            test_restart_result="FAIL"
        else
            # Compare final timestep: continuous backup vs post-restart output
            echo
            echo "Comparing restart outputs (last timestep only)..."
            ${python_exe} ${compare_script} \
                "${hicar_repo}/tests/Test_Cases/output/Repro_Restart_continuous_backup/${OUTPUT_FILENAME}" \
                "${hicar_repo}/tests/Test_Cases/output/Repro_Restart/${OUTPUT_FILENAME}" \
                --tolerance 0.1 \
                --last-timestep-only \
                --figures-dir "${hicar_repo}/tests/Test_Cases/figures/restart"
            if [ $? -eq 0 ]; then
                test_restart_result="PASS"
            else
                test_restart_result="FAIL"
            fi
        fi
    fi
fi

# ===============================================================
# Summary
# ===============================================================
echo
echo "======================================================="
echo "  Reproducibility Test Summary"
echo "======================================================="

any_fail=0

if [ "$test_mode" == "decomposition" ] || [ "$test_mode" == "all" ]; then
    if [ "$test_decomp_result" == "PASS" ]; then
        echo -e "  Domain Decomposition (${decomp_np_half} vs ${decomp_np_full}):   ${GREEN}PASS${NC}"
    elif [ "$test_decomp_result" == "SKIP" ]; then
        echo -e "  Domain Decomposition:              ${BLUE}SKIP (insufficient GPUs)${NC}"
    else
        echo -e "  Domain Decomposition (${decomp_np_half} vs ${decomp_np_full}):   ${RED}${test_decomp_result}${NC}"
        any_fail=1
    fi

    if [ "$test_fm_decomp_result" == "PASS" ]; then
        echo -e "  Fine-Mesh Advection Decomposition: ${GREEN}PASS${NC}"
    elif [ "$test_fm_decomp_result" == "SKIP" ]; then
        echo -e "  Fine-Mesh Advection Decomposition: ${BLUE}SKIP (HICAR-tester not built)${NC}"
    else
        echo -e "  Fine-Mesh Advection Decomposition: ${RED}${test_fm_decomp_result}${NC}"
        any_fail=1
    fi
fi

if [ "$test_mode" == "restart" ] || [ "$test_mode" == "all" ]; then
    if [ "$test_restart_result" == "PASS" ]; then
        echo -e "  Restart Reproducibility:           ${GREEN}PASS${NC}"
    else
        echo -e "  Restart Reproducibility:           ${RED}${test_restart_result}${NC}"
        any_fail=1
    fi
fi

echo "======================================================="

if [ $any_fail -eq 0 ]; then
    echo -e "${GREEN}All reproducibility tests PASSED${NC}"
    exit 0
else
    echo -e "${RED}Some reproducibility tests FAILED${NC}"
    exit 1
fi
