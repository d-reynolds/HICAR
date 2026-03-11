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
if [ ! -f ./input/file_list_TestCase.txt ]; then
    "$hicar_repo/helpers/filelist_script.sh" "forcing/*" input/file_list_TestCase.txt
fi

# Copy supporting files if missing
if [ ! -f ./input/VEGPARM.TBL ]; then
    echo "Copying .TBL files to ./input"
    cp "$hicar_repo/run/"*.TBL ./input/
fi
if [ ! -d ./input/rrtmg_support ]; then
    cp -r "$hicar_repo/run/rrtmg_support" ./input
fi
if [ ! -d ./input/mp_support ]; then
    cp -r "$hicar_repo/run/mp_support" ./input
fi

if [ "$use_gpu" = true ]; then
    if [ ! -d ./input/rrtmgp_support ]; then
        cp -r "$hicar_repo/run/rrtmgp_support" ./input
    fi
fi
# Generate default namelist
default_file=input/default_hicar_options.nml
if [ -f "$default_file" ]; then
    rm "$default_file"
fi
echo "Generating default namelist to ./input"
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
    PY_ENV_PATH=$(pwd)/venv
    echo -e "Creating virtual environment at ${BLUE}${PY_ENV_PATH}${NC}"
    mkdir -p "$PY_ENV_PATH"
    $python_exe -m venv "${PY_ENV_PATH}"
    "${PY_ENV_PATH}/bin/pip" install numpy netCDF4 xarray matplotlib
    python_exe="${PY_ENV_PATH}/bin/python"
fi

export OMP_NUM_THREADS=1

# Cleanup trap for temp files
cleanup() {
    rm -f input/*.bak
}
trap cleanup EXIT

# -------------------------------------------------------
# Helper: generate_standard_nml <output_nml> <output_dir> <restart_dir>
# -------------------------------------------------------
generate_standard_nml() {
    local out_nml="$1"
    local output_dir="$2"
    local restart_dir="$3"

    cp "$default_file" "$out_nml"

    # Apply Standard.sh settings
    cd input
    ../input/nml_gen_scripts/Standard.sh "$(basename "$out_nml")"
    cd ..

    # Override wind, Sx, and output_vars for reproducibility testing
    sed -i'.bak' "s/ wind = 'variational solver'/ wind = 'none'/g" "$out_nml"
    sed -i'.bak' 's/Sx = .True./Sx = .False./g' "$out_nml"
    sed -i'.bak' "s/output_vars = .*$/output_vars = 'all'/g" "$out_nml"

    if [ $use_gpu = true ]; then
        # GPU mode: use all available GPUs and half for decomposition test, so set end time to 10 min for both runs
        sed -i'.bak' "s/rad = 'rrtmg'/rad = 'rrtmgp'/g" "$out_nml"
    fi
    # Shorten run: 10 min instead of 20, output every 5 min for restart checkpoints
    sed -i'.bak' "s/end_date = '2017-02-14 00:20:00'/end_date = '2017-02-14 01:10:00'/g" "$out_nml"
    sed -i'.bak' 's/outputinterval = 600/outputinterval = 300/g' "$out_nml"

    # Override output and restart folders
    sed -i'.bak' "s|output_folder = '../output/Standard/'|output_folder = '${output_dir}'|g" "$out_nml"
    sed -i'.bak' "s|restart_folder = '../restart/Standard/'|restart_folder = '${restart_dir}'|g" "$out_nml"
    rm -f "${out_nml}.bak"
}


# -------------------------------------------------------
# Helper: generate_restart_nml <output_nml> <output_dir> <restart_dir>
# -------------------------------------------------------
generate_restart_nml() {
    local out_nml="$1"
    local output_dir="$2"
    local restart_dir="$3"

    cp "$default_file" "$out_nml"

    # Apply Standard_restart.sh settings
    cd input
    ../input/nml_gen_scripts/Standard_restart.sh "$(basename "$out_nml")"
    cd ..

    # Override wind, Sx, and output_vars for reproducibility testing
    sed -i'.bak' "s/ wind = 'variational solver'/ wind = 'none'/g" "$out_nml"
    sed -i'.bak' 's/Sx = .True./Sx = .False./g' "$out_nml"
    sed -i'.bak' "s/output_vars = .*$/output_vars = 'all'/g" "$out_nml"

    if [ $use_gpu = true ]; then
        # GPU mode: use all available GPUs and half for decomposition test, so set end time to 10 min for both runs
        sed -i'.bak' "s/rad = 'rrtmg'/rad = 'rrtmgp'/g" "$out_nml"
    fi

    # Shorten run: end at 10 min, restart from 5 min checkpoint, output every 5 min
    sed -i'.bak' "s/end_date = '2017-02-14 00:20:00'/end_date = '2017-02-14 01:10:00'/g" "$out_nml"
    sed -i'.bak' "s/restart_date = '2017-02-14 00:10:00'/restart_date = '2017-02-14 01:00:00'/g" "$out_nml"
    #sed -i'.bak' 's/outputinterval = 600/outputinterval = 300/g' "$out_nml"

    # Override output and restart folders
    sed -i'.bak' "s|output_folder = '../output/Standard/'|output_folder = '${output_dir}'|g" "$out_nml"
    sed -i'.bak' "s|restart_folder = '../restart/Standard/'|restart_folder = '${restart_dir}'|g" "$out_nml"
    rm -f "${out_nml}.bak"
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

    return $status
}

# ===============================================================
# Test Results Tracking
# ===============================================================
test_decomp_result=""
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
        decomp_np_full=10
        decomp_np_half=5
    fi

    if [ "$test_decomp_result" != "SKIP" ]; then

    echo
    echo "-------------------------------------------------------"
    echo -e "  ${BLUE}Domain Decomposition Test (${decomp_np_half} vs ${decomp_np_full} ranks)${NC}"
    echo "-------------------------------------------------------"

    # Clear previous figures for this test
    rm -rf "${figures_dir}/decomposition"

    # Create output/restart directories
    for suffix in Repro_MPI_${decomp_np_half} Repro_MPI_${decomp_np_full}; do
        rm -rf "output/${suffix}" "restart/${suffix}"
        mkdir -p "output/${suffix}" "restart/${suffix}"
    done

    # Generate and run with half ranks
    generate_standard_nml "input/Repro_MPI_${decomp_np_half}.nml" "../output/Repro_MPI_${decomp_np_half}/" "../restart/Repro_MPI_${decomp_np_half}/"
    cd input
    run_hicar "Repro_MPI_${decomp_np_half}.nml" "$decomp_np_half" "MPI_${decomp_np_half}ranks"
    mpi_half_status=$?
    cd ..

    if [ $mpi_half_status -ne 0 ]; then
        echo -e "${RED}Domain Decomposition: FAILED (${decomp_np_half}-rank run failed)${NC}"
        test_decomp_result="FAIL"
    else
        # Generate and run with full ranks
        generate_standard_nml "input/Repro_MPI_${decomp_np_full}.nml" "../output/Repro_MPI_${decomp_np_full}/" "../restart/Repro_MPI_${decomp_np_full}/"
        cd input
        run_hicar "Repro_MPI_${decomp_np_full}.nml" "$decomp_np_full" "MPI_${decomp_np_full}ranks"
        mpi_full_status=$?
        cd ..

        if [ $mpi_full_status -ne 0 ]; then
            echo -e "${RED}Domain Decomposition: FAILED (${decomp_np_full}-rank run failed)${NC}"
            test_decomp_result="FAIL"
        else
            # Compare outputs
            echo
            echo "Comparing MPI outputs..."
            $python_exe "$compare_script" \
                "output/Repro_MPI_${decomp_np_half}/${OUTPUT_FILENAME}" \
                "output/Repro_MPI_${decomp_np_full}/${OUTPUT_FILENAME}" \
                --tolerance 0.0 \
                --figures-dir "${figures_dir}/decomposition"
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
    rm -rf "${figures_dir}/restart"

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
    rm -rf "output/Repro_Restart" "restart/Repro_Restart" "output/Repro_Restart_continuous_backup"
    mkdir -p "output/Repro_Restart" "restart/Repro_Restart"

    # Run continuous (full 10 min)
    generate_standard_nml "input/Repro_Restart_cont.nml" "../output/Repro_Restart/" "../restart/Repro_Restart/"
    cd input
    run_hicar "Repro_Restart_cont.nml" "$restart_np" "Restart_continuous"
    cont_status=$?
    cd ..

    if [ $cont_status -ne 0 ]; then
        echo -e "${RED}Restart Reproducibility: FAILED (continuous run failed)${NC}"
        test_restart_result="FAIL"
    else
        # Save backup of continuous output before restart overwrites it
        mkdir -p "output/Repro_Restart_continuous_backup"
        cp "output/Repro_Restart/${OUTPUT_FILENAME}" "output/Repro_Restart_continuous_backup/${OUTPUT_FILENAME}"

        # Run restart (picks up from checkpoint at 5 min)
        # Uses the same output and restart dirs (HICAR appends to existing output file)
        generate_restart_nml "input/Repro_Restart_rst.nml" "../output/Repro_Restart/" "../restart/Repro_Restart/"
        cd input
        run_hicar "Repro_Restart_rst.nml" "$restart_np" "Restart_from_checkpoint"
        rst_status=$?
        cd ..

        if [ $rst_status -ne 0 ]; then
            echo -e "${RED}Restart Reproducibility: FAILED (restart run failed)${NC}"
            test_restart_result="FAIL"
        else
            # Compare final timestep: continuous backup vs post-restart output
            echo
            echo "Comparing restart outputs (last timestep only)..."
            $python_exe "$compare_script" \
                "output/Repro_Restart_continuous_backup/${OUTPUT_FILENAME}" \
                "output/Repro_Restart/${OUTPUT_FILENAME}" \
                --tolerance 0.0 \
                --last-timestep-only \
                --figures-dir "${figures_dir}/restart"
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
