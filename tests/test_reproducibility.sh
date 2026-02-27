#!/bin/bash

# Reproducibility tests for HICAR:
#   decomposition: MPI Reproducibility — run with 5 vs 10 ranks, compare output (must be identical)
#   restart:       Restart Reproducibility — continuous run vs restart run, compare final timestep
#   all:           Run both tests
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

# Get the location of a HICAR executable (debug build required for reproducibility)
if [ -f "$hicar_repo/bin/HICAR_debug" ]; then
    hicar_exe="$hicar_repo/bin/HICAR_debug"
else
    echo -e "${RED}HICAR_debug executable not found in bin directory.${NC}"
    echo -e "${RED}Build with MODE=debug and run 'make install' from the build directory.${NC}"
    exit 1
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

    # Shorten run: 10 min instead of 20, output every 5 min for restart checkpoints
    sed -i'.bak' "s/end_date = '2017-02-14 00:20:00'/end_date = '2017-02-14 00:10:00'/g" "$out_nml"
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

    # Shorten run: end at 10 min, restart from 5 min checkpoint, output every 5 min
    sed -i'.bak' "s/end_date = '2017-02-14 00:20:00'/end_date = '2017-02-14 00:10:00'/g" "$out_nml"
    sed -i'.bak' "s/restart_date = '2017-02-14 00:10:00'/restart_date = '2017-02-14 00:05:00'/g" "$out_nml"
    sed -i'.bak' 's/outputinterval = 600/outputinterval = 300/g' "$out_nml"

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
test_decomp_result="SKIP"
test_restart_result="SKIP"

echo
echo "======================================================="
echo "  HICAR Reproducibility Tests (mode: ${test_mode})"
echo "======================================================="

# ===============================================================
# Test: MPI Domain Decomposition (5 ranks vs 10 ranks)
# ===============================================================
if [ "$test_mode" == "decomposition" ] || [ "$test_mode" == "all" ]; then
    echo
    echo "-------------------------------------------------------"
    echo -e "  ${BLUE}Domain Decomposition Test (5 vs 10 ranks)${NC}"
    echo "-------------------------------------------------------"

    # Clear previous figures for this test
    rm -rf "${figures_dir}/decomposition"

    # Create output/restart directories
    for suffix in Repro_MPI_5 Repro_MPI_10; do
        rm -rf "output/${suffix}" "restart/${suffix}"
        mkdir -p "output/${suffix}" "restart/${suffix}"
    done

    # Generate and run with 5 ranks
    generate_standard_nml "input/Repro_MPI_5.nml" "../output/Repro_MPI_5/" "../restart/Repro_MPI_5/"
    cd input
    run_hicar "Repro_MPI_5.nml" 5 "MPI_5ranks"
    mpi5_status=$?
    cd ..

    if [ $mpi5_status -ne 0 ]; then
        echo -e "${RED}Domain Decomposition: FAILED (5-rank run failed)${NC}"
        test_decomp_result="FAIL"
    else
        # Generate and run with 10 ranks
        generate_standard_nml "input/Repro_MPI_10.nml" "../output/Repro_MPI_10/" "../restart/Repro_MPI_10/"
        cd input
        run_hicar "Repro_MPI_10.nml" 10 "MPI_10ranks"
        mpi10_status=$?
        cd ..

        if [ $mpi10_status -ne 0 ]; then
            echo -e "${RED}Domain Decomposition: FAILED (10-rank run failed)${NC}"
            test_decomp_result="FAIL"
        else
            # Compare outputs
            echo
            echo "Comparing MPI outputs..."
            $python_exe "$compare_script" \
                "output/Repro_MPI_5/${OUTPUT_FILENAME}" \
                "output/Repro_MPI_10/${OUTPUT_FILENAME}" \
                --tolerance 0.0 \
                --figures-dir "${figures_dir}/decomposition"
            if [ $? -eq 0 ]; then
                test_decomp_result="PASS"
            else
                test_decomp_result="FAIL"
            fi
        fi
    fi
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

    # Determine np for restart test (same logic as test_case_runner.sh)
    if [[ "$OSTYPE" == "darwin"* ]]; then
        total_np=$(sysctl -n hw.logicalcpu)
    else
        total_np=$(nproc --all)
    fi
    restart_np=$((total_np / 2))
    restart_np=$((restart_np > 2 ? restart_np : 2))
    restart_np=$((restart_np < 21 ? restart_np : 21))

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
        echo -e "  Domain Decomposition (5 vs 10):   ${GREEN}PASS${NC}"
    else
        echo -e "  Domain Decomposition (5 vs 10):   ${RED}${test_decomp_result}${NC}"
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
