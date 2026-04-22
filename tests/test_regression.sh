#!/bin/bash

# Regression test for HICAR:
#   Builds the executable from a known-good commit, runs both old and new
#   executables with output_vars='all', and compares outputs using compare_outputs.py.
#   Each executable generates its own namelist via --gen-nml so defaults match each version.
#
# Usage: test_regression.sh <hicar_repo> <build_dir> [tolerance]
#   $1 = path to HICAR repo
#   $2 = path to build directory (used to read CMakeCache.txt for build config)
#   $3 = optional tolerance for comparison (default 0.0)

# Define colors for output styling
BLUE='\033[0;36m'
GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m' # No Color

set -euo pipefail

# -------------------------------------------------------
# Parse arguments
# -------------------------------------------------------
if [ $# -lt 2 ]; then
    echo "Usage: $0 <hicar_repo> <build_dir> [tolerance]"
    exit 1
fi

hicar_repo=$(cd "$1" && pwd)
BUILD_DIR=$(cd "$2" && pwd)
TOLERANCE="${3:-0.0}"

REGRESSION_HASH_FILE="$hicar_repo/tests/regression_commit.txt"
COMPARE_SCRIPT="$hicar_repo/tests/compare_outputs.py"
FIGURES_DIR="$hicar_repo/tests/figures/regression"

# Output file name pattern (from Standard.sh: domain=Gaudergrat_250m, start=2017-02-14)
OUTPUT_FILENAME="Gaudergrat_250m_2017-02-14_00-00-00.nc"

# Read the known-good commit hash
if [ ! -f "$REGRESSION_HASH_FILE" ]; then
    echo -e "${RED}regression_commit.txt not found at ${REGRESSION_HASH_FILE}${NC}"
    exit 1
fi
REF_HASH=$(head -n1 "$REGRESSION_HASH_FILE" | tr -d '[:space:]')
if [ -z "$REF_HASH" ]; then
    echo -e "${RED}regression_commit.txt is empty${NC}"
    exit 1
fi

# Validate that REF_HASH is a real commit
if ! git -C "$hicar_repo" cat-file -t "$REF_HASH" &>/dev/null; then
    echo -e "${RED}REF_HASH '${REF_HASH}' is not a valid commit in this repo${NC}"
    exit 1
fi

# Guard against REF_HASH pointing to HEAD (would trivially pass)
CURRENT_HEAD=$(git -C "$hicar_repo" rev-parse HEAD)
REF_FULL=$(git -C "$hicar_repo" rev-parse "$REF_HASH")
if [ "$REF_FULL" = "$CURRENT_HEAD" ]; then
    echo -e "${RED}ERROR: regression_commit.txt points to HEAD — test would trivially pass.${NC}"
    echo -e "${RED}The hash must point to a previous known-good commit.${NC}"
    exit 1
fi

echo
echo "======================================================="
echo "  HICAR Regression Test"
echo "======================================================="
echo -e "  Reference commit: ${BLUE}${REF_HASH}${NC}"
echo -e "  Current HEAD:     ${BLUE}$(git -C "$hicar_repo" rev-parse HEAD)${NC}"
echo -e "  Tolerance:        ${TOLERANCE}"
echo "======================================================="

# -------------------------------------------------------
# Locate current (new) executable
# -------------------------------------------------------
if [ -f "$hicar_repo/bin/HICAR_debug_gpu" ]; then
    new_exe="$hicar_repo/bin/HICAR_debug_gpu"
elif [ -f "$hicar_repo/bin/HICAR_debug" ]; then
    new_exe="$hicar_repo/bin/HICAR_debug"
elif [ -f "$hicar_repo/bin/HICAR" ]; then
    new_exe="$hicar_repo/bin/HICAR"
else
    echo -e "${RED}No HICAR executable found in bin directory.${NC}"
    exit 1
fi
echo -e "  New executable:   ${BLUE}${new_exe}${NC}"

# -------------------------------------------------------
# Build Phase: build old executable from reference commit
# -------------------------------------------------------
OLD_EXE="$hicar_repo/bin/HICAR_old_tester"
OLD_HASH_FILE="$hicar_repo/bin/.regression_hash"

# Worktree cleanup trap
WORKTREE_DIR=""
cleanup_worktree() {
    if [ -n "$WORKTREE_DIR" ] && [ -d "$WORKTREE_DIR" ]; then
        echo "Cleaning up worktree at $WORKTREE_DIR..."
        git -C "$hicar_repo" worktree remove "$WORKTREE_DIR" --force 2>/dev/null || true
        rm -rf "$WORKTREE_DIR" 2>/dev/null || true
    fi
}
trap cleanup_worktree EXIT ERR INT TERM

build_needed=true
if [ -f "$OLD_EXE" ] && [ -f "$OLD_HASH_FILE" ]; then
    cached_hash=$(head -n1 "$OLD_HASH_FILE" | tr -d '[:space:]')
    if [ "$cached_hash" = "$REF_HASH" ]; then
        echo -e "${GREEN}Cached old executable matches reference hash — skipping build${NC}"
        build_needed=false
    fi
fi

if [ "$build_needed" = true ]; then
    echo
    echo "-------------------------------------------------------"
    echo -e "  ${BLUE}Building old executable from ${REF_HASH:0:12}...${NC}"
    echo "-------------------------------------------------------"

    # Ensure the commit exists locally (handles shallow clones)
    git -C "$hicar_repo" fetch origin "$REF_HASH" 2>/dev/null || true

    # Create worktree in /tmp (completely outside the repo)
    WORKTREE_DIR=$(mktemp -d "${TMPDIR:-/tmp}/hicar_regression_XXXXXX")
    # mktemp creates the directory; remove it so git worktree add can use the path
    rmdir "$WORKTREE_DIR"

    # Safety check
    if [ "$WORKTREE_DIR" = "$hicar_repo" ]; then
        echo -e "${RED}FATAL: worktree path equals repo path. Aborting.${NC}"
        exit 1
    fi

    echo "  Worktree location: $WORKTREE_DIR"
    git -C "$hicar_repo" worktree add "$WORKTREE_DIR" "$REF_HASH" --detach

    # Read build configuration from the current build's CMakeCache.txt
    cmake_args=()
    if [ -f "$BUILD_DIR/CMakeCache.txt" ]; then
        # Extract key variables from CMakeCache
        for var in MODE FC PETSC_DIR PETSC_ARCH AMGX_DIR FFTW_DIR NETCDF_DIR MPI_DIR FSM_DIR FSM SNOWPACK_DIR SNOWPACK ASSERTIONS OPENACC; do
            val=$(grep "^${var}:" "$BUILD_DIR/CMakeCache.txt" 2>/dev/null | head -n1 | sed 's/^[^=]*=//')
            if [ -n "$val" ]; then
                cmake_args+=("-D${var}=${val}")
            fi
        done
    else
        echo -e "${RED}Warning: CMakeCache.txt not found at ${BUILD_DIR}/CMakeCache.txt${NC}"
        echo -e "${RED}Building old executable with default settings${NC}"
    fi

    # Build the old executable
    OLD_BUILD_DIR="$WORKTREE_DIR/build_regression"
    mkdir -p "$OLD_BUILD_DIR"

    echo "  Configuring..."
    if ! cmake -S "$WORKTREE_DIR" -B "$OLD_BUILD_DIR" "${cmake_args[@]}" 2>&1 | tail -n 5; then
        echo -e "${RED}CMake configuration failed for old commit. The reference commit may not be buildable.${NC}"
        echo -e "${RED}Consider updating regression_commit.txt to a newer known-good commit.${NC}"
        exit 1
    fi

    echo "  Building..."
    if ! cmake --build "$OLD_BUILD_DIR" --parallel 2>&1 | tail -n 5; then
        echo -e "${RED}Build failed for old commit ${REF_HASH:0:12}.${NC}"
        echo -e "${RED}Consider updating regression_commit.txt to a newer known-good commit.${NC}"
        exit 1
    fi

    echo "  Installing..."
    cmake --install "$OLD_BUILD_DIR" 2>&1 | tail -n 3

    # Find the built executable in the worktree's bin/
    old_built=""
    for name in HICAR_debug_gpu HICAR_debug HICAR HICAR_gpu; do
        if [ -f "$WORKTREE_DIR/bin/$name" ]; then
            old_built="$WORKTREE_DIR/bin/$name"
            break
        fi
    done

    if [ -z "$old_built" ]; then
        echo -e "${RED}No executable found in worktree bin/ after build.${NC}"
        exit 1
    fi

    # Cache the old executable
    mkdir -p "$hicar_repo/bin"
    cp "$old_built" "$OLD_EXE"
    echo "$REF_HASH" > "$OLD_HASH_FILE"
    echo -e "${GREEN}Old executable cached at ${OLD_EXE}${NC}"

    # Cleanup worktree (trap will also try, but do it explicitly)
    git -C "$hicar_repo" worktree remove "$WORKTREE_DIR" --force 2>/dev/null || true
    WORKTREE_DIR=""
fi

echo -e "  Old executable:   ${BLUE}${OLD_EXE}${NC}"

# -------------------------------------------------------
# Setup Phase: copy support files, detect MPI, setup Python
# -------------------------------------------------------

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

# Detect GPU mode
use_gpu=false
if [[ "$new_exe" == *"gpu"* ]]; then
    use_gpu=true
    if [ ! -d ./input/rrtmgp_support ]; then
        cp -r "$hicar_repo/run/rrtmgp_support" ./input
    fi
fi

# Make forcing file list
if [ ! -f ./input/file_list_TestCase.txt ]; then
    "$hicar_repo/helpers/filelist_script.sh" "forcing/*" input/file_list_TestCase.txt
fi

# Find mpiexec (same logic as test_reproducibility.sh)
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

# Determine number of MPI ranks
num_gpus=0
if [ "$use_gpu" = true ]; then
    if command -v nvidia-smi &> /dev/null; then
        num_gpus=$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | wc -l | tr -d ' ')
    fi
    if [ "$num_gpus" -lt 1 ]; then
        echo -e "${RED}GPU mode selected but no GPUs detected via nvidia-smi.${NC}"
        exit 1
    fi
    np=$((num_gpus + 1))
else
    if [[ "$OSTYPE" == "darwin"* ]]; then
        total_np=$(sysctl -n hw.logicalcpu)
    else
        total_np=$(nproc --all)
    fi
    np=$((total_np / 2))
    np=$((np > 2 ? np : 2))
    np=$((np < 21 ? np : 21))
fi
export OMP_NUM_THREADS=1

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

# -------------------------------------------------------
# Helper: run_hicar <exe> <nml_file> <np> <label>
# -------------------------------------------------------
run_hicar() {
    local exe="$1"
    local nml_file="$2"
    local run_np="$3"
    local label="$4"

    echo
    echo -e "Running: ${BLUE}${label}${NC} with ${run_np} MPI ranks"

    local base_label
    base_label=$(echo "$label" | tr ' ' '_')

    if [ -n "$mpiexec_path" ]; then
        $mpiexec_path -np "$run_np" "$exe" "$nml_file" 1>"${base_label}.out" 2>"${base_label}.err" &
        hicar_pid=$!
    elif command -v srun &> /dev/null; then
        srun -N 1 -n "$run_np" $SRUN_FLAGS "$exe" "$nml_file" 1>"${base_label}.out" 2>"${base_label}.err" &
        hicar_pid=$!
    else
        echo -e "${RED}No MPI launcher available${NC}"
        return 1
    fi

    echo -n "Initializing..."

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

# -------------------------------------------------------
# Helper: apply_standard_sed <nml_file> <output_dir> <restart_dir>
#   Applies Standard.sh sed patterns + output_vars='all' overrides
# -------------------------------------------------------
apply_standard_sed() {
    local out_file="$1"
    local output_dir="$2"
    local restart_dir="$3"

    sed -i'.bak' "s/start_date = ''/start_date = '2017-02-14 00:00:00'/g" "$out_file"
    sed -i'.bak' "s/end_date = ''/end_date = '2017-02-14 00:20:00'/g" "$out_file"
    sed -i'.bak' "s/restartinterval = 24/restartinterval = 1/g" "$out_file"
    sed -i'.bak' "s/outputinterval = 3600/outputinterval = 600/g" "$out_file"
    sed -i'.bak' "s|output_folder = '../output/'|output_folder = '${output_dir}'|g" "$out_file"
    sed -i'.bak' "s|restart_folder = '../restart/'|restart_folder = '${restart_dir}'|g" "$out_file"
    sed -i'.bak' 's/dx = 0.0/dx = 250.0/g' "$out_file"
    sed -i'.bak' 's/nz = 500/nz = 40/g' "$out_file"
    sed -i'.bak' "s/ pbl = 'none'/ pbl = 'ysu'/g" "$out_file"
    sed -i'.bak' "s/ lsm = 'none'/ lsm = 'noahmp'/g" "$out_file"
    sed -i'.bak' "s/ sfc = 'none'/ sfc = 'revmm5'/g" "$out_file"
    sed -i'.bak' "s/ water = 'none'/ water = 'simple'/g" "$out_file"
    sed -i'.bak' "s/ mp = 'none'/ mp = 'morrison'/g" "$out_file"
    sed -i'.bak' "s/ rad = 'none'/ rad = 'rrtmg'/g" "$out_file"
    sed -i'.bak' 's/terrain_shading = .False./terrain_shading = .True./g' "$out_file"
    sed -i'.bak' 's/Sx = .False./Sx = .True./g' "$out_file"
    # Force output_vars = 'all' for regression testing
    sed -i'.bak' "s/output_vars = ''/output_vars = 'all'/g" "$out_file"
    sed -i'.bak' 's/inputinterval = 3600/inputinterval = 3600/g' "$out_file"
    sed -i'.bak' "s/LU_Categories = 'MODIFIED_IGBP_MODIS_NOAH'/LU_Categories = 'USGS'/g" "$out_file"
    sed -i'.bak' "s/ time_var = ''/ time_var = 'time'/g" "$out_file"
    sed -i'.bak' "s/ pvar = ''/ pvar = 'P'/g" "$out_file"
    sed -i'.bak' "s/ tvar = ''/ tvar = 'T'/g" "$out_file"
    sed -i'.bak' "s/ qvvar = ''/ qvvar = 'QV'/g" "$out_file"
    sed -i'.bak' "s/ hgtvar = ''/ hgtvar = 'HSURF'/g" "$out_file"
    sed -i'.bak' "s/ zvar = ''/ zvar = 'HFL'/g" "$out_file"
    sed -i'.bak' "s/ latvar = ''/ latvar = 'lat_1'/g" "$out_file"
    sed -i'.bak' "s/ lonvar = ''/ lonvar = 'lon_1'/g" "$out_file"
    sed -i'.bak' "s/ uvar = ''/ uvar = 'U'/g" "$out_file"
    sed -i'.bak' "s/ vvar = ''/ vvar = 'V'/g" "$out_file"
    sed -i'.bak' "s/ lat_hi = ''/ lat_hi = 'lat'/g" "$out_file"
    sed -i'.bak' "s/ lon_hi = ''/ lon_hi = 'lon'/g" "$out_file"
    sed -i'.bak' "s/ hgt_hi = ''/ hgt_hi = 'topo'/g" "$out_file"
    sed -i'.bak' "s/ landvar = ''/ landvar = 'landmask'/g" "$out_file"
    sed -i'.bak' "s/ vegtype_var = ''/ vegtype_var = 'landuse'/g" "$out_file"
    sed -i'.bak' "s/ hlm_var = ''/ hlm_var = 'hlm'/g" "$out_file"
    sed -i'.bak' "s/ svf_var = ''/ svf_var = 'svf'/g" "$out_file"
    sed -i'.bak' "s/ slope_var = ''/ slope_var = 'slope'/g" "$out_file"
    sed -i'.bak' "s/ slope_angle_var = ''/ slope_angle_var = 'slope_rad'/g" "$out_file"
    sed -i'.bak' "s/ aspect_angle_var = ''/ aspect_angle_var = 'aspect_rad'/g" "$out_file"
    sed -i'.bak' 's/smooth_wind_distance = -9999/smooth_wind_distance = 500.0/g' "$out_file"
    sed -i'.bak' "s/init_conditions_file = ''/init_conditions_file = '..\/domains\/Gaudergrat_250m.nc'/g" "$out_file"
    sed -i'.bak' "s/forcing_file_list = ''/forcing_file_list = 'file_list_TestCase.txt'/g" "$out_file"
    sed -i'.bak' 's/qv_is_spec_humidity = .False./qv_is_spec_humidity = .True./g' "$out_file"
    sed -i'.bak' 's/t_is_potential = .True./t_is_potential = .False./g' "$out_file"
    sed -i'.bak' 's/tzone = 0.0/tzone = 1.0/g' "$out_file"
    sed -i'.bak' 's/update_interval_rrtmg = 600/update_interval_rrtmg = 600/g' "$out_file"
    sed -i'.bak' 's/qv_is_spec_humidity = .True./qv_is_spec_humidity = .True./g' "$out_file"
    sed -i'.bak' 's/cfl_reduction_factor = 0.9/cfl_reduction_factor = 1.6/g' "$out_file"
    sed -i'.bak' 's/RK3 = .False./RK3 = .True./g' "$out_file"
    sed -i'.bak' 's/flux_corr = 0/flux_corr = 1/g' "$out_file"
    sed -i'.bak' 's/h_order = 1/h_order = 3/g' "$out_file"
    sed -i'.bak' 's/v_order = 1/v_order = 3/g' "$out_file"
    sed -i'.bak' 's/sleve = .False./sleve = .True./g' "$out_file"
    sed -i'.bak' 's/terrain_smooth_windowsize = 3/terrain_smooth_windowsize = 5/g' "$out_file"
    sed -i'.bak' 's/terrain_smooth_cycles = 5/terrain_smooth_cycles = 100/g' "$out_file"
    sed -i'.bak' 's/sleve_n = 1.2/sleve_n = 1.35/g' "$out_file"
    sed -i'.bak' 's/decay_rate_L_topo = 2/decay_rate_L_topo = 1/g' "$out_file"
    sed -i'.bak' 's/decay_rate_S_topo = 6/decay_rate_S_topo = 3/g' "$out_file"
    sed -i'.bak' 's/dz_levels = 0.0/dz_levels = 23.0793700550371, 25.4054789829865, 27.9546004029731, 30.7456617662627, 33.7986658963874, 37.13462986434, 40.7754878622213, 44.7439504647213, 49.0633117863096, 53.7571952388465, 58.8492279502517, 64.3626335260366, 70.3197328568148, 76.7413432615544, 83.6460676036046, 91.0494673453666, 98.9631170549111, 107.393542879556, 116.341054166015, 125.798485877882, 135.749879770318, 146.169144288172, 157.018746480945, 168.248503199185, 179.79455242322, 191.578597373964, 203.507524314657, 215.473497662955, 227.354631108753, 239.016318976301, 250.313286722599, 261.092382801852, 271.196087038084, 280.46665561911, 288.750764024197, 295.904452287969, 301.798128624738, 306.321354381736, 309.387121454553, 310.935346604612/g' "$out_file"

    if [ "$use_gpu" = true ]; then
        sed -i'.bak' "s/rad = 'rrtmg'/rad = 'rrtmgp'/g" "$out_file"
    fi

    rm -f "${out_file}.bak"
}

# -------------------------------------------------------
# Run Phase
# -------------------------------------------------------
echo
echo "-------------------------------------------------------"
echo -e "  ${BLUE}Generating namelists and running executables${NC}"
echo "-------------------------------------------------------"

# Clear previous regression output/figures
rm -rf output/Regression_old output/Regression_new
rm -rf restart/Regression_old restart/Regression_new
rm -rf "${FIGURES_DIR}"
mkdir -p output/Regression_old output/Regression_new
mkdir -p restart/Regression_old restart/Regression_new

# Generate OLD namelist using old executable's --gen-nml
old_nml="input/Regression_old.nml"
if [ -f "$old_nml" ]; then rm "$old_nml"; fi
echo "Generating old namelist..."
"$OLD_EXE" --gen-nml "$old_nml"
apply_standard_sed "$old_nml" "../output/Regression_old/" "../restart/Regression_old/"

# Generate NEW namelist using new executable's --gen-nml
new_nml="input/Regression_new.nml"
if [ -f "$new_nml" ]; then rm "$new_nml"; fi
echo "Generating new namelist..."
"$new_exe" --gen-nml "$new_nml"
apply_standard_sed "$new_nml" "../output/Regression_new/" "../restart/Regression_new/"

# Run old executable
cd input
run_hicar "$OLD_EXE" "Regression_old.nml" "$np" "Regression_old"
old_status=$?
cd ..

if [ $old_status -ne 0 ]; then
    echo -e "${RED}Old executable run failed. Cannot compare.${NC}"
    exit 1
fi

# Run new executable
cd input
run_hicar "$new_exe" "Regression_new.nml" "$np" "Regression_new"
new_status=$?
cd ..

if [ $new_status -ne 0 ]; then
    echo -e "${RED}New executable run failed. Cannot compare.${NC}"
    exit 1
fi

# -------------------------------------------------------
# Compare Phase
# -------------------------------------------------------
echo
echo "-------------------------------------------------------"
echo -e "  ${BLUE}Comparing outputs${NC}"
echo "-------------------------------------------------------"

$python_exe "$COMPARE_SCRIPT" \
    "output/Regression_old/${OUTPUT_FILENAME}" \
    "output/Regression_new/${OUTPUT_FILENAME}" \
    --tolerance "$TOLERANCE" \
    --figures-dir "${FIGURES_DIR}"
compare_status=$?

# -------------------------------------------------------
# Auto-Update Phase
# -------------------------------------------------------
echo
echo "======================================================="
if [ $compare_status -eq 0 ]; then
    CURRENT_HEAD=$(git -C "$hicar_repo" rev-parse HEAD)
    echo -e "${GREEN}Regression test PASSED.${NC}"
    if [ "${AUTO_UPDATE:-false}" = "true" ]; then
        echo "$CURRENT_HEAD" > "$REGRESSION_HASH_FILE"
        echo -e "${GREEN}Updated regression_commit.txt to ${CURRENT_HEAD}.${NC}"
        echo -e "${GREEN}Stage and commit to persist.${NC}"
    fi
    exit 0
else
    echo -e "${RED}Regression test FAILED.${NC}"
    echo -e "${RED}regression_commit.txt NOT updated.${NC}"
    if [ -d "${FIGURES_DIR}" ]; then
        echo -e "  Comparison figures saved to: ${BLUE}${FIGURES_DIR}${NC}"
    fi
    exit 1
fi
