#!/bin/bash
# Valgrind memcheck of the HICAR unit-test suite (CPU debug build).
#
# Runs HICAR-tester under `valgrind --track-origins` with 4 MPI ranks and FAILS
# when any valgrind error block cites a HICAR Fortran source (<name>.F90:<line>)
# — the real signal, robust to mpich/libc/ld noise — or when the tester itself
# exits non-zero. The full log is written to <build_dir>/valgrind.log.
#
# This is the CPU memcheck CI lane (.github/workflows/valgrind-memcheck.yml) made
# runnable locally. It needs a debug build (-DMODE=debug) and `valgrind` on PATH.
#
# Usage:
#   bash tests/scripts/test_valgrind.sh <hicar_repo> [build_dir]   # build_dir default <repo>/build
#   make test_valgrind                                      # via CMake (debug build)

set -uo pipefail

BLUE='\033[0;36m'; GREEN='\033[0;32m'; RED='\033[0;31m'; NC='\033[0m'

# --- Resolve paths --------------------------------------------------------
if [ -z "${1:-}" ]; then
    echo -e "${RED}Usage: test_valgrind.sh <hicar_repo> [build_dir]${NC}" >&2
    exit 2
fi
hicar_repo=$(cd "$1" 2>/dev/null && pwd) || { echo -e "${RED}repo not found: $1${NC}" >&2; exit 2; }
build_dir=$(cd "${2:-$hicar_repo/build}" 2>/dev/null && pwd) \
    || { echo -e "${RED}build dir not found: ${2:-$hicar_repo/build}${NC}" >&2; exit 2; }
tester="$build_dir/tests/HICAR-tester"
log="$build_dir/valgrind.log"

# --- Preconditions --------------------------------------------------------
if ! command -v valgrind >/dev/null 2>&1; then
    echo -e "${RED}valgrind not found on PATH — install it (e.g. 'apt-get install valgrind') and retry.${NC}" >&2
    exit 2
fi
if [ ! -x "$tester" ]; then
    echo -e "${RED}HICAR-tester not found at $tester — build it first ('make HICAR-tester' or 'make test_valgrind').${NC}" >&2
    exit 2
fi
# Locate a real mpiexec (skip conda/python shims that shadow the MPI launcher).
mpiexec_path=""
IFS=':' read -ra _PD <<< "$PATH"
for d in "${_PD[@]}"; do
    if [ -x "$d/mpiexec" ] && ! [[ "$d" =~ python|conda ]]; then mpiexec_path="$d/mpiexec"; break; fi
done
[ -z "$mpiexec_path" ] && { echo -e "${RED}mpiexec not found on PATH${NC}" >&2; exit 2; }

echo -e "${BLUE}Running HICAR-tester under valgrind (mpiexec -np 4): $tester${NC}"

# --- Run + post-process ---------------------------------------------------
# Capture the tester's exit code and post-process the log ourselves. No
# --error-exitcode: mpich isn't valgrind-clean, so we gate on whether any error
# block cites a HICAR .F90 (the real signal, robust to mpich/libc/ld noise).
# --track-origins reports where the value was created.
"$mpiexec_path" -np 4 valgrind -q \
    --track-origins=yes \
    --leak-check=no --errors-for-leak-kinds=none \
    --num-callers=30 \
    "$tester" 2>&1 | tee "$log"
rc=${PIPESTATUS[0]}

fail=0
echo "----- scanning valgrind output for errors in HICAR Fortran sources -----"
if grep -qE "[A-Za-z0-9_]+\.F90:[0-9]+" "$log"; then
    echo "::error::valgrind reported memory errors citing HICAR Fortran sources (see artifact)"
    grep -nE "depends on uninitialised|Use of uninitialised|Invalid (read|write)|Conditional jump|[A-Za-z0-9_]+\.F90:[0-9]+|was created by" "$log" || true
    fail=1
else
    echo -e "${GREEN}valgrind clean: no error blocks cite a HICAR .F90 source.${NC}"
fi
if [ "$rc" -ne 0 ]; then
    echo "::error::HICAR-tester exited non-zero ($rc) under valgrind (crash or failed unit test)"
    fail=1
fi
exit $fail
