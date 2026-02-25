#!/bin/bash
#
# CLI Options Test Suite for HICAR
#
# Validates that all command-line options in driver.F90 produce the expected
# exit codes and output. Run with:
#   bash tests/test_cli_options.sh /path/to/hicar/repo
#
# Or via cmake:
#   make test_cli

# Colors for output
BLUE='\033[0;36m'
GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m'

# Counters
PASS=0
FAIL=0
TOTAL=0

# Get HICAR repo root (macOS-compatible, no readlink -f)
if [ -n "$1" ]; then
    hicar_repo=$(cd "$1" 2>/dev/null && pwd)
else
    # Default: assume script is in tests/ directory
    hicar_repo=$(cd "$(dirname "$0")/.." 2>/dev/null && pwd)
fi

if [ -z "$hicar_repo" ] || [ ! -d "$hicar_repo" ]; then
    echo -e "${RED}Error: Could not resolve HICAR repo path${NC}"
    echo "Usage: $0 /path/to/hicar/repo"
    exit 1
fi

# Find HICAR executable
if [ -f "$hicar_repo/bin/HICAR" ]; then
    HICAR_EXE="$hicar_repo/bin/HICAR"
elif [ -f "$hicar_repo/bin/HICAR_debug" ]; then
    HICAR_EXE="$hicar_repo/bin/HICAR_debug"
elif [ -f "$hicar_repo/bin/HICAR_gpu" ]; then
    HICAR_EXE="$hicar_repo/bin/HICAR_gpu"
else
    echo -e "${RED}Error: HICAR executable not found in $hicar_repo/bin/${NC}"
    echo "Ensure the project is built and 'make install' was run."
    exit 1
fi

echo "============================================"
echo "  HICAR CLI Options Test Suite"
echo "============================================"
echo "Executable: $HICAR_EXE"
echo ""

# Create temp directory for test artifacts
TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Helper function to run a single test
# Usage: run_test "description" expected_exit_ok check_type check_value [extra_args...]
#   expected_exit_ok: "yes" = must be exit 0, "any" = any exit code is fine
#   check_type: "stdout_contains", "combined_contains", "stdout_nonempty",
#               "stdout_lines_gt", "file_exists_nonempty"
#   check_value: string to search for, line count threshold, or file path
run_test() {
    local description="$1"
    local expected_exit_ok="$2"
    local check_type="$3"
    local check_value="$4"
    shift 4
    local cmd_args=("$@")

    TOTAL=$((TOTAL + 1))
    local stdout_file="$TMPDIR/stdout_${TOTAL}"
    local stderr_file="$TMPDIR/stderr_${TOTAL}"

    # Run the command
    "$HICAR_EXE" "${cmd_args[@]}" >"$stdout_file" 2>"$stderr_file"
    local exit_code=$?

    local failed=0
    local reason=""

    # Check exit code
    if [ "$expected_exit_ok" = "yes" ] && [ $exit_code -ne 0 ]; then
        failed=1
        reason="expected exit 0, got $exit_code"
    fi

    # Check output based on type
    if [ $failed -eq 0 ]; then
        case "$check_type" in
            stdout_contains)
                if ! grep -q "$check_value" "$stdout_file"; then
                    failed=1
                    reason="stdout does not contain '$check_value'"
                fi
                ;;
            combined_contains)
                if ! grep -q "$check_value" "$stdout_file" && ! grep -q "$check_value" "$stderr_file"; then
                    failed=1
                    reason="neither stdout nor stderr contains '$check_value'"
                fi
                ;;
            stdout_nonempty)
                if [ ! -s "$stdout_file" ]; then
                    failed=1
                    reason="stdout is empty"
                fi
                ;;
            stdout_lines_gt)
                local line_count
                line_count=$(wc -l < "$stdout_file" | tr -d ' ')
                if [ "$line_count" -le "$check_value" ]; then
                    failed=1
                    reason="stdout has $line_count lines, expected more than $check_value"
                fi
                ;;
            file_exists_nonempty)
                if [ ! -s "$check_value" ]; then
                    failed=1
                    reason="file '$check_value' does not exist or is empty"
                fi
                ;;
        esac
    fi

    if [ $failed -eq 0 ]; then
        PASS=$((PASS + 1))
        echo -e "  ${GREEN}PASS${NC}  [$TOTAL] $description"
    else
        FAIL=$((FAIL + 1))
        echo -e "  ${RED}FAIL${NC}  [$TOTAL] $description"
        echo -e "        Reason: $reason"
        if [ -s "$stderr_file" ]; then
            echo -e "        stderr: $(head -1 "$stderr_file")"
        fi
    fi
}

# ---- Test Cases ----

echo "Running tests..."
echo ""

# 1. No arguments -> prints usage
run_test "No arguments prints usage" \
    "yes" "stdout_contains" "Usage:" \
    # no args

# 2. -h flag -> prints usage
run_test "-h prints usage" \
    "yes" "stdout_contains" "Usage:" \
    -h

# 3. --help flag -> prints usage
run_test "--help prints usage" \
    "yes" "stdout_contains" "Usage:" \
    --help

# 4. -v mp -> prints variable info
run_test "-v mp produces output" \
    "yes" "stdout_nonempty" "" \
    -v mp

# 5. -v mp pbl -> prints multiple variable info
run_test "-v mp pbl produces output" \
    "yes" "stdout_nonempty" "" \
    -v mp pbl

# 6. -v --all -> produces substantial output (>10 lines)
run_test "-v --all produces substantial output" \
    "yes" "stdout_lines_gt" "10" \
    -v --all

# 7. -v (missing argument) -> error message
run_test "-v with no variable name shows error" \
    "yes" "combined_contains" "ERROR: No variable name" \
    -v

# 8. --out-vars wind -> produces output
run_test "--out-vars wind produces output" \
    "yes" "stdout_nonempty" "" \
    --out-vars wind

# 9. --out-vars (missing argument) -> error message
run_test "--out-vars with no keyword shows error" \
    "yes" "combined_contains" "ERROR: No keyword" \
    --out-vars

# 10. --gen-nml generates a namelist file
GEN_NML_FILE="$TMPDIR/generated.nml"
run_test "--gen-nml creates namelist file" \
    "yes" "file_exists_nonempty" "$GEN_NML_FILE" \
    --gen-nml "$GEN_NML_FILE"

# 11. --gen-nml (missing argument) -> error message
run_test "--gen-nml with no filename shows error" \
    "yes" "combined_contains" "ERROR: No namelist" \
    --gen-nml

# 12. --check-nml with generated namelist -> produces output
run_test "--check-nml produces output" \
    "yes" "stdout_nonempty" "" \
    --check-nml "$GEN_NML_FILE"

# 13. --check-nml (missing argument) -> error message
run_test "--check-nml with no filename shows error" \
    "yes" "combined_contains" "ERROR: No namelist" \
    --check-nml

# 14. Nonexistent file -> error message
run_test "Nonexistent file shows error" \
    "any" "combined_contains" "Options file does not exist" \
    nonexistent_file_that_does_not_exist.nml

# ---- Summary ----

echo ""
echo "============================================"
echo "  Results: $PASS passed, $FAIL failed (of $TOTAL)"
echo "============================================"

if [ $FAIL -gt 0 ]; then
    exit 1
fi
exit 0
