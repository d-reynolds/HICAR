#!/bin/bash
# ===========================================================================
# Auto-generate the example namelists in helpers/example_namelists/.
#
# For every per-example find/replace script in gen_example_nmls/<name>.sh this:
#   1. generates a fresh default namelist from the HICAR executable
#      (./bin/HICAR --gen-nml), so the defaults always track the current build;
#   2. applies <name>.sh (sed find/replace) to set the entries that example
#      customizes;
#   3. reduces the result to ONLY the non-default entries (strip_defaults.py),
#      writing helpers/example_namelists/<name>.nml.
#
# Users never hand-write the <name>.nml files — they add a gen_example_nmls/
# script and the .nml is regenerated here (and in CI on push to main/develop).
#
# Usage:
#   helpers/example_namelists/generate_examples.sh [hicar_exe]
#     hicar_exe   path to the HICAR executable (default: <repo>/bin/HICAR)
# ===========================================================================
set -u

here="$(cd "$(dirname "$0")" && pwd)"
repo="$(cd "$here/../.." && pwd)"
gen_dir="$here/gen_example_nmls"
strip="$here/strip_defaults.py"

hicar_exe="${1:-$repo/bin/HICAR}"
if [ ! -x "$hicar_exe" ]; then
    echo "ERROR: HICAR executable not found/executable: $hicar_exe" >&2
    echo "Build HICAR first, or pass the exe path as the first argument." >&2
    exit 1
fi

python_exe="$(command -v python3 || command -v python || true)"
[ -n "$python_exe" ] || { echo "ERROR: python3 not found" >&2; exit 1; }

shopt -s nullglob
scripts=("$gen_dir"/*.sh)
if [ ${#scripts[@]} -eq 0 ]; then
    echo "No example scripts in $gen_dir — nothing to generate."
    exit 0
fi

# Work in a temp dir. NOTE: `--gen-nml` refuses to overwrite an existing file,
# so the target path must not exist yet (do not point it at an mktemp file).
tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT
default_nml="$tmpdir/default.nml"
"$hicar_exe" --gen-nml "$default_nml" >/dev/null 2>&1
if [ ! -s "$default_nml" ]; then
    echo "ERROR: '$hicar_exe --gen-nml' did not produce a namelist" >&2; exit 1
fi

status=0
for script in "${scripts[@]}"; do
    name="$(basename "$script" .sh)"
    full="$tmpdir/$name.full.nml"
    cp "$default_nml" "$full"

    # Apply the example's find/replace. The scripts call `sed -i'.bak' ... $1`,
    # mirroring tests/Test_Cases/input/nml_gen_scripts; clean up the .bak after.
    if ! bash "$script" "$full"; then
        echo "ERROR: example script failed: $script" >&2
        status=1
    fi
    rm -f "$full.bak"

    "$python_exe" "$strip" "$default_nml" "$full" "$here/$name.nml" \
        --title "HICAR example namelist: $name" || status=1
done

exit $status
