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
# compared over ALL variables against tests/snowpack/tolerances_snowpack.yaml.
#
# Usage:
#   test_snowpack_compare.sh <hicar_repo> <cpp_exe> <fortran_exe> [np]
#   test_snowpack_compare.sh <hicar_repo> --bless [--snowpack-dir DIR] [--reason TEXT] [--force]
#
#   <hicar_repo>  repo root
#   <cpp_exe>     HICAR executable built with the C++ SNOWPACK   (reference)
#   <fortran_exe> HICAR executable built with the Fortran SNOWPACK (candidate)
#   [np]          MPI ranks for each run (default 10; both runs use the same)
#
#   --bless       instead of running the comparison, post the
#                 `hicar-snowpack-parity-blessed=success` commit status to HEAD.
#                 The status DESCRIPTION records the upstream SNOWPACK commit
#                 (snowpack=<sha>) resolved from the fetched checkout, which is
#                 the diff anchor for tests/snowpack/snowpack_divergence_report.sh: when
#                 a later comparison run diverges, `git diff <blessed-sha>..HEAD`
#                 in the SNOWPACK repo lists the upstream C++ changes that still
#                 need porting to the Fortran driver. Bless only commits whose
#                 comparison run actually PASSED.
#   --snowpack-dir DIR  fetched SNOWPACK checkout to read the SHA from
#                       (default: <hicar_repo>/build/external/SNOWPACK)
#   --reason TEXT       short audit note appended to the status description
#
# Exit 0 = within tolerance; 1 = a run failed or outputs diverged.
# ===========================================================================
set -u

BLUE='\033[0;36m'; GREEN='\033[0;32m'; RED='\033[0;31m'; NC='\033[0m'

if [ $# -lt 2 ]; then
    echo "Usage: $0 <hicar_repo> <cpp_exe> <fortran_exe> [np]"
    echo "       $0 <hicar_repo> --bless [--snowpack-dir DIR] [--reason TEXT]"
    exit 2
fi

hicar_repo="$(cd "$1" && pwd)"

# set_var <nml_file> <name> <value> <group>: set a namelist variable by name,
# inserting it into &<group> if absent (replaces the old `sed -i` edits). See
# helpers/example_namelists/set_nml_var.py.
set_var() { "${PYTHON:-python3}" "$hicar_repo/helpers/example_namelists/set_nml_var.py" "$1" "$2" "$3" --group "$4" --insert || exit 1; }

# --- bless mode: post the parity-blessed status (with the SNOWPACK SHA) -----
if [ "$2" = "--bless" ]; then
    shift 2
    snowpack_dir="$hicar_repo/build/external/SNOWPACK"
    reason=""
    force=false
    while [ $# -gt 0 ]; do
        case "$1" in
            --snowpack-dir) snowpack_dir="$2"; shift 2;;
            --reason)       reason="$2"; shift 2;;
            --force)        force=true; shift;;
            *) echo "Unknown bless option: $1"; exit 2;;
        esac
    done
    report_dir="${hicar_repo}/tests/figures/snowpack_compare"
    meta="${report_dir}/parity_meta.txt"

    # SNOWPACK anchor SHA: from the fetched checkout when present (local bless
    # after a build), else from the comparison run's own provenance stamp (the
    # CI bless job runs without a build dir — the evidence artifact carries the
    # SHA the comparison actually used).
    sp_sha=$(git -C "$snowpack_dir" rev-parse HEAD 2>/dev/null || true)
    if [ -z "$sp_sha" ]; then
        sp_sha=$(sed -n 's/^snowpack_commit=//p' "$meta" 2>/dev/null)
    fi
    if [ -z "$sp_sha" ] || [ "$sp_sha" = "unknown" ]; then
        echo -e "${RED}Cannot determine the SNOWPACK anchor SHA: no git checkout at${NC}"
        echo -e "${RED}${snowpack_dir} and no snowpack_commit in ${meta}.${NC}"
        echo "Build SNOWPACK or run the comparison first, or pass --snowpack-dir."
        exit 1
    fi
    head_sha=$(git -C "$hicar_repo" rev-parse HEAD)

    # --- authorization gate: admins/maintainers only -------------------------
    # Blessing moves the parity anchor for everyone, so it is restricted to
    # repo admin/maintain roles (the CI path enforces the same on the workflow
    # actor). NOT bypassable by --force: --force only overrides the evidence
    # checks below, never authorization.
    slug="${GITHUB_REPOSITORY:-$(cd "$hicar_repo" && gh repo view --json nameWithOwner -q .nameWithOwner 2>/dev/null || true)}"
    [ -n "$slug" ] || { echo -e "${RED}Cannot determine owner/repo for gh${NC}"; exit 1; }
    actor=$(gh api user --jq .login 2>/dev/null || true)
    role=$(gh api "/repos/$slug/collaborators/${actor}/permission" \
             --jq '.role_name // .permission' 2>/dev/null || echo unknown)
    case "$role" in
        admin|maintain) : ;;
        *)
            echo -e "${RED}Parity bless is restricted to repo admins/maintainers.${NC}"
            echo "  gh user: ${actor:-unknown}, role on ${slug}: ${role}"
            exit 1
            ;;
    esac

    # --- pass gate: only bless what the comparison actually certified -------
    # Every comparison run stamps tests/figures/snowpack_compare/parity_meta.txt
    # with its PASS/FAIL verdict, the HICAR commit it ran on, and the SNOWPACK
    # SHA it compared. Refuse to bless unless that evidence (a) exists,
    # (b) says PASS, and (c) was produced from THIS commit — so a bless always
    # certifies a real, current, passing comparison. --force overrides (e.g.
    # evidence regenerated on a dirty tree), and is logged in the status text.
    if [ "$force" != true ]; then
        [ -f "$meta" ] || {
            echo -e "${RED}No comparison evidence (${meta}).${NC}"
            echo "Run the comparison first; bless certifies a passing run. (--force to override)"
            exit 1
        }
        meta_status=$(sed -n 's/^compare_status=//p' "$meta")
        meta_hicar=$(sed -n 's/^hicar_commit=//p' "$meta")
        meta_sp=$(sed -n 's/^snowpack_commit=//p' "$meta")
        [ "$meta_status" = "PASS" ] || {
            echo -e "${RED}Refusing to bless: last comparison was ${meta_status:-unknown}, not PASS.${NC} (--force to override)"
            exit 1
        }
        [ "$meta_hicar" = "$head_sha" ] || {
            echo -e "${RED}Refusing to bless: evidence is from commit ${meta_hicar:0:12}, HEAD is ${head_sha:0:12}.${NC}"
            echo "Re-run the comparison on this commit. (--force to override)"
            exit 1
        }
        [ "$meta_sp" = "$sp_sha" ] || {
            echo -e "${RED}Refusing to bless: evidence compared SNOWPACK ${meta_sp:0:12}, checkout now at ${sp_sha:0:12}.${NC}"
            echo "The fetched SNOWPACK changed since the run. (--force to override)"
            exit 1
        }
    fi

    desc="snowpack=${sp_sha}${reason:+; $reason}"
    [ "$force" = true ] && desc="${desc}; FORCED"
    desc="${desc:0:140}"
    gh api -X POST "/repos/$slug/statuses/$head_sha" \
        -f state=success -f context=hicar-snowpack-parity-blessed \
        -f description="$desc" >/dev/null
    echo -e "${GREEN}Blessed ${head_sha}${NC} (hicar-snowpack-parity-blessed=success)"
    echo -e "  SNOWPACK parity anchor: ${BLUE}${sp_sha}${NC}"

    # --- archive the comparison evidence as a permanent GitHub release ------
    # Workflow artifacts expire (<=90 d); a release does not. Each bless gets a
    # `snowpack-parity/<date>-<sha>` release carrying the per-variable stats
    # JSON, the diff maps, and the provenance stamp, so the parity level can be
    # tracked across blesses (compare reports with tests/snowpack/parity_trend.py).
    if [ -f "${report_dir}/parity_report.json" ]; then
        tag="snowpack-parity/$(date -u +%Y%m%d)-${head_sha:0:7}"
        tarball="$(mktemp -d)/parity_diffmaps.tar.gz"
        tar -czf "$tarball" -C "$report_dir" \
            $(cd "$report_dir" && ls -d diffmaps *.png 2>/dev/null) 2>/dev/null || tarball=""
        notes="SNOWPACK C++ vs Fortran parity bless.

* HICAR commit: \`${head_sha}\`
* SNOWPACK (fortran-bindings) anchor: \`${sp_sha}\`
* $(grep -E 'date_utc|ranks|compare_status' "${report_dir}/parity_meta.txt" 2>/dev/null | tr '\n' ' ')
${reason:+* Reason: $reason}

Assets: per-variable comparison stats (parity_report.json), spatial difference
maps (parity_diffmaps.tar.gz). Compare against an earlier bless with
tests/snowpack/parity_trend.py <old.json> <new.json>."
        if gh release create "$tag" --repo "$slug" --target "$head_sha" \
              --title "SNOWPACK parity $(date -u +%Y-%m-%d) (${head_sha:0:7})" \
              --notes "$notes" \
              "${report_dir}/parity_report.json" \
              "${report_dir}/parity_meta.txt" \
              ${tarball:+"$tarball"} 2>/dev/null; then
            echo -e "  Evidence archived: ${BLUE}release ${tag}${NC}"
        else
            echo -e "${RED}  Release archiving failed (no contents:write permission?) — bless status still posted.${NC}"
        fi
    else
        echo -e "${RED}  No parity_report.json found in ${report_dir} — nothing archived.${NC}"
        echo "  (Run the comparison first so the bless can archive its evidence.)"
    fi
    exit 0
fi

if [ $# -lt 3 ]; then
    echo "Usage: $0 <hicar_repo> <cpp_exe> <fortran_exe> [np]"
    exit 2
fi
cpp_exe="$2"
fortran_exe="$3"
# Ranks per run. Default to min(10, available cores) so the run is fast on a
# beefy machine but does NOT oversubscribe a small GitHub CI runner (2-4 cores).
if [ -n "${4:-}" ]; then
    np="$4"
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
tol_spec="${hicar_repo}/tests/snowpack/tolerances_snowpack.yaml"
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
# Record the upstream SNOWPACK commit in the log: this is the SHA that becomes
# the parity anchor if this run is blessed (see --bless), and the audit trail
# for the divergence routine in tests/snowpack/snowpack_divergence_report.sh.
sp_sha=$(git -C "$hicar_repo/build/external/SNOWPACK" rev-parse HEAD 2>/dev/null || echo "unknown")
echo "  SNOWPACK commit:     $sp_sha"

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
# (passing ones included) — at bless time it is archived to a GitHub release so
# the parity level can be consulted over time (see --bless).
report_dir="${hicar_repo}/tests/figures/snowpack_compare"
mkdir -p "$report_dir"
"$python_exe" "$compare_script" "$ref_out" "$cand_out" \
    --mode tolerance --tolerance-spec "$tol_spec" \
    --figures-dir "$report_dir" \
    --report-json "${report_dir}/parity_report.json"
cmp_status=$?

# Spatial difference maps for the key snow/surface fields, generated on PASS as
# well as FAIL (compare_outputs only plots failures): these are the figures the
# bless archives, and the visual baseline for "how did the residual change".
"$python_exe" "${hicar_repo}/tests/snowpack/plot_snowpack_diff.py" \
    "$ref_out" "$cand_out" "${report_dir}/diffmaps" >/dev/null 2>&1 \
    && echo -e "Difference maps: ${BLUE}${report_dir}/diffmaps${NC}" \
    || echo -e "${RED}(diff-map generation failed — non-fatal)${NC}"

# Provenance stamp for the archive.
{
    echo "date_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "hicar_commit=$(git -C "$hicar_repo" rev-parse HEAD 2>/dev/null || echo unknown)"
    echo "snowpack_commit=${sp_sha}"
    echo "ranks=${np}"
    echo "compare_status=$([ $cmp_status -eq 0 ] && echo PASS || echo FAIL)"
} > "${report_dir}/parity_meta.txt"

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
