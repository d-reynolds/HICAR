#!/bin/bash
# ===========================================================================
# Self-resubmitting SLURM job chain for HICAR runs longer than the wall-time.
#
# Submit ONCE from your run/input directory (the same directory you would run
# HICAR from manually, so that relative paths in the namelist resolve):
#
#     sbatch batch_submit_SLURM.sh
#
# How it works (modern rewrite of the ICAR-era chain system):
#   1. Before running, the job queues a successor of itself with
#      --dependency=afternotok:<this job>. If this job dies (typically the
#      scheduler killing it at the wall-time limit), the successor starts;
#      if HICAR finishes its simulation period, the successor is cancelled.
#   2. Each job inspects the restart folder(s) named in the namelist. If
#      restart files exist, it edits the namelist IN PLACE to set
#          restart_run  = .True.
#          restart_date = <datetime of the newest restart file>
#      (the date is parsed from HICAR's restart filename convention,
#      <prefix>_YYYY-MM-DD_HH-MM-SS.nc). No template files, no Python.
#   3. Safety rails: the chain aborts if a successor finds the SAME restart
#      date as the previous attempt (no forward progress — e.g. the model
#      crashes deterministically, or restartinterval writes less often than
#      one wall-time window), and after MAX_CHAIN jobs total.
#
# Requirements:
#   * The namelist must contain `restart_run` and `restart_date` entries in
#     its &restart block (any namelist from `HICAR --gen-nml` does).
#   * restartinterval must be frequent enough that at least one restart file
#     is written per wall-time window, or the chain will (correctly) stop.
#
# Test the restart-detection logic without SLURM:
#     ./batch_submit_SLURM.sh --setup-only my_options.nml
# ===========================================================================
#SBATCH --job-name="HICAR_chain"
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --output=hicar_chain-%j.out
#SBATCH --error=hicar_chain-%j.err
### Uncomment / edit for your system:
##SBATCH --account=YOUR_ACCOUNT
##SBATCH --partition=YOUR_PARTITION

# --------------------------- user configuration ---------------------------
NML="${HICAR_NML:-options.nml}"          # namelist, relative to the submit dir
EXE="${HICAR_EXE:-$HOME/HICAR/bin/HICAR}"
MAX_CHAIN="${MAX_CHAIN:-100}"            # absolute cap on chained jobs
LAUNCH="srun"                            # launcher; ranks come from SBATCH
# ---------------------------------------------------------------------------

set -u

setup_only=false
if [ "${1:-}" = "--setup-only" ]; then
    setup_only=true
    [ -n "${2:-}" ] && NML="$2"
fi

# Run from the submission directory so namelist-relative paths resolve.
cd "${SLURM_SUBMIT_DIR:-$(pwd)}" || exit 1

[ -f "$NML" ] || { echo "ERROR: namelist '$NML' not found in $(pwd)" >&2; exit 1; }
state_file=".hicar_chain_state_$(basename "$NML")"

# --- runaway-chain guard ----------------------------------------------------
depth=${CHAIN_DEPTH:-1}
if [ "$depth" -gt "$MAX_CHAIN" ]; then
    echo "ERROR: chain depth $depth exceeds MAX_CHAIN=$MAX_CHAIN — stopping." >&2
    exit 1
fi
echo "HICAR chain: job $depth (max $MAX_CHAIN), namelist $NML"

# --- find the newest restart file in the folder(s) named in the namelist ----
# restart_folder may be a comma-separated list for nested runs; all nests
# write restarts at the same model time, so the newest datetime is shared.
folders=$(sed -n 's/^[[:space:]]*restart_folder[[:space:]]*=[[:space:]]*//p' "$NML" \
          | head -1 | sed 's/!.*//' | tr ',' '\n' | sed "s/[\"' ]//g" | sed '/^$/d')
[ -n "$folders" ] || folders="../restart/"

newest_stamp=""
for d in $folders; do
    for f in "$d"/*.nc; do
        [ -e "$f" ] || continue
        stamp=$(basename "$f" | grep -o '[0-9]\{4\}-[0-9]\{2\}-[0-9]\{2\}_[0-9]\{2\}-[0-9]\{2\}-[0-9]\{2\}' | head -1)
        [ -n "$stamp" ] || continue
        if [ -z "$newest_stamp" ] || [[ "$stamp" > "$newest_stamp" ]]; then
            newest_stamp="$stamp"
        fi
    done
done

if [ -n "$newest_stamp" ]; then
    # <prefix>_YYYY-MM-DD_HH-MM-SS.nc  ->  'YYYY-MM-DD HH:MM:SS'
    d_part=${newest_stamp%%_*}
    t_part=$(echo "${newest_stamp#*_}" | tr '-' ':')
    restart_date="${d_part} ${t_part}"

    # --- no-progress guard: same restart point as the previous attempt? -----
    if [ -f "$state_file" ] && [ "$(cat "$state_file")" = "$restart_date" ]; then
        echo "ERROR: newest restart ($restart_date) is unchanged since the previous" >&2
        echo "chain job — no forward progress (crash loop, or restartinterval too" >&2
        echo "long for the wall-time). Stopping the chain." >&2
        exit 1
    fi
    echo "$restart_date" > "$state_file"

    echo "Restarting from: $restart_date"
    sed -i.chain_bak \
        -e "s|^\([[:space:]]*\)restart_run[[:space:]]*=.*|\1restart_run = .True.|" \
        -e "s|^\([[:space:]]*\)restart_date[[:space:]]*=.*|\1restart_date = '${restart_date}'|" \
        "$NML"
    grep -q "restart_run = .True." "$NML" || {
        echo "ERROR: failed to set restart options in $NML (no restart_run entry?)" >&2
        exit 1
    }
else
    echo "No restart files found — running $NML from its start date."
    rm -f "$state_file"
fi

if [ "$setup_only" = true ]; then
    echo "--setup-only: namelist prepared; not submitting or running."
    exit 0
fi

# --- queue the successor BEFORE running -------------------------------------
# afternotok: the successor only starts if this job fails or is killed
# (wall-time kill => TIMEOUT counts as not-ok). On success we cancel it.
next_id=$(sbatch --parsable --kill-on-invalid-dep=yes \
                 --dependency=afternotok:"$SLURM_JOB_ID" \
                 --export=ALL,CHAIN_DEPTH=$((depth + 1)),HICAR_NML="$NML",HICAR_EXE="$EXE",MAX_CHAIN="$MAX_CHAIN" \
                 "$0")
echo "Queued successor job $next_id (starts only if this job does not finish)."

# --- run the model -----------------------------------------------------------
if $LAUNCH "$EXE" "$NML"; then
    echo "HICAR completed its simulation period — ending the chain."
    scancel "$next_id" 2>/dev/null
    rm -f "$state_file"
else
    rc=$?
    echo "HICAR exited with status $rc — successor $next_id will take over from the latest restart."
    exit "$rc"
fi
