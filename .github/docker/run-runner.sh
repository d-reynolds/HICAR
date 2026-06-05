#!/usr/bin/env bash
# =============================================================================
# Run the HICAR GPU CI runner with a FRESH container PER JOB.
# =============================================================================
# The runner is --ephemeral (one job, then exits). Because each iteration here
# uses `docker run --rm`, the container is destroyed after every job and the
# next job starts from the pristine image — a job's filesystem writes can never
# persist to a later job. This is defense-in-depth: gpu.yml is configured to run
# ONLY vetted code (no pull_request trigger), but fresh-per-job also prevents
# state leaking between trusted builds.
#
# Contrast with `docker compose up -d` (restart: always), which re-execs the
# runner in the SAME container, so its writable layer persists between jobs.
# Use compose to BUILD the image; use this script to RUN it in production.
#
# Usage:  ./run-runner.sh        (Ctrl-C to stop)
#         IMAGE=hicar-gpu-runner:latest ./run-runner.sh
# Keep it alive across reboots via the systemd unit (hicar-gpu-runner.service).
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"

IMAGE="${IMAGE:-hicar-gpu-runner:latest}"
[ -f runner.env ] || { echo "ERROR: runner.env not found (copy from runner.env.example)"; exit 1; }

echo "Starting HICAR GPU runner loop with image '${IMAGE}' (fresh container per job)."
echo "Stop with Ctrl-C or 'systemctl stop hicar-gpu-runner'."
while true; do
    # GPU access via the NVIDIA Container Runtime (works in CDI mode, where the
    # legacy --gpus flag is disabled). Requires the nvidia runtime to be
    # registered: nvidia-ctk runtime configure --runtime=docker && systemctl restart docker
    docker run --rm \
        --runtime=nvidia \
        -e NVIDIA_VISIBLE_DEVICES=all \
        -e NVIDIA_DRIVER_CAPABILITIES=compute,utility \
        --env-file runner.env \
        --shm-size=2gb \
        --ulimit memlock=-1 --ulimit stack=67108864 \
        --name hicar-gpu-runner \
        "${IMAGE}" || true
    echo "$(date -u +%FT%TZ) runner container exited; starting a fresh one in 5s..."
    sleep 5
done
