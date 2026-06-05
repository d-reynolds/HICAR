#!/usr/bin/env bash
# Register and run an EPHEMERAL GitHub Actions self-hosted runner.
#
# Ephemeral = the runner accepts exactly one job then exits. The container
# orchestrator (docker compose `restart: always`) brings it back and it
# re-registers fresh. This is the recommended pattern for a public repo: each
# job runs on a clean registration and a malicious job cannot persist on the box.
#
# Required environment:
#   REPO_URL    e.g. https://github.com/d-reynolds/HICAR
#   one of:
#     RUNNER_TOKEN  a runner *registration* token (Settings > Actions > Runners > New;
#                   valid ~1h — fine for a long-lived container that registers once)
#     GH_PAT        a PAT with repo admin scope; the script fetches a fresh
#                   registration token via the API on every (re)start
#
# Optional:
#   RUNNER_NAME   default: hicar-gpu-<short-hostname>
#   RUNNER_LABELS default: self-hosted,gpu,nvhpc
set -euo pipefail

: "${REPO_URL:?set REPO_URL to the GitHub repo, e.g. https://github.com/owner/repo}"
RUNNER_NAME="${RUNNER_NAME:-hicar-gpu-$(hostname | cut -c1-8)}"
RUNNER_LABELS="${RUNNER_LABELS:-self-hosted,gpu,nvhpc}"

# Resolve a registration token.
if [ -z "${RUNNER_TOKEN:-}" ]; then
    if [ -n "${GH_PAT:-}" ]; then
        owner_repo="${REPO_URL#https://github.com/}"
        echo "Fetching registration token via API for ${owner_repo}..."
        RUNNER_TOKEN=$(curl -fsSL -X POST \
            -H "Authorization: Bearer ${GH_PAT}" \
            -H "Accept: application/vnd.github+json" \
            "https://api.github.com/repos/${owner_repo}/actions/runners/registration-token" \
            | jq -r .token)
    else
        echo "ERROR: provide RUNNER_TOKEN or GH_PAT" >&2
        exit 1
    fi
fi

cleanup() {
    echo "Removing runner registration..."
    ./config.sh remove --token "${RUNNER_TOKEN}" || true
}
trap cleanup EXIT INT TERM

./config.sh \
    --unattended \
    --url "${REPO_URL}" \
    --token "${RUNNER_TOKEN}" \
    --name "${RUNNER_NAME}" \
    --labels "${RUNNER_LABELS}" \
    --ephemeral \
    --replace

# A quick sanity check that the GPU is visible inside the container.
if command -v nvidia-smi >/dev/null 2>&1; then
    nvidia-smi -L || echo "WARNING: nvidia-smi found no GPUs — was the container started with --gpus all?"
else
    echo "WARNING: nvidia-smi not available in container."
fi

./run.sh
