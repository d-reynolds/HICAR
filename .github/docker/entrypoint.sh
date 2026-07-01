#!/usr/bin/env bash
# Register and run an EPHEMERAL GitHub Actions self-hosted runner.
#
# Ephemeral = the runner accepts exactly one job then exits. The container
# orchestrator (docker compose `restart: always`) brings it back and it
# re-registers fresh. This is the recommended pattern for a public repo: each
# job runs on a clean registration and a malicious job cannot persist on the box.
#
# Required environment:
#   REPO_URL    GitHub repo OR org to register against (scope is auto-detected):
#                 https://github.com/OWNER/REPO -> repository-level runner
#                 https://github.com/ORG        -> organization-level runner (serves
#                                                  every repo in the org)
#   one of:
#     RUNNER_TOKEN  a runner *registration* token (Settings > Actions > Runners > New)
#     GH_PAT        a fine-grained PAT used to mint a fresh registration token on every
#                   (re)start. Grant the LEAST privilege for your scope:
#                     org-level  -> org permission "Self-hosted runners: write" ONLY
#                     repo-level -> repo permission "Administration: write" (no narrower
#                                   repository scope exists for runner registration)
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
        # Auto-detect scope from REPO_URL's path:
        #   "ORG/REPO" (has a slash) -> repository-level registration endpoint
        #   "ORG"      (no slash)    -> organization-level registration endpoint
        url_path="${REPO_URL#https://github.com/}"
        url_path="${url_path%/}"          # tolerate a trailing slash
        if [[ "${url_path}" == */* ]]; then
            reg_api="https://api.github.com/repos/${url_path}/actions/runners/registration-token"
            echo "Fetching repo-level registration token for ${url_path}..."
        else
            reg_api="https://api.github.com/orgs/${url_path}/actions/runners/registration-token"
            echo "Fetching org-level registration token for ${url_path}..."
        fi
        RUNNER_TOKEN=$(curl -fsSL -X POST \
            -H "Authorization: Bearer ${GH_PAT}" \
            -H "Accept: application/vnd.github+json" \
            "${reg_api}" \
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

# SECURITY: scrub the high-value runner-admin PAT from the environment before
# the job runs, so a malicious job step cannot read it from `env`. The runner
# itself no longer needs it; the short-lived RUNNER_TOKEN is kept only so the
# cleanup trap can deregister (ephemeral runners also auto-deregister on exit).
unset GH_PAT

./run.sh
