#!/bin/bash
# Print the most recent commit (an ancestor of HEAD, excluding HEAD itself) that
# carries a `hicar-regression-blessed=success` commit status — i.e. the regression
# reference. This replaces the old tests/regression_commit.txt file.
#
# Requires `gh` and a token with at least read access (GH_TOKEN in CI, or a
# `gh auth login` session locally). Repo slug comes from GITHUB_REPOSITORY (CI) or
# the current repo's remote.
#
# Usage: resolve_blessed_commit.sh [hicar_repo]   # prints the SHA on stdout
set -euo pipefail

repo=$(cd "${1:-.}" && pwd)
slug="${GITHUB_REPOSITORY:-$(gh repo view --json nameWithOwner -q .nameWithOwner 2>/dev/null || true)}"
[ -n "$slug" ] || { echo "resolve_blessed_commit: cannot determine owner/repo" >&2; exit 1; }

head=$(git -C "$repo" rev-parse HEAD)
# Walk recent history newest-first; first blessed ancestor wins.
while read -r sha; do
    [ "$sha" = "$head" ] && continue
    state=$(gh api "/repos/$slug/commits/$sha/status" \
              --jq '[.statuses[] | select(.context=="hicar-regression-blessed")][0].state // "none"' \
              2>/dev/null || echo none)
    if [ "$state" = "success" ]; then
        echo "$sha"
        exit 0
    fi
done < <(git -C "$repo" log -n 200 --format=%H)

echo "resolve_blessed_commit: no hicar-regression-blessed commit found in the last 200 commits" >&2
exit 1
