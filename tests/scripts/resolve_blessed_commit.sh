#!/usr/bin/env bash
# Resolve the most recent ancestor commit carrying a given blessed commit status.
#
# HICAR's blessing mechanism stores no files: a commit is "blessed" by posting a
# GitHub commit status (state=success) in a given context:
#   hicar-regression-blessed       — regression reference (tests/scripts/test_regression.sh)
#   snow-parity                    — last commit where the native-Fortran SNOWPACK
#                                    port matched the C++ build (posted by
#                                    snowpack-compare on PASS). Its status DESCRIPTION
#                                    records the upstream SNOWPACK SHA (snowpack=<sha>),
#                                    the divergence anchor for snowpack_divergence_report.sh.
#
# This script walks the first-parent history from HEAD and prints the SHA of the
# most recent commit whose combined status contains <context>=success.
#
# Usage:
#   resolve_blessed_commit.sh <hicar_repo> [context] [--with-description] [--exclude-head]
#
#   context             status context to search for
#                       (default: hicar-regression-blessed)
#   --with-description  print "<sha>\t<description>" instead of just the SHA
#   --exclude-head      skip HEAD itself (regression semantics: comparing a
#                       commit against itself would trivially pass, so the
#                       reference must be a strict ancestor; the parity anchor
#                       by contrast MAY be HEAD)
#
# Requires: gh (authenticated), git. Exit 1 if none found in the last 200 commits.
set -uo pipefail

repo="${1:?usage: resolve_blessed_commit.sh <hicar_repo> [context] [--with-description] [--exclude-head]}"
context="hicar-regression-blessed"
with_desc=false
exclude_head=false
shift
while [ $# -gt 0 ]; do
    case "$1" in
        --with-description) with_desc=true; shift;;
        --exclude-head)     exclude_head=true; shift;;
        *) context="$1"; shift;;
    esac
done
head_sha=$(git -C "$repo" rev-parse HEAD)

slug="${GITHUB_REPOSITORY:-$(cd "$repo" && gh repo view --json nameWithOwner -q .nameWithOwner 2>/dev/null || true)}"
[ -n "$slug" ] || { echo "resolve_blessed_commit: cannot determine owner/repo for gh" >&2; exit 1; }

# Walk first-parent history (most recent first). One API call per commit, but
# the blessed commit is normally close to HEAD so this terminates quickly.
# The jq emits "FOUND\t<description>" only when the status exists, so a single
# call both tests for the status and retrieves its description.
for sha in $(git -C "$repo" rev-list --first-parent -n 200 HEAD); do
    [ "$exclude_head" = true ] && [ "$sha" = "$head_sha" ] && continue
    line=$(gh api "/repos/$slug/commits/$sha/status" \
        --jq "[.statuses[] | select(.context==\"$context\" and .state==\"success\")][0] | \
              if . == null then empty else \"FOUND\t\(.description // \"\")\" end" 2>/dev/null) || continue
    if [ -n "$line" ]; then
        if [ "$with_desc" = true ]; then
            printf '%s\t%s\n' "$sha" "${line#FOUND	}"
        else
            echo "$sha"
        fi
        exit 0
    fi
done

echo "resolve_blessed_commit: no commit with $context=success found in the last 200 first-parent commits" >&2
exit 1
