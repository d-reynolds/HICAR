#!/bin/bash
# Build a HICAR forcing-file list: one absolute, double-quoted path per line.
#
# Usage:
#   ./filelist_script.sh "<glob pattern>" <output_file>
#   ./filelist_script.sh "../forcing/icar_out_201*" file_list.txt
#
# Quote the glob pattern so it is expanded here, not in your shell's cwd.
# The output file is TRUNCATED first: re-running always yields a list that
# matches the pattern exactly. (Historical behavior was to append, which
# silently produced stale/duplicated lists on re-runs.)

pattern=$1
out_file=$2

if [ -z "$pattern" ] || [ -z "$out_file" ]; then
    echo "Usage: $0 \"<glob pattern>\" <output_file>" >&2
    exit 2
fi

: > "$out_file"
readlink -f $pattern | while read -r line; do
    echo "\"$line\"" >> "$out_file"
done

n=$(wc -l < "$out_file" | tr -d ' ')
if [ "$n" -eq 0 ]; then
    echo "filelist_script: WARNING: no files matched pattern: $pattern" >&2
    exit 1
fi
echo "Wrote $n file path(s) to $out_file"
