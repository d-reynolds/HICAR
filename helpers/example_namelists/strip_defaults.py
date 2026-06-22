#!/usr/bin/env python3
"""Reduce a HICAR namelist to only the entries that differ from the defaults.

Given the freshly generated default namelist and an "example" namelist (the
default with some entries find/replaced to non-default values), emit a minimal
namelist that keeps ONLY the customized entries — so a user reading the example
sees just what that example changes, not all ~340 options.

What is preserved:
  * The full model header (everything the model writes above the first &group
    in `./bin/HICAR --gen-nml`) is reproduced verbatim, so the example keeps the
    model's own documentation of how to set namelist variables.
  * Each kept entry keeps its COMPLETE comment, including multi-line comments
    (the default writer wraps long descriptions / option lists across several
    lines); these continuation lines are carried over, not cut short.
  * All inline comments are re-aligned so the `!` marker starts at a fixed
    column (COMMENT_COL, currently 100); if the value itself runs past that
    column the comment follows one space later.

Variables whose name contains "folder" (e.g. the output / restart folder paths)
are ALWAYS kept, even at their default value, so each example clearly documents
where its output and restart files are written.

Empty sections are KEPT as a bare `&section /` header, never deleted: HICAR
reads its always-present sections (general, restart, domain, forcing, physics,
time_parameters, output) unconditionally, and a *missing* group makes the
namelist read hit end-of-file and `stop` the model. An empty *present* group
reads fine and leaves every value at its in-code default.

Usage:
    strip_defaults.py <default.nml> <example.nml> <out.nml> [--title TEXT]
"""
import argparse
import re

SECTION_RE = re.compile(r"^\s*&(\w+)")
END_RE = re.compile(r"^\s*/\s*$")
VAR_RE = re.compile(r"^\s*([A-Za-z]\w*)\s*=\s*(.*)$")

# 1-indexed column at which the `!` comment marker is placed. The value is
# padded out to this column; if the value is longer the comment follows one
# space after it (so the marker may land beyond COMMENT_COL).
COMMENT_COL = 53


def split_code_comment(line):
    """Split a namelist line into (code, comment) at the first unquoted `!`.

    `code` is the text before the comment (trailing whitespace removed);
    `comment` is the text after the `!` (surrounding whitespace removed), or
    None if the line carries no comment. Quote-aware so a `!` inside a quoted
    string is not mistaken for a comment.
    """
    quote = None
    for i, ch in enumerate(line):
        if quote:
            if ch == quote:
                quote = None
        elif ch in "'\"":
            quote = ch
        elif ch == "!":
            return line[:i].rstrip(), line[i + 1:].strip()
    return line.rstrip(), None


def parse(path):
    """Return an ordered list of (section, [entry]).

    Each entry is a dict with:
      name     - the variable name
      value    - the value text (comment stripped, surrounding ws removed),
                 used to compare against the defaults
      code     - the `    name = value` text (no comment, no trailing ws)
      comments - list of comment fragments for this entry, in order; the first
                 came from the variable's own line, the rest from the wrapped
                 continuation (`<spaces> ! ...`) lines beneath it
    """
    sections = []
    entries = None   # current section's entry list, or None when outside a section
    entry = None     # current entry dict, or None
    with open(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            m = SECTION_RE.match(line)
            if m:
                entries = []
                sections.append((m.group(1), entries))
                entry = None
                continue
            if entries is None:
                continue
            if END_RE.match(line):
                entries = None
                entry = None
                continue
            if line.strip() == "":
                # blank line separates entries -> end the current one
                entry = None
                continue
            code, comment = split_code_comment(line)
            vm = VAR_RE.match(code)
            if vm:
                entry = {
                    "name": vm.group(1),
                    "value": vm.group(2).strip(),
                    "code": code,
                    "comments": [comment] if comment else [],
                }
                entries.append(entry)
            elif comment and entry is not None:
                # a wrapped continuation comment line for the current variable
                entry["comments"].append(comment)
            # else: stray line (no code, no comment) -> ignored
    return sections


def read_header(path):
    """Return the verbatim lines above the first &group (the model header)."""
    header = []
    with open(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if SECTION_RE.match(line):
                break
            header.append(line.rstrip())
    while header and header[0].strip() == "":
        header.pop(0)
    while header and header[-1].strip() == "":
        header.pop()
    return header


def format_entry(code, comments):
    """Render an entry as one or more lines with the comment column aligned."""
    if not comments:
        return [code]
    lines = []
    pad = max(1, COMMENT_COL - 1 - len(code))
    lines.append(code + " " * pad + "! " + comments[0])
    cont = " " * (COMMENT_COL - 1)
    for frag in comments[1:]:
        lines.append(cont + "! " + frag)
    return lines


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("default")
    p.add_argument("example")
    p.add_argument("out")
    p.add_argument("--title", default="HICAR example namelist (auto-generated)")
    args = p.parse_args()

    # default values: {section: {var: value}}
    defaults = {}
    for sec, entries in parse(args.default):
        sec_defaults = defaults.setdefault(sec, {})
        for entry in entries:
            sec_defaults[entry["name"]] = entry["value"]

    example = parse(args.example)

    n_kept = 0
    out_lines = [
        f"! {args.title}",
        "! ---------------------------------------------------------------------------",
        "! Auto-generated by helpers/example_namelists/generate_examples.sh — DO NOT EDIT.",
        "! Only entries that differ from the model defaults are shown (plus the",
        "! output/restart folder paths, always listed); every other option keeps its",
        "! default value. Empty sections are retained on purpose (removing an",
        "! always-read section makes HICAR stop at namelist read).",
        "! Regenerate the full annotated default with:  ./bin/HICAR --gen-nml default.nml",
        "! ---------------------------------------------------------------------------",
        "",
    ]
    # Preserve the model's own header (the documentation block --gen-nml writes
    # above the first group), so the example keeps that guidance verbatim.
    out_lines.extend(read_header(args.default))
    out_lines.append("")

    for sec, entries in example:
        out_lines.append(f" &{sec}")
        sec_defaults = defaults.get(sec, {})
        out_lines.append("")
        for entry in entries:
            # Keep an entry if it is customized (differs from / absent in the
            # defaults) OR if its name contains "folder" — output/restart folder
            # paths are always shown so the example documents where files land,
            # even when left at the default.
            if sec_defaults.get(entry["name"]) != entry["value"] or "folder" in entry["name"].lower():
                out_lines.extend(format_entry(entry["code"], entry["comments"]))
                out_lines.append("")
                n_kept += 1
        out_lines.append(" /")
        out_lines.append("")

    with open(args.out, "w") as fh:
        fh.write("\n".join(out_lines) + "\n")

    print(f"  {args.out}: kept {n_kept} non-default entries across "
          f"{len(example)} sections")


if __name__ == "__main__":
    main()
