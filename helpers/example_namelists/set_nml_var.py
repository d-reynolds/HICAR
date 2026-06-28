#!/usr/bin/env python3
"""Set one namelist variable to a value, in place — independent of its default.

The example generator scripts (gen_example_nmls/*.sh) and the integration-test
generators (tests/Test_Cases/input/nml_gen_scripts/*.sh) call this to assign a
chosen value to a variable in a namelist.

Unlike a `sed 's/name = <default>/name = <value>/'`, this matches the variable
by NAME only: it finds the `name = ...` assignment line (anchored on the whole
identifier, so it never substring-matches another variable — e.g. `start_date`
vs `restart_date`, or `output` vs `output_vars`) and rewrites whatever value is
currently there. Two consequences:

  * The scripts never encode the model's current defaults, so they keep working
    when a default value changes.
  * If the variable is missing, behaviour depends on --group (see below).

The trailing comment on the line is preserved (quote-aware), so the strip step
can still carry the variable's documentation into the example. The exact
comment column does not matter here; strip_defaults.py re-aligns it.

--group NAME documents which namelist group the variable lives in. On its own it
is advisory (used only if an insert is needed). It becomes load-bearing with
--insert.

Missing variables:
  * By default the call fails loudly (exit 1), so a stale script that names a
    renamed/removed variable breaks the build instead of silently leaving the
    default in place. The example generator scripts rely on this.
  * With --insert (which requires --group): if the variable is not already
    present, a fresh `    name = value` line is inserted into the &GROUP
    section. This is for the integration-test scripts, which start from a
    *stripped* example namelist (where variables left at their default —
    outputinterval, restart_run, … — are absent) and still need to set them.
    If GROUP itself is not found, the call fails (exit 1).

Usage:  set_nml_var.py <file.nml> <name> <value> [--group GROUP] [--insert]
"""
import re
import sys

# leading indent (group 1), the whole variable name (group 2), and the
# ` = ` assignment with its surrounding spaces (group 3).
ASSIGN_RE = re.compile(r"^(\s*)([A-Za-z]\w*)(\s*=\s*)")
SECTION_RE = re.compile(r"^\s*&(\w+)")
END_RE = re.compile(r"^\s*/\s*$")


def comment_start(line):
    """Index of the `!` that begins the comment, or -1 if none. Quote-aware so
    a `!` inside a quoted value is not mistaken for the comment."""
    quote = None
    for i, ch in enumerate(line):
        if quote:
            if ch == quote:
                quote = None
        elif ch in "'\"":
            quote = ch
        elif ch == "!":
            return i
    return -1


def replace_in_place(lines, name, value):
    """Rewrite the assignment line for `name`; return True if one was found."""
    for idx, raw in enumerate(lines):
        line = raw.rstrip("\n")
        m = ASSIGN_RE.match(line)
        if not (m and m.group(2) == name):
            continue
        head = m.group(1) + name + m.group(3)        # e.g. "    mp = "
        ci = comment_start(line)
        if ci >= 0:
            lines[idx] = head + value + "   " + line[ci:].rstrip() + "\n"
        else:
            lines[idx] = head + value + "\n"
        return True
    return False


def insert_in_group(lines, name, value, group, path):
    """Insert `    name = value` into the &group section, before its `/`."""
    in_group = False
    for idx, raw in enumerate(lines):
        line = raw.rstrip("\n")
        sm = SECTION_RE.match(line)
        if sm:
            in_group = (sm.group(1).lower() == group.lower())
            continue
        if in_group and END_RE.match(line):
            lines[idx:idx] = [f"    {name} = {value}\n", "\n"]
            return True
    sys.stderr.write(
        f"ERROR: set_nml_var.py: group '&{group}' not found in {path}, "
        f"cannot insert '{name}'\n")
    sys.exit(1)


def main():
    usage = "usage: set_nml_var.py <file.nml> <name> <value> [--group GROUP] [--insert]\n"
    args = sys.argv[1:]
    group = None
    insert = False
    if "--insert" in args:
        insert = True
        args.remove("--insert")
    if "--group" in args:
        gi = args.index("--group")
        try:
            group = args[gi + 1]
        except IndexError:
            sys.stderr.write(usage)
            sys.exit(2)
        del args[gi:gi + 2]
    if len(args) != 3:
        sys.stderr.write(usage)
        sys.exit(2)
    path, name, value = args

    if insert and group is None:
        sys.stderr.write("ERROR: set_nml_var.py: --insert requires --group\n")
        sys.exit(2)

    with open(path) as fh:
        lines = fh.readlines()

    if not replace_in_place(lines, name, value):
        if insert:
            insert_in_group(lines, name, value, group, path)
        else:
            sys.stderr.write(
                f"ERROR: set_nml_var.py: variable '{name}' not found in {path}\n"
                f"       (was it renamed or removed from the model namelist?)\n")
            sys.exit(1)

    with open(path, "w") as fh:
        fh.writelines(lines)


if __name__ == "__main__":
    main()
