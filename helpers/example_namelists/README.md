# Example namelists

Minimal, illustrative HICAR namelists showing how to configure common run
types. **The `*.nml` files here are auto-generated — do not edit them by hand.**

## How it works

Each example is defined by a small script in `gen_example_nmls/`.
`generate_examples.sh` turns each one into a `<name>.nml`:

1. Generate a fresh default namelist from the built executable
   (`./bin/HICAR --gen-nml`), so the defaults always track the current build.
2. Apply `gen_example_nmls/<name>.sh` to set the entries the example
   customizes. Each line is `set_var <name> <value>`, which rewrites that
   variable's line **by name**, whatever its current default. The scripts
   therefore never encode the model defaults: a default value can change and
   the example still gets the value it asked for. If a variable no longer
   exists (renamed/removed) the script fails (exit 1) and the build reports it,
   instead of silently leaving the default in place.
3. Strip every entry still at its default (`strip_defaults.py`), leaving a
   compact namelist that shows only what the example changes. The one exception
   is any variable whose name contains `folder` (the output / restart folder
   paths): these are always kept, even at their default, so each example
   documents where its output and restart files are written.

   The strip step keeps the documentation intact: the model's own header (the
   block `--gen-nml` writes above the first group) is reproduced verbatim, each
   kept entry carries its **complete** comment including any multi-line wrapped
   description, and all inline comments are re-aligned so the `!` marker starts
   at a fixed column (`COMMENT_COL` in `strip_defaults.py`), or just past the
   value if the value runs beyond it.

This regeneration runs automatically in CI on every push to `main` / `develop`
— as the `example-namelists` job in `.github/workflows/hicar-full-test.yml`,
which reuses the full-test release build instead of compiling HICAR a second
time — and commits the refreshed `.nml` files back, so they never drift from
the current option set or defaults. (The job is a sibling of the test jobs: it
runs whether or not the tests pass, and does not gate them.)

## Adding an example

Add `gen_example_nmls/<name>.sh` (copy an existing one). It receives the default
namelist path as `$1`, defines the `set_var` helper, and then lists one
`set_var <name> <value> <group>` per entry the example customizes (`<group>` is
the namelist group the variable lives in, e.g. `physics`, `domain`, `forcing`).
Then run:

```bash
helpers/example_namelists/generate_examples.sh        # uses ./bin/HICAR
# or:  helpers/example_namelists/generate_examples.sh /path/to/HICAR
```

`<name>.nml` appears here. (CI will also regenerate it on push.)

### Writing `set_var` lines

* `set_var <name> <value> <group>`. The value is passed **verbatim** — no `sed`
  escaping. Quote string values (`set_var mp "'morrison'" physics`), leave
  numbers/logicals bare (`set_var nz 40 domain`, `set_var RK3 .True.
  time_parameters`), give per-nest values as comma lists
  (`set_var dx "250.0, 150.0" domain`), and write paths with normal slashes
  (`set_var init_conditions_file "'../domains/dom.nc'" domain`).
* `<group>` is the namelist group the variable belongs to (the `&group` header
  it appears under in `--gen-nml`). For the example scripts it is documentary —
  every variable is already present in the full default, so `set_var` finds and
  rewrites it. (The integration-test generators reuse `set_nml_var.py` with an
  extra `--insert` flag, which uses `<group>` to add a variable that the
  stripped example omitted.)
* Matching is anchored on the whole variable name, so substrings never collide
  (`start_date` is not confused with `restart_date`, nor `output` with
  `output_vars`) — you do not need the old leading-space anchoring trick.
* A `set_var` for a variable the model does not have fails the build, so a typo
  or a stale name is caught rather than silently ignored.

## Why empty sections are kept

Generated examples keep a bare `&section /` header for any section left fully at
defaults rather than deleting it. HICAR reads its always-present sections
(`general`, `restart`, `domain`, `forcing`, `physics`, `time_parameters`,
`output`) unconditionally; a *missing* group makes the namelist read hit
end-of-file and **stop** the model. An empty *present* group reads fine and
leaves every value at its in-code default.

## Full option reference

These examples are intentionally minimal. For the complete, annotated option
list and per-variable docs:

```bash
./bin/HICAR --gen-nml default.nml      # full annotated namelist
./bin/HICAR -v <variable>              # docs for one option
./bin/HICAR -v --all                   # docs for every option
```
