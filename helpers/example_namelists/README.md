# Example namelists

Minimal, illustrative HICAR namelists showing how to configure common run
types. **The `*.nml` files here are auto-generated — do not edit them by hand.**

## How it works

Each example is defined by a small find/replace script in `gen_example_nmls/`.
`generate_examples.sh` turns each one into a `<name>.nml`:

1. Generate a fresh default namelist from the built executable
   (`./bin/HICAR --gen-nml`), so the defaults always track the current build.
2. Apply `gen_example_nmls/<name>.sh` (sed find/replace) to set the entries the
   example customizes — the same pattern the test-case generators use
   (`tests/Test_Cases/input/nml_gen_scripts/`).
3. Strip every entry still at its default (`strip_defaults.py`), leaving a
   compact namelist that shows only what the example changes. The one exception
   is any variable whose name contains `folder` (the output / restart folder
   paths): these are always kept, even at their default, so each example
   documents where its output and restart files are written.

This regeneration runs automatically in CI on every push to `main` / `develop`
— as the `example-namelists` job in `.github/workflows/hicar-full-test.yml`,
which reuses the full-test release build instead of compiling HICAR a second
time — and commits the refreshed `.nml` files back, so they never drift from
the current option set or defaults. (The job is a sibling of the test jobs: it
runs whether or not the tests pass, and does not gate them.)

## Adding an example

Add `gen_example_nmls/<name>.sh` (copy an existing one). It receives the default
namelist path as `$1` and edits it in place; keep using the
`sed -i'.bak' "s/old/new/g" $1` form. Then run:

```bash
helpers/example_namelists/generate_examples.sh        # uses ./bin/HICAR
# or:  helpers/example_namelists/generate_examples.sh /path/to/HICAR
```

`<name>.nml` appears here. (CI will also regenerate it on push.)

### Find/replace gotchas

* **Anchor patterns to avoid substring matches.** `start_date` is a substring of
  `restart_date`, so target it as `" start_date = ''"` (leading space). Check any
  short variable name against the full default namelist before adding it.
* Use `\/` to escape slashes in paths inside the `sed` expression.

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
