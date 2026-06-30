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

This regeneration runs automatically in CI on every push to `main` / `develop`
— as the `example-namelists` job in `.github/workflows/hicar-full-test.yml`

## Full option reference

These examples are intentionally minimal. For the complete, annotated option
list

```bash
./bin/HICAR --gen-nml default.nml      # full annotated namelist
```
