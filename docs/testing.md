# Testing

> For the **automated CI/CD pipeline** (GitHub Actions workflows, the CPU
> full-test gate, the self-hosted GPU lane, regression/baseline handling, and the
> overall design), see [ci_cd_pipeline.md](ci_cd_pipeline.md). This page covers
> running tests locally.

HICAR uses a tiered approach to testing, comprising unit, integration, and end-to-end tests.

Unit tests are small in their scope and seek to test discrete "units" of code. HICAR implements
unit testing via the [test-drive](https://github.com/fortran-lang/test-drive) fortran project.
The fortran code in tests/unit contain these unit tests, grouped by code region. Agruably, 
some of these are not strictly unit tests and verge on integration testing, but if it ain't broke,
don't fix it.

Integration tests use code components together, how well they "integrate" with one another. As mentioned,
some of the "unit" tests are more like integration tests. The tests/scripts directory contains automated
testing, some of which, like the cli tests or decomposition testing, checks that various model components
produce expected results when working together.

Lastly, end-to-end testing checks that a full model run produces expected results. This is obviously difficult
for atmospheric models, where there is no ground truth for a high-resolution, 3D atmospheric state. Instead,
we use end-to-end testing to check that the model does not crash, that we have parity between different setups
which should be similar/identical, and that we do not drift from a known baseline simulation. This last point
is known as regression testing, where we want to ensure that the output does not "regress" relative to the baseline.
Finally, and most expensive, is validation testing, where the model is run for long enough that we can spin up its 
internal state, and compare simulation output to observations. This is currently an outstanding to do in the 
automated test setup, and will be implemented on the GPU self-hosted runner for improved time-to-solution.

## Running the tests (`make` targets)

Every test is exposed as a `make` target in the build directory, so the usual
pattern is just:

```bash
cd build/
make <target>
```

Each target **pulls in its own build and data dependencies** — it (re)builds the
executable it needs and downloads the test data if missing — so the command is
self-contained; you don't have to build or fetch data first.

| Target | What it checks | Extra needs |
|---|---|---|
| `make check` | Unit / invariant suites via `HICAR-tester` (4 MPI ranks) | — |
| `make test_valgrind` | `HICAR-tester` under valgrind `--track-origins` (uninitialised / invalid memory) | valgrind + a debug build |
| `make test_cli` | Every `./bin/HICAR` command-line option behaves | — |
| `make test_cases` | Integration cases run without crash/NaN (Standard, Nested, restarts) | test data (auto) |
| `make test_decomposition` | Output is identical at 5 vs 10 MPI ranks | builds a debug exe |
| `make test_restart` | A restarted run matches an uninterrupted one | builds a debug exe |
| `make test_reproducibility` | Decomposition + restart together | builds a debug exe |
| `make test_regression` | Bit-for-bit vs the blessed reference commit (runs the cases, then diffs) | `gh` + a blessed commit |
| `make test_snowpack` | Native-Fortran SNOWPACK port matches the C++ build (builds both, compares) | builds both CPU exes |
| `make test_gpu` | GPU output matches CPU within tolerance (builds both, compares) | NVHPC + a GPU |

Notes:
- **SLURM:** `make check` / `make test_cases` need `-DSRUN_FLAGS='-A <account>'` at the cmake step (see [Test Cases on SLURM](#test-cases-on-slurm)).
- **`test_valgrind`** runs `make check`'s suite under `valgrind --track-origins` (the CPU memcheck CI lane, runnable locally). Needs a debug build (`-DMODE=debug`) and `valgrind` installed; it writes `build/valgrind.log` and fails on any error block citing a HICAR `.F90` (or a non-zero tester exit).
- **`test_regression`** resolves its reference via `gh` and is soft/skipped until a commit is blessed (see [ci_cd_pipeline.md](ci_cd_pipeline.md)). It reuses existing integration output under `tests/Test_Cases/output` and only re-runs the cases when missing. To pin a reference or use tolerance mode, call the script directly: `bash tests/scripts/test_regression.sh . build "<cases>" --blessed-commit <sha> --mode tolerance`.
- **`test_gpu`** needs an NVHPC-configured build (`-DFC=nvfortran`) and a visible GPU — it can't build the GPU exe on a CPU-only checkout.

The sections below cover the most-used targets in more detail.

## Running Test Cases

To run a standard test case, run

```bash
cd build/
make test_cases
```

The output of the test case will be written to the directory `tests/Test_Cases`

### Test Cases on SLURM

On HPC systems that use SLURM, direct calls to `mpiexec`/`mpirun` are often
restricted — MPI jobs must instead be launched through `srun`, which also needs to
know which account to bill the node-hours to. The two MPI-launching test targets,
**`make check`** (the unit/integration suite) and **`make test_cases`** (the
integration `.nml` cases), support this through a single CMake cache variable,
`SRUN_FLAGS`. Pass it at the cmake step:

```bash
cd build/
cmake ../ -DSRUN_FLAGS='-A 9999'
```

This makes both targets launch with `srun`, billing the account `9999`. Any other
flag you would normally pass to `srun` can be added here too — **with the exception
of `-N` and `-n` (the node and task counts), which the test targets set
themselves.**

`SRUN_FLAGS` only takes effect once CMake's MPI detection has resolved the launcher
(`MPIEXEC_EXECUTABLE`) to `srun`. If `mpiexec`/`mpirun` is on your `PATH` and CMake
selects it instead, point CMake at `srun` explicitly at configure time:

```bash
cd build/
cmake ../ -DMPIEXEC_EXECUTABLE="$(command -v srun)" -DSRUN_FLAGS='-A 9999'
```

`SRUN_FLAGS` is a configure-time setting, so re-run the cmake step (as above) in an
existing build directory to add or change it.

## Testing model components

Some basic integration testing of the model dynamics are implemented, and can be run by calling

```bash
cd build/
make check
```

tests are defined under `tests/unit/`. `test_driver.F90` manages the execution of the different test modules, which are defined as `test_XXXX.F90`. The runnable suite names (as passed to `HICAR-tester`) are: `advection`, `snow_drift`, `control_flow`, `halo_exch`, `geo`, `time`, `utilities` — note the advection suite is `advection`, though its source file is `test_advect.F90`. Run a single suite with `mpiexec -np 2 tests/HICAR-tester <suite>`, or run them all (with no argument) via `make check`.

## SNOWPACK C++ vs Fortran parity (developers)

The native-Fortran SNOWPACK port (`snowpack_driver.F90`, fetched with the
upstream `fortran-bindings` branch) is validated against the C++ SNOWPACK
bindings by `tests/scripts/snowpack/test_snowpack_compare.sh` (run in CI by
`.github/workflows/snowpack-compare.yml`): two HICAR builds differing only in
the snow driver run a 3 h seeded-snowpack case and all outputs are compared
against `tests/tolerances/tolerances_snowpack.yaml`.

If the github workflow action passes, the comparison posts a `snow-parity=success` commit status
whose description records the upstream SNOWPACK SHA used (`snowpack=<sha>`); that
status doubles as the **divergence anchor** (the last upstream state at which the
two implementations matched). If the comparison later diverges (new C++ physics
upstream not yet ported), run

```bash
tests/scripts/snowpack/snowpack_divergence_report.sh <hicar_repo>
```

to get the list of upstream C++ commits/diffs since the last passing parity run —
the porting to-do list. Port them and re-run the comparison; the next passing run
re-records the anchor automatically.

Each comparison run leaves its evidence in `tests/figures/snowpack_compare/`
(per-variable stats in `parity_report.json`, spatial difference maps in
`diffmaps/`), uploaded as the `snowpack-compare-diagnostics` workflow artifact
(retained ~90 days). To compare the residual level between two runs:

```bash
gh run list --workflow=snowpack-compare.yml --branch main      # find the run IDs
gh run download <old-run-id> -n snowpack-compare-diagnostics -D old/
gh run download <new-run-id> -n snowpack-compare-diagnostics -D new/
tests/scripts/snowpack/parity_trend.py old/parity_report.json new/parity_report.json
```

### Reproducing the GH runner (Docker)

To run
in the runner's environment, use the **CI-repro image** — a faithful clone of the
hosted CPU runner (`.github/docker/Dockerfile.ci-repro`):

```bash
# From the repo root. The image BUILD is the CI compile — if HICAR fails to
# build on CI, it fails here too, in the same toolchain. Default = cpu-debug:
docker build -f .github/docker/Dockerfile.ci-repro -t hicar-ci-repro .

```
