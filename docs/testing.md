# Testing

> For the **automated CI/CD pipeline** (GitHub Actions workflows, the CPU
> full-test gate, the self-hosted GPU lane, regression/baseline handling, and the
> overall design), see [ci_cd_pipeline.md](ci_cd_pipeline.md). This page covers
> running tests locally.

## Running Test Cases

To run a standard test case, run

```bash
cd build/
make test_cases
```

The output of the test case will be written to the directory `tests/Test_Cases`

### Test Cases on SLURM

To run the automated test case script on a system using SLURM, a modification to the cmake step is needed. This is because systems running SLURM often restrict direct access to the mpiexec command, and require additional information such as which account to bill node-hours to. To setup the automated test case script to run on a SLURM system, pass the flag `-DSRUN_FLAGS` to the cmake command as such:

```bash
cd build/
cmake ../ -DSRUN_FLAGS='-A 9999'
```

This will setup test cases to run with srun using the user account "9999". Any additional flags which could be passed to srun can also be included here, **with the exception of -N and -n, which are automatically set by the test case script**

## Testing model components (developers)

Some basic integration testing of the model dynamics are implemented, and can be run by calling

```bash
cd build/
make check
```

tests are defined under the `tests/` folder. `test_driver.F90` manages the execution of the different test modules, which are defined as `test_XXXX.F90`. The runnable suite names (as passed to `HICAR-tester`) are: `advection`, `snow_drift`, `control_flow`, `halo_exch`, `geo`, `time`, `utilities` — note the advection suite is `advection`, though its source file is `test_advect.F90`. Run a single suite with `mpiexec -np 2 tests/HICAR-tester <suite>`, or run them all (with no argument) via `make check`.

## SNOWPACK C++ vs Fortran parity (developers)

The native-Fortran SNOWPACK port (`snowpack_driver.F90`, fetched with the
upstream `fortran-bindings` branch) is validated against the C++ SNOWPACK
bindings by `tests/snowpack/test_snowpack_compare.sh` (run in CI by
`.github/workflows/snowpack-compare.yml`): two HICAR builds differing only in
the snow driver run a 3 h seeded-snowpack case and all outputs are compared
against `tests/snowpack/tolerances_snowpack.yaml`.

A passing commit can be **blessed as the parity reference**
(`tests/snowpack/test_snowpack_compare.sh <repo> --bless`, or the workflow dispatch with
`bless=true`). This posts a `hicar-snowpack-parity-blessed=success` commit
status whose description records the upstream SNOWPACK SHA used — the anchor
for tracking upstream drift. Blessing is restricted to repo admins/maintainers
(role-checked on both paths) and requires the comparison evidence from the
current commit to say PASS. When the comparison later diverges (new C++
physics upstream not yet ported), run

```bash
tests/snowpack/snowpack_divergence_report.sh <hicar_repo>
```

to get the list of upstream C++ commits/diffs since the last parity bless —
the porting to-do list. Port them, re-run the comparison, re-bless.

Each comparison run leaves its evidence in `tests/figures/snowpack_compare/`
(per-variable stats in `parity_report.json`, spatial difference maps in
`diffmaps/`, provenance in `parity_meta.txt`). **Blessing archives these as a
permanent GitHub release** (`snowpack-parity/<date>-<sha>`), so the residual
level can be compared across blesses:

```bash
gh release list | grep snowpack-parity
gh release download snowpack-parity/<old> -p parity_report.json -O old.json
gh release download snowpack-parity/<new> -p parity_report.json -O new.json
tests/snowpack/parity_trend.py old.json new.json   # ranks per-variable residual growth
```

## Running the CI suites locally

The workflows are deliberately thin wrappers: every test is a standalone script
taking the repo path, and the build entrypoint CI uses
(`.github/scripts/hicar_install_utils.sh hicar_install`) runs locally too, so
you can reproduce a lane before triggering an Action. What the workflows add on
top (dependency install, caches, artifact upload, commit-status signing, the
bless Environment gate) is GitHub plumbing you don't need locally.

| Workflow job | Local equivalent (from the repo root) |
|---|---|
| full-test / cpu-debug | `HICAR_MODE=debug bash .github/scripts/hicar_install_utils.sh hicar_install`, then `cd build && mpiexec -np 4 tests/HICAR-tester`, `bash tests/Test_Cases/test_case_runner.sh . Standard`, `bash tests/test_reproducibility.sh . all` |
| full-test / cpu-release | `HICAR_MODE=release ... hicar_install`, then `bash tests/Test_Cases/test_case_runner.sh . "Standard,Standard_restart,Nested,Nested_restart"` and `bash tests/test_regression.sh . build "Standard,Nested" --mode exact` (resolves the blessed commit via gh) |
| snowpack-compare | build both exes (`HICAR_MODE=release HICAR_CMAKE_EXTRA="-DSNOWPACK_CPP=ON" ... hicar_install` for the C++ reference; again with `-DOPENACC=OFF` for the default native-Fortran port), then `bash tests/snowpack/test_snowpack_compare.sh . <cpp_exe> <fortran_exe>` |
| gpu | self-hosted only; with NVHPC+GPU locally: `bash tests/compare_cpu_gpu.sh ...` |

Using the same `hicar_install` entrypoint (with the same `HICAR_MODE` /
`HICAR_CMAKE_EXTRA`) is what makes the local run faithful to CI — only the
compiler/OS differ.

### Reproducing the GH runner (Docker)

To run
in the runner's environment, use the **CI-repro image** — a faithful clone of the
hosted CPU runner (`.github/docker/Dockerfile.ci-repro`):

```bash
# From the repo root. The image BUILD *is* the CI compile — if HICAR fails to
# build on CI, it fails here too, in the same toolchain. Default = cpu-debug:
docker build -f .github/docker/Dockerfile.ci-repro -t hicar-ci-repro .

# Reproduce another lane's build via build args:
docker build -f .github/docker/Dockerfile.ci-repro -t hicar-ci-repro-rel \
  --build-arg HICAR_MODE=release .
docker build -f .github/docker/Dockerfile.ci-repro -t hicar-ci-repro-cpp \
  --build-arg HICAR_CMAKE_EXTRA=-DSNOWPACK_CPP=ON .
```
