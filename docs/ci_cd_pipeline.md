# HICAR CI/CD pipeline — design & implementation

This document describes HICAR's automated test infrastructure: the philosophy, the
workflow architecture, the self-hosted GPU runner, and how to operate and extend
it. For *running tests locally*, see [testing.md](testing.md).

## 1. Test pyramid

| Layer | Question | Where it runs |
|---|---|---|
| Build matrix | Does it compile? | every PR (full-test) |
| Unit / invariant | Is a function correct in isolation? | CPU debug (full-test) |
| Smoke / integration | Does a config run N steps without NaN/abort? | CPU release (full-test) |
| Reproducibility | Decomposition- & restart-invariant? | CPU debug (full-test) |
| Regression (bit-for-bit) | Did output change vs a trusted baseline? | CPU release (full-test) |
| Memory safety | Any uninitialised / invalid memory access? | valgrind (**required on `main` PRs**, via the full-test chain) |
| SNOWPACK parity | Does the Fortran port match the C++ reference? | SNOWPACK-compare (**required on `main` PRs**) |
| CPU↔GPU equivalence | Does the GPU match the CPU within tolerance? | GPU lane (**required on `main` PRs**) |
| Validation | Is it physically right? | GPU lane (⬜ stub) |

**Key principle:** localize correctness *down* the pyramid (unit/invariant), use
regression as a catch-all tripwire, and reserve tolerance comparison for things
that cannot be bit-for-bit (CPU↔GPU, GNU↔NVFortran).

**Merge gate:** merging into `main` requires four of these lanes to have *signed*
the PR head commit (`hicar-full-test`, `valgrind`, `gpu-check`, `snow-parity`).

## 2. Workflow architecture (the "signed-gate" model)

Each of the GitHub Actions acts as a gate, signing a given commit with a success
status if it clears that test. This lets us see which failures a given commit
has/had. Merging into `main` requires all four (`hicar-full-test`, `valgrind`, `gpu-check`, `snow-parity`) to be present and green on the PR head commit.

### 2.1 `hicar-full-test.yml` — the CPU correctness gate

GNU / CPU, GitHub-hosted. 

1. Runs unit tests and reproducibility tests using
a debug build. This way, errors in unit tests are caught with full debug
output, and reproducibility is assured by avoiding compiler optimizations that
may break it.
2. Runs the full end-to-end test suite, then compiles the last blessed commit
and performs and integration test using that commit as the baseline.

Triggers: **`pull_request` → `main`** (the gate — runs on
the *test-merge* and signs the PR head), **`push` to `main`/`develop`** (develop
feedback; the `main` push is the **fresh post-merge run** that signs `main`'s own
HEAD), and **workflow_dispatch**.

### 2.2 `gpu.yml` — the GPU test suite

NVFortran, self-hosted. 

1. Compiles two exes, one using OpenACC and one without, both
as release builds using the NVFortran compiler. Then runs the Standard test case
using both exes, and checks that their results are within tolerance
2. **Validation** ⬜ — stubbed (idealized case / conservation checks to come).

Triggers: **`workflow_dispatch` only** — **no `push` and no
`pull_request` trigger**, so a fork PR can never execute on the
self-hosted box.

### 2.3 `bless-baseline.yml` — bless the regression reference commit

Maintainer-only (`workflow_dispatch`), **role-checked** — the dispatching actor must
be a repo `admin`/`maintain` (collaborator-permission API). Posts a **`hicar-regression-blessed=success` commit status** to a (known-good, ideally full-test-passed) commit — the same machinery as the full-test`sign`.

### 2.4 `snowpack-compare.yml` — SNOWPACK C++ vs Fortran parity

HICAR uses a fortran port of snowpack to allow for porting snowpack to GPUs using OpenACC. This
was necessary because the c++ snowpack code uses dynamically sized arrays, which OpenACC does not
handle well. This test is used to ensure that the fortran port of snowpack is still similar to
(in parity with) the c++ code. The test builds HICAR twice on the same gfortran toolchain, where one exe uses the fortran wrappers for a c++ library of the snowpack model, and another uses the native fortran snowpack part. It then runs a 3 h seeded-snowpack comparison (`tests/scripts/snowpack/test_snowpack_compare.sh`)
against `tests/tolerances/tolerances_snowpack.yaml`.

Triggers: **every
`pull_request` → `main`** (no path filter — runs on the test-merge), **`push` to
`main`** (the post-merge run), and manual dispatch.

A **`sign-snow`** job run at the end of the workflow caches the SHA of the upstream
snowpack repo we used to run the test, so we can later see which snowpack commit
was the last one to give parity. Artifacts of the test including the difference log
and difference plots of a few key variables are also archived here, to be viewed as
artifacts of the test.

When the comparison starts failing:

```
tests/scripts/snowpack/snowpack_divergence_report.sh <hicar_repo>
```

Then port the listed changes into the Fortran driver and re-run the comparison
until it passes; the next passing run re-records the anchor automatically.

### 2.5 `valgrind-memcheck.yml` — uninitialised-memory check (required on `main` PRs)

GNU / CPU, GitHub-hosted. Triggers: **`workflow_dispatch`** and **`workflow_run`** —
it fires automatically when **HICAR full-test (CPU)** completes as part of a PR into main. 

A CPU **debug** build runs the unit-test suite (`HICAR-tester`) under 
`valgrind --track-origins=yes`. Scope is deliberately narrow — the debug 
`HICAR-tester` only: integration cases
are far too slow under valgrind, and the GPU/OpenACC build can't run under it at
all.

## 3. Comparison engine

`tests/scripts/compare_outputs.py` is the single comparison engine for all lanes:

- `--mode exact` (bit-for-bit, CPU regression) and `--mode tolerance`
  (CPU↔GPU), driven by per-variable `rtol`/`atol` in `tests/tolerances.yaml`.
- **NaN/Inf introduced by the candidate is always a failure**
- A dropped output variable → failure; an added one → warning.
- On failure it localizes the **worst offender** (variable, `dim=index`, both
  values, abs/rel error), emits GitHub `::error::` annotations, and writes a
  JSON report.

## 4. Operating & extending

- **Add a unit test**: create `tests/test_<name>.F90` with a `collect_<name>_suite`,
  add it to the `tests` list in `tests/CMakeLists.txt` **and** to the `testsuites`
  array in `tests/unit/test_driver.F90`.
- **Tune cross-lane tolerances**: edit `tests/tolerances.yaml` (per-variable
  `rtol`/`atol`); start loose, tighten as the GPU comparison spread is learned.
- **Gate a `main` PR**: full-test, snowpack-parity, and valgrind sign the PR head
  automatically; only GPU is manual — a maintainer reviews the diff and runs
  `gh workflow run gpu.yml --ref <head>` to sign `gpu-check`.
