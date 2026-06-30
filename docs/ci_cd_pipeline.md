# HICAR CI/CD pipeline — design & implementation

This document describes HICAR's automated test infrastructure: the philosophy, the
workflow architecture, the self-hosted GPU runner, and how to operate and extend
it. For *running tests locally*, see [testing.md](testing.md).

> Status legend: ✅ implemented · 🟡 partial / needs validation · ⬜ planned

---

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
the PR head commit (`hicar-full-test`, `valgrind`, `gpu-check`, `snow-parity`) —
see [§2.6](#26-merge-gate--required-checks-on-main).

---

## 2. Workflow architecture (the "signed-gate" model)

Five GitHub Actions workflows: the CPU **full-test** gate, the **GPU** suite, the
**SNOWPACK** parity comparison, the **valgrind** memcheck, and the maintenance
**bless** workflow.

Each gating lane records its result as a **commit status** on the SHA it tested
(`hicar-full-test`, `valgrind`, `gpu-check`, `snow-parity`). Because a status is
attached to a commit — not to a workflow run — it can be posted by any trigger
(push, PR, manual dispatch, or `workflow_run`) and is then honored as a green check
wherever that commit appears. Merging into `main` requires all four to be present
and green on the PR head commit; see [§2.6](#26-merge-gate--required-checks-on-main).

### 2.1 `hicar-full-test.yml` — the CPU correctness gate ✅

GNU / CPU, GitHub-hosted. Triggers: **`pull_request` → `main`** (the gate — runs on
the *test-merge* and signs the PR head), **`push` to `main`/`develop`** (develop
feedback; the `main` push is the **fresh post-merge run** that signs `main`'s own
HEAD), and **workflow_dispatch**.

- **`cpu-debug`** (`MODE=debug`, bounds + `finit-real=nan`): unit/invariant tests
  (`HICAR-tester`) and reproducibility (decomposition + restart).
- **`cpu-release`** (`MODE=release`): runs all integration configs (Standard,
  Nested, restarts) — *run-only* — then a separate **regression** step that
  recompiles the blessed commit and diffs the integration outputs against it
  (see §3). Integration ("did it run?") and regression ("is the output correct?")
  are deliberately split.
- **`sign`** (runs only if both pass): posts a GitHub **commit status**
  `hicar-full-test = success` to the commit SHA using the built-in `GITHUB_TOKEN`
  (no PAT needed for same-repo statuses).

The commit status *is* the "signature" that a commit passed full-test — durable,
queryable, and shown as a green check on the commit.

### 2.2 `gpu.yml` — the GPU test suite ✅ / 🟡

NVFortran, self-hosted. Triggers: **`workflow_dispatch` only** — **no `push` and no
`pull_request` trigger** (see Security, §4), so a fork PR can never execute on the
self-hosted box. A change is GPU-tested *before* merging a PR into main via the manual `gpu-check` gate
(a maintainer dispatches on the PR's head branch — §2.6), and `main` itself is
GPU-tested *after* merge (a maintainer dispatches on `main`).


1. **`gpu-tests`** (self-hosted):
   - Builds a **CPU reference exe** with NVHPC, host (`-DOPENACC=OFF` →
     `bin/HICAR`) and a **GPU exe** with NVHPC (`-DOPENACC=ON` → `bin/HICAR_gpu`).
     Same compiler for both isolates the GPU port from compiler differences.
   - Runs both on the Standard case with an **identical namelist** and compares
     GPU-vs-CPU by tolerance (`tests/scripts/compare_cpu_gpu.sh`).
   - **Validation** ⬜ — stubbed (idealized case / conservation checks to come).
2. **`sign-gpu`** (hosted, runs only if `gpu-tests` passed): posts the
   **`gpu-check = success`** commit status that `main`'s ruleset requires (§2.6).
   The GPU-vs-CPU tolerance comparison is currently *advisory*
   (`continue-on-error`), so `gpu-check` signs on a successful NVHPC **build + run**
   of both exes. **To gate a `main` PR**, dispatch this workflow on the PR's head
   branch — `gh workflow run gpu.yml --ref <your-feature-branch>` — so `github.sha` is the PR head
   commit and the status lands where the ruleset evaluates it.

Only the `HICAR` target is built (`cmake --build --target HICAR` + `cmake
--install`), so the GPU lane does **not** pay the slow `HICAR-tester` device-link.

### 2.3 `bless-baseline.yml` — bless the regression reference commit ✅

Maintainer-only (`workflow_dispatch`), **role-checked** — the dispatching actor must
be a repo `admin`/`maintain` (collaborator-permission API). Posts a **`hicar-regression-blessed=success` commit status** to a
(known-good, ideally full-test-passed) commit — the same machinery as the full-test
`sign`. The regression reference is therefore "the most
recent blessed commit," with no stored file, no PR, no build, no LFS, no token.
Warns if the commit hasn't passed full-test.

### 2.4 `snowpack-compare.yml` — SNOWPACK C++ vs Fortran parity ✅

Builds HICAR twice on the same gfortran toolchain (`-DSNOWPACK_CPP=ON` C++
bindings vs the default native-Fortran port), each into its **own build tree**
(`build_snowpack_cpp` / `build_snowpack_fortran`, via `HICAR_BUILD_DIR`). It then
runs a 3 h seeded-snowpack comparison (`tests/scripts/snowpack/test_snowpack_compare.sh`)
against `tests/tolerances/tolerances_snowpack.yaml`. Triggers: **every
`pull_request` → `main`** (no path filter — runs on the test-merge), **`push` to
`main`** (the post-merge run), and manual dispatch.

" A **`sign-snow`** job posts
the **`snow-parity = success`** commit status that `main`'s ruleset requires (§2.6)
whenever the comparison passes (on a PR on the PR head; on a post-merge `push` on
`main`'s HEAD). Its **description records the upstream SNOWPACK SHA**
(`snowpack=<sha>`) of the fetched `fortran-bindings` checkout, so the success status
*doubles as the divergence anchor* — the last upstream state at which the two
implementations provably matched. 
**Evidence archive (90-day artifact).** Every comparison run writes a
machine-readable per-variable stats report (`parity_report.json`:
max|abs|/max|rel|/violation counts for all variables common to both output files,
passing ones included — the exact count depends on the run's `output_vars` and is
reported at runtime as `n_compared`) and the spatial difference maps (`diffmaps/`)
into `tests/figures/snowpack_compare/`, uploaded as the
`snowpack-compare-diagnostics` workflow artifact (retained ~90 days). To track
residual growth, download the report from two runs (`gh run download <run-id> -n
snowpack-compare-diagnostics`) and analize it with `tests/scripts/snowpack/parity_trend.py
<old.json> <new.json>`. The upstream SNOWPACK commit each run was built against is
recorded in the `snow-parity` status description (`snowpack=<sha>`), not in the
artifact.

**Divergence routine.** SNOWPACK is fetched from the *moving* `fortran-bindings`
branch (cmake/FindSNOWPACK.cmake), so upstream C++ commits land in CI builds
automatically. When the comparison starts failing:

```
tests/scripts/snowpack/snowpack_divergence_report.sh <hicar_repo>
```

resolves the anchor from the most recent `snow-parity=success` status (reading its
`snowpack=<sha>`), deepens the shallow SNOWPACK clone, and prints `git log`/`diff
--stat` between the anchor and the current upstream SHA — split into the C++ core
(`snowpack/`, `meteoio/`: the porting to-do list) and `fortran/` (port-side
changes). Port the listed changes into the Fortran driver and re-run the comparison
until it passes; the next passing run re-records the anchor automatically.

### 2.5 `valgrind-memcheck.yml` — uninitialised-memory check ✅ (required on `main` PRs)

GNU / CPU, GitHub-hosted. Triggers: **`workflow_dispatch`** and **`workflow_run`** —
it fires automatically when **HICAR full-test (CPU)** completes. The `memcheck` job's
`if:` runs only when that full-test **passed** *and* was either a **PR into `main`**
(signs the PR head) or a **post-merge `push` to `main`** (signs `main`'s HEAD);
`develop` pushes are skipped (it reads the triggering run's event + branch and checks
out its `head_sha`). This makes valgrind a
**required merge check** (§2.6) — kept off `develop` pushes (valgrind is ~10–50×
slower than native), but run automatically on every `main`-merge candidate and on
`main` itself.

A CPU **debug** build runs the unit-test suite (`HICAR-tester`) under `valgrind
--track-origins=yes`. This lane exists because the debug build's
`-finit-real/-integer/-logical` flags only poison **local/stack** variables — they
do **not** initialise heap / allocatable memory. An uninitialised `logical`
component of an *allocatable* options object therefore passed every local test
(reading garbage-false) yet failed deterministically on the amd64 CI runner
(reading garbage-true); only valgrind reliably catches that class of bug in
Fortran (gfortran has no MemorySanitizer).

Scope is deliberately narrow — the debug `HICAR-tester` only: integration cases
are far too slow under valgrind, and the GPU/OpenACC build can't run under it at
all. The gate is robust to third-party noise: instead of `--error-exitcode`
(mpich is not valgrind-clean), it **fails only when a valgrind error block cites a
HICAR `.F90` source**, and also fails if the tester itself exits non-zero. The
full valgrind log is uploaded as an artifact. The gate logic lives in
`tests/scripts/test_valgrind.sh` (runnable locally as `make test_valgrind`), which the
workflow's `run:` step invokes.

A **`sign-valgrind`** job posts the outcome as the **`valgrind`** commit status —
`success` on a clean run, `failure` otherwise, so a flagged PR shows a red required
check rather than a perpetually "Expected" one.

### 2.6 Merge gate — required checks on `main`

A PR into `main` can merge only when **four commit-status checks** are present and
green on the PR **head commit**, enforced by the repository ruleset *"Require pull
request for main"* (`required_status_checks`). Each is posted by a `sign` job on the
SHA it tested:

| Status context | Posted by | How it's produced |
|---|---|---|
| `hicar-full-test` | full-test `sign` | **auto** — the PR into `main` runs on the test-merge and passes the CPU suite |
| `valgrind` | valgrind `sign-valgrind` | **auto** — `workflow_run` after the PR's full-test passes |
| `gpu-check` | gpu `sign-gpu` | **manual** — a maintainer dispatches `gpu.yml` on the PR head branch (self-hosted; never auto-triggered — §4) |
| `snow-parity` | snowpack-compare `sign-snow` | **auto** — runs on every `main` PR |

Typical feature-branch → `main` merge:

1. Open a PR from your feature branch into `main`. full-test and snowpack-compare run
   on the **test-merge** and sign the **PR head** (`hicar-full-test`, `snow-parity`);
   valgrind auto-fires (`workflow_run`) after full-test and signs the PR head.
2. A maintainer reviews the diff and dispatches GPU on the PR head branch
   (`gh workflow run gpu.yml --ref <your-feature-branch>`) → signs `gpu-check`.
3. All four green on the PR head commit → the PR merges.
4. The merge **pushes to `main`**, firing a **fresh post-merge run** (full-test +
   snowpack + valgrind on `main`; GPU via a manual dispatch) that
   signs `main`'s **own HEAD**. `main` carries its own statuses rather than
   inheriting the PR head's — statuses are per-SHA and the merge commit is a new SHA.

Statuses are evaluated per-commit, so each must land on the **current** PR head:
push a new commit and you must re-run the manual gate (GPU) on it. Repo admins keep a
ruleset **bypass** (`RepositoryRole: always`) for emergencies.

> **Bootstrapping note.** A lane only becomes triggerable once its workflow file is
> on `main` — both `workflow_dispatch` and `workflow_run` resolve from the default
> branch. The first PR that *introduces* a new required lane therefore cannot run it
> yet; merge that one via admin bypass (or temporarily drop the new context from the
> ruleset), after which every later PR runs it normally.

### 2.7 Flow diagram

```
   PR  feature ─> main   (tests run on the TEST-MERGE, sign the PR HEAD)
     hicar-full-test ─(pass)─> sign  hicar-full-test=success
            └─(success)─> [workflow_run] valgrind ─(clean)─> sign-valgrind  valgrind=success
     snowpack-compare ─(pass)─> sign-snow  snow-parity=success      (auto, every main PR)
     gpu.yml  ── MANUAL dispatch on the PR head branch (maintainer, after review) ──
              ─> gpu-tests ─(pass)─> sign-gpu  gpu-check=success

   merge to main  ⇐  { hicar-full-test, valgrind, gpu-check, snow-parity } all green on the PR HEAD

   push to main (the merge) ─> FRESH post-merge run signs main's OWN head:
     hicar-full-test + snowpack-compare + valgrind     (GPU via manual dispatch)

   push to develop ─> hicar-full-test only (feedback)

   snowpack FAIL ─> snowpack_divergence_report.sh ─> diff anchor..upstream ─> port ─> re-bless
```

---

## 3. Comparison engine ✅

`tests/scripts/compare_outputs.py` is the single comparison engine for all lanes:

- `--mode exact` (bit-for-bit, CPU regression) and `--mode tolerance`
  (CPU↔GPU), driven by per-variable `rtol`/`atol` in `tests/tolerances.yaml`.
- **NaN/Inf introduced by the candidate is always a failure**
- A dropped output variable → failure; an added one → warning.
- On failure it localizes the **worst offender** (variable, `dim=index`, both
  values, abs/rel error), emits GitHub `::error::` annotations, and writes a
  JSON report.

**Regression** (`tests/scripts/test_regression.sh`) uses a **blessed *commit*** as the
reference, not a stored output file. The reference is the most recent commit
carrying a `hicar-regression-blessed=success` status (resolved by
`tests/scripts/resolve_blessed_commit.sh` — branch-aware: the newest blessed *ancestor* of
the commit under test). It recompiles that commit's exe (cached by hash in a git
worktree), runs it on the requested cases to **regenerate** the reference output,
and diffs the *current* integration outputs against it. Both exes are built on the
same runner/toolchain, so the comparison defaults to **bit-for-bit (`--mode
exact`)** — any difference is a real change in model output. Nothing is stored in
the repo; blessing just posts a commit status (§2.3).

**CPU↔GPU** (`tests/scripts/compare_cpu_gpu.sh`) runs the GPU exe and a CPU-host NVHPC exe
(same commit, same compiler) with an identical namelist and compares by tolerance.

---

## 4. Self-hosted GPU runner & security

The GPU lane runs on a containerized self-hosted runner (see
[.github/docker/README.md](https://github.com/HICAR-Model/HICAR/blob/main/.github/docker/README.md)):
NVHPC + CUDA + NCCL with
an NVFortran-built NetCDF/HDF5/PnetCDF stack. HICAR is built per-run (not baked
into the image) because `-gpu=ccnative` needs the device visible at build time.

**Security posture: the runner never executes untrusted PR code.**

- `gpu.yml` has **no `push`/`pull_request` trigger** — only manual
  `workflow_dispatch` by an org maintainer. It is **never auto-triggered** by a PR.
- The runner is **ephemeral** and run as a **fresh container per job**
  (`run-runner.sh`), so no state persists between jobs.
- Runs as non-root, **no Docker socket mounted** (no host escape).
- The GPU job uses no repo secrets; the runner-admin PAT is scrubbed from the
  environment before any job step runs.
- Recommended host hygiene: dedicated box, isolated network with egress firewall,
  `chmod 600 runner.env`.

---

## 5. Operating & extending

- **Add a unit test**: create `tests/test_<name>.F90` with a `collect_<name>_suite`,
  add it to the `tests` list in `tests/CMakeLists.txt` **and** to the `testsuites`
  array in `tests/unit/test_driver.F90`.
- **Tune cross-lane tolerances**: edit `tests/tolerances.yaml` (per-variable
  `rtol`/`atol`); start loose, tighten as the GPU comparison spread is learned.
- **Gate a `main` PR**: full-test, snowpack-parity, and valgrind sign the PR head
  automatically; only GPU is manual — a maintainer reviews the diff and runs
  `gh workflow run gpu.yml --ref <head>` to sign `gpu-check`. See
  [§2.6](#26-merge-gate--required-checks-on-main).
- **The integration runner lives in the data repo**: `tests/Test_Cases/` (including
  `test_case_runner.sh`) is cloned from `HICAR-Model/Test-Data`, not this repo — so
  changes to how integration cases are *run* (e.g. the hang/wall-clock watchdog) go
  there, while the build/compare/regression scripts under `tests/` are here.
