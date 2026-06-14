# HICAR CI/CD pipeline — design & implementation

This document describes HICAR's automated test infrastructure: the philosophy, the
workflow architecture, the self-hosted GPU runner, and how to operate and extend
it. For *running tests locally*, see [testing.md](testing.md).

> Status legend: ✅ implemented · 🟡 partial / needs validation · ⬜ planned

---

## 1. Goals

1. **Don't crash.** Every supported configuration should run without aborting or
   producing NaN/Inf.
2. **Behave as expected.** Detect unintended changes in model output (regression),
   and verify the GPU port matches the CPU.
3. **Informative failures.** A failed test should say *what* broke, *where*, and
   *why* — localized to a variable, grid cell, and tolerance where possible.
4. **No test code in the production hot path.** Tests live in separate
   translation units / build targets; runtime assertions are compile-guarded.
5. **Cover CPU (GNU) and GPU (NVFortran)** without trying to test the full
   Cartesian product of physics options.

Compilers in scope: **GNU** (CPU) and **NVFortran/NVHPC** (CPU host + GPU). Intel
and Cray are deprecated.

---

## 2. Test pyramid

| Layer | Question | Where it runs |
|---|---|---|
| Build matrix | Does it compile? | every PR (full-test) |
| Unit / invariant | Is a function correct in isolation? | CPU debug (full-test) |
| Smoke / integration | Does a config run N steps without NaN/abort? | CPU debug+release |
| Reproducibility | Decomposition- & restart-invariant? | CPU debug (full-test) |
| Regression (bit-for-bit) | Did output change vs a trusted baseline? | CPU release (full-test) |
| CPU↔GPU equivalence | Does the GPU match the CPU within tolerance? | GPU lane |
| Validation | Is it physically right? | GPU lane (⬜ stub) |

**Key principle:** localize correctness *down* the pyramid (unit/invariant), use
regression as a catch-all tripwire, and reserve tolerance comparison for things
that cannot be bit-for-bit (CPU↔GPU, GNU↔NVFortran).

---

## 3. Workflow architecture (the "signed-gate" model)

Two GitHub Actions workflows plus a maintenance workflow.

### 3.1 `hicar-full-test.yml` — the CPU correctness gate ✅

GNU / CPU, GitHub-hosted. Triggers: **pull_request**, **nightly (02:00 UTC)**,
**workflow_dispatch**, and **workflow_call** (reusable — see GPU lane).

- **`cpu-debug`** (`MODE=debug`, bounds + `finit-real=nan`): unit/invariant tests
  (`HICAR-tester`), Standard integration case (run-only), reproducibility
  (decomposition + restart).
- **`cpu-release`** (`MODE=release`): runs all integration configs (Standard,
  Nested, restarts) — *run-only* — then a separate **regression** step that
  recompiles the blessed commit and diffs the integration outputs against it
  (see §4). Integration ("did it run?") and regression ("is the output correct?")
  are deliberately split.
- **`sign`** (runs only if both pass): posts a GitHub **commit status**
  `hicar-full-test = success` to the commit SHA using the built-in `GITHUB_TOKEN`
  (no PAT needed for same-repo statuses).

The commit status *is* the "signature" that a commit passed full-test — durable,
queryable, and shown as a green check on the commit.

### 3.2 `gpu.yml` — the GPU test suite ✅ / 🟡

NVFortran, self-hosted. Triggers: **push to main/develop** (+ a temporary
`ci_pipeline` shakeout trigger), **nightly (04:00 UTC**, after the CPU nightly),
**workflow_dispatch**. **No `pull_request` trigger** (see Security, §5).

1. **`check-signed`** (hosted): queries the commit status; outputs `signed`.
2. **`full-test`** (reusable call to `hicar-full-test.yml`): runs **only if the
   commit is not signed**. If it fails, the GPU job does not run.
3. **`gpu-tests`** (self-hosted, runs if signed OR full-test just passed):
   - Builds a **CPU reference exe** with NVHPC, host (`-DOPENACC=OFF` →
     `bin/HICAR`) and a **GPU exe** with NVHPC (`-DOPENACC=ON` → `bin/HICAR_gpu`).
     Same compiler for both isolates the GPU port from compiler differences.
   - Runs both on the Standard case with an **identical namelist** and compares
     GPU-vs-CPU by tolerance (`tests/compare_cpu_gpu.sh`).
   - **Validation** ⬜ — stubbed (idealized case / conservation checks to come).

Only the `HICAR` target is built (`cmake --build --target HICAR` + `cmake
--install`), so the GPU lane does **not** pay the slow `HICAR-tester` device-link.

### 3.3 `bless-baseline.yml` — bless the regression reference commit ✅

Maintainer-only (`workflow_dispatch`). Posts a **`hicar-regression-blessed=success`
commit status** to a (known-good, ideally full-test-passed) commit — the same
machinery as the full-test `sign`. The regression reference is therefore "the most
recent blessed commit," with no stored file, no PR, no build, no LFS, no token.
Warns if the commit hasn't passed full-test.

### 3.4 `snowpack-compare.yml` — SNOWPACK C++ vs Fortran parity, with bless ✅

Builds HICAR twice on the same gfortran toolchain (`-DSNOWPACK_CPP=ON` C++
bindings vs the default native-Fortran port), runs a 3 h seeded-snowpack
comparison (`tests/snowpack/test_snowpack_compare.sh`) against
`tests/snowpack/tolerances_snowpack.yaml`. Triggers: PRs touching the snow
drivers, nightly, manual.

**Parity blessing.** Analogous to the regression bless, but with a second,
independent status context and one extra payload: a manual dispatch with
`bless=true` posts **`hicar-snowpack-parity-blessed=success`** on the commit,
and the status *description* records the **upstream SNOWPACK commit SHA**
(`snowpack=<sha>`) of the fetched `fortran-bindings` checkout used in the run.
That SHA is the *parity anchor*: the last upstream state at which the Fortran
port provably matched the C++ build.

Blessing is **restricted to repo admins/maintainers** and **gated on a passing
comparison**, with three stacked controls:

1. **Role check** (both paths): the actor's role is queried via the
   collaborator-permission API (`role_name` admin or maintain). The workflow
   fails fast on dispatch with `bless=true` before building; the local
   `--bless` refuses before posting. `--force` never bypasses authorization.
2. **Pass gate**: in CI the bless job runs only if the comparison job
   succeeded (`needs:`); locally `--bless` checks the run's `parity_meta.txt`
   evidence and refuses unless it says PASS, was produced from the current
   HEAD, and matches the SNOWPACK checkout (`--force` overrides the evidence
   checks only, and is recorded in the status text).
3. **Required-reviewer Environment** (platform-enforced): the bless runs in a
   separate job bound to `environment: parity-bless`. With required reviewers
   configured, GitHub pauses that job until a designated reviewer approves it
   in the Actions UI — independent of anything in the workflow file or
   scripts. **One-time setup** (repo admin):
   `Settings > Environments > New environment "parity-bless" > Required
   reviewers > add the owners/maintainers`, or via API:
   `gh api -X PUT repos/<owner>/<repo>/environments/parity-bless --input -`
   with `{"reviewers":[{"type":"User","id":<numeric-uid>}]}`. NB: until the
   reviewers are configured, the auto-created environment does NOT pause the
   job; and environment protection rules require a public repo (or a paid
   plan for private repos).

Residual caveat: GitHub itself gates the statuses API at write access, so a
write-level user could in principle post the status with a raw API call
outside this tooling — controls 1 and 3 guard the sanctioned paths, which is
where accidental or casual blessing happens.

**Evidence archive.** Every comparison run writes a machine-readable
per-variable stats report (`parity_report.json`: max|abs|/max|rel|/violation
counts for ALL 282 variables, passing ones included), the standard spatial
difference maps (`diffmaps/`), and a provenance stamp into
`tests/figures/snowpack_compare/` (uploaded as a workflow artifact, ≤90 d).
**Blessing additionally publishes these as a permanent GitHub release**
tagged `snowpack-parity/<date>-<sha>`, so the parity level is consultable
across blesses indefinitely: list with `gh release list`, download two
reports and rank residual growth with
`tests/snowpack/parity_trend.py <old.json> <new.json>` (flags >3× growth as a
candidate regression and points at the divergence routine below).

**Divergence routine.** SNOWPACK is fetched from the *moving*
`fortran-bindings` branch (cmake/FindSNOWPACK.cmake), so upstream C++ commits
land in CI builds automatically. When the nightly comparison starts failing:

```
tests/snowpack/snowpack_divergence_report.sh <hicar_repo>
```

resolves the anchor from the blessed status, deepens the shallow SNOWPACK
clone, and prints `git log`/`diff --stat` between the anchor and the current
upstream SHA — split into the C++ core (`snowpack/`, `meteoio/`: the porting
to-do list) and `fortran/` (port-side changes). Port the listed changes into
the Fortran driver, re-run the comparison until it passes, then re-bless
(workflow dispatch with `bless=true`, or locally
`tests/snowpack/test_snowpack_compare.sh <repo> --bless --reason '...'`). Every
comparison run also logs the SNOWPACK SHA it used, so CI logs alone can
bracket when a divergence appeared.

### 3.5 Flow diagram

```
PR / nightly ──> hicar-full-test ──(both pass)──> sign commit (hicar-full-test=success)

GPU nightly ──> check-signed ──signed?──┬── yes ──> gpu-tests (CPU ref + GPU, compare)
                                        └── no  ──> full-test ──pass──> gpu-tests

snowpack-compare nightly ──pass──> (manual) bless: hicar-snowpack-parity-blessed
                          ──FAIL──> snowpack_divergence_report.sh
                                      └─> diff anchor..current upstream ─> port ─> re-bless
```

---

## 4. Comparison engine ✅

`tests/compare_outputs.py` is the single comparison engine for all lanes:

- `--mode exact` (bit-for-bit, CPU regression) and `--mode tolerance`
  (CPU↔GPU), driven by per-variable `rtol`/`atol` in `tests/tolerances.yaml`.
- **NaN/Inf introduced by the candidate is always a failure** (the previous
  engine used `nanmax` and silently passed NaNs).
- A dropped output variable → failure; an added one → warning.
- On failure it localizes the **worst offender** (variable, `dim=index`, both
  values, abs/rel error), emits GitHub `::error::` annotations, and writes a
  JSON report.

**Regression** (`tests/test_regression.sh`) uses a **blessed *commit*** as the
reference, not a stored output file. The reference is the most recent commit
carrying a `hicar-regression-blessed=success` status (resolved by
`tests/resolve_blessed_commit.sh` — branch-aware: the newest blessed *ancestor* of
the commit under test). It recompiles that commit's exe (cached by hash in a git
worktree), runs it on the requested cases to **regenerate** the reference output,
and diffs the *current* integration outputs against it. Both exes are built on the
same runner/toolchain, so the comparison defaults to **bit-for-bit (`--mode
exact`)** — any difference is a real change in model output. Nothing is stored in
the repo; blessing just posts a commit status (§3.3).

**CPU↔GPU** (`tests/compare_cpu_gpu.sh`) runs the GPU exe and a CPU-host NVHPC exe
(same commit, same compiler) with an identical namelist and compares by tolerance.

---

## 5. Self-hosted GPU runner & security

The GPU lane runs on a containerized self-hosted runner (see
[.github/docker/README.md](../.github/docker/README.md)): NVHPC + CUDA + NCCL with
an NVFortran-built NetCDF/HDF5/PnetCDF stack. HICAR is built per-run (not baked
into the image) because `-gpu=ccnative` needs the device visible at build time.

**Security posture: the runner never executes untrusted PR code.**

- `gpu.yml` has **no `pull_request` trigger** — only push to main/develop,
  nightly, and dispatch. A fork PR cannot cause execution on the hardware.
- The runner is **ephemeral** and run as a **fresh container per job**
  (`run-runner.sh`), so no state persists between jobs.
- Runs as non-root, **no Docker socket mounted** (no host escape).
- The GPU job uses no repo secrets; the runner-admin PAT is scrubbed from the
  environment before any job step runs.
- Recommended host hygiene: dedicated box, isolated network with egress firewall,
  `chmod 600 runner.env`.

To GPU-test a contributor PR: a maintainer reviews the diff, pulls it onto an
in-repo branch (push to `develop`), which triggers the lane on vetted code.

---

## 6. The configuration-combinatorics problem

We do **not** test the full Cartesian product of physics options. Instead:

1. Push correctness into unit/invariant tests so a module is validated
   independently of the others. ⬜ (mostly planned — see §8)
2. A curated set of representative configs (Standard, Nested, restarts). ✅
3. ⬜ A pairwise-covering config set for nightly (planned, Phase 4).
4. Tier by cadence: fast subset per-PR, broader matrix nightly.

---

## 7. Implementation phases & status

- **Phase 0 — plumbing** ✅: comparison engine, `tolerances.yaml`, stored-baseline
  regression + bless workflow, containerized GPU runner, the workflow set.
- **Phase 1 — localize** ⬜: register the ~15 currently-unbuilt unit tests (in
  `tests/CMakeLists.txt` *and* the `testsuites` array in `test_driver.F90`); add
  `test_conservation.F90` / `test_invariants.F90` (water/mass budgets, limiter
  monotonicity, zero-wind→zero-tendency).
- **Phase 2 — GPU correctness** 🟡: CPU↔GPU tolerance comparison implemented;
  per-kernel equivalence tests planned.
- **Phase 3 — verify** ⬜: idealized/analytic cases (e.g. solid-body cone
  advection) — currently the GPU "validation" stub.
- **Phase 4 — breadth** ⬜: pairwise config matrix, tiered PR vs nightly.
- **Phase 5 — performance** ⬜: per-kernel timing regression.

---

## 8. Operating & extending

- **Bless the regression reference**: run `bless-baseline.yml` (optionally pass a
  commit; blank = default-branch HEAD). It posts a `hicar-regression-blessed`
  commit status — no PR, no file. Bless a commit that passed full-test (it warns
  if it hasn't). Until at least one ancestor is blessed, the regression step is
  skipped/`continue-on-error`.
- **Split of run vs check**: `test_case_runner.sh` only *runs* integration cases;
  `test_regression.sh` does the comparison against the blessed commit.
- **Trigger the GPU lane manually**: `gh workflow run gpu.yml` (once on the default
  branch), or push to `develop`.
- **Add a unit test**: create `tests/test_<name>.F90` with a `collect_<name>_suite`,
  add it to the `tests` list in `tests/CMakeLists.txt` **and** to the `testsuites`
  array in `tests/test_driver.F90`.
- **Tune cross-lane tolerances**: edit `tests/tolerances.yaml` (per-variable
  `rtol`/`atol`); start loose, tighten as the nightly GPU spread is learned.
- **Build-time note**: HICAR's NVHPC GPU build is dominated by the relocatable
  device-code link (`acc routine` across files); ~7 min is the floor and is not
  reducible by optimization level. Build test/CPU-reference exes at `-O0` /
  without `deepcopy`.

---

## 9. Known issues / assumptions to validate

- 🟡 The CPU reference build (`-DOPENACC=OFF`, NVFortran host) must compile and run
  HICAR with GPU-only paths (NCCL/cuFFT) gated off.
- 🟡 `compare_cpu_gpu.sh` runs both exes with `RAD=rrtmgp`; if rrtmgp won't run on
  the host build, set `RAD=rrtmg` (keep it identical for both exes).
- 🟡 At least one ancestor commit must be blessed (`hicar-regression-blessed`
  status) for the regression step to run; until then it's skipped. The regression
  rebuilds the blessed commit when it changes (cached by hash in
  `bin/HICAR_blessed`), which adds a HICAR build to the release job.
- The temporary `ci_pipeline` push trigger in `gpu.yml` must be removed before
  merging to main.
