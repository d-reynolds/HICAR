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
| Memory safety | Any uninitialised / invalid memory access? | valgrind (nightly **+ required on `main` PRs**) |
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

**Nightly no-op skip.** The three lanes that run on the nightly `schedule` (GPU,
SNOWPACK-compare, valgrind) share a small `changed` gate job: it compares the
commit under test against the head SHA of that workflow's most recent completed
run on `main`, and skips the whole run when `main` has not advanced since — there
is no point re-testing an idle `main` every night. Only the scheduled trigger is
gated; manual dispatch and event-driven (PR/push) runs always run.

### 2.1 `hicar-full-test.yml` — the CPU correctness gate ✅

GNU / CPU, GitHub-hosted. Triggers: **push** and **pull_request** (`main`/`develop`),
**workflow_dispatch**, and **workflow_call** (reusable — see GPU lane).

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

NVFortran, self-hosted. Triggers: **nightly (04:00 UTC)** and
**workflow_dispatch** — **no `push` and no `pull_request` trigger** (see Security,
§4), so a fork PR can never execute on the self-hosted box. The nightly runs
against `main`; `develop` is GPU-tested *after* it merges to `main`, and *before*
merge via the manual `gpu-check` gate (dispatch on the head branch — §2.6).

0. **`changed`** (hosted): on the nightly schedule, compares the commit under
   test against the head SHA of the most recent completed run of this workflow on
   `main`. If they match (main hasn't advanced since the last GPU run), the whole
   pipeline is skipped — no point re-testing an unchanged `main`. Manual dispatch
   always runs.
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
4. **`sign-gpu`** (hosted, runs only if `gpu-tests` passed): posts the
   **`gpu-check = success`** commit status that `main`'s ruleset requires (§2.6).
   The GPU-vs-CPU tolerance comparison is currently *advisory*
   (`continue-on-error`), so `gpu-check` signs on a successful NVHPC **build + run**
   of both exes. **To gate a `main` PR**, dispatch this workflow on the PR's head
   branch — `gh workflow run gpu.yml --ref develop` — so `github.sha` is the PR head
   commit and the status lands where the ruleset evaluates it.

Only the `HICAR` target is built (`cmake --build --target HICAR` + `cmake
--install`), so the GPU lane does **not** pay the slow `HICAR-tester` device-link.

### 2.3 `bless-baseline.yml` — bless the regression reference commit ✅

Maintainer-only (`workflow_dispatch`). Posts a **`hicar-regression-blessed=success`
commit status** to a (known-good, ideally full-test-passed) commit — the same
machinery as the full-test `sign`. The regression reference is therefore "the most
recent blessed commit," with no stored file, no PR, no build, no LFS, no token.
Warns if the commit hasn't passed full-test.

### 2.4 `snowpack-compare.yml` — SNOWPACK C++ vs Fortran parity, with bless ✅

Builds HICAR twice on the same gfortran toolchain (`-DSNOWPACK_CPP=ON` C++
bindings vs the default native-Fortran port), each into its **own build tree**
(`build_cpp` / `build_fortran`, via `HICAR_BUILD_DIR`). It then runs a 3 h seeded-snowpack comparison
(`tests/snowpack/test_snowpack_compare.sh`) against
`tests/snowpack/tolerances_snowpack.yaml`. Triggers: PRs touching the snow drivers
(paths-filtered), nightly, manual.

A **`sign-snow`** job posts the **`snow-parity = success`** commit status that
`main`'s ruleset requires (§2.6) whenever the comparison passes. It signs
automatically on a PR that touches the snow drivers; for a `main` PR that does
**not** touch them, dispatch the workflow on the head branch
(`gh workflow run snowpack-compare.yml --ref develop`) to post the status.

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
   with `{"reviewers":[{"type":"User","id":<numeric-uid>}]}`.

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

### 2.5 `valgrind-memcheck.yml` — uninitialised-memory check ✅ (required on `main` PRs)

GNU / CPU, GitHub-hosted. Triggers: **nightly (03:30 UTC)**, **workflow_dispatch**,
and **`workflow_run`** — it fires automatically when **HICAR full-test (CPU)**
completes. Its `changed` gate runs the memcheck only when that full-test **passed**
*and* the tested commit has an **open PR to `main`** (it checks out the full-test's
`head_sha`). This makes valgrind a **required merge check** (§2.6) — kept off every
push (valgrind is ~10–50× slower than native), but run automatically once a commit
is actually a `main`-merge candidate.

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
full valgrind log is uploaded as an artifact.

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
| `hicar-full-test` | full-test `sign` | **auto** — push to `develop` (or a PR) passes the CPU suite |
| `valgrind` | valgrind `sign-valgrind` | **auto** — `workflow_run` after full-test passes, when an open `main` PR exists |
| `gpu-check` | gpu `sign-gpu` | **manual** — dispatch `gpu.yml` on the head branch |
| `snow-parity` | snowpack-compare `sign-snow` | **auto** on snow-driver PRs; **manual** dispatch otherwise |

Typical `develop → main` merge:

1. Push to `develop` → full-test runs and signs `hicar-full-test`; valgrind then
   auto-fires (`workflow_run`) and signs `valgrind`.
2. `gh workflow run gpu.yml --ref develop` → signs `gpu-check`.
3. `gh workflow run snowpack-compare.yml --ref develop` (only if the PR didn't
   already auto-run it) → signs `snow-parity`.
4. All four green on the head commit → the PR merges.

Statuses are evaluated per-commit, so each must land on the **current** PR head:
push a new commit and you must re-run the manual gates on it. Repo admins keep a
ruleset **bypass** (`RepositoryRole: always`) for emergencies.

> **Bootstrapping note.** A lane only becomes triggerable once its workflow file is
> on `main` — both `workflow_dispatch` and `workflow_run` resolve from the default
> branch. The first PR that *introduces* a new required lane therefore cannot run it
> yet; merge that one via admin bypass (or temporarily drop the new context from the
> ruleset), after which every later PR runs it normally.

### 2.7 Flow diagram

```
                 push develop / PR
                        │
                        ▼
                 hicar-full-test ──(both pass)──> sign  hicar-full-test=success
                        │
                        └──(success)──> [workflow_run] valgrind ──> changed? (open main PR)
                                                          └──(clean)──> sign-valgrind  valgrind=success

   manual dispatch on the PR head branch (--ref develop):
     gpu.yml          ─> changed ─> check-signed ─┬─ signed ──> gpu-tests ─(pass)─> sign-gpu  gpu-check=success
                                                  └─ unsigned ─> full-test ─pass─> gpu-tests ─> sign-gpu
     snowpack-compare ─> compare ─(pass)─> sign-snow  snow-parity=success   (auto on snow-driver PRs)

   merge to main  ⇐  { hicar-full-test, valgrind, gpu-check, snow-parity } all green on the PR head

   ── nightly (each gated by `changed`: skip if main unchanged) ──
     snowpack-compare 03:00 · valgrind 03:30 · gpu 04:00
     snowpack FAIL ─> snowpack_divergence_report.sh ─> diff anchor..upstream ─> port ─> re-bless
```

---

## 3. Comparison engine ✅

`tests/compare_outputs.py` is the single comparison engine for all lanes:

- `--mode exact` (bit-for-bit, CPU regression) and `--mode tolerance`
  (CPU↔GPU), driven by per-variable `rtol`/`atol` in `tests/tolerances.yaml`.
- **NaN/Inf introduced by the candidate is always a failure**
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
the repo; blessing just posts a commit status (§2.3).

**CPU↔GPU** (`tests/compare_cpu_gpu.sh`) runs the GPU exe and a CPU-host NVHPC exe
(same commit, same compiler) with an identical namelist and compares by tolerance.

---

## 4. Self-hosted GPU runner & security

The GPU lane runs on a containerized self-hosted runner (see
[.github/docker/README.md](https://github.com/HICAR-Model/HICAR/blob/main/.github/docker/README.md)):
NVHPC + CUDA + NCCL with
an NVFortran-built NetCDF/HDF5/PnetCDF stack. HICAR is built per-run (not baked
into the image) because `-gpu=ccnative` needs the device visible at build time.

**Security posture: the runner never executes untrusted PR code.**

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
  array in `tests/test_driver.F90`.
- **Tune cross-lane tolerances**: edit `tests/tolerances.yaml` (per-variable
  `rtol`/`atol`); start loose, tighten as the nightly GPU spread is learned.
- **Gate a `main` PR** (sign the head commit with the manual checks): after
  full-test + valgrind have auto-signed, run `gh workflow run gpu.yml --ref <head>`
  and (if the PR didn't auto-run it) `gh workflow run snowpack-compare.yml --ref
  <head>`. See [§2.6](#26-merge-gate--required-checks-on-main).
- **The integration runner lives in the data repo**: `tests/Test_Cases/` (including
  `test_case_runner.sh`) is cloned from `HICAR-Model/Test-Data`, not this repo — so
  changes to how integration cases are *run* (e.g. the hang/wall-clock watchdog) go
  there, while the build/compare/regression scripts under `tests/` are here.
