# HICAR self-hosted GPU CI runner

A containerized GitHub Actions runner that provides the NVHPC + CUDA + NCCL
toolchain, so the
`gpu.yml` workflow can build and run HICAR on a real GPU.

## Why a container

* **Reproducible & redeployable** — `docker compose up -d` on any GPU box gives
  an identical toolchain; no hand-installed dependencies to drift.
* **Ephemeral & safe for a public repo** — each job runs on a fresh registration
  and the runner exits afterward (`--ephemeral`), so a malicious PR cannot
  persist on the hardware. Combined with `gpu.yml` only triggering on the
  `run-gpu` label (maintainers) or the nightly schedule, untrusted code never
  auto-runs here.

## One-time host setup

1. Install the NVIDIA driver and the **NVIDIA Container Toolkit** so Docker can
   expose GPUs (`docker run --rm --gpus all nvcr.io/nvidia/nvhpc:... nvidia-smi`).
2. Log in to NGC if needed to pull the NVHPC base image
   (`docker login nvcr.io`).

## Deploy

```bash
cd .github/docker
cp runner.env.example runner.env       # fill in REPO_URL and GH_PAT (chmod 600 it)
docker compose build                   # build the image (deps w/ nvfortran, ~30–60 min, once)

# Production run: FRESH container per job (recommended).
./run-runner.sh                        # foreground; Ctrl-C to stop
# ...or keep it alive across reboots via the systemd unit:
#   edit paths in hicar-gpu-runner.service, then
#   sudo cp hicar-gpu-runner.service /etc/systemd/system/
#   sudo systemctl enable --now hicar-gpu-runner
#   journalctl -u hicar-gpu-runner -f
```

The runner appears under **Settings → Actions → Runners** with labels
`self-hosted, gpu, nvhpc`.

> `docker compose up -d` also works for quick **dev** runs, but it re-execs the
> runner in the *same* container (writable layer persists between jobs).
> `run-runner.sh` uses `docker run --rm` so every job starts from the pristine
> image — prefer it in production.

## Security model

Self-hosted runners on a public repo are a known risk: a triggered workflow runs
the code at that ref on your hardware. This setup is hardened so that **only code
a maintainer has already vetted ever runs**:

* **No untrusted PR code, ever.** `gpu.yml` has **no `pull_request` trigger** — it
  runs only on push to `master`/`develop`, the nightly schedule, and manual
  `workflow_dispatch`. A fork PR cannot cause execution on the box. Contributor
  PRs get GPU validation only after review + merge (or a maintainer dispatching
  the workflow against a ref they've read).
* **Ephemeral + fresh container per job** — `--ephemeral` runner + `docker run
  --rm` loop: one job per container, then destroyed. No persistence between jobs.
* **Non-root, no Docker socket** — jobs run as user `runner`; the host Docker
  socket is never mounted, so there's no trivial container→host escape.
* **No secrets in the GPU job** — `gpu.yml` sets `permissions: contents: read`
  and uses no repo secrets. The runner-admin PAT is **scrubbed from the
  environment** (`unset GH_PAT` in entrypoint.sh) before any job step runs, so a
  job cannot read it from `env`.
* **Recommended host hygiene** — run on a dedicated box, ideally an isolated
  network segment with an **egress firewall** (limits exfiltration / lateral
  movement in the unlikely event a trusted build is compromised). Keep
  `runner.env` mode `600`; it holds the registration PAT.

Also set, in **Settings → Actions → General → Fork pull request workflows**:
"Require approval for all outside collaborators" (defense in depth even though no
PR trigger exists here).

### Validating a contributor PR on GPU (the cost of this posture)

Because there's no PR trigger, GPU-testing an external PR is a deliberate,
maintainer-initiated act **after reading the diff** (including `.github/` and
build scripts). Pull the PR onto a branch *in this repo* — which only a
maintainer can push to — and the `push` trigger runs it on the GPU box:

```bash
gh pr checkout <PR#>                    # review the FULL diff first!
git push origin HEAD:develop           # push to develop -> triggers gpu.yml
# (or push to a throwaway in-repo branch and: gh workflow run gpu.yml --ref <that-branch>)
```

Never point the runner at the fork's ref directly — staging it in-repo is the
review gate.

## Values to confirm for YOUR GPU box

The defaults target **NVHPC 26.3** `cuda_multi` (the nvfortran that builds HICAR
cleanly — 26.1 has an internal compiler error on `output_obj.F90`). Match this to
the `nvfortran --version` of your known-good build. Adjust in `docker-compose.yml`
(the args propagate to the Dockerfile):

| Arg          | What it controls                                   | How to check |
|--------------|----------------------------------------------------|--------------|
| `NVHPC_TAG`  | NVHPC base image (cuda_multi picks CUDA ≤ driver)  | [NGC nvhpc tags](https://catalog.ngc.nvidia.com/orgs/nvidia/containers/nvhpc/tags); match your `nvfortran --version` |
| `NVHPC_VER`  | Path `…/Linux_x86_64/<VER>` inside the image       | `ls /opt/nvidia/hpc_sdk/Linux_x86_64` in the image |

NCCL and MPI use NVHPC's version-agnostic `comm_libs/{nccl,mpi}` symlinks, so
there's no separate CUDA arg to track. Confirm `NVHPC_VER` once the image is pulled:
```bash
docker run --rm --entrypoint bash nvcr.io/nvidia/nvhpc:<TAG> -c \
  'ls /opt/nvidia/hpc_sdk/Linux_x86_64'
```

`-gpu=ccnative` means HICAR is compiled for the GPU present at build time — the
build happens inside the running container (which has GPU access), so no
compute-capability arg is needed.

## Known risks to validate on first build

These can't be verified without an NVHPC + GPU host; expect to iterate once:

* **NetCDF C stack with `gcc` vs `nvc`** — built with `OMPI_CC=gcc` for safety.
  If linking HICAR complains about ABI, rebuild netcdf-c with `nvc`.
* **NCCL path** — `NCCL_ROOT` points at `comm_libs/nccl`; HICAR's CMake
  auto-detects it. If NCCL isn't found, the build falls back to GPU-aware MPI.

## Updating

```bash
docker compose build --pull     # rebuild on dependency or NVHPC bumps
docker compose up -d
```
