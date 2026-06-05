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
cp runner.env.example runner.env       # fill in REPO_URL and GH_PAT
docker compose build                   # builds deps with nvfortran (~30–60 min, once)
docker compose up -d
docker compose logs -f                 # watch it register; should print nvidia-smi -L
```

The runner appears under **Settings → Actions → Runners** with labels
`self-hosted, gpu, nvhpc`.

## Values to confirm for YOUR GPU box

The defaults target the latest **NVHPC 26.x** `cuda_multi` image (best guess —
confirm the exact tag string on NGC). Adjust in `docker-compose.yml` (the args
propagate to the Dockerfile):

| Arg          | What it controls                                   | How to check |
|--------------|----------------------------------------------------|--------------|
| `NVHPC_TAG`  | NVHPC base image (cuda_multi picks CUDA ≤ driver)  | [NGC nvhpc tags](https://catalog.ngc.nvidia.com/orgs/nvidia/containers/nvhpc/tags); `nvidia-smi` for driver CUDA |
| `NVHPC_VER`  | Path `…/Linux_x86_64/<VER>` inside the image       | `ls /opt/nvidia/hpc_sdk/Linux_x86_64` in the image |
| `NVHPC_CUDA` | Path `…/comm_libs/<CUDA>/nccl`                      | `ls .../comm_libs` in the image |

Quickest way to confirm VER/CUDA once the image is pulled:
```bash
docker run --rm --entrypoint bash nvcr.io/nvidia/nvhpc:<TAG> -c \
  'ls /opt/nvidia/hpc_sdk/Linux_x86_64 && ls /opt/nvidia/hpc_sdk/Linux_x86_64/*/comm_libs'
```

`-gpu=ccnative` means HICAR is compiled for the GPU present at build time — the
build happens inside the running container (which has `--gpus all`), so no
compute-capability arg is needed.

## Known risks to validate on first build

These can't be verified without an NVHPC + GPU host; expect to iterate once:

* **NetCDF C stack with `gcc` vs `nvc`** — built with `OMPI_CC=gcc` for safety.
  If linking HICAR complains about ABI, rebuild netcdf-c with `nvc`.
* **NCCL path** — set via `NVHPC_CUDA`; HICAR's CMake auto-detects via
  `NCCL_ROOT`. If NCCL isn't found, the build falls back to GPU-aware MPI.

## Updating

```bash
docker compose build --pull     # rebuild on dependency or NVHPC bumps
docker compose up -d
```
