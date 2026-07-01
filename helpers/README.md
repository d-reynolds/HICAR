# helpers/

Utility scripts for setting up and running HICAR.

| Script | Purpose | Docs |
|---|---|---|
| `gen_HICAR_dir.sh` | Create a working-directory tree (input/output/restart/forcing/domains) and populate the supporting files for a run. | [tutorial](../docs/tutorial.md) |
| `filelist_script.sh` | Build a forcing-file list (one absolute quoted path per line) from a glob. | [forcing_data](../docs/forcing_data.md) |
| `batch_submit_SLURM.sh` | Self-resubmitting SLURM job chain for runs longer than the wall-time: each job restarts HICAR from the newest restart file (`restart_run`/`restart_date` set automatically) and queues a successor that only fires if the job is killed. Test locally with `--setup-only`. | header comments |

## domains/ — domain-generation toolchain

| File | Purpose |
|---|---|
| `gen_HICAR_dom.py` | Main driver: generate a HICAR domain file from a DEM + land-use data (terrain descriptors, horizon matrix / sky-view factor via HORAYZON). Edit the paths at the top, then run. |
| `HICAR_Domain.py`, `ProjHelpers.py` | Modules used by `gen_HICAR_dom.py`. |
| `HICAR_dom.yml` | Conda environment for the toolchain (`conda env create -f HICAR_dom.yml`). |
| `ds_cutter.py` | Cut/subset an existing domain dataset (handles staggered variables). |
| `merge_fields.py` | Merge numbered field variants (`name1..nameN`) in a dataset into combined variables. |
| `Regrid_script.sh` | `ESMF_Regrid` wrapper to regrid variables onto a destination grid. |

Setup and usage: [domain_generation](../docs/domain_generation.md). The
HORAYZON dependency is cloned/installed per those instructions (a local
`domains/HORAYZON/` checkout is gitignored).
