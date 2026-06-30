# The High-resolution Intermediate Complexity Atmospheric Research Model (HICAR)

<!-- The four merge-gate lanes that sign main's HEAD: hicar-full-test, valgrind, gpu-check, snow-parity -->
[![full-test](https://github.com/HICAR-Model/HICAR/actions/workflows/hicar-full-test.yml/badge.svg?branch=main)](https://github.com/HICAR-Model/HICAR/actions/workflows/hicar-full-test.yml)
[![valgrind](https://github.com/HICAR-Model/HICAR/actions/workflows/valgrind-memcheck.yml/badge.svg?branch=main)](https://github.com/HICAR-Model/HICAR/actions/workflows/valgrind-memcheck.yml)
[![GPU](https://github.com/HICAR-Model/HICAR/actions/workflows/gpu.yml/badge.svg?branch=main)](https://github.com/HICAR-Model/HICAR/actions/workflows/gpu.yml)
[![SNOWPACK parity](https://github.com/HICAR-Model/HICAR/actions/workflows/snowpack-compare.yml/badge.svg?branch=main)](https://github.com/HICAR-Model/HICAR/actions/workflows/snowpack-compare.yml)

[![License: GPLv3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)
[![version](https://img.shields.io/github/v/tag/HICAR-Model/HICAR?filter=v*&label=version)](https://github.com/HICAR-Model/HICAR/tags)
[![DOI](https://zenodo.org/badge/638935780.svg)](https://zenodo.org/badge/latestdoi/638935780)
[![Paper: GMD](https://img.shields.io/badge/paper-GMD-1f7a8c.svg)](https://doi.org/10.5194/gmd-2023-16)
[![Docs](https://img.shields.io/badge/docs-mkdocs-blue.svg)](docs/index.md)

HICAR is a variant of the Intermediate Complexity Atmospheric Research (ICAR) model developed for sub-kilometer resolutions. The model is developed for downscaling of kilometer-scale NWP model output to resolutions used for land-surface simulations. HICAR features physics parameterizations shared by traditional weather models such as WRF, but with massively simplified dynamics which enable for run times nearly 600x faster than WRF.

More information about the model features enabling this can be found in Reynolds et al., 2023.

## Getting Started

To start working with the model, follow the instructions under docs/compiling to compile the model. 

Once the model has been compiled, a working directory file structure, supporting files, and test case can be setup and downloaded with the following script:

```bash
./helpers/gen_HICAR_dir.sh /path/to/working/directory /path/to/HICAR/repo
```

This script creates the working-directory tree (`input/`, `output/`, `restart/`, `forcing/`, `domains/`) and copies the supporting files (the NoahMP `.TBL` file and the `mp_support`/`rrtmg_support`/`rrtmgp_support` directories) from the repo's `run/` directory into `input/`. It runs non-interactively and does not overwrite existing directories; the supporting files must already have been downloaded by the cmake configure step.

#### Compilation Requirements
While being fast to run compared to traditional weather models, HICAR has still been developed and intended for use on High Performance Computing (HPC) machines, although it can also be run on a local machine. HICAR thus uses a few package requirements which are common to HPC environments. They are:

- Parallel NetCDF4
- MPI
- FFTW

HICAR currently supports either the GNU Fortran or NVFortran compilers.

#### Static data requirements

HICAR uses a domain file which defines land-surface variables and some terrain-descriptors. To generate a HICAR domain file using an existing DEM and land use data, the python script gen_HICAR_dom.py, located in helpers/domains/, can be used. See below for more details on how to run this script.

Example static data for running a 1-day simulation can be found under [HICAR-model/Test-data](https://github.com/HICAR-Model/Test-Data)

#### Forcing data requirements

HICAR requires at least the following 5 fields to run, all of which are contained within one netCDF file:
- U winds
- V winds
- Pressure
- Temperature (potential or normal temperature)
- Humidity (mixing ratio or specific humidity)

When using the variational wind solver, providing W winds from forcing data can lead to better estimates of wind speeds.

To perform nested runs, HICAR can be forced with the output from a previous HICAR simulation. Thus HICAR also supports the forcing of all hydrometeors and all of their moments as according to the microphysics scheme chosen. It is recommended to specify these forcing variables when forcing HICAR with output from coarser resolution HICAR runs.

HICAR reads in forcing data from a forcing file list supplied to the model in the namelist. A shell script for generating a forcing file list from a given directory is found within helpers/filelist_script.sh

Example forcing data for running a 1-day simulation can be found under: [HICAR-model/Test-data](https://github.com/HICAR-Model/Test-Data)

#### Namelist

Example namelists for a range of configurations are provided under `helpers/example_namelists/`. The complete set of namelist options, with default values and inline comments describing their function, can be generated from the compiled executable with `./bin/HICAR --gen-nml my_options.nml` (see [docs/namelist_options.md](docs/namelist_options.md)).

#### Supplementary data

Supplementary data for running HICAR, for example look-up tables needed for the ISHMAEL microphysics scheme, can be found in the following repo: https://github.com/NCAR/icar_supporting_files


#### Generating Static Data

HICAR relies on pre-computed static data to speed up some of it’s online calculations. To generate a HICAR domain file, an existing netCDF file with lat, lon, and a DEM is required; landuse categories and a land mask are optional (the land mask is derived from landuse when present, and the domain otherwise defaults to all-land). The lat and lon variables must be named **lat** and **lon**, and the terrain variable must be named **topo**. Additionally, a larger extent DEM of the same resolution is needed to generate parameters for terrain-shading of radiation. I.e., if you have a 50m resolution domain, a larger DEM with an extent ~20km beyond the boundaries of the target domain is also needed.

For information on how to generate the rest of the variables used by HICAR, namely those for calculating terrain shading, see docs/domain_generation.md

#### Reference

Reynolds, D. S., Gutmann, E., Kruyt, B., Haugeneder, M., Jonas, T., Gerber, F., Lehning, M., and Mott, R.: The High-resolution Intermediate Complexity Atmospheric Research (HICAR v1.0) Model Enables Fast Dynamic Downscaling to the Hectometer Scale, Geosci. Model Dev. Discuss. [preprint], https://doi.org/10.5194/gmd-2023-16, in review, 2023. 

Gutmann, E. D., I. Barstad, M. P. Clark, J. R. Arnold, and R. M. Rasmussen (2016), *The Intermediate Complexity Atmospheric Research Model*, J. Hydrometeor, doi:[10.1175/JHM-D-15-0155.1](http://dx.doi.org/10.1175/JHM-D-15-0155.1).

Steger, C. R., Steger, B. and Schär, C. (2022): HORAYZON v1.2: an efficient and flexible ray-tracing algorithm to compute horizon and sky view factor, Geosci. Model Dev., 15, 6817–6840, https://doi.org/10.5194/gmd-15-6817-2022
