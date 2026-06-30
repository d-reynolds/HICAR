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

#### Reference

Reynolds, D. S., Gutmann, E., Kruyt, B., Haugeneder, M., Jonas, T., Gerber, F., Lehning, M., and Mott, R.: The High-resolution Intermediate Complexity Atmospheric Research (HICAR v1.0) Model Enables Fast Dynamic Downscaling to the Hectometer Scale, Geosci. Model Dev. Discuss. [preprint], https://doi.org/10.5194/gmd-2023-16, in review, 2023. 

Gutmann, E. D., I. Barstad, M. P. Clark, J. R. Arnold, and R. M. Rasmussen (2016), *The Intermediate Complexity Atmospheric Research Model*, J. Hydrometeor, doi:[10.1175/JHM-D-15-0155.1](http://dx.doi.org/10.1175/JHM-D-15-0155.1).
