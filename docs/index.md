# HICAR — The High-resolution Intermediate Complexity Atmospheric Research Model

HICAR is an atmospheric model for **dynamic downscaling at the hectometer
(sub-kilometer) scale**. It takes the output of a coarser atmospheric model — a
kilometer-scale NWP forecast, a reanalysis, or a climate model — and produces a
high-resolution atmospheric state suitable for driving land-surface, snow, and
hydrological simulations over complex (e.g. alpine) terrain.

HICAR is a variant of the [Intermediate Complexity Atmospheric Research (ICAR)
model](https://github.com/NCAR/icar). Like ICAR, it pairs physics
parameterizations shared with full weather models such as WRF (microphysics,
radiation, PBL, land surface, cumulus) with **massively simplified dynamics**.
Instead of solving the full non-hydrostatic equations of motion, HICAR computes a
high-resolution wind field — either from analytical linear theory or from an
iterative variational wind solver — and uses it to advect heat, moisture, and
hydrometeors around the domain. This makes HICAR roughly two to three orders of
magnitude faster than a full NWP model at comparable resolution.

What distinguishes HICAR from ICAR is the move to **sub-kilometer resolution and
HPC-scale parallelism**. HICAR is parallelized with **MPI** (ICAR v2 used Fortran
coarrays) and is additionally **GPU-accelerated with OpenACC**, so it can run on
both CPU clusters and GPU nodes. It adds one-way grid nesting, asynchronous
parallel I/O, terrain-shaded radiation, and a choice of snow models (SNOWPACK,
FSM2, and simpler schemes).

## Requirements

To run HICAR you need: time-varying 3-D atmospheric forcing (U/V winds,
pressure, temperature, and humidity — at minimum — in a single netCDF file, from
WRF output, reanalysis, a GCM, or a coarser HICAR run), a high-resolution static
**domain file** describing the terrain and land surface, and a namelist
specifying the model options. See [Forcing data](forcing_data.md) and
[Domain generation](domain_generation.md) for how to prepare these inputs.

## Getting started

- **[Tutorial](tutorial.md)** — set up a working directory and run the provided
  test case end to end.
- **[Compiling](compiling.md)** — dependencies, CMake build, and the GPU build.
- **[Running](running.md)** — invoking the model, MPI rank layout, and a Slurm
  template.
- **[Namelist options](namelist_options.md)** — the self-documenting executable
  and nested-run configuration.

For the code, see the **[Code overview](code_overview.md)** and
**[Developing](developing.md)** pages.

## References

If you use HICAR, please cite:

> Reynolds, D. S., Gutmann, E., Kruyt, B., Haugeneder, M., Jonas, T., Gerber, F.,
> Lehning, M., and Mott, R.: *The High-resolution Intermediate Complexity
> Atmospheric Research (HICAR v1.0) Model Enables Fast Dynamic Downscaling to the
> Hectometer Scale*, Geosci. Model Dev., 2023.
> <https://doi.org/10.5194/gmd-2023-16>

<!-- Maintainer note: confirm the final GMD volume/page citation (vs. the
     gmd-2023-16 discussion DOI) before a documentation release. -->

The underlying ICAR model is described in:

> Gutmann, E. D., I. Barstad, M. P. Clark, J. R. Arnold, and R. M. Rasmussen
> (2016): *The Intermediate Complexity Atmospheric Research Model*, J.
> Hydrometeor, 18, 957–973.
> <https://doi.org/10.1175/JHM-D-15-0155.1>
