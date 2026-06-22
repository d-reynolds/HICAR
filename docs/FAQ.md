# FAQ

## Does HICAR run on GPUs?

Yes. HICAR is GPU-accelerated with **OpenACC**. Build it with `-DOPENACC=ON`
using the NVHPC / NVFortran compiler (see [Compiling](compiling.md)); the same
source also builds and runs on CPUs with GNU Fortran. On GPU systems each compute
MPI rank drives one GPU, and halo exchanges can use NCCL when it is available.

## Which compilers are supported?

**GNU Fortran** (CPU) and **NVFortran / NVHPC** (CPU host and GPU). Intel and
Cray support is deprecated and no longer tested.

## Which snow model should I use?

HICAR always builds with **SNOWPACK**. By default it uses the native-Fortran port
(`snowpack_driver.F90`); passing `-DSNOWPACK_CPP=ON` selects the C++
Alpine3D/SNOWPACK wrapper instead. **FSM2** is available as an optional model that
must be compiled separately and linked with `-DFSM=ON` (see
[Compiling](compiling.md)).

## Does HICAR support nested domains?

Yes — **one-way nesting**. Each domain has exactly one parent, and the forcing
data must encompass all domains. Nesting is configured with the `nests` and
`parent_nest` options, and per-domain settings are given as ordered lists. See
[Namelist options](namelist_options.md#nested-runs) for the details.

## What forcing data does HICAR need?

At a minimum: **U and V winds, pressure, temperature, and humidity**, all in a
single netCDF file. Supplying **W winds** improves the iterative variational wind
solver. When forcing HICAR from the output of a coarser HICAR run, it is
recommended to also provide the hydrometeor fields. See
[Forcing data](forcing_data.md).

## Does HICAR include cumulus, PBL, land-surface, and radiation schemes?

Yes. HICAR carries the same WRF-derived physics as ICAR, generally in their
**full** form (the *dynamics* are what is simplified, not the physics). Available
options include cumulus (Tiedtke, Kain-Fritsch, BMJ, NSAS), the full YSU PBL
scheme, the Noah-MP land-surface model, RRTMG/RRTMGP radiation, and a range of
microphysics schemes (Thompson, Morrison, WSM, ISHMAEL, simple). The active
schemes are selected in the `&physics` namelist.

## When linear winds are used, is the look-up table for dry or moist N²?

This only applies when running the optional linear-theory wind solver (the
iterative variational solver is the main solver). The look-up table itself knows
nothing about dry vs. moist: the atmosphere responds in roughly the same way for a
given N² either way. At run time the model computes the dry or moist N² from the
current atmospheric conditions and uses that to look up the correct perturbation.
