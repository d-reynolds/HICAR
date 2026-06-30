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

HICAR always builds with [SNOWPACK](https://code.wsl.ch/snow-models/snowpack), 
and this configuration is recommended. By default it uses a native-Fortran port
(`snowpack_driver.F90`); passing `-DSNOWPACK_CPP=ON` selects the C++
SNOWPACK wrapper instead. **FSM2** is available as an optional model that
must be compiled separately and linked with `-DFSM=ON` (see
[Compiling](compiling.md)). Unfortunately, the FSM2 model that HICAR has been
developed with is closed source, so contact the developers to get a working repo.

## Does HICAR support nested domains?

Yes — **one-way nesting**. Each domain has exactly one parent, and the forcing
data must encompass all domains. Nesting is configured with the `nests` and
`parent_nest` options, and per-domain settings are given as ordered lists. See
[Namelist options](namelist_options.md#nested-runs) for the details. HICAR supports
up to 10 domains, but your RAM will probably run out before that point.

## What forcing data does HICAR need?

At a minimum: **U and V winds, pressure, temperature, and humidity**, all in a
single netCDF file, all 3D. Supplying **W winds** improves the iterative variational wind
solver. See [Forcing data](forcing_data.md).

## What static data does HICAR need?

At a minimum: elevation, lat, and lon data on a 2D grid. Landuse is recommended, especially
if using any of the land surface options. To derive additional static data, such as
sky view fraction and other variables needed for terrain-shaded radiation, see [Domain Generation](domain_generation.md).