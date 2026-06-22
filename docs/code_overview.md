# Code overview

This page is a map of the HICAR source tree to help you find the right file.
It is a high-level orientation, not an API reference — the authoritative
description of any routine is the code itself. For how to *add* to the code
(new variables, new physics), see [Developing](developing.md).

All Fortran source lives under `src/`, grouped into `main/`, `objects/`,
`physics/`, `io/`, `utilities/`, and `constants/`.

## A note on the `_h` / `_obj` split

Most of HICAR's derived types are split across two files using **Fortran
submodules**: a `*_h.F90` file declares the `module` — the type definition and
the interfaces of its procedures — and a `*_obj.F90` file is the
`submodule(...)` that contains the implementations. For example
`objects/domain_h.F90` declares `module domain_interface` and the `domain_t`
type, while `objects/domain_obj.F90` is `submodule(domain_interface)
domain_implementation` and implements its methods. When looking for *what* a
type can do, read the `_h`; for *how* it does it, read the `_obj`.

## `main/` — driver and time integration

- **`driver.F90`** — the main program. It initializes MPI, reads the namelist
  (`init_options`), splits the MPI ranks into compute and I/O groups
  (`split_processes`), initializes each nest (`component_init`), then runs the
  main loop (`component_loop`) and shuts down (`component_program_end`).
- **`flow_events.F90`** — orchestrates the per-nest event loop: initializing,
  waking nested domains, reading new forcing, stepping, and writing output
  (`component_init`, `component_loop`, `component_read`, `component_write`,
  `component_program_end`). This is the high-level control flow.
- **`nest_manager.F90`** — manages one-way grid nesting: which nest runs next,
  when child nests should be updated from their parent, and switching the active
  nest context (`switch_nest_context`). Note that multiple nests run in one
  process, so physics-module init routines are re-invoked on each context
  switch.
- **`time_step.F90`** — `step()` integrates a single nest forward over one
  forcing interval, sub-stepping under the advection CFL limit and calling the
  physics drivers (advection, microphysics, PBL, radiation, land surface,
  cumulus) in sequence.
- **`init.F90`** — model and option initialization (`init_model`,
  `init_options`, `split_processes`).

## `objects/` — the model state

The model state and infrastructure are derived types, each as an `_h`/`_obj`
pair:

- **`domain`** — the central object: the model variables, the grid
  decomposition, the halos, and the methods to allocate, exchange, and update
  them.
- **`boundary`** — the forcing/boundary-condition state read from the forcing
  files.
- **`halo`** — halo (ghost-cell) exchange between MPI ranks. The backend is MPI,
  or NCCL on GPU builds when available.
- **`grid`** — grid geometry and the domain decomposition (`ims:ime` memory,
  `its:ite` tile/interior, `ids:ide` global — see Developing).
- **`options`** — the parsed namelist options, including which variables to
  allocate and output.
- **`variable` / `variable_dict`** — the `variable_t` container (1-D…4-D arrays
  plus metadata) and a dictionary of them keyed by the `kVARS` indices.
- **`meta_data`** — variable metadata (names, units, attributes, valid ranges).
- **`flow`** — the per-nest "component" object stepped by the flow-event loop.

## `constants/` — indices and physical constants

- **`icar_constants.F90`** — defines the `kVARS` structure (every model
  variable has an integer index `kVARS%name`), the storage-size constants, and
  the integer codes that select physics schemes, wind solvers, snow models,
  etc.
- **`wrf_constants.F90`** — physical constants shared with WRF-derived physics.

## `physics/` — parameterizations and dynamics

Physics files follow a **prefix naming convention**, and each category has a
`*_driver.F90` that is the entry point called from `time_step.F90`:

| Prefix / driver | Category |
|---|---|
| `mp_*`, `mp_driver.F90` | Microphysics (Thompson, Morrison, WSM, ISHMAEL, simple) |
| `advect.F90`, `adv_*`, `advection_driver.F90` | Advection (incl. flux-corrected transport) |
| `pbl_*`, `pbl_driver.F90` | Planetary boundary layer (YSU, simple, diagnostic) |
| `ra_*`, `ra_driver.F90` | Radiation (RRTMG, RRTMGP, simple) |
| `lsm_driver.F90`, Noah-MP (`noahmp/`) | Land surface |
| `cu_*`, `cu_driver.F90` | Cumulus (Tiedtke, Kain-Fritsch, BMJ, NSAS) |
| `sfc_*`, `sfc_driver.F90` | Surface-layer scheme |

**Dynamics / wind solvers.** Unlike a full NWP model, HICAR's "dynamics" are the
construction of the high-resolution wind field:

- `linear_winds.F90` — analytical linear-theory perturbation (with a look-up
  table).
- `wind_iterative.F90` — the iterative variational wind solver (mass-conserving
  wind balancing).
- `wind.F90`, `wind_surf.F90`, `wind_thermal.F90` — wind balancing, surface, and
  thermal-wind adjustments.

**Snow / surface-mass models.** `sm_driver.F90` dispatches to the configured
snow model: `snowpack_driver.F90` (the native-Fortran SNOWPACK port, default) or
`sm_SNOWPACK.F90` (the C++ Alpine3D/SNOWPACK wrapper, with `-DSNOWPACK_CPP=ON`),
or `sm_FSM.F90` (FSM2, optional). `snow_drift.F90` and `snowslide.F90` add
blowing-snow drift and gravitational snow transport. `water_simple.F90`,
`water_lake.F90`, and `water_flake.F90` (`flake_core.F90`) are the lake/water
options.

## `io/` — asynchronous parallel I/O

HICAR overlaps reading and writing with computation by dedicating some MPI ranks
to I/O (see [Running](running.md)):

- **`ioserver` (`ioserver_h/obj`)** — the dedicated I/O ranks that gather data
  and read/write the netCDF files.
- **`ioclient` (`ioclient_h/obj`)** — the compute-rank side that hands data to,
  and receives data from, the I/O servers.
- **`output` (`output_h/obj`)**, **`reader` (`reader_h/obj`)** — output and
  forcing-file readers.
- **`default_output_metadata.F90`** — the `get_varmeta` function: the single
  source of truth mapping each `kVARS` index to its output name, units,
  dimensions, and valid range.
- **`io_routines.F90`**, **`lt_lut_io.F90`** (linear-theory look-up table),
  **`time_io.F90`** — supporting read/write helpers.

## `utilities/`

General-purpose helpers: calendar/time handling (`time_obj.F90`,
`time_delta_obj.F90`, `time_io.F90`), geographic/horizontal interpolation
(`geo_reader.F90`), vertical interpolation (`vinterp.F90`), atmospheric helper
functions (`atm_utilities.F90`), string handling (`string.F90`), FFTs
(`fftshift.F90`, `fftw.F90`, `cufft_interface.F90`), MPI helpers
(`mpi_utils.F90`), the NCCL interface (`nccl_interface.F90`, `nccl_wrapper.c`),
namelist parsing (`namelist_utilities.F90`), array helpers
(`array_utilities.F90`), assertions (`assertions.F90`), and timers
(`timer_obj.F90`).

## Tests

Component tests live under `tests/` as [test-drive] suites: `test_driver.F90`
dispatches to the individual `test_<name>.F90` files. The runnable suite names are
`advection`, `snow_drift`, `control_flow`, `halo_exch`, `geo`, `time`, and
`utilities` (note the advection suite is `advection`, though its file is
`test_advect.F90`). Run them all with `make check`, or a single suite with
`mpiexec -np 2 tests/HICAR-tester <suite>`.
See [Testing](testing.md) for local runs and
[CI/CD pipeline](ci_cd_pipeline.md) for how they fit into the automated gates.

[test-drive]: https://github.com/fortran-lang/test-drive
