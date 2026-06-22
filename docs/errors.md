# Common errors

This page collects the failure modes users hit most often, with their fixes.

## 1. Segmentation fault on startup

Often this is the shell's stack-size limit, especially for larger domains. Raise
it before running:

```text
bash:  ulimit -s unlimited
csh:   unlimit stacksize
```

Some systems cap the stack at 64 MB (`ulimit -s 65532`); that may or may not be
enough depending on the domain size. If raising the stack does not help, rebuild
HICAR in debug mode (`-DMODE=debug`, which adds bounds checking and
`-finit-real=nan`) to localize the fault.

## 2. "Too few MPI ranks" / no compute tasks

HICAR reserves **one rank per node for asynchronous I/O**, so it must be run with
**at least two MPI ranks** (`-np 2`). Running with `-np 1` leaves zero compute
ranks and the model cannot run. See [Running](running.md) for how the I/O and
compute ranks are divided.

## 3. Output file already exists with a different grid

If output files from a previous run are present but were written on a different
grid (different number of levels, latitudes, etc.) or with a different variable
set, HICAR tries to write into the existing file and stops because the dimensions
or variables don't match. Typical messages:

```text
NetCDF: Variable not found
setup_varids: Searching for variable in existing file: nsq
```

```text
NetCDF: Start+count exceeds dimension bound
output/hicar_1990_01_01_00-00.nc: qv
```

The fix is to delete or move the existing output files, or change your options so
the output is consistent with the existing files.

## 4. netCDF file-type vs. parallel-I/O backend mismatch

HICAR reads and writes netCDF in parallel, and the file *format* must match the
backend the netCDF library was built against:

- **PnetCDF** handles only **classic** netCDF formats.
- **Parallel HDF5** handles only the **netCDF-4** format.

If your forcing/domain files are netCDF-4 but the library was built only with
PnetCDF (or vice versa), reads or writes will fail. Make sure the netCDF library
HICAR was compiled with supports the file type you are using (`nc-config
--has-nc4` should report `yes` for netCDF-4 support). A symptom of building
without netCDF-4 support is a compile-time error like:

```text
Error: Symbol 'nf90_netcdf4' at (1) has no IMPLICIT type
```

Recompile the netCDF stack with HDF5 (and PnetCDF) support — the
`hicar_install_utils.sh` helper does this; see [Compiling](compiling.md).

## 5. Slow parallel I/O with OpenMPI

When using **OpenMPI**, its default OMPIO component is slow for parallel HDF5
writes. Force the older ROMIO component instead:

```bash
export OMPI_MCA_io=romio321
```

This can dramatically speed up output on OpenMPI systems. (It is not needed with
MPICH/Cray MPI.)

## 6. Changing a namelist option has no effect

Two common causes:

- **Namelist group not read.** The `lt`, `mp`, `adv`,
  `sm`, `lsm`, `cu`, `rad`, `pbl`, `sfc`, and `wind` groups are only read when the
  corresponding `use_<x>_options` flag in the `&general` group is set to `.true.` (the default behavior).
  This lets those groups be omitted for shorter namelists. Set the relevant
  `use_*_options = .true.` to have the group read.
- **The option was rejected as invalid.** Run `./HICAR --check-nml your.nml`
  first to validate the namelist; HICAR will report unknown or mis-set options
  without running. See [Namelist options](namelist_options.md).

## 7. Floating-point errors / NaNs

Everyone's nightmare. Check your forcing data and static data
first. Then, try running the `make test_cases` check locally from your build directory to confirm that the executable works properly. A debug build (`-DMODE=debug`) may detect the first NaN
or error and help with debugging. If none of these checks turns up a solution, please contact the developers.
