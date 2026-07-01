# Common errors

This page details some common errors users may hit, and their known fixes/explanations.

## 1. "Too few MPI ranks" / no compute tasks

HICAR reserves **one rank per node for asynchronous I/O**, so it must be run with
at least two MPI ranks (`-np 2`). Running with `-np 1` leaves zero compute
ranks and the model cannot run. See [Running](running.md) for how the I/O and
compute ranks are divided.

## 2. Output file already exists with a different grid

If output files from a previous run are present but were written on a different
grid (different number of levels, latitudes, etc.) or with a different variable
set, HICAR tries to write into the existing file and stops because the dimensions
or variables don't match. Typical messages:

```text
NetCDF: Start+count exceeds dimension bound
saving:qv
```

The fix is to delete or move the existing output files, or change your options so
the output is consistent with the existing files.

## 3. netCDF file-type vs. parallel-I/O backend mismatch

HICAR reads and writes netCDF in parallel, and the file *format* must match the
backend the netCDF library was built against:

- **PnetCDF** handles only **classic** netCDF formats.
- **Parallel HDF5** handles only the **netCDF-4** format.

If your forcing/domain files are netCDF-4 but the library was built only with
PnetCDF (or vice versa), I/O operations will fail. Make sure that the netCDF library
HICAR was compiled with supports the file type you are using (the generic netCDF
query `nc-config --has-nc4` should report `yes` for netCDF-4 support — note HICAR's
own build instead probes `nc-config --has-parallel` to detect parallel-I/O capability). A symptom of building
without netCDF-4 support is a compile-time error like:

```text
Error: Symbol 'nf90_netcdf4' at (1) has no IMPLICIT type
```

Recompile the netCDF stack with HDF5 (and PnetCDF) support — the
`hicar_install_utils.sh` helper does this; see [Compiling](compiling.md).

## 4. Slow parallel I/O with OpenMPI

When using **OpenMPI**, its default OMPIO component can be slow for parallel HDF5
writes; forcing the older ROMIO component can help. HICAR does not set this for you —
it is external OpenMPI runtime tuning. The ROMIO component name is version-specific,
so discover it first:

```bash
ompi_info | grep -i romio     # e.g. romio321 on OpenMPI <=4.x, romio341 on 5.x
export OMPI_MCA_io=<name>      # use the component name your build reports (i.e. romio321 on OpenMPI <=4.x)
```

This can dramatically speed up output on OpenMPI systems. (It is not needed with
MPICH/Cray MPI.)

## 5. Changing a namelist option has no effect

Two common causes:

- **Namelist group not read.** The `lt`, `mp`, `adv`,
  `sm`, `lsm`, `cu`, `rad`, `pbl`, `sfc`, and `wind` groups are only read when the
  corresponding `use_<x>_options` flag in the `&general` group is set to `.true.` (the default behavior).
  This lets those groups be omitted for shorter namelists. Set the relevant
  `use_*_options = .true.` to have the group read.
- **The option was rejected as invalid.** Run `./HICAR --check-nml your.nml`
  first to validate the namelist; HICAR will report unknown or mis-set options
  without running. This may be higher up in the model output, so be sure to read all of the
  output. See [Namelist options](namelist_options.md) for more information.

## 6. Floating-point errors / NaNs

Everyone's nightmare. See [It Crashed](it_crashed.md) 

## 7. Segmentation fault on startup

This may be the shell's stack-size limit, especially for larger domains. Raise
it before running:

```bash
ulimit -s unlimited
```

