
# Dependencies

Before compiling HICAR, a few dependencies must first be installed on your system. **If you are using an HPC which supports modules, first check to see if the relevant dependencies are available on your system as modules (See section "Example (HPC systems)" below for details).** If certain dependencies are not available as modules, then proceed by installing the relevant libraries as described below.

HICAR relies upon three external libraries:

- MPI
- FFTW
- Parallel NetCDF

GPU builds additionally use **NCCL** (for halo exchange) and **cuFFT**, which
ship with the NVHPC toolkit; neither needs to be installed separately.

The easiest way to install all dependencies on a Linux system is to call:

```bash
HICAR/.github/scripts/hicar_install_utils.sh hicar_dependencies
```

and follow the prompts. The following subsections detail commands to install each of the individual depencencies.

## MPI

HICAR has been tested using packed MPICH distributions. This can be installed on a linux system with:

```bash
sudo apt-get install mpich
```

## FFTW

HICAR has been tested using packed FFTW3 distributions. This can be installed on a linux system with:

```bash
sudo apt-get install libfftw3-dev
```

## NetCDF

HICAR uses NetCDF-fortran to read and write NetCDF files in parallel. If you are not using a module for parallel NetCDF, compiling NetCDF-fortran from source is recommended. NetCDF can be configured to work with classic NetCDF file formats, NetCDF-4 file formats, or both. To work with classic NetCDF files, PnetCDF is required. To work with NetCDF-4 files, HDF5 is required.

On a linux system, the script HICAR/.github/scripts/hicar_install_utils.sh can be used to install parallel NetCDF for classic and NetCDF-4 file formats using the following commands:

```bash
HICAR/.github/scripts/hicar_install_utils.sh install_zlib
HICAR/.github/scripts/hicar_install_utils.sh install_hdf5
HICAR/.github/scripts/hicar_install_utils.sh install_PnetCDF
HICAR/.github/scripts/hicar_install_utils.sh install_netcdf_c
HICAR/.github/scripts/hicar_install_utils.sh install_netcdf_fortran
```

# Compiling HICAR

The easiest way to compile HICAR is using the install utility command:

```bash
HICAR/.github/scripts/hicar_install_utils.sh hicar_install
```

Instructions for a more manual install of the model follows:

HICAR is compiled using cmake to first generate a makefile. The full steps to generate a makefile are as follows

```bash
cd ~/
git clone https://github.com/HICAR-Model/HICAR.git
cd HICAR
mkdir build
cd build
cmake ../
```

The cmake file will attempt to find any fortran comiplers already set in your environment. If the correct one is not found, it can be set with

```bash
cmake ../ -DFC=gfortran
```

HICAR is compiled and tested with the **GNU Fortran** (CPU) and **NVFortran /
NVHPC** (CPU host and GPU) compilers. NVFortran/NVHPC is required for GPU builds
(`-DOPENACC=ON`). **Intel and Cray compiler support is deprecated and no longer
tested.**

The generated makefile can then be run with the standard

```bash
make -j 4
make install
```

This will then install the executable HICAR in the directory HICAR/bin/

For cmake options available when generating the makefile, see the section "Options" below.

## Example

This example shows how to compile the model after installing dependencies, assuming a unix environment.

Starting from the root reposititory (HICAR/):

```bash
mkdir build
cd build
export NETCDF_DIR=/path/to/netcdf/root
export FFTW_DIR=/path/to/fftw/root
export PATH=/path/to/netcdf/bin:${PATH}                               # This line is needed to find the nc-config command installed with NetCDF, which
                                                                      # is used to determine the correct libraries to link to.
export LD_LIBRARY_PATH=/path/to/netcdf/root:${LD_LIBRARY_PATH}      # This line is quite important when NetCDF has been manually installed and linked
cmake ../
make -j 4
make install
```

## Example (HPC systems)

This example shows how to compile the model on an HPC running linux with the dependencies added as modules.
These modules automatically set the search paths necesarry to find the packages.

Starting from the root reposititory (HICAR/):

```bash
mkdir build
cd build

# exact module names will vary, these are relevant for the CSCS HPC Daint
module load daint-mc                     
module load CMake
module load fftw                    # Load FFTW
module load netcdf-hdf5parallel     # Load Parallel NetCDF

cmake ../                                # Generate the makefile
make -j 4                                
make install
```

This will then install the executable HICAR in the directory HICAR/bin/

## Example (GPU build)

To build HICAR with GPU acceleration, use the NVHPC / NVFortran compiler and pass
`-DOPENACC=ON`. NCCL and cuFFT (shipped with NVHPC) are picked up automatically
when available.

```bash
mkdir build
cd build

# Load the NVHPC toolchain and an NVFortran-built NetCDF/HDF5/PnetCDF stack.
# Exact module names vary by system.
module load nvhpc
module load netcdf

cmake ../ -DFC=nvfortran -DOPENACC=ON
make -j 4
make install
```

The GPU build maps one GPU to each compute MPI rank at run time (see
[Running](running.md)). The relocatable device-code link is the dominant cost of
a GPU build and takes several minutes regardless of optimization level.

## Options

The following are options which can be passed to the cmake command using `-DOPTION`. For example, to generate a makefile for compiling HICAR as a debug run without the FSM snow model linked, you could run:

```bash
cmake ../ -DFSM=OFF -DMODE=debug
```

Full list of user options, not including standard cmake options, are:

```bash
    MODE=              # Set compilation flags, useful for debugging
        release        # full optimization, slower compile (DEFAULT)
        debug          # debug compile with optimizations
        debugslow      # debug compile w/o optimizations
        profile        # set profiling options for gnu or intel compilers

    FC=                # Set the fortran compiler to use (can be auto-detected for most cases)

    OPENACC=           # Enable OpenACC GPU acceleration (requires the NVHPC/NVFortran compiler)
        OFF            # (DEFAULT) CPU build
        ON             # GPU build

    NCCL=              # Link to NCCL for GPU halo exchange, if available
        ON             # (DEFAULT) used only when OPENACC=ON and NCCL is found
        OFF            #

    FSM=               # Option to link HICAR to optional FSM code libraries compiled separately
        OFF            # (DEFAULT)
        ON             # If no libraries are found, HICAR is not linked to FSM

    SNOWPACK_CPP=      # Snow-model implementation. SNOWPACK is always compiled in;
                       #   this flag only selects which implementation is used.
        OFF            # (DEFAULT) native-Fortran SNOWPACK port (snowpack_driver.F90);
                       #   fetches the SNOWPACK fortran-bindings repo at configure time
        ON             # C++ Alpine3D/SNOWPACK wrapper (sm_SNOWPACK.F90); builds C++ SNOWPACK + MeteoIO

    ASSERTIONS=        # Check for logical assertions at runtime. Used sparingly, little effect.
        ON             # (DEFAULT)
        OFF            #

    SRUN_FLAGS=        # srun flags for running the automated test cases on SLURM systems.
                       # See the Testing page for details.
```

> The `SRUN_FLAGS` option only affects how the bundled test cases are launched on
> SLURM systems; see [Testing](testing.md#test-cases-on-slurm).

> Note: SNOWPACK (native-Fortran port) is now always compiled into HICAR; there
> is no flag to build without it. `-DSNOWPACK_CPP=ON` selects the C++ wrapper
> instead. `-DSNOWPACK_FORTRAN=ON` is deprecated (the port is the default) and
> `-DSNOWPACK=OFF` is no longer supported. Re-configuring an existing build
> directory keeps previously cached values; use a fresh build directory.

The following flags are used to help cmake locate installed dependencies. Though not explicitly necesarry, they help cmake find dependencies, especially if they have been installed to unusual locations on the system.

```bash
    FFTW_DIR=          # Path to FFTW installation. If not set, defaults to environment variable, if set

    NETCDF_DIR=        # Path to NETCDF installation. If not set, defaults to environment variable, if set

    MPI_DIR=           # Path to MPI installation. If not set, defaults to environment variable, if set

    FSM_DIR=           # Path to the compiled FSM installation (root directory of the lib and build directories)
                       # Defaults to HICAR/FSM2trans/, assuming that the cmake routine for FSM2trans has been run without modification
```

 If using modules, these directories (besides `FSM_DIR`) should be set automatically. Still, they can be checked for by running

```bash
module show MODULE_NAME
```

# Compiling FSM

If the user wants to use the snowmodel [FSM2](https://egusphere.copernicus.org/preprints/2023/egusphere-2023-2071/), it must also be compiled prior to compiling HICAR. The process for compiling FSM2 is simple, and similar to that for HICAR. 

To couple with FSM2, the user must first get the FSM2 distribution. To get just the source code (recommended), setup a sparse checkout as follows:

```bash
git clone --no-checkout --filter=blob:none https://github.com/oshd-slf/jim_operational.git FSM2
cd FSM2
git sparse-checkout set --cone
git checkout 64221dc358bc22026d08b61fc45f282f12089b9f
git sparse-checkout set FSM_SOURCE_CODE/code/
```

Note that FSM2 development is not gaurenteed to always dovetail nicely with the HICAR interface. The last tested FSM2 commit hash was:

```bash
64221dc358bc22026d08b61fc45f282f12089b9f
```

To compile FSM2 as a library which can be called from HICAR, the steps below can be followed. In the below block, `FSM2_Dir` refers to the root directory of the FSM2 git repo from the previous step.

```bash
cp HICAR/src/physics/FSM2_interface/FSM2_CMakeLists.txt FSM2_Dir/FSM_SOURCE_CODE/CMakeLists.txt  # Copy the FSM2 CMake file to the FSM2 distribution
cd FSM2_Dir/FSM_SOURCE_CODE
mkdir build                # Make the build directory
cd build
cmake ../                  # Generate FSM2 makefile
make -j 4
make install               # FSM2 library is now installed to ../lib
```

FSM2 is now installed. To link FSM2 to HICAR when installing HICAR, set the `FSM_DIR` variable when calling cmake for HICAR. In this way, the previous cmake command under the section "Compiling HICAR" can be changed to:

```bash
cmake ../ -DFC=gfortran -DFSM_DIR=FSM2_Dir/FSM_SOURCE_CODE
```
