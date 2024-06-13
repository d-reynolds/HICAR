
# Dependencies

Before compiling HICAR, a few dependencies must first be installed on your system. **If you are using an HPC which supports modules, first check to see if the relevant dependencies are available on your system as modules (See section "Example (HPC systems)" below for details).** If certain dependencies are not available as modules, then proceed by installing  as described below.

HICAR relies upon four external libraries:

- MPI
- FFTW
- PETSc
- Parallel NetCDF

For an example of how to install these dependencies on a Linux system, see the script HICAR/.github/scripts/hicar_install_utils.sh

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

## PETSc

HICAR has been tested using both packaged PETSc distributions, and by building PETSc from source. Packages can be installed on a linux system with:

```bash
sudo apt-get install petsc-dev
```

While information on installing PETSc from source (or other packaged distributions) can be found at: <https://petsc.org/release/install/>

When installing PETSc from source, be aware that the default installation option is to install the library in debug mode. Information on how to install optimized PETSc is found under the section <https://petsc.org/release/install/install/#compilers>

After installing PETSc, either from source or as a package, read the information about setting the environment variables `PETSC_DIR` and `PETSC_ARCH` here: <https://petsc.org/release/install/multibuild/>

## NetCDF

HICAR uses NetCDF-fortran to read and write NetCDF files in parallel. If you are not using a module for parallel NetCDF, compiling NetCDF-fortran from source is recommended. NetCDF can be configured to work with classic NetCDF file formats, NetCDF-4 file formats, or both. To work with classic NetCDF files, PnetCDF is required. To work with NetCDF-4 files, HDF5 is required.

On a linux system, the script HICAR/.github/scripts/hicar_install_utils.sh can be used to install parallel NetCDF for classic and NetCDF-4 file formats using the following commands:

```bash
HICAR/.github/scripts/hicar_install_utils.sh install_zlib
HICAR/.github/scripts/hicar_install_utils.sh install_hdf5
HICAR/.github/scripts/hicar_install_utils.sh install_pnetcdf
HICAR/.github/scripts/hicar_install_utils.sh install_netcdf_c
HICAR/.github/scripts/hicar_install_utils.sh install_netcdf_fortran
```

# Compiling HICAR

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

Currently, HICAR has been succesfully compiled and tested using the Cray and GNU fortran compilers. **Intel compiler support is still included as legacy, but its function is not insured.**

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
export NETCDF_DIR=/path/to/netcdf/intall
export FFTW_DIR=/path/to/fftw/intall
export PETSC_DIR=/path/to/petsc/intall
export PATH=/path/to/netcdf/bin:${PATH}                               # This line is needed to find the nc-config command installed with NetCDF, which
                                                                      # is used to determine the correct libraries to link to.
export LD_LIBRARY_PATH=/path/to/netcdf/intall:${LD_LIBRARY_PATH}      # This line is quite important when NetCDF has been manually installed and linked
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
module load cray-petsc                   # Load PETSc
module load cray-tpsl                    # Load scientific libraries needed for PETSc
module load cray-fftw                    # Load FFTW
module load cray-netcdf-hdf5parallel     # Load Parallel NetCDF

cmake ../                                # Generate the makefile
make -j 4                                
make install
```

This will then install the executable HICAR in the directory HICAR/bin/

## Options

The following are options which can be passed to the cmake command using `-DOPTION`. For example, to generate a makefile for compiling HICAR as a debug run with the FSM snow model linked, you could run:

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

    FSM=               # Option to link HICAR to optional FSM code libraries compiled separately
        ON             # (DEFAULT; If no libraries are found, HICAR is not linked to FSM)
        OFF            #

    CAF=               # Flag for using Coarray-Fortran for halo exchanges if compiling using the cray fortran compiler
        OFF            # (DEFAULT)
        ON             #

    ASSERTIONS=        # Check for logical assertions at runtime. Used sparingly, little effect.
        ON             # (DEFAULT)
        OFF            #
```

The following flags are used to help cmake locate installed dependencies. Though not explicitly necesarry, they help cmake find dependencies, especially if they have been installed to unusual locations on the system.

```bash
    PETSC_DIR=         # Path to PETSC installation. If not set, defaults to environment variable, if set

    PETSC_ARCH=        # System architecture for PETSC installation. If not set, defaults to environment variable, if set

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

If the user wants to use the snowmodel [FSM2trans](https://egusphere.copernicus.org/preprints/2023/egusphere-2023-2071/), it must also be compiled prior to compiling HICAR. The process for compiling FSM2trans is simple, and similar to that for HICAR:

```bash
cd HICAR/FSM2trans/build   # Navigate to the FSM2trans folder
mkdir build                # Make the build directory
cd build
cmake ../                  # Generate makefile
make -j 4
make install               # FSM2trans is now installed
```

FSM2trans is now installed and will be automatically linked when compiling HICAR
