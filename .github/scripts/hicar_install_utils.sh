#!/usr/bin/env bash

# see link for size of runner
# https://docs.github.com/en/actions/using-github-hosted-runners/about-github-hosted-runners#supported-runners-and-hardware-resources

# if no arguments passed to the function, or if the only argument is "-h" or "--help" then list 
# the available functions
if [ $# -eq 0 ] || [ $# -eq 1 -a \( "$1" == "-h" -o "$1" == "--help" \) ]; then
    bold=$(tput bold)
    normal=$(tput sgr0)
    echo "Usage: hicar_install_utils.sh [function1] [function2] ..."
    echo ""
    echo "Available functions:"
    echo "  ${bold}hicar_dependencies${normal}: install dependencies for HICAR"
    echo "  ${bold}hicar_install${normal}: build and install HICAR"
    echo "  ${bold}gen_test_run_data${normal}: build test executable and download test data"
    echo "  ${bold}execute_test_run${normal}: run CLI, unit, and integration tests"
    echo "  ${bold}install_zlib${normal}: install zlib"
    echo "  ${bold}install_hdf5${normal}: install hdf5"
    echo "  ${bold}install_PnetCDF${normal}: install PnetCDF"
    echo "  ${bold}install_netcdf_c${normal}: install netcdf-c"
    echo "  ${bold}install_netcdf_fortran${normal}: install netcdf-fortran"
    exit 0
fi

# check if the GITHUB_WORKSPACE variable is set
if [ -z "$GITHUB_WORKSPACE" ]; then
    export HUMAN_RUN=1
    echo "------------------------------------------"
    echo "        HICAR installation script"
    echo "------------------------------------------"
    # this is not a github action, so prompt the user for the location of the HICAR directory 
    # and the install and work directories
    GITHUB_WORKSPACE=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
    # now move two directories up to get to the root of the HICAR directory
    GITHUB_WORKSPACE=$(dirname $(dirname $GITHUB_WORKSPACE))
    echo "Installing HICAR found in repository: $GITHUB_WORKSPACE"
    echo ""
    if [ -z "$WORKDIR" ]; then 
        echo "Please enter the location of the working directory where library files will be downloaded (blank for $HOME/workdir)"
        read WORKDIR
        echo "Using working directory: $WORKDIR"
        echo ""
    fi
    if [ -z "$INSTALLDIR" ]; then
        echo "Please enter the location of the install directory where libraries will be installed (blank for $HOME/installdir)"
        read INSTALLDIR
        echo "Using install directory: $INSTALLDIR"
    fi
else
    export HUMAN_RUN=0
fi   

set -e
set -x

export JN=-j8

if [ -z "$WORKDIR" ]; then
    export WORKDIR=$HOME/workdir
    if [ ! -d "$WORKDIR" ]; then
        mkdir -p $WORKDIR
    fi
fi

if [ -z "$INSTALLDIR" ]; then
    export INSTALLDIR=$HOME/installdir
    if [ ! -d "$INSTALLDIR" ]; then
        mkdir -p $INSTALLDIR
    fi
fi
export LD_LIBRARY_PATH=${INSTALLDIR}/lib:${LD_LIBRARY_PATH}

# Under GitHub Actions, each step runs in a fresh shell, so an `export` here is
# lost by the next step. Persist the runtime library path to $GITHUB_ENV so the
# later steps that RUN the HICAR executable (CLI test, unit tests, integration)
# can resolve the dependency libs (HDF5/NetCDF/FFTW) in $INSTALLDIR/lib. Guarded
# so the path is added once, not accumulated across the deps/build/data steps.
if [ -n "$GITHUB_ENV" ] && [[ ":${PERSISTED_LD_LIBRARY_PATH:-}:" != *":${INSTALLDIR}/lib:"* ]]; then
    echo "LD_LIBRARY_PATH=${LD_LIBRARY_PATH}" >> "$GITHUB_ENV"   # already has ${INSTALLDIR}/lib prefixed above
    echo "PERSISTED_LD_LIBRARY_PATH=${INSTALLDIR}/lib" >> "$GITHUB_ENV"
fi

function install_zlib {
    echo install_zlib
    cd $WORKDIR

    if [ ! -d "$WORKDIR/zlib-1.3.1" ]; then
        # zlib.net only hosts the current release at its top-level URL (older
        # versions move to /fossils and the old URL 404s). Use the permanent
        # GitHub release asset instead.
        wget --no-check-certificate -q https://github.com/madler/zlib/releases/download/v1.3.1/zlib-1.3.1.tar.gz
        tar -xvzf zlib-1.3.1.tar.gz
    fi

    cd zlib-1.3.1/

    # Always run configure for zlib (its tarball ships with a Makefile
    # that lacks install targets until configure sets the prefix)
    ./configure --prefix=$INSTALLDIR &> config.log

    make -j 8 &> make.log
    make install
    make check
}

function install_hdf5 {
    echo install_hdf5
    cd $WORKDIR

    export CPPFLAGS=-I$INSTALLDIR/include
    export LDFLAGS=-L$INSTALLDIR/lib
    export CC=${CC:-mpicc}

    if [ ! -d "$WORKDIR/hdf5-1.14.3" ]; then
        wget --no-check-certificate -q https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.14/hdf5-1.14.3/src/hdf5-1.14.3.tar.gz
        tar -xzf hdf5-1.14.3.tar.gz
    fi
    
    cd hdf5-1.14.3

    # ALWAYS run configure. The HDF5 release tarball ships a top-level Makefile
    # whose default prefix is <builddir>/hdf5, so a `[ ! -f Makefile ]` guard would
    # skip configure and silently install to the wrong place (downstream netcdf-c
    # then fails with "cannot find -lhdf5"). Re-running configure forces our prefix.
    ./configure --prefix=$INSTALLDIR --enable-parallel --with-zlib=$INSTALLDIR

    make -j 8
    make install

    # Fail loudly if the library did not land where downstream builds expect it.
    if ! ls "$INSTALLDIR"/lib/libhdf5.* >/dev/null 2>&1; then
        echo "ERROR: libhdf5 not found in $INSTALLDIR/lib after install" >&2
        exit 1
    fi

    export HDF5=$INSTALLDIR
    export HDF5_DIR=$INSTALLDIR
    export LD_LIBRARY_PATH=$INSTALLDIR/lib:$LD_LIBRARY_PATH
}

function install_PnetCDF {
    echo install_PnetCDF
    cd $WORKDIR

    export CPPFLAGS=-I$INSTALLDIR/include 
    export LDFLAGS=-L$INSTALLDIR/lib

    if [ ! -d "$WORKDIR/pnetcdf-1.13.0" ]; then
        wget --no-check-certificate -q https://parallel-netcdf.github.io/Release/pnetcdf-1.13.0.tar.gz
        tar -xzf pnetcdf-1.13.0.tar.gz
    fi

    cd pnetcdf-1.13.0

    # check if make file exists and if not, run configure
    if [ ! -f "Makefile" ]; then
        ./configure --prefix=${INSTALLDIR}
    fi

    make -j 8
    make install
}

function install_netcdf_c {
    echo install_netcdf_c
    cd $WORKDIR

    export CPPFLAGS=-I$INSTALLDIR/include
    export LDFLAGS=-L$INSTALLDIR/lib
    export CC=${CC:-mpicc}
    export LIBS=-ldl

    if [ ! -d "$WORKDIR/netcdf-c-4.9.2" ]; then
        wget --no-check-certificate -q https://github.com/Unidata/netcdf-c/archive/refs/tags/v4.9.2.tar.gz
        tar -xzf v4.9.2.tar.gz

    fi
    
    cd netcdf-c-4.9.2

    rm -f "$INSTALLDIR"/lib/*.la

    # check if make file exists and if not, run configure
    if [ ! -f "Makefile" ]; then
        ./configure --prefix=${INSTALLDIR} --disable-shared --enable-pnetcdf --enable-parallel-tests
    fi

    make -j 8
    make install
}

function install_netcdf_fortran {
    cd $WORKDIR

    export CPPFLAGS=-I$INSTALLDIR/include 
    export LDFLAGS=-L$INSTALLDIR/lib
    export PATH=$INSTALLDIR/bin:$PATH

    export LD_LIBRARY_PATH=${INSTALLDIR}/lib:${LD_LIBRARY_PATH}
    export LIBS=$(nc-config --libs)

    if [ ! -d "$WORKDIR/netcdf-fortran-4.6.1" ]; then
        wget --no-check-certificate -q https://github.com/Unidata/netcdf-fortran/archive/refs/tags/v4.6.1.tar.gz

        tar -xzf v4.6.1.tar.gz
    fi

    cd netcdf-fortran-4.6.1

    rm -f "$INSTALLDIR"/lib/*.la

    # check if make file exists and if not, run configure
    if [ ! -f "Makefile" ]; then
        CC=${CC:-mpicc} FC=${FC:-mpif90} F77=${F77:-mpif77} ./configure --prefix=${INSTALLDIR} --disable-shared
    fi

    make -j 8
    make install
}



function hicar_dependencies {
    echo hicar_dependencies
    sudo apt-get update
    sudo apt-get install mpich
    sudo apt-get install libcurl4-gnutls-dev
    sudo apt-get install libfftw3-dev

    install_zlib
    install_hdf5
    install_PnetCDF
    install_netcdf_c
    install_netcdf_fortran

    # put installed bin directory in PATH
    export PATH=${INSTALLDIR}/bin:$PATH
}

function hicar_install {
    echo hicar_install
    pwd
    cd ${GITHUB_WORKSPACE}
    if [ ! -d "$GITHUB_WORKSPACE/build" ]; then
        mkdir build
        cd build
    else
        cd build
        rm -rf *
    fi
    export NETCDF_DIR=${INSTALLDIR}
    export FFTW_DIR=/usr
    export PATH=${INSTALLDIR}/bin:$PATH
    export LD_LIBRARY_PATH=${INSTALLDIR}/lib:${LD_LIBRARY_PATH}
    cmake ../ -DFSM=OFF -DMODE=${HICAR_MODE:-debug} ${HICAR_CMAKE_EXTRA:-}
    make ${JN}
    make install
    
    echo "hicar install succeeded"

}

function gen_test_run_data {
    echo gen_test_run_data
    cd ${GITHUB_WORKSPACE}/build
    make download_test_data
    make HICAR-tester
}

function execute_test_run {
    echo execute_test_run
    cd ${GITHUB_WORKSPACE}

    # 1. CLI options test (fast, no MPI needed)
    echo "--- Running CLI options tests ---"
    bash tests/test_cli_options.sh ${GITHUB_WORKSPACE}

    # 2. Unit tests (MPI)
    echo "--- Running unit tests ---"
    cd build
    mpiexec -np 4 tests/HICAR-tester
    cd ..

    # 3. Integration tests (MPI + test data)
    echo "--- Running integration tests ---"
    cd tests/Test_Cases
    bash test_case_runner.sh ${GITHUB_WORKSPACE} Standard
    cd ../..
}

for func in "$@"
do
    case $func in
        hicar_dependencies)
            hicar_dependencies;;
        hicar_install)
            hicar_install;;
        gen_test_run_data)
            gen_test_run_data;;
        execute_test_run)
            execute_test_run;;
        install_zlib)
            install_zlib;;
        install_hdf5)
            install_hdf5;;
        install_PnetCDF)
            install_PnetCDF;;
        install_netcdf_c)
            install_netcdf_c;;
        install_netcdf_fortran)
            install_netcdf_fortran;;
        *)
            echo "$func unknown";;
    esac
done
