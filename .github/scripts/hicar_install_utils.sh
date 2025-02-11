#!/usr/bin/env bash

# see link for size of runner
# https://docs.github.com/en/actions/using-github-hosted-runners/about-github-hosted-runners#supported-runners-and-hardware-resources

# check if the GITHUB_WORKSPACE variable is set
if [ -z "$GITHUB_WORKSPACE" ]; then

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
    echo "Please enter the location of the working directory where library files will be downloaded (blank for $HOME/workdir)"
    read WORKDIR
    if [ -z "$WORKDIR" ]; then
        export WORKDIR=$HOME/workdir
    fi
    echo "Using working directory: $WORKDIR"
    echo ""
    echo "Please enter the location of the install directory where libraries will be installed (blank for $HOME/installdir)"
    if [ -z "$INSTALLDIR" ]; then
        export INSTALLDIR=$HOME/installdir
    fi
    read INSTALLDIR
    echo "Using install directory: $INSTALLDIR"
fi   

set -e
set -x

export JN=-j8

if [ -z "$WORKDIR" ]; then
    export WORKDIR=$HOME/workdir
    mkdir -p $WORKDIR
fi

if [ -z "$INSTALLDIR" ]; then
    export INSTALLDIR=$HOME/installdir
    mkdir -p $INSTALLDIR
fi
export LD_LIBRARY_PATH=${INSTALLDIR}/lib:${LD_LIBRARY_PATH}

function install_zlib {
    echo install_zlib
    cd $WORKDIR

    if [ ! -d "$WORKDIR/zlib-1.3.1" ]; then
        wget --no-check-certificate -q https://www.zlib.net/zlib-1.3.1.tar.gz
        tar -xvzf zlib-1.3.1.tar.gz
    fi

    cd zlib-1.3.1/

    # check if make file exists and if not, run configure
    if [ ! -f "Makefile" ]; then
        ./configure --prefix=$INSTALLDIR &> config.log
    fi

    make -j 8 &> make.log
    make check
    make install
}

function install_hdf5 {
    echo install_hdf5
    cd $WORKDIR

    export CPPFLAGS=-I$INSTALLDIR/include 
    export LDFLAGS=-L$INSTALLDIR/lib
    export CC=mpicc

    if [ ! -d "$WORKDIR/hdf5-1.14.3" ]; then
        wget --no-check-certificate -q https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.14/hdf5-1.14.3/src/hdf5-1.14.3.tar.gz
        tar -xzf hdf5-1.14.3.tar.gz
    fi
    
    cd hdf5-1.14.3

    # check if make file exists and if not, run configure
    if [ ! -f "Makefile" ]; then
        ./configure --prefix=$INSTALLDIR --enable-parallel --with-zlib=$INSTALLDIR #&> config.log
    fi

    make -j 8
    make install

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
    export CC=mpicc
    export LIBS=-ldl

    if [ ! -d "$WORKDIR/netcdf-c-4.9.2" ]; then
        wget --no-check-certificate -q https://github.com/Unidata/netcdf-c/archive/refs/tags/v4.9.2.tar.gz
        tar -xzf v4.9.2.tar.gz

    fi
    
    cd netcdf-c-4.9.2

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

    # check if make file exists and if not, run configure
    if [ ! -f "Makefile" ]; then
        CC=mpicc FC=mpif90 F77=mpif77 ./configure --prefix=${INSTALLDIR} --disable-shared
    fi

    make -j 8
    make install
}

function install_petsc {
    cd $WORKDIR

    if [ ! -d "$WORKDIR/petsc" ]; then
        git clone -b release https://gitlab.com/petsc/petsc.git petsc

        cd petsc
    fi

    cd petsc
    git pull # obtain new release fixes (since a prior clone or pull)

    # check if make file exists and if not, run configure
    if [ ! -f "Makefile" ]; then
        ./configure --prefix=$INSTALLDIR --with-debugging=0 #&> config.log
    fi

    make -j 8
    make check
    make install

}


function hicar_dependencies {
    echo hicar_dependencies
    sudo apt-get update
    sudo apt-get install mpich
    sudo apt-get install libcurl4-gnutls-dev
    sudo apt-get install libfftw3-dev
    # sudo apt-get install petsc-dev

    install_petsc
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
    export PETSC_DIR=/usr #${INSTALLDIR}
    export PATH=${INSTALLDIR}/bin:$PATH
    export LD_LIBRARY_PATH=${INSTALLDIR}/lib:${LD_LIBRARY_PATH}
    cmake ../ -DFSM=OFF
    make ${JN}
    make install
    
    echo "hicar install succeeded"

}

function gen_test_run_data {
    cd ${GITHUB_WORKSPACE}/helpers
    mkdir ${GITHUB_WORKSPACE}/Model_runs/
    printf 'y\ny\n' | ./gen_HICAR_dir.sh ${GITHUB_WORKSPACE}/Model_runs/ ${GITHUB_WORKSPACE}
}

function execute_test_run {
    cd ${GITHUB_WORKSPACE}/Model_runs/HICAR/input
    export LD_LIBRARY_PATH=${INSTALLDIR}/lib:${LD_LIBRARY_PATH}
    echo "Starting HICAR run"
    mpirun -np 2 ${GITHUB_WORKSPACE}/bin/HICAR HICAR_Test_Case.nml

    time_dim=$(ncdump -v ../output/*.nc | grep "time = UNLIMITED" | sed 's/[^0-9]*//g')

    if [[ ${time_dim} == "1" ]]; then
	echo "FAILURE: HICAR output time dimension should not be equal to one, it was ${time_dim}"
	exit 1
    else
	echo "SUCCESS: time dimension is equal to ${time_dim}"
	exit 0
    fi
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
        execute_test_run)
            execute_test_run;;
        *)
            echo "$func unknown";;
    esac
done
