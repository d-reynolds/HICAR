#!/usr/bin/env bash

set -e
set -x

export FC=gfortran-9
# see link for size of runner
# https://docs.github.com/en/actions/using-github-hosted-runners/about-github-hosted-runners#supported-runners-and-hardware-resources
export JN=-j4

if [ -z "$WORKDIR" ]; then
    export WORKDIR=$HOME/workdir
    mkdir -p $WORKDIR
fi

if [ -z "$INSTALLDIR" ]; then
    export INSTALLDIR=$HOME/installdir
    mkdir -p $INSTALLDIR
fi
export LD_LIBRARY_PATH=${INSTALLDIR}/lib:${LD_LIBRARY_PATH}


function install_szip {
    echo install_szip
    cd $WORKDIR
    wget --no-check-certificate -q http://www.hdfgroup.org/ftp/lib-external/szip/2.1.1/src/szip-2.1.1.tar.gz
    tar -xzf szip-2.1.1.tar.gz
    cd szip-2.1.1
    ./configure --prefix=$INSTALLDIR &> config.log
    make -j 4 &> make.log
    make install
}

function install_zlib {
    echo install_zlib
    cd $WORKDIR
    wget --no-check-certificate -q https://www.zlib.net/zlib-1.3.1.tar.gz
    tar -xvzf zlib-1.3.1.tar.gz
    cd zlib-1.3.1/
    ./configure --prefix=$INSTALLDIR &> config.log
    make -j 4 &> make.log
    make check
    make install
}

function install_hdf5 {
    echo install_hdf5
    cd $WORKDIR
    wget --no-check-certificate -q https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.7/src/hdf5-1.10.7.tar.gz
    tar -xzf hdf5-1.10.7.tar.gz
    cd hdf5-1.10.7
    # FCFLAGS="-DH5_USE_110_API" ./configure --prefix=$INSTALLDIR &> config.log
    CC=mpicc ./configure --prefix=$INSTALLDIR --enable-parallel --with-zlib=$INSTALLDIR #&> config.log
    make -j 4
    make install
    # CFLAGS=-DH5_USE_110_API make
    # (CFLAGS=-DH5_USE_110_API make | awk 'NR%100 == 0')
    export HDF5=$INSTALLDIR
    export LD_LIBRARY_PATH=$INSTALLDIR/lib:$LD_LIBRARY_PATH
}

function install_netcdf_c {
    echo install_netcdf_c
    cd $WORKDIR
    wget --no-check-certificate -q ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-c-4.8.0.tar.gz
    tar -xzf netcdf-c-4.8.0.tar.gz
    cd netcdf-c-4.8.0
    ./configure --prefix=$INSTALLDIR #&> config.log
    make #&> make.log
    make install
}

function install_netcdf_fortran {
    export SRCNCDF=${GITHUB_WORKSPACE}/srcNETCDF
    mkdir $SRCNCDF
    cd $SRCNCDF

    wget --no-check-certificate -q https://github.com/Unidata/netcdf-c/archive/refs/tags/v4.8.0.tar.gz
    wget --no-check-certificate -q https://github.com/Unidata/netcdf-fortran/archive/refs/tags/v4.6.1.tar.gz

    tar -xzf v4.8.0.tar.gz
    tar -xzf v4.6.1.tar.gz

    cd netcdf-c-4.8.0
    export CPPFLAGS=-I$INSTALLDIR/include 
    export LDFLAGS=-L$INSTALLDIR/lib
    cmake ./ -D"NETCDF_ENABLE_PARALLEL4=ON" -D"CMAKE_INSTALL_PREFIX=${INSTALLDIR}"
    make all check install

    export PATH=$INSTALLDIR/bin:$PATH
    export NETCDF=$INSTALLDIR

    cd ../netcdf-fortran-4.6.1
    export CPPFLAGS=$CPPFLAGS" -I${INSTALLDIR}/include"
    export LDFLAGS=$LDFLAGS" -L${INSTALLDIR}/lib"
    export LD_LIBRARY_PATH=${INSTALLDIR}/lib:${LD_LIBRARY_PATH}
    cmake ./—prefix=${INSTALLDIR}
    make check
    Make install
}



function hicar_dependencies {
    echo hicar_dependencies
    sudo apt-get update
    sudo apt-get install mpich
    sudo apt-get install libcurl4-gnutls-dev
    sudo apt-get install libfftw3-dev
    # install_zlib
    # install_hdf5
    sudo apt-get install libhdf5-mpi-dev
    install_netcdf_fortran
    sudo apt-get install petsc-dev

    # Installing HDF5 currently not working for NetCDF
    # sudo apt-get install libhdf5-dev libhdf5-openmpi-dev

    # export CPPFLAGS="$CPPFLAGS -I${INSTALLDIR}/include"
    # export LDFLAGS="$LDFLAGS -L${INSTALLDIR}/lib"

    # # Install szip (used by hdf5)
    # install_szip
    # # Install HDF5
    # install_hdf5

    # # Install NetCDF-C
    # install_netcdf_c
    # # Install NetCDF fortran
    # install_netcdf_fortran

    # put installed bin directory in PATH
    export PATH=${INSTALLDIR}/bin:$PATH
}

function hicar_install {
    echo hicar_install
    pwd
    cd ${GITHUB_WORKSPACE}
    mkdir build
    cd build
    cmake ../ -DFSM=OFF
    make ${JN}
    make install
    
    export NETCDF=${INSTALLDIR}
    #export FFTW=/usr
    #export JN=-j4

    
    # test build
    # make -C src clean; VERBOSE=1 make -C src ${JN} MODE=debugslow

    echo "hicar install succeeded"

}

function hicar_script {

    cd ${GITHUB_WORKSPACE}
    cd ./src

    cp ../run/complete_icar_options.nml ./icar_options.nml
    export OMP_NUM_THREADS=2

    mpiexec -n 4 ./icar
    # for i in *_test; do
    #     echo $i
    #     ./${i}
    # done

    echo "hicar script succeeded"
    cd ..
}

function gen_test_run_data {
    cd ${GITHUB_WORKSPACE}/helpers
    mkdir ${GITHUB_WORKSPACE}/Model_runs/
    echo "y y" | ./gen_HICAR_dir.sh ${GITHUB_WORKSPACE}/Model_runs/ ${GITHUB_WORKSPACE}
    }

function execute_test_run {
    cd ${GITHUB_WORKSPACE}/Model_runs/input

    echo "Starting HICAR run"
    ./bin/HICAR HICAR_Test_Case.nml
    time_dim=$(ncdump -v ../output/*.nc | grep "time = UNLIMITED" | sed 's/[^0-9]*//g')

    if [[ ${time_dim} == "1" ]]; then
	echo "FAILURE: HICAR output time dimension should not be equal to one, it was ${time_dim}"
	exit 1
    else
	echo "SUCCESS: time dimension is equal to ${time_dim}"
	exit 0
    fi
}

function icar_after_success {
  echo icar_after_success
  echo "icar build succeeded"
}

function icar_after_failure {
  echo icar_after_failure
  echo "icar build failed"
}

for func in "$@"
do
    case $func in
        hicar_dependencies)
            hicar_dependencies;;
        hicar_install)
            hicar_install;;
	hicar_script)
	    hicar_script;;
	gen_test_run_data)
	    gen_test_run_data;;
	execute_test_run)
	    execute_test_run;;
        *)
            echo "$func unknown"
    esac
done
