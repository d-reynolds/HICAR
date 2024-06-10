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
    make &> make.log
    make install ${JN}
}

function install_hdf5 {
    echo install_hdf5
    cd $WORKDIR
    wget --no-check-certificate -q https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.7/src/hdf5-1.10.7.tar.gz
    tar -xzf hdf5-1.10.7.tar.gz
    cd hdf5-1.10.7
    # FCFLAGS="-DH5_USE_110_API" ./configure --prefix=$INSTALLDIR &> config.log
    ./configure --prefix=$INSTALLDIR #&> config.log
    make
    # CFLAGS=-DH5_USE_110_API make
    # (CFLAGS=-DH5_USE_110_API make | awk 'NR%100 == 0')
    make install ${JN}
    export LIBDIR=${INSTALLDIR}/lib
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
    export INSNCDF=${GITHUB_WORKSPACE}/NETCDF_INST
    mkdir $SRCNCDF
    mkdir $INSNCDF
    cd $SRCNCDF

    Git clone https://github.com/Unidata/netcdf-c.git
    Git clone https://github.com/Unidata/netcdf-fortran.git

    cd netcdf-c/
    cmake ./ -D"NETCDF_ENABLE_PARALLEL4=ON" -D"CMAKE_INSTALL_PREFIX=${INSNCDF}"
    make all check install

    cd ../netcdf-fortran
    export NCDIR=${INSNCDF}
    export NFDIR=${INSNCDF}
    export CPPFLAGS=$CPPFLAGS" -I${NCDIR}/include"
    export LDFLAGS=$LDFLAGS" -L${NCDIR}/lib"
    export LD_LIBRARY_PATH=${NCDIR}/lib:${LD_LIBRARY_PATH}
    cmake ./—prefix=${NFDIR}
    make check
    Make install
}



function hicar_dependencies {
    echo hicar_dependencies

    sudo apt-get update
    sudo apt-get install mpich
    sudo apt-get install libcurl4-gnutls-dev
    sudo apt-get install libfftw3-dev
    install_szip;;
    install_hdf5;;
    install_netcdf_fortran;;
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
    echo "y y" | ./gen_HICAR_dir.sh ~/ ${GITHUB_WORKSPACE}
    }

function execute_test_run {
    cd ~/HICAR/input

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
