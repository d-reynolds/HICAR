# ---------------------------------------------------------------------------
# HICAR test / CI "make" targets.
#
# Moved out of the root CMakeLists.txt (it was a ~160-line block at the end) so
# the test-runner logic lives under tests/. This file is include()d from the root
# CMakeLists.txt at the very end, AFTER find_package(MPI), the HICAR / hicar_lib
# executables, and add_subdirectory(tests) (which creates HICAR-tester +
# download_test_data). include() runs in the ROOT scope, so PROJECT_SOURCE_DIR is
# the repo root and every target referenced below already exists.
#
# Do NOT turn this into tests/CMakeLists.txt or add_subdirectory() it: that scope
# redefines PROJECT_SOURCE_DIR (its own project() call) and is processed before
# HICAR exists, which would silently turn the DEPENDS into broken file deps.
#
# Targets: check, test_valgrind, test_cli, test_cases, test_decomposition,
# test_restart, test_reproducibility, test_regression, test_gpu, test_snowpack
# (+ helpers HICAR_installed, HICAR_debug). See docs/testing.md.
# ---------------------------------------------------------------------------

# Add a custom target for running tests with MPI

separate_arguments(SRUN_FLAGS_LIST UNIX_COMMAND "${SRUN_FLAGS}")
SET(SRUN_FLAGS_LIST_TEST_CASE "-SRUN_FLAGS='${SRUN_FLAGS}'")

set(TEST_CASES "Standard,Nested,Standard_restart,Nested_restart")

# Build and install HICAR to bin/ before the single-exe test targets run: the
# test scripts invoke bin/HICAR, so depending on this guarantees they use a
# freshly compiled exe rather than a stale install.
add_custom_target(HICAR_installed
    COMMAND ${CMAKE_COMMAND} --install ${CMAKE_BINARY_DIR}
    DEPENDS HICAR
    COMMENT "Building + installing HICAR to bin/"
)

#check if the string "srun" appears in the MPIEXEC_EXECUTABLE variable
if (MPIEXEC_EXECUTABLE MATCHES "srun")
    #see if the variable "account" is not set
    if (NOT DEFINED SRUN_FLAGS)
        message(ERROR "No account specified for srun, please pass your SLURM account name to cmake using -DSRUN_FLAGS=<SRUN_FLAGS>")
    endif()

    add_custom_target(check
        COMMAND ${CMAKE_COMMAND} -E echo "Running code integration testes with srun..."
        COMMAND "${MPIEXEC_EXECUTABLE}" ${MPIEXEC_NUMPROC_FLAG} 4 ${SRUN_FLAGS_LIST} ${CMAKE_BINARY_DIR}/tests/${PROJECT_NAME}-tester
        DEPENDS ${PROJECT_NAME}-tester
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        COMMENT "Running integration test suite with 4 MPI processes"
        VERBATIM
    )

    add_custom_target(test_cases
        COMMAND ${CMAKE_COMMAND} -E echo "Running integration cases: ${TEST_CASES}"
        COMMAND "${PROJECT_SOURCE_DIR}/tests/Test_Cases/test_case_runner.sh" ${PROJECT_SOURCE_DIR} "${TEST_CASES}" ${SRUN_FLAGS_LIST_TEST_CASE}
        DEPENDS HICAR_installed download_test_data
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests/Test_Cases
        COMMENT "Running integration test case suite"
        VERBATIM
    )

else()
    add_custom_target(check
        COMMAND ${CMAKE_COMMAND} -E echo "Running code integration testes with MPI..."
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${CMAKE_BINARY_DIR}/tests/${PROJECT_NAME}-tester
        DEPENDS ${PROJECT_NAME}-tester
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        COMMENT "Running integration test suite with 4 MPI processes"
        VERBATIM
    )

    add_custom_target(test_cases
        COMMAND ${CMAKE_COMMAND} -E echo "Running integration cases: ${TEST_CASES}"
        COMMAND ${PROJECT_SOURCE_DIR}/tests/Test_Cases/test_case_runner.sh ${PROJECT_SOURCE_DIR} "${TEST_CASES}"
        DEPENDS HICAR_installed download_test_data
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests/Test_Cases
        COMMENT "Running integration test case suite"
        VERBATIM
    )

endif()

# Unit tests under valgrind memcheck (tests/test_valgrind.sh): the same
# HICAR-tester `check` runs (mpiexec -np 4), wrapped in valgrind --track-origins.
# This is the CPU  memcheck CI lane (valgrind-memcheck.yml) made runnable locally; needs a debug
# build and `valgrind` installed (the script errors clearly if it is missing).
add_custom_target(test_valgrind
    COMMAND ${CMAKE_COMMAND} -E echo "Running unit tests under valgrind memcheck (mpiexec -np 4)..."
    COMMAND ${PROJECT_SOURCE_DIR}/tests/test_valgrind.sh ${PROJECT_SOURCE_DIR} ${CMAKE_BINARY_DIR}
    DEPENDS ${PROJECT_NAME}-tester
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Running HICAR unit tests under valgrind memcheck"
    VERBATIM
)

add_custom_target(test_cli
    COMMAND ${CMAKE_COMMAND} -E echo "Running CLI option tests..."
    COMMAND ${PROJECT_SOURCE_DIR}/tests/test_cli_options.sh ${PROJECT_SOURCE_DIR}
    DEPENDS HICAR_installed
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests
    COMMENT "Running CLI options test suite"
    VERBATIM
)

# HICAR_debug target: build debug executable for reproducibility testing.
# If already in debug mode, just depend on the main HICAR target.
# Otherwise, spawn a separate debug sub-build that installs HICAR_debug to bin/.
if(MODE STREQUAL "debug" OR MODE STREQUAL "debugslow")
    add_custom_target(HICAR_debug DEPENDS HICAR)
else()
    set(_DEBUG_BUILD_DIR "${CMAKE_BINARY_DIR}/debug_build")
    add_custom_target(HICAR_debug
        COMMAND ${CMAKE_COMMAND} -S ${PROJECT_SOURCE_DIR} -B ${_DEBUG_BUILD_DIR}
            -DMODE=debug
            -DFC=${FC}
            -DFFTW_DIR=${FFTW_DIR}
            -DNETCDF_DIR=${NETCDF_DIR}
            -DMPI_DIR=${MPI_DIR}
            -DFSM_DIR=${FSM_DIR}
            -DFSM=${FSM}
            -DSNOWPACK_DIR=${SNOWPACK_DIR}
            -DSNOWPACK_CPP=${SNOWPACK_CPP}
            -DASSERTIONS=${ASSERTIONS}
            -DOPENACC=${OPENACC}
        COMMAND ${CMAKE_COMMAND} --build ${_DEBUG_BUILD_DIR} --parallel
        COMMAND ${CMAKE_COMMAND} --install ${_DEBUG_BUILD_DIR}
        COMMENT "Building HICAR_debug executable for reproducibility testing"
    )
endif()

add_custom_target(test_decomposition
    COMMAND ${PROJECT_SOURCE_DIR}/tests/test_reproducibility.sh ${PROJECT_SOURCE_DIR} decomposition
    DEPENDS HICAR_debug download_test_data
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests/Test_Cases
    COMMENT "Running MPI domain decomposition reproducibility test (5 vs 10 ranks)"
    VERBATIM
)

add_custom_target(test_restart
    COMMAND ${PROJECT_SOURCE_DIR}/tests/test_reproducibility.sh ${PROJECT_SOURCE_DIR} restart
    DEPENDS HICAR_debug download_test_data
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests/Test_Cases
    COMMENT "Running restart reproducibility test"
    VERBATIM
)

add_custom_target(test_reproducibility
    COMMAND ${PROJECT_SOURCE_DIR}/tests/test_reproducibility.sh ${PROJECT_SOURCE_DIR} all
    DEPENDS HICAR_debug download_test_data
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests/Test_Cases
    COMMENT "Running all reproducibility tests"
    VERBATIM
)

# Regression test: diff current integration output against a "blessed" reference
# commit (tests/test_regression.sh, which resolves the blessed commit via gh). The
# script decides for itself whether to (re)run the integration cases — it reuses
# existing output under tests/Test_Cases/output and only runs `make test_cases`
# when it is missing
add_custom_target(test_regression
    COMMAND ${CMAKE_COMMAND} -E echo "Regression: diff vs the blessed reference commit (mode exact)..."
    COMMAND ${PROJECT_SOURCE_DIR}/tests/test_regression.sh
            ${PROJECT_SOURCE_DIR} ${CMAKE_BINARY_DIR} "${TEST_CASES}" --mode exact
    DEPENDS HICAR_installed download_test_data
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests/Test_Cases
    COMMENT "Regression: diff integration output vs blessed commit"
    VERBATIM
)

# GPU-vs-CPU comparison (tests/compare_cpu_gpu.sh). `make test_gpu` compiles BOTH
# exes fresh (CPU=OPENACC=OFF, GPU=OPENACC=ON, same compiler) into
# build_gpu_compare_cpu / build_gpu_compare_gpu, compares, then removes those
# trees so the next run rebuilds fresh. Needs an NVHPC-configured build
# (FC=nvfortran) for the GPU exe.
add_custom_target(test_gpu
    COMMAND ${CMAKE_COMMAND} -E echo "Building CPU + GPU exes, then comparing (tolerance)..."
    COMMAND ${PROJECT_SOURCE_DIR}/tests/compare_cpu_gpu.sh ${PROJECT_SOURCE_DIR} --build ${CMAKE_BINARY_DIR}
    DEPENDS download_test_data
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests/Test_Cases
    COMMENT "Building CPU/GPU exes then running GPU-vs-CPU comparison"
    VERBATIM
)

# SNOWPACK C++ vs native-Fortran parity (tests/snowpack/test_snowpack_compare.sh).
# `make test_snowpack` compiles BOTH exes fresh (-DSNOWPACK_CPP=ON vs the default
# port, both CPU) into build_snowpack_cpp / build_snowpack_fortran, compares, then
# removes those trees so the next run rebuilds fresh. The --build arg is the build
# directory whose CMakeCache supplies the compiler/library config.
add_custom_target(test_snowpack
    COMMAND ${CMAKE_COMMAND} -E echo "Building C++ + Fortran SNOWPACK exes, then comparing..."
    COMMAND ${PROJECT_SOURCE_DIR}/tests/snowpack/test_snowpack_compare.sh ${PROJECT_SOURCE_DIR} --build ${CMAKE_BINARY_DIR}
    DEPENDS download_test_data
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests/Test_Cases
    COMMENT "Building SNOWPACK C++/Fortran exes then running parity comparison"
    VERBATIM
)
