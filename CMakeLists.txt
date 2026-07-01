# ============================================================================
# HICAR root build — organized as dependency-ordered phases:
#   A Project & toolchain
#   B Options & config
#   C Sources & metadata
#   D Compiler identity & flags
#   E Dependencies
#   F Targets
#   G Finalize
# ============================================================================

# ============================================================================
# PHASE A — Project & toolchain
# ============================================================================
# Always set the cmake min version.
cmake_minimum_required(VERSION 3.20)

project(HICAR)
set(CMAKE_CXX_STANDARD 14)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

#Do not use any conda packages
if(DEFINED ENV{CONDA_PREFIX})
  set(CMAKE_IGNORE_PATH $ENV{CONDA_PREFIX} ${CMAKE_IGNORE_PATH})
  message(STATUS "Ignoring Conda environment path: $ENV{CONDA_PREFIX}")
endif()

# Resolve FC to an absolute path BEFORE enable_language(Fortran) so FetchContent
# subprojects (rte-rrtmgp) get an absolute CMAKE_Fortran_COMPILER.
set(_fc_request "")
if (DEFINED FC AND NOT FC STREQUAL "")
	set(_fc_request "${FC}")
elseif (DEFINED ENV{FC} AND NOT "$ENV{FC}" STREQUAL "")
	set(_fc_request "$ENV{FC}")
endif()
if (NOT _fc_request STREQUAL "")
	find_program(FC_FULL_PATH NAMES "${_fc_request}")
	if (FC_FULL_PATH)
		set(_fc_request "${FC_FULL_PATH}")   # bare name -> absolute path
	endif()
	set(CMAKE_Fortran_COMPILER "${_fc_request}" CACHE FILEPATH "Fortran compiler" FORCE)
	unset(FC_FULL_PATH CACHE)
endif()
unset(_fc_request)

enable_language (Fortran)

# Put all fetched external dependencies (noah_mp, rte-rrtmgp, rrtmgp_data,
# icar_supporting_files, test-drive, ...) under build/external — alongside
# SNOWPACK/MeteoIO — instead of CMake's default build/_deps. Must be set before
# any FetchContent use (incl. phase E and the phase-G tests test-drive finder).
set(FETCHCONTENT_BASE_DIR "${CMAKE_BINARY_DIR}/external")

# ============================================================================
# PHASE B — Build options & configuration
# ============================================================================
####### Set up the user variables #######
set(MODE "release" CACHE STRING "Set compilation flags, useful for debugging")
set(FC "" CACHE STRING "Set the fortran compiler to use (can be auto-detected for most cases)")

if (NOT DEFINED FFTW_DIR OR FFTW_DIR STREQUAL "")
	set(FFTW_DIR $ENV{FFTW_DIR} CACHE PATH "Path to FFTW installation")
endif()
if (NOT DEFINED NETCDF_DIR OR NETCDF_DIR STREQUAL "")
	set(NETCDF_DIR $ENV{NETCDF_DIR} CACHE PATH "Path to NETCDF installation")
endif()
if (NOT DEFINED MPI_DIR OR MPI_DIR STREQUAL "")
	set(MPI_DIR $ENV{MPICH_DIR} CACHE PATH "Path to MPI installation")
endif()
if (NOT DEFINED FSM_DIR OR FSM_DIR STREQUAL "")
	set(FSM_DIR $ENV{FSM_DIR} CACHE PATH "Path to the compiled FSM installation (root directory of the lib and build directories)")
endif()

if (NOT DEFINED SNOWPACK_DIR OR SNOWPACK_DIR STREQUAL "")
	set(SNOWPACK_DIR $ENV{SNOWPACK_DIR} CACHE PATH "Path to the compiled SNOWPACK installation (root directory of the lib and include directories)")
endif()
if (NOT DEFINED METEOIO_DIR OR METEOIO_DIR STREQUAL "")
	set(METEOIO_DIR $ENV{METEOIO_DIR} CACHE PATH "Path to the compiled MeteoIO installation (root directory of the lib and include directories)")
endif()

if (NOT "$ENV{NOAHMP_LIBRARY}" STREQUAL "")
	set(NOAHMP_LIBRARY $ENV{NOAHMP_LIBRARY} CACHE PATH "Path to the Noah-MP source directory")
endif()

if (NOT DEFINED MPICH_DIR AND DEFINED ENV{MPICH_DIR})
    set(MPI_DIR $ENV{MPICH_DIR} CACHE STRING "Path to MPI directory" FORCE)
endif()


# Track if user explicitly set OPENACC flag before option() creates cache entry
if(DEFINED OPENACC)
	set(OPENACC_SET_BY_USER TRUE)
else()
	set(OPENACC_SET_BY_USER FALSE)
endif()

option(ASSERTIONS "Check logical assertions at runtime" ON)
option(FSM "Link to the FSM2 snow model, requires separate compilation of the FSM2 model" OFF)
option(SNOWPACK_CPP "Use the C++ SNOWPACK wrapper instead of the native-Fortran port" OFF)
option(NCCL "Link to the NCCL library, if available" ON)
option(OPENACC "Enable OpenACC support for GPU acceleration" OFF)
set(SRUN_FLAGS "" CACHE STRING "srun flags to use for running the test cases")

# SNOWPACK is always compiled into HICAR; SNOWPACK_CPP is the single switch that
# selects the implementation (the whole build keys off it directly — there is no
# internal SNOWPACK / SNOWPACK_FORTRAN variable):
#   default            -> native-Fortran SNOWPACK port (snowpack_driver.F90)
#   -DSNOWPACK_CPP=ON  -> C++ SNOWPACK wrapper (sm_SNOWPACK.F90, -DSNOWPACK_CPP)
#
# Removed/deprecated flags (warn, then honor the closest modern equivalent):
if(DEFINED SNOWPACK AND NOT SNOWPACK)
    message(WARNING "-DSNOWPACK=OFF is no longer supported: SNOWPACK is always compiled "
        "(native-Fortran port by default). Pass -DSNOWPACK_CPP=ON for the C++ wrapper.")
endif()
if(DEFINED SNOWPACK_FORTRAN)
    message(DEPRECATION "SNOWPACK_FORTRAN is deprecated: the native-Fortran SNOWPACK port is "
        "now the build default. Pass -DSNOWPACK_CPP=ON to select the C++ wrapper instead.")
    if(NOT SNOWPACK_FORTRAN)
        set(SNOWPACK_CPP ON CACHE BOOL "Use the C++ SNOWPACK wrapper" FORCE)
    endif()
    unset(SNOWPACK_FORTRAN CACHE)
endif()
#########################################

add_subdirectory(${PROJECT_SOURCE_DIR}/tests EXCLUDE_FROM_ALL)

# ============================================================================
# PHASE C — Sources & build metadata
# ============================================================================
########## Set up the sources #########
FILE(GLOB_RECURSE SOURCE_FILES  ${PROJECT_SOURCE_DIR}/src/*.F90)

FILE(GLOB_RECURSE ra_mo_files  ${PROJECT_SOURCE_DIR}/src/physics/ra_mo_*.F90)
list(REMOVE_ITEM SOURCE_FILES ${ra_mo_files})


macro(Disable_FSM)
	set(FSM OFF)
	FILE(GLOB_RECURSE FSM2_INTERFACE  ${PROJECT_SOURCE_DIR}/src/physics/FSM2_interface/*)
	list(REMOVE_ITEM SOURCE_FILES  ${PROJECT_SOURCE_DIR}/src/physics/sm_FSM.F90)
	list(REMOVE_ITEM SOURCE_FILES ${FSM2_INTERFACE})
endmacro()

#-DFSM -- controls compilation of FSM2 driver and linking to the FSM2 library
IF (FSM AND IS_DIRECTORY ${FSM_DIR}/lib AND
			IS_DIRECTORY ${FSM_DIR}/build)
elseif(FSM AND NOT (IS_DIRECTORY ${FSM_DIR}/lib AND
					IS_DIRECTORY ${FSM_DIR}/build))
	message(WARNING "FSM2 not found, disabling FSM2 support.\nTo enable, set the FSM_DIR variable to the root directory of the FSM2 install (i.e. the directory which contains the 'lib' and 'build' subdirectories).")
	Disable_FSM()
else()
	Disable_FSM()
endif()

find_package(Git QUIET)
if(GIT_EXECUTABLE)
	execute_process(
		COMMAND ${GIT_EXECUTABLE} describe --long --dirty --all --always
		WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
		OUTPUT_VARIABLE GIT_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
	string(REGEX REPLACE "^heads/" "" GIT_VERSION "${GIT_VERSION}")

	# GIT_COMMIT = exact short commit; GIT_TAG = nearest `v*` release tag + offset/dirty
	# (--tags so lightweight tags count, --match v* to skip dev tags; bare hash if none).
	execute_process(
		COMMAND ${GIT_EXECUTABLE} rev-parse --short=12 HEAD
		WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
		OUTPUT_VARIABLE GIT_COMMIT_HASH OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
	execute_process(
		COMMAND ${GIT_EXECUTABLE} describe --tags --match v* --always --dirty
		WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
		OUTPUT_VARIABLE GIT_DESCRIBE_TAG OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
endif()
if(NOT GIT_VERSION)
	set(GIT_VERSION "unknown")
endif()
if(NOT GIT_COMMIT_HASH)
	set(GIT_COMMIT_HASH "unknown")
endif()
if(NOT GIT_DESCRIBE_TAG)
	set(GIT_DESCRIBE_TAG "untagged")
endif()

# ============================================================================
# PHASE D — Compiler identity & flags
# ============================================================================
############## Auto-detect  ################
message(STATUS "")
message(STATUS "Trying to auto-detect the environment...")
list(APPEND CMAKE_MESSAGE_INDENT "  ")


# enable_language(Fortran) above resolved CMAKE_Fortran_COMPILER (the FC we forced, or
# CMake's auto-detection). Classify by executable basename — a directory component
# containing e.g. "gfortran" can't misclassify.
message(STATUS "Using Fortran compiler: ${CMAKE_Fortran_COMPILER}")
get_filename_component(FC_NAME "${CMAKE_Fortran_COMPILER}" NAME)
if (FC_NAME MATCHES "gfortran")
	set(COMPILER "gnu")
elseif (FC_NAME MATCHES "mpif90")
	set(COMPILER "gnu")
elseif (FC_NAME MATCHES "mpifort")
	set(COMPILER "gnu")
elseif (FC_NAME MATCHES "nvfortran")
	set(COMPILER "nvfortran")
	# Auto-enable OpenACC for NVFortran if not explicitly set
	if (NOT OPENACC_SET_BY_USER)
		set(OPENACC ON CACHE BOOL "Enable OpenACC support for GPU acceleration" FORCE)
		message(STATUS "Auto-enabling OpenACC support for NVFortran compiler")
	endif()

else()
	message(FATAL_ERROR "Unknown Fortran compiler (FC) specified")
endif()
unset(FC_NAME)

string( TOLOWER "${MODE}" MODE )
string( TOLOWER "${COMPILER}" COMPILER )

list(POP_BACK CMAKE_MESSAGE_INDENT)
message(STATUS "...done auto-detecting the environment")



###########################################################
############### Define Compiler Flags Once ###############
message(STATUS "")
message(STATUS "Setting up compiler flags...")

# Define compiler flags based on MODE and COMPILER
# These will be used consistently throughout the project
# Define as a list (semicolon-separated) for proper handling in CMake
# Optimization flag tail shared by the gnu release + profile builds.
set(_gnu_opt_tail -ffast-math -ffree-line-length-none -ftree-vectorize -fimplicit-none -funroll-loops -march=native -fno-protect-parens)
if(MODE STREQUAL "debug")
	set(CMAKE_BUILD_TYPE Debug)

	if (COMPILER STREQUAL "gnu")
		# -finit-* poison uninitialised scalars so reads of them are deterministic
		# (and usually fatal in tests) rather than passing by luck on one platform
		# and failing on another. -finit-real=nan alone misses logicals/integers:
		# an uninitialised `restart` logical caused a CI-only control_flow failure.
		set(HICAR_Fortran_FLAGS -O0 -g -fopenmp -fbounds-check -fcheck=all -Wall -fbacktrace -finit-real=nan -finit-integer=-99999999 -finit-logical=true -finit-derived -ffree-line-length-none -Wconversion)
	elseif (COMPILER STREQUAL "nvfortran")
		if(OPENACC)
			set(HICAR_Fortran_FLAGS -O0 -g -traceback -Mbounds -Minfo=accel -Mpreprocess -march=native -acc=gpu -gpu=mem:separate,ccnative)
		else()
			set(HICAR_Fortran_FLAGS -O0 -g -traceback -mp -Mbounds -Mpreprocess -march=native)
		endif()
	endif()
elseif(MODE STREQUAL "debugslow")
	set(CMAKE_BUILD_TYPE Debug)

	if (COMPILER STREQUAL "gnu")
		set(HICAR_Fortran_FLAGS -g -ffree-line-length-none -Wconversion)
	endif()
elseif(MODE STREQUAL "profile")
	set(CMAKE_BUILD_TYPE Profile)

	if (COMPILER STREQUAL "gnu")
		set(HICAR_Fortran_FLAGS -g -O3 -fopenmp -fbounds-check -fbacktrace -finit-real=nan -fopt-info-missed=opt_info.txt ${_gnu_opt_tail})
	endif()
elseif(MODE STREQUAL "release")
	set(CMAKE_BUILD_TYPE Release)
	if (COMPILER STREQUAL "gnu")
		set(HICAR_Fortran_FLAGS -O3 -fopenmp ${_gnu_opt_tail})
	elseif (COMPILER STREQUAL "nvfortran")
		if(OPENACC)
			set(HICAR_Fortran_FLAGS -O3 -Minfo=accel -march=native -acc=gpu -gpu=mem:separate,ccnative)
		else()
			set(HICAR_Fortran_FLAGS -O3 -mp -march=native)
		endif()
	endif()
endif()

# Fail loudly instead of silently building with NO flags: unknown MODE, or a
# MODE×COMPILER combo with no branch above (e.g. nvfortran + debugslow/profile).
if(NOT DEFINED HICAR_Fortran_FLAGS)
	message(FATAL_ERROR "No compiler flags for MODE='${MODE}' with COMPILER='${COMPILER}'. "
		"Supported: debug/release (gnu + nvfortran); debugslow/profile (gnu only).")
endif()

# Convert list to space-separated string for display and external projects
string(REPLACE ";" " " HICAR_Fortran_FLAGS_STR "${HICAR_Fortran_FLAGS}")
message(STATUS "Fortran compiler flags: ${HICAR_Fortran_FLAGS_STR}")
message(STATUS "...done setting up compiler flags")

# ============================================================================
# PHASE E — Dependencies (find_package + FetchContent)
# ============================================================================
############### Find or install dependencies ##############
message(STATUS "")
message(STATUS "Searching for dependencies...")
list(APPEND CMAKE_MESSAGE_INDENT "  ")

include(GNUInstallDirs)


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake" ${FFTW_DIR} ${NETCDF_DIR}  ${MPI_DIR} ${FSM_DIR} ${SNOWPACK_DIR})
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${FFTW_DIR} ${NETCDF_DIR} ${MPI_DIR} ${FSM_DIR} ${SNOWPACK_DIR})

find_package(PkgConfig REQUIRED)
include(FetchContent)

# 6. Link to the package libraries
set(NETCDF_F90 "YES")
find_package(MPI REQUIRED)
find_package(OpenMP)
find_package(FFTW REQUIRED)
find_package(NetCDF REQUIRED)
if (NOT SNOWPACK_CPP)
	# Native-Fortran port (default): compile the pure-Fortran snowpack modules
	# directly, no C++ SNOWPACK needed.
	if(NOT SNOWPACK_DIR OR SNOWPACK_DIR STREQUAL "")
		find_package(SNOWPACK REQUIRED)
	endif()
	set(SNOWPACK_FORTRAN_SRC_DIR "${SNOWPACK_DIR}/fortran")
	file(GLOB SNOWPACK_FORTRAN_SOURCES "${SNOWPACK_FORTRAN_SRC_DIR}/snowpack_*.F90")
	if(NOT SNOWPACK_FORTRAN_SOURCES)
		message(FATAL_ERROR "No snowpack_*.F90 files found in ${SNOWPACK_FORTRAN_SRC_DIR}")
	endif()
	list(APPEND SOURCE_FILES ${SNOWPACK_FORTRAN_SOURCES})
	# Remove the C++ SNOWPACK Fortran driver (replaced by snowpack_driver)
	list(REMOVE_ITEM SOURCE_FILES ${PROJECT_SOURCE_DIR}/src/physics/sm_SNOWPACK.F90)
	message(STATUS "SNOWPACK sources (native-Fortran port): ${SNOWPACK_FORTRAN_SOURCES}")
else()
	# C++ SNOWPACK wrapper: keep sm_SNOWPACK.F90, link the C++ libs.
	find_package(SNOWPACK REQUIRED)
endif()

#########################################################
############# Declare all content used ##################
#########################################################

# rte-rrtmgp library for OpenACC support
FetchContent_Declare(
	rte-rrtmgp
	GIT_REPOSITORY    https://github.com/earth-system-radiation/rte-rrtmgp.git
	GIT_TAG           v1.9.2
	GIT_SHALLOW       TRUE
)

# rrtgmp supporting data for rrtmgp scheme
FetchContent_Declare(
	rrtmgp_data
	GIT_REPOSITORY    https://github.com/earth-system-radiation/rrtmgp-data.git
	GIT_TAG           main
	GIT_SHALLOW       TRUE
)

# Noah-MP library
FetchContent_Declare(
	noah_mp
	GIT_REPOSITORY    https://github.com/d-reynolds/noahmp.git
	GIT_TAG           OpenACC
	GIT_SHALLOW       TRUE
)

# Fetch specific folders from icar_supporting_files repo
FetchContent_Declare(
	icar_supporting_files
	GIT_REPOSITORY    https://github.com/NCAR/icar_supporting_files.git
	GIT_TAG           main
)


#########################################################
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${HICAR_Fortran_FLAGS_STR}")
# Handle OpenACC flag and library detection
if(OPENACC)
	find_package(OpenACC)
	if(NOT OpenACC_Fortran_FOUND)
		message(WARNING "OpenACC support requested but OpenACC library not found. Disabling OpenACC support.")
		set(OPENACC OFF CACHE BOOL "Enable OpenACC support for GPU acceleration" FORCE)
	else()
		message(STATUS "OpenACC support enabled")
	endif()
endif()

# Fetch and build rte-rrtmgp (used for both CPU and GPU radiation)
message(STATUS "Fetching rte-rrtmgp library...")

# Configure build options for rte-rrtmgp
set(BUILD_TESTING OFF)
set(BUILD_C_HEADERS OFF)
# Select the rte-rrtmgp kernel implementation:
#   accel   -> OpenACC/OpenMP-offload kernels for GPU builds
#   default -> portable Fortran kernels for CPU (multiprocessor) builds
if(OPENACC)
	set(KERNEL_MODE "accel" CACHE STRING "Select the kernel mode for rte-rrtmgp" FORCE)
else()
	set(KERNEL_MODE "default" CACHE STRING "Select the kernel mode for rte-rrtmgp" FORCE)
endif()
message(STATUS "rte-rrtmgp KERNEL_MODE: ${KERNEL_MODE}")
set(BUILD_SHARED_LIBS OFF)

# rte-rrtmgp's frontend carries OpenMP *target-offload* (`!$omp target ... map`)
# directives next to its OpenACC ones. The CPU (default) kernels are plain
# Fortran and need no OpenMP, but HICAR's gfortran flags include -fopenmp (for the
# host-threaded snow models). Inheriting -fopenmp makes gfortran try to compile
# those offload directives and fail ("Unclassifiable OpenMP directive"). Build
# rte-rrtmgp with -fopenmp stripped; HICAR and the other deps keep it. (No-op for
# the nvfortran GPU build, which drives OpenACC with -acc and has no -fopenmp, so
# the accel kernels are unaffected.) The add_subdirectory that creates the
# rte-rrtmgp targets runs inside FetchContent_MakeAvailable, so strip around it.
set(_hicar_saved_fortran_flags "${CMAKE_Fortran_FLAGS}")

if (COMPILER STREQUAL "gnu")
	string(REPLACE "-fopenmp" "" CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
elseif (COMPILER STREQUAL "nvfortran")
	string(REPLACE "-mp" "" CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
endif()

FetchContent_MakeAvailable(rte-rrtmgp)
	
# Use the predefined HICAR compiler flags, with openmp stripped out			
message(STATUS "Using compiler flags for rte-rrtmgp: ${CMAKE_Fortran_FLAGS}")

if(NOT TARGET rte-rrtmgp::rte AND NOT TARGET rte-rrtmgp::rrtmgp)
	# Add the rte-rrtmgp subdirectory with configured options
	add_subdirectory(${rte-rrtmgp_SOURCE_DIR} ${rte-rrtmgp_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()

# Restore HICAR's full flags (incl. -fopenmp) for HICAR and the remaining deps.
set(CMAKE_Fortran_FLAGS "${_hicar_saved_fortran_flags}")

list(APPEND SOURCE_FILES ${ra_mo_files})

message(STATUS "rte-rrtmgp library configured successfully")
		

if (NOT DEFINED NOAHMP_LIBRARY)
	# Now fetch Noah-MP repo
	message(STATUS "Fetching External Noah-MP library...")

	FetchContent_MakeAvailable(noah_mp)

	# Use the predefined HICAR compiler flags
	message(STATUS "Using compiler flags for Noah-MP: ${CMAKE_Fortran_FLAGS}")

	# Add the noah_mp subdirectory with configured options
	if(NOT TARGET noahmp)
		add_subdirectory(${noah_mp_SOURCE_DIR}/drivers/hicar ${noah_mp_BINARY_DIR} EXCLUDE_FROM_ALL)
	endif()

	# copy NoahmpTable.TBL from the ${noah_mp_SOURCE_DIR}/parameters into the run directory
	file(COPY ${noah_mp_SOURCE_DIR}/parameters/NoahmpTable.TBL DESTINATION ${PROJECT_SOURCE_DIR}/run)
	message(STATUS "Noah-MP library configured successfully")
else()
	message(STATUS "Using precompiled Noah-MP library specified by NOAHMP_LIBRARY environment variable: ${NOAHMP_LIBRARY}")
endif()

############### Download supporting run files ##############
# (data-only fetches; grouped with the other FetchContent here in phase E — they
# only need their Declare + FETCHCONTENT_BASE_DIR, no flags or targets.)

# Download the repo but don't build it
FetchContent_MakeAvailable(icar_supporting_files)

# Copy only specific folders to the run directory
file(COPY "${icar_supporting_files_SOURCE_DIR}/mp_support" DESTINATION "${PROJECT_SOURCE_DIR}/run")
file(COPY "${icar_supporting_files_SOURCE_DIR}/rrtmg_support" DESTINATION "${PROJECT_SOURCE_DIR}/run")

# Clean up downloaded repo to save space
file(REMOVE_RECURSE "${icar_supporting_files_SOURCE_DIR}")
file(REMOVE_RECURSE "${FETCHCONTENT_BASE_DIR}/icar_supporting_files-subbuild")
file(REMOVE_RECURSE "${FETCHCONTENT_BASE_DIR}/icar_supporting_files-build")

# if we are using RTE_RRTMGP, copy supporting files from the git repository: https://github.com/earth-system-radiation/rrtmgp-data
message(STATUS "Fetching rte-rrtmgp data files...")

# Download the repo but don't build it
FetchContent_MakeAvailable(rrtmgp_data)

# Copy rrtmgp data files to the run directory
file(COPY "${rrtmgp_data_SOURCE_DIR}/" DESTINATION "${PROJECT_SOURCE_DIR}/run/rrtmgp_support")

# Clean up downloaded repo to save space
file(REMOVE_RECURSE "${rrtmgp_data_SOURCE_DIR}")
file(REMOVE_RECURSE "${FETCHCONTENT_BASE_DIR}/rrtmgp_data-subbuild")
file(REMOVE_RECURSE "${FETCHCONTENT_BASE_DIR}/rrtmgp_data-build")

message(STATUS "RTE-RRTMGP data files copied to run/rrtmgp_support")

list(POP_BACK CMAKE_MESSAGE_INDENT)
message(STATUS "...done searching for dependencies")

# ============================================================================
# PHASE F — Targets (hicar_lib + HICAR: definitions, links, options)
# ============================================================================
# SOURCE_FILES is final here (globs + FSM removes + SNOWPACK-native append + ra_mo
# re-append all happened in phases C/E); create the targets from it.
list(REMOVE_ITEM SOURCE_FILES ${PROJECT_SOURCE_DIR}/src/main/driver.F90)
add_library(hicar_lib STATIC ${SOURCE_FILES})
add_executable(HICAR ${PROJECT_SOURCE_DIR}/src/main/driver.F90)

# HICAR-only compile definitions (target-scoped so they don't leak into the
# FetchContent subprojects). PUBLIC so they also reach the HICAR executable.
if(ASSERTIONS)
	target_compile_definitions(hicar_lib PUBLIC USE_ASSERTIONS=.true.)
else()
	target_compile_definitions(hicar_lib PUBLIC USE_ASSERTIONS=.false.)
endif()
target_compile_definitions(hicar_lib PUBLIC
	VERSION="${GIT_VERSION}"
	GIT_COMMIT="${GIT_COMMIT_HASH}"
	GIT_TAG="${GIT_DESCRIBE_TAG}"
	USE_RTE_RRTMGP)

if(OpenMP_CXX_FOUND)
        target_link_libraries(hicar_lib PUBLIC OpenMP::OpenMP_CXX)
elseif(OpenMP_Fortran_FOUND)
	target_link_libraries(hicar_lib PUBLIC OpenMP::OpenMP_Fortran)
endif()

if (FFTW_FOUND)
	target_link_libraries(hicar_lib PUBLIC ${FFTW_LIBRARIES})
	target_include_directories(hicar_lib PUBLIC ${FFTW_INCLUDES})
	if (OPENACC AND OpenACC_Fortran_FOUND)
		# cuFFT for direct GPU FFT calls in linear wind LUT computation
		find_library(CUFFT_LIBRARY NAMES cufft
			HINTS ENV NVHPC_ROOT ENV CUDA_HOME ENV CUDA_PATH
			PATH_SUFFIXES math_libs/lib64 lib64 lib)
		if (NOT CUFFT_LIBRARY)
			find_library(CUFFT_LIBRARY NAMES cufft)
		endif()
		if (CUFFT_LIBRARY)
			message(STATUS "  cuFFT found: ${CUFFT_LIBRARY}")
			target_link_libraries(hicar_lib PUBLIC ${CUFFT_LIBRARY})
		else()
			message(WARNING "cuFFT not found — GPU-accelerated FFT in linear winds will not be available")
		endif()
	endif()
endif()

if (NETCDF_FOUND)
    # Check for parallel NetCDF support
    execute_process(COMMAND nc-config --has-parallel
                    OUTPUT_VARIABLE NETCDF_HAS_PARALLEL
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                    ERROR_QUIET)
    
    if(NETCDF_HAS_PARALLEL STREQUAL "yes")
        message(STATUS "NetCDF built with parallel support")
        target_compile_definitions(hicar_lib PUBLIC NETCDF_PARALLEL)
    else()
        message(STATUS "NetCDF built WITHOUT parallel support")
    endif()
    
	find_program(NF_CONFIG nf-config)
	if(NF_CONFIG)
		execute_process(COMMAND ${NF_CONFIG} --flibs 
						OUTPUT_VARIABLE NETCDF_LIBRARIES 
						OUTPUT_STRIP_TRAILING_WHITESPACE
						ERROR_QUIET)
		separate_arguments(NETCDF_LIBRARIES)

		# Clean up any duplicate libraries
		list(REMOVE_DUPLICATES NETCDF_LIBRARIES)

		# nf-config --flibs emits a bare "-lnetcdf" with no matching -L on split-prefix
		# installs (netcdf-c and netcdf-fortran in different roots). Resolve to an absolute
		# path so the linker doesn't have to search, and so static linking works.
		#
		# Prefer the netcdf-c library FindNetCDF already resolved from NETCDF_DIR with
		# NO_DEFAULT_PATH. A plain find_library here (no NO_DEFAULT_PATH) searches the
		# system default paths too, so a stray homebrew/conda libnetcdf.dylib can shadow
		# the intended (often parallel) netcdf and silently break parallel I/O at runtime
		# ("Parallel operation on file opened for non-parallel access"). Fall back to a
		# NETCDF_DIR-scoped search, then the broad search, only if that is unavailable.
		if(NETCDF_LIBRARIES_C)
			set(NETCDF_C_LIBRARY "${NETCDF_LIBRARIES_C}")
		else()
			find_library(NETCDF_C_LIBRARY NAMES netcdf
				HINTS "${NETCDF_DIR}/lib" "$ENV{NETCDF_DIR}/lib" NO_DEFAULT_PATH)
			if(NOT NETCDF_C_LIBRARY)
				find_library(NETCDF_C_LIBRARY NAMES netcdf
					HINTS ${NETCDF_DIR} ENV NETCDF_DIR ENV NETCDF_C_DIR ENV NETCDF
					PATH_SUFFIXES lib lib64)
			endif()
		endif()
		# Replace -lnetcdf with the absolute path at the *same position* so that
		# libnetcdf.a stays before -lpnetcdf and -lhdf5 in link order.  APPEND
		# would push it past those libraries, breaking single-pass archive resolution.
		list(FIND NETCDF_LIBRARIES -lnetcdf _netcdf_c_pos)
		list(REMOVE_ITEM NETCDF_LIBRARIES -lnetcdf)
		if(NETCDF_C_LIBRARY)
			if(_netcdf_c_pos GREATER -1)
				list(INSERT NETCDF_LIBRARIES ${_netcdf_c_pos} ${NETCDF_C_LIBRARY})
			else()
				list(APPEND NETCDF_LIBRARIES ${NETCDF_C_LIBRARY})
			endif()
			message(STATUS "Resolved netcdf-c library: ${NETCDF_C_LIBRARY}")
		else()
			message(FATAL_ERROR
				"netcdf-c library not found. nf-config --flibs requested -lnetcdf but the "
				"linker cannot locate it. Set NETCDF_C_DIR (or NETCDF_DIR) to the netcdf-c "
				"install prefix.")
		endif()

		target_link_libraries(hicar_lib PUBLIC ${NETCDF_LIBRARIES})
		message(STATUS "NetCDF libraries (with dependencies): ${NETCDF_LIBRARIES}")
	else()
		message(WARNING "nf-config not found, falling back to manual dependency resolution")
		target_link_libraries(hicar_lib PUBLIC ${NETCDF_LIBRARIES})
		find_package(HDF5 COMPONENTS C HL)
		if(HDF5_FOUND)
			target_link_libraries(hicar_lib PUBLIC ${HDF5_LIBRARIES})
		endif()
	endif()
    
    target_include_directories(hicar_lib PUBLIC ${NETCDF_INCLUDES})
endif()


if (FSM)
	find_library(FSM2_LIB NAMES FSM2 REQUIRED)
	target_include_directories(hicar_lib PUBLIC ${FSM_DIR}/build)
	target_link_libraries(hicar_lib PUBLIC ${FSM2_LIB})
	target_compile_definitions(hicar_lib PUBLIC FSM)
endif()

if(NOT SNOWPACK_CPP)
    # Native-Fortran SNOWPACK port (default): no C++ libraries needed and no
    # compile define — the sources select the port path with #ifndef SNOWPACK_CPP.
    message(STATUS "SNOWPACK support enabled (native-Fortran port, no C++ SNOWPACK/MeteoIO needed)")
elseif(SNOWPACK_FOUND)
    # If SNOWPACK should be built from source, add dependency
    if(BUILD_SNOWPACK_FROM_SOURCE)
        add_dependencies(hicar_lib SNOWPACK)
    endif()
    
	# Set RPATH to find SNOWPACK and MeteoIO libraries at runtime
	# Use BUILD_WITH_INSTALL_RPATH to avoid install_name_tool errors on macOS
	set_target_properties(HICAR PROPERTIES
		INSTALL_RPATH "${SNOWPACK_DIR}/lib;${METEOIO_DIR}/lib"
		BUILD_WITH_INSTALL_RPATH TRUE
	)

    # Link SNOWPACK and MeteoIO libraries
    target_link_libraries(hicar_lib PUBLIC ${SNOWPACK_FORTRAN_LIBRARY} ${SNOWPACK_LIBRARY} ${METEOIO_LIBRARY})
    target_include_directories(hicar_lib PUBLIC ${SNOWPACK_INCLUDE_DIR})

    # The C++ SNOWPACK Fortran bindings (static libsnowpack_fortran.a) reference
    # C++ standard-library symbols (operator new, __cxa_throw, std::logic_error,
    # ...), but the final HICAR executable is linked by the Fortran driver, which
    # does not pull in the C++ runtime. Link it explicitly: libc++ on macOS/clang,
    # libstdc++ with GNU/Linux. CXX_STDLIB is set by FindSNOWPACK when SNOWPACK is
    # built from source; fall back to a platform default when it is pre-installed.
    if(NOT CXX_STDLIB)
        if(APPLE)
            set(CXX_STDLIB "c++")
        else()
            set(CXX_STDLIB "stdc++")
        endif()
    endif()
    target_link_libraries(hicar_lib PUBLIC ${CXX_STDLIB})
    
	# C++ SNOWPACK wrapper: -DSNOWPACK_CPP selects the C++ wrapper path
	# in the sources (#ifdef SNOWPACK_CPP).
	target_compile_definitions(hicar_lib PUBLIC SNOWPACK_CPP)
    message(STATUS "SNOWPACK support enabled (C++ SNOWPACK wrapper)")
endif()
if (OpenACC_Fortran_FOUND AND OPENACC)
	target_link_libraries(hicar_lib PUBLIC OpenACC::OpenACC_Fortran)

	# Search for NCCL for GPU halo exchange acceleration
	find_path(NCCL_INCLUDE_DIR nccl.h
		HINTS $ENV{NVHPC_ROOT}/comm_libs/nccl/include
		      $ENV{NCCL_DIR}/include
		      $ENV{NCCL_ROOT}/include)
	find_library(NCCL_LIBRARY nccl
		HINTS $ENV{NVHPC_ROOT}/comm_libs/nccl/lib
		      $ENV{NCCL_DIR}/lib
		      $ENV{NCCL_ROOT}/lib)
	if(NCCL_INCLUDE_DIR AND NCCL_LIBRARY AND NCCL)
		message(STATUS "Found NCCL: ${NCCL_LIBRARY}")
		target_compile_definitions(hicar_lib PUBLIC USE_NCCL)
		# NCCL requires CUDA runtime (cudaSetDevice, cudaFree, cudaStream*, etc.)
		if(NOT TARGET CUDA::cudart)
			enable_language(CUDA)
			find_package(CUDAToolkit REQUIRED)
		endif()
		target_include_directories(hicar_lib PRIVATE ${NCCL_INCLUDE_DIR})
		target_link_libraries(hicar_lib PUBLIC ${NCCL_LIBRARY} CUDA::cudart)
		target_sources(hicar_lib PRIVATE
			${PROJECT_SOURCE_DIR}/src/utilities/nccl_wrapper.c
			${PROJECT_SOURCE_DIR}/src/utilities/nccl_interface.F90)
	else()
		message(STATUS "NCCL not found -- using GPU-aware MPI for halo exchange")
	endif()
endif()

# Link rte-rrtmgp libraries
if(TARGET rte-rrtmgp::rte AND TARGET rte-rrtmgp::rrtmgp)
	target_link_libraries(hicar_lib PUBLIC rte-rrtmgp::rte rte-rrtmgp::rrtmgp)
	message(STATUS "Linked rte-rrtmgp libraries to HICAR")
endif()

target_link_libraries(hicar_lib PUBLIC MPI::MPI_Fortran)
target_include_directories(hicar_lib PUBLIC ${MPI_INCLUDE_DIR})

# Link Noah-MP(Required) to hicar_lib 
if(NOT DEFINED NOAHMP_LIBRARY)
	target_link_libraries(hicar_lib PUBLIC noahmp)
else()
	target_link_libraries(hicar_lib PUBLIC ${NOAHMP_LIBRARY})
	target_include_directories(hicar_lib PUBLIC $ENV{NOAHMP_INCLUDES})
endif()

# link the hicar_lib library to the HICAR executable
target_link_libraries(HICAR PUBLIC hicar_lib)
###########################################################

##########################################################
############## Set up the compiler options ###############
target_compile_options(hicar_lib INTERFACE ${HICAR_Fortran_FLAGS})

# Executable name: debug/debugslow -> HICAR_debug, profile -> HICAR_profile,
# release -> HICAR; nvfortran GPU (OPENACC) builds get a _gpu suffix (not for profile).

if(MODE MATCHES "debug")
	set(_hicar_name "HICAR_debug")
elseif(MODE STREQUAL "profile")
	set(_hicar_name "HICAR_profile")
else()
	set(_hicar_name "HICAR")
endif()
if(COMPILER STREQUAL "nvfortran" AND OPENACC AND NOT MODE STREQUAL "profile")
	set(_hicar_name "${_hicar_name}_gpu")
endif()
set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "${_hicar_name}")
unset(_hicar_name)


# sm_* files (snow model drivers) need OpenMP but not OpenACC (for CPU snow models).
# With the native-Fortran SNOWPACK port, sm_driver.F90 calls GPU kernels and needs OpenACC.
if(COMPILER STREQUAL "nvfortran" AND OPENACC)
	if(SNOWPACK_CPP)
		# CPU snow models (C++ SNOWPACK wrapper): disable OpenACC on all sm_* files
		file(GLOB SM_SOURCES ${PROJECT_SOURCE_DIR}/src/physics/sm_*.F90)
		foreach(source ${SM_SOURCES})
			set_source_files_properties(${source} PROPERTIES COMPILE_OPTIONS "-mp")
		endforeach()
		# sm_SNOWPACK.F90 has OpenMP parallel regions but no OpenACC directives.
		# Disable OpenACC to prevent CUDA context initialization for OpenMP threads.
		if(EXISTS ${PROJECT_SOURCE_DIR}/src/physics/sm_SNOWPACK.F90)
			set_source_files_properties(
				${PROJECT_SOURCE_DIR}/src/physics/sm_SNOWPACK.F90
				PROPERTIES COMPILE_OPTIONS "-mp;-noacc")
		endif()
	else()
		# GPU SNOWPACK: sm_driver.F90 needs OpenACC (calls snowpack_step with !$acc present),
		# but sm_FSM.F90 still needs OpenMP-only
		file(GLOB SM_FSM_SOURCES ${PROJECT_SOURCE_DIR}/src/physics/sm_FSM*.F90)
		foreach(source ${SM_FSM_SOURCES})
			set_source_files_properties(${source} PROPERTIES COMPILE_OPTIONS "-mp")
		endforeach()
	endif()
endif()

# ============================================================================
# PHASE G — Finalize (banner, install, tooling)
# ============================================================================
############## Print User Information #################
set(MAIN_CFLAGS "${HICAR_Fortran_FLAGS_STR}")

# if FFTW_DIR is empty
if (NOT DEFINED FFTW_DIR OR FFTW_DIR STREQUAL "")
	set(FFTW_DIR ${FFTW_INCLUDES})
endif()

message(STATUS "")
message(STATUS "--------------------------------------------------------------------------------------------------------------------")
message(STATUS "")
message(STATUS "The fortran compiler ${CMAKE_Fortran_COMPILER} will be used to compile HICAR with the following flags:")
message(STATUS "")
message(STATUS "  ${CMAKE_Fortran_COMPILER} ${MAIN_CFLAGS}")
message(STATUS "")
message(STATUS "The following libraries are linked:")
message(STATUS "  NetCDF: ${NETCDF_LIBRARIES}")
message(STATUS "  FFTW:   ${FFTW_LIBRARIES}")
if (CUFFT_LIBRARY)
	message(STATUS "  cuFFT:  ${CUFFT_LIBRARY}")
endif()
message(STATUS "  MPI:    ${MPI_LIBRARIES}")
if (FSM)
	message(STATUS "  FSM:    ${FSM_DIR}")
else()
	message(STATUS "  FSM:    Not Used")
endif()
if (SNOWPACK_CPP)
	message(STATUS "  SNOWPACK:    C++ wrapper (${SNOWPACK_DIR})")
else()
	message(STATUS "  SNOWPACK:    native-Fortran port (${SNOWPACK_DIR})")
endif()
if (OPENACC AND OpenACC_Fortran_FOUND)
	message(STATUS "  OpenACC: Enabled")
else()
	message(STATUS "  OpenACC: Not Used")
endif()
message(STATUS "")
#######################################################

#This sets the location for library installation -- i.e., ${CMAKE_INSTALL_LIBDIR}=/usr/local/lib 
# and also the header for library if it has, in this case we do not have. ${CMAKE_INSTALL_INCLUDEDIR}=/usr/local/include


# Set the install prefix to the project source directory so that
# "make install" works without root permissions
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR} CACHE PATH "Install path prefix" FORCE)
endif()
set(CMAKE_INSTALL_EXEDIR ${PROJECT_SOURCE_DIR}/bin)
set(CMAKE_Fortran_MODULE_DIRECTORY  ${CMAKE_BINARY_DIR})

# this is to copy the library to the diseired directory
install(TARGETS HICAR  DESTINATION  ${CMAKE_INSTALL_EXEDIR})

ADD_CUSTOM_TARGET( distclean make clean
				   COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_BINARY_DIR}/CMakeFiles
				   COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_BINARY_DIR}/cmake_install.cmake ${CMAKE_BINARY_DIR}/install_manifest.txt ${CMAKE_BINARY_DIR}/Makefile ${CMAKE_BINARY_DIR}/CMakeCache.txt
					)

# Test / CI "make" targets (check, test_*, helpers) live under tests/.
include(${CMAKE_CURRENT_SOURCE_DIR}/tests/test_targets.cmake)
