# Always set the cmake min version.
cmake_minimum_required(VERSION 3.20)

project(HICAR)
set(CMAKE_CXX_STANDARD 14)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

enable_language (Fortran)

#Do not use any conda packages
if(DEFINED ENV{CONDA_PREFIX})
  set(CMAKE_IGNORE_PATH $ENV{CONDA_PREFIX} ${CMAKE_IGNORE_PATH})
  message(STATUS "Ignoring Conda environment path: $ENV{CONDA_PREFIX}")
endif()

#########################################
####### Set up the user variables #######
set(MODE "release" CACHE STRING "Set compilation flags, useful for debugging")
set(FC "" CACHE STRING "Set the fortran compiler to use (can be auto-detected for most cases)")

if (NOT DEFINED PETSC_DIR OR PETSC_DIR STREQUAL "")
	set(PETSC_DIR $ENV{PETSC_DIR} CACHE PATH "Path to PETSC installation")
endif()
if (NOT DEFINED PETSC_ARCH OR PETSC_ARCH STREQUAL "")
	set(PETSC_ARCH $ENV{PETSC_ARCH} CACHE PATH "System architecture for PETSC installation")
endif()
if (NOT DEFINED AMGX_DIR OR AMGX_DIR STREQUAL "")
	set(AMGX_DIR $ENV{AMGX_DIR} CACHE PATH "Path to NVIDIA AmgX installation")
endif()
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

if (NOT $ENV{NOAHMP_LIBRARY} STREQUAL "")
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
option(SNOWPACK "Link to the SNOWPACK snow model, requires separate compilation of the SNOWPACK model" OFF)
option(SNOWPACK_GPU "Use GPU-accelerated pure-Fortran SNOWPACK (replaces C++ SNOWPACK, requires OpenACC)" OFF)
option(NCCL "Link to the NCCL library, if available" ON)
option(OPENACC "Enable OpenACC support for GPU acceleration" OFF)
option(SRUN_FLAGS "srun flags to use for running the test cases" "")

# SNOWPACK_GPU implies both SNOWPACK and OPENACC
if(SNOWPACK_GPU)
    set(SNOWPACK ON CACHE BOOL "Enable SNOWPACK support in HICAR" FORCE)
    if(NOT OPENACC_SET_BY_USER)
        set(OPENACC ON CACHE BOOL "Enable OpenACC support for GPU acceleration" FORCE)
    endif()
endif()
#########################################

#Set constants for project
add_subdirectory(${PROJECT_SOURCE_DIR}/tests EXCLUDE_FROM_ALL)
#######################################


#########################################################
######## Set up compiler options for the build ##########

#######################################
########## Set up the sources #########
FILE(GLOB_RECURSE SOURCE_FILES  ${PROJECT_SOURCE_DIR}/src/*.F90)

# remove the two wind solvers, which will later be conditionally added back in
# depending on library support
list(REMOVE_ITEM SOURCE_FILES ${PROJECT_SOURCE_DIR}/src/physics/wind_iterative_petsc.F90)
list(REMOVE_ITEM SOURCE_FILES ${PROJECT_SOURCE_DIR}/src/physics/wind_iterative_amgx.F90)
FILE(GLOB_RECURSE ra_mo_files  ${PROJECT_SOURCE_DIR}/src/physics/ra_mo_*.F90)
list(REMOVE_ITEM SOURCE_FILES ${ra_mo_files})


macro(Disable_FSM)
	set(FSM OFF)
	FILE(GLOB_RECURSE FSM2_INTERFACE  ${PROJECT_SOURCE_DIR}/src/physics/FSM2_interface/*)
	list(REMOVE_ITEM SOURCE_FILES  ${PROJECT_SOURCE_DIR}/src/physics/sm_FSM.F90)
	list(REMOVE_ITEM SOURCE_FILES ${FSM2_INTERFACE})
endmacro()

#-DUSE_ASSERTIONS -- controls whether to use assertions in the code
if(ASSERTIONS)
  	add_definitions(-DUSE_ASSERTIONS=.true.)
else()
  	add_definitions(-DUSE_ASSERTIONS=.false.)
endif()

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

#-DVERSION -- sets the version of the code (used in output) to the current git version
set(GIT_VERSION $(shell git describe --long --dirty --all --always | sed -e's/heads\///'))
add_definitions(-DVERSION=\"$(GIT_VERSION)\")


############################################
############## Auto-detect  ################
message(STATUS "")
message(STATUS "Trying to auto-detect the environment...")
list(APPEND CMAKE_MESSAGE_INDENT "  ")

# Environment

# Detect if we are running on a cray system, where ftn is always use as the fortran compiler
# and the compiler type is determined by the compiler flags passed to ftn
if (DEFINED ENV{CRAYPE_VERSION})
	message(STATUS "Auto-detected Cray system, setting compiler to ftn")

	set(CMAKE_Fortran_COMPILER_ID "Cray")
	set(CMAKE_Fortran_COMPILER_LOADED TRUE)
	if (NOT DEFINED FC OR FC STREQUAL "")
		#if (DEFINED ENV{FC})
			#message(STATUS "Using Fortran compiler flags for the compiler supplied by the environment variable FC")
			#set(COMPILER $ENV{FC})
			#else()
			if ("$ENV{PE_ENV}" STREQUAL "GNU")
				message(STATUS "Auto-detected GNU Programming Environment, using gnu compiler flags")
				set(COMPILER "gnu")
			elseif("$ENV{PE_ENV}" STREQUAL "INTEL")
				message(STATUS "Auto-detected Intel Programming Environment, using intel compiler flags")
				set(COMPILER "intel")
			elseif ("$ENV{PE_ENV}" STREQUAL "CRAY")
				message(STATUS "Auto-detected Cray Programming Environment, using cray compiler flags")
				set(COMPILER "cray")
			else()
				message(FATAL_ERROR "Unknown Programming Environment")
			endif()
	endif()
	set(FC "ftn")
else()
	# If we are not on a cray system, we can use the compiler specified by the user
	if (NOT DEFINED FC OR FC STREQUAL "")
		if (DEFINED ENV{FC})
			message(STATUS "Using Fortran compiler supplied by the environment variable FC: $ENV{FC}")
			set(FC $ENV{FC})
		elseif(DEFINED CMAKE_Fortran_COMPILER)
			message(STATUS "Using Auto-detected Fortran compiler: ${CMAKE_Fortran_COMPILER}")
			set(FC ${CMAKE_Fortran_COMPILER})
		else()
			message(FATAL_ERROR "No Fortran compiler specified")
		endif()
	elseif(DEFINED FC)
		message(STATUS "Using Fortran compiler specified by the user in variable FC: ${FC}")
	endif()
	if (FC MATCHES "gfortran")
		set(CMAKE_Fortran_COMPILER "${FC}")
		set(CMAKE_Fortran_COMPILER_ID "GNU")
		set(CMAKE_Fortran_COMPILER_LOADED TRUE)
		set(COMPILER "gnu")
	elseif (FC MATCHES "mpif90")
		#set(CMAKE_Fortran_COMPILER "mpif90")
		#set(CMAKE_Fortran_COMPILER_ID "mpif90")
		#set(CMAKE_Fortran_COMPILER_LOADED TRUE)
		set(COMPILER "gnu")
        elseif (FC MATCHES "mpifort")
                #set(CMAKE_Fortran_COMPILER "mpifort")
                #set(CMAKE_Fortran_COMPILER_ID "mpifort")
                #set(CMAKE_Fortran_COMPILER_LOADED TRUE)
                set(COMPILER "gnu")
	elseif (FC MATCHES "nvfortran")
		set(CMAKE_Fortran_COMPILER "${FC}")
		set(CMAKE_Fortran_COMPILER_ID "NVHPC")
		set(CMAKE_Fortran_COMPILER_LOADED TRUE)
		set(COMPILER "nvfortran")
		# Auto-enable OpenACC for NVFortran if not explicitly set
		if (NOT OPENACC_SET_BY_USER)
			set(OPENACC ON CACHE BOOL "Enable OpenACC support for GPU acceleration" FORCE)
			message(STATUS "Auto-enabling OpenACC support for NVFortran compiler")
		endif()

	else()
		message(FATAL_ERROR "Unknown Fortran compiler (FC) specified")
	endif()
endif()

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
if(MODE STREQUAL "debug")
	set(CMAKE_BUILD_TYPE Debug)

	if (COMPILER STREQUAL "gnu")
		set(HICAR_Fortran_FLAGS -O0 -g -fopenmp -fbounds-check -fcheck=all -Wall -fbacktrace -finit-real=nan -ffree-line-length-none -Wconversion)
	elseif (COMPILER STREQUAL "intel")
		set(HICAR_Fortran_FLAGS -debug -O2 -u -check all -check noarg_temp_created -traceback -fpe0 -fast-transcendentals -xhost)
	elseif (COMPILER STREQUAL "cray")
		set(HICAR_Fortran_FLAGS -h noomp -eI -g -O0)
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
	elseif (COMPILER STREQUAL "intel")
		set(HICAR_Fortran_FLAGS -debug -debug-parameters all -traceback -ftrapuv -g -fpe0 -u -check all -check noarg_temp_created -CB)
	elseif (COMPILER STREQUAL "cray")
		set(HICAR_Fortran_FLAGS -h noomp -eI -m 2 -O0 -g -R bcs)
	endif()
elseif(MODE STREQUAL "profile")
	set(CMAKE_BUILD_TYPE Profile)

	if (COMPILER STREQUAL "gnu")
		set(HICAR_Fortran_FLAGS -g -O3 -fopenmp -fbounds-check -fbacktrace -finit-real=nan -fopt-info-missed=opt_info.txt -ffast-math -ffree-line-length-none -ftree-vectorize -fimplicit-none -funroll-loops -march=native -fno-protect-parens)
	elseif (COMPILER STREQUAL "intel")
		set(HICAR_Fortran_FLAGS -pg -debug inline-debug-info -shared-intel)
	elseif (COMPILER STREQUAL "cray")
		set(HICAR_Fortran_FLAGS -h noomp -eI -hfp3 -O2 -g -R bcs)
	endif()
elseif(MODE STREQUAL "release")
	set(CMAKE_BUILD_TYPE Release)
	if (COMPILER STREQUAL "gnu")
		set(HICAR_Fortran_FLAGS -O3 -fopenmp -ffast-math -ffree-line-length-none -ftree-vectorize -fimplicit-none -funroll-loops -march=native -fno-protect-parens)
	elseif (COMPILER STREQUAL "intel")
		set(HICAR_Fortran_FLAGS -u -O3 -xHost -ftz -fpe0)
	elseif (COMPILER STREQUAL "cray")
		set(HICAR_Fortran_FLAGS -h omp -hfp3 -O2 -e0 -eI -eZ)
	elseif (COMPILER STREQUAL "nvfortran")
		if(OPENACC)
			set(HICAR_Fortran_FLAGS -O3 -Minfo=accel -march=native -acc=gpu -gpu=mem:separate,ccnative)
		else()
			set(HICAR_Fortran_FLAGS -O3 -mp -march=native)
		endif()
	endif()
endif()

# Convert list to space-separated string for display and external projects
string(REPLACE ";" " " HICAR_Fortran_FLAGS_STR "${HICAR_Fortran_FLAGS}")
message(STATUS "Fortran compiler flags: ${HICAR_Fortran_FLAGS_STR}")
message(STATUS "...done setting up compiler flags")
###########################################################




###########################################################
############### Find or install dependencies ##############
message(STATUS "")
message(STATUS "Searching for dependencies...")
list(APPEND CMAKE_MESSAGE_INDENT "  ")

include(GNUInstallDirs)

# set root of location to find PETSc's pkg-config
set(PETSC "${PETSC_DIR}/${PETSC_ARCH}")
set(ENV{PKG_CONFIG_PATH} "$ENV{PKG_CONFIG_PATH}:${PETSC}/lib/pkgconfig")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake" ${PETSC} ${AMGX_DIR} ${FFTW_DIR} ${NETCDF_DIR}  ${MPI_DIR} ${FSM_DIR} ${SNOWPACK_DIR})
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} "${PETSC}/lib/pkgconfig/" ${AMGX_DIR} ${FFTW_DIR} ${NETCDF_DIR} ${MPI_DIR} ${MPI_DIR} ${FSM_DIR} ${SNOWPACK_DIR})

find_package(PkgConfig REQUIRED)
include(FetchContent)

# 1. Check if not already included via Cray Module and ftn
# 2. If not included via Cray module, link to the package libraries manually

# 3. If not found, try to find package using find library
# 4. If not found, try to find package using find package
# 5. If not on machine, try to fetch and install package
# 5. Enforce requirement, check that package was found/installed

# 6. Link to the package libraries
set(NETCDF_F90 "YES")
find_package(MPI REQUIRED)
find_package(OpenMP)
find_package(FFTW REQUIRED)
find_package(NetCDF REQUIRED)
if (SNOWPACK_GPU)
	# GPU SNOWPACK: compile pure-Fortran GPU modules directly, no C++ SNOWPACK needed
	if(NOT SNOWPACK_DIR OR SNOWPACK_DIR STREQUAL "")
		find_package(SNOWPACK REQUIRED)
	endif()
	set(SNOWPACK_GPU_SRC_DIR "${SNOWPACK_DIR}/fortran")
	file(GLOB SNOWPACK_GPU_SOURCES "${SNOWPACK_GPU_SRC_DIR}/snowpack_gpu_*.F90")
	if(NOT SNOWPACK_GPU_SOURCES)
		message(FATAL_ERROR "No snowpack_gpu_*.F90 files found in ${SNOWPACK_GPU_SRC_DIR}")
	endif()
	list(APPEND SOURCE_FILES ${SNOWPACK_GPU_SOURCES})
	# Remove the C++ SNOWPACK Fortran driver (replaced by snowpack_gpu_driver)
	list(REMOVE_ITEM SOURCE_FILES ${PROJECT_SOURCE_DIR}/src/physics/sm_SNOWPACK.F90)
	message(STATUS "SNOWPACK GPU sources: ${SNOWPACK_GPU_SOURCES}")
elseif (SNOWPACK)
	find_package(SNOWPACK REQUIRED)
else()
	# remove sm_SNOWPACK.F90, which contains the SNOWPACK driver, from the source files if SNOWPACK support is disabled
	FILE(GLOB_RECURSE SNOWPACK_INTERFACE  ${PROJECT_SOURCE_DIR}/src/physics/sm_SNOWPACK.F90)
	list(REMOVE_ITEM SOURCE_FILES ${SNOWPACK_INTERFACE})
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
		
		# Fetch and build rte-rrtmgp when OpenACC is enabled
		message(STATUS "Fetching rte-rrtmgp library for OpenACC support...")

		# Configure build options for rte-rrtmgp
		set(BUILD_TESTING OFF)
		set(BUILD_C_HEADERS OFF)
		set(KERNEL_MODE "accel" CACHE STRING "Select the kernel mode for rte-rrtmgp" FORCE)
		set(BUILD_SHARED_LIBS OFF)

		FetchContent_MakeAvailable(rte-rrtmgp)
			
		# Use the predefined HICAR compiler flags			
		message(STATUS "Using compiler flags for rte-rrtmgp: ${CMAKE_Fortran_FLAGS}")

		if(NOT TARGET rte-rrtmgp::rte AND NOT TARGET rte-rrtmgp::rrtmgp)
			# Add the rte-rrtmgp subdirectory with configured options
			add_subdirectory(${rte-rrtmgp_SOURCE_DIR} ${rte-rrtmgp_BINARY_DIR} EXCLUDE_FROM_ALL)
		endif()
		
		# Add preprocessor definition to indicate rte-rrtmgp is available
		add_definitions(-DUSE_RTE_RRTMGP)
		list(APPEND SOURCE_FILES ${ra_mo_files})

		message(STATUS "rte-rrtmgp library configured successfully")


		# Link AmgX library if available
		if(DEFINED AMGX_DIR AND NOT AMGX_DIR STREQUAL "")
			find_library(AMGX_LIBRARY NAMES amgx PATHS ${AMGX_DIR}/lib ${AMGX_DIR}/lib64 NO_DEFAULT_PATH)
			if(AMGX_LIBRARY)
				message(STATUS "Found AmgX library: ${AMGX_LIBRARY}")
				# Store the actual include directory path
				set(AMGX_INCLUDE_DIR "${AMGX_DIR}/include" CACHE PATH "Path to AmgX include directory")
			else()
				message(WARNING "AmgX library not found in ${AMGX_DIR}. Variational wind solver will not be available for this GPU build.")
			endif()
		else()
			# Try to auto-detect AmgX from NVIDIA HPC SDK
			if(DEFINED ENV{NVHPC_ROOT})
				set(NVHPC_AMGX_PATH "$ENV{NVHPC_ROOT}/math_libs")
				find_library(AMGX_LIBRARY NAMES amgx 
					PATHS 
						${NVHPC_AMGX_PATH}/lib64
						${NVHPC_AMGX_PATH}/lib
					NO_DEFAULT_PATH)
				if(AMGX_LIBRARY)
					message(STATUS "Auto-detected AmgX from NVIDIA HPC SDK: ${AMGX_LIBRARY}")
					set(AMGX_DIR ${NVHPC_AMGX_PATH} CACHE PATH "Path to NVIDIA AmgX installation (auto-detected)" FORCE)
					set(AMGX_INCLUDE_DIR "${NVHPC_AMGX_PATH}/include" CACHE PATH "Path to AmgX include directory")
				else()
					message(WARNING "AMGX_DIR not set and auto-detection failed. Variational wind solver will not be available.")
					message(WARNING "  To enable AmgX, set AMGX_DIR or ensure NVHPC_ROOT is set correctly.")
				endif()
			else()
				message(WARNING "AMGX_DIR not set. Variational wind solver will not be available for this GPU build. Set AMGX_DIR to enable.")
			endif()
		endif()
		
		# If AmgX was found, add source files and set preprocessor definition
		if(AMGX_LIBRARY)
			add_definitions(-DUSE_AMGX)
			list(APPEND SOURCE_FILES ${PROJECT_SOURCE_DIR}/src/physics/wind_iterative_amgx.F90)
			list(APPEND SOURCE_FILES ${PROJECT_SOURCE_DIR}/src/physics/mpi_comm_wrapper.c)
			message(STATUS "AmgX support enabled - added wind_iterative_amgx.F90 and mpi_comm_wrapper.c")
		endif()
	endif()
endif()
#########################################################

if(NOT AMGX_LIBRARY)
	find_package(PETSc)
	if (PETSC_FOUND)
		if (${PETSC_VERSION} VERSION_GREATER_EQUAL 3.23)
			message(FATAL_ERROR "PETSc version ${PETSC_VERSION} is greater than or equal to 3.23, code currently only supports versions less than 3.23.")
		endif()
		add_definitions(-DUSE_PETSC)
		list(APPEND SOURCE_FILES ${PROJECT_SOURCE_DIR}/src/physics/wind_iterative_petsc.F90)
	endif()
endif()
#########################################################

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


###########################################################


list(REMOVE_ITEM SOURCE_FILES ${PROJECT_SOURCE_DIR}/src/main/driver.F90)
add_library(hicar_lib STATIC ${SOURCE_FILES})
add_executable(HICAR ${PROJECT_SOURCE_DIR}/src/main/driver.F90)

set_property(TARGET hicar_lib PROPERTY EXCLUDE_FROM_DEFAULT_BUILD FALSE)
set_property(TARGET HICAR PROPERTY EXCLUDE_FROM_DEFAULT_BUILD FALSE)


if(AMGX_LIBRARY)
	# Enable C, C++, and CUDA languages for AmgX support
	enable_language(C CXX CUDA)
	
	# Find CUDA libraries required by AmgX
	find_package(CUDAToolkit REQUIRED)
	
	# Link AmgX library, C++ standard library, and CUDA libraries (all required for AmgX)
	# AmgX requires: cudart (runtime), cublas (linear algebra), cusparse (sparse matrices), cusolver (dense linear solvers)
	target_link_libraries(hicar_lib PUBLIC ${AMGX_LIBRARY} stdc++ CUDA::cudart CUDA::cublas CUDA::cusparse CUDA::cusolver)
	target_include_directories(hicar_lib PUBLIC ${AMGX_INCLUDE_DIR})
	
	message(STATUS "Linked AmgX and CUDA libraries to hicar_lib")
elseif(PETSC_FOUND)
	target_link_libraries(hicar_lib PUBLIC ${PETSC_LIBRARIES})
	target_include_directories(hicar_lib PUBLIC ${PETSC_INCLUDES})
endif()

if(OpenMP_CXX_FOUND)
        target_link_libraries(hicar_lib PUBLIC OpenMP::OpenMP_CXX)
elseif(OpenMP_Fortran_FOUND)
	target_link_libraries(hicar_lib PUBLIC OpenMP::OpenMP_Fortran)
endif()

if (FFTW_FOUND)
	target_link_libraries(hicar_lib PUBLIC ${FFTW_LIBRARIES})
	target_include_directories(hicar_lib PUBLIC ${FFTW_INCLUDES})
endif()

if (NETCDF_FOUND)
    # Check for parallel NetCDF support
    execute_process(COMMAND nc-config --has-parallel
                    OUTPUT_VARIABLE NETCDF_HAS_PARALLEL
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                    ERROR_QUIET)
    
    if(NETCDF_HAS_PARALLEL STREQUAL "yes")
        message(STATUS "NetCDF built with parallel support")
        add_definitions(-DNETCDF_PARALLEL)
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

		# Clean up any duplicate libraries and references to -lnetcdf
		list(REMOVE_DUPLICATES NETCDF_LIBRARIES)
		# clean references to base netcdf-c library
		# list(REMOVE_ITEM NETCDF_LIBRARIES -lnetcdf)

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
	add_definitions(-DFSM)
endif()

if(SNOWPACK_GPU)
    # GPU SNOWPACK: pure Fortran, no C++ libraries needed
    add_definitions(-DSNOWPACK -DSNOWPACK_GPU)
    message(STATUS "SNOWPACK GPU support enabled (pure Fortran OpenACC, no C++ SNOWPACK/MeteoIO needed)")
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
    
	add_definitions(-DSNOWPACK)
    message(STATUS "SNOWPACK support enabled")
endif()
if (OpenACC_Fortran_FOUND AND OPENACC)
	target_link_libraries(hicar_lib PUBLIC OpenACC::OpenACC_Fortran)
	
	# Link rte-rrtmgp libraries when OpenACC is enabled
	if(TARGET rte-rrtmgp::rte AND TARGET rte-rrtmgp::rrtmgp)
		target_link_libraries(hicar_lib PUBLIC rte-rrtmgp::rte rte-rrtmgp::rrtmgp)
		message(STATUS "Linked rte-rrtmgp libraries to HICAR")
	endif()

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
		add_definitions(-DUSE_NCCL)
		# NCCL requires CUDA runtime (cudaSetDevice, cudaFree, cudaStream*, etc.)
		# This may already be available if AMGX is linked, but we need it regardless.
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

target_link_libraries(hicar_lib PUBLIC MPI::MPI_Fortran)
target_include_directories(hicar_lib PUBLIC ${MPI_INCLUDE_DIR})

# add include directory for the model itself
target_include_directories(hicar_lib PUBLIC ${PROJECT_SOURCE_DIR}/include)

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

if(MODE STREQUAL "debug" OR MODE STREQUAL "debugslow")
	set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "HICAR_debug")
	if (COMPILER STREQUAL "nvfortran" AND OPENACC)
		set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "HICAR_debug_gpu")
	endif()
elseif(MODE STREQUAL "profile")
	set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "HICAR_profile")
elseif(MODE STREQUAL "release")
	set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "HICAR")
	if (COMPILER STREQUAL "cray")
		set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "HICAR_cray")
	elseif (COMPILER STREQUAL "nvfortran" AND OPENACC)
		set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "HICAR_gpu")
	endif()
endif()

# Use the predefined HICAR compiler flags (already a list)
target_compile_options(HICAR PUBLIC ${HICAR_Fortran_FLAGS})


if (COMPILER STREQUAL "cray" AND FSM)
	#Set compiler flags for FSM_interface files to be the same as those for target HICAR, with the addition of -e0
	# Set compile options for each file
	file(GLOB FORTRAN_SOURCES ${PROJECT_SOURCE_DIR}/src/physics/FSM2_interface/*.F90)

	foreach(source ${FORTRAN_SOURCES})
		set_property(SOURCE ${source} APPEND PROPERTY COMPILE_OPTIONS -e0)
	endforeach()
endif()

# sm_* files (snow model drivers) need OpenMP but not OpenACC (for CPU snow models).
# With SNOWPACK_GPU, sm_driver.F90 calls GPU kernels and needs OpenACC.
if(COMPILER STREQUAL "nvfortran" AND OPENACC)
	if(NOT SNOWPACK_GPU)
		# CPU snow models: disable OpenACC on all sm_* files
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
		# GPU SNOWPACK: sm_driver.F90 needs OpenACC (calls snowpack_gpu_step with !$acc present),
		# but sm_FSM.F90 still needs OpenMP-only
		file(GLOB SM_FSM_SOURCES ${PROJECT_SOURCE_DIR}/src/physics/sm_FSM*.F90)
		foreach(source ${SM_FSM_SOURCES})
			set_source_files_properties(${source} PROPERTIES COMPILE_OPTIONS "-mp")
		endforeach()
	endif()
endif()

# Apply the same compile options to hicar_lib as are set on HICAR
get_target_property(COMPILER_OPTIONS HICAR COMPILE_OPTIONS)
if(COMPILER_OPTIONS)
	target_compile_options(hicar_lib PUBLIC ${COMPILER_OPTIONS})
endif()

#set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LFLAGS}")

##########################################################


###########################################################
############### Download supporting run files ##############

# Download the repo but don't build it
FetchContent_MakeAvailable(icar_supporting_files)

# Copy only specific folders to the run directory
file(COPY "${icar_supporting_files_SOURCE_DIR}/mp_support" DESTINATION "${PROJECT_SOURCE_DIR}/run")
file(COPY "${icar_supporting_files_SOURCE_DIR}/rrtmg_support" DESTINATION "${PROJECT_SOURCE_DIR}/run")

# Clean up downloaded repo to save space
file(REMOVE_RECURSE "${icar_supporting_files_SOURCE_DIR}")
file(REMOVE_RECURSE "${CMAKE_BINARY_DIR}/_deps/icar_supporting_files-subbuild")
file(REMOVE_RECURSE "${CMAKE_BINARY_DIR}/_deps/icar_supporting_files-build")

# if we are using RTE_RRTMGP, copy supporting files from the git repository: https://github.com/earth-system-radiation/rrtmgp-data
if(OPENACC AND OpenACC_Fortran_FOUND)
	message(STATUS "Fetching rte-rrtmgp data files...")
	
	# Download the repo but don't build it
	FetchContent_MakeAvailable(rrtmgp_data)
	
	# Copy rrtmgp data files to the run directory
	file(COPY "${rrtmgp_data_SOURCE_DIR}/" DESTINATION "${PROJECT_SOURCE_DIR}/run/rrtmgp_support")
		
	# Clean up downloaded repo to save space
	file(REMOVE_RECURSE "${rrtmgp_data_SOURCE_DIR}")
	file(REMOVE_RECURSE "${CMAKE_BINARY_DIR}/_deps/rrtmgp_data-subbuild")
	file(REMOVE_RECURSE "${CMAKE_BINARY_DIR}/_deps/rrtmgp_data-build")
	
	message(STATUS "RTE-RRTMGP data files copied to run/rrtmgp_support")
endif()


list(POP_BACK CMAKE_MESSAGE_INDENT)
message(STATUS "...done searching for dependencies")

###########################################################






######################################################k#
############## Print User Information #################
get_target_property(MAIN_CFLAGS hicar_lib COMPILE_OPTIONS)
string(REGEX REPLACE ";" " " MAIN_CFLAGS "${MAIN_CFLAGS}")

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
message(STATUS "  MPI:    ${MPI_LIBRARIES}")
if(DEFINED AMGX_DIR AND NOT AMGX_DIR STREQUAL "" AND AMGX_LIBRARY)
	message(STATUS "  AmgX:   ${AMGX_LIBRARY}")
elseif(DEFINED PETSC_DIR AND NOT PETSC_DIR STREQUAL "" AND PETSC_LIBRARIES)
	message(STATUS "  PETSc:  ${PETSC_LIBRARIES}")
endif()
if (FSM)
	message(STATUS "  FSM:    ${FSM_DIR}")
else()
	message(STATUS "  FSM:    Not Used")
endif()
if (SNOWPACK)
	message(STATUS "  SNOWPACK:    ${SNOWPACK_DIR}")
else()
	message(STATUS "  SNOWPACK:    Not Used")
endif()
if (OPENACC AND OpenACC_Fortran_FOUND)
	message(STATUS "  OpenACC: Enabled")
	if(TARGET rte-rrtmgp::rte AND TARGET rte-rrtmgp::rrtmgp)
		message(STATUS "  RTE-RRTMGP: Enabled (linked for OpenACC support)")
	endif()
else()
	message(STATUS "  OpenACC: Not Used")
	message(STATUS "  RTE-RRTMGP: Not Used")
endif()
message(STATUS "")
#######################################################

#This sets the location for library installation -- i.e., ${CMAKE_INSTALL_LIBDIR}=/usr/local/lib 
# and also the header for library if it has, in this case we do not have. ${CMAKE_INSTALL_INCLUDEDIR}=/usr/local/include


# Create a variable for the include directory
#set(INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)

# Add the include directory to the compiler flags
#target_compile_options(HICAR PRIVATE -I${INCLUDE_DIR})

# get_target_property(hicar_lib HICAR LINK_LIBRARIES)
# message(STATUS "Libraries linked against HICAR target: ${HICAR_LIBS}")

# Set the install prefix to the project source directory so that
# "make install" works without root permissions
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR} CACHE PATH "Install path prefix" FORCE)
endif()
set(CMAKE_INSTALL_EXEDIR ${PROJECT_SOURCE_DIR}/bin)
set(CMAKE_Fortran_MODULE_DIRECTORY  ${CMAKE_BINARY_DIR})

# this is to copy the library to the diseired directory
install(TARGETS HICAR  DESTINATION  ${CMAKE_INSTALL_EXEDIR})


# this is to copy build folder including .mod files in $HOME/include directory
# note that CMAKE_Fortran_MODULE_DIRECTORY is ouput path for modules or anything
#install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/build DESTINATION ${CMAKE_Fortran_MODULE_DIRECTORY}) 
#install(PROGRAMS ${modeFiles} DESTINATION ${CMAKE_Fortran_MODULE_DIRECTORY}/sec) # not working
    

ADD_CUSTOM_TARGET( distclean make clean
				   COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_BINARY_DIR}/CMakeFiles
				   COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_BINARY_DIR}/cmake_install.cmake ${CMAKE_BINARY_DIR}/install_manifest.txt ${CMAKE_BINARY_DIR}/Makefile ${CMAKE_BINARY_DIR}/CMakeCache.txt
					)
# Add a custom target for running tests with MPI

separate_arguments(SRUN_FLAGS_LIST UNIX_COMMAND "${SRUN_FLAGS}")
SET(SRUN_FLAGS_LIST_TEST_CASE "-SRUN_FLAGS='${SRUN_FLAGS}'")

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
	COMMAND ${CMAKE_COMMAND} -E echo "Running test cases..."
	COMMAND "${PROJECT_SOURCE_DIR}/tests/Test_Cases/test_case_runner.sh" ${PROJECT_SOURCE_DIR} "Standard,Standard_nested" ${SRUN_FLAGS_LIST_TEST_CASE} 
	DEPENDS download_test_data
	WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests/Test_Cases
	COMMENT "Running test case suite"
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
	COMMAND ${CMAKE_COMMAND} -E echo "Running test cases..."
	COMMAND ${PROJECT_SOURCE_DIR}/tests/Test_Cases/test_case_runner.sh ${PROJECT_SOURCE_DIR} Standard,Standard_restart
	DEPENDS download_test_data
	WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests/Test_Cases
	COMMENT "Running test case suite"
	VERBATIM
	)

endif()

add_custom_target(test_cli
    COMMAND ${CMAKE_COMMAND} -E echo "Running CLI option tests..."
    COMMAND ${PROJECT_SOURCE_DIR}/tests/test_cli_options.sh ${PROJECT_SOURCE_DIR}
    DEPENDS HICAR
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
            -DPETSC_DIR=${PETSC_DIR}
            -DPETSC_ARCH=${PETSC_ARCH}
            -DAMGX_DIR=${AMGX_DIR}
            -DFFTW_DIR=${FFTW_DIR}
            -DNETCDF_DIR=${NETCDF_DIR}
            -DMPI_DIR=${MPI_DIR}
            -DFSM_DIR=${FSM_DIR}
            -DFSM=${FSM}
			-DSNOWPACK_DIR=${SNOWPACK_DIR}
			-DSNOWPACK=${SNOWPACK}
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

add_custom_target(test_regression
    COMMAND ${PROJECT_SOURCE_DIR}/tests/test_regression.sh
        ${PROJECT_SOURCE_DIR}
        ${CMAKE_BINARY_DIR}
    DEPENDS HICAR download_test_data
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests/Test_Cases
    COMMENT "Running regression test against known-good commit"
    VERBATIM
)
