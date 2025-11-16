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
option(FSM "Link to the FSM2 snow model, requires separate compilation of the FSM2 model" ON)
option(OPENACC "Enable OpenACC support for GPU acceleration" OFF)
option(SRUN_FLAGS "srun flags to use for running the test cases" "")
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


macro(Disable_FSM)
	set(FSM OFF)
	list(REMOVE_ITEM SOURCE_FILES ${PROJECT_SOURCE_DIR}/src/physics/lsm_FSMdrv.F90)
	FILE(GLOB_RECURSE FSM2_INTERFACE  ${PROJECT_SOURCE_DIR}/src/physics/FSM2_interface/*)
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
	add_definitions(-DFSM)
elseif(FSM AND NOT (IS_DIRECTORY ${FSM_DIR}/lib AND
					IS_DIRECTORY ${FSM_DIR}/build))
	message(WARNING "FSM2 not found, disabling FSM2 support.\nTo enable, set the FSM_DIR variable to the root directory of the FSM2 install (i.e. the directory which contains the 'lib' and 'build' subdirectories).")
	Disable_FSM()
elseif(NOT IS_DIRECTORY FSM_DIR AND NOT FSM)
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
		set(CMAKE_Fortran_COMPILER "gfortran")
		set(CMAKE_Fortran_COMPILER_ID "GNU")
		set(CMAKE_Fortran_COMPILER_LOADED TRUE)
		set(COMPILER "gnu")
	elseif (FC MATCHES "ifort")
		set(CMAKE_Fortran_COMPILER "ifort")
		set(CMAKE_Fortran_COMPILER_ID "Intel")
		set(CMAKE_Fortran_COMPILER_LOADED TRUE)
		set(COMPILER "intel")
	elseif (FC MATCHES "mpif90")
		set(CMAKE_Fortran_COMPILER "mpif90")
		set(CMAKE_Fortran_COMPILER_ID "mpif90")
		set(CMAKE_Fortran_COMPILER_LOADED TRUE)
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
############### Find or install dependencies ##############
message(STATUS "")
message(STATUS "Searching for dependencies...")
list(APPEND CMAKE_MESSAGE_INDENT "  ")

include(GNUInstallDirs)

# set root of location to find PETSc's pkg-config
set(PETSC "${PETSC_DIR}/${PETSC_ARCH}")
set(ENV{PKG_CONFIG_PATH} "$ENV{PKG_CONFIG_PATH}:${PETSC}/lib/pkgconfig")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake" ${PETSC} ${AMGX_DIR} ${FFTW_DIR} ${NETCDF_DIR}  ${MPI_DIR} ${FSM_DIR})
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} "${PETSC}/lib/pkgconfig/" ${AMGX_DIR} ${FFTW_DIR} ${NETCDF_DIR} ${MPI_DIR} ${MPI_DIR} ${FSM_DIR})

find_package(PkgConfig REQUIRED)
include(FetchContent)

# 1. Check if not already included via Cray Module and ftn
# 2. If not included via Cray module, link to the package libraries manually

# 3. If not found, try to find package using find library
# 4. If not found, try to find package using find package
# 5. If not on machine, try to fetch and install package
# 5. Enforce requirement, check that package was found/installed

# 6. Link to the package libraries
set (NETCDF_F90 "YES")
find_package(MPI REQUIRED)
find_package(OpenMP)
find_package(FFTW REQUIRED)
find_package(NetCDF REQUIRED)

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
		
		# Get the content without building it yet
		FetchContent_Declare(
			rte-rrtmgp
			GIT_REPOSITORY    https://github.com/earth-system-radiation/rte-rrtmgp.git
			GIT_TAG           v1.9.2
			GIT_SHALLOW       TRUE
		)

		FetchContent_GetProperties(rte-rrtmgp)
		if(NOT rte-rrtmgp_POPULATED)
			FetchContent_Populate(rte-rrtmgp)

			# Configure build options for rte-rrtmgp
			set(BUILD_TESTING OFF)
			set(BUILD_C_HEADERS OFF)
			set(KERNEL_MODE "accel" CACHE STRING "Select the kernel mode for rte-rrtmgp" FORCE)
			set(BUILD_SHARED_LIBS OFF)
			
			# Set compiler flags for rte-rrtmgp to match HICAR's flags
			# Save current CMAKE_Fortran_FLAGS
			set(SAVED_CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
			
			# Determine and set appropriate compiler flags based on MODE and COMPILER
			if(MODE STREQUAL "debug")
				if (COMPILER STREQUAL "gnu")
					set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -c -O0 -g -fopenmp -fbounds-check -fcheck=all -Wall -fbacktrace -finit-real=nan -ffree-line-length-none -Wconversion")
				elseif (COMPILER STREQUAL "intel")
					set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -debug -c -O2 -u -check all -check noarg_temp_created -traceback -fpe0 -fast-transcendentals -xhost")
				elseif (COMPILER STREQUAL "cray")
					set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -h noomp -c -eI -g -O0")
				elseif (COMPILER STREQUAL "nvfortran")
					if(OPENACC)
						set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O0 -g -traceback -Mnoopenmp -Mbounds -Minfo=accel -Mpreprocess -march=native -acc=gpu -gpu=mem:separate,ccnative")
					else()
						set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O0 -g -traceback -Mnoopenmp -Mbounds -Mpreprocess -march=native")
					endif()
				endif()
			elseif(MODE STREQUAL "release")
				if (COMPILER STREQUAL "gnu")
					set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -fopenmp -c -ffast-math -ffree-line-length-none -floop-block -floop-strip-mine -ftree-vectorize -fimplicit-none -funroll-loops -floop-nest-optimize -march=native -fno-protect-parens")
				elseif (COMPILER STREQUAL "intel")
					set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -c -u -O3 -xHost -ftz -fpe0")
				elseif (COMPILER STREQUAL "cray")
					set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -h omp -hfp3 -O2 -c -eI -eZ")
				elseif (COMPILER STREQUAL "nvfortran")
					if(OPENACC)
						set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -Minfo=accel -march=native -acc=gpu -gpu=mem:separate,ccnative")
					else()
						set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -march=native")
					endif()
				endif()
			endif()
			
			message(STATUS "Using compiler flags for rte-rrtmgp: ${CMAKE_Fortran_FLAGS}")

			# Add the rte-rrtmgp subdirectory with configured options
			add_subdirectory(${rte-rrtmgp_SOURCE_DIR} ${rte-rrtmgp_BINARY_DIR} EXCLUDE_FROM_ALL)
			
			# Restore original CMAKE_Fortran_FLAGS
			set(CMAKE_Fortran_FLAGS "${SAVED_CMAKE_Fortran_FLAGS}")
		endif()
		
		# Add preprocessor definition to indicate rte-rrtmgp is available
		add_definitions(-DUSE_RTE_RRTMGP)
		
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
else()
	find_package(PETSc)
	if (PETSC_FOUND)
		if (${PETSC_VERSION} VERSION_GREATER_EQUAL 3.23)
			message(FATAL_ERROR "PETSc version ${PETSC_VERSION} is greater than or equal to 3.23, code currently only supports versions less than 3.23.")
		endif()
		add_definitions(-DUSE_PETSC)
		list(APPEND SOURCE_FILES ${PROJECT_SOURCE_DIR}/src/physics/wind_iterative_petsc.F90)
		target_link_libraries(hicar_lib PUBLIC ${PETSC_LIBRARIES})
		target_include_directories(hicar_lib PUBLIC ${PETSC_INCLUDES})
	endif()
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
	#This is the most robust way to find the NetCDF flags -- covering both custom installs, static or shared, and cray PEs
	execute_process(COMMAND nf-config --flibs OUTPUT_VARIABLE NETCDF_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
	#if the nf-config method worked...
	if(${NETCDF_FLAGS} MATCHES "-lnetcdff")
		target_link_libraries(hicar_lib PUBLIC ${NETCDF_FLAGS})
	else()
	#else hope that the library flags as found by find_package are enough (should be as long as netCDF is a shared library)
		target_link_libraries(hicar_lib PUBLIC ${NETCDF_LIBRARIES})
	endif()
	target_include_directories(hicar_lib PUBLIC ${NETCDF_INCLUDES})
endif()


if (FSM)
	find_library(FSM2_LIB NAMES FSM2 REQUIRED)
	target_include_directories(hicar_lib PUBLIC ${FSM_DIR}/build)
	target_link_libraries(hicar_lib PUBLIC ${FSM2_LIB})
endif()

if (OpenACC_Fortran_FOUND AND OPENACC)
	target_link_libraries(hicar_lib PUBLIC OpenACC::OpenACC_Fortran)
	
	# Link rte-rrtmgp libraries when OpenACC is enabled
	if(TARGET rte-rrtmgp::rte AND TARGET rte-rrtmgp::rrtmgp)
		target_link_libraries(hicar_lib PUBLIC rte-rrtmgp::rte rte-rrtmgp::rrtmgp)
		message(STATUS "Linked rte-rrtmgp libraries to HICAR")
	endif()
endif()

target_link_libraries(hicar_lib PUBLIC MPI::MPI_Fortran)
target_include_directories(hicar_lib PUBLIC ${MPI_INCLUDE_DIR})

# add include directory for the model itself
target_include_directories(hicar_lib PUBLIC ${PROJECT_SOURCE_DIR}/include)


# link the hicar_lib library to the HICAR executable
target_link_libraries(HICAR PUBLIC hicar_lib)
###########################################################

##########################################################
############## Set up the compiler options ###############
# Print default flags for different build types
message(STATUS "CMAKE_Fortran_FLAGS_DEBUG: ${CMAKE_Fortran_FLAGS_DEBUG}")
message(STATUS "CMAKE_Fortran_FLAGS_RELEASE: ${CMAKE_Fortran_FLAGS_RELEASE}")
message(STATUS "CMAKE_Fortran_FLAGS_RELWITHDEBINFO: ${CMAKE_Fortran_FLAGS_RELWITHDEBINFO}")
message(STATUS "CMAKE_Fortran_FLAGS_MINSIZEREL: ${CMAKE_Fortran_FLAGS_MINSIZEREL}")
message(STATUS "CMAKE_Fortran_FLAGS: ${CMAKE_Fortran_FLAGS}")

if(MODE STREQUAL "debug")
	set(CMAKE_BUILD_TYPE Debug)
	set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "HICAR_debug")

    if (COMPILER STREQUAL "gnu")
	    target_compile_options(HICAR PUBLIC -c -O0 -g -fopenmp -fbounds-check -fcheck=all -Wall -fbacktrace -finit-real=nan -ffree-line-length-none -Wconversion)
	elseif (COMPILER STREQUAL "intel")
		target_compile_options(HICAR PUBLIC -debug -c -O2 -u -check all -check noarg_temp_created -traceback -fpe0 -fast-transcendentals -xhost)
	elseif (COMPILER STREQUAL "cray")
		target_compile_options(HICAR PUBLIC -h noomp -c -eI -g -O0) #-R cpbs --cray-print-opts
	elseif (COMPILER STREQUAL "nvfortran")
		if(OPENACC)
			target_compile_options(HICAR PUBLIC -O0 -g -traceback -Mnoopenmp -Mbounds -Minfo=accel -Mpreprocess -march=native -acc=gpu,multicore -gpu=mem:separate,ccnative)
		else()
			target_compile_options(HICAR PUBLIC -O0 -g -traceback -Mnoopenmp -Mbounds -Mpreprocess -march=native)
		endif()
	endif()
elseif(MODE STREQUAL "debugslow")
	set(CMAKE_BUILD_TYPE Debug)
	set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "HICAR_debug")

	if (COMPILER STREQUAL "gnu")
		target_compile_options(HICAR PUBLIC -c -g -ffree-line-length-none -Wconversion)
	elseif (COMPILER STREQUAL "intel")
		target_compile_options(HICAR PUBLIC -debug -debug-parameters all -traceback -ftrapuv -g -fpe0 -c -u -check all -check noarg_temp_created -CB)
	elseif (COMPILER STREQUAL "cray")
		target_compile_options(HICAR PUBLIC -c -h noomp -eI -m 2 -O0 -g -R bcs)
	endif()
elseif(MODE STREQUAL "profile")
	set(CMAKE_BUILD_TYPE Profile)

	if (COMPILER STREQUAL "gnu")
		target_compile_options(HICAR PUBLIC -g -O3 -fopenmp -fbounds-check -fbacktrace -finit-real=nan -fopt-info-missed=opt_info.txt -c -ffast-math -ffree-line-length-none -floop-block -floop-strip-mine -ftree-vectorize -fimplicit-none -funroll-loops -floop-nest-optimize -march=native  -fno-protect-parens)
		target_link_options(HICAR PUBLIC -pg)
	elseif (COMPILER STREQUAL "intel")
		target_compile_options(HICAR PUBLIC -pg -debug inline-debug-info -shared-intel)
	elseif (COMPILER STREQUAL "cray")
		target_compile_options(HICAR PUBLIC -c -h noomp -eI -hfp3 -O2 -g -R bcs)
	endif()
elseif(MODE STREQUAL "release")
	set(CMAKE_BUILD_TYPE Release)
	set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "HICAR")

	if (COMPILER STREQUAL "gnu")
		target_compile_options(HICAR PUBLIC -O3 -fopenmp -c -ffast-math -ffree-line-length-none -floop-block -floop-strip-mine -ftree-vectorize -fimplicit-none -funroll-loops -floop-nest-optimize -march=native  -fno-protect-parens)  #-fcheck=all -Wall  # -ffast-math #-flto #
	elseif (COMPILER STREQUAL "intel")
		target_compile_options(HICAR PUBLIC -c -u -O3 -xHost -ftz -fpe0)# -check stack,bounds -fp-stack-check
	elseif (COMPILER STREQUAL "cray")
		target_compile_options(HICAR PUBLIC -h omp -hfp3 -O2 -c -eI -eZ)
		set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "HICAR_cray")
	elseif (COMPILER STREQUAL "nvfortran")
		if(OPENACC)
			set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "HICAR_gpu")
			target_compile_options(HICAR PUBLIC -O3 -Mnoopenmp -Minfo=accel -march=native -acc=gpu,multicore -gpu=mem:separate,ccnative)
		else()
			target_compile_options(HICAR PUBLIC -O3 -Mnoopenmp -march=native)
		endif()
	endif()
endif()

if (COMPILER STREQUAL "cray" AND FSM)
	#Set compiler flags for FSM_interface files to be the same as those for target HICAR, with the addition of -e0
	# Set compile options for each file
	file(GLOB FORTRAN_SOURCES ${PROJECT_SOURCE_DIR}/src/physics/FSM2_interface/*.F90)

	foreach(source ${FORTRAN_SOURCES})
		set_property(SOURCE ${source} APPEND PROPERTY COMPILE_OPTIONS -e0)
	endforeach()
endif()

# Apply the same compile options to hicar_lib as are set on HICAR
get_target_property(COMPILER_OPTIONS HICAR COMPILE_OPTIONS)
if(COMPILER_OPTIONS)
	target_compile_options(hicar_lib PUBLIC ${COMPILER_OPTIONS})
endif()

#set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FFLAGS}")
#set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LFLAGS}")

##########################################################


###########################################################
############### Download supporting run files ##############

# Fetch specific folders from icar_supporting_files repo
FetchContent_Declare(
	icar_supporting_files
	GIT_REPOSITORY    https://github.com/NCAR/icar_supporting_files.git
	GIT_TAG           main
)

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
	
	FetchContent_Declare(
		rrtmgp_data
		GIT_REPOSITORY    https://github.com/earth-system-radiation/rrtmgp-data.git
		GIT_TAG           main
		GIT_SHALLOW       TRUE
	)
	
	# Download the repo but don't build it
	FetchContent_MakeAvailable(rrtmgp_data)
	
	# Copy rrtmgp data files to the run directory
	file(COPY "${rrtmgp_data_SOURCE_DIR}" DESTINATION "${PROJECT_SOURCE_DIR}/run/rrtmgp_support")
		
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

# Set the location for library installation -- i.e., /usr/lib in this case
# not really necessary in this example. Use "sudo make install" to apply 
set(CMAKE_INSTALL_EXEDIR ${PROJECT_SOURCE_DIR}/bin)
set(CMAKE_Fortran_MODULE_DIRECTORY  ${CMAKE_BINARY_DIR})

# this is to copy the library to the diseired directory
install(TARGETS HICAR  DESTINATION  ${CMAKE_INSTALL_EXEDIR})


# this is to copy build folder including .mod files in $HOME/include directory
# note that CMAKE_Fortran_MODULE_DIRECTORY is ouput path for modules or anything
#install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/build DESTINATION ${CMAKE_Fortran_MODULE_DIRECTORY}) 
#install(PROGRAMS ${modeFiles} DESTINATION ${CMAKE_Fortran_MODULE_DIRECTORY}/sec) # not working
    
#add_executable(test_FSM_interface test/test_FSM_interface.f90)
#target_link_libraries(test_FSM_interface FSM_interface ${PROJECT_LINK_LIBS})

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
