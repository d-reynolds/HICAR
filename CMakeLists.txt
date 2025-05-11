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


option(ASSERTIONS "Check logical assertions at runtime" ON)
option(FSM "Link to the FSM2 snow model, requires separate compilation of the FSM2 model" ON)
#########################################

#Set constants for project
add_subdirectory(${PROJECT_SOURCE_DIR}/tests EXCLUDE_FROM_ALL)
#######################################


#########################################################
######## Set up compiler options for the build ##########

#######################################
########## Set up the sources #########
FILE(GLOB_RECURSE SOURCE_FILES  ${PROJECT_SOURCE_DIR}/src/*.F90)

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

#########################################################


list(REMOVE_ITEM SOURCE_FILES ${PROJECT_SOURCE_DIR}/src/main/driver.F90)
add_library(hicar_lib STATIC ${SOURCE_FILES})
add_executable(HICAR ${PROJECT_SOURCE_DIR}/src/main/driver.F90)

set_property(TARGET hicar_lib PROPERTY EXCLUDE_FROM_DEFAULT_BUILD FALSE)
set_property(TARGET HICAR PROPERTY EXCLUDE_FROM_DEFAULT_BUILD FALSE)

###########################################################
############### Find or install dependencies ##############
message(STATUS "")
message(STATUS "Searching for dependencies...")
list(APPEND CMAKE_MESSAGE_INDENT "  ")

include(GNUInstallDirs)

# set root of location to find PETSc's pkg-config
set(PETSC "${PETSC_DIR}/${PETSC_ARCH}")
set(ENV{PKG_CONFIG_PATH} "$ENV{PKG_CONFIG_PATH}:${PETSC}/lib/pkgconfig")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake" ${PETSC} ${FFTW_DIR} ${NETCDF_DIR}  ${MPI_DIR} ${FSM_DIR})
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} "${PETSC}/lib/pkgconfig/" ${FFTW_DIR} ${NETCDF_DIR} ${MPI_DIR} ${MPI_DIR} ${FSM_DIR})

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
find_package(PETSc REQUIRED)

if (PETSC_FOUND)
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

target_link_libraries(hicar_lib PUBLIC MPI::MPI_Fortran)
target_include_directories(hicar_lib PUBLIC ${MPI_INCLUDE_DIR})

# add include directory for the model itself
target_include_directories(hicar_lib PUBLIC ${PROJECT_SOURCE_DIR}/include)


# link the hicar_lib library to the HICAR executable
target_link_libraries(HICAR PUBLIC hicar_lib)
###########################################################


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

list(POP_BACK CMAKE_MESSAGE_INDENT)
message(STATUS "...done searching for dependencies")

###########################################################


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
	else()
		message(FATAL_ERROR "Unknown Fortran compiler (FC) specified")
	endif()
endif()

string( TOLOWER "${MODE}" MODE )
string( TOLOWER "${COMPILER}" COMPILER )

list(POP_BACK CMAKE_MESSAGE_INDENT)
message(STATUS "...done auto-detecting the environment")

################################################


##########################################################
############## Set up the compiler options ###############

if(MODE STREQUAL "debug")
	set(CMAKE_BUILD_TYPE Debug)
	set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "HICAR_debug")

    if (COMPILER STREQUAL "gnu")
		target_compile_options(HICAR PUBLIC -c -O0 -g -fopenmp -fbounds-check -fbacktrace -finit-real=nan -ffree-line-length-none -Wconversion)
	elseif (COMPILER STREQUAL "intel")
		target_compile_options(HICAR PUBLIC -debug -c -O2 -u -check all -check noarg_temp_created -traceback -fpe0 -fast-transcendentals -xhost)
	elseif (COMPILER STREQUAL "cray")
		target_compile_options(HICAR PUBLIC -h noomp -c -eI -g -O0) #-R cpbs --cray-print-opts
	endif()
elseif(MODE STREQUAL "debugslow")
	set(CMAKE_BUILD_TYPE Debug)

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
		target_compile_options(HICAR PUBLIC -g -O3 -fopenmp -fbounds-check -fbacktrace -fopt-info-missed=opt_info.txt -c -ffast-math -ffree-line-length-none -floop-block -floop-strip-mine -ftree-vectorize -fimplicit-none -funroll-loops -floop-nest-optimize -march=native  -fno-protect-parens)
		target_link_options(HICAR PUBLIC -pg)
	elseif (COMPILER STREQUAL "intel")
		target_compile_options(HICAR PUBLIC -pg -debug inline-debug-info -shared-intel)
	elseif (COMPILER STREQUAL "cray")
		target_compile_options(HICAR PUBLIC -c -h noomp -eI -m 2 -O0 -hpl=my_program_debug.pl -g -R bcs)
	endif()
elseif(MODE STREQUAL "release")
	set(CMAKE_BUILD_TYPE Release)
	set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "HICAR")

	if (COMPILER STREQUAL "gnu")
		target_compile_options(HICAR PUBLIC -O3 -fopenmp -c -ffast-math -ffree-line-length-none -floop-block -floop-strip-mine -ftree-vectorize -fimplicit-none -funroll-loops -floop-nest-optimize -march=native  -fno-protect-parens)  #-fcheck=all -Wall  # -ffast-math #-flto #
	elseif (COMPILER STREQUAL "intel")
		target_compile_options(HICAR PUBLIC -c -u -O3 -xHost -ftz -fpe0)# -check stack,bounds -fp-stack-check
	elseif (COMPILER STREQUAL "cray")
		target_compile_options(HICAR PUBLIC -h omp -hfp3 -O2 -c -eI)
		set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "HICAR_cray")

		# set(LINK "-h noomp")
	endif()
endif()

# Apply the same compile options to hicar_lib as are set on HICAR
get_target_property(COMPILER_OPTIONS HICAR COMPILE_OPTIONS)
if(COMPILER_OPTIONS)
	target_compile_options(hicar_lib PUBLIC ${COMPILER_OPTIONS})
endif()

#set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FFLAGS}")
#set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LFLAGS}")

##########################################################


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
message(STATUS "Libraries were searched for in the following directories:")
message(STATUS "  NetCDF: ${NETCDF_DIR}")
message(STATUS "  PETSc:  ${PETSC_DIR}/${PETSC_ARCH}")
message(STATUS "  FFTW:   ${FFTW_DIR}")
message(STATUS "  MPI:    ${MPI_DIR}")
if (FSM)
	message(STATUS "  FSM:    ${FSM_DIR}")
else()
	message(STATUS "  FSM:    Not Used")
endif()
message(STATUS "")
#######################################################

#This sets the location for library installation -- i.e., ${CMAKE_INSTALL_LIBDIR}=/usr/local/lib 
# and also the header for library if it has, in this case we do not have. ${CMAKE_INSTALL_INCLUDEDIR}=/usr/local/include


# Create a variable for the include directory
#set(INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)

# Add the include directory to the compiler flags
#target_compile_options(HICAR PRIVATE -I${INCLUDE_DIR})

get_target_property(hicar_lib HICAR LINK_LIBRARIES)
message(STATUS "Libraries linked against HICAR target: ${HICAR_LIBS}")

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
		message(ERROR "No account specified for srun, please pass your SLURM account name to cmake using -SRUN_FLAGS=<SRUN_FLAGS>")
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
	DEPENDS ${PROJECT_NAME}-tester 
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
	DEPENDS ${PROJECT_NAME}-tester 
	WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests/Test_Cases
	COMMENT "Running test case suite"
	VERBATIM
	)  

endif()
