# Always set the cmake min version.
cmake_minimum_required(VERSION 3.20)

project(HICAR)
set(CMAKE_CXX_STANDARD 14)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

enable_language (Fortran)

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
option(CAF "Turn on Coarray-Fortran (only available with Cray Compiler)" OFF)
option(FSM "Link to the FSM2 snow model, requires separate compilation of the FSM2 model" ON)
#########################################

#Set constants for project
set(BUILD ${PROJECT_SOURCE_DIR}/build)
add_subdirectory(${PROJECT_SOURCE_DIR}/tests)
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

#-DCRAY_PE -- controls usage of Coarray-Fortran over MPI
if(CAF AND COMPILER STREQUAL "cray")
  	add_definitions(-DCRAY_PE)
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


add_executable(HICAR ${SOURCE_FILES})
add_library(hicar_lib STATIC ${SOURCE_FILES})

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

# The installed petsc package on cray machines (piz daint) somehow does not
# automatically link to some dependencies of the petsc library, so we need to
# manually link to them. This could be done by compiling the whole model with the
# cray compiler flag --static, but this would also link to the static version of the
# other libraries, which is not what we want. So we manually link to the petsc library
# dependencies here.
message(STATUS "Looking for PETSc...")
list(APPEND CMAKE_MESSAGE_INDENT "  ")

#First, try to find Petsc with just find_package
if(PKG_CONFIG_FOUND)
	message(STATUS "1. Trying to find PETSc with pkg-config...")

	pkg_search_module(PETSC IMPORTED_TARGET PETSc)
	
	if(PETSC_FOUND AND NOT (DEFINED ENV{CRAYPE_VERSION} AND (PETSC MATCHES ".*cray.*")))
		message(STATUS "PETSc found with pkg-config")
		target_link_libraries(HICAR PUBLIC PkgConfig::PETSC)
		target_link_libraries(hicar_lib PUBLIC PkgConfig::PETSC)
	elseif (EXISTS "${PETSC}/lib/petsc/conf/petscvariables")
		message(STATUS "2. Trying to find PETSc by mimicking the cray-petsc makefile setup...")

		# try to find PETSc by mimicking the makefile setup
		# this is recommended on the PETSc site for makefiles:
		# https://petsc.org/main/manual/getting_started/#writing-c-c-or-fortran-applications
		# and the cray_PE currently does not robustly support using CMake for PETSc
		# issue: (https://lists.mcs.anl.gov/pipermail/petsc-users/2020-March/040491.html)
		file(READ "${PETSC}/lib/petsc/conf/petscvariables" MAKEFILE_CONTENTS)

		# Convert the MAKEFILE_CONTENTS to a list of lines
		string(REGEX REPLACE ";\\ ?" ";" MAKEFILE_CONTENTS "${MAKEFILE_CONTENTS}")
		string(REPLACE "\n" ";" MAKEFILE_CONTENTS "${MAKEFILE_CONTENTS}")
		list(REMOVE_DUPLICATES MAKEFILE_CONTENTS)
		
		# Initialize the PETSC_WITH_EXTERNAL_LIB variable
		set(PETSC_WITH_EXTERNAL_LIB "")
		
		# Iterate through the list of lines
		foreach(LINE ${MAKEFILE_CONTENTS})
			# Try to match the PETSC_WITH_EXTERNAL_LIB definition
			string(REGEX MATCH "PETSC_WITH_EXTERNAL_LIB = (.+)" MATCH_RESULT "${LINE}")
		
			# If a match is found, extract the value and break out of the loop
			if(MATCH_RESULT)
				string(REGEX REPLACE "PETSC_WITH_EXTERNAL_LIB = " "" PETSC_WITH_EXTERNAL_LIB "${MATCH_RESULT}")
				message(STATUS "PETSC_WITH_EXTERNAL_LIB: ${PETSC_WITH_EXTERNAL_LIB}")
				target_link_libraries(HICAR  PUBLIC ${PETSC_WITH_EXTERNAL_LIB})
				break()
			endif()
		endforeach()

		#If this doesn't work, and pkg-config was found, then just force linking to static libraries
		if(PETSC_WITH_EXTERNAL_LIB STREQUAL "" AND PKG_CONFIG_FOUND)
			message(WARNING "Failed to extract PETSC_WITH_EXTERNAL_LIB from Makefile")
			message(STATUS "3. Trying to find PETSc by linking to static libraries...")
			pkg_search_module(PETSC REQUIRED IMPORTED_TARGET PETSc)
			if(PETSC_FOUND)
				# Remove the petsc library from the list of libraries to link to
				# the ftn wrapper will automatically link to the correct petsc library
				list(REMOVE_ITEM PETSC_STATIC_LDFLAGS "-lpetsc")
				# Link to the petsc library static dependencies
				target_link_libraries(HICAR  PUBLIC ${PETSC_STATIC_LDFLAGS})
			else()
				message(FATAL_ERROR "PETSc not found")
			endif()
		endif()
	else()
                message(FATAL_ERROR "PETSc not found")
	endif()
endif()
list(POP_BACK CMAKE_MESSAGE_INDENT)

target_link_libraries(HICAR PUBLIC MPI::MPI_Fortran)
if(OpenMP_CXX_FOUND)
    target_link_libraries(HICAR PUBLIC OpenMP::OpenMP_CXX)
endif()
target_link_libraries(HICAR PUBLIC ${FFTW_LIBRARIES})

if (NETCDF_FOUND)
	#This is the most robust way to find the NetCDF flags -- covering both custom installs, static or shared, and cray PEs
	execute_process(COMMAND nf-config --flibs OUTPUT_VARIABLE NETCDF_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
	#if the nf-config method worked...
	if(${NETCDF_FLAGS} MATCHES "-lnetcdff")
		target_link_libraries(HICAR PUBLIC ${NETCDF_FLAGS})
	else()
	#else hope that the library flags as found by find_package are enough (should be as long as netCDF is a shared library)
		target_link_libraries(HICAR PUBLIC ${NETCDF_LIBRARIES})
	endif()
endif()


if (FSM)
	find_library(FSM2_LIB NAMES FSM2 REQUIRED)
	target_include_directories(HICAR PUBLIC ${FSM_DIR}/build)
	target_link_libraries(HICAR PUBLIC ${FSM2_LIB})

	target_include_directories(hicar_lib PUBLIC ${FSM_DIR}/build)

endif()

target_include_directories(HICAR PUBLIC ${MPI_INCLUDE_DIR})
target_include_directories(HICAR PUBLIC ${FFTW_INCLUDES})
target_include_directories(HICAR PUBLIC ${NETCDF_INCLUDES})

# Add the include directory to the compiler flags
target_include_directories(HICAR PUBLIC ${PROJECT_SOURCE_DIR}/include)

target_include_directories(hicar_lib PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories(hicar_lib PUBLIC ${MPI_INCLUDE_DIR})
target_include_directories(hicar_lib PUBLIC ${FFTW_INCLUDES})
target_include_directories(hicar_lib PUBLIC ${NETCDF_INCLUDES})

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
		target_compile_options(HICAR PUBLIC -h noomp -c -eI -G0 -m2 -rl -R bcdsp) #-R cpbs --cray-print-opts
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
		target_compile_options(HICAR PUBLIC -g -O3 -fopenmp -fopt-info-missed=opt_info.txt -mfma -c -ffast-math -ffree-line-length-none -floop-block -floop-strip-mine -ftree-vectorize -fimplicit-none -funroll-loops -floop-nest-optimize -march=native  -fno-protect-parens)
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
		target_compile_options(HICAR PUBLIC -O3 -fopenmp -finline-limit=10000 -mfma -c -ffast-math -ffree-line-length-none -floop-block -floop-strip-mine -ftree-vectorize -fimplicit-none -funroll-loops -floop-nest-optimize -march=native  -fno-protect-parens)  #-fcheck=all -Wall  # -ffast-math #-flto #
	elseif (COMPILER STREQUAL "intel")
		target_compile_options(HICAR PUBLIC -c -u -O3 -xHost -ftz -fpe0)# -check stack,bounds -fp-stack-check
	elseif (COMPILER STREQUAL "cray")
		target_compile_options(HICAR PUBLIC -h omp -hfp3 -O2 -c -eI)
		set(LINK "-h noomp")
	endif()
endif()

#set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FFLAGS}")
#set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LFLAGS}")

##########################################################


######################################################k#
############## Print User Information #################
get_target_property(MAIN_CFLAGS HICAR COMPILE_OPTIONS)
string(REGEX REPLACE ";" " " MAIN_CFLAGS "${MAIN_CFLAGS}")

message(STATUS "")
message(STATUS "--------------------------------------------------------------------------------------------------------------------")
message(STATUS "")
message(STATUS "The fortran compiler ${CMAKE_Fortran_COMPILER} will be used to compile HICAR with the following flags:")
message(STATUS "")
message(STATUS "  ${CMAKE_Fortran_COMPILER} ${MAIN_CFLAGS}")
message(STATUS "")
message(STATUS "Libraries were searched for in the following directories:")
message(STATUS "  NetCDF: ${NETCDF_DIR}")
message(STATUS "  PETSc:  ${PETSC}")
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

get_target_property(HICAR_LIBS HICAR LINK_LIBRARIES)
message(STATUS "Libraries linked against HICAR target: ${HICAR_LIBS}")

# Set the location for library installation -- i.e., /usr/lib in this case
# not really necessary in this example. Use "sudo make install" to apply 
set(CMAKE_INSTALL_EXEDIR ${PROJECT_SOURCE_DIR}/bin)
set(CMAKE_Fortran_MODULE_DIRECTORY  ${PROJECT_SOURCE_DIR}/build)

# this is to copy the library to the diseired directory
install(TARGETS HICAR  DESTINATION  ${CMAKE_INSTALL_EXEDIR})


# this is to copy build folder including .mod files in $HOME/include directory
# note that CMAKE_Fortran_MODULE_DIRECTORY is ouput path for modules or anything
#install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/build DESTINATION ${CMAKE_Fortran_MODULE_DIRECTORY}) 
#install(PROGRAMS ${modeFiles} DESTINATION ${CMAKE_Fortran_MODULE_DIRECTORY}/sec) # not working
    
#add_executable(test_FSM_interface test/test_FSM_interface.f90)
#target_link_libraries(test_FSM_interface FSM_interface ${PROJECT_LINK_LIBS})

ADD_CUSTOM_TARGET( distclean make clean
				   COMMAND ${CMAKE_COMMAND} -E remove_directory ${PROJECT_SOURCE_DIR}/build/CMakeFiles
				   COMMAND ${CMAKE_COMMAND} -E remove ${PROJECT_SOURCE_DIR}/build/cmake_install.cmake ${PROJECT_SOURCE_DIR}/build/install_manifest.txt ${PROJECT_SOURCE_DIR}/build/Makefile ${PROJECT_SOURCE_DIR}/build/CMakeCache.txt
					)

