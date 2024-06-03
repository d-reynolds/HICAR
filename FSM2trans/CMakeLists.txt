# Always set the cmake min version.
cmake_minimum_required(VERSION 3.0)

project(FSM_interface)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake")

enable_language (Fortran)

#########################################
####### Set up the user variables #######
set(FC "" CACHE STRING "Set the fortran compiler to use (can be auto-detected for most cases)")
#########################################

FILE(GLOB_RECURSE SOURCE_FILES  code/*.F90)

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
	set(CMAKE_Fortran_COMPILER "ftn")
	set(CMAKE_Fortran_COMPILER_LOADED TRUE)
	if (NOT DEFINED FC OR FC STREQUAL "")
		if (DEFINED $ENV{FC})
			message(STATUS "Using Fortran compiler supplied by the environment variable FC")
			set(COMPILER $ENV{FC})
		else()
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
	endif()
else()
	# If we are not on a cray system, we can use the compiler specified by the user
	if (NOT DEFINED FC OR FC STREQUAL "")
		if (DEFINED $ENV{FC})
			message(STATUS "Using Fortran compiler supplied by the environment variable FC: ${FC}")
			set(FC $ENV{FC})
		else()
			message(FATAL_ERROR "No Fortran compiler specified")
		endif()
	elseif(DEFINED FC)
		message(STATUS "Using Fortran compiler specified by the user in variable FC: ${FC}")
	endif()
	if (FC STREQUAL "gfortran")
		set(CMAKE_Fortran_COMPILER "gfortran")
		set(CMAKE_Fortran_COMPILER_ID "GNU")
		set(CMAKE_Fortran_COMPILER_LOADED TRUE)
		set(COMPILER "gnu")
	elseif (FC STREQUAL "ifort")
		set(CMAKE_Fortran_COMPILER "ifort")
		set(CMAKE_Fortran_COMPILER_ID "Intel")
		set(CMAKE_Fortran_COMPILER_LOADED TRUE)
		set(COMPILER "intel")
	else()
		message(FATAL_ERROR "Unknown Fortran compiler (FC) specified")
	endif()
endif()

string( TOLOWER "${MODE}" MODE )
string( TOLOWER "${COMPILER}" COMPILER )

list(POP_BACK CMAKE_MESSAGE_INDENT)
message(STATUS "...done auto-detecting the environment")

############################################

set(CMAKE_BUILD_TYPE Release)

if (COMPILER STREQUAL "gnu")
	set(CMAKE_Fortran_FLAGS "-fopenmp -lgomp -O2 -mfma -c -ffree-line-length-none -ftree-vectorize -fimplicit-none -funroll-loops -march=native  -fno-protect-parens")  #-fcheck=all -Wall  # -ffast-math #-flto #
elseif (COMPILER STREQUAL "intel")
	set(CMAKE_Fortran_FLAGS "-c -u -qopenmp -liomp5 -O3 -xHost -ftz -fpe0")# -check stack,bounds -fp-stack-check
elseif (COMPILER STREQUAL "cray")
	set(CMAKE_Fortran_FLAGS "-h noomp -hfp3 -O2 -c -eI")
endif()

#This sets the location for library installation -- i.e., ${CMAKE_INSTALL_LIBDIR}=/usr/local/lib 
# and also the header for library if it has, in this case we do not have. ${CMAKE_INSTALL_INCLUDEDIR}=/usr/local/include
include(GNUInstallDirs)

    
add_library(FSM_interface STATIC ${SOURCE_FILES})


# Set the location for library installation -- i.e., /usr/lib in this case
# not really necessary in this example. Use "sudo make install" to apply 
set(CMAKE_INSTALL_LIBDIR ${PROJECT_SOURCE_DIR}/lib/)
set(CMAKE_Fortran_MODULE_DIRECTORY  ${PROJECT_SOURCE_DIR}/build/)

# this is to copy the library to the diseired directory
install(TARGETS FSM_interface  DESTINATION  ${CMAKE_INSTALL_LIBDIR})


# this is to copy build folder including .mod files in $HOME/include directory
# note that CMAKE_Fortran_MODULE_DIRECTORY is ouput path for modules or anything
#install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/build DESTINATION ${CMAKE_Fortran_MODULE_DIRECTORY}) 
#install(PROGRAMS ${modeFiles} DESTINATION ${CMAKE_Fortran_MODULE_DIRECTORY}/sec) # not working
    
#add_executable(test_FSM_interface test/test_FSM_interface.f90)
#target_link_libraries(test_FSM_interface FSM_interface ${PROJECT_LINK_LIBS})

ADD_CUSTOM_TARGET( distclean make clean
				   COMMAND ${CMAKE_COMMAND} -E remove_directory CMakeFiles
				   COMMAND ${CMAKE_COMMAND} -E remove *~ cmake_install.cmake install_manifest.txt Makefile CMakeCache.txt
				   COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_Fortran_MODULE_DIRECTORY}
				   COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_INSTALL_LIBDIR}
					)
