#[=[
        FindSNOWPACK.cmake - CMake module to locate or build the SNOWPACK library

        This module defines the following variables:
            SNOWPACK_FOUND        - True if SNOWPACK was found or built
            SNOWPACK_INCLUDE_DIR - Include directories for SNOWPACK headers
            SNOWPACK_LIBRARY      - Library to link against for SNOWPACK
            SNOWPACK_FORTRAN_LIBRARY - Fortran bindings library for SNOWPACK
            METEOIO_INCLUDE_DIR - Include directories for MeteoIO headers
            METEOIO_LIBRARY      - Library to link against for MeteoIO
            BUILD_SNOWPACK_FROM_SOURCE - True if SNOWPACK was built from source

        Behavior:
            1. First checks if required variables are already set in the cache
            2. If not, look for SNOWPACK and meteoIO in SNOWPACK_DIR and METEOIO_DIR
                2a. If SNOWPACK_DIR and METEOIO_DIR are not set, default to build directories under CMAKE_BINARY_DIR/external
            3. If not found, builds MeteoIO and SNOWPACK from source using ExternalProject

        Usage:
            find_package(SNOWPACK [REQUIRED] [QUIET])
            if(SNOWPACK_FOUND)
                target_include_directories(mytarget PRIVATE ${SNOWPACK_INCLUDE_DIR})
                target_link_libraries(mytarget ${SNOWPACK_LIBRARY})
            endif()
#]=]

macro(FIND_LIB lib_name var_name hint_path)
    #First try to find libs at exactly the path specified by hint_path
    find_library (${var_name} NAMES ${lib_name} HINTS "${hint_path}/lib" NO_DEFAULT_PATH)
    #...if that didn't work, then search more generally
    if (NOT ${var_name})
      find_library (${var_name} NAMES ${lib_name})
    endif()
endmacro(FIND_LIB)

macro(FIND_HEADER_PATH header_name var_name hint_path)
    #First try to find headers at exactly the path specified by hint_path
    find_path (${var_name} NAMES ${header_name} HINTS "${hint_path}/include" NO_DEFAULT_PATH)
    #...if that didn't work, then search more generally
    if (NOT ${var_name})
      find_path (${var_name} NAMES ${header_name})
    endif()
endmacro(FIND_HEADER_PATH)

include(ExternalProject)

if (SNOWPACK_INCLUDE_DIR AND SNOWPACK_LIBRARY)
  # Already in cache, be silent
  set (SNOWPACK_FIND_QUIETLY TRUE)
else()
    message(STATUS "")
    message(STATUS "--------------------------------------------------------------------------------------------------------------------")
endif (SNOWPACK_INCLUDE_DIR AND SNOWPACK_LIBRARY)

if (NOT SNOWPACK_DIR)
    set(SNOWPACK_DIR "${CMAKE_BINARY_DIR}/external/SNOWPACK" CACHE PATH "Directory where SNOWPACK is installed" FORCE)
endif()

if (NOT METEOIO_DIR)
    set(METEOIO_DIR "${CMAKE_BINARY_DIR}/external/meteoio" CACHE PATH "Directory where MeteoIO is installed" FORCE)
endif()


#First try to find includes at exactly the path specified by SNOWPACK_DIR and METEOIO_DIR
FIND_HEADER_PATH(snowpack_mod.mod SNOWPACK_INCLUDE_DIR "${SNOWPACK_DIR}")
FIND_HEADER_PATH(meteoio.h METEOIO_INCLUDE_DIR "${METEOIO_DIR}")

#First try to find libs at exactly the path specified by SNOWPACK_DIR and METEOIO_DIR
FIND_LIB(snowpack SNOWPACK_LIBRARY "${SNOWPACK_DIR}")
FIND_LIB(snowpack_fortran SNOWPACK_FORTRAN_LIBRARY "${SNOWPACK_DIR}")
FIND_LIB(meteoio METEOIO_LIBRARY "${METEOIO_DIR}")

if (SNOWPACK_GPU) 
    # If GPU SNOWPACK is enabled, also look for the GPU Fortran bindings library
    message(STATUS "SNOWPACK_GPU is enabled, only fetching snowpack repo...")

    set(SNOWPACK_DIR "${CMAKE_BINARY_DIR}/external/SNOWPACK" CACHE PATH "Directory where SNOWPACK is installed" FORCE)
    set(SNOWPACK_BUILD "${CMAKE_BINARY_DIR}/external/SNOWPACK-build" CACHE PATH "Directory where SNOWPACK is installed" FORCE)
    set(SNOWPACK_STAMPS "${CMAKE_BINARY_DIR}/external/SNOWPACK-stamps" CACHE PATH "Directory where SNOWPACK is installed" FORCE)

    # Fetch SNOWPACK
    FetchContent_Declare(SNOWPACK
        GIT_REPOSITORY    https://git.wsl.ch/snow-models/snowpack.git
        GIT_TAG           fortran-bindings
        GIT_SHALLOW       TRUE
        PREFIX            "${SNOWPACK_STAMPS}"
        SOURCE_DIR        "${SNOWPACK_DIR}"
        BINARY_DIR        "${SNOWPACK_BUILD}"
    )
	FetchContent_Populate(SNOWPACK)

    set(SNOWPACK_FOUND TRUE)

elseif (SNOWPACK_LIBRARY AND SNOWPACK_INCLUDE_DIR AND SNOWPACK_FORTRAN_LIBRARY AND METEOIO_LIBRARY AND METEOIO_INCLUDE_DIR)
    message(STATUS "Found pre-installed SNOWPACK library: ${SNOWPACK_LIBRARY}")
    message(STATUS "Using SNOWPACK include directory: ${SNOWPACK_INCLUDE_DIR}")
    set(SNOWPACK_FOUND TRUE)
    set(SNOWPACK TRUE CACHE BOOL "Enable SNOWPACK support in HICAR")
    # Mark as built from source so the main project can add dependencies
    set(BUILD_SNOWPACK_FROM_SOURCE FALSE CACHE BOOL "SNOWPACK should be built from source" FORCE)
else()
    # If still not found, but user wants it, build from source using ExternalProject
    message(STATUS "SNOWPACK not found. Building MeteoIO and SNOWPACK from source...")
    
    # Determine the correct C++ standard library for the platform
    if(APPLE)
        set(CXX_STDLIB "c++")
        # On macOS, the Makefile-based Fortran bindings build needs the SDK sysroot
        # so that Apple Clang can find C++ standard library headers (e.g. <iostream>)
        execute_process(
            COMMAND xcrun --show-sdk-path
            OUTPUT_VARIABLE MACOS_SDK_PATH
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        if(MACOS_SDK_PATH)
            set(SNOWPACK_SYSROOT_FLAG "-isysroot ${MACOS_SDK_PATH}")
        endif()
    else()
        set(CXX_STDLIB "stdc++")
    endif()

    # Determine the correct Fortran free-form flag for SNOWPACK bindings
    # Shroud-generated .f files are free-form despite the .f extension
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
        set(SNOWPACK_FORTRAN_FREEFORM_FLAG "-ffree-form -ffree-line-length-none")
    elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
        set(SNOWPACK_FORTRAN_FREEFORM_FLAG "-Mfreeform")
    elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" OR CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
        set(SNOWPACK_FORTRAN_FREEFORM_FLAG "-free")
    elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
        set(SNOWPACK_FORTRAN_FREEFORM_FLAG "-f free")
    else()
        message(WARNING "Unknown Fortran compiler '${CMAKE_Fortran_COMPILER_ID}' for SNOWPACK bindings, assuming GNU-style free-form flag")
        set(SNOWPACK_FORTRAN_FREEFORM_FLAG "-ffree-form")
    endif()

    set(SNOWPACK_FORTRAN_FFLAGS "${SNOWPACK_FORTRAN_FREEFORM_FLAG} -O3 -fPIC")

    # SNOWPACK's SetCompilerOptions.cmake has no NVHPC/nvc++ support.
    # Without this, OPTIM is empty and CMAKE_CXX_FLAGS_RELEASE is force-set to "",
    # causing the entire SNOWPACK C++ library to compile with no optimization.
    # USER_COMPILER_OPTIONS is added to EXTRA → CMAKE_CXX_FLAGS, bypassing the empty OPTIM.
    if(CMAKE_CXX_COMPILER_ID STREQUAL "NVHPC")
        set(SNOWPACK_USER_COMPILER_OPTIONS "-O3 -DNDEBUG -DNOSAFECHECKS")
    else()
        set(SNOWPACK_USER_COMPILER_OPTIONS "")
    endif()

    set(METEOIO_DIR "${CMAKE_BINARY_DIR}/external/meteoio" CACHE PATH "MeteoIO source directory" FORCE)
    set(METEOIO_BUILD "${CMAKE_BINARY_DIR}/external/METEOIO-build" CACHE PATH "MeteoIO build directory" FORCE)
    set(METEOIO_STAMPS "${CMAKE_BINARY_DIR}/external/METEOIO-stamps" CACHE PATH "MeteoIO build directory" FORCE)

    # Build MeteoIO first (SNOWPACK dependency)
    ExternalProject_Add(MeteoIO
        GIT_REPOSITORY    https://git.wsl.ch/snow-models/meteoio.git
        GIT_TAG           master
        GIT_SHALLOW       TRUE
        UPDATE_DISCONNECTED TRUE
        PREFIX            "${METEOIO_STAMPS}"
        SOURCE_DIR        "${METEOIO_DIR}"
        BINARY_DIR        "${METEOIO_BUILD}"
        CMAKE_ARGS
            -DCMAKE_INSTALL_PREFIX=${METEOIO_DIR}
            -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
            -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
            -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
            -DUSER_COMPILER_OPTIONS=${SNOWPACK_USER_COMPILER_OPTIONS}
            -DCMAKE_LIBRARY_OUTPUT_DIRECTORY=${METEOIO_DIR}/lib
            -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY=${METEOIO_DIR}/lib
            -DCMAKE_RUNTIME_OUTPUT_DIRECTORY=${METEOIO_DIR}/lib
        BUILD_COMMAND     ${CMAKE_COMMAND} --build . --parallel
        INSTALL_COMMAND   ${CMAKE_COMMAND} --install . --prefix ${METEOIO_DIR} --component headers
            COMMAND ${CMAKE_COMMAND} -E make_directory ${METEOIO_DIR}/lib
            COMMAND sh -c "cp -P ${METEOIO_BUILD}/meteoio/${CMAKE_SHARED_LIBRARY_PREFIX}meteoio${CMAKE_SHARED_LIBRARY_SUFFIX}* ${METEOIO_DIR}/lib/ 2>/dev/null || true"
            COMMAND sh -c "cp -P ${METEOIO_BUILD}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}meteoio${CMAKE_SHARED_LIBRARY_SUFFIX}* ${METEOIO_DIR}/lib/ 2>/dev/null || true"
        LOG_DOWNLOAD      ON
        LOG_CONFIGURE     ON
        LOG_BUILD         ON
        LOG_INSTALL       ON
    )
    
    set(METEOIO_LIBRARY ${METEOIO_DIR}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}meteoio${CMAKE_SHARED_LIBRARY_SUFFIX} CACHE FILEPATH "MeteoIO library" FORCE) 
    set(METEOIO_INCLUDE_DIR ${METEOIO_DIR}/include CACHE PATH "MeteoIO include directory" FORCE)


    set(SNOWPACK_DIR "${CMAKE_BINARY_DIR}/external/SNOWPACK" CACHE PATH "Directory where SNOWPACK is installed" FORCE)
    set(SNOWPACK_BUILD "${CMAKE_BINARY_DIR}/external/SNOWPACK-build" CACHE PATH "Directory where SNOWPACK is installed" FORCE)
    set(SNOWPACK_STAMPS "${CMAKE_BINARY_DIR}/external/SNOWPACK-stamps" CACHE PATH "Directory where SNOWPACK is installed" FORCE)

    set(SNOWPACK_FORTRAN_CXXFLAGS "${SNOWPACK_SYSROOT_FLAG} -O3 -fPIC -I${SNOWPACK_DIR} -I${METEOIO_DIR}")

    # Build SNOWPACK (depends on MeteoIO)
    ExternalProject_Add(SNOWPACK
        DEPENDS           MeteoIO
        GIT_REPOSITORY    https://git.wsl.ch/snow-models/snowpack.git
        GIT_TAG           fortran-bindings
        GIT_SHALLOW       TRUE
        UPDATE_DISCONNECTED TRUE
        PREFIX            "${SNOWPACK_STAMPS}"
        SOURCE_DIR        "${SNOWPACK_DIR}"
        BINARY_DIR        "${SNOWPACK_BUILD}"
        CMAKE_ARGS
            -DCMAKE_INSTALL_PREFIX=${SNOWPACK_DIR}
            -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
            -DMETEOIO_INCLUDE_DIR=${METEOIO_INCLUDE_DIR}
            -DMETEOIO_LIBRARY=${METEOIO_LIBRARY}
            -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
            -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
            -DUSER_COMPILER_OPTIONS=${SNOWPACK_USER_COMPILER_OPTIONS}
            -DCMAKE_LIBRARY_OUTPUT_DIRECTORY=${SNOWPACK_DIR}/lib
            -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY=${SNOWPACK_DIR}/lib
            -DCMAKE_RUNTIME_OUTPUT_DIRECTORY=${SNOWPACK_DIR}/lib
        BUILD_COMMAND     ${CMAKE_COMMAND} --build . --parallel
        INSTALL_COMMAND   ${CMAKE_COMMAND} --install . --prefix ${SNOWPACK_DIR} --component headers
            COMMAND ${CMAKE_COMMAND} -E make_directory ${SNOWPACK_DIR}/lib
            COMMAND sh -c "cp -P ${SNOWPACK_BUILD}/snowpack/${CMAKE_SHARED_LIBRARY_PREFIX}snowpack${CMAKE_SHARED_LIBRARY_SUFFIX}* ${SNOWPACK_DIR}/lib/ 2>/dev/null || true"
            COMMAND sh -c "cp -P ${SNOWPACK_BUILD}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}snowpack${CMAKE_SHARED_LIBRARY_SUFFIX}* ${SNOWPACK_DIR}/lib/ 2>/dev/null || true"
            COMMAND sh -c "python -m pip install llnl-shroud 2>/dev/null || pip install llnl-shroud"
            COMMAND ${CMAKE_COMMAND} -E chdir ${SNOWPACK_DIR}/fortran make
                "FC=${CMAKE_Fortran_COMPILER}"
                "CXX=${CMAKE_CXX_COMPILER}"
                "FFLAGS=${SNOWPACK_FORTRAN_FFLAGS}"
                "CXXFLAGS=${SNOWPACK_FORTRAN_CXXFLAGS}"
                "LDFLAGS=-L../lib -lsnowpack -L../../meteoio/lib -lmeteoio -l${CXX_STDLIB} -Wl,-rpath,../lib -Wl,-rpath,../../meteoio/lib"
            COMMAND ${CMAKE_COMMAND} -E chdir ${SNOWPACK_DIR}/fortran make install
        LOG_DOWNLOAD      ON
        LOG_CONFIGURE     ON
        LOG_BUILD         ON
        LOG_INSTALL       ON
    )
    
    # Set the paths where SNOWPACK will be installed
    set(SNOWPACK_LIBRARY "${SNOWPACK_DIR}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}snowpack${CMAKE_SHARED_LIBRARY_SUFFIX}" CACHE FILEPATH "SNOWPACK library" FORCE)
    set(SNOWPACK_FORTRAN_LIBRARY "${SNOWPACK_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}snowpack_fortran${CMAKE_STATIC_LIBRARY_SUFFIX}" CACHE FILEPATH "SNOWPACK Fortran bindings library" FORCE)
    set(SNOWPACK_INCLUDE_DIR "${SNOWPACK_DIR}/include/snowpack" CACHE PATH "SNOWPACK include directory" FORCE)
    
    message(STATUS "SNOWPACK will be built and installed to: ${SNOWPACK_DIR}")
    message(STATUS "SNOWPACK library will be: ${SNOWPACK_LIBRARY}")
    message(STATUS "SNOWPACK include dir will be: ${SNOWPACK_INCLUDE_DIR}")
    
    set(SNOWPACK_FOUND TRUE)
    set(SNOWPACK TRUE CACHE BOOL "Enable SNOWPACK support in HICAR" FORCE)

    # Mark as built from source so the main project can add dependencies
    set(BUILD_SNOWPACK_FROM_SOURCE TRUE CACHE BOOL "SNOWPACK should be built from source" FORCE)
endif()

if( (SNOWPACK OR SNOWPACK_GPU) AND NOT SNOWPACK_FOUND)
	message(WARNING "SNOWPACK not found, disabling SNOWPACK support.\nTo enable, set the SNOWPACK_DIR variable to the root directory of the SNOWPACK install (i.e. the directory which contains the 'lib' and 'build' subdirectories).")
	set(SNOWPACK OFF)
endif()

if (NOT SNOWPACK_FIND_QUIETLY)
    message(STATUS "--------------------------------------------------------------------------------------------------------------------")
    message(STATUS "")
endif()