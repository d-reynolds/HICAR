# Always set the cmake min version.
cmake_minimum_required(VERSION 3.20)

project(HICAR LANGUAGES Fortran)

enable_language (Fortran)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/../cmake")

find_package(PkgConfig REQUIRED)
include(ExternalProject)


# Define test data directory
set(TEST_DATA_DIR "${PROJECT_SOURCE_DIR}/Test_Cases")

# Setup the external project for test data
ExternalProject_Add(hicar_test_data
  GIT_REPOSITORY    https://github.com/HICAR-Model/Test-Data.git
  GIT_TAG           main
  SOURCE_DIR        "${TEST_DATA_DIR}"
  CONFIGURE_COMMAND ""  # No configuration needed for data-only repo
  BUILD_COMMAND     ""  # No build needed
  INSTALL_COMMAND   ""  # No install needed
  TEST_COMMAND      ""  # No tests needed
  UPDATE_DISCONNECTED ${FETCHCONTENT_UPDATES_DISCONNECTED}
)

# 6. Link to the package libraries
#set (NETCDF_F90 "YES")
find_package(MPI REQUIRED)
#find_package(OpenMP)
#find_package(FFTW REQUIRED)
#find_package(NetCDF REQUIRED)

if(NOT TARGET "test-drive::test-drive")
  find_package("test-drive" REQUIRED)
endif()

# Unit testing
set(
  tests
  advect
  halo_exch
)
set(
  SOURCE_FILES
  "test_driver.F90"
)
foreach(t IN LISTS tests)
  string(MAKE_C_IDENTIFIER ${t} t) 
  list(APPEND SOURCE_FILES "test_${t}.F90")
endforeach()

add_executable(
  "${PROJECT_NAME}-tester"
  "${SOURCE_FILES}"
)
target_link_libraries(
  "${PROJECT_NAME}-tester"
  PRIVATE
  "hicar_lib"
  "test-drive::test-drive"
)

target_link_libraries("${PROJECT_NAME}-tester" PUBLIC MPI::MPI_Fortran)

target_include_directories("${PROJECT_NAME}-tester" PUBLIC ${PROJECT_SOURCE_DIR}/../build)

#if(OpenMP_CXX_FOUND)
#    target_link_libraries("${PROJECT_NAME}-tester" PUBLIC OpenMP::OpenMP_CXX)
#endif()

foreach(t IN LISTS tests)
  add_test("test_${t}.F90" "${PROJECT_NAME}-tester" "${t}")
endforeach()

# set(CMAKE_BUILD_TYPE Debug)

if (COMPILER STREQUAL "gnu")
  target_compile_options("${PROJECT_NAME}-tester" PUBLIC -c -O0 -g -fopenmp -fbounds-check -fbacktrace -finit-real=nan -ffree-line-length-none -Wconversion)
elseif (COMPILER STREQUAL "intel")
  target_compile_options("${PROJECT_NAME}-tester" PUBLIC -debug -c -O2 -u -check all -check noarg_temp_created -traceback -fpe0 -fast-transcendentals -xhost)
elseif (COMPILER STREQUAL "cray")
  target_compile_options("${PROJECT_NAME}-tester" PUBLIC -h noomp -c -eI -G0 -m2 -rl -R bcdsp) #-R cpbs --cray-print-opts
endif()

# Create a custom target for downloading test data
add_custom_target(download_test_data
  DEPENDS hicar_test_data
  COMMENT "Downloading HICAR test data from GitHub"
)

# Make the HICAR-tester target depend on the test data
add_dependencies("${PROJECT_NAME}-tester" download_test_data)

# Add test data directory to compile definitions
target_compile_definitions(
  "${PROJECT_NAME}-tester"
  PRIVATE
  TEST_DATA_DIR="${TEST_DATA_DIR}"
)