cmake_minimum_required(VERSION 3.10)
project(HICAR LANGUAGES Fortran)

enable_language (Fortran)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/../cmake")

find_package(PkgConfig REQUIRED)

# 6. Link to the package libraries
#set (NETCDF_F90 "YES")
#find_package(MPI REQUIRED)
#find_package(OpenMP)
#find_package(FFTW REQUIRED)
#find_package(NetCDF REQUIRED)

if(NOT TARGET "test-drive::test-drive")
  find_package("test-drive" REQUIRED)
endif()

# Unit testing
set(
  tests
  halo_exch
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
target_include_directories("${PROJECT_NAME}-tester" PUBLIC ${PROJECT_SOURCE_DIR}/../build)

#target_link_libraries("${PROJECT_NAME}-tester" PUBLIC MPI::MPI_Fortran)
#if(OpenMP_CXX_FOUND)
#    target_link_libraries("${PROJECT_NAME}-tester" PUBLIC OpenMP::OpenMP_CXX)
#endif()

foreach(t IN LISTS tests)
  add_test("test_${t}.F90" "${PROJECT_NAME}-tester" "${t}")
endforeach()

set(CMAKE_BUILD_TYPE Debug)
set_target_properties("${PROJECT_NAME}-tester" PROPERTIES OUTPUT_NAME "test")

if (COMPILER STREQUAL "gnu")
  target_compile_options("${PROJECT_NAME}-tester" PUBLIC -c -O0 -g -fopenmp -fbounds-check -fbacktrace -finit-real=nan -ffree-line-length-none -Wconversion)
elseif (COMPILER STREQUAL "intel")
  target_compile_options("${PROJECT_NAME}-tester" PUBLIC -debug -c -O2 -u -check all -check noarg_temp_created -traceback -fpe0 -fast-transcendentals -xhost)
elseif (COMPILER STREQUAL "cray")
  target_compile_options("${PROJECT_NAME}-tester" PUBLIC -h noomp -c -eI -G0 -m2 -rl -R bcdsp) #-R cpbs --cray-print-opts
endif()
