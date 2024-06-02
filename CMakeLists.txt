# Always set the cmake min version.
cmake_minimum_required(VERSION 3.0)

project(HICAR)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake")

enable_language (Fortran)


#FILE(GLOB_RECURSE modeFiles build/*.mod) #not working
FILE(GLOB_RECURSE SOURCE_FILES  src/*.F90)

#This sets the location for library installation -- i.e., ${CMAKE_INSTALL_LIBDIR}=/usr/local/lib 
# and also the header for library if it has, in this case we do not have. ${CMAKE_INSTALL_INCLUDEDIR}=/usr/local/include
include(GNUInstallDirs)

    
add_library(HICAR STATIC ${SOURCE_FILES})



# Set the location for library installation -- i.e., /usr/lib in this case
# not really necessary in this example. Use "sudo make install" to apply 
set(CMAKE_INSTALL_LIBDIR /users/dreynold/lib/FSM_interface)
set(CMAKE_Fortran_MODULE_DIRECTORY  /users/dreynold/include/FSM_interface)

# this is to copy the library to the diseired directory
install(TARGETS HICAR  DESTINATION  ${CMAKE_INSTALL_LIBDIR})


# this is to copy build folder including .mod files in $HOME/include directory
# note that CMAKE_Fortran_MODULE_DIRECTORY is ouput path for modules or anything
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/build DESTINATION ${CMAKE_Fortran_MODULE_DIRECTORY}) 
#install(PROGRAMS ${modeFiles} DESTINATION ${CMAKE_Fortran_MODULE_DIRECTORY}/sec) # not working
    
#add_executable(test_FSM_interface test/test_FSM_interface.f90)
#target_link_libraries(test_FSM_interface FSM_interface ${PROJECT_LINK_LIBS})

# this is for gnu
SET(CMAKE_Fortran_FLAGS "-fcoarray=lib -mfma -ffree-line-length-none -ftree-vectorize -fimplicit-none -funroll-loops -march=native  -fno-protect-parens ")

# this is for cray
#SET(CMAKE_Fortran_FLAGS "-h omp vector2 -O3 -c -hfp2 -eI")

ADD_CUSTOM_TARGET( distclean make clean
				   COMMAND ${CMAKE_COMMAND} -E remove_directory CMakeFiles
				   COMMAND ${CMAKE_COMMAND} -E remove *~ cmake_install.cmake install_manifest.txt Makefile CMakeCache.txt
				   COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_Fortran_MODULE_DIRECTORY}
				   COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_INSTALL_LIBDIR}
					)
