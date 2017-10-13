#configure for fotran use
cmake_minimum_required (VERSION 2.6)
enable_language (Fortran)

#project name
project(2DpotWell)
#version info
set (cmakeTest_VERSION_MAJOR 1)
set (cmakeTest_VERSION_MINOR 0)

#set some directions
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

#use all files provided in src, i.e. my .f90 files
file(GLOB_RECURSE sources src/*.f90 src/*.F90)

#include_directories("${PROJECT-SOURCE_DIR}/src")
#add_subdirectory(src)

include_directories("${PROJECT_SOURCE_DIR}/thirdparty")
add_subdirectory(thirdparty)

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
	message("only intel compiler is supported")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  # ifort (untested)
  	#set (CMAKE_Fortran_FLAGS "-mkl -qopenmp -g -check all -fpe0 -warn -traceback -debug extended")
  	set (CMAKE_Fortran_FLAGS "-mkl -O3 -qopenmp")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")


#set (CMAKE_Fortran_FLAGS "-mkl -O3 -g -Wall -qopenmp")
#linker flags
set(COVERAGE_LINK_FLAGS "-mkl -g  -qopenmp")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${COVERAGE_LINK_FLAGS}")

# generate executable and link library
add_executable(main.exe ${sources})
target_link_libraries(main.exe config_fortran quadpack)

