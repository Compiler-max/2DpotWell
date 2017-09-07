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



include_directories("${PROJECT_SOURCE_DIR}/thirdparty")
add_subdirectory(thirdparty)



#complie flags
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qopenmp -O3")
#linker flags
set(COVERAGE_LINK_FLAGS "-mkl ${OPENMP_LINK_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${COVERAGE_LINK_FLAGS}")

# generate executable and link library
set(EXEC main.exe)
add_executable(${EXEC} ${sources})
target_link_libraries(${EXEC} config_fortran quadpack)


