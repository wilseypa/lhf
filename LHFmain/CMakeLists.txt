# CMakeList.txt : CMake project for LHF, include source and define
# project specific logic here.
#

#IF(DEFINED MPI_CXX_COMPILER)
#  SET(CMAKE_CXX_COMPILER "${MPI_CXX_COMPILER}" CACHE FILEPATH "CXX compiler overridden with MPI C++ wrapper")
#ENDIF(DEFINED MPI_CXX_COMPILER)

cmake_minimum_required(VERSION 3.11)



# Use external project for qhull integration
include(FetchContent)

FetchContent_Declare(qhullext
    GIT_REPOSITORY https://github.com/qhull/qhull.git
)

FetchContent_MakeAvailable(qhullext)
include_directories(./_deps/qhullext-src/src/)

find_package(MPI REQUIRED)
add_definitions(-DOMPI_SKIP_MPICXX)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

find_package(OpenMP)

set (CMAKE_CXX_STANDARD 14)
set (CMAKE_CXX_STANDARD_LIBRARIES "-lpthread")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -fopenmp")

include_directories(../Complex)
link_directories(../Complex)

include_directories(../Pipes)
link_directories(../Pipes)

include_directories(../Preprocessing)
link_directories(../Preprocessing)

include_directories(../Utils)
link_directories(../Utils)

#Add Source to this project's executable
add_library(LHFlib SHARED "LHF.cpp" "LHF.hpp")
set_target_properties(LHFlib PROPERTIES VERSION ${PROJECT_VERSION})
set_target_properties(LHFlib PROPERTIES SOVERSION 1)

target_link_libraries(LHFlib simplexBase pipePacket utils readInput writeOutput basePipe preprocessor OpenMP::OpenMP_CXX MPI::MPI_C qhullcpp qhullstatic_r)

INSTALL (TARGETS LHFlib
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib)

add_executable(LHF "main.cpp")
target_link_libraries(LHF PRIVATE LHFlib)

MESSAGE (STATUS "CMAKE_INCLUDE_PATH: " ${CMAKE_INCLUDE_PATH})
MESSAGE (STATUS "FLAGS: " ${CMAKE_CXX_FLAGS})
MESSAGE (STATUS "Ext Dir: " ${EXTERNAL_INSTALL_LOCATION})
