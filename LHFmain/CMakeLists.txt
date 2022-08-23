# CMakeList.txt : CMake project for LHF, include source and define
# project specific logic here.
#

#IF(DEFINED MPI_CXX_COMPILER)
#  SET(CMAKE_CXX_COMPILER "${MPI_CXX_COMPILER}" CACHE FILEPATH "CXX compiler overridden with MPI C++ wrapper")
#ENDIF(DEFINED MPI_CXX_COMPILER)

cmake_minimum_required(VERSION 3.5)

find_package(MPI REQUIRED)
add_definitions(-DOMPI_SKIP_MPICXX)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

find_package(OpenMP)

set (CMAKE_CXX_STANDARD 14)
set (CMAKE_CXX_STANDARD_LIBRARIES "-lpthread -lmpi")
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
add_executable(LHF "main.cpp")

target_link_libraries(LHF PRIVATE simplexBase pipePacket utils basePipe preprocessor OpenMP::OpenMP_CXX)



# TODO: Add tests and install targets if needed.

MESSAGE (STATUS "CMAKE_INCLUDE_PATH: " ${CMAKE_INCLUDE_PATH})
MESSAGE (STATUS "FLAGS: " ${CMAKE_CXX_FLAGS})
