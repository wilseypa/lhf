# CMakeList.txt : CMake project for LHF, include source and define
# project specific logic here.
#

#IF(DEFINED MPI_CXX_COMPILER)
#  SET(CMAKE_CXX_COMPILER "${MPI_CXX_COMPILER}" CACHE FILEPATH "CXX compiler overridden with MPI C++ wrapper")
#ENDIF(DEFINED MPI_CXX_COMPILER)

cmake_minimum_required(VERSION 3.11)

find_package(MPI REQUIRED)
add_definitions(-DOMPI_SKIP_MPICXX)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

find_package(OpenMP REQUIRED)

#Add Source to this project's executable
add_library(LHFlib SHARED "LHF.cpp" "LHF.hpp")
set_target_properties(LHFlib PROPERTIES VERSION ${PROJECT_VERSION})
set_target_properties(LHFlib PROPERTIES SOVERSION 1)
set_target_properties(LHFlib PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_link_libraries(LHFlib PUBLIC readInput argParser basePipe writeOutput pipePacket preprocessor utils OpenMP::OpenMP_CXX MPI::MPI_CXX)

INSTALL (TARGETS LHFlib
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib)

add_executable(LHF "main.cpp")
MESSAGE (STATUS "CMAKE_INCLUDE_PATH: " ${PROJECT_SOURCE_DIR})
target_link_libraries(LHF PRIVATE LHFlib OpenMP::OpenMP_CXX MPI::MPI_CXX)

MESSAGE (STATUS "CMAKE_INCLUDE_PATH: " ${CMAKE_INCLUDE_PATH})
MESSAGE (STATUS "FLAGS: " ${CMAKE_CXX_FLAGS})
MESSAGE (STATUS "Ext Dir: " ${EXTERNAL_INSTALL_LOCATION})