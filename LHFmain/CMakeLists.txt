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

set (CMAKE_CXX_STANDARD 14)
set (CMAKE_CXX_STANDARD_LIBRARIES "-lpthread")

if(UNIX AND NOT APPLE)
        set(LINUX TRUE)
endif()
# if(NOT LINUX) should work, too, if you need that
if(LINUX) 
    message(STATUS ">>> Linux")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -fopenmp")
else()
    message(STATUS ">>> Not Linux")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
endif()

include_directories(${PROJECT_SOURCE_DIR}/Complex)
link_directories(${PROJECT_SOURCE_DIR}/Complex)

include_directories(${PROJECT_SOURCE_DIR}/Pipes)
link_directories(${PROJECT_SOURCE_DIR}/Pipes)

include_directories(${PROJECT_SOURCE_DIR}/Preprocessing)
link_directories(${PROJECT_SOURCE_DIR}/Preprocessing)

include_directories(${PROJECT_SOURCE_DIR}/Utils)
link_directories(${PROJECT_SOURCE_DIR}/Utils)

include_directories(${PROJECT_SOURCE_DIR}/_deps/qhullext-src/src/)
link_directories(${PROJECT_SOURCE_DIR}/_deps/qhullext-build/)

#Add Source to this project's executable
add_library(LHFlib SHARED "LHF.cpp" "LHF.hpp")
set_target_properties(LHFlib PROPERTIES VERSION ${PROJECT_VERSION})
set_target_properties(LHFlib PROPERTIES SOVERSION 1)
set_target_properties(LHFlib PROPERTIES POSITION_INDEPENDENT_CODE ON)

#target_link_libraries(LHFlib pipePacket utils readInput writeOutput simplexBase basePipe preprocessor OpenMP::OpenMP_CXX MPI::MPI_CXX qhullcpp qhull_r delaunayPipe)
target_link_libraries(LHFlib pipePacket utils readInput writeOutput simplexBase basePipe preprocessor OpenMP::OpenMP_CXX MPI::MPI_CXX qhullcpp qhull_r)

INSTALL (TARGETS LHFlib
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib)

add_executable(LHF "main.cpp")
MESSAGE (STATUS "CMAKE_INCLUDE_PATH: " ${PROJECT_SOURCE_DIR})
target_link_libraries(LHF PRIVATE LHFlib)

MESSAGE (STATUS "CMAKE_INCLUDE_PATH: " ${CMAKE_INCLUDE_PATH})
MESSAGE (STATUS "FLAGS: " ${CMAKE_CXX_FLAGS})
MESSAGE (STATUS "Ext Dir: " ${EXTERNAL_INSTALL_LOCATION})
