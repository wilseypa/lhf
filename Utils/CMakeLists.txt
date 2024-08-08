# CMakeList.txt : CMake project for LHF, include source and define
# project specific logic here.
#

cmake_minimum_required(VERSION 3.5)

# Find required packages
find_package(Eigen3 REQUIRED)
find_package(OpenMP REQUIRED)

add_library(utils STATIC utils.cpp utils.hpp)
target_link_libraries(utils PUBLIC kdTree Eigen3::Eigen OpenMP::OpenMP_CXX)
target_include_directories(utils PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${EIGEN3_INCLUDE_DIRS} ${OpenMP_CXX_INCLUDE_DIRS})

add_library(unionFind STATIC unionFind.cpp unionFind.hpp)
target_include_directories(unionFind PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library(readInput STATIC readInput.cpp readInput.hpp)
target_include_directories(readInput PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library(argParser STATIC argParser.cpp argParser.hpp)
target_include_directories(argParser PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library(writeOutput STATIC writeOutput.cpp writeOutput.hpp)
target_link_libraries(writeOutput PUBLIC utils)
target_include_directories(writeOutput PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

INSTALL (TARGETS utils
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib)

# TODO: Add tests and install targets if needed.

MESSAGE (STATUS "CMAKE_INCLUDE_PATH: " ${CMAKE_INCLUDE_PATH})
MESSAGE (STATUS "FLAGS: " ${CMAKE_CXX_FLAGS})