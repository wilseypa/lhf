# CMakeList.txt : CMake project for LHF, include source and define
# project specific logic here.
#

cmake_minimum_required(VERSION 3.5)

add_library(argParser argParser.cpp argParser.hpp)
add_library(readInput readInput.cpp readInput.hpp)
add_library(writeOutput writeOutput.cpp writeOutput.hpp)
add_library(cluster cluster.cpp cluster.hpp)
add_library(utils utils.cpp utils.hpp)

include_directories(../Pipes)
link_directories(../Pipes)

include_directories(../Complex)
link_directories(../Complex)

# target_link_libraries(LHF PipeFunctions)
target_link_libraries(utils argParser readInput writeOutput cluster)

INSTALL (TARGETS utils
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib)

# TODO: Add tests and install targets if needed.

MESSAGE (STATUS "CMAKE_INCLUDE_PATH: " ${CMAKE_INCLUDE_PATH})
MESSAGE (STATUS "FLAGS: " ${CMAKE_CXX_FLAGS})
