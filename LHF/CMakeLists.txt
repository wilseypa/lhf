# CMakeList.txt : CMake project for LHF, include source and define
# project specific logic here.
#

cmake_minimum_required(VERSION 3.5)

set (CMAKE_CXX_STANDARD 14)
set (CMAKE_CXX_STANDARD_LIBRARIES "-lpthread")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")

include_directories(../Pipes)
link_directories(../Pipes)

include_directories(../Complex)
link_directories(../Complex)

include_directories(../Utils)
link_directories(../Utils)

#Add Source to this project's executable
add_executable(LHF "main.cpp")

add_library(combinedLib INTERFACE)

#target_link_libraries(combinedLib INTERFACE basePipe simplexBase utils readInput writeOutput argParser pipePacket bettiPipe ripsPipe neighGraphPipe simplexTree simplexArrayList distMatrixPipe)
target_link_libraries(combinedLib INTERFACE basePipe simplexBase utils readInput writeOutput argParser pipePacket simplexBase simplexTree simplexArrayList distMatrixPipe neighGraphPipe bettiPipe ripsPipe)

target_link_libraries(LHF combinedLib)

# TODO: Add tests and install targets if needed.

MESSAGE (STATUS "CMAKE_INCLUDE_PATH: " ${CMAKE_INCLUDE_PATH})
MESSAGE (STATUS "FLAGS: " ${CMAKE_CXX_FLAGS})
