# CMakeList.txt : CMake project for LHF, include source and define
# project specific logic here.
#

cmake_minimum_required(VERSION 3.5)

add_library(pipePacket "pipePacket.cpp" "pipePacket.hpp")
add_library(basePipe "basePipe.cpp" "basePipe.hpp")
add_library(distMatrixPipe "distMatrixPipe.cpp" "distMatrixPipe.hpp")
add_library(neighGraphPipe "neighGraphPipe.cpp" "neighGraphPipe.hpp")
add_library(ripsPipe "ripsPipe.cpp" "ripsPipe.hpp")
add_library(upscalePipe "upscalePipe.cpp" "upscalePipe.hpp")
add_library(persistencePairs "persistencePairs.cpp" "persistencePairs.hpp")
add_library(slidingWindow "slidingWindow.cpp" "slidingWindow.hpp")
add_library(fastPersistence "fastPersistence.cpp" "fastPersistence.hpp")
add_library(incrementalPersistence "incrementalPersistence.cpp" "incrementalPersistence.hpp")
add_library(involutedPersistence "involutedPersistence.cpp" "involutedPersistence.hpp")
add_library(naiveWindow "naiveWindow.cpp" "naiveWindow.hpp")

include_directories(../Complex)
link_directories(../Complex)

include_directories(../Pipes)
link_directories(../Pipes)

include_directories(../Utils)
link_directories(../Utils)

# target_link_libraries(LHF PipeFunctions)
target_link_libraries(basePipe utils simplexBase pipePacket distMatrixPipe neighGraphPipe ripsPipe upscalePipe persistencePairs slidingWindow fastPersistence incrementalPersistence naiveWindow involutedPersistence)

INSTALL (TARGETS basePipe
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib)

# TODO: Add tests and install targets if needed.

MESSAGE (STATUS "CMAKE_INCLUDE_PATH: " ${CMAKE_INCLUDE_PATH})
MESSAGE (STATUS "FLAGS: " ${CMAKE_CXX_FLAGS})
