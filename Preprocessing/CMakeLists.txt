# CMakeList.txt : CMake project for LHF, include source and define
# project specific logic here.
#

cmake_minimum_required(VERSION 3.5)

#Source folder libraries to be compiled into libPreprocessor.so
add_library(preprocessor "preprocessor.cpp" "preprocessor.hpp")
add_library(kMeansPlusPlus "kMeansPlusPlus.cpp" "kMeansPlusPlus.hpp")
add_library(streamingKmeans "streamingKmeans.cpp" "streamingKmeans.hpp")
add_library(streamingUtils "streamingUtils.cpp" "streamingUtils.hpp")
add_library(denStream "denStream.cpp" "denStream.hpp")
add_library(dbscan "dbscan.cpp" "dbscan.hpp")
add_library(kdTree "kdTree.cpp" "kdTree.hpp")

#For pipePacket
include_directories(${PROJECT_SOURCE_DIR}/Pipes)
link_directories(${PROJECT_SOURCE_DIR}/Pipes)

#PipePacket requires simplexBase (complexes)
include_directories(${PROJECT_SOURCE_DIR}/Complex)
link_directories(${PROJECT_SOURCE_DIR}/Complex)

#For utility functions
include_directories(${PROJECT_SOURCE_DIR}/Utils)
link_directories(${PROJECT_SOURCE_DIR}/Utils)

target_link_libraries(preprocessor kMeansPlusPlus streamingKmeans streamingUtils denStream dbscan kdTree utils)

INSTALL (TARGETS preprocessor
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib)

# TODO: Add tests and install targets if needed.
MESSAGE (STATUS "CMAKE_INCLUDE_PATH: " ${CMAKE_INCLUDE_PATH})
MESSAGE (STATUS "FLAGS: " ${CMAKE_CXX_FLAGS})
