# CMakeList.txt : CMake project for LHF, include source and define
# project specific logic here.
#

cmake_minimum_required(VERSION 3.5)

add_library(preprocessor STATIC preprocessor.cpp preprocessor.hpp)
target_link_libraries(preprocessor PUBLIC kMeansPlusPlus pipePacket)
target_include_directories(preprocessor PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library(kMeansPlusPlus STATIC kMeansPlusPlus.cpp kMeansPlusPlus.hpp)
target_link_libraries(kMeansPlusPlus PUBLIC utils preprocessor cluster)
target_include_directories(kMeansPlusPlus PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library(streamingKmeans STATIC streamingKmeans.cpp streamingKmeans.hpp)
target_link_libraries(streamingKmeans PUBLIC utils streamingUtils preprocessor)
target_include_directories(streamingKmeans PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library(streamingUtils STATIC streamingUtils.cpp streamingUtils.hpp)
target_link_libraries(streamingUtils PUBLIC utils preprocessor)
target_include_directories(streamingUtils PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library(denStream STATIC denStream.cpp denStream.hpp)
target_link_libraries(denStream PUBLIC utils dbscan preprocessor)
target_include_directories(denStream PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library(dbscan STATIC dbscan.cpp dbscan.hpp)
target_link_libraries(dbscan PUBLIC utils kdTree kDTree)
target_include_directories(dbscan PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library(kdTree STATIC kdTree.cpp kdTree.hpp)
target_link_libraries(kdTree PUBLIC utils)
target_include_directories(kdTree PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

INSTALL (TARGETS preprocessor
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib)

# TODO: Add tests and install targets if needed.
MESSAGE (STATUS "CMAKE_INCLUDE_PATH: " ${CMAKE_INCLUDE_PATH})
MESSAGE (STATUS "FLAGS: " ${CMAKE_CXX_FLAGS})