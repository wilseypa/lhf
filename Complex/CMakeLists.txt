# CMakeList.txt : CMake project for LHF, include source and define
# project specific logic here.
#

cmake_minimum_required(VERSION 3.10)

find_package(OpenMP REQUIRED)

add_library(simplexBase STATIC simplexBase.cpp)
target_link_libraries(simplexBase PUBLIC utils pipePacket simplexTree simplexArrayList alphaComplex witnessComplex betaComplex)
target_include_directories(simplexBase PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/Utils)

add_library(simplexArrayList STATIC simplexArrayList.cpp)
target_link_libraries(simplexArrayList PUBLIC simplexBase)
target_include_directories(simplexArrayList PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library(alphaComplex STATIC alphaComplex.cpp)
target_link_libraries(alphaComplex PUBLIC utils kdTree simplexArrayList OpenMP::OpenMP_CXX)
target_include_directories(alphaComplex PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library(witnessComplex STATIC witnessComplex.cpp)
target_link_libraries(witnessComplex PUBLIC simplexArrayList)
target_include_directories(witnessComplex PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library(simplexTree STATIC simplexTree.cpp)
target_link_libraries(simplexTree PUBLIC utils kdTree simplexBase)
target_include_directories(simplexTree PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library(betaComplex STATIC betaComplex.cpp)
target_link_libraries(betaComplex PUBLIC alphaComplex)
target_include_directories(betaComplex PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

INSTALL (TARGETS simplexBase
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib)

# TODO: Add tests and install targets if needed.

MESSAGE (STATUS "CMAKE_INCLUDE_PATH: " ${CMAKE_INCLUDE_PATH})
MESSAGE (STATUS "FLAGS: " ${CMAKE_CXX_FLAGS})