# CMakeList.txt : CMake project for LHF, include source and define
# project specific logic here.
#

cmake_minimum_required(VERSION 3.10)

add_library(simplexBase STATIC simplexBase.cpp simplexBase.hpp)
target_link_libraries(simplexBase PUBLIC utils pipePacket simplexArrayList alphaComplex witnessComplex betaComplex)
target_include_directories(simplexBase PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/Utils)

add_library(simplexArrayList STATIC simplexArrayList.cpp simplexArrayList.hpp)
target_link_libraries(simplexArrayList PUBLIC simplexBase)
target_include_directories(simplexArrayList PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library(alphaComplex STATIC alphaComplex.cpp alphaComplex.hpp)
target_link_libraries(alphaComplex PUBLIC utils kdTree simplexArrayList simplexBase)
target_include_directories(alphaComplex PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library(witnessComplex STATIC witnessComplex.cpp witnessComplex.hpp)
target_link_libraries(witnessComplex PUBLIC simplexArrayList simplexBase)
target_include_directories(witnessComplex PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
#add_library(simplexTree "simplexTree.cpp" "simplexTree.hpp")

install(TARGETS simplexBase simplexArrayList alphaComplex witnessComplex DESTINATION lib)
install(FILES simplexBase.hpp simplexArrayList.hpp alphaComplex.hpp witnessComplex.hpp DESTINATION include)

# TODO: Add tests and install targets if needed.

MESSAGE (STATUS "CMAKE_INCLUDE_PATH: " ${CMAKE_INCLUDE_PATH})
MESSAGE (STATUS "FLAGS: " ${CMAKE_CXX_FLAGS})