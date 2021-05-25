# CMakeList.txt : CMake project for LHF, include source and define
# project specific logic here.
#

cmake_minimum_required(VERSION 3.5)

add_library(simplexBase "simplexBase.cpp" "simplexBase.hpp")
add_library(simplexArrayList "simplexArrayList.cpp" "simplexArrayList.hpp")
#add_library(simplexTree "simplexTree.cpp" "simplexTree.hpp")

include_directories(${PROJECT_SOURCE_DIR}/Utils)
link_directories(${PROJECT_SOURCE_DIR}/Utils)

include_directories(../Preprocessing)
link_directories(../Preprocessing)

# target_link_libraries(LHF PipeFunctions)
target_link_libraries(simplexBase utils pipePacket simplexArrayList) #simplexTree)

INSTALL (TARGETS simplexBase
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib)

# TODO: Add tests and install targets if needed.

MESSAGE (STATUS "CMAKE_INCLUDE_PATH: " ${CMAKE_INCLUDE_PATH})
MESSAGE (STATUS "FLAGS: " ${CMAKE_CXX_FLAGS})
