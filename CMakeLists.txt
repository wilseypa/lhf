# CMakeList.txt : Top-level CMake project file, do global configuration 
# and include sub-projects here.
#
enable_testing()

cmake_minimum_required(VERSION 3.11)
set (CMAKE_CXX_STANDARD_LIBRARIES "-lpthread")
set (CMAKE_CXX_STANDARD 14)

if(UNIX AND NOT APPLE)
        set(LINUX TRUE)
endif()
# if(NOT LINUX) should work, too, if you need that

SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS} -O1")

if(LINUX) 
    message(STATUS ">>> Linux")
else()
    message(STATUS ">>> Not Linux")
    link_libraries("-undefined dynamic_lookup")
endif()

SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}")
SET(CMAKE_POSITION_INDEPENDENT_CODE ON)

project("LHF" VERSION 1.0.0 DESCRIPTION "Lightweight Homology Framework")

# Include sub-projects
add_subdirectory("LHFmain")
add_subdirectory("Complex")
add_subdirectory("Pipes")
add_subdirectory("Utils")
add_subdirectory("Preprocessing")

add_subdirectory(cmake_tests)
