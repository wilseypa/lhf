# CMakeList.txt : Top-level CMake project file, do global configuration 
# and include sub-projects here.
#
enable_testing()

cmake_minimum_required(VERSION 3.5)
set (CMAKE_CXX_STANDARD_LIBRARIES "-lpthread")
set (CMAKE_CXX_STANDARD 14)

project("LHF")

# Include sub-projects
add_subdirectory("LHFmain")
add_subdirectory("Complex")
add_subdirectory("Pipes")
add_subdirectory("Utils")
add_subdirectory("Preprocessing")

add_subdirectory(test)
