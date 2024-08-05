# CMakeList.txt : Top-level CMake project file, do global configuration 
# and include sub-projects here.

cmake_minimum_required(VERSION 3.11)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_LIBRARIES "-lpthread")

enable_testing()

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")

if(UNIX AND NOT APPLE) 
    set(LINUX TRUE)
    message(STATUS ">>> Linux")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
elseif(APPLE)
    set(MACOS TRUE)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xpreprocessor -fopenmp")
    set (OpenMP_C_FLAGS "-Xpreprocessor -fopenmp")
    set (OpenMP_C_LIB_NAMES "omp")
    set (OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp")
    set (OpenMP_CXX_LIB_NAMES "omp")
    set (OpenMP_omp_LIBRARY "/opt/homebrew/opt/libomp/lib/libomp.dylib")
    include_directories(/opt/homebrew/opt/libomp/include)
    include_directories(/opt/homebrew/include)
    link_libraries("-undefined dynamic_lookup")
    message(STATUS ">>> MacOS")
else()
    set(WINDOWS TRUE)
    message(STATUS ">>> Windows")
endif()

SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}")
SET(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Project information
project("LHF" VERSION 1.0.0 DESCRIPTION "Lightweight Homology Framework")

# Include sub-projects
add_subdirectory("Complex")
add_subdirectory("Pipes")
add_subdirectory("Utils")
add_subdirectory("Preprocessing")
add_subdirectory("LHFmain")

#add_subdirectory(cmake_tests)