# CMakeList.txt : CMake project for LHF, include source and define
# project specific logic here.
#

cmake_minimum_required(VERSION 3.11)

add_library(pipePacket "pipePacket.cpp" "pipePacket.hpp")
add_library(basePipe "basePipe.cpp" "basePipe.hpp")
add_library(distMatrixPipe "distMatrixPipe.cpp" "distMatrixPipe.hpp")
add_library(neighGraphPipe "neighGraphPipe.cpp" "neighGraphPipe.hpp")
add_library(ripsPipe "ripsPipe.cpp" "ripsPipe.hpp")
add_library(betaSkeletonBasedComplex "betaSkeletonBasedComplex.cpp" "betaSkeletonBasedComplex.hpp")
add_library(betaSubSkeletonComplex "betaSubSkeletonComplex.cpp" "betaSubSkeletonComplex.hpp")
add_library(upscalePipe "upscalePipe.cpp" "upscalePipe.hpp")
add_library(slidingWindow "slidingWindow.cpp" "slidingWindow.hpp")
add_library(fastPersistence "fastPersistence.cpp" "fastPersistence.hpp")
add_library(incrementalPersistence "incrementalPersistence.cpp" "incrementalPersistence.hpp")
add_library(naiveWindow "naiveWindow.cpp" "naiveWindow.hpp")
add_library(qhullPipe "qhullPipe.cpp" "qhullPipe.hpp")
add_library(delaunayPipe "delaunayPipe.cpp" "delaunayPipe.hpp")


include_directories(${PROJECT_SOURCE_DIR}/Complex)
link_directories(${PROJECT_SOURCE_DIR}/Complex)

include_directories(${PROJECT_SOURCE_DIR}/Preprocessing)
link_directories(${PROJECT_SOURCE_DIR}/Preprocessing)

include_directories(${PROJECT_SOURCE_DIR}/Utils)
link_directories(${PROJECT_SOURCE_DIR}/Utils)

# Use external project for qhull integration
include(FetchContent)

FetchContent_Declare(qhullext
    GIT_REPOSITORY "https://github.com/qhull/qhull.git"
    #BUILD_STATIC_LIBS OFF
    #BUILD_SHARED_LIBS ON
)

if(${CMAKE_VERSION} VERSION_LESS 3.14)
    FetchContent_GetProperties(qhullext)
    if(NOT qhullext_POPULATED)
        FetchContent_Populate(qhullext)
        add_subdirectory(${qhullext_SOURCE_DIR} ${qhullext_BINARY_DIR})
    endif()
else()
    FetchContent_MakeAvailable(qhullext)
endif()
    
include_directories(${PROJECT_SOURCE_DIR}/_deps/qhullext-src/src/)
link_directories(${PROJECT_SOURCE_DIR}/_deps/qhullext-build/)


#Integrating with CGAL library
FetchContent_Declare(cgalext
    GIT_REPOSITORY "https://github.com/CGAL/cgal.git"
    #BUILD_STATIC_LIBS OFF
    #BUILD_SHARED_LIBS ON
)

if(${CMAKE_VERSION} VERSION_LESS 3.14)
    FetchContent_GetProperties(cgalext)
    if(NOT cgalext_POPULATED)
        FetchContent_Populate(cgalext)
        add_subdirectory(${cgalext_SOURCE_DIR} ${cgalext_BINARY_DIR})
    endif()
else()
    FetchContent_MakeAvailable(cgalext)
endif()

set(CGAL_DIR ${cgalext_BINARY_DIR})

find_package(CGAL REQUIRED)
find_package(OpenMP)
find_package(Eigen3)
include_directories(${EIGEN3_INCLUDE_DIR})
target_link_libraries(qhullPipe qhullcpp qhull_r)

#target_link_libraries(basePipe utils simplexBase pipePacket distMatrixPipe neighGraphPipe incrementalPersistence fastPersistence ripsPipe naiveWindow betaSkeletonBasedComplex betaSubSkeletonComplex upscalePipe qhullcpp qhull_r qhullPipe slidingWindow)
target_link_libraries(basePipe utils simplexBase pipePacket distMatrixPipe neighGraphPipe incrementalPersistence fastPersistence ripsPipe naiveWindow upscalePipe slidingWindow betaSkeletonBasedComplex betaSubSkeletonComplex delaunayPipe CGAL::CGAL OpenMP::OpenMP_CXX qhullPipe)

INSTALL (TARGETS basePipe
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib)

# TODO: Add tests and install targets if needed.

MESSAGE (STATUS "CMAKE_INCLUDE_PATH: " ${CMAKE_INCLUDE_PATH})
MESSAGE (STATUS "FLAGS: " ${CMAKE_CXX_FLAGS})
