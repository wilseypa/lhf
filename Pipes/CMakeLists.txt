# CMakeList.txt : CMake project for LHF, include source and define
# project specific logic here.
#

cmake_minimum_required(VERSION 3.11)

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

find_package(CGAL)
#Integrating with CGAL library
if(NOT CGAL_FOUND)
    FetchContent_Declare(cgalext
        GIT_REPOSITORY "https://github.com/CGAL/cgal.git"
        GIT_TAG "v5.6.1"
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
endif()

find_package(OpenMP REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(TBB)

add_library(pipePacket STATIC "pipePacket.cpp" "pipePacket.hpp")
target_link_libraries(pipePacket PUBLIC simplexBase)
target_include_directories(pipePacket PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/Complex)

add_library(basePipe STATIC "basePipe.cpp" "basePipe.hpp")
target_link_libraries(basePipe PUBLIC
 distMatrixPipe neighGraphPipe 
 ripsPipe betaSkeletonBasedComplex betaSubSkeletonComplex
 upscalePipe slidingWindow fastPersistence 
 incrementalPersistence naiveWindow qhullPipe delaunayPipe incrementalPipe)
target_include_directories(basePipe PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library(distMatrixPipe STATIC "distMatrixPipe.cpp" "distMatrixPipe.hpp")
target_link_libraries(distMatrixPipe PUBLIC utils basePipe)
target_include_directories(distMatrixPipe PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library(neighGraphPipe STATIC "neighGraphPipe.cpp" "neighGraphPipe.hpp")
target_link_libraries(neighGraphPipe PUBLIC basePipe)
target_include_directories(neighGraphPipe PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library(ripsPipe STATIC "ripsPipe.cpp" "ripsPipe.hpp")
target_link_libraries(ripsPipe PUBLIC utils basePipe)
target_include_directories(ripsPipe PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library(betaSkeletonBasedComplex STATIC "betaSkeletonBasedComplex.cpp" "betaSkeletonBasedComplex.hpp")
target_link_libraries(betaSkeletonBasedComplex PUBLIC basePipe alphaComplex qhullPipe readInput utils kdTree qhullcpp qhull_r)
target_include_directories(betaSkeletonBasedComplex PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${qhullext_SOURCE_DIR}/src)

add_library(betaSubSkeletonComplex STATIC "betaSubSkeletonComplex.cpp" "betaSubSkeletonComplex.hpp")
target_link_libraries(betaSubSkeletonComplex PUBLIC basePipe alphaComplex readInput utils kdTree)
target_include_directories(betaSubSkeletonComplex PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library(upscalePipe STATIC "upscalePipe.cpp" "upscalePipe.hpp")
target_link_libraries(upscalePipe PUBLIC basePipe utils)
target_include_directories(upscalePipe PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library(slidingWindow STATIC "slidingWindow.cpp" "slidingWindow.hpp")
target_link_libraries(slidingWindow PUBLIC basePipe readInput utils)
target_include_directories(slidingWindow PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library(fastPersistence STATIC "fastPersistence.cpp" "fastPersistence.hpp")
target_link_libraries(fastPersistence PUBLIC basePipe utils unionFind)
target_include_directories(fastPersistence PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library(incrementalPersistence STATIC "incrementalPersistence.cpp" "incrementalPersistence.hpp")
target_link_libraries(incrementalPersistence PUBLIC basePipe simplexBase simplexArrayList utils unionFind)
target_include_directories(incrementalPersistence PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library(naiveWindow STATIC "naiveWindow.cpp" "naiveWindow.hpp")
target_link_libraries(naiveWindow PUBLIC basePipe readInput utils)
target_include_directories(naiveWindow PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library(qhullPipe STATIC "qhullPipe.cpp" "qhullPipe.hpp")
target_link_libraries(qhullPipe PUBLIC basePipe utils alphaComplex qhullcpp qhull_r)
target_include_directories(qhullPipe PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${qhullext_SOURCE_DIR}/src)

add_library(delaunayPipe STATIC "delaunayPipe.cpp" "delaunayPipe.hpp")
target_link_libraries(delaunayPipe PUBLIC basePipe alphaComplex utils CGAL::CGAL OpenMP::OpenMP_CXX)
target_include_directories(delaunayPipe PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

if(APPLE)
    add_compile_definitions(NO_PARALLEL_ALGORITHMS=1)
endif()

add_library(incrementalPipe STATIC "incrementalPipe.cpp" "incrementalPipe.hpp")
if(TBB_FOUND)
target_link_libraries(incrementalPipe PUBLIC basePipe utils OpenMP::OpenMP_CXX TBB::tbb)
else()
target_link_libraries(incrementalPipe PUBLIC basePipe utils OpenMP::OpenMP_CXX)
endif()
target_include_directories(incrementalPipe PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${EIGEN3_INCLUDE_DIR})

INSTALL (TARGETS basePipe
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib)

# TODO: Add tests and install targets if needed.

MESSAGE (STATUS "CMAKE_INCLUDE_PATH: " ${CMAKE_INCLUDE_PATH})
MESSAGE (STATUS "FLAGS: " ${CMAKE_CXX_FLAGS})