cmake_minimum_required(VERSION 3.5)

set (CMAKE_CXX_STANDARD 14)
set (CMAKE_CXX_STANDARD_LIBRARIES "-lpthread")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")

include_directories(../Complex)
link_directories(../Complex)

include_directories(../Pipes)
link_directories(../Pipes)

include_directories(../Preprocessing)
link_directories(../Preprocessing)

include_directories(../Utils)
link_directories(../Utils)

add_executable(test_functions testExample.cpp)
add_executable(test_complex complexTests.cpp)
add_executable(test_pipes pipeTests.cpp)
add_executable(test_write writeTests.cpp)
add_executable(test_read readTests.cpp)
add_executable(test_preproc preprocTests.cpp)

target_link_libraries(test_complex LHFlib)
target_link_libraries(test_pipes LHFlib)
target_link_libraries(test_write LHFlib)
target_link_libraries(test_read LHFlib)
target_link_libraries(test_preproc LHFlib)

add_test(FunctionTest test_functions)
add_test(ComplexTests test_complex)
add_test(PipeTests test_pipes)
add_test(WriteTests test_write)
add_test(ReadTests test_read)
add_test(PreprocTests test_preproc)
