
include_directories(../Utils)
link_directories(../Utils)

include_directories(../Complex)
link_directories(../Complex)

add_executable(test_functions testExample.cpp)
add_executable(test_complex complexTests.cpp)


target_link_libraries(test_complex utils simplexBase simplexArrayList simplexTree indSimplexTree)


add_test(FunctionTest test_functions)
add_test(ComplexTests test_complex)
