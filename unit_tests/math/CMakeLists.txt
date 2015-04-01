add_executable(test_histogram test_histogram.cpp)   
target_link_libraries(test_histogram ${GTEST_BOTH_LIBRARIES} calin_math calin_proto)  
add_test(TestHistogram test_histogram) 

add_executable(test_hessian test_hessian.cpp ../calib/karkar_data.cpp)   
target_link_libraries(test_hessian ${GTEST_BOTH_LIBRARIES} calin_calib
calin_math
calin_proto calin_nlopt)  
add_test(TestHessian test_hessian) 

add_executable(test_function test_function.cpp)   
target_link_libraries(test_function ${GTEST_BOTH_LIBRARIES} calin_math)
add_test(TestFunction test_function) 
