add_executable(test_histogram test_histogram.cpp)  
target_link_libraries(test_histogram ${GTEST_BOTH_LIBRARIES} calin_math calin_proto)
add_test(TestHistogram test_histogram)


