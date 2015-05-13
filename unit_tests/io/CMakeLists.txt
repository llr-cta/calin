add_executable(test_log test_log.cpp)  
target_link_libraries(test_log ${GTEST_BOTH_LIBRARIES} calin_io ${EXTRA_LIBS})
add_test(TestLog test_log)
