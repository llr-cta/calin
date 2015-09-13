add_executable(test_log test_log.cpp)  
target_link_libraries(test_log ${GTEST_BOTH_LIBRARIES} calin_io ${EXTRA_LIBS})
add_test(TestLog test_log)

add_executable(test_sql_transceiver test_sql_transceiver.cpp)  
target_link_libraries(test_sql_transceiver ${GTEST_BOTH_LIBRARIES} calin_io calin_proto ${EXTRA_LIBS})
add_test(TestSQLTransceiver test_sql_transceiver)
