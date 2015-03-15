add_executable(test_spe_fit test_spe_fit.cpp karkar_data.cpp)  
target_link_libraries(test_spe_fit ${GTEST_BOTH_LIBRARIES} calin_calib
calin_math calin_proto calin_nlopt calin_f2c ${EXTRA_LIBS})
add_test(TestSpeFit test_spe_fit)
