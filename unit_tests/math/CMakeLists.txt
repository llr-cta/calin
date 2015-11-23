add_executable(test_histogram test_histogram.cpp)   
target_link_libraries(test_histogram ${GTEST_BOTH_LIBRARIES} calin_math calin_proto)  
add_test(TestHistogram test_histogram) 

add_executable(test_hessian test_hessian.cpp ../calib/karkar_data.cpp)   
target_link_libraries(test_hessian ${GTEST_BOTH_LIBRARIES} calin_calib
				   calin_math calin_proto calin_nlopt)  
add_test(TestHessian test_hessian) 

add_executable(test_pdf_1d test_pdf_1d.cpp)   
target_link_libraries(test_pdf_1d ${GTEST_BOTH_LIBRARIES} calin_math ${EXTRA_LIBS})
add_test(TestFunction test_pdf_1d) 

add_executable(test_hex_array test_hex_array.cpp)   
target_link_libraries(test_hex_array ${GTEST_BOTH_LIBRARIES} calin_math ${EXTRA_LIBS})
add_test(TestHexArray test_hex_array) 

add_executable(test_vs_vec3d test_vs_vec3d.cpp)   
target_link_libraries(test_vs_vec3d ${GTEST_BOTH_LIBRARIES} calin_math ${EXTRA_LIBS})
add_test(TestVSVec3D test_vs_vec3d) 

add_executable(test_rng test_rng.cpp)   
target_link_libraries(test_rng ${GTEST_BOTH_LIBRARIES} calin_math ${EXTRA_LIBS})
add_test(TestRNG test_rng) 
