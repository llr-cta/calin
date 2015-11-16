if(GEANT4_FOUND)
	add_executable(test_geant4 test_geant4.cpp)   
	target_link_libraries(test_geant4 ${GTEST_BOTH_LIBRARIES}
				     calin_simulation calin_math calin_proto)
	add_test(TestGeant4 test_geant4)
endif()

add_executable(test_vs_optics test_vs_optics.cpp)   
target_link_libraries(test_vs_optics ${GTEST_BOTH_LIBRARIES}
				     calin_simulation calin_io calin_proto)
add_test(TestVSOptics test_vs_optics)
