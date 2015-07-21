add_executable(test_geant4 test_geant4.cpp)   
target_link_libraries(test_geant4 ${GTEST_BOTH_LIBRARIES}
				     calin_air_shower calin_math calin_proto)
add_test(TestGeant4 test_geant4) 
