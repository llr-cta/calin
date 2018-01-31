# calin/unit_tests/simulation/CMakeLists.txt -- Stephen Fegan
#
# Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
# LLR, Ecole Polytechnique, CNRS/IN2P3
#
# This file is part of "calin"
#
# "calin" is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License version 2 or later, as published by
# the Free Software Foundation.
#
# "calin" is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

if(GEANT4_FOUND)
	add_executable(test_geant4 test_geant4.cpp)
	target_link_libraries(test_geant4 ${GTEST_BOTH_LIBRARIES}
				     calin_simulation calin_math calin_proto)
	add_test(TestGeant4 test_geant4)
	install(TARGETS test_geant4
	  DESTINATION ${CALIN_UNIT_TEST_INSTALL_DIR})
endif()

#add_executable(test_vs_optics test_vs_optics.cpp)
#target_link_libraries(test_vs_optics ${GTEST_BOTH_LIBRARIES}
#				     calin_simulation calin_io calin_proto)
#add_test(TestVSOptics test_vs_optics)
