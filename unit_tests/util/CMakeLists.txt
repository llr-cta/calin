# calin/unit_tests/util/CMakeLists.txt -- Stephen Fegan
#
# Copyright 2018, Stephen Fegan <sfegan@llr.in2p3.fr>
# Laboratoire Leprince-Ringuet, CNRS/IN2P3, Ecole Polytechnique, Institut Polytechnique de Paris
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

add_executable(test_vcl test_vcl.cpp)
target_link_libraries(test_vcl ${GTEST_BOTH_LIBRARIES}
	calin_math calin_util)
#add_test(TestVCL test_vcl)
install(TARGETS test_vcl
  DESTINATION ${CALIN_UNIT_TEST_INSTALL_DIR})
