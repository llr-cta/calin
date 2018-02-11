# calin/unit_tests/io/CMakeLists.txt -- Stephen Fegan
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

add_executable(test_log test_log.cpp)
target_link_libraries(test_log ${GTEST_BOTH_LIBRARIES} calin_io ${EXTRA_LIBS})
add_test(TestLog test_log)

add_executable(test_sql_transceiver test_sql_transceiver.cpp)
target_link_libraries(test_sql_transceiver ${GTEST_BOTH_LIBRARIES} calin_io calin_proto ${EXTRA_LIBS})
add_test(TestSQLTransceiver test_sql_transceiver)

add_executable(test_data_source test_data_source.cpp)
target_link_libraries(test_data_source ${GTEST_BOTH_LIBRARIES} calin_io calin_proto ${EXTRA_LIBS})
add_test(TestDataSource test_data_source)

install(TARGETS test_log test_sql_transceiver test_data_source
  DESTINATION ${CALIN_UNIT_TEST_INSTALL_DIR})
