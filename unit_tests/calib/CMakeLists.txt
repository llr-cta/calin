# calin/unit_tests/calib/CMakeLists.txt -- Stephen Fegan
#
# Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
# LLR, Ecole polytechnique, CNRS/IN2P3, Universite Paris-Saclay
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

add_executable(test_spe_fit test_spe_fit.cpp karkar_data.cpp)
target_link_libraries(test_spe_fit ${GTEST_BOTH_LIBRARIES} calin_calib
calin_math calin_proto calin_nlopt calin_f2c ${EXTRA_LIBS})
add_test(TestSpeFit test_spe_fit)

add_executable(test_pmt_ses_models test_pmt_ses_models.cpp)
target_link_libraries(test_pmt_ses_models ${GTEST_BOTH_LIBRARIES} calin_calib
  calin_math ${EXTRA_LIBS})
add_test(TestPMTSESModels test_pmt_ses_models)
