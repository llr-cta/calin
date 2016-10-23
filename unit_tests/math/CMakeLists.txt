# calin/unit_tests/math/CMakeLists.txt -- Stephen Fegan
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

#add_executable(test_vs_vec3d_obsolete test_vs_vec3d_obsolete.cpp)
#target_link_libraries(test_vs_vec3d_obsolete ${GTEST_BOTH_LIBRARIES} calin_math ${EXTRA_LIBS})
#add_test(TestVSVec3DObsolete test_vs_vec3d_obsolete)

#add_executable(test_vs_vec3d test_vs_vec3d.cpp)
#target_link_libraries(test_vs_vec3d ${GTEST_BOTH_LIBRARIES} calin_math ${EXTRA_LIBS})
#add_test(TestVSVec3D test_vs_vec3d)

add_executable(test_rng test_rng.cpp)
target_link_libraries(test_rng ${GTEST_BOTH_LIBRARIES} calin_math ${EXTRA_LIBS})
add_test(TestRNG test_rng)
