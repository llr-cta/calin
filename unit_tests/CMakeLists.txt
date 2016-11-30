# calin/unit_tests/CMakeLists.txt -- Stephen Fegan
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

add_subdirectory(gtest-1.8.0)
include_directories(gtest-1.8.0/include)
set(GTEST_BOTH_LIBRARIES gtest gtest_main)

add_subdirectory(simulation)
add_subdirectory(math)
add_subdirectory(calib)
add_subdirectory(io)
