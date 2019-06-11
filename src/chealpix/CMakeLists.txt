# calin/src/cminpack/CMakeLists.txt -- Stephen Fegan
#
# Copyright 2017, Stephen Fegan <sfegan@llr.in2p3.fr>
# LLR, CNRS, Ecole Polytechnique, Institut Polytechnique de Paris
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

set (chealpix_srcs chealpix.c)

include_directories(${CMAKE_SOURCE_DIR}/include/chealpix)
add_library (calin_chealpix SHARED ${chealpix_srcs})

install(TARGETS calin_chealpix DESTINATION ${CALIN_LIB_INSTALL_DIR})
