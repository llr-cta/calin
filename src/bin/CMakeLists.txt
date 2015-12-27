# calin/src/bin/CMakeLists.txt -- Stephen Fegan
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

add_executable(do_spe_fit do_spe_fit)
target_link_libraries(do_spe_fit calin_calib calin_math calin_proto calin_nlopt ${EXTRA_LIBS})

add_executable(do_spe_expgauss_fit do_spe_expgauss_fit)
target_link_libraries(do_spe_expgauss_fit calin_calib calin_math calin_proto calin_nlopt ${EXTRA_LIBS})

install(TARGETS do_spe_fit do_spe_expgauss_fit DESTINATION bin)
