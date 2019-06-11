# calin/src/calib/CMakeLists.txt -- Stephen Fegan
#
# Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
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

add_library(calin_calib SHARED spe_fit.cpp spe_fit_poisson_gaussian.cpp
  spe_fit_general_poisson.cpp spe_robust.cpp pmt_model_pg.cpp
  spe_fit_fast_single_value_general_poisson.cpp
  pmt_ses_models.cpp)
target_link_libraries(calin_calib calin_math calin_proto_calib)
target_link_libraries(calin_calib ${FFTW_LIBRARIES} ${EXTRA_LIBS})

install(TARGETS calin_calib DESTINATION ${CALIN_LIB_INSTALL_DIR})
