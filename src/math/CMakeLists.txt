# calin/src/math/CMakeLists.txt -- Stephen Fegan
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

set(CALIN_TARGET_LIBRARY calin_math)
set(CALIN_SOURCES accumulator.cpp histogram.cpp function.cpp
	pdf_1d.cpp optimizer.cpp nlopt_optimizer.cpp cminpack_optimizer.cpp
	minuit75_optimizer.cpp hessian.cpp brent.cpp
	log_quadratic_spline_pdf_1d.cpp covariance_calc.cpp
	special.cpp hex_array.cpp hex_array_vvv.cpp vs_vec3d.cpp
	vs_vec4d.cpp vs_particle.cpp rng.cpp fftw_util.cpp)

add_library(${CALIN_TARGET_LIBRARY} SHARED ${CALIN_SOURCES})
target_link_libraries(${CALIN_TARGET_LIBRARY}
	calin_nlopt calin_cminpack calin_f2c calin_io calin_proto_math ${EXTRA_LIBS})

install(TARGETS ${CALIN_TARGET_LIBRARY} DESTINATION lib)
