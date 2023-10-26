# calin/src/math/CMakeLists.txt -- Stephen Fegan
#
# Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
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

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --save-temps")

set(CALIN_TARGET_LIBRARY calin_math)
set(CALIN_SOURCES accumulator.cpp histogram.cpp function.cpp
	pdf_1d.cpp pdf_1d_gaussian.cpp pdf_1d_limited_exponential.cpp
	pdf_1d_two_component.cpp pdf_1d_binned_gaussian.cpp
  pdf_1d_two_component_constrained.cpp
	hessian.cpp brent.cpp
	log_quadratic_spline_pdf_1d.cpp covariance_calc.cpp
	special.cpp hex_array.cpp hex_array_vvv.cpp square_array.cpp
	geometry.cpp vector3d_util.cpp ray.cpp
	rng.cpp fftw_util.cpp healpix_array.cpp position_generator.cpp
	direction_generator.cpp ray_generator.cpp data_modeling.cpp m_estimate.cpp
	regular_grid.cpp lomb_scargle.cpp	b_spline.cpp
	spline_interpolation.cpp spline_interpolation_2d.cpp least_squares.cpp nspace.cpp
	rng_gaussian_ziggurat.cpp rng_exponential_ziggurat.cpp
	rng_x_exp_minus_x_squared_ziggurat.cpp)

add_library(${CALIN_TARGET_LIBRARY} SHARED ${CALIN_SOURCES})
target_link_libraries(${CALIN_TARGET_LIBRARY} calin_util
	calin_chealpix calin_io calin_provenance
	calin_proto_math ${EXTRA_LIBS}) # -lblosc)

install(TARGETS ${CALIN_TARGET_LIBRARY} DESTINATION ${CALIN_LIB_INSTALL_DIR})
