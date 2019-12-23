# calin/swig/math/CMakeLists.txt -- Stephen Fegan
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

set(CALIN_SWIG_FILES function.i pdf_1d.i histogram.i optimizer.i
	hessian.i special.i hex_array.i rng.i ray.i vector3d_util.i # vs_physics.i
	interpolation_1d.i healpix_array.i ray_generator.i geometry.i moments_calc.i
	data_modeling.i fftw_util.i regular_grid.i lomb_scargle.i
	# simd.i fft_simd.i
	b_spline.i spline_interpolation.i least_squares.i)
set(CALIN_SWIG_TARGET math)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
make_swig_modules(calin_${CALIN_SWIG_TARGET}
	${CALIN_SWIG_INSTALL_DIR}/${CALIN_SWIG_TARGET} ${CALIN_SWIG_FILES})
