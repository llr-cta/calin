# calin/swig/diagnostics/CMakeLists.txt -- Stephen Fegan -- 2016-03-04
#
# Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
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

set(CALIN_SWIG_FILES counter.i module.i event_number.i waveform.i functional.i
	value_capture.i delta_t.i)
set(CALIN_SWIG_TARGET diagnostics)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
make_swig_modules(calin_${CALIN_SWIG_TARGET}
	${CALIN_SWIG_INSTALL_DIR}/${CALIN_SWIG_TARGET} ${CALIN_SWIG_FILES})
