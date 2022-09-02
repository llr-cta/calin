# calin/swig/simulation/CMakeLists.txt -- Stephen Fegan
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

set(CALIN_SWIG_FILES atmosphere.i pmt.i vs_optics.i ray_processor.i tracker.i
	detector_efficiency.i world_magnetic_model.i vcl_iact.i sct_optics.i
	ray_propagator.i)

if(GEANT4_FOUND)
	set(CALIN_SWIG_FILES ${CALIN_SWIG_FILES} geant4_shower_generator.i)
endif()

set(CALIN_SWIG_TARGET simulation)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
make_swig_modules(calin_${CALIN_SWIG_TARGET}
	${CALIN_SWIG_INSTALL_DIR}/${CALIN_SWIG_TARGET} ${CALIN_SWIG_FILES})
