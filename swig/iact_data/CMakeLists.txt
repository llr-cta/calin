# calin/swig/io/CMakeLists.txt -- Stephen Fegan
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

set(CALIN_SWIG_FILES algorithms.i # telescope_data_source.i
	event_visitor.i event_dispatcher.i
	instrument_layout.i nectarcam_layout.i lstcam_layout.i
	waveform_treatment_event_visitor.i nectarcam_ancillary_data.i
	raw_acada_event_data_source_l0.i raw_acada_event_data_source_r1v0.i
	raw_acada_event_data_source_r1v1.i)
set(CALIN_SWIG_TARGET iact_data)

include_directories(${CTA_ADH_APIS_INCLUDE_DIRS})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
make_swig_modules(calin_${CALIN_SWIG_TARGET}
	${CALIN_SWIG_INSTALL_DIR}/${CALIN_SWIG_TARGET} ${CALIN_SWIG_FILES})
