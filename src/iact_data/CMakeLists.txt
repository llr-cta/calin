# calin/src/io/CMakeLists.txt -- Stephen Fegan
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

set(CALIN_TARGET_LIBRARY calin_iact_data)
set(CALIN_SOURCES algorithms.cpp
	telescope_data_source.cpp event_visitor.cpp
	parallel_event_dispatcher.cpp
	instrument_layout.cpp
	nectarcam_layout.cpp nectarcam_configuration.cpp
	lstcam_layout.cpp
	waveform_treatment_event_visitor.cpp waveform_treatment_event_visitor_impl.cpp
	nectarcam_ancillary_data.cpp
	acada_data_source.cpp zfits_acada_data_source.cpp 
	acada_event_decoder.cpp
	nectarcam_acada_event_decoder_l0.cpp 
	nectarcam_acada_event_decoder_r1v0.cpp
	lstcam_acada_event_decoder_r1v0.cpp
	unified_acada_event_decoder_r1v1.cpp
	cta_acada_event_decoder_r1v0.cpp
	zfits_data_source.cpp
	nectarcam_data_source.cpp
	lstcam_data_source.cpp
	cta_data_source.cpp)


if(CTA_CAMERASTOACTL_FOUND)
	set(CALIN_SOURCES ${CALIN_SOURCES}
		
		)
endif()

add_library(${CALIN_TARGET_LIBRARY} SHARED ${CALIN_SOURCES})
target_link_libraries(${CALIN_TARGET_LIBRARY}
	calin_proto_io calin_proto_iact_data calin_io calin_math calin_util)

include_directories(${CTA_ADH_APIS_INCLUDE_DIRS})
target_link_libraries(${CALIN_TARGET_LIBRARY} ${CTA_ADH_APIS_LIBRARIES})

install(TARGETS ${CALIN_TARGET_LIBRARY} DESTINATION ${CALIN_LIB_INSTALL_DIR})
