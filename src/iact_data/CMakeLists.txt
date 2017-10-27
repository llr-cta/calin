# calin/src/io/CMakeLists.txt -- Stephen Fegan
#
# Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
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

set(CALIN_TARGET_LIBRARY calin_iact_data)
set(CALIN_SOURCES telescope_data_source.cpp event_visitor.cpp
	event_dispatcher.cpp functional_event_visitor.cpp instrument_layout.cpp 
	nectarcam_layout.cpp)
if(CTA_CAMERASTOACTL_FOUND)
	set(CALIN_SOURCES ${CALIN_SOURCES} zfits_actl_data_source.cpp
		zfits_data_source.cpp nectarcam_data_source.cpp)
endif()

add_library(${CALIN_TARGET_LIBRARY} SHARED ${CALIN_SOURCES})
target_link_libraries(${CALIN_TARGET_LIBRARY}
	calin_proto_io calin_proto_iact_data calin_io calin_math)

if(CTA_CAMERASTOACTL_FOUND)
	include_directories(${CTA_CAMERASTOACTL_INCLUDE_DIRS})
	target_link_libraries(${CALIN_TARGET_LIBRARY} ${CTA_CAMERASTOACTL_LIBRARIES})
endif()

install(TARGETS ${CALIN_TARGET_LIBRARY} DESTINATION ${CALIN_LIB_INSTALL_DIR})
