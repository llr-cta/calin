# calin/proto/iact_data/CMakeLists.txt -- Stephen Fegan
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

set(CALIN_SWIG_PROTO instrument_layout.proto telescope_event.proto
	zfits_data_source.proto nectarcam_data_source.proto)
set(CALIN_PACKAGE iact_data)
set(CALIN_TARGET_LIBRARY calin_proto_${CALIN_PACKAGE})

PROTOBUF_GENERATE_CPP(ProtoSources ProtoHeaders ${CALIN_SWIG_PROTO})
add_library(${CALIN_TARGET_LIBRARY} SHARED ${ProtoSources})
include_directories(${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(${CALIN_TARGET_LIBRARY} calin_proto ${PROTOBUF_LIBRARY})
install(TARGETS ${CALIN_TARGET_LIBRARY} DESTINATION lib)
install(FILES __init__.py DESTINATION ${CALIN_SWIG_INSTALL_DIR}/${CALIN_PACKAGE})
make_proto_swig_modules(${CALIN_TARGET_LIBRARY}
		${CALIN_SWIG_INSTALL_DIR}/${CALIN_PACKAGE}
		${CALIN_PROTO_INSTALL_DIR}/${CALIN_PACKAGE} ${CALIN_SWIG_PROTO})
