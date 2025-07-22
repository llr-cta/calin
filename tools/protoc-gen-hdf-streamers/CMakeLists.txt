# calin/tools/protoc-gen-hdf-streaners/CMakeLists.txt -- Stephen Fegan
#																														2017-03-14
#
# Copyright 2025, Stephen Fegan <sfegan@llr.in2p3.fr>
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

add_executable(protoc-gen-hdf-streamers protoc-gen-hdf-streamers.cpp
	hdf_streamers_generator.cpp)

include_directories(${CMAKE_CURRENT_BUILD_DIR})

target_link_libraries(protoc-gen-hdf-streamers calin_util_native
	${PROTOBUF_LIBRARY} calin_proto_base_native ${PROTOBUF_PROTOC_LIBRARIES} -lpthread)

target_compile_options(protoc-gen-hdf-streamers PUBLIC -march=native)

install(TARGETS protoc-gen-hdf-streamers DESTINATION ${CALIN_BIN_INSTALL_DIR})
