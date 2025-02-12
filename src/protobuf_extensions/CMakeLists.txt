# calin/src/proto/CMakeLists.txt -- Stephen Fegan
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

set(CALIN_TARGET_LIBRARY protobuf_extensions)
set(CALIN_SOURCES hdf_stream_reader.cpp hdf_stream_writer.cpp)

add_library(${CALIN_TARGET_LIBRARY} SHARED ${CALIN_SOURCES})
target_link_libraries(${CALIN_TARGET_LIBRARY} ${PROTOBUF_LIBRARY} calin_proto
  ${HDF5_LIBRARIES})

install(TARGETS ${CALIN_TARGET_LIBRARY} DESTINATION ${CALIN_LIB_INSTALL_DIR})
