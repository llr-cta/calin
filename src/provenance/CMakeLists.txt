# calin/src/provenance/CMakeLists.txt -- Stephen Fegan
#
# Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
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

set(CALIN_TARGET_LIBRARY calin_provenance)
set(CALIN_SOURCES system_info.cpp)

add_library(${CALIN_TARGET_LIBRARY} SHARED ${CALIN_SOURCES})
target_link_libraries(${CALIN_TARGET_LIBRARY} ${PROTOBUF_LIBRARY}
  calin_proto_provenance calin_proto)

install(TARGETS ${CALIN_TARGET_LIBRARY} DESTINATION lib)
