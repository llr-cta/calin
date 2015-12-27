# calin/src/io/CMakeLists.txt -- Stephen Fegan
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

set(CALIN_TARGET_LIBRARY calin_io)
set(CALIN_SOURCES log.cpp sql_statement.cpp sql_transceiver.cpp
		  sqlite3_statement.cpp sqlite3_transceiver.cpp util.cpp)

include_directories(${PYTHON_INCLUDE_PATH})
add_library(${CALIN_TARGET_LIBRARY} SHARED ${CALIN_SOURCES})
target_link_libraries(${CALIN_TARGET_LIBRARY} calin_proto_io
${SQLITE3_LIBRARY} ${PYTHON_LIBRARIES})

install(TARGETS ${CALIN_TARGET_LIBRARY} DESTINATION lib)
