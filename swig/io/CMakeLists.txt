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

set(CALIN_SWIG_FILES sql_serializer.i zmq_inproc_push_pull.i data_source.i
	options_processor.i)
set(CALIN_SWIG_TARGET io)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
make_swig_modules(calin_${CALIN_SWIG_TARGET}
	${CALIN_SWIG_INSTALL_DIR}/${CALIN_SWIG_TARGET} ${CALIN_SWIG_FILES})
