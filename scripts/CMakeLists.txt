# calin/scripts/CMakeLists.txt -- Stephen Fegan
#
# Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
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

install(PROGRAMS cta_stage1.py DESTINATION ${CALIN_BIN_INSTALL_DIR})
install(PROGRAMS render_stage1_results.py DESTINATION ${CALIN_BIN_INSTALL_DIR})
install(PROGRAMS scriptSPEFit.py DESTINATION ${CALIN_BIN_INSTALL_DIR})
