# calin/src/bin/CMakeLists.txt -- Stephen Fegan
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

add_executable(write_system_info write_system_info.cpp)
target_link_libraries(write_system_info calin_provenance calin_util)

install(TARGETS write_system_info DESTINATION ${CALIN_BIN_INSTALL_DIR})

if(CTA_CAMERASTOACTL_FOUND)
	include_directories(${CTA_CAMERASTOACTL_INCLUDE_DIRS})
	add_executable(dump_nectar_zfits dump_nectar_zfits.cpp)
  target_link_libraries(dump_nectar_zfits calin_iact_data calin_util)
  install(TARGETS dump_nectar_zfits DESTINATION ${CALIN_BIN_INSTALL_DIR})

	add_executable(write_nectar_zfits_as_raw write_nectar_zfits_as_raw.cpp)
  target_link_libraries(write_nectar_zfits_as_raw calin_iact_data calin_util)
  install(TARGETS write_nectar_zfits_as_raw DESTINATION ${CALIN_BIN_INSTALL_DIR})
endif()
