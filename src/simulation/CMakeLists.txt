# calin/src/simulation/CMakeLists.txt -- Stephen Fegan
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

set(CALIN_TARGET_LIBRARY calin_simulation)
set(CALIN_SOURCES atmosphere.cpp pmt.cpp
	vso_pixel.cpp vso_mirror.cpp vso_obscuration.cpp vso_telescope.cpp
  vso_array.cpp vso_raytracer.cpp
	pe_processor.cpp ray_processor.cpp vso_ray_processor.cpp
	air_cherenkov_tracker.cpp iact_array_tracker.cpp
	tracker.cpp misc_trackers.cpp
  straight_track_generator.cpp quadrature_iact_array_integration.cpp
 	bfield_track_generator.cpp vso_quadrature_iact_array_integration.cpp
	detector_efficiency.cpp world_magnetic_model.cpp
	GeomagnetismLibrary.c)

if(GEANT4_FOUND)
	set(CALIN_SOURCES ${CALIN_SOURCES}
		geant4_shower_generator.cpp
		geant4_shower_generator_internals.cpp)
endif()

add_library(${CALIN_TARGET_LIBRARY} SHARED ${CALIN_SOURCES})
target_link_libraries(${CALIN_TARGET_LIBRARY} calin_io calin_math
	calin_proto_simulation calin_iact_data calin_proto_iact_data calin_provenance)

if(GEANT4_FOUND)
	target_link_libraries(${CALIN_TARGET_LIBRARY} ${GEANT4_LIBRARIES})
endif()

install(TARGETS ${CALIN_TARGET_LIBRARY} DESTINATION ${CALIN_LIB_INSTALL_DIR})
