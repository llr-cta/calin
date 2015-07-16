find_package(Geant4 REQUIRED)

include(${Geant4_USE_FILE})

set(CALIN_TARGET_LIBRARY calin_air_shower)
set(CALIN_SOURCES atmosphere.cpp geant4_shower_generator.cpp
		  geant4_shower_generator_internals.cpp)

add_library(${CALIN_TARGET_LIBRARY} SHARED ${CALIN_SOURCES})
target_link_libraries(${CALIN_TARGET_LIBRARY} calin_io ${Geant4_LIBRARIES})

install(TARGETS ${CALIN_TARGET_LIBRARY} DESTINATION lib)
