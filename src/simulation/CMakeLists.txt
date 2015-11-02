set(CALIN_TARGET_LIBRARY calin_simulation)
set(CALIN_SOURCES atmosphere.cpp)

if(GEANT4_FOUND)
	set(CALIN_SOURCES ${CALIN_SOURCES} tracker.cpp
			  geant4_shower_generator.cpp
		  	  geant4_shower_generator_internals.cpp)
endif()

add_library(${CALIN_TARGET_LIBRARY} SHARED ${CALIN_SOURCES})
target_link_libraries(${CALIN_TARGET_LIBRARY} calin_io)

if(GEANT4_FOUND)
	target_link_libraries(${CALIN_TARGET_LIBRARY} ${GEANT4_LIBRARIES})
endif()

install(TARGETS ${CALIN_TARGET_LIBRARY} DESTINATION lib)
