set(CALIN_TARGET_LIBRARY calin_io)
set(CALIN_SOURCES log.cpp sql_transceiver.cpp sqlite3_transceiver.cpp)

include_directories(${PYTHON_INCLUDE_PATH}) 
add_library(${CALIN_TARGET_LIBRARY} SHARED ${CALIN_SOURCES})
target_link_libraries(${CALIN_TARGET_LIBRARY} calin_proto
${SQLITE3_LIBRARY} ${PYTHON_LIBRARIES})

install(TARGETS ${CALIN_TARGET_LIBRARY} DESTINATION lib)
