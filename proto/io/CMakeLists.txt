set(CALIN_PACKAGE io)
set(CALIN_TARGET_LIBRARY calin_proto_${CALIN_PACKAGE})

PROTOBUF_GENERATE_CPP(ProtoSources ProtoHeaders
					 sql_transceiver.proto)
add_library(${CALIN_TARGET_LIBRARY} SHARED ${ProtoSources})
include_directories(${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(${CALIN_TARGET_LIBRARY} calin_proto ${PROTOBUF_LIBRARY})
install(TARGETS ${CALIN_TARGET_LIBRARY} DESTINATION lib)

install(FILES sql_transceiver.proto
	DESTINATION ${CALIN_PROTO_INSTALL_DIR}/${CALIN_PACKAGE})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/sql_transceiver.pb.h
	DESTINATION ${CALIN_PROTO_INSTALL_DIR}/${CALIN_PACKAGE})
