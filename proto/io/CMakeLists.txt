PROTOBUF_GENERATE_CPP(ProtoSources ProtoHeaders
					 sql_transceiver.proto)

set(CALIN_TARGET_LIBRARY calin_proto_io)

add_library(${CALIN_TARGET_LIBRARY} SHARED ${ProtoSources}) 
include_directories(${CMAKE_CURRENT_BINARY_DIR}) 
target_link_libraries(${CALIN_TARGET_LIBRARY} calin_proto ${PROTOBUF_LIBRARY}) 
install(TARGETS ${CALIN_TARGET_LIBRARY} DESTINATION lib)

