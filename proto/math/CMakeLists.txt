PROTOBUF_GENERATE_CPP(ProtoSources ProtoHeaders
					 histogram.proto rng.proto)

set(CALIN_TARGET_LIBRARY calin_proto_math)

add_library(${CALIN_TARGET_LIBRARY} SHARED ${ProtoSources}) 
include_directories(${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(${CALIN_TARGET_LIBRARY} calin_proto ${PROTOBUF_LIBRARY}) 
install(TARGETS ${CALIN_TARGET_LIBRARY} DESTINATION lib)

