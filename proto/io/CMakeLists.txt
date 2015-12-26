PROTOBUF_GENERATE_CPP(ProtoSources ProtoHeaders
					 sql_transceiver.proto)

set(CALIN_PACKAGE io)
set(CALIN_TARGET_LIBRARY calin_proto_${CALIN_PACKAGE})
set(CALIN_SWIG_INSTALL_DIR python/calin/ix/${CALIN_PACKAGE})

add_library(${CALIN_TARGET_LIBRARY} SHARED ${ProtoSources})
include_directories(${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(${CALIN_TARGET_LIBRARY} calin_proto ${PROTOBUF_LIBRARY})
install(TARGETS ${CALIN_TARGET_LIBRARY} DESTINATION lib)
install(FILES sql_transceiver.proto DESTINATION proto/${CALIN_PACKAGE})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/sql_transceiver.pb.h DESTINATION proto/${CALIN_PACKAGE})
