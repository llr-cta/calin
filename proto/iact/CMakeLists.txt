set(CALIN_SWIG_PROTO instrument_layout.proto telescope_event.proto)
set(CALIN_PACKAGE iact)
set(CALIN_TARGET_LIBRARY calin_proto_${CALIN_PACKAGE})

PROTOBUF_GENERATE_CPP(ProtoSources ProtoHeaders ${CALIN_SWIG_PROTO})
add_library(${CALIN_TARGET_LIBRARY} SHARED ${ProtoSources})
include_directories(${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(${CALIN_TARGET_LIBRARY} calin_proto ${PROTOBUF_LIBRARY})
install(TARGETS ${CALIN_TARGET_LIBRARY} DESTINATION lib)
install(FILES __init__.py DESTINATION ${CALIN_SWIG_INSTALL_DIR}/${CALIN_PACKAGE})
make_proto_swig_modules(${CALIN_TARGET_LIBRARY}
		${CALIN_SWIG_INSTALL_DIR}/${CALIN_PACKAGE}
		${CALIN_PROTO_INSTALL_DIR}/${CALIN_PACKAGE} ${CALIN_SWIG_PROTO})
