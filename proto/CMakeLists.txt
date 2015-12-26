set(CALIN_SWIG_PROTO common_types.proto unittest.proto)

PROTOBUF_GENERATE_CPP(ProtoSources ProtoHeaders
					 calin.proto ${CALIN_SWIG_PROTO})

message(STATUS "PROTOBUF_IMPORT_DIRS: " ${PROTOBUF_IMPORT_DIRS})
message(STATUS "PROTOBUF_INCLUDE_DIRS: " ${PROTOBUF_INCLUDE_DIRS})

set(CALIN_TARGET_LIBRARY calin_proto)
set(CALIN_SWIG_INSTALL_DIR python/calin/ix)

add_library(${CALIN_TARGET_LIBRARY} SHARED ${ProtoSources})
include_directories(${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(${CALIN_TARGET_LIBRARY} ${PROTOBUF_LIBRARY})
install(TARGETS ${CALIN_TARGET_LIBRARY} DESTINATION lib)
install(FILES calin.proto DESTINATION proto)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/calin.pb.h DESTINATION proto/${CALIN_PACKAGE})
install(FILES __init__.py DESTINATION ${CALIN_SWIG_INSTALL_DIR})

######################################################################
# GENERATE SWIG INTERFACE FILES
######################################################################

foreach (proto_interface_file ${CALIN_SWIG_PROTO})
	string(REPLACE ".proto" "" output_file ${proto_interface_file})
        set(interface_file ${output_file}.pb.i)
	add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${interface_file}
				  COMMAND protoc --plugin=protoc-gen-swig=${CMAKE_BINARY_DIR}/tools/protoc-gen-swig/protoc-gen-swig --swig_out=${CMAKE_BINARY_DIR}/proto -I${PROTOBUF_IMPORT_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/${proto_interface_file}
		  		  DEPENDS ${proto_interface_file} ${CMAKE_BINARY_DIR}/tools/protoc-gen-swig/protoc-gen-swig)
	add_custom_command(OUTPUT ${output_file}.pb_wrap.cxx
				  COMMAND ${SWIG_EXECUTABLE} ${CMAKE_SWIG_FLAGS} -interface _${output_file}_pb  -c++ -python ${interface_file}
				  DEPENDS ${interface_file})
	add_library(_${output_file}_pb MODULE ${output_file}.pb_wrap.cxx)
	set_target_properties(_${output_file}_pb PROPERTIES NO_SONAME ON)
	set_target_properties(_${output_file}_pb PROPERTIES PREFIX "")
	target_link_libraries(_${output_file}_pb calin_proto ${CALIN_TARGET_LIBRARY} ${PYTHON_LIBRARIES} ${PROTOBUF_LIBRARY})
	install(TARGETS _${output_file}_pb  DESTINATION ${CALIN_SWIG_INSTALL_DIR})
	install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${output_file}.py  DESTINATION ${CALIN_SWIG_INSTALL_DIR})
	install(FILES ${proto_interface_file} DESTINATION proto)
	install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${output_file}.pb.h DESTINATION proto/${CALIN_PACKAGE})
endforeach (proto_interface_file)

######################################################################
# SUB DIRECTORIES
######################################################################

add_subdirectory(math)
add_subdirectory(simulation)
add_subdirectory(io)
add_subdirectory(iact)
