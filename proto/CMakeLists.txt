set(CALIN_SWIG_PROTO common_types.proto unittest.proto)

PROTOBUF_GENERATE_CPP(ProtoSources ProtoHeaders
					 calin.proto ${CALIN_SWIG_PROTO})

message(STATUS "PROTOBUF_IMPORT_DIRS: " ${PROTOBUF_IMPORT_DIRS}) 
message(STATUS "PROTOBUF_INCLUDE_DIRS: " ${PROTOBUF_INCLUDE_DIRS}) 

set(CALIN_TARGET_LIBRARY calin_proto)

add_library(${CALIN_TARGET_LIBRARY} SHARED ${ProtoSources}) 
include_directories(${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(${CALIN_TARGET_LIBRARY} ${PROTOBUF_LIBRARY}) 
install(TARGETS ${CALIN_TARGET_LIBRARY} DESTINATION lib)

######################################################################
# GENERATE SWIG INTERFACE FILES
###################################################################### 

foreach (proto_interface_file ${CALIN_SWIG_PROTO})
	string(REPLACE ".proto" "" output_file ${proto_interface_file})
        set(interface_file ${output_file}.pb.i) 
	add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${interface_file}
				  COMMAND protoc --plugin=protoc-gen-swig=${CMAKE_BINARY_DIR}/tools/protoc-gen-swig/protoc-gen-swig --swig_out=${CMAKE_CURRENT_BINARY_DIR} -I${PROTOBUF_IMPORT_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/${proto_interface_file}
		  		  DEPENDS ${proto_interface_file} ${CMAKE_BINARY_DIR}/tools/protoc-gen-swig/protoc-gen-swig)
       set_source_files_properties(${interface_file} PROPERTIES CPLUSPLUS ON) 
       swig_add_module(${output_file} python ${CMAKE_CURRENT_BINARY_DIR}/${interface_file})
      swig_link_libraries(${output_file} ${CALIN_TARGET_LIBRARY} ${PYTHON_LIBRARIES} ${PROTOBUF_LIBRARIES})
     install(TARGETS ${SWIG_MODULE_${output_file}_REAL_NAME}  DESTINATION python/calin/ix)
        install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${output_file}.py  DESTINATION python/calin/ix)
endforeach (proto_interface_file)


add_subdirectory(math)
add_subdirectory(simulation)
add_subdirectory(io)
