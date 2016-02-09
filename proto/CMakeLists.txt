# calin/proto/CMakeLists.txt -- Stephen Fegan
#
# Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
# LLR, Ecole polytechnique, CNRS/IN2P3, Universite Paris-Saclay
#
# This file is part of "calin"
#
# "calin" is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License version 2 or later, as published by
# the Free Software Foundation.
#
# "calin" is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

set(CALIN_SWIG_PROTO common_types.proto unittest.proto)
set(CALIN_TARGET_LIBRARY calin_proto)
set(CALIN_SWIG_INSTALL_DIR python/calin/ix)
set(CALIN_PROTO_INSTALL_DIR proto)

message(STATUS "PROTOBUF_IMPORT_DIRS: " ${PROTOBUF_IMPORT_DIRS})
message(STATUS "PROTOBUF_INCLUDE_DIRS: " ${PROTOBUF_INCLUDE_DIRS})

PROTOBUF_GENERATE_CPP(ProtoSources ProtoHeaders
					 calin.proto ${CALIN_SWIG_PROTO})

add_library(${CALIN_TARGET_LIBRARY} SHARED ${ProtoSources})
include_directories(${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(${CALIN_TARGET_LIBRARY} ${PROTOBUF_LIBRARY})
install(TARGETS ${CALIN_TARGET_LIBRARY} DESTINATION lib)
install(FILES calin.proto DESTINATION proto)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/calin.pb.h DESTINATION proto/${CALIN_PACKAGE})
install(FILES __init__.py DESTINATION ${CALIN_SWIG_INSTALL_DIR})

######################################################################
# MACRO TO GENERATE SWIG INTERFACE FILES
######################################################################

macro(make_proto_modules arg_proto_library arg_proto_install_dir)
	set(proto_cxx "")
	foreach (proto_interface_file ${ARGN})
		string(REPLACE ".proto" "" base_name ${proto_interface_file})
		add_custom_command(OUTPUT ${base_name}.pb.cc
			COMMAND ${PROTOBUF_PROTOC_EXECUTABLE} -I${PROTOBUF_IMPORT_DIRS} --cpp_out=${CMAKE_BINARY_DIR}/proto ${CMAKE_CURRENT_SOURCE_DIR}/${proto_interface_file}
			DEPENDS ${proto_interface_file})
		set(proto_cxx ${proto_cxx} ${base_name}.pb.cc)
	endforeach (proto_interface_file)
	add_library(${arg_proto_library} SHARED ${proto_cxx})
	include_directories(${CMAKE_CURRENT_BINARY_DIR})
	target_link_libraries(${arg_proto_library} calin_proto ${PROTOBUF_LIBRARY})
	install(TARGETS ${arg_proto_library} DESTINATION lib)
	install(FILES ${proto_interface_file} DESTINATION ${arg_proto_install_dir})
	install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${base_name}.pb.h DESTINATION ${arg_proto_install_dir})
endmacro()

macro(make_proto_swig_modules arg_proto_library arg_swig_install_dir arg_proto_install_dir)
	foreach (proto_interface_file ${ARGN})
		string(REPLACE ".proto" "" base_name ${proto_interface_file})
  	set(interface_file ${base_name}.pb.i)
		add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${interface_file}
			COMMAND ${PROTOBUF_PROTOC_EXECUTABLE} --plugin=protoc-gen-swig=${CMAKE_BINARY_DIR}/tools/protoc-gen-swig/protoc-gen-swig --swig_out=${CMAKE_BINARY_DIR}/proto -I${PROTOBUF_IMPORT_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/${proto_interface_file}
		  DEPENDS ${proto_interface_file} ${CMAKE_BINARY_DIR}/tools/protoc-gen-swig/protoc-gen-swig)
		add_custom_command(OUTPUT ${base_name}.pb_wrap.cxx
		  COMMAND ${SWIG_EXECUTABLE} ${CMAKE_SWIG_FLAGS} -interface _${base_name}_pb  -c++ -python ${interface_file}
		  DEPENDS ${interface_file})
		include_directories(${CMAKE_CURRENT_BINARY_DIR})
		add_library(_${base_name}_pb MODULE ${base_name}.pb_wrap.cxx)
		set_target_properties(_${base_name}_pb PROPERTIES NO_SONAME ON)
		set_target_properties(_${base_name}_pb PROPERTIES PREFIX "")
		target_link_libraries(_${base_name}_pb calin_proto ${arg_proto_library} ${PYTHON_LIBRARIES} ${PROTOBUF_LIBRARY})
		install(TARGETS _${base_name}_pb  DESTINATION ${arg_swig_install_dir})
		install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${base_name}.py  DESTINATION ${arg_swig_install_dir})
		install(FILES ${proto_interface_file} DESTINATION ${arg_proto_install_dir})
		install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${base_name}.pb.h DESTINATION ${arg_proto_install_dir})
	endforeach (proto_interface_file)
endmacro()

######################################################################
# GENERATE SWIG INTERFACE FILES
######################################################################

make_proto_swig_modules(${CALIN_TARGET_LIBRARY} ${CALIN_SWIG_INSTALL_DIR}
		${CALIN_PROTO_INSTALL_DIR} ${CALIN_SWIG_PROTO})

######################################################################
# SUB DIRECTORIES
######################################################################

add_subdirectory(math)
add_subdirectory(simulation)
add_subdirectory(io)
add_subdirectory(iact_data)
add_subdirectory(diagnostics)
