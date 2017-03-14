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
set(CALIN_PROTO_SWIG_INSTALL_DIR ${CALIN_PYTHON_INSTALL_DIR}/ix)

message(STATUS "PROTOBUF_IMPORT_DIRS: " ${PROTOBUF_IMPORT_DIRS})
message(STATUS "PROTOBUF_INCLUDE_DIRS: " ${PROTOBUF_INCLUDE_DIRS})
message(STATUS "PROTOBUF_LIBRARY: " ${PROTOBUF_LIBRARY})

######################################################################
# MACRO TO GENERATE PROTOBUF C++ AND SWIG FILES AND COMPILE THEM
######################################################################

string(REPLACE " " ";" CXX_FLAGS_LIST ${CMAKE_CXX_FLAGS})

execute_process(COMMAND ${CMAKE_CXX_COMPILER} ${CXX_FLAGS_LIST}
	-I ${PROTOBUF_INCLUDE_DIRS}
	-o ${CMAKE_BINARY_DIR}/protoc-gen-swig_dep
	${CMAKE_SOURCE_DIR}/tools/protoc-gen-swig_dep/protoc-gen-swig_dep.cpp
	${PROTOBUF_LIBRARY} ${PROTOBUF_PROTOC_LIBRARY})

macro(make_proto_modules arg_proto_library arg_package)
	if("${arg_package}" STREQUAL "")
		set(proto_install_dir ${CALIN_PROTO_INSTALL_DIR})
		set(proto_header_install_dir ${CALIN_PROTO_HEADER_INSTALL_DIR})
	else()
		set(proto_install_dir ${CALIN_PROTO_INSTALL_DIR}/${arg_package})
		set(proto_header_install_dir ${CALIN_PROTO_HEADER_INSTALL_DIR}/${arg_package})
	endif()
	set(proto_cxx "")
	set(proto_proto "")
	foreach (proto_interface_file ${ARGN})
		get_filename_component(file_ext ${proto_interface_file} EXT)
		if("${file_ext}" STREQUAL ".proto")
			string(REPLACE ".proto" "" base_name ${proto_interface_file})
			add_custom_command(OUTPUT ${base_name}.pb.cc ${base_name}.pb.h
				COMMAND ${PROTOBUF_PROTOC_EXECUTABLE} -I${PROTOBUF_IMPORT_DIRS} --cpp_out=${CMAKE_BINARY_DIR}/proto --plugin=protoc-gen-calin-cpp-extensions=${CMAKE_BINARY_DIR}/tools/protoc-gen-calin-cpp-extensions/protoc-gen-calin-cpp-extensions --calin-cpp-extensions_out=${CMAKE_BINARY_DIR}/proto ${CMAKE_CURRENT_SOURCE_DIR}/${proto_interface_file}
				DEPENDS ${proto_interface_file} ${CMAKE_BINARY_DIR}/tools/protoc-gen-calin-cpp-extensions/protoc-gen-calin-cpp-extensions)
			set(proto_cxx ${proto_cxx} ${base_name}.pb.cc)
			set(proto_proto ${proto_proto} ${proto_interface_file})
			install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${base_name}.pb.h DESTINATION ${proto_header_install_dir})
		else()
			set(proto_cxx ${proto_cxx} ${proto_interface_file})
		endif()
	endforeach (proto_interface_file)
	add_library(${arg_proto_library} SHARED ${proto_cxx})
	include_directories(${CMAKE_CURRENT_BINARY_DIR})
	if(${arg_proto_library} STREQUAL "calin_proto")
		target_link_libraries(${arg_proto_library} calin_proto_base ${PROTOBUF_LIBRARY})
	else()
		target_link_libraries(${arg_proto_library} calin_proto ${PROTOBUF_LIBRARY})
	endif()
	install(TARGETS ${arg_proto_library} DESTINATION ${CALIN_LIB_INSTALL_DIR})
	install(FILES ${proto_proto} DESTINATION ${proto_install_dir})
endmacro()

macro(make_proto_swig_modules arg_proto_library arg_swig_install_dir)
	foreach (proto_interface_file ${ARGN})
		string(REPLACE ".proto" "" base_name ${proto_interface_file})
  	set(interface_file ${base_name}.pb.i)
		add_custom_command(OUTPUT ${interface_file}
			COMMAND ${PROTOBUF_PROTOC_EXECUTABLE} --plugin=protoc-gen-swig=${CMAKE_BINARY_DIR}/tools/protoc-gen-swig/protoc-gen-swig --swig_out=${CMAKE_BINARY_DIR}/proto -I${PROTOBUF_IMPORT_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/${proto_interface_file}
		  DEPENDS ${proto_interface_file} ${CMAKE_BINARY_DIR}/tools/protoc-gen-swig/protoc-gen-swig)
	endforeach (proto_interface_file)
	foreach (proto_interface_file ${ARGN})
		string(REPLACE ".proto" "" base_name ${proto_interface_file})
		set(interface_file ${base_name}.pb.i)
		set(interface_dep ${CMAKE_CURRENT_BINARY_DIR}/${base_name}.pb.swig_dep)
		execute_process(COMMAND ${PROTOBUF_PROTOC_EXECUTABLE}
			--plugin=protoc-gen-swig_dep=${CMAKE_BINARY_DIR}/protoc-gen-swig_dep
			--swig_dep_out=${CMAKE_SOURCE_DIR},${CMAKE_BINARY_DIR}:${CMAKE_BINARY_DIR}/proto -I${PROTOBUF_IMPORT_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/${proto_interface_file})
		file(READ ${interface_dep} SWIG_DEPS)
		string(STRIP "${SWIG_DEPS}" SWIG_DEPS)
#		string(REGEX REPLACE "([^;]+)" "${CMAKE_BINARY_DIR}/proto/\\1" SWIG_DEPS "${SWIG_DEPS}")
#		MESSAGE(STATUS "${SWIG_DEPS}")
		add_custom_command(OUTPUT ${base_name}.pb_wrap.cxx
		  COMMAND ${SWIG_EXECUTABLE} ${CMAKE_SWIG_FLAGS} -interface _${base_name}_pb  -c++ -python ${interface_file}
			DEPENDS "${interface_file};${CMAKE_BINARY_DIR}/swig/calin_config.i;${SWIG_DEPS}")
		include_directories(${CMAKE_CURRENT_BINARY_DIR})
		add_library(_${base_name}_pb MODULE ${base_name}.pb_wrap.cxx)
		set_target_properties(_${base_name}_pb PROPERTIES NO_SONAME ON)
		set_target_properties(_${base_name}_pb PROPERTIES PREFIX "")
		target_link_libraries(_${base_name}_pb calin_proto ${arg_proto_library} ${PYTHON_LIBRARIES} ${PROTOBUF_LIBRARY})
		install(TARGETS _${base_name}_pb  DESTINATION ${arg_swig_install_dir})
		install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${base_name}.py  DESTINATION ${arg_swig_install_dir})
		install(FILES __init__.py DESTINATION ${arg_swig_install_dir})
	endforeach (proto_interface_file)
endmacro()

######################################################################
# GENERATE SWIG INTERFACE FILES
######################################################################

add_custom_command(OUTPUT calin.pb.cc calin.pb.h
	COMMAND ${PROTOBUF_PROTOC_EXECUTABLE} -I${PROTOBUF_IMPORT_DIRS} --cpp_out=${CMAKE_BINARY_DIR}/proto ${CMAKE_CURRENT_SOURCE_DIR}/calin.proto
	DEPENDS calin.proto)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/calin.pb.h DESTINATION ${CALIN_PROTO_HEADER_INSTALL_DIR})
install(FILES calin.proto DESTINATION ${CALIN_PROTO_INSTALL_DIR})
add_library(calin_proto_base SHARED calin.pb.cc)
include_directories(${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(calin_proto_base ${PROTOBUF_LIBRARY})
install(TARGETS calin_proto_base DESTINATION ${CALIN_LIB_INSTALL_DIR})

make_proto_modules(${CALIN_TARGET_LIBRARY} ""
	${CALIN_SWIG_PROTO})
make_proto_swig_modules(${CALIN_TARGET_LIBRARY} ${CALIN_PROTO_SWIG_INSTALL_DIR}
	${CALIN_SWIG_PROTO})

######################################################################
# SUB DIRECTORIES
######################################################################

add_subdirectory(math)
add_subdirectory(simulation)
add_subdirectory(io)
add_subdirectory(provenance)
add_subdirectory(iact_data)
add_subdirectory(diagnostics)
add_subdirectory(calib)
add_subdirectory(scripts)
