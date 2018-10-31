# calin/swig/CMakeLists.txt -- Stephen Fegan
#
# Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
# LLR, Ecole Polytechnique, CNRS/IN2P3
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

INCLUDE (CheckCXXSourceRuns)

CHECK_CXX_SOURCE_RUNS("#include<cstdint>
#include<cstdlib>
template<typename T> int is_long_int() { return EXIT_FAILURE; }
template<> int is_long_int<long int>() { return EXIT_SUCCESS; }
int main() { return is_long_int<int64_t>(); }" CALIN_INT64_IS_LONGINT)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/calin_config.i.in
  ${CMAKE_CURRENT_BINARY_DIR}/calin_config.i)

set(CALIN_SWIG_FILES google_protobuf.i)
#set(CALIN_SWIG_FILES calin_global_definitions.i google_protobuf.i)
set(CALIN_SWIG_INSTALL_DIR ${CALIN_PYTHON_INSTALL_DIR})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
#include_directories(${CMAKE_CURRENT_BINARY_DIR}) # For calin_config.i

execute_process(COMMAND ${CMAKE_CXX_COMPILER} ${CXX_FLAGS_LIST}
	-o ${CMAKE_BINARY_DIR}/gen-swig_dep
	${CMAKE_SOURCE_DIR}/tools/gen-swig_dep/gen-swig_dep.cpp)

macro(make_swig_modules arg_calin_library arg_swig_install_dir)
  foreach (interface_file ${ARGN})
    string(REPLACE ".i" "" base_name ${interface_file})
		set(interface_dep ${CMAKE_CURRENT_BINARY_DIR}/${base_name}.swig_dep)
    execute_process(COMMAND ${CMAKE_BINARY_DIR}/gen-swig_dep
			${CMAKE_CURRENT_SOURCE_DIR}/${interface_file} ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR})
		file(READ ${interface_dep} SWIG_DEPS)
		string(STRIP "${SWIG_DEPS}" SWIG_DEPS)
		add_custom_command(OUTPUT ${base_name}_wrap.cxx
		  COMMAND ${SWIG_EXECUTABLE} ${CMAKE_SWIG_FLAGS} -outcurrentdir -interface _${base_name} -c++ -python ${CMAKE_CURRENT_SOURCE_DIR}/${interface_file}
			DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/${interface_file};${SWIG_DEPS}"
      COMMENT "swig: ${interface_file}")
		include_directories(${CMAKE_CURRENT_BINARY_DIR})
		add_library(_${base_name} MODULE ${base_name}_wrap.cxx)
		set_target_properties(_${base_name} PROPERTIES NO_SONAME ON)
		set_target_properties(_${base_name} PROPERTIES PREFIX "")
		target_link_libraries(_${base_name} calin_proto ${arg_calin_library} ${PYTHON_LIBRARIES} ${PROTOBUF_LIBRARY})
		install(TARGETS _${base_name}  DESTINATION ${arg_swig_install_dir})
		install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${base_name}.py DESTINATION ${arg_swig_install_dir})
  endforeach (interface_file)

  set(CALIN_PY_INIT_SWIG_LIST FALSE)
  foreach (interface_file ${ARGN})
    string(REPLACE ".i" "" base_name ${interface_file})
    if(CALIN_PY_INIT_SWIG_LIST)
      set(CALIN_PY_INIT_SWIG_LIST "${CALIN_PY_INIT_SWIG_LIST}, \"${base_name}\"")
    else(CALIN_PY_INIT_SWIG_LIST)
      set(CALIN_PY_INIT_SWIG_LIST "\"${base_name}\"")
    endif(CALIN_PY_INIT_SWIG_LIST)
  endforeach (interface_file)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/__init__.py.in
    ${CMAKE_CURRENT_BINARY_DIR}/__init__.py)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/__init__.py
    DESTINATION ${arg_swig_install_dir})
endmacro()

make_swig_modules("" ${CALIN_SWIG_INSTALL_DIR} ${CALIN_SWIG_FILES})

add_subdirectory(math)
add_subdirectory(io)
add_subdirectory(iact_data)
add_subdirectory(diagnostics)
add_subdirectory(calib)
add_subdirectory(simulation)
add_subdirectory(provenance)
add_subdirectory(util)
