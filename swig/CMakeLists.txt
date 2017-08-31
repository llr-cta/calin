# calin/swig/CMakeLists.txt -- Stephen Fegan
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

INCLUDE (CheckCXXSourceRuns)

CHECK_CXX_SOURCE_RUNS("#include<cstdint>
#include<cstdlib>
template<typename T> int is_long_int() { return EXIT_FAILURE; }
template<> int is_long_int<long int>() { return EXIT_SUCCESS; }
int main() { return is_long_int<int64_t>(); }" CALIN_INT64_IS_LONGINT)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/calin_config.i.in
  ${CMAKE_CURRENT_BINARY_DIR}/calin_config.i)

set(CALIN_SWIG_FILES calin_global_definitions.i google_protobuf.i)
set(CALIN_SWIG_INSTALL_DIR ${CALIN_PYTHON_INSTALL_DIR})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
#include_directories(${CMAKE_CURRENT_BINARY_DIR}) # For calin_config.i

execute_process(COMMAND ${CMAKE_CXX_COMPILER} ${CXX_FLAGS_LIST}
	-o ${CMAKE_BINARY_DIR}/gen-swig_dep
	${CMAKE_SOURCE_DIR}/tools/gen-swig_dep/gen-swig_dep.cpp)

macro(make_swig_modules arg_calin_library arg_swig_install_dir)
  foreach (interface_file ${ARGN})
    set_source_files_properties(${interface_file} PROPERTIES CPLUSPLUS ON)
    string(REPLACE ".i" "" base_name ${interface_file})
    set(interface_dep ${CMAKE_CURRENT_BINARY_DIR}/${base_name}.swig_dep)
    execute_process(COMMAND ${CMAKE_BINARY_DIR}/gen-swig_dep
			${CMAKE_CURRENT_SOURCE_DIR}/${interface_file} ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR})
		file(READ ${interface_dep} SWIG_MODULE_${base_name}_EXTRA_DEPS)
    # This could eventually be removed by advancing the minimum required CMake version
    if (${CMAKE_VERSION} VERSION_LESS "3.8.0")
      swig_add_module(${base_name} python ${interface_file})
    else()
      swig_add_library(${base_name} LANGUAGE python TYPE MODULE SOURCES ${interface_file})
    endif()
    swig_link_libraries(${base_name} ${arg_calin_library} ${PYTHON_LIBRARIES}
      ${PROTOBUF_LIBRARIES})
    install(TARGETS ${SWIG_MODULE_${base_name}_REAL_NAME}
      DESTINATION ${arg_swig_install_dir})
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${base_name}.py
      DESTINATION ${arg_swig_install_dir})
  endforeach (interface_file)
endmacro()

make_swig_modules("" ${CALIN_SWIG_INSTALL_DIR} ${CALIN_SWIG_FILES})
install(FILES __init__.py DESTINATION ${CALIN_SWIG_INSTALL_DIR})

add_subdirectory(math)
add_subdirectory(io)
add_subdirectory(iact_data)
add_subdirectory(diagnostics)
add_subdirectory(calib)
add_subdirectory(simulation)
add_subdirectory(provenance)
add_subdirectory(util)
