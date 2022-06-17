# calin/tools/protoc-gen-calin-cpp-extensions/CMakeLists.txt -- Stephen Fegan
#																														2017-03-14
#
# Copyright 2017, Stephen Fegan <sfegan@llr.in2p3.fr>
# Laboratoire Leprince-Ringuet, CNRS/IN2P3, Ecole Polytechnique, Institut Polytechnique de Paris
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

add_executable(protoc-gen-calin-cpp-extensions protoc-gen-calin-cpp-extensions.cpp
	calin_cpp_extensions_generator.cpp)

include_directories(${CMAKE_CURRENT_BUILD_DIR})

target_link_libraries(protoc-gen-calin-cpp-extensions calin_util_lite
	${PROTOBUF_LIBRARY} calin_proto_base ${PROTOBUF_PROTOC_LIBRARIES} -lpthread)

target_compile_options(protoc-gen-calin-cpp-extensions PUBLIC -march=)

install(TARGETS protoc-gen-calin-cpp-extensions DESTINATION ${CALIN_BIN_INSTALL_DIR})
