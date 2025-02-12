# calin/CMakeLists.txt -- Stephen Fegan
#
# Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
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

cmake_minimum_required(VERSION 3.15)
project(calin)
enable_testing()

# ============================================================= #
# Find the various libraries and packages that calin depends on #
# ============================================================= #

## Use CMake policy to specify complete library path
cmake_policy(SET CMP0060 NEW)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/CMakeModules/")

include_directories(include ${CMAKE_BINARY_DIR}/include)

if(NOT CALIN_BUILD_ARCH)
	set(CALIN_BUILD_ARCH "native")
endif()
message(STATUS "Tuning build for architecture: " ${CALIN_BUILD_ARCH})
if(${CALIN_BUILD_ARCH} STREQUAL "generic")
	set(CALIN_BUILD_ARCH_OPT "-mtune=${CALIN_BUILD_ARCH}")
else()
	set(CALIN_BUILD_ARCH_OPT "-march=${CALIN_BUILD_ARCH}")
endif()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CALIN_BUILD_ARCH_OPT} -Wall")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CALIN_BUILD_ARCH_OPT} -std=c++17 -Wall")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} -O3")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -O3")

if(NOT CALIN_PYTHON_SUB_DIR)
	set(CALIN_PYTHON_SUB_DIR "python")
endif()

if(NOT CALIN_SHARE_SUB_DIR)
	set(CALIN_SHARE_SUB_DIR "share")
endif()

if(NOT CALIN_INCLUDE_SUB_DIR)
	set(CALIN_INCLUDE_SUB_DIR "include")
endif()

if(NOT CALIN_LIB_SUB_DIR)
	set(CALIN_LIB_SUB_DIR "lib")
endif()

if(NOT CALIN_BIN_SUB_DIR)
	set(CALIN_BIN_SUB_DIR "bin")
endif()

set(CALIN_PROTO_INSTALL_DIR ${CALIN_SHARE_SUB_DIR}/calin/proto)
set(CALIN_PROTO_HEADER_INSTALL_DIR ${CALIN_INCLUDE_SUB_DIR}/calin/proto)
set(CALIN_DATA_INSTALL_DIR ${CALIN_SHARE_SUB_DIR}/calin/data)
set(CALIN_EXAMPLES_INSTALL_DIR ${CALIN_SHARE_SUB_DIR}/calin/examples)
set(CALIN_LIB_INSTALL_DIR ${CALIN_LIB_SUB_DIR})
set(CALIN_BIN_INSTALL_DIR ${CALIN_BIN_SUB_DIR})
set(CALIN_HEADER_INSTALL_DIR ${CALIN_INCLUDE_SUB_DIR}/calin)
set(CALIN_PYTHON_INSTALL_DIR ${CALIN_PYTHON_SUB_DIR}/calin)
set(CALIN_UNIT_TEST_INSTALL_DIR ${CALIN_SHARE_SUB_DIR}/calin/unit_tests)

## Find Protobuf libraries and set protobuf base path
set(PROTOBUF_IMPORT_DIRS ${PROJECT_SOURCE_DIR}/proto)
find_package(Protobuf REQUIRED)
include_directories(${PROTOBUF_INCLUDE_DIRS} ${CMAKE_BINARY_DIR}/proto)
set(EXTRA_LIBS ${EXTRA_LIBS} ${PROTOBUF_LIBRARY})

## Find FFTW libraries
find_package(FFTW REQUIRED)
include_directories(${FFTW_INCLUDE_DIRS})
set(EXTRA_LIBS ${EXTRA_LIBS} ${FFTW_LIBRARIES} )

## Find GSL libraries
find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
set(EXTRA_LIBS ${EXTRA_LIBS} ${GSL_LIBRARIES} )

## Find Sqlite3 libraries
find_package(SQLite3 REQUIRED)
include_directories(${SQLITE3_INCLUDE_DIR})
set(EXTRA_LIBS ${EXTRA_LIBS} ${SQLITE3_LIBRARY} )

## Find HDF5
find_package(HDF5 REQUIRED)
include_directories(${HDF5_INCLUDE_DIRS})

## Find Python interpreter and libraries
find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)
include_directories(${Python3_INCLUDE_DIRS} ${Python3_NumPy_INCLUDE_DIRS})
message(STATUS "Python3_LIBRARIES:" ${Python3_LIBRARIES})

## Find SWIG and set SWIG include paths to "swig", "include" and "proto"
find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})
set(CMAKE_SWIG_FLAGS -I${CMAKE_SOURCE_DIR}/swig -I${CMAKE_BINARY_DIR}/swig -I${CMAKE_SOURCE_DIR}/include -I${CMAKE_BINARY_DIR}/include -I${CMAKE_BINARY_DIR}/proto)

## Find ZeroMQ
find_package(ZMQ REQUIRED)
include_directories(${ZMQ_INCLUDE_DIRS})

## Find Xerces-C
find_package(XercesC)
set(CALIN_HAVE_XERCESC FALSE)
if(XercesC_FOUND)
	include_directories(${XercesC_INCLUDE_DIRS})
	set(EXTRA_LIBS ${EXTRA_LIBS} ${XercesC_LIBRARIES} )
	set(CALIN_HAVE_XERCESC TRUE)
endif()

## Find the CTA ADH_APIS package
find_package(CTA_ADH_APIS)

## Find Geant4 if requested
SET(CALIN_HAVE_GEANT4 FALSE)
if(NO_GEANT4)
	message(STATUS "Skipping GEANT4 as requested")
	set(GEANT4_FOUND FALSE)
else()
	find_package(Geant4)
	if(GEANT4_FOUND)
		message(STATUS "GEANT4_INCLUDE_DIR: " ${GEANT4_INCLUDE_DIR})
		message(STATUS "GEANT4_LIBRARY_DIR: " ${GEANT4_LIBRARY_DIR})
		include_directories(${GEANT4_INCLUDE_DIR})
		link_directories(${GEANT4_LIBRARY_DIR})
		set(CALIN_HAVE_GEANT4 TRUE)
	endif()
endif()

# ========================================== #
# Add the subdirectories to the build system #
# ========================================== #

add_subdirectory(include)
add_subdirectory(scripts)
add_subdirectory(tools)
add_subdirectory(proto)
add_subdirectory(src)
add_subdirectory(swig)
add_subdirectory(unit_tests)

install(DIRECTORY include/ DESTINATION ${CALIN_HEADER_INSTALL_DIR})
install(DIRECTORY data/ DESTINATION ${CALIN_DATA_INSTALL_DIR})
install(DIRECTORY examples/ DESTINATION ${CALIN_EXAMPLES_INSTALL_DIR})
install(DIRECTORY python/ DESTINATION ${CALIN_PYTHON_INSTALL_DIR})
