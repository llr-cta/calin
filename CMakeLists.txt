# calin/CMakeLists.txt -- Stephen Fegan
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

cmake_minimum_required(VERSION 3.3)
project(calin)
enable_testing()

# ============================================================= #
# Find the various libraries and packages that calin depends on #
# ============================================================= #

cmake_policy(SET CMP0060 NEW)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/CMakeModules/")

include_directories(include ${CMAKE_BINARY_DIR}/include)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall")

### Find Protobuf libraries and set protobuf base path
set(PROTOBUF_IMPORT_DIRS ${PROJECT_SOURCE_DIR}/proto)
find_package(Protobuf REQUIRED)
include_directories(${PROTOBUF_INCLUDE_DIRS} ${CMAKE_BINARY_DIR}/proto)
set(EXTRA_LIBS ${EXTRA_LIBS} ${PROTOBUF_LIBRARY})

## Find FFTW libraries
find_package(FFTW REQUIRED)
include_directories(${FFTW_INCLUDES})
set(EXTRA_LIBS ${EXTRA_LIBS} ${FFTW_LIBRARIES} )

## Find GSL libraries
find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
set(EXTRA_LIBS ${EXTRA_LIBS} ${GSL_LIBRARIES} )

## Find Sqlite3 libraries
find_package(SQLite3 REQUIRED)
include_directories(${SQLITE3_INCLUDE_DIR})
set(EXTRA_LIBS ${EXTRA_LIBS} ${SQLITE3_LIBRARY} )

## Find Python interpreter and libraries
find_package(PythonInterp)
find_package(PythonLibs)
include_directories(${PYTHON_INCLUDE_PATH})

## Find Numpy since we use numpy arrays here
find_package(Numpy)
include_directories(${PYTHON_NUMPY_INCLUDE_DIR})

## Find SWIG and set SWIG include paths to "swig", "include" and "proto"
find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})
set(CMAKE_SWIG_FLAGS -py3 -I${CMAKE_SOURCE_DIR}/swig -I${CMAKE_SOURCE_DIR}/include -I${CMAKE_BINARY_DIR}/proto)

## Find ZeroMQ
find_package(ZMQ REQUIRED)
include_directories(${ZMQ_INCLUDE_DIRS})

## Find the CTA CamerasToACTL package
SET(CALIN_HAVE_CTA_CAMERASTOACTL FALSE)
find_package(CTA_CamerasToACTL)
if(CTA_CAMERASTOACTL_FOUND)
	SET(CALIN_HAVE_CTA_CAMERASTOACTL TRUE)
endif()

## Find Geant4 if requested
SET(CALIN_HAVE_GEANT4 FALSE)
if(NO_GEANT4)
	message(STATUS "Skipping GEANT4 as requested")
	SET(GEANT4_FOUND FALSE)
else()
	find_package(Geant4)
	if(GEANT4_FOUND)
		message(STATUS "GEANT4_INCLUDE_DIR: " ${GEANT4_INCLUDE_DIR})
		message(STATUS "GEANT4_LIBRARY_DIR: " ${GEANT4_LIBRARY_DIR})
		include_directories(${GEANT4_INCLUDE_DIR})
		link_directories(${GEANT4_LIBRARY_DIR})
		SET(CALIN_HAVE_GEANT4 TRUE)
	endif()
endif()

# ================================================ #
# See: https://cmake.org/Wiki/CMake_RPATH_handling #
# ================================================ #

# use, i.e. don't skip the full RPATH for the build tree
#SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
#SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

#SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
#SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

#MESSAGE(${CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES})
# the RPATH to be used when installing, but only if it's not a system directory
#LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
#IF("${isSystemDir}" STREQUAL "-1")
#	SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
#ENDIF("${isSystemDir}" STREQUAL "-1")

# ========================================== #
# Add the subdirectories to the build system #
# ========================================== #

add_subdirectory(include)
add_subdirectory(tools)
add_subdirectory(proto)
add_subdirectory(src)
add_subdirectory(swig)
add_subdirectory(unit_tests)

install(DIRECTORY include/ DESTINATION include/calin)
