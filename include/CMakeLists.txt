# calin/include/CMakeLists.txt -- Stephen Fegan
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

CMAKE_HOST_SYSTEM_INFORMATION(RESULT CALIN_SYSTEM_FQDN QUERY FQDN)
STRING(TIMESTAMP CALIN_BUILD_DATE)

FIND_PACKAGE(Git)
IF(GIT_FOUND)
  EXECUTE_PROCESS(
    COMMAND ${GIT_EXECUTABLE} remote get-url origin
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE CALIN_GIT_ORIGIN_URL
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  EXECUTE_PROCESS(
    COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE CALIN_GIT_BRANCH
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  EXECUTE_PROCESS(
    COMMAND ${GIT_EXECUTABLE} rev-parse -q HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE CALIN_GIT_COMMIT_SHA1
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  EXECUTE_PROCESS(
    COMMAND ${GIT_EXECUTABLE} show -s --format=%cI HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE CALIN_GIT_COMMIT_DATE
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  EXECUTE_PROCESS(
    COMMAND ${GIT_EXECUTABLE} diff-index --quiet HEAD
		WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
    RESULT_VARIABLE CALIN_GIT_IS_DIRTY
  )
  IF(CALIN_GIT_IS_DIRTY)
    SET(CALIN_GIT_REPO_STATUS "dirty")
  ELSE()
    SET(CALIN_GIT_REPO_STATUS "clean")
  ENDIF()
ELSE()
  MESSAGE("Git executable not found : repository provenace will be unavailable")
ENDIF()

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/calin_global_config.hpp.in
  ${CMAKE_CURRENT_BINARY_DIR}/calin_global_config.hpp)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/calin_global_config.hpp
  DESTINATION include/calin)
