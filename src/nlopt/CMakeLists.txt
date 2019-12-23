# calin/src/nlopt/CMakeLists.txt -- Stephen Fegan
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

include_directories(${CMAKE_SOURCE_DIR}/include/nlopt ${CMAKE_CURRENT_BINARY_DIR})

INCLUDE (CheckFunctionExists)
INCLUDE (CheckIncludeFiles)
INCLUDE (CheckLibraryExists)
INCLUDE (CheckSymbolExists)
INCLUDE (CheckTypeSize)
INCLUDE (CheckCSourceCompiles)

SET(NLOPT_MAJOR_VERSION 2)
SET(NLOPT_MINOR_VERSION 4)
SET(NLOPT_BUGFIX_VERSION 2)

CHECK_FUNCTION_EXISTS (copysign HAVE_COPYSIGN)
CHECK_INCLUDE_FILES (dlfcn.h HAVE_DLFCN_H)
CHECK_INCLUDE_FILES (getopt.h HAVE_GETOPT_H)
CHECK_FUNCTION_EXISTS (getpid HAVE_GETPID)
CHECK_C_SOURCE_COMPILES("#include<sys/syscall.h>
#include<unistd.h>
int main() { syscall(SYS_gettid); return 0; }" HAVE_GETTID_SYSCALL)
CHECK_FUNCTION_EXISTS (gettimeofday HAVE_GETTIMEOFDAY)
CHECK_INCLUDE_FILES (inttypes.h HAVE_INTTYPES_H)
CHECK_FUNCTION_EXISTS (isinf HAVE_ISINF)
CHECK_FUNCTION_EXISTS (isnan HAVE_ISNAN)
CHECK_LIBRARY_EXISTS(m sin "" HAVE_LIBM)
CHECK_INCLUDE_FILES (memory.h HAVE_MEMORY_H)
CHECK_FUNCTION_EXISTS (qsort_r HAVE_QSORT_R)
CHECK_INCLUDE_FILES (stdint.h HAVE_STDINT_H)
CHECK_INCLUDE_FILES (stdlib.h HAVE_STDLIB_H)
CHECK_INCLUDE_FILES (strings.h HAVE_STRINGS_H)
CHECK_INCLUDE_FILES (string.h HAVE_STRING_H)
CHECK_INCLUDE_FILES (sys/stat.h HAVE_SYS_STAT_H)
CHECK_INCLUDE_FILES (sys/types.h HAVE_SYS_TYPES_H)
CHECK_FUNCTION_EXISTS (time HAVE_TIME)
CHECK_TYPE_SIZE(uint32_t HAVE_UINT32_T)
CHECK_TYPE_SIZE("unsigned int" SIZEOF_UNSIGNED_INT)
CHECK_TYPE_SIZE("unsigned long" SIZEOF_UNSIGNED_LONG)
CHECK_C_SOURCE_COMPILES("__thread int i; int main() { return 0; }" HAVE_THREADLOCAL___THREAD)
CHECK_C_SOURCE_COMPILES("#include<time.h>
#include<sys/time.h>
int main() { return 0; }" TIME_WITH_SYS_TIME)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/config.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/config.h)

# Relax some of the -Wall options here
add_compile_options(-Wno-unknown-warning-option)
add_compile_options(-Wno-deprecated-declarations)
add_compile_options(-Wno-self-assign)
add_compile_options(-Wno-unused-const-variable)
add_compile_options(-Wno-maybe-uninitialized)
add_compile_options(-Wno-unused-but-set-variable)
add_compile_options(-Wno-array-bounds)

add_library(calin_nlopt SHARED DIRect.c DIRserial.c DIRsubrout.c auglag.c bobyqa.c ccsa_quadratic.c cdirect.c cobyla.c crs.c deprecated.c direct_wrap.c esch.c f77api.c general.c global.cc hybrid.c isres.c linalg.cc local.cc mlsl.c mma.c mssubs.c mt19937ar.c newuoa.c nldrmd.c optimize.c options.c plip.c plis.c pnet.c praxis.c pssubs.c qsort_r.c redblack.c redblack_test.c rescale.c sbplx.c slsqp.c sobolseq.c stogo.cc stop.c timer.c tools.cc)

install(TARGETS calin_nlopt DESTINATION ${CALIN_LIB_INSTALL_DIR})
