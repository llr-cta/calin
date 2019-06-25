# calin/src/cminpack/CMakeLists.txt -- Stephen Fegan
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

set (cminpack_srcs
#  cminpack.h
  cminpackP.h
  chkder.c  enorm.c   hybrd1.c  hybrj.c   lmdif1.c  lmstr1.c  qrfac.c   r1updt.c
  dogleg.c  fdjac1.c  hybrd.c   lmder1.c  lmdif.c   lmstr.c   qrsolv.c  rwupdt.c
  dpmpar.c  fdjac2.c  hybrj1.c  lmder.c   lmpar.c   qform.c   r1mpyq.c  covar.c covar1.c
#  minpack.h
#  chkder_.c enorm_.c  hybrd1_.c hybrj_.c  lmdif1_.c lmstr1_.c qrfac_.c  r1updt_.c
# dogleg_.c fdjac1_.c hybrd_.c  lmder1_.c lmdif_.c  lmstr_.c  qrsolv_.c rwupdt_.c
# dpmpar_.c fdjac2_.c hybrj1_.c lmder_.c  lmpar_.c  qform_.c  r1mpyq_.c covar_.c
  )
set (cminpack_hdrs
    cminpack.h minpack.h)

include_directories(${CMAKE_SOURCE_DIR}/include/cminpack)
add_library (calin_cminpack SHARED ${cminpack_srcs})

install(TARGETS calin_cminpack DESTINATION ${CALIN_LIB_INSTALL_DIR})

#install (TARGETS cminpack
 #  LIBRARY DESTINATION ${CMINPACK_LIB_INSTALL_DIR} COMPONENT library
  # ARCHIVE DESTINATION ${CMINPACK_LIB_INSTALL_DIR} COMPONENT library
  # RUNTIME DESTINATION ${CMINPACK_LIB_INSTALL_DIR} COMPONENT library)
#install (FILES ${cminpack_hdrs} DESTINATION ${CMINPACK_INCLUDE_INSTALL_DIR}
 #   COMPONENT cminpack_hdrs)
