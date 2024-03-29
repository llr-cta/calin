# calin/src/f2c/CMakeLists.txt -- Stephen Fegan
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

include_directories(${CMAKE_SOURCE_DIR}/include/f2c)

add_compile_options(-Wno-unknown-warning-option
  -Wno-shift-op-parentheses -Wno-bitwise-op-parentheses
  -Wno-maybe-uninitialized -Wno-parentheses -Wno-format-extra-args)

add_library(calin_f2c SHARED f77vers.c i77vers.c s_rnge.c abort_.c exit_.c
        getenv_.c signal_.c s_stop.c s_paus.c
        system_.c cabs.c ctype.c derf_.c derfc_.c erf_.c erfc_.c
        sig_die.c uninit.c pow_ci.c pow_dd.c pow_di.c pow_hh.c
        pow_ii.c pow_ri.c pow_zi.c pow_zz.c c_abs.c c_cos.c c_div.c
        c_exp.c c_log.c c_sin.c c_sqrt.c z_abs.c z_cos.c z_div.c
        z_exp.c z_log.c z_sin.c z_sqrt.c r_abs.c r_acos.c r_asin.c
        r_atan.c r_atn2.c r_cnjg.c r_cos.c r_cosh.c r_dim.c r_exp.c
        r_imag.c r_int.c r_lg10.c r_log.c r_mod.c r_nint.c r_sign.c
        r_sin.c r_sinh.c r_sqrt.c r_tan.c r_tanh.c d_abs.c d_acos.c
        d_asin.c d_atan.c d_atn2.c d_cnjg.c d_cos.c d_cosh.c d_dim.c
        d_exp.c d_imag.c d_int.c d_lg10.c d_log.c d_mod.c d_nint.c
        d_prod.c d_sign.c d_sin.c d_sinh.c d_sqrt.c d_tan.c d_tanh.c
        i_abs.c i_dim.c i_dnnt.c i_indx.c i_len.c i_mod.c i_nint.c
        i_sign.c lbitbits.c lbitshft.c h_abs.c h_dim.c h_dnnt.c
        h_indx.c h_len.c h_mod.c h_nint.c h_sign.c l_ge.c l_gt.c
        l_le.c l_lt.c hl_ge.c hl_gt.c hl_le.c hl_lt.c ef1asc_.c
        ef1cmc_.c f77_aloc.c s_cat.c s_cmp.c s_copy.c backspac.c
        close.c dfe.c dolio.c due.c endfile.c err.c fmt.c fmtlib.c
        ftell_.c iio.c ilnw.c inquire.c lread.c lwrite.c open.c
        rdfmt.c rewind.c rsfe.c rsli.c rsne.c sfe.c sue.c typesize.c
        uio.c util.c wref.c wrtfmt.c wsfe.c wsle.c wsne.c xwsne.c
        dtime_.c etime_.c)

#// main.c getarg_.c iargc_.c pow_qq.c qbitbits.c qbitshft.c ftell64_.c

install(TARGETS calin_f2c DESTINATION ${CALIN_LIB_INSTALL_DIR})
