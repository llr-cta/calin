#!/bin/bash

# calin/math/fftw_dft_codelets/generate_codelets.sh -- Stephen Fegan
#
# Use the FFTW genfft package to generate DFT codelets of various sizes
#
# Copyright 2018, Stephen Fegan <sfegan@llr.in2p3.fr>
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

CODELET_SIZES="8 12 15 16 18 20 24 30 32 36 40 48 56 60 64"

GENFFTDIR=/Users/sfegan/GitHub/fftw3/genfft
INDENT="indent -nbad -bap -nbc -br -brs -c33 -cd33 -ncdb -ce -ci4 -cli0 -d0 -di1 -nfc1 -ip0 -l75 -lp -npcs -npsl -nsc -nut -i2"

cat > define_all_dft_codelets.cpp_incl <<EOF
// Automatically generated file
EOF

cat > declare_all_dft_codelets.hpp_incl <<EOF
// Automatically generated file

void dft_r2c(int N, float_type* r_in, float_type* c_out, int rs = 1, int cs = 1);
void dft_c2r(int N, float_type* r_in, float_type* c_out, int rs = 1, int cs = 1);
void dft_r2hc(int N, float_type* r_in, float_type* c_out, int rs = 1, int cs = 1);
void dft_hc2r(int N, float_type* r_in, float_type* c_out, int rs = 1, int cs = 1);

template<int N> void dft_r2c(float_type* r_in, float_type* c_out, int rs = 1, int cs = 1);
template<int N> void dft_c2r(float_type* r_in, float_type* c_out, int rs = 1, int cs = 1);
template<int N> void dft_r2hc(float_type* r_in, float_type* c_out, int rs = 1, int cs = 1);
template<int N> void dft_hc2r(float_type* r_in, float_type* c_out, int rs = 1, int cs = 1);
EOF

for N in $CODELET_SIZES
do
  COMMENT=""
  if test "$((N-2*(N/2)))" == "1"
  then
    COMMENT="// "
  fi

  cat >> define_all_dft_codelets.cpp_incl <<EOF

// ********************************* N = ${N} *********************************

#include "dft_r2cf_${N}.c"
#include "dft_r2cb_${N}.c"
#include "dft_c2c_${N}.c"

template<> void dft_r2c<${N}>(float_type* r_in, float_type* c_out, int rs, int cs)
{
  c_out[1] = ZERO();
  ${COMMENT}c_out[$((2*(N/2)))*cs+1] = ZERO();
  dft_codelet_r2cf_${N}(r_in, r_in+rs, c_out, c_out+1, 2*rs, 2*cs, 2*cs, 1, 0, 0);
}

template<> void dft_c2r<${N}>(float_type* r_out, float_type* c_in, int rs, int cs)
{
  dft_codelet_r2cb_${N}(r_out, r_out+rs, c_in, c_in+1, 2*rs, 2*cs, 2*cs, 1, 0, 0);
}

template<> void dft_r2hc<${N}>(float_type* r_in, float_type* c_out, int rs, int cs)
{
  dft_codelet_r2cf_${N}(r_in, r_in+rs, c_out, c_out+${N}*cs, 2*rs, cs, -cs, 1, 0, 0);
}

template<> void dft_hc2r<${N}>(float_type* r_out, float_type* c_in, int rs, int cs)
{
  dft_codelet_r2cb_${N}(r_out, r_out+rs, c_in, c_in+${N}*cs, 2*rs, cs, -cs, 1, 0, 0);
}

EOF

  cat >> declare_all_dft_codelets.hpp_incl <<EOF

// ********************************* N = ${N} *********************************

template<> void dft_r2c<${N}>(float_type* r_in, float_type* c_out, int rs, int cs);
template<> void dft_c2r<${N}>(float_type* r_out, float_type* c_in, int rs, int cs);
template<> void dft_r2hc<${N}>(float_type* r_in, float_type* c_out, int rs, int cs);
template<> void dft_hc2r<${N}>(float_type* r_out, float_type* c_in, int rs, int cs);
EOF

  # echo Generating r2cf codelet for n=$N
  # $GENFFTDIR/gen_r2cf.native -n ${N} -standalone -fma -generic-arith -compact -name dft_codelet_r2cf_${N} | $INDENT > dft_r2cf_${N}.c
  # echo Generating r2cb codelet for n=$N
  # $GENFFTDIR/gen_r2cb.native -n ${N} -standalone -fma -generic-arith -compact -name dft_codelet_r2cb_${N} | $INDENT > dft_r2cb_${N}.c
  # echo Generating c2c codelet for n=$N
  # $GENFFTDIR/gen_notw.native -n ${N} -standalone -fma -generic-arith -compact -name dft_codelet_c2c_${N} | $INDENT > dft_c2c_${N}.c
done

cat >> define_all_dft_codelets.cpp_incl <<EOF

// ******************************** Dispatcher ********************************

void dft_r2c(int N, R* r_in, R* c_out, int rs, int cs)
{
  switch(N) {
EOF
for N in $CODELET_SIZES
do
  echo "  case $N: return dft_r2c<$N>(r_in, c_out, rs, cs);" >> define_all_dft_codelets.cpp_incl
done
cat >> define_all_dft_codelets.cpp_incl <<EOF
  default: throw std::runtime_error("DFT codelet for size " + std::to_string(N) +
    " not available at compile time");
  }
}

void dft_c2r(int N, R* r_in, R* c_out, int rs, int cs)
{
  switch(N) {
EOF
for N in $CODELET_SIZES
do
  echo "  case $N: return dft_c2r<$N>(r_in, c_out, rs, cs);" >> define_all_dft_codelets.cpp_incl
done
cat >> define_all_dft_codelets.cpp_incl <<EOF
  default: throw std::runtime_error("DFT codelet for size " +
    std::to_string(N) + " not available at compile time");
  }
}

void dft_r2hc(int N, R* r_in, R* c_out, int rs, int cs)
{
  switch(N) {
EOF
for N in $CODELET_SIZES
do
  echo "  case $N: return dft_r2hc<$N>(r_in, c_out, rs, cs);" >> define_all_dft_codelets.cpp_incl
done
cat >> define_all_dft_codelets.cpp_incl <<EOF
  default: throw std::runtime_error("DFT codelet for size " +
    std::to_string(N) + " not available at compile time");
  }
}

void dft_hc2r(int N, R* r_in, R* c_out, int rs, int cs)
{
  switch(N) {
EOF
for N in $CODELET_SIZES
do
  echo "  case $N: return dft_hc2r<$N>(r_in, c_out, rs, cs);" >> define_all_dft_codelets.cpp_incl
done
cat >> define_all_dft_codelets.cpp_incl <<EOF
  default: throw std::runtime_error("DFT codelet for size " +
    std::to_string(N) + " not available at compile time");
  }
}
EOF

cp declare_all_dft_codelets.hpp_incl ../../../include/math/genfft_codelets
