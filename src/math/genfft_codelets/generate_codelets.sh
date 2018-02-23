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

CODELET_SIZES="8 12 15 16 18 20 24 28 30 32 36 40 48 56 60 64"
CODELET_CSV=`echo $CODELET_SIZES | tr ' ' ','`

GENFFTDIR=/Users/sfegan/GitHub/fftw3/genfft
INDENT="indent -nbad -bap -nbc -br -brs -c33 -cd33 -ncdb -ce -ci4 -cli0 -d0 -di1 -nfc1 -ip0 -l75 -lp -npcs -npsl -nsc -nut -i2"

cat > define_all_dft_codelets.cpp_incl <<EOF
// Automatically generated file

template<int N> inline void dft_r2c(float_type* r_in, float_type* c_out, int rs = 1, int cs = 1);
template<int N> inline void dft_c2r(float_type* r_in, float_type* c_out, int rs = 1, int cs = 1);
template<int N> inline void dft_r2hc(float_type* r_in, float_type* c_out, int rs = 1, int cs = 1);
template<int N> inline void dft_hc2r(float_type* r_in, float_type* c_out, int rs = 1, int cs = 1);
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
//#include "dft_c2c_${N}.c"

template<> inline void dft_r2c<${N}>(float_type* r_in, float_type* c_out, int rs, int cs)
{
  c_out[1] = ZERO();
  ${COMMENT}c_out[$((2*(N/2)))*cs+1] = ZERO();
  dft_codelet_r2cf_${N}(r_in, r_in+rs, c_out, c_out+1, 2*rs, 2*cs, 2*cs, 1, 0, 0);
}

template<> inline void dft_c2r<${N}>(float_type* r_out, float_type* c_in, int rs, int cs)
{
  dft_codelet_r2cb_${N}(r_out, r_out+rs, c_in, c_in+1, 2*rs, 2*cs, 2*cs, 1, 0, 0);
}

template<> inline void dft_r2hc<${N}>(float_type* r_in, float_type* c_out, int rs, int cs)
{
  dft_codelet_r2cf_${N}(r_in, r_in+rs, c_out, c_out+${N}*cs, 2*rs, cs, -cs, 1, 0, 0);
}

template<> inline void dft_hc2r<${N}>(float_type* r_out, float_type* c_in, int rs, int cs)
{
  dft_codelet_r2cb_${N}(r_out, r_out+rs, c_in, c_in+${N}*cs, 2*rs, cs, -cs, 1, 0, 0);
}

EOF

  # echo Generating r2cf codelet for n=$N
  # $GENFFTDIR/gen_r2cf.native -n ${N} -standalone -fma -generic-arith -compact -name dft_codelet_r2cf_${N} | $INDENT > dft_r2cf_${N}.c
  # echo Generating r2cb codelet for n=$N
  # $GENFFTDIR/gen_r2cb.native -n ${N} -sign 1 -standalone -fma -generic-arith -compact -name dft_codelet_r2cb_${N} | $INDENT > dft_r2cb_${N}.c
  ## echo Generating c2c codelet for n=$N
  ## $GENFFTDIR/gen_notw.native -n ${N} -standalone -fma -generic-arith -compact -name dft_codelet_c2c_${N} | $INDENT > dft_c2c_${N}.c
done

cat >> define_all_dft_codelets.cpp_incl <<EOF

// ******************************** Dispatcher ********************************

template<int N> class Codelet_FixedSizeRealToComplexDFT:
  public FixedSizeRealToComplexDFT<R>
{
public:
  Codelet_FixedSizeRealToComplexDFT(unsigned real_stride = 1, unsigned complex_stride = 1):
    FixedSizeRealToComplexDFT<R>(N, real_stride, complex_stride) { }
  virtual ~Codelet_FixedSizeRealToComplexDFT() { }
  void r2c(R* r_in, R* c_out) override {
    dft_r2c<N>(r_in, c_out, FixedSizeRealToComplexDFT<R>::rs_, FixedSizeRealToComplexDFT<R>::cs_);
  }
  void c2r(R* r_out, R* c_in) override {
    dft_c2r<N>(r_out, c_in, FixedSizeRealToComplexDFT<R>::rs_, FixedSizeRealToComplexDFT<R>::cs_);
  }
};

FixedSizeRealToComplexDFT<R>* new_codelet_r2c_dft(unsigned N,
  unsigned real_stride = 1, unsigned complex_stride = 1)
{
  switch(N) {
EOF
for N in $CODELET_SIZES
do
  echo "  case $N: return new Codelet_FixedSizeRealToComplexDFT<$N>(real_stride, complex_stride);" >> define_all_dft_codelets.cpp_incl
done
cat >> define_all_dft_codelets.cpp_incl <<EOF
  default: return nullptr;
  }
}

std::vector<unsigned> list_available_codelets()
{
  return { $CODELET_CSV };
}
EOF
