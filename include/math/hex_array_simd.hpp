/*

   calin/math/hex_array_simd.hpp -- Stephen Fegan -- 2018-01-05

   Collection of functions which translate between hexagonal and Cartesian
   geometries, and provide other useful calculations for hex grids.

   AVX SIMD version.

   Copyright 2018, Stephen Fegan <sfegan@llr.in2p3.fr>
   Laboratoire Leprince-Ringuet, CNRS/IN2P3, Ecole Polytechnique, Institut Polytechnique de Paris

   This file is part of "calin"

   "calin" is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 or
   later, as published by the Free Software Foundation.

   "calin" is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

*/

#pragma once

#include <math/simd.hpp>
#include <math/hex_array.hpp>

namespace calin { namespace math { namespace hex_array {

#if defined(__AVX2__) and defined(__FMA__)

// *****************************************************************************
//
// General functions for the hex array building to UV <-> HexID and neighbors
//
// *****************************************************************************

CALIN_MM256PS_CONST(_c_m256_one_half, 0.5f);
#define _c_m256_vx _c_m256_one_half
CALIN_MM256PS_CONST(_c_m256_vy, 0.5f*CALIN_HEX_ARRAY_SQRT3);
CALIN_MM256PS_CONST(_c_m256_vy_inv, 1.0f/(0.5f*CALIN_HEX_ARRAY_SQRT3));
CALIN_MM256PS_CONST(_c_m256_one_third, 1.0f/3.0f);
CALIN_MM256PS_CONST(_c_m256_four_thirds, 4.0f/3.0f);

inline __m256i avx2_positive_hexid_to_ringid_loop(const __m256i hexid)
{
  // This algorithm is relatively slow in comparisson to the scalar version
  // but still faster overall conidering we compute 8 rigids in one go
  const __m256i six = _mm256_set1_epi32(6);
  const __m256i one = _mm256_set1_epi32(1);
  __m256i ringid = _mm256_setzero_si256();
  __m256i nsites = one;
  __m256i nring = _mm256_setzero_si256();
  __m256i mask = _mm256_cmpgt_epi32(nsites, hexid);
  while(~_mm256_movemask_epi8(mask)) {
    ringid = _mm256_blendv_epi8(_mm256_add_epi32(ringid, one), ringid, mask);
    nring = _mm256_add_epi32(nring, six);
    nsites = _mm256_add_epi32(nsites, nring);
    mask = _mm256_cmpgt_epi32(nsites, hexid);
  }
  return ringid;
}

inline __m256i avx2_positive_hexid_to_ringid_root(const __m256i hexid)
{
  // The following algorithm works until hexid=12,589,056
  // const unsigned iarg = 1+4*(hexid-1)/3;
  // return (unsigned(std::sqrt(float(iarg)))+1)/2;
  __m256 arg = _mm256_cvtepi32_ps(hexid);
  arg = _mm256_fmsub_ps(arg, calin::math::simd::c_m256(_c_m256_four_thirds),
    calin::math::simd::c_m256(_c_m256_one_third));
  arg = _mm256_sqrt_ps(arg);
  arg = _mm256_fmadd_ps(arg, calin::math::simd::c_m256(_c_m256_one_half),
    calin::math::simd::c_m256(_c_m256_one_half));
  arg = _mm256_floor_ps(arg);
  return _mm256_cvtps_epi32(arg);
}

inline __m256i avx2_positive_hexid_to_ringid(const __m256i hexid)
{
  return avx2_positive_hexid_to_ringid_root(hexid);
}

inline __m256i avx2_hexid_to_ringid(const __m256i hexid)
{
  const __m256i mask = _mm256_cmpeq_epi32(hexid, _mm256_setzero_si256());
  return _mm256_andnot_si256(mask, avx2_positive_hexid_to_ringid(hexid));
}

inline __m256i avx2_ringid_to_nsites_contained(const __m256i ringid)
{
  // return 3*ringid*(ringid+1)+1;
  const __m256i one = _mm256_set1_epi32(1);
  __m256i nsites = _mm256_add_epi32(ringid, one);
  nsites = _mm256_mullo_epi32(ringid, nsites);
  nsites = _mm256_sub_epi32(_mm256_slli_epi32(nsites, 2), nsites);
  nsites = _mm256_add_epi32(nsites, one);
  return nsites;
}

inline void avx2_positive_hexid_to_ringid_segid_runid(
  const __m256i hexid, __m256i& ringid, __m256i& segid, __m256i& runid)
{
  // ringid = positive_hexid_to_ringid(hexid);
  // unsigned iring = hexid - ringid_to_nsites_contained(ringid-1);
  // segid = int(iring/ringid);
  // runid = iring - segid*ringid;
  const __m256i one = _mm256_set1_epi32(1);
  ringid = avx2_positive_hexid_to_ringid(hexid);
  runid = _mm256_sub_epi32(hexid,
    avx2_ringid_to_nsites_contained(_mm256_sub_epi32(ringid,one)));
  segid = _mm256_setzero_si256();

  const __m256i ringid_minus_one = _mm256_sub_epi32(ringid, one);

  __m256i mask = _mm256_cmpgt_epi32(runid, ringid_minus_one);
  runid = _mm256_sub_epi32(runid, _mm256_and_si256(mask, ringid));
  segid = _mm256_add_epi32(segid, _mm256_and_si256(mask, one));

  mask = _mm256_cmpgt_epi32(runid, ringid_minus_one);
  runid = _mm256_sub_epi32(runid, _mm256_and_si256(mask, ringid));
  segid = _mm256_add_epi32(segid, _mm256_and_si256(mask, one));

  mask = _mm256_cmpgt_epi32(runid, ringid_minus_one);
  runid = _mm256_sub_epi32(runid, _mm256_and_si256(mask, ringid));
  segid = _mm256_add_epi32(segid, _mm256_and_si256(mask, one));

  mask = _mm256_cmpgt_epi32(runid, ringid_minus_one);
  runid = _mm256_sub_epi32(runid, _mm256_and_si256(mask, ringid));
  segid = _mm256_add_epi32(segid, _mm256_and_si256(mask, one));

  mask = _mm256_cmpgt_epi32(runid, ringid_minus_one);
  runid = _mm256_sub_epi32(runid, _mm256_and_si256(mask, ringid));
  segid = _mm256_add_epi32(segid, _mm256_and_si256(mask, one));
}

inline void avx2_hexid_to_ringid_segid_runid(
  const __m256i hexid, __m256i& ringid, __m256i& segid, __m256i& runid)
{
  // if(hexid==0) { ringid = segid = runid = 0; return; }
  // return positive_hexid_to_ringid_segid_runid(hexid, ringid, segid, runid);
  avx2_positive_hexid_to_ringid_segid_runid(hexid, ringid, segid, runid);
  const __m256i mask = _mm256_cmpeq_epi32(hexid, _mm256_setzero_si256());
  ringid = _mm256_andnot_si256(mask, ringid);
  segid = _mm256_andnot_si256(mask, segid);
  runid = _mm256_andnot_si256(mask, runid);
}

inline __m256i avx2_positive_ringid_segid_runid_to_hexid(
  const __m256i ringid, const __m256i segid, const __m256i runid)
{
  // return ringid_to_nsites_contained(ringid-1)+segid*ringid+runid;
  const __m256i one = _mm256_set1_epi32(1);
  __m256i nsites = avx2_ringid_to_nsites_contained(_mm256_sub_epi32(ringid, one));
  nsites = _mm256_add_epi32(nsites, _mm256_mullo_epi32(segid, ringid));
  nsites = _mm256_add_epi32(nsites, runid);
  return nsites;
}

inline __m256i avx2_ringid_segid_runid_to_hexid(
  const __m256i ringid, const __m256i segid, const __m256i runid)
{
  // return (ringid==0) ? 0 :
  //     positive_ringid_segid_runid_to_hexid(ringid, segid, runid);
  const __m256i mask = _mm256_cmpeq_epi32(ringid, _mm256_setzero_si256());
  return _mm256_andnot_si256(mask,
    avx2_positive_ringid_segid_runid_to_hexid(ringid, segid, runid));
}

inline __m256i avx2_uv_to_ringid(const __m256i u, const __m256i v)
{
  // return static_cast<unsigned>(std::max({std::abs(u), std::abs(v),
  //           std::abs(u+v)}));
  __m256i ringid = _mm256_abs_epi32(u);
  ringid = _mm256_max_epu32(ringid, _mm256_abs_epi32(v));
  ringid = _mm256_max_epu32(ringid, _mm256_abs_epi32(_mm256_add_epi32(u,v)));
  return ringid;
}

inline void avx2_hexid_to_uv_ccw(const __m256i hexid, __m256i& u, __m256i& v)
{
  // if(hexid==0) { u = v = 0; return; }
  // unsigned ringid;
  // unsigned segid;
  // unsigned runid;
  // positive_hexid_to_ringid_segid_runid(hexid, ringid, segid, runid);
  // switch(segid)
  // {
  //   case 0: u = ringid-runid; v = runid;        break;
  //   case 1: u = -runid;       v = ringid;       break;
  //   case 2: u = -ringid;      v = ringid-runid; break;
  //   case 3: u = runid-ringid; v = -runid;       break;
  //   case 4: u = runid;        v = -ringid;      break;
  //   case 5: u = ringid;       v = runid-ringid; break;
  //   default: assert(0);
  // }
  const __m256i one = _mm256_set1_epi32(1);
  __m256i ringid = avx2_positive_hexid_to_ringid(hexid);
  __m256i iring = _mm256_sub_epi32(hexid,
    avx2_ringid_to_nsites_contained(_mm256_sub_epi32(ringid,one)));

  u = ringid;
  v = _mm256_setzero_si256();

  __m256i irun = _mm256_min_epu32(iring, ringid);
  u = _mm256_sub_epi32(u, irun);
  v = _mm256_add_epi32(v, irun);
  iring = _mm256_sub_epi32(iring, irun);

  irun = _mm256_min_epu32(iring, ringid);
  u = _mm256_sub_epi32(u, irun);
  iring = _mm256_sub_epi32(iring, irun);

  irun = _mm256_min_epu32(iring, ringid);
  v = _mm256_sub_epi32(v, irun);
  iring = _mm256_sub_epi32(iring, irun);

  irun = _mm256_min_epu32(iring, ringid);
  u = _mm256_add_epi32(u, irun);
  v = _mm256_sub_epi32(v, irun);
  iring = _mm256_sub_epi32(iring, irun);

  irun = _mm256_min_epu32(iring, ringid);
  u = _mm256_add_epi32(u, irun);
  iring = _mm256_sub_epi32(iring, irun);

  v = _mm256_add_epi32(v, iring);

  const __m256i mask = _mm256_cmpeq_epi32(hexid, _mm256_setzero_si256());
  u = _mm256_andnot_si256(mask, u);
  v = _mm256_andnot_si256(mask, v);
}

inline void avx2_hexid_to_uv_cw(const __m256i hexid, __m256i& u, __m256i& v)
{
#if 0 // This code is correct but it's not worth maintaining two versions
  const __m256i one = _mm256_set1_epi32(1);
  __m256i ringid = avx2_positive_hexid_to_ringid(hexid);
  __m256i iring = _mm256_sub_epi32(hexid,
    avx2_ringid_to_nsites_contained(_mm256_sub_epi32(ringid,one)));

  u = ringid;
  v = _mm256_setzero_si256();

  __m256i irun = _mm256_min_epu32(iring, ringid);
  v = _mm256_sub_epi32(v, irun);
  iring = _mm256_sub_epi32(iring, irun);

  irun = _mm256_min_epu32(iring, ringid);
  u = _mm256_sub_epi32(u, irun);
  iring = _mm256_sub_epi32(iring, irun);

  irun = _mm256_min_epu32(iring, ringid);
  u = _mm256_sub_epi32(u, irun);
  v = _mm256_add_epi32(v, irun);
  iring = _mm256_sub_epi32(iring, irun);

  irun = _mm256_min_epu32(iring, ringid);
  v = _mm256_add_epi32(v, irun);
  iring = _mm256_sub_epi32(iring, irun);

  irun = _mm256_min_epu32(iring, ringid);
  u = _mm256_add_epi32(u, irun);
  iring = _mm256_sub_epi32(iring, irun);

  u = _mm256_add_epi32(u, irun);
  v = _mm256_add_epi32(v, iring);

  const __m256i mask = _mm256_cmpeq_epi32(hexid, _mm256_setzero_si256());
  u = _mm256_andnot_si256(mask, u);
  v = _mm256_andnot_si256(mask, v);
#else
  // hexid_to_uv_ccw(hexid, u, v);
  // u += v;
  // v = -v;
  avx2_hexid_to_uv_ccw(hexid, u, v);
  u = _mm256_add_epi32(u, v);
  v = _mm256_sign_epi32(v, _mm256_cmpeq_epi32(v, v));
#endif
}

inline __m256i avx2_uv_to_hexid_ccw(const __m256i u, const __m256i v)
{
  // if(u==0 and v==0)return 0;
  // int ringid = uv_to_ringid(u,v);
  // unsigned segid;
  // int runid;
  // int upv = u+v;
  // if(upv==ringid and v!=ringid)         { segid=0; runid=v; }
  // else if(v==ringid and u!=-ringid)     { segid=1; runid=-u; }
  // else if(u==-ringid and upv!=-ringid)  { segid=2; runid=ringid-v; }
  // else if(u+v==-ringid and v!=-ringid)  { segid=3; runid=-v; }
  // else if(v==-ringid and u!=ringid)     { segid=4; runid=u; }
  // else /*if(u==ringid and upv!=ringid)*/{ segid=5; runid=ringid+v; }
  // return positive_ringid_segid_runid_to_hexid(ringid, segid, runid);
  const __m256i one = _mm256_set1_epi32(1);
  const __m256i minus_one = _mm256_set1_epi32(-1);
  const __m256i ringid = avx2_uv_to_ringid(u,v);
  const __m256i minus_ringid = _mm256_sign_epi32(ringid, minus_one);
  const __m256i upv = _mm256_add_epi32(u, v);
  __m256i not_found_mask = minus_one;
  __m256i hexid = avx2_ringid_to_nsites_contained(_mm256_sub_epi32(ringid, one));

  // Seg ID = 0
  // if(upv==ringid and v!=ringid)         { segid=0; runid=v; }
  __m256i here_mask = _mm256_cmpeq_epi32(upv, ringid);
  hexid = _mm256_add_epi32(hexid, _mm256_and_si256(not_found_mask,
    _mm256_blendv_epi8(ringid, v, here_mask)));
  not_found_mask = _mm256_andnot_si256(here_mask, not_found_mask);
  // hexid = _mm256_add_epi32(hexid, _mm256_or_si256(
  //   _mm256_and_si256(here_mask, v),
  //   _mm256_and_si256(not_found_mask, ringid)));

  // Seg ID = 1
  // else if(v==ringid and u!=-ringid)     { segid=1; runid=-u; }
  here_mask = _mm256_cmpeq_epi32(v, ringid);
  hexid = _mm256_sub_epi32(hexid, _mm256_and_si256(not_found_mask,
    _mm256_blendv_epi8(minus_ringid, u, here_mask)));
  not_found_mask = _mm256_andnot_si256(here_mask, not_found_mask);
  // hexid = _mm256_sub_epi32(hexid, _mm256_or_si256(
  //   _mm256_and_si256(here_mask, u),
  //   _mm256_and_si256(not_found_mask, minus_ringid)));

  // Seg ID = 2
  // else if(u==-ringid and upv!=-ringid)  { segid=2; runid=ringid-v; }
  here_mask = _mm256_cmpeq_epi32(u, minus_ringid);
  hexid = _mm256_sub_epi32(hexid, _mm256_and_si256(not_found_mask,
    _mm256_blendv_epi8(minus_ringid, upv, here_mask)));
  not_found_mask = _mm256_andnot_si256(here_mask, not_found_mask);
  // hexid = _mm256_sub_epi32(hexid, _mm256_or_si256(
  //   _mm256_and_si256(here_mask, upv),
  //   _mm256_and_si256(not_found_mask, minus_ringid)));

  // Seg ID = 3
  // else if(u+v==-ringid and v!=-ringid)  { segid=3; runid=-v; }
  here_mask = _mm256_cmpeq_epi32(upv, minus_ringid);
  hexid = _mm256_sub_epi32(hexid, _mm256_and_si256(not_found_mask,
    _mm256_blendv_epi8(minus_ringid, v, here_mask)));
  not_found_mask = _mm256_andnot_si256(here_mask, not_found_mask);
  // hexid = _mm256_sub_epi32(hexid, _mm256_or_si256(
  //   _mm256_and_si256(here_mask, v),
  //   _mm256_and_si256(not_found_mask, minus_ringid)));

  // Seg ID = 4
  // else if(v==-ringid and u!=ringid)     { segid=4; runid=u; }
  here_mask = _mm256_cmpeq_epi32(v, minus_ringid);
  hexid = _mm256_add_epi32(hexid, _mm256_and_si256(not_found_mask,
    _mm256_blendv_epi8(ringid, u, here_mask)));
  not_found_mask = _mm256_andnot_si256(here_mask, not_found_mask);
  // hexid = _mm256_add_epi32(hexid, _mm256_or_si256(
  //   _mm256_and_si256(here_mask, u),
  //   _mm256_and_si256(not_found_mask, ringid)));

  // Seg ID = 5
  // else /*if(u==ringid and upv!=ringid)*/{ segid=5; runid=ringid+v; }
  hexid = _mm256_add_epi32(hexid, _mm256_and_si256(not_found_mask, upv));

  const __m256i mask = _mm256_cmpeq_epi32(ringid, _mm256_setzero_si256());
  hexid = _mm256_andnot_si256(mask, hexid);
  return hexid;
}

inline __m256i avx2_uv_to_hexid_cw(__m256i u, __m256i v)
{
  // u += v;
  // v = -v;
  // return uv_to_hexid_ccw(u, v);
  u = _mm256_add_epi32(u, v);
  v = _mm256_sign_epi32(v, _mm256_cmpeq_epi32(v, v));
  return avx2_uv_to_hexid_ccw(u, v);
}

inline void avx2_hexid_to_uv(__m256i hexid, __m256i& u, __m256i& v)
{
  // hexid_to_uv_ccw(hexid, u, v);
  avx2_hexid_to_uv_ccw(hexid, u, v);
}

inline __m256i avx2_uv_to_hexid(__m256i u, __m256i v)
{
  // return uv_to_hexid_ccw(u,v);
  return avx2_uv_to_hexid_ccw(u, v);
}

inline void avx2_hexid_to_uv(__m256i hexid, __m256i& u, __m256i& v, bool clockwise)
{
  if(clockwise)avx2_hexid_to_uv_cw(hexid, u, v);
  else avx2_hexid_to_uv_ccw(hexid, u, v);
}

inline __m256i avx2_uv_to_hexid(__m256i u, __m256i v, bool clockwise)
{
  if(clockwise)return avx2_uv_to_hexid_cw(u, v);
  else return avx2_uv_to_hexid_ccw(u, v);
}

// *****************************************************************************
//
// XY <-> UV
//
// *****************************************************************************

// XY <-> UV without rotation matrix

inline void avx2_uv_to_xy_f(const __m256i u, const __m256i v, __m256& x, __m256& y)
{
  y = _mm256_cvtepi32_ps(v);
  x = _mm256_cvtepi32_ps(u);
  x = _mm256_fmadd_ps(y, calin::math::simd::c_m256(_c_m256_vx), x);
  y = _mm256_mul_ps(y, calin::math::simd::c_m256(_c_m256_vy));
}

inline void avx2_xy_to_uv_f(__m256 x, __m256 y, __m256i& u, __m256i& v)
{
  // Convert X,Y first into U,V space then round to nearest
  // integer. That gets us close to correct answer, mapping XY to a
  // lozenge-shaped space rather than hexagonal. We then correct the
  // four regions that lie outside the hexagonal cell assigning them
  // to their correct neighboring cell.
  // Writer's note: see ~/Google Drive/Work/calin

  // double dv = y*c_vy_inv;
  // double du = x-dv*c_vx;
  // u = std::lround(du);
  // v = std::lround(dv);
  // du -= u;
  // dv -= v;

  y = _mm256_mul_ps(y, calin::math::simd::c_m256(_c_m256_vy_inv));
  x = _mm256_fnmadd_ps(y, calin::math::simd::c_m256(_c_m256_vx), x);
  u = _mm256_cvtps_epi32(x);
  v = _mm256_cvtps_epi32(y);
  x = _mm256_sub_ps(x, _mm256_cvtepi32_ps(u));
  y = _mm256_sub_ps(y, _mm256_cvtepi32_ps(v));

  // double c3 = dv-du;
  const __m256i c3 = _mm256_castps_si256(_mm256_sub_ps(y, x));

  __m256i uvshift;
  __m256i mask;

  // double c1 = du+0.5*dv;
  // double c2 = dv+0.5*du;
  // if(c3<0) {
  //   if(c1>=1) u++;
  //   else if(c2<-1) v--;
  // } else {
  //   if(c2>=1) v++;
  //   else if(c1<-1) u--;
  // }

  uvshift = _mm256_cvtps_epi32(_mm256_fmadd_ps(y, calin::math::simd::c_m256(_c_m256_one_half), x));
  mask = _mm256_srai_epi32(_mm256_xor_si256(uvshift, c3), 31);
  u = _mm256_blendv_epi8(u, _mm256_add_epi32(u, uvshift), mask);

  uvshift = _mm256_cvtps_epi32(_mm256_fmadd_ps(x, calin::math::simd::c_m256(_c_m256_one_half), y));
  mask = _mm256_srai_epi32(_mm256_xor_si256(uvshift, c3), 31);
  v = _mm256_blendv_epi8(_mm256_add_epi32(v, uvshift), v, mask);
}

inline void avx2_xy_to_uv_with_remainder_f(
  __m256& x_in_dx_out, __m256& y_in_dy_out, __m256i& u, __m256i& v)
{
  avx2_xy_to_uv_f(x_in_dx_out, y_in_dy_out, u, v);
  x_in_dx_out = _mm256_sub_ps(x_in_dx_out, _mm256_cvtepi32_ps(u));
  __m256 vf = _mm256_cvtepi32_ps(v);
  x_in_dx_out = _mm256_fnmadd_ps(vf, calin::math::simd::c_m256(_c_m256_vx), x_in_dx_out);
  y_in_dy_out = _mm256_fnmadd_ps(vf, calin::math::simd::c_m256(_c_m256_vy), y_in_dy_out);
}

// XY <-> UV with (shear-free) affine transformation

inline void avx2_uv_to_xy_trans_f(const __m256i u, const __m256i v, __m256& x, __m256& y,
  const float crot, const float srot, const float scale, const float dx = 0, const float dy = 0)
{
  // uv_to_xy(u,v,x,y);
  // double xx = x*crot - y*srot;
  // y = scale * (y*crot + x*srot) + dy;
  // x = scale * xx + dx;
  avx2_uv_to_xy_f(u,v,x,y);
  const __m256 vsrot = _mm256_set1_ps(srot);
  const __m256 vcrot = _mm256_set1_ps(crot);
  const __m256 vscale = _mm256_set1_ps(scale);
  __m256 xx = _mm256_mul_ps(y, vsrot);
  xx = _mm256_fmsub_ps(x, vcrot, xx);
  y = _mm256_mul_ps(y, vcrot);
  y = _mm256_fmadd_ps(x, vsrot, y);
  x = xx;
  x = _mm256_fmadd_ps(x, vscale, _mm256_set1_ps(dx));
  y = _mm256_fmadd_ps(y, vscale, _mm256_set1_ps(dy));
}

inline void avx2_xy_trans_to_uv_f(__m256 x, __m256 y, __m256i& u, __m256i& v,
  const float crot, const float srot, const float scale, const float dx = 0, const float dy = 0)
{
  // x = (x - dx)/scale;
  // y = (y - dy)/scale;
  // double xx = x*crot + y*srot;
  // y = y*crot - x*srot;
  // xy_to_uv(xx,y,u,v);
  const __m256 vsrot = _mm256_set1_ps(srot);
  const __m256 vcrot = _mm256_set1_ps(crot);
  const __m256 vscale = _mm256_set1_ps(1.0f/scale);
  x = _mm256_mul_ps(_mm256_sub_ps(x, _mm256_set1_ps(dx)), vscale);
  y = _mm256_mul_ps(_mm256_sub_ps(y, _mm256_set1_ps(dy)), vscale);
  __m256 yy = _mm256_mul_ps(x, vsrot);
  yy = _mm256_fmsub_ps(y, vcrot, yy);
  x = _mm256_mul_ps(x, vcrot);
  x = _mm256_fmadd_ps(y, vsrot, x);
  avx2_xy_to_uv_f(x, yy, u, v);
}

inline void avx2_xy_trans_to_uv_with_remainder_f(
  __m256& x_in_dx_out, __m256& y_in_dy_out, __m256i& u, __m256i& v,
  const float crot, const float srot, const float scale, const float dx = 0, const float dy = 0)
{
  const __m256 vsrot = _mm256_set1_ps(srot);
  const __m256 vcrot = _mm256_set1_ps(crot);
  __m256 vscale = _mm256_set1_ps(1.0f/scale);
  x_in_dx_out = _mm256_mul_ps(_mm256_sub_ps(x_in_dx_out, _mm256_set1_ps(dx)), vscale);
  y_in_dy_out = _mm256_mul_ps(_mm256_sub_ps(y_in_dy_out, _mm256_set1_ps(dy)), vscale);
  __m256 yy = _mm256_mul_ps(x_in_dx_out, vsrot);
  yy = _mm256_fmsub_ps(y_in_dy_out, vcrot, yy);
  x_in_dx_out = _mm256_mul_ps(x_in_dx_out, vcrot);
  x_in_dx_out = _mm256_fmadd_ps(y_in_dy_out, vsrot, x_in_dx_out);
  avx2_xy_to_uv_with_remainder_f(x_in_dx_out, yy, u, v);

  vscale = _mm256_set1_ps(scale);
  y_in_dy_out = _mm256_mul_ps(yy, vcrot);
  y_in_dy_out = _mm256_fmadd_ps(x_in_dx_out, vsrot, y_in_dy_out);
  y_in_dy_out = _mm256_mul_ps(y_in_dy_out, vscale);

  x_in_dx_out = _mm256_mul_ps(x_in_dx_out, vcrot);
  x_in_dx_out = _mm256_fnmadd_ps(yy, vsrot, x_in_dx_out);
  x_in_dx_out = _mm256_mul_ps(x_in_dx_out, vscale);
}

// XY <-> HEXID without rotation

inline __m256i avx2_xy_to_hexid_f(const __m256 x, const __m256 y)
{
  __m256i u;
  __m256i v;
  avx2_xy_to_uv_f(x, y, u, v);
  return avx2_uv_to_hexid(u, v);
}

inline __m256i avx2_xy_to_hexid_with_remainder_f(__m256& x_in_dx_out, __m256& y_in_dy_out)
{
  __m256i u;
  __m256i v;
  avx2_xy_to_uv_with_remainder_f(x_in_dx_out, y_in_dy_out, u, v);
  return avx2_uv_to_hexid(u, v);
}

inline __m256i avx2_xy_to_hexid_with_remainder_f(__m256& x_in_dx_out, __m256& y_in_dy_out, bool clockwise)
{
  __m256i u;
  __m256i v;
  avx2_xy_to_uv_with_remainder_f(x_in_dx_out, y_in_dy_out, u, v);
  return avx2_uv_to_hexid(u, v, clockwise);
}

inline void avx2_hexid_to_xy_f(__m256i hexid, __m256& x, __m256& y)
{
  __m256i u;
  __m256i v;
  avx2_hexid_to_uv(hexid, u, v);
  avx2_uv_to_xy_f(u,v,x,y);
}

inline void avx2_hexid_to_xy_f(__m256i hexid, __m256& x, __m256& y, bool clockwise)
{
  __m256i u;
  __m256i v;
  avx2_hexid_to_uv(hexid, u, v, clockwise);
  avx2_uv_to_xy_f(u,v,x,y);
}

// XY <-> HEXID with (shear-free) affine transformation

inline __m256i avx2_xy_trans_to_hexid_f(__m256 x, __m256 y,
  const float crot, const float srot, const float scale, const float dx = 0, const float dy = 0)
{
  __m256i u;
  __m256i v;
  avx2_xy_trans_to_uv_f(x, y, u, v, crot, srot, scale, dx, dy);
  return avx2_uv_to_hexid(u, v);
}

inline __m256i avx2_xy_trans_to_hexid_f(__m256 x, __m256 y, bool clockwise,
  const float crot, const float srot, const float scale, const float dx = 0, const float dy = 0)
{
  __m256i u;
  __m256i v;
  avx2_xy_trans_to_uv_f(x, y, u, v, crot, srot, scale, dx, dy);
  return avx2_uv_to_hexid(u, v, clockwise);
}

inline __m256i avx2_xy_trans_to_hexid_with_remainder_f(
  __m256& x_in_dx_out, __m256& y_in_dy_out,
  const float crot, const float srot, const float scale, const float dx = 0, const float dy = 0)
{
  __m256i u;
  __m256i v;
  avx2_xy_trans_to_uv_with_remainder_f(x_in_dx_out, y_in_dy_out, u, v,
    crot, srot, scale, dx, dy);
  return avx2_uv_to_hexid(u, v);
}

inline __m256i avx2_xy_trans_to_hexid_with_remainder_f(
  __m256& x_in_dx_out, __m256& y_in_dy_out, bool clockwise,
  const float crot, const float srot, const float scale, const float dx = 0, const float dy = 0)
{
  __m256i u;
  __m256i v;
  avx2_xy_trans_to_uv_with_remainder_f(x_in_dx_out, y_in_dy_out, u, v,
    crot, srot, scale, dx, dy);
  return avx2_uv_to_hexid(u, v, clockwise);
}

inline void avx2_hexid_to_xy_trans_f(__m256i hexid, __m256& x, __m256& y,
  const float crot, const float srot, const float scale, const float dx = 0, const float dy = 0)
{
  __m256i u;
  __m256i v;
  avx2_hexid_to_uv(hexid, u, v);
  avx2_uv_to_xy_trans_f(u, v, x, y, crot, srot, scale, dx, dy);
}

inline void avx2_hexid_to_xy_trans_f(__m256i hexid, __m256& x, __m256& y, bool clockwise,
  const float crot, const float srot, const float scale, const float dx = 0, const float dy = 0)
{
  __m256i u;
  __m256i v;
  avx2_hexid_to_uv(hexid, u, v, clockwise);
  avx2_uv_to_xy_trans_f(u, v, x, y, crot, srot, scale, dx, dy);
}

#undef _c_m256_vx

#endif // defined(__AVX2__) and defined(__FMA__)

unsigned test_avx2_positive_hexid_to_ringid_loop(unsigned hexid);
unsigned test_avx2_positive_hexid_to_ringid_root(unsigned hexid);
unsigned test_avx2_hexid_to_ringid(unsigned hexid);
unsigned test_avx2_ringid_to_nsites_contained(unsigned ringid);
void test_avx2_positive_hexid_to_ringid_segid_runid(unsigned hexid,
  unsigned& ringid, unsigned& segid, unsigned& runid);
void test_avx2_hexid_to_ringid_segid_runid(unsigned hexid,
  unsigned& ringid, unsigned& segid, unsigned& runid);
unsigned test_avx2_positive_ringid_segid_runid_to_hexid(
  unsigned ringid, unsigned segid, unsigned runid);
unsigned test_avx2_ringid_segid_runid_to_hexid(
  unsigned ringid, unsigned segid, unsigned runid);
unsigned test_avx2_uv_to_ringid(int u, int v);
void test_avx2_hexid_to_uv_ccw(unsigned hexid, int& u, int& v);
void test_avx2_hexid_to_uv_cw(unsigned hexid, int& u, int& v);
unsigned test_avx2_uv_to_hexid_ccw(int u, int v);
unsigned test_avx2_uv_to_hexid_cw(int u, int v);


void test_avx2_uv_to_xy_f(int u, int v, float& x, float& y);
void test_avx2_xy_to_uv_f(float x, float y, int& u, int& v);
void test_avx2_xy_to_uv_with_remainder_f(float& x_in_dx_out, float& y_in_dy_out,
  int& u, int& v);
void test_avx2_uv_to_xy_trans_f(int u, int v, float& x, float& y,
  float crot, float srot, float scale, float dx = 0, float dy = 0);
void test_avx2_xy_trans_to_uv_f(float x, float y, int& u, int& v,
  float crot, float srot, float scale, float dx = 0, float dy = 0);
void test_avx2_xy_trans_to_uv_with_remainder_f(float& x_in_dx_out, float& y_in_dy_out,
  int& u, int& v,
  float crot, float srot, float scale, float dx = 0, float dy = 0);

unsigned test_avx2_xy_to_hexid_f(float x, float y);
unsigned test_avx2_xy_to_hexid_with_remainder_f(float& x_in_dx_out, float& y_in_dy_out);
unsigned test_avx2_xy_to_hexid_with_remainder_f(float& x_in_dx_out, float& y_in_dy_out, bool clockwise);
void test_avx2_hexid_to_xy_f(unsigned hexid, float& x, float& y);
void test_avx2_hexid_to_xy_f(unsigned hexid, float& x, float& y, bool clockwise);

unsigned test_avx2_xy_trans_to_hexid_f(float x, float y,
  float crot, float srot, float scale, float dx = 0, float dy = 0);
unsigned test_avx2_xy_trans_to_hexid_f(float x, float y, bool clockwise,
  float crot, float srot, float scale, float dx = 0, float dy = 0);
unsigned test_avx2_xy_trans_to_hexid_with_remainder_f(float& x_in_dx_out,
  float& y_in_dy_out,
  float crot, float srot, float scale, float dx = 0, float dy = 0);
unsigned test_avx2_xy_trans_to_hexid_with_remainder_f(float& x_in_dx_out,
  float& y_in_dy_out, bool clockwise,
  float crot, float srot, float scale, float dx = 0, float dy = 0);
void test_avx2_hexid_to_xy_trans_f(unsigned hexid, float& x, float& y,
  float crot, float srot, float scale, float dx = 0, float dy = 0);
void test_avx2_hexid_to_xy_trans_f(unsigned hexid, float& x, float& y, bool clockwise,
  float crot, float srot, float scale, float dx = 0, float dy = 0);


} } } // namespace calin::math::hex_array
