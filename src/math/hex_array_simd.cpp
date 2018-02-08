/*

   calin/io/hex_array.cpp -- Stephen Fegan -- 2015-10-21

   Collection of functions which translate between hexagonal and Cartesian
   geometries, and provide other useful calculations for hex grids.

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, Ecole Polytechnique, CNRS/IN2P3

   This file is part of "calin"

   "calin" is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 or
   later, as published by the Free Software Foundation.

   "calin" is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

*/

#include <stdexcept>

#include <math/simd.hpp>
#include <math/hex_array_simd.hpp>

unsigned calin::math::hex_array::test_avx2_positive_hexid_to_ringid_loop(unsigned hexid)
{
#if defined(__AVX2__) and defined(__FMA__)
  __m256i vhexid = _mm256_set1_epi32(hexid);
  __m256i vringid = avx2_positive_hexid_to_ringid_loop(vhexid);
  return vringid[0];
#else
  throw std::runtime_error("AVX2 and FMA not available at compile time");
#endif
}

unsigned calin::math::hex_array::test_avx2_positive_hexid_to_ringid_root(unsigned hexid)
{
#if defined(__AVX2__) and defined(__FMA__)
  __m256i vhexid = _mm256_set1_epi32(hexid);
  __m256i vringid = avx2_positive_hexid_to_ringid_root(vhexid);
  return vringid[0];
#else
  throw std::runtime_error("AVX2 and FMA not available at compile time");
#endif
}

unsigned calin::math::hex_array::test_avx2_hexid_to_ringid(unsigned hexid)
{
#if defined(__AVX2__) and defined(__FMA__)
  __m256i vhexid = _mm256_set1_epi32(hexid);
  __m256i vringid = avx2_hexid_to_ringid(vhexid);
  return vringid[0];
#else
  throw std::runtime_error("AVX2 and FMA not available at compile time");
#endif
}

unsigned calin::math::hex_array::test_avx2_ringid_to_nsites_contained(unsigned ringid)
{
#if defined(__AVX2__) and defined(__FMA__)
  __m256i vringid = _mm256_set1_epi32(ringid);
  __m256i vnsites = avx2_ringid_to_nsites_contained(vringid);
  return vnsites[0];
#else
  throw std::runtime_error("AVX2 and FMA not available at compile time");
#endif
}

void calin::math::hex_array::test_avx2_positive_hexid_to_ringid_segid_runid(unsigned hexid,
  unsigned& ringid, unsigned& segid, unsigned& runid)
{
#if defined(__AVX2__) and defined(__FMA__)
  __m256i vhexid = _mm256_set1_epi32(hexid);
  __m256i vringid;
  __m256i vsegid;
  __m256i vrunid;
  avx2_positive_hexid_to_ringid_segid_runid(vhexid, vringid, vsegid, vrunid);
  ringid = vringid[0];
  segid = vsegid[0];
  runid = vrunid[0];
#else
  throw std::runtime_error("AVX2 and FMA not available at compile time");
#endif
}

void calin::math::hex_array::test_avx2_hexid_to_ringid_segid_runid(unsigned hexid,
  unsigned& ringid, unsigned& segid, unsigned& runid)
{
#if defined(__AVX2__) and defined(__FMA__)
  __m256i vhexid = _mm256_set1_epi32(hexid);
  __m256i vringid;
  __m256i vsegid;
  __m256i vrunid;
  avx2_hexid_to_ringid_segid_runid(vhexid, vringid, vsegid, vrunid);
  ringid = vringid[0];
  segid = vsegid[0];
  runid = vrunid[0];
#else
  throw std::runtime_error("AVX2 and FMA not available at compile time");
#endif
}

unsigned calin::math::hex_array::test_avx2_positive_ringid_segid_runid_to_hexid(
  unsigned ringid, unsigned segid, unsigned runid)
{
#if defined(__AVX2__) and defined(__FMA__)
  __m256i vringid = _mm256_set1_epi32(ringid);
  __m256i vsegid = _mm256_set1_epi32(segid);
  __m256i vrunid = _mm256_set1_epi32(runid);
  __m256i vhexid = avx2_positive_ringid_segid_runid_to_hexid(vringid, vsegid, vrunid);
  return vhexid[0];
#else
  throw std::runtime_error("AVX2 and FMA not available at compile time");
#endif
}

unsigned calin::math::hex_array::test_avx2_ringid_segid_runid_to_hexid(
  unsigned ringid, unsigned segid, unsigned runid)
{
#if defined(__AVX2__) and defined(__FMA__)
  __m256i vringid = _mm256_set1_epi32(ringid);
  __m256i vsegid = _mm256_set1_epi32(segid);
  __m256i vrunid = _mm256_set1_epi32(runid);
  __m256i vhexid = avx2_ringid_segid_runid_to_hexid(vringid, vsegid, vrunid);
  return vhexid[0];
#else
  throw std::runtime_error("AVX2 and FMA not available at compile time");
#endif
}

unsigned calin::math::hex_array::test_avx2_uv_to_ringid(int u, int v)
{
#if defined(__AVX2__) and defined(__FMA__)
  __m256i vu = _mm256_set1_epi32(u);
  __m256i vv = _mm256_set1_epi32(v);
  __m256i vringid = avx2_uv_to_ringid(vu, vv);
  return vringid[0];
#else
  throw std::runtime_error("AVX2 and FMA not available at compile time");
#endif
}

void calin::math::hex_array::test_avx2_hexid_to_uv_ccw(unsigned hexid, int& u, int& v)
{
#if defined(__AVX2__) and defined(__FMA__)
  __m256i vhexid = _mm256_set1_epi32(hexid);
  __m256i vu;
  __m256i vv;
  avx2_hexid_to_uv_ccw(vhexid, vu, vv);
  u = vu[0];
  v = vv[0];
#else
  throw std::runtime_error("AVX2 and FMA not available at compile time");
#endif
}

void calin::math::hex_array::test_avx2_hexid_to_uv_cw(unsigned hexid, int& u, int& v)
{
#if defined(__AVX2__) and defined(__FMA__)
  __m256i vhexid = _mm256_set1_epi32(hexid);
  __m256i vu;
  __m256i vv;
  avx2_hexid_to_uv_cw(vhexid, vu, vv);
  u = vu[0];
  v = vv[0];
#else
  throw std::runtime_error("AVX2 and FMA not available at compile time");
#endif
}

unsigned calin::math::hex_array::test_avx2_uv_to_hexid_ccw(int u, int v)
{
#if defined(__AVX2__) and defined(__FMA__)
  __m256i vu = _mm256_set1_epi32(u);
  __m256i vv = _mm256_set1_epi32(v);
  __m256i vhexid = avx2_uv_to_hexid_ccw(vu, vv);
  return vhexid[0];
#else
  throw std::runtime_error("AVX2 and FMA not available at compile time");
#endif
}

unsigned calin::math::hex_array::test_avx2_uv_to_hexid_cw(int u, int v)
{
#if defined(__AVX2__) and defined(__FMA__)
  __m256i vu = _mm256_set1_epi32(u);
  __m256i vv = _mm256_set1_epi32(v);
  __m256i vhexid = avx2_uv_to_hexid_cw(vu, vv);
  return vhexid[0];
#else
  throw std::runtime_error("AVX2 and FMA not available at compile time");
#endif
}






void calin::math::hex_array::test_avx2_uv_to_xy_f(int u, int v, float& x, float& y)
{
#if defined(__AVX2__) and defined(__FMA__)
  __m256i vu = _mm256_set1_epi32(u);
  __m256i vv = _mm256_set1_epi32(v);
  __m256 vx;
  __m256 vy;
  avx2_uv_to_xy_f(vu, vv, vx, vy);
  x = vx[0];
  y = vy[0];
#else
  throw std::runtime_error("AVX2 and FMA not available at compile time");
#endif
}

void calin::math::hex_array::test_avx2_xy_to_uv_f(float x, float y, int& u, int& v)
{
#if defined(__AVX2__) and defined(__FMA__)
  __m256 vx = _mm256_set1_ps(x);
  __m256 vy = _mm256_set1_ps(y);
  __m256i vu;
  __m256i vv;
  avx2_xy_to_uv_f(vx, vy, vu, vv);
  u = vu[0];
  v = vv[0];
#else
  throw std::runtime_error("AVX2 and FMA not available at compile time");
#endif
}

void calin::math::hex_array::
test_avx2_xy_to_uv_with_remainder_f(float& x_in_dx_out, float& y_in_dy_out, int& u, int& v)
{
  #if defined(__AVX2__) and defined(__FMA__)
  __m256 vx_in_dx_out = _mm256_set1_ps(x_in_dx_out);
  __m256 vy_in_dy_out = _mm256_set1_ps(y_in_dy_out);
  __m256i vu;
  __m256i vv;
  avx2_xy_to_uv_f(vx_in_dx_out, vy_in_dy_out, vu, vv);
  u = vu[0];
  v = vv[0];
  x_in_dx_out = vx_in_dx_out[0];
  y_in_dy_out = vy_in_dy_out[0];
#else
  throw std::runtime_error("AVX2 and FMA not available at compile time");
#endif
}

void calin::math::hex_array::
test_avx2_uv_to_xy_trans_f(int u, int v, float& x, float& y,
  float crot, float srot, float scale, float dx, float dy)
{
#if defined(__AVX2__) and defined(__FMA__)
  __m256i vu = _mm256_set1_epi32(u);
  __m256i vv = _mm256_set1_epi32(v);
  __m256 vx;
  __m256 vy;
  avx2_uv_to_xy_trans_f(vu, vv, vx, vy, crot, srot, scale, dx, dy);
  x = vx[0];
  y = vy[0];
#else
  throw std::runtime_error("AVX2 and FMA not available at compile time");
#endif
}

void calin::math::hex_array::
test_avx2_xy_trans_to_uv_f(float x, float y, int& u, int& v,
  float crot, float srot, float scale, float dx, float dy)
{
#if defined(__AVX2__) and defined(__FMA__)
  __m256 vx = _mm256_set1_ps(x);
  __m256 vy = _mm256_set1_ps(y);
  __m256i vu;
  __m256i vv;
  avx2_xy_trans_to_uv_f(vx, vy, vu, vv, crot, srot, scale, dx, dy);
  u = vu[0];
  v = vv[0];
#else
  throw std::runtime_error("AVX2 and FMA not available at compile time");
#endif
}

void calin::math::hex_array::
test_avx2_xy_trans_to_uv_with_remainder_f(float& x_in_dx_out, float& y_in_dy_out,
  int& u, int& v, float crot, float srot, float scale, float dx, float dy)
{
#if defined(__AVX2__) and defined(__FMA__)
  __m256 vx_in_dx_out = _mm256_set1_ps(x_in_dx_out);
  __m256 vy_in_dy_out = _mm256_set1_ps(y_in_dy_out);
  __m256i vu;
  __m256i vv;
  avx2_xy_trans_to_uv_with_remainder_f(vx_in_dx_out, vy_in_dy_out, vu, vv, crot, srot, scale, dx, dy);
  u = vu[0];
  v = vv[0];
  x_in_dx_out = vx_in_dx_out[0];
  y_in_dy_out = vy_in_dy_out[0];
#else
  throw std::runtime_error("AVX2 and FMA not available at compile time");
#endif
}

unsigned calin::math::hex_array::test_avx2_xy_to_hexid_f(float x, float y)
{
#if defined(__AVX2__) and defined(__FMA__)
  __m256 vx = _mm256_set1_ps(x);
  __m256 vy = _mm256_set1_ps(y);
  __m256i vhexid;
  vhexid = avx2_xy_to_hexid_f(vx, vy);
  return vhexid[0];
#else
  throw std::runtime_error("AVX2 and FMA not available at compile time");
#endif
}

unsigned calin::math::hex_array::
test_avx2_xy_to_hexid_with_remainder_f(float& x_in_dx_out, float& y_in_dy_out)
{
#if defined(__AVX2__) and defined(__FMA__)
  __m256 vx_in_dx_out = _mm256_set1_ps(x_in_dx_out);
  __m256 vy_in_dy_out = _mm256_set1_ps(y_in_dy_out);
  __m256i vhexid;
  vhexid = avx2_xy_to_hexid_with_remainder_f(vx_in_dx_out, vy_in_dy_out);
  x_in_dx_out = vx_in_dx_out[0];
  y_in_dy_out = vy_in_dy_out[0];
  return vhexid[0];
#else
  throw std::runtime_error("AVX2 and FMA not available at compile time");
#endif
}

unsigned calin::math::hex_array::
test_avx2_xy_to_hexid_with_remainder_f(float& x_in_dx_out, float& y_in_dy_out, bool clockwise)
{
#if defined(__AVX2__) and defined(__FMA__)
  __m256 vx_in_dx_out = _mm256_set1_ps(x_in_dx_out);
  __m256 vy_in_dy_out = _mm256_set1_ps(y_in_dy_out);
  __m256i vhexid;
  vhexid = avx2_xy_to_hexid_with_remainder_f(vx_in_dx_out, vy_in_dy_out, clockwise);
  x_in_dx_out = vx_in_dx_out[0];
  y_in_dy_out = vy_in_dy_out[0];
  return vhexid[0];
#else
  throw std::runtime_error("AVX2 and FMA not available at compile time");
#endif
}

void calin::math::hex_array::
test_avx2_hexid_to_xy_f(unsigned hexid, float& x, float& y)
{
#if defined(__AVX2__) and defined(__FMA__)
__m256i vhexid = _mm256_set1_epi32(hexid);
  __m256 vx;
  __m256 vy;
  avx2_hexid_to_xy_f(vhexid, vx, vy);
  x = vx[0];
  y = vy[0];
#else
  throw std::runtime_error("AVX2 and FMA not available at compile time");
#endif
}

void calin::math::hex_array::
test_avx2_hexid_to_xy_f(unsigned hexid, float& x, float& y, bool clockwise)
{
#if defined(__AVX2__) and defined(__FMA__)
__m256i vhexid = _mm256_set1_epi32(hexid);
  __m256 vx;
  __m256 vy;
  avx2_hexid_to_xy_f(vhexid, vx, vy, clockwise);
  x = vx[0];
  y = vy[0];
#else
  throw std::runtime_error("AVX2 and FMA not available at compile time");
#endif
}
