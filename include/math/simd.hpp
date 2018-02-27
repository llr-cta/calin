/*

   calin/math/simd.hpp -- Stephen Fegan -- 2018-01-04

   SIMD functions

   Copyright 2018, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#pragma once

#if defined(__AVX__)
#include <immintrin.h>
#endif // defined(__AVX__)

namespace calin { namespace math { namespace simd {

#if defined(__AVX__)
#define _CALIN_DO_ONE_AVX_SWIZZLE_U16(i) \
  tmp =  _mm256_unpackhi_epi16(*(data+i), *(data+i+1)); \
  *(data+i) =  _mm256_unpacklo_epi16(*(data+i), *(data+i+1)); \
  *(data+i+1) = tmp

#define _CALIN_DO_ONE_AVX_SWIZZLE_U32(i) \
  tmp =  _mm256_unpackhi_epi32(*(data+i), *(data+i+2)); \
  *(data+i) =  _mm256_unpacklo_epi32(*(data+i), *(data+i+2)); \
  *(data+i+2) = tmp

#define _CALIN_DO_ONE_AVX_SWIZZLE_U64(i) \
  tmp =  _mm256_unpackhi_epi64(*(data+i), *(data+i+4)); \
  *(data+i) =  _mm256_unpacklo_epi64(*(data+i), *(data+i+4)); \
  *(data+i+4) = tmp

#define _CALIN_DO_ONE_AVX_SWIZZLE_U128(i) \
  tmp =  _mm256_permute2x128_si256(*(data+i), *(data+i+8), 0x31); \
  *(data+i) =  _mm256_permute2x128_si256(*(data+i), *(data+i+8), 0x20); \
  *(data+i+8) = tmp

// Swizzle a array of U16 structs (aos) to struct of U16 arrays - AVX version
inline void avx_m256_swizzle_u16(__m256i* data)
{
  __m256i tmp;
  _CALIN_DO_ONE_AVX_SWIZZLE_U16(0);
  _CALIN_DO_ONE_AVX_SWIZZLE_U16(2);
  _CALIN_DO_ONE_AVX_SWIZZLE_U16(4);
  _CALIN_DO_ONE_AVX_SWIZZLE_U16(6);
  _CALIN_DO_ONE_AVX_SWIZZLE_U16(8);
  _CALIN_DO_ONE_AVX_SWIZZLE_U16(10);
  _CALIN_DO_ONE_AVX_SWIZZLE_U16(12);
  _CALIN_DO_ONE_AVX_SWIZZLE_U16(14);

  _CALIN_DO_ONE_AVX_SWIZZLE_U32(0);
  _CALIN_DO_ONE_AVX_SWIZZLE_U32(1);
  _CALIN_DO_ONE_AVX_SWIZZLE_U32(4);
  _CALIN_DO_ONE_AVX_SWIZZLE_U32(5);
  _CALIN_DO_ONE_AVX_SWIZZLE_U32(8);
  _CALIN_DO_ONE_AVX_SWIZZLE_U32(9);
  _CALIN_DO_ONE_AVX_SWIZZLE_U32(12);
  _CALIN_DO_ONE_AVX_SWIZZLE_U32(13);

  _CALIN_DO_ONE_AVX_SWIZZLE_U64(0);
  _CALIN_DO_ONE_AVX_SWIZZLE_U64(1);
  _CALIN_DO_ONE_AVX_SWIZZLE_U64(2);
  _CALIN_DO_ONE_AVX_SWIZZLE_U64(3);
  _CALIN_DO_ONE_AVX_SWIZZLE_U64(8);
  _CALIN_DO_ONE_AVX_SWIZZLE_U64(9);
  _CALIN_DO_ONE_AVX_SWIZZLE_U64(10);
  _CALIN_DO_ONE_AVX_SWIZZLE_U64(11);

  _CALIN_DO_ONE_AVX_SWIZZLE_U128(0);
  _CALIN_DO_ONE_AVX_SWIZZLE_U128(1);
  _CALIN_DO_ONE_AVX_SWIZZLE_U128(2);
  _CALIN_DO_ONE_AVX_SWIZZLE_U128(3);
  _CALIN_DO_ONE_AVX_SWIZZLE_U128(4);
  _CALIN_DO_ONE_AVX_SWIZZLE_U128(5);
  _CALIN_DO_ONE_AVX_SWIZZLE_U128(6);
  _CALIN_DO_ONE_AVX_SWIZZLE_U128(7);

  std::swap(data[1], data[4]);
  std::swap(data[3], data[6]);
  std::swap(data[9], data[12]);
  std::swap(data[11], data[14]);
}
#endif // defined(__AVX__)

#if defined(__AVX2__)
#define _CALIN_DO_ONE_AVX2_SWIZZLE_U128(i) \
  tmp =  _mm256_permute2x128_si256(*(data+i), *(data+i+8), 0x31); \
  *(data+i) =  _mm256_permute2x128_si256(*(data+i), *(data+i+8), 0x20); \
  *(data+i+8) = tmp

// Swizzle a array of U16 structs (aos) to struct of U16 arrays - AVX2 version
inline void avx2_m256_swizzle_u16(__m256i* data)
{
  __m256i tmp;
  _CALIN_DO_ONE_AVX_SWIZZLE_U16(0);
  _CALIN_DO_ONE_AVX_SWIZZLE_U16(2);
  _CALIN_DO_ONE_AVX_SWIZZLE_U16(4);
  _CALIN_DO_ONE_AVX_SWIZZLE_U16(6);
  _CALIN_DO_ONE_AVX_SWIZZLE_U16(8);
  _CALIN_DO_ONE_AVX_SWIZZLE_U16(10);
  _CALIN_DO_ONE_AVX_SWIZZLE_U16(12);
  _CALIN_DO_ONE_AVX_SWIZZLE_U16(14);

  _CALIN_DO_ONE_AVX_SWIZZLE_U32(0);
  _CALIN_DO_ONE_AVX_SWIZZLE_U32(1);
  _CALIN_DO_ONE_AVX_SWIZZLE_U32(4);
  _CALIN_DO_ONE_AVX_SWIZZLE_U32(5);
  _CALIN_DO_ONE_AVX_SWIZZLE_U32(8);
  _CALIN_DO_ONE_AVX_SWIZZLE_U32(9);
  _CALIN_DO_ONE_AVX_SWIZZLE_U32(12);
  _CALIN_DO_ONE_AVX_SWIZZLE_U32(13);

  _CALIN_DO_ONE_AVX_SWIZZLE_U64(0);
  _CALIN_DO_ONE_AVX_SWIZZLE_U64(1);
  _CALIN_DO_ONE_AVX_SWIZZLE_U64(2);
  _CALIN_DO_ONE_AVX_SWIZZLE_U64(3);
  _CALIN_DO_ONE_AVX_SWIZZLE_U64(8);
  _CALIN_DO_ONE_AVX_SWIZZLE_U64(9);
  _CALIN_DO_ONE_AVX_SWIZZLE_U64(10);
  _CALIN_DO_ONE_AVX_SWIZZLE_U64(11);

  _CALIN_DO_ONE_AVX2_SWIZZLE_U128(0);
  _CALIN_DO_ONE_AVX2_SWIZZLE_U128(1);
  _CALIN_DO_ONE_AVX2_SWIZZLE_U128(2);
  _CALIN_DO_ONE_AVX2_SWIZZLE_U128(3);
  _CALIN_DO_ONE_AVX2_SWIZZLE_U128(4);
  _CALIN_DO_ONE_AVX2_SWIZZLE_U128(5);
  _CALIN_DO_ONE_AVX2_SWIZZLE_U128(6);
  _CALIN_DO_ONE_AVX2_SWIZZLE_U128(7);

  std::swap(data[1], data[4]);
  std::swap(data[3], data[6]);
  std::swap(data[9], data[12]);
  std::swap(data[11], data[14]);
}
#endif // defined(__AVX2__)

#if defined(__AVX__)

inline void _CALIN_DO_ONE_AVX_SWIZZLE_FLT32(__m256& x, __m256& y)
{
  __m256 t;
  t = _mm256_unpackhi_ps(x, y);
  x = _mm256_unpacklo_ps(x, y);
  y = t;
}

inline void _CALIN_DO_ONE_AVX_SWIZZLE_FLT64(__m256& x, __m256& y)
{
  __m256d t;
  t = _mm256_unpackhi_pd(_mm256_castps_pd(x), _mm256_castps_pd(y));
  x = _mm256_castpd_ps(_mm256_unpacklo_pd(_mm256_castps_pd(x), _mm256_castps_pd(y)));
  y = _mm256_castpd_ps(t);
}

inline void _CALIN_DO_ONE_AVX_SWIZZLE_FLT128_HILO(__m256& x, __m256& y)
{
  __m128 t;
  t = _mm256_extractf128_ps(x, 1);
  x = _mm256_insertf128_ps(x, _mm256_extractf128_ps(y, 0), 1);
  y = _mm256_insertf128_ps(y, t, 0);
}

inline void _CALIN_DO_ONE_AVX_SWIZZLE_FLT128_LOLO(__m256& x, __m256& y)
{
  __m128 t;
  t = _mm256_extractf128_ps(x, 0);
  x = _mm256_insertf128_ps(x, _mm256_extractf128_ps(y, 0), 0);
  y = _mm256_insertf128_ps(y, t, 0);
}

inline void _CALIN_DO_ONE_AVX_SWIZZLE_FLT128_HIHI(__m256& x, __m256& y)
{
  __m128 t;
  t = _mm256_extractf128_ps(x, 1);
  x = _mm256_insertf128_ps(x, _mm256_extractf128_ps(y, 1), 1);
  y = _mm256_insertf128_ps(y, t, 1);
}

inline void transpose(__m256& a, __m256& b, __m256& c, __m256& d, __m256& e, __m256& f, __m256& g, __m256& h)
{
  _CALIN_DO_ONE_AVX_SWIZZLE_FLT32(a, b);
  _CALIN_DO_ONE_AVX_SWIZZLE_FLT32(c, d);
  _CALIN_DO_ONE_AVX_SWIZZLE_FLT32(e, f);
  _CALIN_DO_ONE_AVX_SWIZZLE_FLT32(g, h);

  _CALIN_DO_ONE_AVX_SWIZZLE_FLT64(a, c);
  _CALIN_DO_ONE_AVX_SWIZZLE_FLT64(b, d);
  _CALIN_DO_ONE_AVX_SWIZZLE_FLT64(e, g);
  _CALIN_DO_ONE_AVX_SWIZZLE_FLT64(f, h);

  _CALIN_DO_ONE_AVX_SWIZZLE_FLT128_HILO(a, e);
  _CALIN_DO_ONE_AVX_SWIZZLE_FLT128_HILO(b, g);
  _CALIN_DO_ONE_AVX_SWIZZLE_FLT128_HILO(c, f);
  _CALIN_DO_ONE_AVX_SWIZZLE_FLT128_HILO(d, h);
  _CALIN_DO_ONE_AVX_SWIZZLE_FLT128_LOLO(b, c);
  _CALIN_DO_ONE_AVX_SWIZZLE_FLT128_HIHI(f, g);
}

inline void _CALIN_DO_ONE_AVX_SWIZZLE_DBL64(__m256d& x, __m256d& y)
{
  __m256d t;
  t = _mm256_unpackhi_pd(x, y);
  x = _mm256_unpacklo_pd(x, y);
  y = t;
}

inline void _CALIN_DO_ONE_AVX_SWIZZLE_DBL128(__m256d& x, __m256d& y)
{
  __m128d t;
  t = _mm256_extractf128_pd(x, 1);
  x = _mm256_insertf128_pd(x, _mm256_extractf128_pd(y, 0), 1);
  y = _mm256_insertf128_pd(y, t, 0);
}

inline void transpose(__m256d& a, __m256d& b, __m256d& c, __m256d& d)
{
  _CALIN_DO_ONE_AVX_SWIZZLE_DBL64(a, b);
  _CALIN_DO_ONE_AVX_SWIZZLE_DBL64(c, d);
  _CALIN_DO_ONE_AVX_SWIZZLE_DBL128(a, c);
  _CALIN_DO_ONE_AVX_SWIZZLE_DBL128(b, d);
}

#endif

// First fit direct
// constexpr float _ps0=7.8539816171e-01;
// constexpr float _ps1=-8.0745451510e-02;
// constexpr float _ps2=2.4900604274e-03;
// constexpr float _ps3=-3.5995584999e-05;

// constexpr float _pc0=-3.0842513735e-01;
// constexpr float _pc1=1.5854339023e-02;
// constexpr float _pc2=-3.2596461752e-04;
// constexpr float _pc3=3.5445940383e-06;

// Freeze sin coefficient at pi/4
// constexpr float _ps0=7.8539816340e-01;
// constexpr float _ps1=-8.0745505976e-02;
// constexpr float _ps2=2.4902556507e-03;
// constexpr float _ps3=-3.6153601377e-05;

// MINMAX
constexpr float _ps0 = 7.8539816098e-01;
constexpr float _ps1 = -8.0745434831e-02;
constexpr float _ps2 = 2.4900071298e-03;
constexpr float _ps3 = -3.5954509280e-05;

constexpr float _pc0 = -3.0842513734e-01;
constexpr float _pc1 = 1.5854338156e-02;
constexpr float _pc2 = -3.2596140215e-04;
constexpr float _pc3 = 3.5419672354e-06;

#define CALIN_MM256PS_CONST(name, val) \
static const float name[8] __attribute__((aligned(32))) = {(val),(val),(val),(val),(val),(val),(val),(val)}

#if defined(__AVX__)
inline __m256 c_m256(const float (&c)[8])
{
  return *(const __m256 *)&c;
}
#endif

#if defined(__AVX2__) && defined(__FMA__)

// AVX2 implementation of single precision floating point sin & cos functions
// valid over the domain -pi/4 (x=-1.0) to pi/4 (x=+1.0). Different polynomial
// coefficients can be given if desired (for testing!)
inline void avx2_sincosf_domain_pi_4_poly3(const __m256 x, __m256& s, __m256& c,
  const float ps0=_ps0, const float ps1=_ps1, const float ps2=_ps2, const float ps3=_ps3,
  const float pc0=_pc0, const float pc1=_pc1, const float pc2=_pc2, const float pc3=_pc3)
{
  __m256 xx = _mm256_mul_ps(x,x);

  s = _mm256_set1_ps(ps3);
  s = _mm256_fmadd_ps(xx, s, _mm256_set1_ps(ps2));
  s = _mm256_fmadd_ps(xx, s, _mm256_set1_ps(ps1));
  s = _mm256_fmadd_ps(xx, s, _mm256_set1_ps(ps0));
  s = _mm256_mul_ps(x, s);

  c = _mm256_set1_ps(pc3);
  c = _mm256_fmadd_ps(xx, c, _mm256_set1_ps(pc2));
  c = _mm256_fmadd_ps(xx, c, _mm256_set1_ps(pc1));
  c = _mm256_fmadd_ps(xx, c, _mm256_set1_ps(pc0));
  c = _mm256_fmadd_ps(xx, c, _mm256_set1_ps(1.0));
}

CALIN_MM256PS_CONST(_c_m256_ps0, _ps0);
CALIN_MM256PS_CONST(_c_m256_ps1, _ps1);
CALIN_MM256PS_CONST(_c_m256_ps2, _ps2);
CALIN_MM256PS_CONST(_c_m256_ps3, _ps3);
CALIN_MM256PS_CONST(_c_m256_pc0, _pc0);
CALIN_MM256PS_CONST(_c_m256_pc1, _pc1);
CALIN_MM256PS_CONST(_c_m256_pc2, _pc2);
CALIN_MM256PS_CONST(_c_m256_pc3, _pc3);
CALIN_MM256PS_CONST(_c_m256_one, 1.0);

inline void avx2_sincosf_domain_pi_4_poly3_def(const __m256 x, __m256& s, __m256& c)
{
  __m256 xx = _mm256_mul_ps(x,x);

  s = c_m256(_c_m256_ps3);
  s = _mm256_fmadd_ps(xx, s, c_m256(_c_m256_ps2));
  s = _mm256_fmadd_ps(xx, s, c_m256(_c_m256_ps1));
  s = _mm256_fmadd_ps(xx, s, c_m256(_c_m256_ps0));
  s = _mm256_mul_ps(x, s);

  c = c_m256(_c_m256_pc3);
  c = _mm256_fmadd_ps(xx, c, c_m256(_c_m256_pc2));
  c = _mm256_fmadd_ps(xx, c, c_m256(_c_m256_pc1));
  c = _mm256_fmadd_ps(xx, c, c_m256(_c_m256_pc0));
  c = _mm256_fmadd_ps(xx, c, c_m256(_c_m256_one));
}


// AVX2 implementation of single precision floating point sin & cos functions
// valid over the domain -pi (x=-4.0) to pi (x=+4.0). No checking on the domain
// is performed
inline void avx2_sincosf_domain_pi_poly3(const __m256 x, __m256& s, __m256& c)
{
  // xq gives integer quadrent id of argument as either -4,-2,0,2,4
  // (the final quadrent therefore has a split id of -4 and +4)
  __m256i xq = _mm256_cvtps_epi32(_mm256_add_ps(x,_mm256_set1_ps(0.5)));
  xq = _mm256_slli_epi32(_mm256_srli_epi32(xq,1),1);

  // Calculate argument to give to pi/4 domain sincos function above - residual
  // after subtraction of quadrant id - this puts argument in domain -1.0 to 1.0
  __m256 xr = _mm256_sub_ps(x, _mm256_cvtepi32_ps(xq));

  // Call reduced domain sincos function with default polynomial arguments
  avx2_sincosf_domain_pi_4_poly3(xr, s, c);

  // Blend components calculated with the sine and cosine polynomials based
  // on the quadrant -- we swap sin and cos based on bit#2 - swapping
  // thereby only quadrants +2 and -2
  __m256i mask = _mm256_slli_epi32(xq,30);
  xr = _mm256_blendv_ps(s,c,_mm256_castsi256_ps(mask));
  c = _mm256_blendv_ps(c,s,_mm256_castsi256_ps(mask));

  // For sin we must now negate for quadtants -4, -2 and 4 - this arranges
  // for bit#3 of the quadrant id to be shifted into bit#32 where we use it
  // to flip (XOR) the sign of the float
  mask = _mm256_slli_epi32(_mm256_srli_epi32(xq,2),31);
  s = _mm256_xor_ps(xr, _mm256_castsi256_ps(mask));

  // For cos we must now negate for quadrants -4, 2 and 4 - this arranges
  // for bit#3 of the quadrant id to be XORd with bit#2, and then shifted into
  // bit#32 where we use it to flip (xor) the sign of the float
  mask = _mm256_xor_si256(_mm256_slli_epi32(xq,30), mask);
  c = _mm256_xor_ps(c, _mm256_castsi256_ps(mask));
}

#endif // defined __AVX2__

float test_sin_avx2_sincosf_domain_pi_4_poly3(const float x,
  const float ps0=_ps0, const float ps1=_ps1, const float ps2=_ps2, const float ps3=_ps3);

float test_cos_avx2_sincosf_domain_pi_4_poly3(const float x,
  const float pc0=_pc0, const float pc1=_pc1, const float pc2=_pc2, const float pc3=_pc3);

float test_sin_avx2_sincosf_domain_pi_poly3(const float x);

float test_cos_avx2_sincosf_domain_pi_poly3(const float x);

} } } // namespace calin::math::simd
