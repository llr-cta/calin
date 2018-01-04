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

#if defined(__AVX2__)
#include <immintrin.h>
#endif // defined(__AVX2__)

namespace calin { namespace math { namespace simd {

constexpr float _ps0=7.8539816171e-01;
constexpr float _ps1=-8.0745451510e-02;
constexpr float _ps2=2.4900604274e-03;
constexpr float _ps3=-3.5995584999e-05;

constexpr float _pc0=-3.0842513735e-01;
constexpr float _pc1=1.5854339023e-02;
constexpr float _pc2=-3.2596461752e-04;
constexpr float _pc3=3.5445940383e-06;

#if defined(__AVX2__) && defined(__FMA__)

// AVX2 implementation of single precision floating point sin & cos functions
// valid over the domain -pi/4 (x=-1.0) to pi/4 (x=+1.0)
void avx2_sincosf_domain_pi_4_poly3(const __m256 x, __m256& s, __m256& c,
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

// AVX2 implementation of single precision floating point sin & cos functions
// valid over the domain -pi (x=-4.0) to pi (x=+4.0). No checking on the domain
// is
void avx2_sincosf_domain_pi_poly3(const __m256 x, __m256& s, __m256& c)
{
  // xq gives integer quadrent id of argument as either -4,-2,0,2,4
  // (the final quadrent therefore has id +/-4)
  __m256i xq = _mm256_cvtps_epi32(_mm256_add_ps(x,_mm256_set1_ps(0.5)));
  xq = _mm256_slli_epi32(_mm256_srli_epi32(xq,1),1);

  // Argument to pi/4 domain function is residual after subtraction of quadrant
  // id - this puts argument in domain -1.0 to 1.0
  __m256 xr = _mm256_sub_ps(x, _mm256_cvtepi32_ps(xq));

  // Call reduced domain sincos function with default polynomial arguments
  avx2_sincosf_domain_pi_4_poly3(xr, s, c);

  // Blend compnents calculated with the sine and cosine polynomials based
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

  // For cos we must now negate for quadtants -4, 2 and 4 - this arranges
  // for bit#3 of the quadran id to be XORd with bit#2, and then shifted into
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
