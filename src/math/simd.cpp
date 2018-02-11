/*

   calin/math/simd.cpp -- Stephen Fegan -- 2018-01-04

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

#include <stdexcept>

#include <math/simd.hpp>

float calin::math::simd::test_sin_avx2_sincosf_domain_pi_4_poly3(const float x,
  const float ps0, const float ps1, const float ps2, const float ps3)
{
#if defined(__AVX2__) && defined(__FMA__)
  __m256 vx = _mm256_set1_ps(x);
  __m256 vs;
  __m256 vc;
  avx2_sincosf_domain_pi_4_poly3(vx, vs, vc, ps0, ps1, ps2, ps3, _pc0, _pc1, _pc2, _pc3);
  float s[8];
  _mm256_store_ps(s, vs);
  return s[0];
#else
  throw std::runtime_error("AVX2 not available at compile time");
#endif
}

float calin::math::simd::test_cos_avx2_sincosf_domain_pi_4_poly3(const float x,
  const float pc0, const float pc1, const float pc2, const float pc3)
{
#if defined(__AVX2__) && defined(__FMA__)
  __m256 vx = _mm256_set1_ps(x);
  __m256 vs;
  __m256 vc;
  avx2_sincosf_domain_pi_4_poly3(vx, vs, vc, _ps0, _ps1, _ps2, _ps3, pc0, pc1, pc2, pc3);
  float c[8];
  _mm256_store_ps(c, vc);
  return c[0];
#else
  throw std::runtime_error("AVX2 not available at compile time");
#endif
}

float calin::math::simd::test_sin_avx2_sincosf_domain_pi_poly3(const float x)
{
#if defined(__AVX2__) && defined(__FMA__)
  __m256 vx = _mm256_set1_ps(x);
  __m256 vs;
  __m256 vc;
  avx2_sincosf_domain_pi_poly3(vx, vs, vc);
  float s[8];
  _mm256_store_ps(s, vs);
  return s[0];
#else
  throw std::runtime_error("AVX2 not available at compile time");
#endif
}

float calin::math::simd::test_cos_avx2_sincosf_domain_pi_poly3(const float x)
{
#if defined(__AVX2__) && defined(__FMA__)
  __m256 vx = _mm256_set1_ps(x);
  __m256 vs;
  __m256 vc;
  avx2_sincosf_domain_pi_poly3(vx, vs, vc);
  float c[8];
  _mm256_store_ps(c, vc);
  return c[0];
#else
  throw std::runtime_error("AVX2 not available at compile time");
#endif
}
