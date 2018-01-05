/*

   calin/unit_tests/math/test_simd.cpp -- Stephen Fegan -- 2018-01-05

   Unit tests for SIMD

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

#include <iostream>
#include <iomanip>
#include <gtest/gtest.h>
#include <vector>
#include <tuple>
#include <cmath>

#include <math/rng.hpp>
#include <math/simd.hpp>

using calin::math::rng::NR3_AVX2_RNGCore;
using calin::math::rng::NR3RNGCore;
using namespace calin::math::simd;

constexpr unsigned N100M = 100000000;

TEST(TestSIMD, SpeedTest100M_Random64SinCos32)
{
  float sum_s = 0;
  float sum_c = 0;
  uint64_t seed = calin::math::rng::RNG::uint64_from_random_device();
  NR3RNGCore core(seed);
  for(unsigned i=0;i<N100M/2;i++) {
    uint64_t x64 = core.uniform_uint64();
    float x = 2.328306437e-10 * 2.0 * M_PI * unsigned(x64);
    sum_s += ::sinf(x);
    sum_c += ::cosf(x);
    x = 2.328306437e-10 * 2.0 * M_PI * unsigned(x64>>32);
    sum_s += ::sinf(x);
    sum_c += ::cosf(x);
  }

  EXPECT_GE(sum_s, -float(N100M));
  EXPECT_LE(sum_s, float(N100M));
  EXPECT_GE(sum_c, -float(N100M));
  EXPECT_LE(sum_c, float(N100M));
}

#ifdef CALIN_HAS_NR3_AVX2_RNGCORE
TEST(TestSIMD, SpeedTest100M_Random256SinCos32)
{
  float sum_s = 0;
  float sum_c = 0;
  uint64_t seed = calin::math::rng::RNG::uint64_from_random_device();
  NR3_AVX2_RNGCore core(seed);
  for(unsigned i=0;i<N100M/8;i++) {
    float x[8] __attribute__((aligned(32)));
    __m256 vx = core.uniform_zc_psvec256(2*M_PI);
    _mm256_store_ps(x, vx);
    sum_s += ::sinf(x[0]);
    sum_c += ::cosf(x[0]);
    sum_s += ::sinf(x[1]);
    sum_c += ::cosf(x[1]);
    sum_s += ::sinf(x[2]);
    sum_c += ::cosf(x[2]);
    sum_s += ::sinf(x[3]);
    sum_c += ::cosf(x[3]);
    sum_s += ::sinf(x[4]);
    sum_c += ::cosf(x[4]);
    sum_s += ::sinf(x[5]);
    sum_c += ::cosf(x[5]);
    sum_s += ::sinf(x[6]);
    sum_c += ::cosf(x[6]);
    sum_s += ::sinf(x[7]);
    sum_c += ::cosf(x[7]);
  }

  EXPECT_GE(sum_s, -float(N100M));
  EXPECT_LE(sum_s, float(N100M));
  EXPECT_GE(sum_c, -float(N100M));
  EXPECT_LE(sum_c, float(N100M));
}

TEST(TestSIMD, SpeedTest100M_Random256SinCos256)
{
  __m256 sum_s = _mm256_setzero_ps();
  __m256 sum_c = _mm256_setzero_ps();
  uint64_t seed = calin::math::rng::RNG::uint64_from_random_device();
  NR3_AVX2_RNGCore core(seed);
  for(unsigned i=0;i<N100M/8;i++) {
    __m256 x = core.uniform_zc_psvec256(8.0);
    __m256 s;
    __m256 c;
    avx2_sincosf_domain_pi_poly3(x, s, c);
    sum_s = _mm256_add_ps(sum_s, s);
    sum_c = _mm256_add_ps(sum_c, c);
  }
  EXPECT_GE(reinterpret_cast<float*>(&sum_s)[0], -float(N100M));
  EXPECT_LE(reinterpret_cast<float*>(&sum_s)[0], float(N100M));
  EXPECT_GE(reinterpret_cast<float*>(&sum_c)[0], -float(N100M));
  EXPECT_LE(reinterpret_cast<float*>(&sum_c)[0], float(N100M));
}
#endif // defined CALIN_HAS_NR3_AVX2_RNGCORE

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
