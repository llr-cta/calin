/*

   calin/unit_tests/math/test_rng_vcl.cpp -- Stephen Fegan -- 2018-08-14

   Unit tests for VCL RNG

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

#include <iostream>
#include <iomanip>
#include <gtest/gtest.h>

#include <math/accumulator.hpp>
#include <math/rng_vcl.hpp>
#include <util/vcl.hpp>

using namespace calin::math::accumulator;
using namespace calin::math::rng;
using namespace calin::util::vcl;

TEST(NR3SpeedTestRNG, Scalar_1G_UInt64)
{
  uint64_t sum = 0ULL;
  uint64_t seed = RNG::uint64_from_random_device();
  NR3RNGCore core(seed);
  for(unsigned i=0;i<1000000000;i++)
    sum += core.uniform_uint64();
  EXPECT_GE(sum, 0ULL);
}

#ifdef CALIN_HAS_NR3_AVX2_RNGCORE
TEST(NR3SpeedTestRNG, AVX2_1G_UInt64)
{
  uint64_t sum = 0ULL;
  uint64_t seed = RNG::uint64_from_random_device();
  NR3_AVX2_RNGCore core(seed);
  for(unsigned i=0;i<1000000000;i++)
    sum += core.uniform_uint64();
  EXPECT_GE(sum, 0ULL);
}

TEST(NR3SpeedTestRNG, AVX2_VEC_1G_UInt64)
{
  __m256i sum = _mm256_setzero_si256();
  uint64_t seed = RNG::uint64_from_random_device();
  NR3_AVX2_RNGCore core(seed);
  for(unsigned i=0;i<250000000;i++)
    sum = _mm256_add_epi64(sum, core.uniform_m256i());
  EXPECT_GE(sum[0], 0ULL);
  EXPECT_GE(sum[1], 0ULL);
  EXPECT_GE(sum[2], 0ULL);
  EXPECT_GE(sum[3], 0ULL);
}
#endif // defined CALIN_HAS_NR3_AVX2_RNGCORE

TEST(NR3SpeedTestRNG, Scalar_1G_Float)
{
  float sum = 0;
  uint64_t seed = RNG::uint64_from_random_device();
  NR3RNGCore core(seed);
  for(unsigned i=0;i<1000000000;i++)
    sum += 2.328306437e-10 * float(unsigned(core.uniform_uint64()));
  EXPECT_GE(sum,0.0);
}

TEST(NR3SpeedTestRNG, Scalar_1G_Double)
{
  double sum = 0;
  uint64_t seed = RNG::uint64_from_random_device();
  NR3RNGCore core(seed);
  for(unsigned i=0;i<1000000000;i++)
    sum += 5.42101086242752217E-20 * double(core.uniform_uint64());
  EXPECT_GE(sum,0.0);
}

#ifdef CALIN_HAS_NR3_AVX2_RNGCORE
TEST(NR3SpeedTestRNG, AVX2_VEC_1G_Float)
{
  __m256 sum = _mm256_setzero_ps();
  uint64_t seed = RNG::uint64_from_random_device();
  NR3_AVX2_RNGCore core(seed);
  for(unsigned i=0;i<125000000;i++)
    sum = _mm256_add_ps(sum, core.uniform_m256());
  EXPECT_GE(sum[0],0.0);
  EXPECT_GE(sum[1],0.0);
  EXPECT_GE(sum[2],0.0);
  EXPECT_GE(sum[3],0.0);
  EXPECT_GE(sum[4],0.0);
  EXPECT_GE(sum[5],0.0);
  EXPECT_GE(sum[6],0.0);
  EXPECT_GE(sum[7],0.0);
}
#endif // defined CALIN_HAS_NR3_AVX2_RNGCORE

template<typename VCLArchitecture> class VCLSpeedTestRNG :
  public testing::Test
{
public:
};

using ArchTypes = ::testing::Types<VCL128Architecture, VCL256Architecture, VCL512Architecture>;

TYPED_TEST_CASE(VCLSpeedTestRNG, ArchTypes);

TYPED_TEST(VCLSpeedTestRNG, VEC_1G_Uint64)
{
  uint64_t seed = RNG::uint64_from_random_device();
  NR3_VCLRNGCore<TypeParam> core(seed, __PRETTY_FUNCTION__, "core");
  const unsigned N = unsigned(UINT64_C(64000000000)/TypeParam::vec_bits);
  typename TypeParam::uint64_vt sum(0);
  for(unsigned i=0;i<N;i++) {
    sum += core.uniform_uint64();
  }
  EXPECT_TRUE(horizontal_and(sum >= UINT64_C(0)));
}

TYPED_TEST(VCLSpeedTestRNG, VEC_1G_Float)
{
  uint64_t seed = RNG::uint64_from_random_device();
  VCLRNG<TypeParam> core(seed, __PRETTY_FUNCTION__, "core");
  const unsigned N = unsigned(UINT64_C(32000000000)/TypeParam::vec_bits);
  typename TypeParam::float_vt sum(0);
  for(unsigned i=0;i<N;i++) {
    sum = core.uniform_float();
  }
  EXPECT_TRUE(horizontal_and(sum >= 0.0));
}

TYPED_TEST(VCLSpeedTestRNG, VEC_1G_Double)
{
  uint64_t seed = RNG::uint64_from_random_device();
  VCLRNG<TypeParam> core(seed, __PRETTY_FUNCTION__, "core");
  const unsigned N = unsigned(UINT64_C(64000000000)/TypeParam::vec_bits);
  typename TypeParam::double_vt sum(0);
  for(unsigned i=0;i<N;i++) {
    sum = core.uniform_double();
  }
  EXPECT_TRUE(horizontal_and(sum >= 0.0));
}

TYPED_TEST(VCLSpeedTestRNG, VEC_1G_Double_Alt)
{
  uint64_t seed = RNG::uint64_from_random_device();
  VCLRNG<TypeParam> core(seed, __PRETTY_FUNCTION__, "core");
  const unsigned N = unsigned(UINT64_C(64000000000)/TypeParam::vec_bits);
  typename TypeParam::double_vt sum(0);
  for(unsigned i=0;i<N;i++) {
    sum = core.uniform_double_alt();
  }
  EXPECT_TRUE(horizontal_and(sum >= 0.0));
}

TYPED_TEST(VCLSpeedTestRNG, VEC_1G_Double_53bit)
{
  uint64_t seed = RNG::uint64_from_random_device();
  VCLRNG<TypeParam> core(seed, __PRETTY_FUNCTION__, "core");
  const unsigned N = unsigned(UINT64_C(64000000000)/TypeParam::vec_bits);
  typename TypeParam::double_vt sum(0);
  for(unsigned i=0;i<N;i++) {
    sum = core.uniform_double_53bit();
  }
  EXPECT_TRUE(horizontal_and(sum >= 0.0));
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
