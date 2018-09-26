/*

   calin/unit_tests/math/test_rng_vcl.cpp -- Stephen Fegan -- 2018-08-14

   Unit tests for VCL RNG

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

#include <math/accumulator.hpp>
#include <math/rng_vcl.hpp>
#include <util/vcl.hpp>

using namespace calin::math::accumulator;
using namespace calin::math::rng;
using namespace calin::util::vcl;

TEST(NR3SpeedTestRNG, Scalar_UInt64)
{
  uint64_t sum = 0ULL;
  uint64_t seed = RNG::uint64_from_random_device();
  NR3RNGCore core(seed);
  for(unsigned i=0;i<1000000000;i++)
    sum += core.uniform_uint64();
  EXPECT_GE(sum, 0ULL);
}

#ifdef CALIN_HAS_NR3_AVX2_RNGCORE
TEST(NR3SpeedTestRNG, AVX2_UInt64)
{
  uint64_t sum = 0ULL;
  uint64_t seed = RNG::uint64_from_random_device();
  NR3_AVX2_RNGCore core(seed);
  for(unsigned i=0;i<1000000000;i++)
    sum += core.uniform_uint64();
  EXPECT_GE(sum, 0ULL);
}

TEST(NR3SpeedTestRNG, AVX2_VEC_UInt64)
{
  __m256i sum = _mm256_setzero_si256();
  uint64_t seed = RNG::uint64_from_random_device();
  NR3_AVX2_RNGCore core(seed);
  for(unsigned i=0;i<250000000;i++)
    sum = _mm256_add_epi64(sum, core.uniform_m256i());
  EXPECT_GE(core.uniform_uint64(), 0ULL);
}
#endif // defined CALIN_HAS_NR3_AVX2_RNGCORE

TEST(NR3SpeedTestRNG, Scalar_Float)
{
  double sum = 0;
  uint64_t seed = RNG::uint64_from_random_device();
  NR3RNGCore core(seed);
  for(unsigned i=0;i<2000000000;i++)
    sum += 2.328306437e-10 * float(unsigned(core.uniform_uint64()));
  EXPECT_GE(sum,0.0);
}

#ifdef CALIN_HAS_NR3_AVX2_RNGCORE
TEST(NR3SpeedTestRNG, AVX2_VEC_Float)
{
  __m256 sum = _mm256_setzero_ps();
  uint64_t seed = RNG::uint64_from_random_device();
  NR3_AVX2_RNGCore core(seed);
  for(unsigned i=0;i<250000000;i++)
    sum = _mm256_add_ps(sum, core.uniform_m256());
  EXPECT_GE(reinterpret_cast<float*>(&sum)[0],0.0);
}
#endif // defined CALIN_HAS_NR3_AVX2_RNGCORE

template<typename VCLArchitecture> class VCLSpeedTestRNG :
  public testing::Test
{
public:
};

using ArchTypes = ::testing::Types<VCL128Architecture, VCL256Architecture, VCL512Architecture>;

TYPED_TEST_CASE(VCLSpeedTestRNG, ArchTypes);

TYPED_TEST(VCLSpeedTestRNG, VEC_Uint64)
{
  uint64_t seed = RNG::uint64_from_random_device();
  NR3_VCLRNGCore<TypeParam> core(seed, __PRETTY_FUNCTION__, "core");
  const unsigned N = unsigned(UINT64_C(64000000000)/TypeParam::vec_bits);
  typename TypeParam::uint64_vt x(0);
  for(unsigned i=0;i<N;i++) {
    x = core.uniform_uint64();
  }
  EXPECT_TRUE(horizontal_and(x >= UINT64_C(0)));
}

TYPED_TEST(VCLSpeedTestRNG, VEC_Float)
{
  uint64_t seed = RNG::uint64_from_random_device();
  VCLRNG<TypeParam> core(seed, __PRETTY_FUNCTION__, "core");
  const unsigned N = unsigned(UINT64_C(64000000000)/TypeParam::vec_bits);
  typename TypeParam::float_vt x;
  for(unsigned i=0;i<N;i++) {
    x = core.uniform_float_zc();
  }
  EXPECT_TRUE(horizontal_and(x <= 0.5));
  EXPECT_TRUE(horizontal_and(x >= -0.5));
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
