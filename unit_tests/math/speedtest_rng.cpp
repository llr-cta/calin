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

TEST(NR3SpeedTestRNG, Scalar_1G_Normal_BM)
{
  double sum = 0;
  uint64_t seed = RNG::uint64_from_random_device();
  RNG rng(seed);
  for(unsigned i=0;i<1000000000;i++)
    sum += rng.normal();
  EXPECT_NE(sum,0.0);
}

TEST(NR3SpeedTestRNG, Scalar_1G_Normal_Ziggurat)
{
  double sum = 0;
  uint64_t seed = RNG::uint64_from_random_device();
  RNG rng(seed);
  for(unsigned i=0;i<1000000000;i++)
    sum += rng.normal_ziggurat();
  EXPECT_NE(sum,0.0);
}

TEST(NR3SpeedTestRNG, Scalar_1G_Exponential_Log)
{
  double sum = 0;
  uint64_t seed = RNG::uint64_from_random_device();
  RNG rng(seed);
  for(unsigned i=0;i<1000000000;i++)
    sum += rng.exponential();
  EXPECT_GE(sum,0.0);
}

TEST(NR3SpeedTestRNG, Scalar_1G_Exponential_Ziggurat)
{
  double sum = 0;
  uint64_t seed = RNG::uint64_from_random_device();
  RNG rng(seed);
  for(unsigned i=0;i<1000000000;i++)
    sum += rng.exponential_ziggurat();
  EXPECT_GE(sum,0.0);
}

TEST(NR3SpeedTestRNG, Scalar_1G_XExpMinusXSquare)
{
  double sum = 0;
  uint64_t seed = RNG::uint64_from_random_device();
  RNG rng(seed);
  for(unsigned i=0;i<1000000000;i++)
    sum += std::sqrt(-std::log(rng.uniform_double()));
  EXPECT_GE(sum,0.0);
}

TEST(NR3SpeedTestRNG, Scalar_1G_XExpMinusXSquare_Ziggurat)
{
  double sum = 0;
  uint64_t seed = RNG::uint64_from_random_device();
  RNG rng(seed);
  for(unsigned i=0;i<1000000000;i++)
    sum += rng.x_exp_minus_x_squared_ziggurat();
  EXPECT_GE(sum,0.0);
}


template<typename VCLArchitecture> class alignas(VCLArchitecture::vec_bytes) VCLSpeedTestRNG :
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
