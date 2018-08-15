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

template<typename VCLArchitecture> class NR3_VCLCoreTests :
  public testing::Test
{
public:
};

using ArchTypes = ::testing::Types<VCL128Architecture, VCL256Architecture, VCL512Architecture>;

TYPED_TEST_CASE(NR3_VCLCoreTests, ArchTypes);

TYPED_TEST(NR3_VCLCoreTests, SameSeedsSameDeviates)
{
  uint64_t seed = RNG::uint64_from_random_device();
  NR3_VCLRNGCore<TypeParam> core(seed);
  NR3_VCLRNGCore<TypeParam> core2(seed);
  EXPECT_EQ(core.seed(), core2.seed());
  EXPECT_EQ(core.calls(), core2.calls());
  EXPECT_TRUE(horizontal_and(core.sequence_seeds() == core2.sequence_seeds()));
  for(unsigned i=0;i<1000;i++) {
    EXPECT_TRUE(horizontal_and(core.uniform_uint64() == core2.uniform_uint64()));
  }
}

TYPED_TEST(NR3_VCLCoreTests, DifferentSeedsDifferentDeviates)
{
  uint64_t seed = RNG::uint64_from_random_device();
  NR3_VCLRNGCore<TypeParam> core(seed);
  NR3_VCLRNGCore<TypeParam> core2(seed+1);
  EXPECT_NE(core.seed(), core2.seed());
  EXPECT_TRUE(horizontal_and(core.sequence_seeds() != core2.sequence_seeds()));
  for(unsigned i=0;i<1000;i++) {
    EXPECT_TRUE(horizontal_and(core.uniform_uint64() != core2.uniform_uint64()));
  }
}

TYPED_TEST(NR3_VCLCoreTests, SaveAndRestoreState)
{
  NR3_VCLRNGCore<TypeParam> core;
  for(unsigned i=0;i<1000;i++) {
    core.uniform_uint64();
  }
  auto* state = core.as_proto();
  NR3_VCLRNGCore<TypeParam> core2(NR3_VCLRNGCore<TypeParam>::core_data(*state), /* restore_state = */ true);
  EXPECT_EQ(core.seed(), core2.seed());
  EXPECT_EQ(core.calls(), core2.calls());
  EXPECT_TRUE(horizontal_and(core.sequence_seeds() == core2.sequence_seeds()));
  for(unsigned i=0;i<1000;i++) {
    EXPECT_TRUE(horizontal_and(core.uniform_uint64() == core2.uniform_uint64()));
  }
  delete state;
}

TYPED_TEST(NR3_VCLCoreTests, SaveAndRestoreState_Factory)
{
  NR3_VCLRNGCore<TypeParam> core;
  for(unsigned i=0;i<1000;i++) {
    core.uniform_uint64();
  }
  auto* state = core.as_proto();
  VCLRNGCore<TypeParam>* core2 =
    VCLRNGCore<TypeParam>::create_from_proto(*state, /* restore_state = */ true);
  for(unsigned i=0;i<1000;i++) {
    EXPECT_TRUE(horizontal_and(core.uniform_uint64() == core2->uniform_uint64()));
  }
  delete core2;
  delete state;
}

TYPED_TEST(NR3_VCLCoreTests, EqualsScalarNR3Implentation)
{
  NR3_VCLRNGCore<TypeParam> core;
  NR3RNGCore* scalar_cores[TypeParam::num_uint64];
  for(unsigned j=0;j<TypeParam::num_uint64;j++)
    scalar_cores[j] = new NR3RNGCore(core.sequence_seeds()[j]);
  for(unsigned i=0;i<1000000;i++) {
    typename TypeParam::uint64_vt x = core.uniform_uint64();
    for(unsigned j=0;j<TypeParam::num_uint64;j++)
      EXPECT_EQ(x[j], scalar_cores[j]->uniform_uint64());
  }
  for(unsigned j=0;j<TypeParam::num_uint64;j++)
    delete scalar_cores[j];
}

TYPED_TEST(NR3_VCLCoreTests, 64GBitSpeedTest)
{
  uint64_t seed = RNG::uint64_from_random_device();
  NR3_VCLRNGCore<TypeParam> core(seed);
  const unsigned N = unsigned(UINT64_C(64000000000)/TypeParam::vec_bits);
  typename TypeParam::uint64_vt x(0);
  for(unsigned i=0;i<N;i++) {
    x = core.uniform_uint64();
  }
  EXPECT_TRUE(horizontal_and(x >= UINT64_C(0)));
}

// TYPED_TEST(NR3_VCLCoreTests, 64GBitSpeedTest_Cmp)
// {
//   uint64_t seed = RNG::uint64_from_random_device();
//   NR3_VCLRNGCore<TypeParam> core(seed);
//   const unsigned N = unsigned(UINT64_C(64000000000)/TypeParam::vec_bits);
//   typename TypeParam::uint64_vt x(0);
//   for(unsigned i=0;i<N;i++) {
//     x = core.uniform_uint64_reference();
//   }
//   EXPECT_TRUE(horizontal_and(x >= UINT64_C(0)));
// }


TYPED_TEST(NR3_VCLCoreTests, 64GBitSpeedTest_Float)
{
  uint64_t seed = RNG::uint64_from_random_device();
  VCLRNG<TypeParam> core(seed);
  const unsigned N = unsigned(UINT64_C(64000000000)/TypeParam::vec_bits);
  typename TypeParam::float_vt x;
  for(unsigned i=0;i<N;i++) {
    x = core.uniform_float_zc();
  }
  EXPECT_TRUE(horizontal_and(x <= 0.5));
  EXPECT_TRUE(horizontal_and(x >= -0.5));
}

template<typename float_vt> void verify_float_range(const std::string& what,
  const float_vt& x, float xmin, float xmax)
{
  EXPECT_FALSE(horizontal_or(x>xmax))
    << "Out of range: " << what << " above maximum : " << x << " > " << xmax;
  EXPECT_FALSE(horizontal_or(x<xmin))
    << "Out of range: " << what << " below minimum : " << x << " < " << xmin;
}

TYPED_TEST(NR3_VCLCoreTests, UniformMoments)
{
  uint64_t seed = RNG::std_test_seed; //RNG::uint64_from_random_device();
  VCLRNG<TypeParam> core(seed);
  const unsigned N = 1000000;
  BasicKahanAccumulator<typename TypeParam::float_vt> sumx;
  BasicKahanAccumulator<typename TypeParam::float_vt> sumxx(0.0);
  BasicKahanAccumulator<typename TypeParam::float_vt> sumxxx(0.0);
  typename TypeParam::float_vt xmax(-1.0);
  typename TypeParam::float_vt xmin(1.0);
  typename TypeParam::float_vt x;
  for(unsigned i=0;i<N;i++) {
    x = core.uniform_float_zc();
    xmax = max(x, xmax);
    xmin = min(x, xmin);
    sumx.accumulate(x);
    sumxx.accumulate(x*x);
    sumxxx.accumulate(x*x*x);
  }

  verify_float_range("xmax", xmax, 0.5*0.9999, 0.5);
  verify_float_range("xmin", xmin, -0.5, -0.5*0.9999);

  typename TypeParam::float_vt m1 = sumx.total()/double(N);
  typename TypeParam::float_vt m2 = sumxx.total()/double(N);
  typename TypeParam::float_vt m3 = sumxxx.total()/double(N);

  verify_float_range("m1", m1, -0.001, 0.001);
  verify_float_range("m2", m2, 1/12.0-0.001, 1/12.0+0.001);
  verify_float_range("m3", m3, -0.0001, 0.0001);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
