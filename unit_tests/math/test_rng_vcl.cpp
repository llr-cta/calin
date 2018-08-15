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

#include <math/rng_vcl.hpp>
#include <util/vcl.hpp>

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
  typename TypeParam::float_vt sum(0);
  for(unsigned i=0;i<N;i++) {
    sum += core.uniform_float();
  }
  EXPECT_TRUE(horizontal_and(sum >= 0.0));
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
