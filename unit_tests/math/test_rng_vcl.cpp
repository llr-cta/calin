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

template<typename VCLArchitecture> class alignas(VCLArchitecture::vec_bytes) NR3_VCLRNGCoreTests :
  public testing::Test
{
public:
};

using ArchTypes = ::testing::Types<VCL128Architecture, VCL256Architecture, VCL512Architecture>;

TYPED_TEST_CASE(NR3_VCLRNGCoreTests, ArchTypes);

TYPED_TEST(NR3_VCLRNGCoreTests, SameSeedsSameDeviates)
{
  uint64_t seed = RNG::uint64_from_random_device();
  NR3_VCLRNGCore<TypeParam> core(seed, __PRETTY_FUNCTION__, "core");
  NR3_VCLRNGCore<TypeParam> core2(seed, __PRETTY_FUNCTION__, "core2");
  EXPECT_EQ(core.seed(), core2.seed());
  EXPECT_EQ(core.calls(), core2.calls());
  EXPECT_TRUE(horizontal_and(core.sequence_seeds() == core2.sequence_seeds()));
  for(unsigned i=0;i<1000;i++) {
    EXPECT_TRUE(horizontal_and(core.uniform_uint64() == core2.uniform_uint64()));
  }
}

TYPED_TEST(NR3_VCLRNGCoreTests, DifferentSeedsDifferentDeviates)
{
  uint64_t seed = RNG::uint64_from_random_device();
  NR3_VCLRNGCore<TypeParam> core(seed, __PRETTY_FUNCTION__, "core");
  NR3_VCLRNGCore<TypeParam> core2(seed+1, __PRETTY_FUNCTION__, "core2");
  EXPECT_NE(core.seed(), core2.seed());
  EXPECT_TRUE(horizontal_and(core.sequence_seeds() != core2.sequence_seeds()));
  for(unsigned i=0;i<1000;i++) {
    EXPECT_TRUE(horizontal_and(core.uniform_uint64() != core2.uniform_uint64()));
  }
}

TYPED_TEST(NR3_VCLRNGCoreTests, SaveAndRestoreState)
{
  NR3_VCLRNGCore<TypeParam> core(__PRETTY_FUNCTION__, "core");
  for(unsigned i=0;i<1000;i++) {
    core.uniform_uint64();
  }
  auto* state = core.as_proto();
  NR3_VCLRNGCore<TypeParam> core2(NR3_VCLRNGCore<TypeParam>::core_data(*state),
    /* restore_state = */ true, __PRETTY_FUNCTION__, "core2");
  EXPECT_EQ(core.seed(), core2.seed());
  EXPECT_EQ(core.calls(), core2.calls());
  EXPECT_TRUE(horizontal_and(core.sequence_seeds() == core2.sequence_seeds()));
  for(unsigned i=0;i<1000;i++) {
    EXPECT_TRUE(horizontal_and(core.uniform_uint64() == core2.uniform_uint64()));
  }
  delete state;
}

TYPED_TEST(NR3_VCLRNGCoreTests, SaveAndRestoreState_Factory)
{
  NR3_VCLRNGCore<TypeParam> core(__PRETTY_FUNCTION__, "core");
  for(unsigned i=0;i<1000;i++) {
    core.uniform_uint64();
  }
  auto* state = core.as_proto();
  VCLRNGCore<TypeParam>* core2 =
    VCLRNGCore<TypeParam>::create_from_proto(*state, /* restore_state = */ true,
      __PRETTY_FUNCTION__, "core");
  for(unsigned i=0;i<1000;i++) {
    EXPECT_TRUE(horizontal_and(core.uniform_uint64() == core2->uniform_uint64()));
  }
  delete core2;
  delete state;
}

TYPED_TEST(NR3_VCLRNGCoreTests, EqualsScalarNR3Implentation)
{
  NR3_VCLRNGCore<TypeParam> core(__PRETTY_FUNCTION__, "core");
  NR3RNGCore* scalar_cores[TypeParam::num_uint64];
  for(unsigned j=0;j<TypeParam::num_uint64;j++)
    scalar_cores[j] = new NR3RNGCore(core.sequence_seeds()[j],
      __PRETTY_FUNCTION__, "scalar_cores["+std::to_string(j)+"]");
  for(unsigned i=0;i<1000000;i++) {
    typename TypeParam::uint64_vt x = core.uniform_uint64();
    for(unsigned j=0;j<TypeParam::num_uint64;j++)
      EXPECT_EQ(x[j], scalar_cores[j]->uniform_uint64());
  }
  for(unsigned j=0;j<TypeParam::num_uint64;j++)
    delete scalar_cores[j];
}

template<typename VCLArchitecture> class alignas(VCLArchitecture::vec_bytes) VCLRNGTests :
  public testing::Test
{
public:
};

using ArchTypes = ::testing::Types<VCL128Architecture, VCL256Architecture, VCL512Architecture>;

TYPED_TEST_CASE(VCLRNGTests, ArchTypes);

template<typename float_vt> void verify_float_range(const std::string& what,
  const float_vt& x, float xmin, float xmax)
{
  EXPECT_FALSE(horizontal_or(x>xmax))
    << "Out of range: " << what << " above maximum : " << x << " > " << xmax;
  EXPECT_FALSE(horizontal_or(x<xmin))
    << "Out of range: " << what << " below minimum : " << x << " < " << xmin;
}

template<typename double_vt> void verify_double_range(const std::string& what,
  const double_vt& x, double xmin, double xmax)
{
  EXPECT_FALSE(horizontal_or(x>xmax))
    << "Out of range: " << what << " above maximum : " << x << " > " << xmax;
  EXPECT_FALSE(horizontal_or(x<xmin))
    << "Out of range: " << what << " below minimum : " << x << " < " << xmin;
}

TYPED_TEST(VCLRNGTests, UniformFloatZCMoments)
{
  uint64_t seed = RNG::std_test_seed; //RNG::uint64_from_random_device();
  VCLRNG<TypeParam> core(seed, __PRETTY_FUNCTION__, "core");
  const unsigned N = 1000000;
  BasicKahanAccumulator<typename TypeParam::float_vt> sumx;
  BasicKahanAccumulator<typename TypeParam::float_vt> sumxx;
  BasicKahanAccumulator<typename TypeParam::float_vt> sumxxx;
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

  typename TypeParam::float_vt m1 = sumx.total()/float(N);
  typename TypeParam::float_vt m2 = sumxx.total()/float(N);
  typename TypeParam::float_vt m3 = sumxxx.total()/float(N);

  verify_float_range("m1", m1, -0.001, 0.001);
  verify_float_range("m2", m2, 1/12.0-0.001, 1/12.0+0.001);
  verify_float_range("m3", m3, -0.0001, 0.0001);
}

TYPED_TEST(VCLRNGTests, UniformDoubleZCMoments)
{
  uint64_t seed = RNG::std_test_seed; //RNG::uint64_from_random_device();
  VCLRNG<TypeParam> core(seed, __PRETTY_FUNCTION__, "core");
  const unsigned N = 1000000;
  BasicKahanAccumulator<typename TypeParam::double_vt> sumx;
  BasicKahanAccumulator<typename TypeParam::double_vt> sumxx;
  BasicKahanAccumulator<typename TypeParam::double_vt> sumxxx;
  typename TypeParam::double_vt xmax(-1.0);
  typename TypeParam::double_vt xmin(1.0);
  typename TypeParam::double_vt x;
  for(unsigned i=0;i<N;i++) {
    x = core.uniform_double_zc();
    xmax = max(x, xmax);
    xmin = min(x, xmin);
    sumx.accumulate(x);
    sumxx.accumulate(x*x);
    sumxxx.accumulate(x*x*x);
  }

  verify_double_range("xmax", xmax, 0.5*0.9999, 0.5);
  verify_double_range("xmin", xmin, -0.5, -0.5*0.9999);

  typename TypeParam::double_vt m1 = sumx.total()/double(N);
  typename TypeParam::double_vt m2 = sumxx.total()/double(N);
  typename TypeParam::double_vt m3 = sumxxx.total()/double(N);

  verify_double_range("m1", m1, -0.001, 0.001);
  verify_double_range("m2", m2, 1/12.0-0.001, 1/12.0+0.001);
  verify_double_range("m3", m3, -0.0001, 0.0001);
}

TYPED_TEST(VCLRNGTests, ExponentialFloatMoments)
{
  uint64_t seed = RNG::std_test_seed; //RNG::uint64_from_random_device();
  VCLRNG<TypeParam> core(seed, __PRETTY_FUNCTION__, "core");
  const unsigned N = 1000000;
  BasicKahanAccumulator<typename TypeParam::float_vt> sumx;
  BasicKahanAccumulator<typename TypeParam::float_vt> sumxx;
  BasicKahanAccumulator<typename TypeParam::float_vt> sumxxx;
  typename TypeParam::float_vt x;
  for(unsigned i=0;i<N;i++) {
    x = core.exponential_float();
    sumx.accumulate(x);
    sumxx.accumulate(x*x);
    sumxxx.accumulate(x*x*x);
  }

  typename TypeParam::float_vt m1 = sumx.total()/float(N);
  typename TypeParam::float_vt m2 = sumxx.total()/float(N);
  typename TypeParam::float_vt m3 = sumxxx.total()/float(N);

  verify_float_range("m1", m1, 1.0-0.003, 1.0+0.003);
  verify_float_range("m2", m2, 2.0-0.01,  2.0+0.01);
  verify_float_range("m3", m3, 6.0-0.1,   6.0+0.1);
}

TYPED_TEST(VCLRNGTests, ExponentialDoubleMoments)
{
  uint64_t seed = RNG::std_test_seed; //RNG::uint64_from_random_device();
  VCLRNG<TypeParam> core(seed, __PRETTY_FUNCTION__, "core");
  const unsigned N = 1000000;
  BasicKahanAccumulator<typename TypeParam::double_vt> sumx;
  BasicKahanAccumulator<typename TypeParam::double_vt> sumxx;
  BasicKahanAccumulator<typename TypeParam::double_vt> sumxxx;
  typename TypeParam::double_vt x;
  for(unsigned i=0;i<N;i++) {
    x = core.exponential_double();
    sumx.accumulate(x);
    sumxx.accumulate(x*x);
    sumxxx.accumulate(x*x*x);
  }

  typename TypeParam::double_vt m1 = sumx.total()/double(N);
  typename TypeParam::double_vt m2 = sumxx.total()/double(N);
  typename TypeParam::double_vt m3 = sumxxx.total()/double(N);

  verify_double_range("m1", m1, 1.0-0.003, 1.0+0.003);
  verify_double_range("m2", m2, 2.0-0.012, 2.0+0.012);
  verify_double_range("m3", m3, 6.0-0.1,   6.0+0.1);
}

TYPED_TEST(VCLRNGTests, SinCosFloatMoments)
{
  uint64_t seed = RNG::std_test_seed; //RNG::uint64_from_random_device();
  VCLRNG<TypeParam> core(seed, __PRETTY_FUNCTION__, "core");
  const unsigned N = 1000000;
  BasicKahanAccumulator<typename TypeParam::float_vt> sums;
  BasicKahanAccumulator<typename TypeParam::float_vt> sumss;
  BasicKahanAccumulator<typename TypeParam::float_vt> sumsss;
  BasicKahanAccumulator<typename TypeParam::float_vt> sumc;
  BasicKahanAccumulator<typename TypeParam::float_vt> sumcc;
  BasicKahanAccumulator<typename TypeParam::float_vt> sumccc;
  BasicKahanAccumulator<typename TypeParam::float_vt> sumsc;
  typename TypeParam::float_vt s;
  typename TypeParam::float_vt c;
  for(unsigned i=0;i<N;i++) {
    core.sincos_float(s,c);
    sums.accumulate(s);
    sumss.accumulate(s*s);
    sumsss.accumulate(s*s*s);
    sumc.accumulate(c);
    sumcc.accumulate(c*c);
    sumccc.accumulate(c*c*c);
    sumsc.accumulate(s*c);
  }

  typename TypeParam::float_vt m1s = sums.total()/float(N);
  typename TypeParam::float_vt m2s = sumss.total()/float(N);
  typename TypeParam::float_vt m3s = sumsss.total()/float(N);
  typename TypeParam::float_vt m1c = sumc.total()/float(N);
  typename TypeParam::float_vt m2c = sumcc.total()/float(N);
  typename TypeParam::float_vt m3c = sumccc.total()/float(N);
  typename TypeParam::float_vt m2sc = sumsc.total()/float(N);

  verify_float_range("m1s", m1s, -0.005, 0.005);
  verify_float_range("m2s", m2s, 0.5-0.005,  0.5+0.005);
  verify_float_range("m3s", m3s, -0.005, 0.005);
  verify_float_range("m1c", m1c, -0.005, 0.005);
  verify_float_range("m2c", m2c, 0.5-0.005,  0.5+0.005);
  verify_float_range("m3c", m3c, -0.005, 0.005);
  verify_float_range("m2sc", m2sc, -0.002, 0.002);
}

TYPED_TEST(VCLRNGTests, SinCosDoubleMoments)
{
  uint64_t seed = RNG::std_test_seed; //RNG::uint64_from_random_device();
  VCLRNG<TypeParam> core(seed, __PRETTY_FUNCTION__, "core");
  const unsigned N = 1000000;
  BasicKahanAccumulator<typename TypeParam::double_vt> sums;
  BasicKahanAccumulator<typename TypeParam::double_vt> sumss;
  BasicKahanAccumulator<typename TypeParam::double_vt> sumsss;
  BasicKahanAccumulator<typename TypeParam::double_vt> sumc;
  BasicKahanAccumulator<typename TypeParam::double_vt> sumcc;
  BasicKahanAccumulator<typename TypeParam::double_vt> sumccc;
  BasicKahanAccumulator<typename TypeParam::double_vt> sumsc;
  typename TypeParam::double_vt s;
  typename TypeParam::double_vt c;
  for(unsigned i=0;i<N;i++) {
    core.sincos_double(s,c);
    sums.accumulate(s);
    sumss.accumulate(s*s);
    sumsss.accumulate(s*s*s);
    sumc.accumulate(c);
    sumcc.accumulate(c*c);
    sumccc.accumulate(c*c*c);
    sumsc.accumulate(s*c);
  }

  typename TypeParam::double_vt m1s = sums.total()/double(N);
  typename TypeParam::double_vt m2s = sumss.total()/double(N);
  typename TypeParam::double_vt m3s = sumsss.total()/double(N);
  typename TypeParam::double_vt m1c = sumc.total()/double(N);
  typename TypeParam::double_vt m2c = sumcc.total()/double(N);
  typename TypeParam::double_vt m3c = sumccc.total()/double(N);
  typename TypeParam::double_vt m2sc = sumsc.total()/double(N);

  verify_double_range("m1s", m1s, -0.005, 0.005);
  verify_double_range("m2s", m2s, 0.5-0.005,  0.5+0.005);
  verify_double_range("m3s", m3s, -0.005, 0.005);
  verify_double_range("m1c", m1c, -0.005, 0.005);
  verify_double_range("m2c", m2c, 0.5-0.005,  0.5+0.005);
  verify_double_range("m3c", m3c, -0.005, 0.005);
  verify_double_range("m2sc", m2sc, -0.002, 0.002);
}

TYPED_TEST(VCLRNGTests, NormalFloatBMMoments)
{
  uint64_t seed = RNG::std_test_seed; //RNG::uint64_from_random_device();
  VCLRNG<TypeParam> core(seed, __PRETTY_FUNCTION__, "core");
  const unsigned N = 1000000;
  BasicKahanAccumulator<typename TypeParam::float_vt> sumx;
  BasicKahanAccumulator<typename TypeParam::float_vt> sumxx;
  BasicKahanAccumulator<typename TypeParam::float_vt> sumxxx;
  BasicKahanAccumulator<typename TypeParam::float_vt> sumy;
  BasicKahanAccumulator<typename TypeParam::float_vt> sumyy;
  BasicKahanAccumulator<typename TypeParam::float_vt> sumyyy;
  BasicKahanAccumulator<typename TypeParam::float_vt> sumxy;
  typename TypeParam::float_vt x;
  typename TypeParam::float_vt y;
  for(unsigned i=0;i<N;i++) {
    core.normal_pair_float_bm(x,y);
    sumx.accumulate(x);
    sumxx.accumulate(x*x);
    sumxxx.accumulate(x*x*x);
    sumy.accumulate(y);
    sumyy.accumulate(y*y);
    sumyyy.accumulate(y*y*y);
    sumxy.accumulate(x*y);
  }

  typename TypeParam::float_vt m1x = sumx.total()/float(N);
  typename TypeParam::float_vt m2x = sumxx.total()/float(N);
  typename TypeParam::float_vt m3x = sumxxx.total()/float(N);
  typename TypeParam::float_vt m1y = sumy.total()/float(N);
  typename TypeParam::float_vt m2y = sumyy.total()/float(N);
  typename TypeParam::float_vt m3y = sumyyy.total()/float(N);
  typename TypeParam::float_vt m2xy = sumxy.total()/float(N);

  verify_float_range("m1x", m1x, -0.005, 0.005);
  verify_float_range("m2x", m2x, 1.0-0.005,  1.0+0.005);
  verify_float_range("m3x", m3x, -0.02, 0.02);
  verify_float_range("m1y", m1y, -0.005, 0.005);
  verify_float_range("m2y", m2y, 1.0-0.005,  1.0+0.005);
  verify_float_range("m3y", m3y, -0.02, 0.02);
  verify_float_range("m2xy", m2xy, -0.005, 0.005);
}

TYPED_TEST(VCLRNGTests, NormalDoubleBMMoments)
{
  uint64_t seed = RNG::std_test_seed; //RNG::uint64_from_random_device();
  VCLRNG<TypeParam> core(seed, __PRETTY_FUNCTION__, "core");
  const unsigned N = 1000000;
  BasicKahanAccumulator<typename TypeParam::double_vt> sumx;
  BasicKahanAccumulator<typename TypeParam::double_vt> sumxx;
  BasicKahanAccumulator<typename TypeParam::double_vt> sumxxx;
  BasicKahanAccumulator<typename TypeParam::double_vt> sumy;
  BasicKahanAccumulator<typename TypeParam::double_vt> sumyy;
  BasicKahanAccumulator<typename TypeParam::double_vt> sumyyy;
  BasicKahanAccumulator<typename TypeParam::double_vt> sumxy;
  typename TypeParam::double_vt x;
  typename TypeParam::double_vt y;
  for(unsigned i=0;i<N;i++) {
    core.normal_pair_double_bm(x,y);
    sumx.accumulate(x);
    sumxx.accumulate(x*x);
    sumxxx.accumulate(x*x*x);
    sumy.accumulate(y);
    sumyy.accumulate(y*y);
    sumyyy.accumulate(y*y*y);
    sumxy.accumulate(x*y);
  }

  typename TypeParam::double_vt m1x = sumx.total()/double(N);
  typename TypeParam::double_vt m2x = sumxx.total()/double(N);
  typename TypeParam::double_vt m3x = sumxxx.total()/double(N);
  typename TypeParam::double_vt m1y = sumy.total()/double(N);
  typename TypeParam::double_vt m2y = sumyy.total()/double(N);
  typename TypeParam::double_vt m3y = sumyyy.total()/double(N);
  typename TypeParam::double_vt m2xy = sumxy.total()/double(N);

  verify_double_range("m1x", m1x, -0.005, 0.005);
  verify_double_range("m2x", m2x, 1.0-0.005,  1.0+0.005);
  verify_double_range("m3x", m3x, -0.02, 0.02);
  verify_double_range("m1y", m1y, -0.005, 0.005);
  verify_double_range("m2y", m2y, 1.0-0.005,  1.0+0.005);
  verify_double_range("m3y", m3y, -0.02, 0.02);
  verify_double_range("m2xy", m2xy, -0.005, 0.005);
}

TYPED_TEST(VCLRNGTests, NormalFloatZigguratMoments)
{
  uint64_t seed = RNG::std_test_seed; //RNG::uint64_from_random_device();
  VCLRNG<TypeParam> core(seed, __PRETTY_FUNCTION__, "core");
  const unsigned N = 1000000;
  BasicKahanAccumulator<typename TypeParam::float_vt> sumx;
  BasicKahanAccumulator<typename TypeParam::float_vt> sumxx;
  BasicKahanAccumulator<typename TypeParam::float_vt> sumxxx;
  typename TypeParam::float_vt x;
  for(unsigned i=0;i<N;i++) {
    x = core.normal_float_ziggurat();
    sumx.accumulate(x);
    sumxx.accumulate(x*x);
    sumxxx.accumulate(x*x*x);
  }

  typename TypeParam::float_vt m1x = sumx.total()/float(N);
  typename TypeParam::float_vt m2x = sumxx.total()/float(N);
  typename TypeParam::float_vt m3x = sumxxx.total()/float(N);

  verify_float_range("m1x", m1x, -0.005, 0.005);
  verify_float_range("m2x", m2x, 1.0-0.005,  1.0+0.005);
  verify_float_range("m3x", m3x, -0.02, 0.02);
}

TYPED_TEST(VCLRNGTests, NormalDoubleZigguratMoments)
{
  uint64_t seed = RNG::std_test_seed; //RNG::uint64_from_random_device();
  VCLRNG<TypeParam> core(seed, __PRETTY_FUNCTION__, "core");
  const unsigned N = 1000000;
  BasicKahanAccumulator<typename TypeParam::double_vt> sumx;
  BasicKahanAccumulator<typename TypeParam::double_vt> sumxx;
  BasicKahanAccumulator<typename TypeParam::double_vt> sumxxx;
  typename TypeParam::double_vt x;
  for(unsigned i=0;i<N;i++) {
    x = core.normal_double_ziggurat();
    sumx.accumulate(x);
    sumxx.accumulate(x*x);
    sumxxx.accumulate(x*x*x);
  }

  typename TypeParam::double_vt m1x = sumx.total()/double(N);
  typename TypeParam::double_vt m2x = sumxx.total()/double(N);
  typename TypeParam::double_vt m3x = sumxxx.total()/double(N);

  verify_double_range("m1x", m1x, -0.005, 0.005);
  verify_double_range("m2x", m2x, 1.0-0.005,  1.0+0.005);
  verify_double_range("m3x", m3x, -0.02, 0.02);
}

TYPED_TEST(VCLRNGTests, CDFUniformFloatZCMoments)
{
  uint64_t seed = RNG::std_test_seed; //RNG::uint64_from_random_device();
  VCLRNG<TypeParam> core(seed, __PRETTY_FUNCTION__, "core");
  std::vector<std::pair<double, double>> cdf;
  cdf.emplace_back(-0.5,0.0);
  cdf.emplace_back(0.5,1.0);
  std::vector<float> inverse_cdf = core.generate_inverse_cdf_float(cdf, 0);

  const unsigned N = 1000000;
  BasicKahanAccumulator<typename TypeParam::float_vt> sumx;
  BasicKahanAccumulator<typename TypeParam::float_vt> sumxx;
  BasicKahanAccumulator<typename TypeParam::float_vt> sumxxx;
  typename TypeParam::float_vt xmax(-1.0);
  typename TypeParam::float_vt xmin(1.0);
  typename TypeParam::float_vt x;
  for(unsigned i=0;i<N;i++) {
    x = core.from_inverse_cdf_float(inverse_cdf);
    xmax = max(x, xmax);
    xmin = min(x, xmin);
    sumx.accumulate(x);
    sumxx.accumulate(x*x);
    sumxxx.accumulate(x*x*x);
  }

  verify_float_range("xmax", xmax, 0.5*0.9999, 0.5);
  verify_float_range("xmin", xmin, -0.5, -0.5*0.9999);

  typename TypeParam::float_vt m1 = sumx.total()/float(N);
  typename TypeParam::float_vt m2 = sumxx.total()/float(N);
  typename TypeParam::float_vt m3 = sumxxx.total()/float(N);

  verify_float_range("m1", m1, -0.001, 0.001);
  verify_float_range("m2", m2, 1/12.0-0.001, 1/12.0+0.001);
  verify_float_range("m3", m3, -0.0001, 0.0001);
}

TYPED_TEST(VCLRNGTests, CDFUniformDoubleZCMoments)
{
  uint64_t seed = RNG::std_test_seed; //RNG::uint64_from_random_device();
  VCLRNG<TypeParam> core(seed, __PRETTY_FUNCTION__, "core");
  std::vector<std::pair<double, double>> cdf;
  cdf.emplace_back(-0.5,0.0);
  cdf.emplace_back(0.5,1.0);
  std::vector<double> inverse_cdf = core.generate_inverse_cdf_double(cdf, 0);

  const unsigned N = 1000000;
  BasicKahanAccumulator<typename TypeParam::double_vt> sumx;
  BasicKahanAccumulator<typename TypeParam::double_vt> sumxx;
  BasicKahanAccumulator<typename TypeParam::double_vt> sumxxx;
  typename TypeParam::double_vt xmax(-1.0);
  typename TypeParam::double_vt xmin(1.0);
  typename TypeParam::double_vt x;
  for(unsigned i=0;i<N;i++) {
    x = core.from_inverse_cdf_double(inverse_cdf);
    xmax = max(x, xmax);
    xmin = min(x, xmin);
    sumx.accumulate(x);
    sumxx.accumulate(x*x);
    sumxxx.accumulate(x*x*x);
  }

  verify_double_range("xmax", xmax, 0.5*0.9999, 0.5);
  verify_double_range("xmin", xmin, -0.5, -0.5*0.9999);

  typename TypeParam::double_vt m1 = sumx.total()/double(N);
  typename TypeParam::double_vt m2 = sumxx.total()/double(N);
  typename TypeParam::double_vt m3 = sumxxx.total()/double(N);

  verify_double_range("m1", m1, -0.001, 0.001);
  verify_double_range("m2", m2, 1/12.0-0.001, 1/12.0+0.001);
  verify_double_range("m3", m3, -0.0001, 0.0001);
}

TYPED_TEST(VCLRNGTests, CDFNormalFloatZCMoments)
{
  uint64_t seed = RNG::std_test_seed; //RNG::uint64_from_random_device();
  VCLRNG<TypeParam> core(seed, __PRETTY_FUNCTION__, "core");
  std::vector<std::pair<double, double>> cdf;
  cdf.emplace_back(-20.0, 0.0);
  for(double x=-10.0; x<=10.0; x+=0.1) {
    cdf.emplace_back(x, 0.5*(std::erf(x)+1.0));
  }
  cdf.emplace_back(20.0, 1.0);
  std::vector<float> inverse_cdf = core.generate_inverse_cdf_float(cdf, 4096);

  //for(auto x: inverse_cdf) { std::cout << x << '\n'; }

  const unsigned N = 1000000;
  BasicKahanAccumulator<typename TypeParam::float_vt> sumx;
  BasicKahanAccumulator<typename TypeParam::float_vt> sumxx;
  BasicKahanAccumulator<typename TypeParam::float_vt> sumxxx;
  typename TypeParam::float_vt x;
  for(unsigned i=0;i<N;i++) {
    x = core.from_inverse_cdf_float(inverse_cdf);
    sumx.accumulate(x);
    sumxx.accumulate(x*x);
    sumxxx.accumulate(x*x*x);
  }

  typename TypeParam::float_vt m1 = sumx.total()/float(N);
  typename TypeParam::float_vt m2 = sumxx.total()/float(N);
  typename TypeParam::float_vt m3 = sumxxx.total()/float(N);

  verify_float_range("m1", m1, -0.005, 0.005);
  verify_float_range("m2", m2, 0.5-0.01,  0.5+0.01);
  verify_float_range("m3", m3, -0.02, 0.02);
}

TYPED_TEST(VCLRNGTests, PoissonInt64Moments_ExpNegLambda_0_1)
{
  uint64_t seed = RNG::std_test_seed;
  VCLRNG<TypeParam> core(seed, __PRETTY_FUNCTION__, "core");
  const unsigned N = 1000000;

  BasicKahanAccumulator<typename TypeParam::double_vt> sumx;
  BasicKahanAccumulator<typename TypeParam::double_vt> sumxx;
  BasicKahanAccumulator<typename TypeParam::double_vt> sumxxx;
  typename TypeParam::double_vt exp_neg_lambda = std::exp(-0.1);

  for(unsigned i = 0; i < N; ++i) {
    typename TypeParam::double_vt x = to_double(core.poisson_small_exp_neg_lambda_int64(exp_neg_lambda));
    sumx.accumulate(x);
    sumxx.accumulate(x*x);
    sumxxx.accumulate(x*x*x);
  }

  typename TypeParam::double_vt m1 = sumx.total()/double(N);
  typename TypeParam::double_vt m2 = sumxx.total()/double(N);
  typename TypeParam::double_vt m3 = sumxxx.total()/double(N);

  double m1_exp = 0.1;
  double m1_tol = 0.003;
  double m2_exp = 0.11;
  double m2_tol = 0.003;
  double m3_exp = 0.131;
  double m3_tol = 0.003;

  std::string tag = "lambda=0.1 (exp neg lambda) ";
  verify_double_range(tag + "m1", m1, m1_exp - m1_tol, m1_exp + m1_tol);
  verify_double_range(tag + "m2", m2, m2_exp - m2_tol, m2_exp + m2_tol);
  verify_double_range(tag + "m3", m3, m3_exp - m3_tol, m3_exp + m3_tol);
}

template<typename TypeParam>
void test_poisson_int64_moments(double lambda_val,
  double m1_exp, double m1_tol,
  double m2_exp, double m2_tol,
  double m3_exp, double m3_tol,
  double small_lambda_max = 10.0)
{
  uint64_t seed = RNG::std_test_seed;
  VCLRNG<TypeParam> core(seed, __PRETTY_FUNCTION__, "core");
  const unsigned N = 1000000;

  BasicKahanAccumulator<typename TypeParam::double_vt> sumx;
  BasicKahanAccumulator<typename TypeParam::double_vt> sumxx;
  BasicKahanAccumulator<typename TypeParam::double_vt> sumxxx;
  typename TypeParam::double_vt lambda(lambda_val);

  for(unsigned i = 0; i < N; ++i) {
    typename TypeParam::double_vt x = to_double(core.poisson_int64(lambda, small_lambda_max));
    sumx.accumulate(x);
    sumxx.accumulate(x*x);
    sumxxx.accumulate(x*x*x);
  }

  typename TypeParam::double_vt m1 = sumx.total()/double(N);
  typename TypeParam::double_vt m2 = sumxx.total()/double(N);
  typename TypeParam::double_vt m3 = sumxxx.total()/double(N);

  std::string tag = "lambda=" + std::to_string(lambda_val) + " (small_max=" + std::to_string(small_lambda_max) + ") ";
  verify_double_range(tag + "m1", m1, m1_exp - m1_tol, m1_exp + m1_tol);
  verify_double_range(tag + "m2", m2, m2_exp - m2_tol, m2_exp + m2_tol);
  verify_double_range(tag + "m3", m3, m3_exp - m3_tol, m3_exp + m3_tol);
}

TYPED_TEST(VCLRNGTests, PoissonInt64Moments_Lambda0_1)
{
  test_poisson_int64_moments<TypeParam>(0.1, 0.1, 0.003, 0.11, 0.003, 0.131, 0.003);
}

TYPED_TEST(VCLRNGTests, PoissonInt64Moments_Lambda1)
{
  test_poisson_int64_moments<TypeParam>(1.0, 1.0, 0.008, 2.0, 0.03, 5.0, 0.10);
}

TYPED_TEST(VCLRNGTests, PoissonInt64Moments_Lambda10)
{
  test_poisson_int64_moments<TypeParam>(10.0, 10.0, 0.03, 110.0, 0.6, 1310.0, 6.0);
}

TYPED_TEST(VCLRNGTests, PoissonInt64Moments_Lambda10_SmallAlgo)
{
  test_poisson_int64_moments<TypeParam>(10.0, 10.0, 0.03, 110.0, 0.6, 1310.0, 6.0, /* small_lambda_max = */ 10.5);
}

TYPED_TEST(VCLRNGTests, PoissonInt64Moments_Lambda100)
{
  test_poisson_int64_moments<TypeParam>(100.0, 100.0, 0.1, 10100.0, 20.0, 1030100.0, 2500.0);
}

TYPED_TEST(VCLRNGTests, BorelTannerMoments)
{
  uint64_t seed = RNG::std_test_seed; //RNG::uint64_from_random_device();
  VCLRNG<TypeParam> core(seed, __PRETTY_FUNCTION__, "core");

  typename TypeParam::int64_vt k = 1;
  typename TypeParam::double_vt mu = 0.1;

  const unsigned N = 1000000;
  BasicKahanAccumulator<typename TypeParam::double_vt> sumx;
  BasicKahanAccumulator<typename TypeParam::double_vt> sumxx;
  BasicKahanAccumulator<typename TypeParam::double_vt> sumxxx;
  typename TypeParam::double_vt x;
  for(unsigned i=0;i<N;i++) {
    x = to_double(core.borel_tanner_int64(k, mu));
    sumx.accumulate(x);
    sumxx.accumulate(x*x);
    sumxxx.accumulate(x*x*x);
  }

  typename TypeParam::double_vt m1 = sumx.total()/double(N);
  typename TypeParam::double_vt m2 = sumxx.total()/double(N);
  typename TypeParam::double_vt m3 = sumxxx.total()/double(N);

  verify_double_range("m1", m1, 1.1100, 1.1120); // Expected mean = 1/(1-0.1) = 1.111111...
  verify_double_range("m2", m2, 1.3617, 1.3817); // Expected second moment = (1+0.1)/(1-0.1)^3 = 1.371742...
  verify_double_range("m3", m3, 2.0121, 2.0521); // Expected third moment = (1+4*0.1+0.1^2)/(1-0.1)^5 = 2.032098...
}

TYPED_TEST(VCLRNGTests, BorelTannerK1Moments)
{
  uint64_t seed = RNG::std_test_seed; //RNG::uint64_from_random_device();
  VCLRNG<TypeParam> core(seed, __PRETTY_FUNCTION__, "core");

  typename TypeParam::double_vt exp_neg_mu = std::exp(-0.1);

  const unsigned N = 1000000;
  BasicKahanAccumulator<typename TypeParam::double_vt> sumx;
  BasicKahanAccumulator<typename TypeParam::double_vt> sumxx;
  BasicKahanAccumulator<typename TypeParam::double_vt> sumxxx;
  typename TypeParam::double_vt x;
  for(unsigned i=0;i<N;i++) {
    x = to_double(core.borel_tanner_k1_exp_neg_mu_int64(exp_neg_mu));
    sumx.accumulate(x);
    sumxx.accumulate(x*x);
    sumxxx.accumulate(x*x*x);
  }

  typename TypeParam::double_vt m1 = sumx.total()/double(N);
  typename TypeParam::double_vt m2 = sumxx.total()/double(N);
  typename TypeParam::double_vt m3 = sumxxx.total()/double(N);

  verify_double_range("m1", m1, 1.1100, 1.1120); // Expected mean = 1/(1-0.1) = 1.111111...
  verify_double_range("m2", m2, 1.3617, 1.3817); // Expected second moment = (1+0.1)/(1-0.1)^3 = 1.371742...
  verify_double_range("m3", m3, 2.0121, 2.0521); // Expected third moment = (1+4*0.1+0.1^2)/(1-0.1)^5 = 2.032098...
}


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
