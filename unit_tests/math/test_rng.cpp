/*

   calin/unit_tests/math/test_rng.cpp -- Stephen Fegan -- 2015-11-22

   Unit tests for RNG

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, Ecole polytechnique, CNRS/IN2P3, Universite Paris-Saclay

   This file is part of "calin"

   "calin" is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 or
   later, as published by the Free Software Foundation.

   "calin" is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

*/

// TODO : Come up with more systematic approach for bounds on moments. For M1
// and M2 it shouldn't be a problem (see for example http://goo.gl/qP91oD).
// M3 is probably more tricky - I assume this would require calculating the
// 6th order moments.

#include <iostream>
#include <iomanip>
#include <gtest/gtest.h>
#include <vector>
#include <tuple>

#include <google/protobuf/text_format.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>

#include <math/rng.hpp>

using namespace calin::math::rng;

TEST(TestRNG, RandomDeviceFillsAllBits64) {
  std::vector<unsigned> count(64,0);
  for(unsigned i=0;i<100;i++)
  {
    uint64_t x = RNG::uint64_from_random_device();
    for(unsigned j=0;j<64;j++) {
      if(x&1)count[j]++;
      x>>=1; }
  }
  for(unsigned j=0;j<64;j++) {
    EXPECT_GE(count[j], 20U); EXPECT_LE(count[j], 80U); }
}

TEST(TestRNG, RandomDeviceFillsAllBits32) {
  std::vector<unsigned> count(32,0);
  for(unsigned i=0;i<100;i++)
  {
    uint32_t x = RNG::uint32_from_random_device();
    for(unsigned j=0;j<32;j++) {
      if(x&1)count[j]++;
      x>>=1; }
  }
  for(unsigned j=0;j<32;j++) {
    EXPECT_GE(count[j], 20U); EXPECT_LE(count[j], 80U); }
}

template<typename CORE> class CoreTests : public testing::Test {
 public:
};

TYPED_TEST_CASE_P(CoreTests);

TYPED_TEST_P(CoreTests, FillsAllBits64)
{
  TypeParam core(RNG::uint64_from_random_device());
  std::vector<unsigned> count(64,0);
  for(unsigned i=0;i<1000000;i++)
  {
    uint64_t x = core.uniform_uint64();
    for(unsigned j=0;j<64;j++) {
      if(x&1)count[j]++;
      x>>=1; }
  }
  for(unsigned j=0;j<64;j++) {
    EXPECT_GE(count[j], 490000U)
        << "uint64_t bit#" << j << " set too few times";
    EXPECT_LE(count[j], 510000U)
        << "uint64_t bit#" << j << " set too many times";
  }
}

TYPED_TEST_P(CoreTests, FillsAllBits32)
{
  TypeParam core(RNG::uint64_from_random_device());
  RNG rng(&core, false);
  std::vector<unsigned> count(32,0);
  for(unsigned i=0;i<1000000;i++)
  {
    uint32_t x = rng.uniform_uint32();
    for(unsigned j=0;j<32;j++) {
      if(x&1)count[j]++;
      x>>=1; }
  }
  for(unsigned j=0;j<32;j++) {
    EXPECT_GE(count[j], 490000U)
        << "uint32_t bit#" << j << " set too few times";
    EXPECT_LE(count[j], 510000U)
        << "uint32_t bit#" << j << " set too many times";
  }
}

TYPED_TEST_P(CoreTests, DifferentSeeds)
{
  TypeParam core1(1ULL);
  TypeParam core2(2ULL);
  unsigned nsame = 0;
  for(unsigned i=0;i<10000;i++)
    if(core1.uniform_uint64() == core2.uniform_uint64())nsame++;
  EXPECT_LE(nsame,1U);
}

TYPED_TEST_P(CoreTests, SameSeeds)
{
  TypeParam core1(1ULL);
  TypeParam core2(1ULL);
  unsigned ndiffer = 0;
  for(unsigned i=0;i<10000;i++)
    if(core1.uniform_uint64() != core2.uniform_uint64())ndiffer++;
  EXPECT_EQ(ndiffer,0U);
}

TYPED_TEST_P(CoreTests, SaveAndRestoreState64)
{
  for(unsigned N=1;N<100;N++)
  {
    TypeParam core(RNG::uint64_from_random_device());
    std::vector<uint64_t> data0;
    std::vector<uint64_t> data1;
    for(unsigned i=0;i<N;i++)
      data0.push_back(core.uniform_uint64());
    calin::ix::math::rng::RNGData proto;
    core.save_to_proto(&proto);
#if 0
    google::protobuf::io::OstreamOutputStream stream(&std::cout);
    if(N==99)google::protobuf::TextFormat::Print(proto, &stream);
#endif
    for(unsigned i=0;i<N;i++)
      data1.push_back(core.uniform_uint64());
    RNGCore* core2 = RNGCore::create_from_proto(proto);
    for(unsigned i=0;i<N;i++)
      EXPECT_EQ(data0[i], core2->uniform_uint64())
          << "Core2 mistmatch - data0[" << i << ']';
    for(unsigned i=0;i<N;i++)
      EXPECT_EQ(data1[i], core2->uniform_uint64())
          << "Core2 mistmatch - data1[" << i << ']';
    RNGCore* core3 = RNGCore::create_from_proto(proto, true);
    for(unsigned i=0;i<N;i++)
      EXPECT_EQ(data1[i], core3->uniform_uint64())
          << "Core3 mistmatch - data1[" << i << ']';
    delete core2;
    delete core3;
  }
}

TYPED_TEST_P(CoreTests, RestoreFromSeedOnlyProto)
{
  uint64_t seed = RNG::uint64_from_random_device();
  TypeParam core(seed);
  calin::ix::math::rng::RNGData proto;
  TypeParam::mutable_core_data(&proto)->set_seed(seed);
#if 0
  google::protobuf::io::OstreamOutputStream stream(&std::cout);
  google::protobuf::TextFormat::Print(proto, &stream);
#endif
  TypeParam core2(TypeParam::core_data(proto));
  for(unsigned i=0;i<1000;i++)
    EXPECT_EQ(core.uniform_uint64(), core2.uniform_uint64());
}

REGISTER_TYPED_TEST_CASE_P(CoreTests, FillsAllBits64, FillsAllBits32,
                           SaveAndRestoreState64, RestoreFromSeedOnlyProto,
                           DifferentSeeds, SameSeeds);
typedef ::testing::Types<NR3RNGCore, Ranlux48RNGCore, MT19937RNGCore> CoreTypes;
INSTANTIATE_TYPED_TEST_CASE_P(TestRNG, CoreTests, CoreTypes);

TEST(TestRNG, SaveAndRestoreStateU32)
{
  for(unsigned N=1;N<100;N++)
  {
    RNG rng;
    std::vector<uint32_t> data0;
    std::vector<uint32_t> data1;
    for(unsigned i=0;i<N;i++)
      data0.push_back(rng.uniform_uint32());
    calin::ix::math::rng::RNGData proto;
    rng.save_to_proto(&proto);
#if 0
    google::protobuf::io::OstreamOutputStream stream(&std::cout);
    if(N==99)google::protobuf::TextFormat::Print(proto, &stream);
#endif
    for(unsigned i=0;i<N;i++)
      data1.push_back(rng.uniform_uint32());
    RNG rng2 = RNG(proto);
    for(unsigned i=0;i<N;i++)
      EXPECT_EQ(data0[i], rng2.uniform_uint32())
          << "RNG2 mistmatch - data0[" << i << ']';
    for(unsigned i=0;i<N;i++)
      EXPECT_EQ(data1[i], rng2.uniform_uint32())
          << "RNG2 mistmatch - data1[" << i << ']';
    RNG rng3 = RNG(proto, true);
    for(unsigned i=0;i<N;i++)
      EXPECT_EQ(data1[i], rng3.uniform_uint32())
          << "RNG3 mistmatch - data1[" << i << ']';
  }
}

TEST(TestRNG, SaveAndRestoreStateBM)
{
  for(unsigned N=1;N<100;N++)
  {
    RNG rng;
    std::vector<double> normal_data0;
    std::vector<double> normal_data1;
    for(unsigned i=0;i<N;i++)
      normal_data0.push_back(rng.normal());
    calin::ix::math::rng::RNGData proto;
    rng.save_to_proto(&proto);
    for(unsigned i=0;i<N;i++)
      normal_data1.push_back(rng.normal());
    RNG rng2(proto);
    for(unsigned i=0;i<N;i++)
      EXPECT_EQ(normal_data0[i], rng2.normal())
          << "rng2 mistmatch - normal_data0[" << i << ']';
    for(unsigned i=0;i<N;i++)
      EXPECT_EQ(normal_data1[i], rng2.normal())
          << "rng2 mistmatch - normal_data1[" << i << ']';
    RNG rng3(proto, true);
    for(unsigned i=0;i<N;i++)
      EXPECT_EQ(normal_data1[i], rng3.normal())
          << "rng3 mistmatch - data1[" << i << ']';
  }
}

TEST(TestRNG, RestoreFromEmptyProto)
{
  calin::ix::math::rng::RNGData proto1;
  RNG rng1(proto1);
  calin::ix::math::rng::RNGData proto2;
  rng1.save_to_proto(&proto2);
#if 0
  google::protobuf::io::OstreamOutputStream stream(&std::cout);
  google::protobuf::TextFormat::Print(proto2, &stream);
#endif
  RNG rng2(proto2);
  for(unsigned i=0;i<1000;i++)
    EXPECT_EQ(rng1.normal(), rng2.normal());
}

TEST(TestRNG, RestoreFromSeedOnlyProto)
{
  uint64_t seed = RNG::uint64_from_random_device();
  RNG rng(seed);
  calin::ix::math::rng::RNGData proto;
  proto.mutable_nr3_core()->set_seed(seed);
  RNG rng2(proto);
  for(unsigned i=0;i<1000;i++)
    EXPECT_EQ(rng.normal(), rng2.normal());
}

template<typename T> std::tuple<double, double, double>
calc_moments(const T& generator, bool print = false,
             RNGCore* core = new NR3RNGCore(RNG::std_test_seed /*RNG::uint64_from_random_device()*/),
             unsigned N = 1000000)
{
  RNG rng(core, true);
  double sum_x = 0;
  double sum_xx = 0;
  double sum_xxx = 0;
  for(unsigned i=0; i<N; i++) {
    double x;
    generator(rng,x);
    sum_x += x;
    sum_xx += x*x;
    sum_xxx += x*x*x;
  }
  double m1 = sum_x/double(N);
  double m2 = sum_xx/double(N) - m1*m1;
  double m3 = sum_xxx/double(N) - 3*sum_xx/double(N)*m1 + 2*m1*m1*m1;
  if(print)std::cout << m1 << ' ' << m2 << ' ' << m3 << '\n';
  return std::make_tuple(m1,m2,m3);
}

TEST(TestRNG, UniformMoments) {
  double m1, m2, m3;
  std::tie(m1,m2,m3) = calc_moments([](RNG& rng,double& x){ x=rng.uniform(); });
  EXPECT_NEAR(m1, 0.5, 0.01);
  EXPECT_NEAR(m2, 1.0/12.0, 0.01);
  EXPECT_NEAR(m3, 0.0, 0.001);
}

TEST(TestRNG, FloatUniformMoments) {
  double m1, m2, m3;
  std::tie(m1,m2,m3) = calc_moments([](RNG& rng,double& x){
      x = rng.uniform_float (); });
  EXPECT_NEAR(m1, 0.5, 0.01);
  EXPECT_NEAR(m2, 1.0/12.0, 0.01);
  EXPECT_NEAR(m3, 0.0, 0.001);
}

TEST(TestRNG, ExponentialMoments) {
  double m1, m2, m3;
  std::tie(m1,m2,m3) = calc_moments([](RNG& rng,double& x){
      x=rng.exponential(); });
  EXPECT_NEAR(m1, 1.0, 0.01);
  EXPECT_NEAR(m2, 1.0, 0.01);
  EXPECT_NEAR(m3, 2.0, 0.05);
}

class ExponentialWithMean : public testing::TestWithParam<double> {
 public:
};

TEST_P(ExponentialWithMean, Moments) {
  double mean = GetParam();
  double m1, m2, m3;
  std::tie(m1,m2,m3) = calc_moments([mean](RNG& rng,double& x){
      x=rng.exponential(mean); });
  EXPECT_NEAR(m1, mean, 0.01*mean);
  EXPECT_NEAR(m2, mean*mean, 0.01*mean*mean);
  EXPECT_NEAR(m3, 2.0*mean*mean*mean, 0.05*mean*mean*mean);
}

INSTANTIATE_TEST_CASE_P(TestRNG,
                        ExponentialWithMean,
                        ::testing::Values(1.0, 2.0, 3.0));

TEST(TestRNG, NormalMoments) {
  double m1, m2, m3;
  std::tie(m1,m2,m3) = calc_moments([](RNG& rng,double& x){
      x=rng.normal(); });
  EXPECT_NEAR(m1, 0.0, 0.01);
  EXPECT_NEAR(m2, 1.0, 0.01);
  EXPECT_NEAR(m3, 0.0, 0.05);
}

class NormalWithMeanAndSigma :
    public testing::TestWithParam<std::pair<double,double>> {
 public:
};

TEST_P(NormalWithMeanAndSigma, Moments) {
  double mean = GetParam().first;
  double sigma = GetParam().second;
  double m1, m2, m3;
  std::tie(m1,m2,m3) = calc_moments([mean,sigma](RNG& rng,double& x){
      x=rng.normal(mean,sigma); });
  EXPECT_NEAR(m1, mean, 0.01*std::abs(mean));
  EXPECT_NEAR(m2, sigma*sigma, 0.01*sigma*sigma);
  EXPECT_NEAR(m3, 0, 0.05*sigma*sigma*sigma);
}

INSTANTIATE_TEST_CASE_P(TestRNG,
                        NormalWithMeanAndSigma,
                        ::testing::Values(std::make_pair(1.0, 2.0),
                                          std::make_pair(-1.0, 2.0)));

class PoissonWithMean : public testing::TestWithParam<double> {
 public:
};

TEST_P(PoissonWithMean, Moments) {
  double mean = GetParam();
  double m1, m2, m3;
  std::tie(m1,m2,m3) = calc_moments([mean](RNG& rng,double& x){
      x=double(rng.poisson(mean)); });
  EXPECT_NEAR(m1, mean, 0.01*mean);
  EXPECT_NEAR(m2, mean, 0.01*mean);
  EXPECT_NEAR(m3, mean, 0.1*mean);
}

// Make sure to test values <5 and >5 since the generator uses
// different algorithms for these cases
INSTANTIATE_TEST_CASE_P(TestRNG,
                        PoissonWithMean,
                        ::testing::Values(1.0, 2.0, 3.0, 4.99999, 5.0,
                                          5.00001, 10.0, 100.0));

class GammaByAlphaAndBeta :
    public testing::TestWithParam<std::pair<double,double>> {
 public:
};

TEST_P(GammaByAlphaAndBeta, Moments) {
  double alpha = GetParam().first;
  double beta = GetParam().second;
  double m1, m2, m3;
  std::tie(m1,m2,m3) = calc_moments([alpha,beta](RNG& rng,double& x){
      x=rng.gamma_by_alpha_and_beta(alpha,beta); });
  EXPECT_NEAR(m1, alpha/beta, 0.01*alpha/beta/beta);
  EXPECT_NEAR(m2, alpha/(beta*beta), 0.05*alpha/(beta*beta*beta));
  EXPECT_NEAR(m3, 2.0*alpha/(beta*beta*beta), 0.1*alpha/(beta*beta*beta));
}

// Test alpha <1 and alpha>1
INSTANTIATE_TEST_CASE_P(TestRNG,
                        GammaByAlphaAndBeta,
                        ::testing::Values(std::make_pair(0.1, 1.0),
                                          std::make_pair(0.99999, 1.0),
                                          std::make_pair(1.0, 1.0),
                                          std::make_pair(1.00001, 1.0),
                                          std::make_pair(1.0, 2.0),
                                          std::make_pair(2.0, 1.0)));

class GammaByMeanAndSigma :
    public testing::TestWithParam<std::pair<double,double>> {
 public:
};

TEST_P(GammaByMeanAndSigma, Moments) {
  double mean = GetParam().first;
  double sigma = GetParam().second;
  double m1, m2, m3;
  std::tie(m1,m2,m3) = calc_moments([mean,sigma](RNG& rng,double& x){
      x=rng.gamma_by_mean_and_sigma(mean,sigma); });
  EXPECT_NEAR(m1, mean, 0.01*mean);
  EXPECT_NEAR(m2, sigma*sigma, 0.01*sigma*sigma);
  EXPECT_NEAR(m3, 2.0*(sigma*sigma*sigma*sigma)/mean,
              0.1*(sigma*sigma*sigma*sigma)/mean);
}

INSTANTIATE_TEST_CASE_P(TestRNG,
                        GammaByMeanAndSigma,
                        ::testing::Values(std::make_pair(1.0, 0.1),
                                          std::make_pair(1.0, 1.0),
                                          std::make_pair(2.0, 1.0)));

class PolyaByMeanAndExcessSigma :
    public testing::TestWithParam<std::pair<double,double>> {
 public:
};

TEST_P(PolyaByMeanAndExcessSigma, Moments) {
  double mean = GetParam().first;
  double xs_sigma = GetParam().second;
  double m1, m2, m3;
  std::tie(m1,m2,m3) = calc_moments([mean,xs_sigma](RNG& rng,double& x){
      x=rng.polya(mean,xs_sigma); });
  EXPECT_NEAR(m1, mean, 0.01*mean);
  EXPECT_NEAR(m2, mean+xs_sigma*xs_sigma, 0.01*(mean+xs_sigma*xs_sigma));
  EXPECT_NEAR(m3, (mean+2.0*xs_sigma*xs_sigma)*(mean+xs_sigma*xs_sigma)/mean,
              0.05*(mean+2.0*xs_sigma*xs_sigma)*(mean+xs_sigma*xs_sigma)/mean);
}

INSTANTIATE_TEST_CASE_P(TestRNG,
                        PolyaByMeanAndExcessSigma,
                        ::testing::Values(std::make_pair(1.0, 0.1),
                                          std::make_pair(1.0, 1.0),
                                          std::make_pair(2.0, 1.0),
                                          std::make_pair(100.0, 10.0)));

class Binomial :
    public testing::TestWithParam<std::pair<double,int>> {
 public:
};

TEST_P(Binomial, Moments) {
  double p = GetParam().first;
  int n = GetParam().second;
  double np = double(n)*p;
  double m1, m2, m3;
  std::tie(m1,m2,m3) = calc_moments([p,n](RNG& rng,double& x) {
      x=rng.binomial(p,n); });
  EXPECT_NEAR(m1, np, 0.01*sqrt(np*(1-p)));
  EXPECT_NEAR(m2, np*(1-p), 0.05*(np*(1-p)));
  EXPECT_NEAR(m3, (1-2*p)*np*(1-p),  0.05*sqrt(np*(1-p))*(np*(1-p)));
}

INSTANTIATE_TEST_CASE_P(TestRNG,
                        Binomial,
                        ::testing::Values(std::make_pair(0.1, 1),
                                          std::make_pair(0.9, 1),
                                          std::make_pair(0.1, 24),
                                          std::make_pair(0.1, 25),
                                          std::make_pair(0.001, 100),
                                          std::make_pair(0.0005, 1000),
                                          std::make_pair(0.1, 1000),
                                          std::make_pair(0.1, 1000000)));


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
