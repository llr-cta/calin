/*

   calin/unit_tests/util/test_vcl.cpp -- Stephen Fegan -- 2018-08-08

   Unit tests for vcl

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

#include <util/vcl.hpp>
#include <math/rng_vcl.hpp>
#include <provenance/chronicle.hpp>

using namespace calin::util::vcl;
using namespace calin::math::rng;

TEST(TestVCL, Print) {
  Vec8f v(1.234);
  std::cout << v << '\n';
  Vec8i vi(-101);
  std::cout << vi << '\n';
  std::cout << to_float(vi) << '\n';
  std::cout << to_float(vi) + 0.5 << '\n';
}

TEST(TestVCL, RNG128) {
  std::string fixture = ::testing::UnitTest::GetInstance()->current_test_info()->name();
  uint64_t seeds[2] = { 1, 2 };
  NR3_VCLRNGCore<VCL128Architecture> rng(seeds,fixture,"rng");
  NR3RNGCore rng0(seeds[0],fixture,"rng0");
  NR3RNGCore rng1(seeds[1],fixture,"rng1");
  std::cout << rng.uniform_uint64() << ' '
    << rng0.uniform_uint64() << ' '
    << rng1.uniform_uint64() << '\n';

  std::cout << rng.uniform_uint64() << ' '
    << rng0.uniform_uint64() << ' '
    << rng1.uniform_uint64() << '\n';
}

TEST(TestVCL, RNG256) {
  std::string fixture = ::testing::UnitTest::GetInstance()->current_test_info()->name();
  uint64_t seeds[4] = { 1, 2, 3, 4 };
  NR3_VCLRNGCore<VCL256Architecture> rng(seeds,fixture,"rng");
  NR3RNGCore rng0(seeds[0],fixture,"rng0");
  NR3RNGCore rng1(seeds[1],fixture,"rng1");
  NR3RNGCore rng2(seeds[2],fixture,"rng2");
  NR3RNGCore rng3(seeds[3],fixture,"rng3");
  std::cout << rng.uniform_uint64() << ' '
    << rng0.uniform_uint64() << ' '
    << rng1.uniform_uint64() << ' '
    << rng2.uniform_uint64() << ' '
    << rng3.uniform_uint64() << '\n';

  std::cout << rng.uniform_uint64() << ' '
    << rng0.uniform_uint64() << ' '
    << rng1.uniform_uint64() << ' '
    << rng2.uniform_uint64() << ' '
    << rng3.uniform_uint64() << '\n';
}

TEST(TestVCL, RNG512) {
  std::string fixture = ::testing::UnitTest::GetInstance()->current_test_info()->name();
  uint64_t seeds[8] = { 1, 2, 3, 4, 5, 6, 7, 8 };
  NR3_VCLRNGCore<VCL512Architecture> rng(seeds,fixture,"rng");
  NR3RNGCore rng0(seeds[0],fixture,"rng0");
  NR3RNGCore rng1(seeds[1],fixture,"rng1");
  NR3RNGCore rng2(seeds[2],fixture,"rng2");
  NR3RNGCore rng3(seeds[3],fixture,"rng3");
  NR3RNGCore rng4(seeds[4],fixture,"rng4");
  NR3RNGCore rng5(seeds[5],fixture,"rng5");
  NR3RNGCore rng6(seeds[6],fixture,"rng6");
  NR3RNGCore rng7(seeds[7],fixture,"rng7");
  std::cout << rng.uniform_uint64() << ' '
    << rng0.uniform_uint64() << ' '
    << rng1.uniform_uint64() << ' '
    << rng2.uniform_uint64() << ' '
    << rng3.uniform_uint64() << ' '
    << rng4.uniform_uint64() << ' '
    << rng5.uniform_uint64() << ' '
    << rng6.uniform_uint64() << ' '
    << rng7.uniform_uint64() << '\n';

  std::cout << rng.uniform_uint64() << ' '
    << rng0.uniform_uint64() << ' '
    << rng1.uniform_uint64() << ' '
    << rng2.uniform_uint64() << ' '
    << rng3.uniform_uint64() << ' '
    << rng4.uniform_uint64() << ' '
    << rng5.uniform_uint64() << ' '
    << rng6.uniform_uint64() << ' '
    << rng7.uniform_uint64() << '\n';
}

TEST(TestVCL, RNG256_FLT) {
  std::string fixture = ::testing::UnitTest::GetInstance()->current_test_info()->name();
  VCLRNG<VCL256Architecture> rng(fixture,"rng");
  std::cout << rng.uniform_float() << '\n';
  std::cout << rng.uniform_float() << '\n';

  unsigned N = 10000000;
  Vec8f x;
  Vec8f sum_x(0);
  Vec8f sum_xx(0);
  Vec8f max_x(-100);
  Vec8f min_x(100);
  for(unsigned i=0;i<N;i++)
  {
    x = rng.uniform_float();
    sum_x += x;
    sum_xx += x*x;
    max_x = max(x, max_x);
    min_x = min(x, min_x);
  }
  std::cout << sum_x << " -- " << sum_xx << '\n';
  std::cout << min_x << " -- " << max_x << " -- " << 1-max_x << '\n';
  float mean_x = horizontal_add(sum_x)/(8*N);
  float var_x = horizontal_add(sum_xx)/(8*N) - mean_x*mean_x;
  std::cout << mean_x << ' ' << var_x << '\n';
}

TEST(TestVCL, RNG256_FLT_EXP) {
  std::string fixture = ::testing::UnitTest::GetInstance()->current_test_info()->name();
  VCLRNG<VCL256Architecture> rng(fixture,"rng");
  std::cout << rng.exponential_float() << '\n';
  std::cout << rng.exponential_float() << '\n';
}

TEST(TestVCL, RNG256_FLT_THETA) {
  std::string fixture = ::testing::UnitTest::GetInstance()->current_test_info()->name();
  VCLRNG<VCL256Architecture> rng(fixture,"rng");
  Vec8f t, s, c;
  t = rng.uniform_float_zc(2.0*M_PI); s=sincos(&c,t);
  std::cout << t << " -- " << s << " -- " << c << " -- " << s*s+c*c << '\n';
  t = rng.uniform_float_zc(2.0*M_PI); s=sincos(&c,t);
  std::cout << t << " -- " << s << " -- " << c << " -- " << s*s+c*c << '\n';
  unsigned N = 10000000;
  Vec8f sum_s(0);
  Vec8f sum_ss(0);
  Vec8f sum_c(0);
  Vec8f sum_cc(0);
  for(unsigned i=0;i<N;i++)
  {
    t = rng.uniform_float_zc(2.0*M_PI);
    s=sincos(&c,t);
    sum_s += s;
    sum_ss += s*s;
    sum_c += c;
    sum_cc += c*c;
  }
  std::cout << sum_s << " -- " << sum_ss << '\n';
  std::cout << sum_c << " -- " << sum_cc << '\n';
  float mean_s = horizontal_add(sum_s)/(8*N);
  float var_s = horizontal_add(sum_ss)/(8*N) - mean_s*mean_s;
  float mean_c = horizontal_add(sum_c)/(8*N);
  float var_c = horizontal_add(sum_cc)/(8*N) - mean_c*mean_c;
  std::cout << mean_s << ' ' << mean_c << ' ' << var_s << ' ' << var_c << '\n';
}

TEST(TestVCL, RNG256_FLT_SINCOS) {
  std::string fixture = ::testing::UnitTest::GetInstance()->current_test_info()->name();
  VCLRNG<VCL256Architecture> rng(fixture,"rng");
  Vec8f s, c;
  rng.sincos_float(s,c);
  std::cout << s << " -- " << c << " -- " << s*s+c*c << '\n';
  rng.sincos_float(s,c);
  std::cout << s << " -- " << c << " -- " << s*s+c*c << '\n';

  unsigned N = 10000000;
  Vec8f sum_s(0);
  Vec8f sum_ss(0);
  Vec8f sum_c(0);
  Vec8f sum_cc(0);
  for(unsigned i=0;i<N;i++)
  {
    rng.sincos_float(s,c);
    sum_s += s;
    sum_ss += s*s;
    sum_c += c;
    sum_cc += c*c;
  }
  std::cout << sum_s << " -- " << sum_ss << '\n';
  std::cout << sum_c << " -- " << sum_cc << '\n';
  float mean_s = horizontal_add(sum_s)/(8*N);
  float var_s = horizontal_add(sum_ss)/(8*N) - mean_s*mean_s;
  float mean_c = horizontal_add(sum_c)/(8*N);
  float var_c = horizontal_add(sum_cc)/(8*N) - mean_c*mean_c;
  std::cout << mean_s << ' ' << mean_c << ' ' << var_s << ' ' << var_c << '\n';
}

TEST(TestVCL, RNG256_DBL_SINCOS) {
  std::string fixture = ::testing::UnitTest::GetInstance()->current_test_info()->name();
  VCLRNG<VCL256Architecture> rng(fixture,"rng");
  Vec4d s, c;
  rng.sincos_double(s,c);
  std::cout << s << " -- " << c << " -- " << s*s+c*c << '\n';
  rng.sincos_double(s,c);
  std::cout << s << " -- " << c << " -- " << s*s+c*c << '\n';

  unsigned N = 10000000;
  Vec4d sum_s(0);
  Vec4d sum_ss(0);
  Vec4d sum_c(0);
  Vec4d sum_cc(0);
  for(unsigned i=0;i<N;i++)
  {
    rng.sincos_double(s,c);
    sum_s += s;
    sum_ss += s*s;
    sum_c += c;
    sum_cc += c*c;
  }
  std::cout << sum_s << " -- " << sum_ss << '\n';
  std::cout << sum_c << " -- " << sum_cc << '\n';
  double mean_s = horizontal_add(sum_s)/(4*N);
  double var_s = horizontal_add(sum_ss)/(4*N) - mean_s*mean_s;
  double mean_c = horizontal_add(sum_c)/(4*N);
  double var_c = horizontal_add(sum_cc)/(4*N) - mean_c*mean_c;
  std::cout << mean_s << ' ' << mean_c << ' ' << var_s << ' ' << var_c << '\n';
}

TEST(TestVCL, RNG256_FLT_NORMAL_BM) {
  std::string fixture = ::testing::UnitTest::GetInstance()->current_test_info()->name();
  VCLRNG<VCL256Architecture> rng(fixture,"rng");
  Vec8f x, y;
  rng.normal_two_float_bm(x,y);
  std::cout << x << " -- " << y << '\n';
  rng.normal_two_float_bm(x,y);
  std::cout << x << " -- " << y << '\n';

  unsigned N = 10000000;
  Vec8f sum_x(0);
  Vec8f sum_xx(0);
  Vec8f sum_y(0);
  Vec8f sum_yy(0);
  for(unsigned i=0;i<N;i++)
  {
    rng.normal_two_float_bm(x,y);
    sum_x += x;
    sum_xx += x*x;
    sum_y += y;
    sum_yy += y*y;
  }
  std::cout << sum_x << " -- " << sum_xx << '\n';
  std::cout << sum_y << " -- " << sum_yy << '\n';
  float mean_x = horizontal_add(sum_x)/(8*N);
  float var_x = horizontal_add(sum_xx)/(8*N) - mean_x*mean_x;
  float mean_y = horizontal_add(sum_y)/(8*N);
  float var_y = horizontal_add(sum_yy)/(8*N) - mean_y*mean_y;
  std::cout << mean_x << ' ' << mean_y << ' ' << var_x << ' ' << var_y << '\n';
}

TEST(TestVCL, RNG256_DBL_NORMAL_BM) {
  std::string fixture = ::testing::UnitTest::GetInstance()->current_test_info()->name();
  VCLRNG<VCL256Architecture> rng(fixture,"rng");
  Vec4d x, y;
  rng.normal_two_double_bm(x,y);
  std::cout << x << " -- " << y << '\n';
  rng.normal_two_double_bm(x,y);
  std::cout << x << " -- " << y << '\n';

  unsigned N = 10000000;
  Vec4d sum_x(0);
  Vec4d sum_xx(0);
  Vec4d sum_y(0);
  Vec4d sum_yy(0);
  for(unsigned i=0;i<N;i++)
  {
    rng.normal_two_double_bm(x,y);
    sum_x += x;
    sum_xx += x*x;
    sum_y += y;
    sum_yy += y*y;
  }
  std::cout << sum_x << " -- " << sum_xx << '\n';
  std::cout << sum_y << " -- " << sum_yy << '\n';
  double mean_x = horizontal_add(sum_x)/(4*N);
  double var_x = horizontal_add(sum_xx)/(4*N) - mean_x*mean_x;
  double mean_y = horizontal_add(sum_y)/(4*N);
  double var_y = horizontal_add(sum_yy)/(4*N) - mean_y*mean_y;
  std::cout << mean_x << ' ' << mean_y << ' ' << var_x << ' ' << var_y << '\n';
}

TEST(TestVCL, ZZ_DumpChronicle) {
  calin::ix::provenance::chronicle::Chronicle* c =
    calin::provenance::chronicle::copy_the_chronicle();
  //std::cout << c->DebugString();
  delete c;
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
