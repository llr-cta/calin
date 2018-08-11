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
}

TEST(TestVCL, RNG256_FLT_EXP) {
  std::string fixture = ::testing::UnitTest::GetInstance()->current_test_info()->name();
  VCLRNG<VCL256Architecture> rng(fixture,"rng");
  std::cout << rng.exponential_float() << '\n';
  std::cout << rng.exponential_float() << '\n';
}

TEST(TestVCL, RNG256_FLT_SINCOS) {
  std::string fixture = ::testing::UnitTest::GetInstance()->current_test_info()->name();
  VCLRNG<VCL256Architecture> rng(fixture,"rng");
  Vec8f s, c;
  rng.sincos_float(s,c);
  std::cout << s << " -- " << c << " -- " << s*s+c*c << '\n';
  rng.sincos_float(s,c);
  std::cout << s << " -- " << c << " -- " << s*s+c*c << '\n';
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
