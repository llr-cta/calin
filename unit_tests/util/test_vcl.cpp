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
#include <math/accumulator.hpp>
#include <provenance/chronicle.hpp>

using namespace calin::util::vcl;
using namespace calin::math::rng;
using namespace calin::math::accumulator;

TEST(TestVCL, Print) {
  Vec8f v(1.234);
  std::cout << v << '\n';
  Vec8i vi(-101);
  std::cout << vi << '\n';
  std::cout << to_float(vi) << '\n';
  std::cout << to_float(vi) + 0.5 << '\n';
}

TEST(TestVCL, Transpose128_U16) {
  Vec8us x[8];
  for(unsigned i=0;i<8;i++)
    x[i] = Vec8us(i*8+0,i*8+1,i*8+2,i*8+3,i*8+4,i*8+5,i*8+6,i*8+7);
  transpose(x);
  for(unsigned j=0;j<8;j++)
    for(unsigned i=0;i<8;i++)
      EXPECT_EQ(x[j][i], i*8+j);
}

TEST(TestVCL, Transpose256_U16) {
  Vec16us x[16];
  for(unsigned i=0;i<16;i++)
    x[i] = Vec16us(i*16+0,i*16+1,i*16+2,i*16+3,i*16+4,i*16+5,i*16+6,i*16+7,
      i*16+8,i*16+9,i*16+10,i*16+11,i*16+12,i*16+13,i*16+14,i*16+15);
  transpose(x);
  for(unsigned j=0;j<16;j++)
    for(unsigned i=0;i<16;i++)
      EXPECT_EQ(x[j][i], i*16+j);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
