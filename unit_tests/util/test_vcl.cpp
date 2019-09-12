/*

   calin/unit_tests/util/test_vcl.cpp -- Stephen Fegan -- 2018-08-08

   Unit tests for vcl

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

TEST(TestVCL, Transpose128_U32) {
  Vec4ui x[4];
  for(unsigned i=0;i<4;i++)
    x[i] = Vec4ui(i*4+0,i*4+1,i*4+2,i*4+3);
  transpose(x);
  for(unsigned j=0;j<4;j++)
    for(unsigned i=0;i<4;i++)
      EXPECT_EQ(x[j][i], i*4+j);
}

TEST(TestVCL, Transpose128_U64) {
  Vec2uq x[2];
  for(unsigned i=0;i<2;i++)
    x[i] = Vec2uq(i*2+0,i*2+1);
  transpose(x);
  for(unsigned j=0;j<2;j++)
    for(unsigned i=0;i<2;i++)
      EXPECT_EQ(x[j][i], i*2+j);
}

TEST(TestVCL, Transpose128_FLT) {
  Vec4f x[4];
  for(unsigned i=0;i<4;i++)
    x[i] = Vec4f(i+0,i+0.1,i+0.2,i+0.3);
  transpose(x);
  for(unsigned j=0;j<4;j++)
    for(unsigned i=0;i<4;i++)
      EXPECT_NEAR(x[j][i], i+0.1*j, 0.01);
}

TEST(TestVCL, Transpose128_DBL) {
  Vec2d x[2];
  for(unsigned i=0;i<2;i++)
    x[i] = Vec2d(i*2+0,i*2+0.1);
  transpose(x);
  for(unsigned j=0;j<2;j++)
    for(unsigned i=0;i<2;i++)
      EXPECT_NEAR(x[j][i], i*2+0.1*j, 0.01);
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

TEST(TestVCL, Transpose256_U32) {
  Vec8ui x[8];
  for(unsigned i=0;i<8;i++)
    x[i] = Vec8ui(i*8+0,i*8+1,i*8+2,i*8+3,i*8+4,i*8+5,i*8+6,i*8+7);
  transpose(x);
  for(unsigned j=0;j<8;j++)
    for(unsigned i=0;i<8;i++)
      EXPECT_EQ(x[j][i], i*8+j);
}

TEST(TestVCL, Transpose256_U64) {
  Vec4uq x[4];
  for(unsigned i=0;i<4;i++)
    x[i] = Vec4uq(i*4+0,i*4+1,i*4+2,i*4+3);
  transpose(x);
  for(unsigned j=0;j<4;j++)
    for(unsigned i=0;i<4;i++)
      EXPECT_EQ(x[j][i], i*4+j);
}

TEST(TestVCL, Transpos256_FLT) {
  Vec8f x[8];
  for(unsigned i=0;i<8;i++)
    x[i] = Vec8f(i+0.0,i+0.1,i+0.2,i+0.3,i+0.4,i+0.5,i+0.6,i+0.7);
  transpose(x);
  for(unsigned j=0;j<8;j++)
    for(unsigned i=0;i<8;i++)
      EXPECT_NEAR(x[j][i], i+0.1*j, 0.01);
}

TEST(TestVCL, Transpose256_DBL) {
  Vec4d x[4];
  for(unsigned i=0;i<4;i++)
    x[i] = Vec4d(i+0,i+0.1,i+0.2,i+0.3);
  transpose(x);
  for(unsigned j=0;j<4;j++)
    for(unsigned i=0;i<4;i++)
      EXPECT_NEAR(x[j][i], i+0.1*j, 0.01);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
