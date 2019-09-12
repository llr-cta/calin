/*

   calin/unit_tests/math/test_accumulator.cpp -- Stephen Fegan -- 2018-10-13

   Unit tests for Accumulator class

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
#include <cmath>

#include <math/accumulator.hpp>
#include <util/vcl.hpp>

using namespace calin::math::accumulator;
using namespace calin::util::vcl;

TEST(TestAccumulator, Double_OnePlusOne) {
  KahanAccumulator acc;
  acc.accumulate(1.0);
  acc.accumulate(1.0);
  EXPECT_EQ(acc.total(), 2.0);
}

TEST(TestAccumulator, Double_OnePlusEpsilon) {
  KahanAccumulator acc;
  acc.accumulate(1.0);
  acc.accumulate(DBL_EPSILON);
  EXPECT_EQ(acc.total(), 1.0+DBL_EPSILON);
  EXPECT_EQ(acc.correction(), 0.0);
}

TEST(TestAccumulator, Double_OnePlusHalfEpsilon) {
  KahanAccumulator acc;
  acc.accumulate(1.0);
  acc.accumulate(0.5*DBL_EPSILON);
  EXPECT_EQ(acc.total(), 1.0);
  EXPECT_EQ(acc.correction(), 0.5*DBL_EPSILON);
}

TEST(TestAccumulator, Float_OnePlusOne) {
  BasicKahanAccumulator<float> acc;
  acc.accumulate(1.0f);
  acc.accumulate(1.0f);
  EXPECT_EQ(acc.total(), 2.0f);
}

TEST(TestAccumulator, Float_OnePlusEpsilon) {
  BasicKahanAccumulator<float> acc;
  acc.accumulate(1.0f);
  acc.accumulate(FLT_EPSILON);
  EXPECT_EQ(acc.total(), 1.0f+FLT_EPSILON);
  EXPECT_EQ(acc.correction(), 0.0f);
}

TEST(TestAccumulator, Float_OnePlusHalfEpsilon) {
  BasicKahanAccumulator<float> acc;
  acc.accumulate(1.0f);
  acc.accumulate(0.5f*FLT_EPSILON);
  EXPECT_EQ(acc.total(), 1.0f);
  EXPECT_EQ(acc.correction(), 0.5f*FLT_EPSILON);
}

TEST(TestAccumulator, Vec4d_OnePlusOne) {
  BasicKahanAccumulator<Vec4d> acc;
  acc.accumulate(1.0);
  acc.accumulate(1.0);
  EXPECT_TRUE(horizontal_and(acc.total() == Vec4d(2.0)));
}

TEST(TestAccumulator, Vec4d_OnePlusEpsilon) {
  BasicKahanAccumulator<Vec4d> acc;
  acc.accumulate(1.0);
  acc.accumulate(DBL_EPSILON);
  EXPECT_TRUE(horizontal_and(acc.total() == Vec4d(1.0+DBL_EPSILON)));
  EXPECT_TRUE(horizontal_and(acc.correction() == Vec4d(0.0)));
}

TEST(TestAccumulator, Vec4d_OnePlusHalfEpsilon) {
  BasicKahanAccumulator<Vec4d> acc;
  acc.accumulate(1.0);
  acc.accumulate(0.5*DBL_EPSILON);
  EXPECT_TRUE(horizontal_and(acc.total() == Vec4d(1.0)));
  EXPECT_TRUE(horizontal_and(acc.correction() == Vec4d(0.5*DBL_EPSILON)));
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
