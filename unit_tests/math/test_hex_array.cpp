/*

   calin/unit_tests/math/test_hex_array.cpp -- Stephen Fegan -- 2015-10-21

   Unit tests for hex array classes

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
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
#include <vector>

#include "math/hex_array.hpp"
#include "math/hex_array_simd.hpp"

using namespace calin::math::hex_array;
using namespace calin::math::hex_array::vvv;

constexpr unsigned NLOOP_SPEED_TEST_50RINGS = 10000;

TEST(TestHexArray, HexIDToRingIDLoop_SpeedTest50Rings) {
  for(unsigned iloop = 0; iloop<NLOOP_SPEED_TEST_50RINGS; iloop++)
  {
    unsigned hexid = 1;
    for(unsigned iring=1;iring<50;iring++)
      for(unsigned ichan=0;ichan<6*iring;ichan++)
      {
        EXPECT_EQ(iring, positive_hexid_to_ringid_loop(hexid));
        hexid++;
      }
  }
}

TEST(TestHexArray, HexIDToRingIDRoot_SpeedTest50Rings) {
  for(unsigned iloop = 0; iloop<NLOOP_SPEED_TEST_50RINGS; iloop++)
  {
    unsigned hexid = 1;
    for(unsigned iring=1;iring<50;iring++)
      for(unsigned ichan=0;ichan<6*iring;ichan++)
      {
        EXPECT_EQ(iring, positive_hexid_to_ringid_root(hexid));
        hexid++;
      }
  }
}

#if defined(__AVX2__) and defined(__FMA__)
TEST(TestHexArray, AVX2_HexIDToRingIDRoot_SpeedTest50Rings) {
  for(unsigned iloop = 0; iloop<NLOOP_SPEED_TEST_50RINGS/8; iloop++)
  {
    unsigned hexid = 1;
    for(unsigned iring=1;iring<50;iring++)
      for(unsigned ichan=0;ichan<6*iring;ichan++)
      {
        EXPECT_EQ(iring, test_avx2_positive_hexid_to_ringid_root(hexid));
        hexid++;
      }
  }
}
#endif

TEST(TestHexArray, HexIDToRingIDRoot_2000Rings) {
  unsigned hexid = 1;
  for(unsigned iring=1;iring<2000;iring++)
    for(unsigned ichan=0;ichan<6*iring;ichan++)
    {
      ASSERT_EQ(iring, positive_hexid_to_ringid_root(hexid))
        << "With hexid=" << hexid;
      hexid++;
    }
}

#if defined(__AVX2__) and defined(__FMA__)
TEST(TestHexArray, AVX2_HexIDToRingIDRoot_2000Rings) {
  unsigned hexid = 1;
  for(unsigned iring=1;iring<2000;iring++)
    for(unsigned ichan=0;ichan<6*iring;ichan++)
    {
      ASSERT_EQ(iring, test_avx2_positive_hexid_to_ringid_root(hexid))
        << "With hexid=" << hexid;
      hexid++;
    }
}
#endif

TEST(TestHexArray, SomeNeighbors) {
  EXPECT_EQ(hexid_to_neighbor_hexids(0),
            std::vector<unsigned>({1,2,3,4,5,6}));
  EXPECT_EQ(hexid_to_neighbor_hexids(1),
            std::vector<unsigned>({7,8,2,0,6,18}));
  EXPECT_EQ(hexid_to_neighbor_hexids(7),
            std::vector<unsigned>({19,20,8,1,18,36}));
  EXPECT_EQ(hexid_to_neighbor_hexids(8),
            std::vector<unsigned>({20,21,9,2,1,7}));
  EXPECT_EQ(hexid_to_neighbor_hexids(113),
            std::vector<unsigned>({80,79,112,152,153,114}));
}

TEST(TestHexArray, HexIDToXY_NewCodeSpeedTest) {
  for(unsigned iloop=0;iloop<1000;iloop++)
  {
    unsigned hexid = 1;
    for(unsigned iring=1;iring<50;iring++)
      for(unsigned ichan=0;ichan<6*iring;ichan++)
      {
        double x1,y1;
        //double x2,y2;
        hexid_to_xy(hexid, x1, y1, true);
        //int vvv_hexid = hexid+1;
        //nh_to_xy(&vvv_hexid, &x2, &y2);
        //EXPECT_NEAR(x1,x2,1e-6);
        hexid++;
      }
  }
}

TEST(TestHexArray, HexIDToXY_VVVCodeSpeedTest) {
  for(unsigned iloop=0;iloop<1000;iloop++)
  {
    unsigned hexid = 1;
    for(unsigned iring=1;iring<50;iring++)
      for(unsigned ichan=0;ichan<6*iring;ichan++)
      {
        //double x1,y1;
        double x2,y2;
        //hexid_to_xy(hexid, x1, y1, true);
        int vvv_hexid = hexid+1;
        nh_to_xy(&vvv_hexid, &x2, &y2);
        //EXPECT_NEAR(x1,x2,1e-6);
        hexid++;
      }
  }
}

TEST(TestHexArray, HexIDToXY_ComparisonWithVVVCode) {
  unsigned hexid = 1;
  for(unsigned iring=1;iring<500;iring++)
    for(unsigned ichan=0;ichan<6*iring;ichan++)
    {
      double x1,y1;
      double x2,y2;
      hexid_to_xy(hexid, x1, y1, true);
      int vvv_hexid = hexid+1;
      nh_to_xy(&vvv_hexid, &x2, &y2);
      EXPECT_NEAR(x1,x2,1e-6);
      EXPECT_NEAR(y1,y2,1e-6);
      hexid++;
    }
}

TEST(TestHexArray, XYToHexID_NewCodeSpeedTest) {
  for(double x=-10.005; x<10.015; x+=0.01)
    for(double y=-10.005; y<10.015; y+=0.01)
    {
      double xx1 = x;
      double yy1 = y;
      /* unsigned hexid = */ xy_to_hexid_with_remainder(xx1, yy1, true);
      //double xx2 = x;
      //double yy2 = y;
      //int hexid2;
      //xy_to_nh(&xx2,&yy2,&hexid2);
    }
}

TEST(TestHexArray, XYToHexID_VVVCodeSpeedTest) {
  for(double x=-10.005; x<10.015; x+=0.01)
    for(double y=-10.005; y<10.015; y+=0.01)
    {
      //double xx1 = x;
      //double yy1 = y;
      //unsigned hexid = xy_to_hexid_with_remainder(xx1, yy1, true);
      double xx2 = x;
      double yy2 = y;
      int hexid2;
      xy_to_nh(&xx2,&yy2,&hexid2);
    }
}

TEST(TestHexArray, XYToHexID_ComparisonWithVVVCode) {
  for(double x=-10.005; x<10.015; x+=0.02)
    for(double y=-10.005; y<10.015; y+=0.02)
    {
      double xx1 = x;
      double yy1 = y;
      unsigned hexid = xy_to_hexid_with_remainder(xx1, yy1, true);
      double xx2 = x;
      double yy2 = y;
      int hexid2;
      xy_to_nh(&xx2,&yy2,&hexid2);
      EXPECT_EQ(hexid, (unsigned)hexid2-1);
      EXPECT_NEAR(xx1,xx2,1e-6);
      EXPECT_NEAR(yy1,yy2,1e-6);
    }
}

#if defined(__AVX2__) and defined(__FMA__)
TEST(TestHexArray, AVX2_HexIDToFromRingSegRun_APriori) {
  unsigned hexid=1;
  for(unsigned iring=1;iring<100;iring++)
    for(unsigned iseg=0;iseg<6;iseg++)
      for(unsigned irun=0;irun<iring;irun++)
      {
        ASSERT_EQ(hexid,
          test_avx2_positive_ringid_segid_runid_to_hexid(iring, iseg, irun)) <<
            iring << ' ' << iseg << ' ' << irun;
        unsigned ringid = 0;
        unsigned segid = 0;
        unsigned runid = 0;
        test_avx2_positive_hexid_to_ringid_segid_runid(hexid, ringid, segid, runid);
        ASSERT_EQ(iring, ringid);
        ASSERT_EQ(iseg, segid);
        ASSERT_EQ(irun, runid);
        hexid++;
      }
}
#endif

TEST(TestHexArray, HexIDToFromRingSegRun_APriori) {
  unsigned hexid=1;
  for(unsigned iring=1;iring<100;iring++)
    for(unsigned iseg=0;iseg<6;iseg++)
      for(unsigned irun=0;irun<iring;irun++)
      {
        ASSERT_EQ(hexid,
          positive_ringid_segid_runid_to_hexid(iring, iseg, irun)) <<
            iring << ' ' << iseg << ' ' << irun;
        unsigned ringid = 0;
        unsigned segid = 0;
        unsigned runid = 0;
        positive_hexid_to_ringid_segid_runid(hexid, ringid, segid, runid);
        ASSERT_EQ(iring, ringid);
        ASSERT_EQ(iseg, segid);
        ASSERT_EQ(irun, runid);
        hexid++;
      }
}

TEST(TestHexArray, HexIDToFromRingSegRun_EQ) {
  for(unsigned hexid=1;hexid<100000;hexid++)
  {
    unsigned ringid = 0;
    unsigned segid = 0;
    unsigned runid = 0;
    positive_hexid_to_ringid_segid_runid(hexid, ringid, segid, runid);
    ASSERT_EQ(hexid,
      positive_ringid_segid_runid_to_hexid(ringid, segid, runid))
        << ringid << ' ' << segid << ' ' << runid;
  }
}

#if defined(__AVX2__) and defined(__FMA__)
TEST(TestHexArray, AVX2_HexIDToFromRingSegRun_EQ) {
  for(unsigned hexid=1;hexid<100000;hexid++)
  {
    unsigned ringid = 0;
    unsigned segid = 0;
    unsigned runid = 0;
    test_avx2_positive_hexid_to_ringid_segid_runid(hexid, ringid, segid, runid);
    ASSERT_EQ(hexid,
      test_avx2_positive_ringid_segid_runid_to_hexid(ringid, segid, runid))
        << ringid << ' ' << segid << ' ' << runid;
  }
}
#endif

TEST(TestHexArray, HexIDToFromUV_CW_EQ) {
  for(unsigned hexid=1;hexid<100000;hexid++)
  {
    int u = 0;
    int v = 0;
    hexid_to_uv_cw(hexid, u, v);
    ASSERT_EQ(hexid, uv_to_hexid_cw(u, v))
      << hexid << ' ' << u << ' ' << v;
  }
}

TEST(TestHexArray, HexIDToFromUV_CCW_EQ) {
  for(unsigned hexid=1;hexid<100000;hexid++)
  {
    int u = 0;
    int v = 0;
    hexid_to_uv_ccw(hexid, u, v);
    ASSERT_EQ(hexid, uv_to_hexid_ccw(u, v))
      << hexid << ' ' << u << ' ' << v;
  }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
