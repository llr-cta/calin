/*

   calin/unit_tests/math/test_hex_array.cpp -- Stephen Fegan -- 2015-10-21

   Unit tests for hex array classes

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <random>
#include <iostream>
#include <iomanip>
#include <gtest/gtest.h>
#include <vector>

#include "math/hex_array.hpp"
// #include "math/hex_array_simd.hpp"

using namespace calin::math::hex_array;
using namespace calin::math::hex_array::vvv;

constexpr unsigned NLOOP_SPEED_TEST_50RINGS = 10000;
constexpr unsigned NLOOP_HEXID_TOFROM_UV = 1000;

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
#if 0
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

#if 0
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
        hexid_to_xy(hexid, x1, y1, true);
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

#if 0
#if defined(__AVX2__) and defined(__FMA__)
TEST(TestHexArray, AVX2_HexIDToXY_Equals_Scalar) {
  unsigned hexid = 1;
  for(unsigned iring=1;iring<500;iring++)
    for(unsigned ichan=0;ichan<6*iring;ichan++)
    {
      double x1,y1;
      float x2,y2;
      hexid_to_xy(hexid, x1, y1);
      test_avx2_hexid_to_xy_f(hexid, x2, y2);
      EXPECT_NEAR(float(x1),x2,std::abs(x2)*1e-6);
      EXPECT_NEAR(float(y1),y2,std::abs(y2)*1e-6);
      hexid++;
    }
}

TEST(TestHexArray, AVX2_HexIDToXY_CW_Equals_Scalar) {
  unsigned hexid = 1;
  for(unsigned iring=1;iring<500;iring++)
    for(unsigned ichan=0;ichan<6*iring;ichan++)
    {
      double x1,y1;
      float x2,y2;
      hexid_to_xy(hexid, x1, y1, true);
      test_avx2_hexid_to_xy_f(hexid, x2, y2, true);
      EXPECT_NEAR(float(x1),x2,std::abs(x2)*1e-6);
      EXPECT_NEAR(float(y1),y2,std::abs(y2)*1e-6);
      hexid++;
    }
}
#endif
#endif

TEST(TestHexArray, XYToHexID_NewCodeSpeedTest) {
  double dx = 0.005;
  for(double x=-10.005; x<10.015; x+=dx)
    for(double y=-10.005; y<10.015; y+=dx)
    {
      double xx1 = x;
      double yy1 = y;
      xy_to_hexid_with_remainder(xx1, yy1, true);
    }
}

TEST(TestHexArray, XYToHexID_VVVCodeSpeedTest) {
  double dx = 0.005;
  for(double x=-10.005; x<10.015; x+=dx)
    for(double y=-10.005; y<10.015; y+=dx)
    {
      double xx2 = x;
      double yy2 = y;
      int hexid2;
      xy_to_nh(&xx2,&yy2,&hexid2);
    }
}

#if 0
#if defined(__AVX2__) and defined(__FMA__)
TEST(TestHexArray, XYToHexID_AVX2CodeSpeedTest) {
  volatile __m256i hexid;
  float dx = 0.005;
  for(float x=-10.005; x<10.015; x+=dx) {
    __m256 vx = _mm256_set1_ps(x);
    for(float y=-10.005; y<10.015; y+=8*dx)
    {
      __m256 vy = _mm256_set_ps(y+7*dx,y+6*dx,y+5*dx,y+4*dx,y+3*dx,y+2*dx,y+dx,y);
      hexid = avx2_xy_to_hexid_with_remainder_f(vx,vy);
    }
  }
  EXPECT_GE(int32_t(hexid[0]), 0);
}
#endif
#endif

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

#if 0
#if defined(__AVX2__) and defined(__FMA__)
TEST(TestHexArray, AVX2_XYToHexID_Equals_Scalar) {
  for(double x=-10.005; x<10.015; x+=0.02)
    for(double y=-10.005; y<10.015; y+=0.02)
    {
      double xx1 = x;
      double yy1 = y;
      unsigned hexid = xy_to_hexid_with_remainder(xx1, yy1);
      float xx2 = x;
      float yy2 = y;
      unsigned hexid2 = test_avx2_xy_to_hexid_with_remainder_f(xx2, yy2);
      EXPECT_EQ(hexid, hexid2);
      EXPECT_NEAR(xx1,xx2,1e-6);
      EXPECT_NEAR(yy1,yy2,1e-6);
    }
}

TEST(TestHexArray, AVX2_XYToHexID_CW_Equals_Scalar) {
  for(double x=-10.005; x<10.015; x+=0.02)
    for(double y=-10.005; y<10.015; y+=0.02)
    {
      double xx1 = x;
      double yy1 = y;
      unsigned hexid = xy_to_hexid_with_remainder(xx1, yy1, true);
      float xx2 = x;
      float yy2 = y;
      unsigned hexid2 = test_avx2_xy_to_hexid_with_remainder_f(xx2, yy2, true);
      EXPECT_EQ(hexid, hexid2);
      EXPECT_NEAR(xx1,xx2,1e-6);
      EXPECT_NEAR(yy1,yy2,1e-6);
    }
}

TEST(TestHexArray, AVX2_XYToHexID_Trans_Equals_Scalar) {
  float dx = 0.32426534147237f;
  float dy = 1.25432437427634f;
  float theta = 12.21752765/180.0*M_PI;
  float ctheta = std::cos(theta);
  float stheta = std::sin(theta);
  float scale = 0.9326575752f;
  for(double x=-10.005; x<10.015; x+=0.02)
    for(double y=-10.005; y<10.015; y+=0.02)
    {
      double xx1 = x;
      double yy1 = y;
      unsigned hexid = xy_trans_to_hexid_with_remainder(xx1, yy1, ctheta, stheta, scale, dx, dy);
      float xx2 = x;
      float yy2 = y;
      unsigned hexid2 = test_avx2_xy_trans_to_hexid_with_remainder_f(xx2, yy2, ctheta, stheta, scale, dx, dy);
      EXPECT_EQ(hexid, hexid2);
      EXPECT_NEAR(xx1,xx2,1e-5);
      EXPECT_NEAR(yy1,yy2,1e-5);
    }
}

TEST(TestHexArray, AVX2_XYToHexID_CW_Trans_Equals_Scalar) {
  float dx = 0.32426534147237f;
  float dy = 1.25432437427634f;
  float theta = 12.21752765/180.0*M_PI;
  float ctheta = std::cos(theta);
  float stheta = std::sin(theta);
  float scale = 0.9326575752f;
  for(double x=-10.005; x<10.015; x+=0.02)
    for(double y=-10.005; y<10.015; y+=0.02)
    {
      double xx1 = x;
      double yy1 = y;
      unsigned hexid = xy_trans_to_hexid_with_remainder(xx1, yy1, true, ctheta, stheta, scale, dx, dy);
      float xx2 = x;
      float yy2 = y;
      unsigned hexid2 = test_avx2_xy_trans_to_hexid_with_remainder_f(xx2, yy2, true, ctheta, stheta, scale, dx, dy);
      EXPECT_EQ(hexid, hexid2);
      EXPECT_NEAR(xx1,xx2,1e-5);
      EXPECT_NEAR(yy1,yy2,1e-5);
    }
}
#endif

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

#if 0
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
#endif

TEST(TestHexArray, RotateHexID_CCW) {
  ASSERT_EQ(rotate_ccw1_hexid_ccw(0), 0U);
  ASSERT_EQ(rotate_ccw1_hexid_ccw(1), 2U);
  ASSERT_EQ(rotate_ccw1_hexid_ccw(2), 3U);
  ASSERT_EQ(rotate_ccw1_hexid_ccw(3), 4U);
  ASSERT_EQ(rotate_ccw1_hexid_ccw(4), 5U);
  ASSERT_EQ(rotate_ccw1_hexid_ccw(5), 6U);
  ASSERT_EQ(rotate_ccw1_hexid_ccw(6), 1U);

  ASSERT_EQ(rotate_ccw2_hexid_ccw(0), 0U);
  ASSERT_EQ(rotate_ccw2_hexid_ccw(1), 3U);
  ASSERT_EQ(rotate_ccw2_hexid_ccw(2), 4U);
  ASSERT_EQ(rotate_ccw2_hexid_ccw(3), 5U);
  ASSERT_EQ(rotate_ccw2_hexid_ccw(4), 6U);
  ASSERT_EQ(rotate_ccw2_hexid_ccw(5), 1U);
  ASSERT_EQ(rotate_ccw2_hexid_ccw(6), 2U);

  ASSERT_EQ(rotate_ccw3_hexid_ccw(0), 0U);
  ASSERT_EQ(rotate_ccw3_hexid_ccw(1), 4U);
  ASSERT_EQ(rotate_ccw3_hexid_ccw(2), 5U);
  ASSERT_EQ(rotate_ccw3_hexid_ccw(3), 6U);
  ASSERT_EQ(rotate_ccw3_hexid_ccw(4), 1U);
  ASSERT_EQ(rotate_ccw3_hexid_ccw(5), 2U);
  ASSERT_EQ(rotate_ccw3_hexid_ccw(6), 3U);

  ASSERT_EQ(rotate_ccw4_hexid_ccw(0), 0U);
  ASSERT_EQ(rotate_ccw4_hexid_ccw(1), 5U);
  ASSERT_EQ(rotate_ccw4_hexid_ccw(2), 6U);
  ASSERT_EQ(rotate_ccw4_hexid_ccw(3), 1U);
  ASSERT_EQ(rotate_ccw4_hexid_ccw(4), 2U);
  ASSERT_EQ(rotate_ccw4_hexid_ccw(5), 3U);
  ASSERT_EQ(rotate_ccw4_hexid_ccw(6), 4U);

  ASSERT_EQ(rotate_ccw5_hexid_ccw(0), 0U);
  ASSERT_EQ(rotate_ccw5_hexid_ccw(2), 1U);
  ASSERT_EQ(rotate_ccw5_hexid_ccw(3), 2U);
  ASSERT_EQ(rotate_ccw5_hexid_ccw(4), 3U);
  ASSERT_EQ(rotate_ccw5_hexid_ccw(5), 4U);
  ASSERT_EQ(rotate_ccw5_hexid_ccw(6), 5U);
}

TEST(TestHexArray, RotateHexID_CW) {
  ASSERT_EQ(rotate_ccw5_hexid_cw(0), 0U);
  ASSERT_EQ(rotate_ccw5_hexid_cw(1), 2U);
  ASSERT_EQ(rotate_ccw5_hexid_cw(2), 3U);
  ASSERT_EQ(rotate_ccw5_hexid_cw(3), 4U);
  ASSERT_EQ(rotate_ccw5_hexid_cw(4), 5U);
  ASSERT_EQ(rotate_ccw5_hexid_cw(5), 6U);
  ASSERT_EQ(rotate_ccw5_hexid_cw(6), 1U);

  ASSERT_EQ(rotate_ccw4_hexid_cw(0), 0U);
  ASSERT_EQ(rotate_ccw4_hexid_cw(1), 3U);
  ASSERT_EQ(rotate_ccw4_hexid_cw(2), 4U);
  ASSERT_EQ(rotate_ccw4_hexid_cw(3), 5U);
  ASSERT_EQ(rotate_ccw4_hexid_cw(4), 6U);
  ASSERT_EQ(rotate_ccw4_hexid_cw(5), 1U);
  ASSERT_EQ(rotate_ccw4_hexid_cw(6), 2U);

  ASSERT_EQ(rotate_ccw3_hexid_cw(0), 0U);
  ASSERT_EQ(rotate_ccw3_hexid_cw(1), 4U);
  ASSERT_EQ(rotate_ccw3_hexid_cw(2), 5U);
  ASSERT_EQ(rotate_ccw3_hexid_cw(3), 6U);
  ASSERT_EQ(rotate_ccw3_hexid_cw(4), 1U);
  ASSERT_EQ(rotate_ccw3_hexid_cw(5), 2U);
  ASSERT_EQ(rotate_ccw3_hexid_cw(6), 3U);

  ASSERT_EQ(rotate_ccw2_hexid_cw(0), 0U);
  ASSERT_EQ(rotate_ccw2_hexid_cw(1), 5U);
  ASSERT_EQ(rotate_ccw2_hexid_cw(2), 6U);
  ASSERT_EQ(rotate_ccw2_hexid_cw(3), 1U);
  ASSERT_EQ(rotate_ccw2_hexid_cw(4), 2U);
  ASSERT_EQ(rotate_ccw2_hexid_cw(5), 3U);
  ASSERT_EQ(rotate_ccw2_hexid_cw(6), 4U);

  ASSERT_EQ(rotate_ccw1_hexid_cw(0), 0U);
  ASSERT_EQ(rotate_ccw1_hexid_cw(1), 6U);
  ASSERT_EQ(rotate_ccw1_hexid_cw(2), 1U);
  ASSERT_EQ(rotate_ccw1_hexid_cw(3), 2U);
  ASSERT_EQ(rotate_ccw1_hexid_cw(4), 3U);
  ASSERT_EQ(rotate_ccw1_hexid_cw(5), 4U);
  ASSERT_EQ(rotate_ccw1_hexid_cw(6), 5U);
}

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
  for(unsigned iloop=0; iloop<NLOOP_HEXID_TOFROM_UV; iloop++)
  for(unsigned hexid=1;hexid<100000;hexid++)
  {
    int u = 0;
    int v = 0;
    hexid_to_uv_ccw(hexid, u, v);
    ASSERT_EQ(hexid, uv_to_hexid_ccw(u, v))
      << hexid << ' ' << u << ' ' << v;
  }
}

TEST(TestHexArray, RandHexIDToFromUV_CCW_EQ) {
  for(unsigned iloop=0; iloop<NLOOP_HEXID_TOFROM_UV; iloop++) {
    std::mt19937 gen(15939); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> dis(0, 100000);
    for(unsigned idev=0;idev<100000;idev++)
    {
      unsigned hexid = dis(gen);
      int u = 0;
      int v = 0;
      hexid_to_uv_ccw(hexid, u, v);
      ASSERT_EQ(hexid, uv_to_hexid_ccw(u, v))
        << hexid << ' ' << u << ' ' << v;
    }
  }
}

#if 0
#if defined(__AVX2__) and defined(__FMA__)
TEST(TestHexArray, AVX2_HexIDToFromUV_CW_EQ) {
  for(unsigned hexid=1;hexid<100000;hexid++)
  {
    int u = 0;
    int v = 0;
    test_avx2_hexid_to_uv_cw(hexid, u, v);
    ASSERT_EQ(hexid, test_avx2_uv_to_hexid_cw(u, v))
      << hexid << ' ' << u << ' ' << v;
  }
}

TEST(TestHexArray, AVX2_HexIDToFromUV_CW_Equals_Scalar) {
  for(unsigned hexid=1;hexid<100000;hexid++)
  {
    int u = 0;
    int v = 0;
    hexid_to_uv_cw(hexid, u, v);
    ASSERT_EQ(hexid, test_avx2_uv_to_hexid_cw(u, v))
      << hexid << ' ' << u << ' ' << v;
  }
}

TEST(TestHexArray, AVX2_HexIDToFromUV_CCW_EQ) {
  for(unsigned iloop=0; iloop<NLOOP_HEXID_TOFROM_UV/8; iloop++)
  for(unsigned hexid=1;hexid<100000;hexid++)
  {
    int u = 0;
    int v = 0;
    test_avx2_hexid_to_uv_ccw(hexid, u, v);
    ASSERT_EQ(hexid, test_avx2_uv_to_hexid_ccw(u, v))
      << hexid << ' ' << u << ' ' << v;
  }
}

TEST(TestHexArray, AVX2_HexIDToFromUV_CCW_Equals_Scalar) {
  for(unsigned hexid=1;hexid<100000;hexid++)
  {
    int u = 0;
    int v = 0;
    hexid_to_uv_ccw(hexid, u, v);
    ASSERT_EQ(hexid, test_avx2_uv_to_hexid_ccw(u, v))
      << hexid << ' ' << u << ' ' << v;
  }
}

TEST(TestHexArray, AVX2_RandHexIDToFromUV_CCW_EQ) {
  for(unsigned iloop=0; iloop<NLOOP_HEXID_TOFROM_UV/8; iloop++) {
    std::mt19937 gen(15939); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> dis(0, 100000);
    for(unsigned idev=0;idev<100000;idev++)
    {
      unsigned hexid = dis(gen);
      int u = 0;
      int v = 0;
      test_avx2_hexid_to_uv_ccw(hexid, u, v);
      ASSERT_EQ(hexid, test_avx2_uv_to_hexid_ccw(u, v))
        << hexid << ' ' << u << ' ' << v;
    }
  }
}
#endif
#endif

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
