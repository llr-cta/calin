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

#include <random>
#include <iostream>
#include <iomanip>
#include <gtest/gtest.h>
#include <vector>

#include <util/vcl.hpp>
#include <math/hex_array.hpp>
#include <math/hex_array_vcl.hpp>

using namespace calin::util::vcl;
using namespace calin::math::hex_array;

template<typename VCLArchitecture> class VCLHexArrayTest :
  public VCLArchitecture, public testing::Test
{
public:
  using typename VCLArchitecture::uint32_vt;
  using typename VCLArchitecture::int32_vt;
  using typename VCLArchitecture::uint64_vt;
  using typename VCLArchitecture::int64_vt;
  using typename VCLArchitecture::float_vt;
  using typename VCLArchitecture::double_vt;

  using typename VCLArchitecture::uint32_bvt;
  using typename VCLArchitecture::int32_bvt;
  using typename VCLArchitecture::uint64_bvt;
  using typename VCLArchitecture::int64_bvt;
  using typename VCLArchitecture::float_bvt;
  using typename VCLArchitecture::double_bvt;
};

using ArchTypes = ::testing::Types<VCL128Architecture, VCL256Architecture, VCL512Architecture>;

TYPED_TEST_CASE(VCLHexArrayTest, ArchTypes);

TYPED_TEST(VCLHexArrayTest, HexIDToRingIDLoop_50Rings) {
  typename TypeParam::int32_vt hexid = 1;
  for(unsigned iring=1;iring<50;iring++)
    for(unsigned ichan=0;ichan<6*iring;ichan++)
    {
      typename TypeParam::int32_vt ringid = VCLHexArray<TypeParam>::positive_hexid_to_ringid_loop(hexid);
      EXPECT_TRUE(horizontal_and(ringid == iring))
        << "With hexid=" << hexid << " giving " << ringid;
      hexid++;
    }
}

TYPED_TEST(VCLHexArrayTest, HexIDToRingIDRoot_5000Rings) {
  typename TypeParam::int32_vt hexid = 1;
  for(unsigned iring=1;iring<5000;iring++)
    for(unsigned ichan=0;ichan<6*iring;ichan++)
    {
      typename TypeParam::int32_vt ringid = VCLHexArray<TypeParam>::positive_hexid_to_ringid_root(hexid);
      if(!horizontal_and(ringid == iring))
      {
        std::cout << "Mismatch for hexid=" << hexid[0] << " (ring=" << iring <<
          ") giving " << ringid[0] << '\n';
        EXPECT_GE(iring,2000);
        goto finish_test;
      }
      hexid++;
    }
finish_test:
  EXPECT_TRUE(true);
}

TYPED_TEST(VCLHexArrayTest, HexIDToFromRingSegRun_APriori) {
  typename TypeParam::int32_vt hexid = 1;
  for(unsigned iring=1;iring<100;iring++)
    for(unsigned iseg=0;iseg<6;iseg++)
      for(unsigned irun=0;irun<iring;irun++)
      {
        ASSERT_TRUE(horizontal_and(VCLHexArray<TypeParam>::positive_ringid_segid_runid_to_hexid(iring, iseg, irun) == hexid))
          << iring << ' ' << iseg << ' ' << irun;
        typename TypeParam::int32_vt ringid = 0;
        typename TypeParam::int32_vt segid = 0;
        typename TypeParam::int32_vt runid = 0;
        VCLHexArray<TypeParam>::positive_hexid_to_ringid_segid_runid(hexid, ringid, segid, runid);
        ASSERT_TRUE(horizontal_and(ringid == iring)) << iring << ' ' << iseg << ' ' << irun;
        ASSERT_TRUE(horizontal_and(segid == iseg)) << iring << ' ' << iseg << ' ' << irun;
        ASSERT_TRUE(horizontal_and(runid == irun)) << iring << ' ' << iseg << ' ' << irun;
        hexid++;
      }
}

TYPED_TEST(VCLHexArrayTest, HexIDToFromRingSegRun_EQ) {
  for(unsigned ihex=1;ihex<100000;ihex++)
  {
    typename TypeParam::int32_vt ringid = 0;
    typename TypeParam::int32_vt segid = 0;
    typename TypeParam::int32_vt runid = 0;
    VCLHexArray<TypeParam>::positive_hexid_to_ringid_segid_runid(ihex, ringid, segid, runid);
    typename TypeParam::int32_vt hexid =
      VCLHexArray<TypeParam>::positive_ringid_segid_runid_to_hexid(ringid, segid, runid);
    ASSERT_TRUE(horizontal_and(hexid == ihex))
      << ringid << ' ' << segid << ' ' << runid;
  }
}


// TEST(VCLHexArrayTest, HexIDToXY_Equals_Scalar) {
//   unsigned hexid = 1;
//   for(unsigned iring=1;iring<500;iring++)
//     for(unsigned ichan=0;ichan<6*iring;ichan++)
//     {
//       double x1,y1;
//       float x2,y2;
//       hexid_to_xy(hexid, x1, y1);
//       test_avx2_hexid_to_xy_f(hexid, x2, y2);
//       EXPECT_NEAR(float(x1),x2,abs(x2)*1e-6);
//       EXPECT_NEAR(float(y1),y2,abs(y2)*1e-6);
//       hexid++;
//     }
// }


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
      EXPECT_NEAR(float(x1),x2,abs(x2)*1e-6);
      EXPECT_NEAR(float(y1),y2,abs(y2)*1e-6);
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
      EXPECT_NEAR(float(x1),x2,abs(x2)*1e-6);
      EXPECT_NEAR(float(y1),y2,abs(y2)*1e-6);
      hexid++;
    }
}
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
