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

TYPED_TEST(VCLHexArrayTest, HexIDToFromUV_CW_Equals_Scalar) {
  for(unsigned ihex=1;ihex<100000;ihex++)
  {
    int u = 0;
    int v = 0;
    hexid_to_uv_cw(ihex, u, v);
    typename TypeParam::int32_vt hexid = uv_to_hexid_cw(u, v);
    ASSERT_TRUE(horizontal_and(hexid == ihex))
      << ihex << ' ' << hexid << ' ' << u << ' ' << v;
  }
}

TYPED_TEST(VCLHexArrayTest, HexIDToFromUV_CCW_Equals_Scalar) {
  for(unsigned ihex=1;ihex<100000;ihex++)
  {
    int u = 0;
    int v = 0;
    hexid_to_uv_ccw(ihex, u, v);
    typename TypeParam::int32_vt hexid = uv_to_hexid_ccw(u, v);
    ASSERT_TRUE(horizontal_and(hexid == ihex))
      << ihex << ' ' << hexid << ' ' << u << ' ' << v;
  }
}

TYPED_TEST(VCLHexArrayTest, HexIDToFromUV_CW_EQ) {
  for(unsigned ihex=1;ihex<100000;ihex++)
  {
    typename TypeParam::int32_vt u = 0;
    typename TypeParam::int32_vt v = 0;
    VCLHexArray<TypeParam>::hexid_to_uv_cw(ihex, u, v);
    typename TypeParam::int32_vt hexid = VCLHexArray<TypeParam>::uv_to_hexid_cw(u, v);
    ASSERT_TRUE(horizontal_and(hexid == ihex))
      << ihex << ' ' << u << ' ' << v << ' ' << hexid;
  }
}

TYPED_TEST(VCLHexArrayTest, HexIDToFromUV_CCW_EQ) {
  for(unsigned ihex=1;ihex<100000;ihex++)
  {
    typename TypeParam::int32_vt u = 0;
    typename TypeParam::int32_vt v = 0;
    VCLHexArray<TypeParam>::hexid_to_uv_ccw(ihex, u, v);
    typename TypeParam::int32_vt hexid = VCLHexArray<TypeParam>::uv_to_hexid_ccw(u, v);
    ASSERT_TRUE(horizontal_and(hexid == ihex))
      << ihex << ' ' << u << ' ' << v << ' ' << hexid;
  }
}

TYPED_TEST(VCLHexArrayTest, HexIDToXY_Equals_Scalar) {
  unsigned hexid = 1;
  for(unsigned iring=1;iring<500;iring++)
    for(unsigned ichan=0;ichan<6*iring;ichan++)
    {
      double x1,y1;
      typename TypeParam::float_vt x2,y2;
      hexid_to_xy(hexid, x1, y1);
      VCLHexArray<TypeParam>::hexid_to_xy_f(hexid, x2, y2);
      ASSERT_NEAR(float(x1),x2[0],abs(x1)*1e-6);
      ASSERT_NEAR(float(y1),y2[0],abs(y1)*1e-6);
      hexid++;
    }
}

TYPED_TEST(VCLHexArrayTest, HexIDToXY_CW_Equals_Scalar) {
  unsigned hexid = 1;
  for(unsigned iring=1;iring<500;iring++)
    for(unsigned ichan=0;ichan<6*iring;ichan++)
    {
      double x1,y1;
      typename TypeParam::float_vt x2,y2;
      hexid_to_xy(hexid, x1, y1, true);
      VCLHexArray<TypeParam>::hexid_to_xy_f(hexid, x2, y2, true);
      ASSERT_NEAR(float(x1),x2[0],abs(x1)*1e-6);
      ASSERT_NEAR(float(y1),y2[0],abs(y1)*1e-6);
      hexid++;
    }
}

TYPED_TEST(VCLHexArrayTest, XYToHexID_Equals_Scalar) {
  for(double x=-10.005; x<10.015; x+=0.02)
    for(double y=-10.005; y<10.015; y+=0.02)
    {
      double xx1 = x;
      double yy1 = y;
      unsigned hexid = xy_to_hexid_with_remainder(xx1, yy1);
      typename TypeParam::float_vt xx2 = x;
      typename TypeParam::float_vt yy2 = y;
      typename TypeParam::int32_vt hexid2 =
        VCLHexArray<TypeParam>::xy_to_hexid_with_remainder_f(xx2, yy2);
      ASSERT_EQ(hexid, hexid2[0]);
      ASSERT_NEAR(xx1,xx2[0],1e-6);
      ASSERT_NEAR(yy1,yy2[0],1e-6);
    }
}

TYPED_TEST(VCLHexArrayTest, XYToHexID_CW_Equals_Scalar) {
  for(double x=-10.005; x<10.015; x+=0.02)
    for(double y=-10.005; y<10.015; y+=0.02)
    {
      double xx1 = x;
      double yy1 = y;
      unsigned hexid = xy_to_hexid_with_remainder(xx1, yy1, true);
      typename TypeParam::float_vt xx2 = x;
      typename TypeParam::float_vt yy2 = y;
      typename TypeParam::int32_vt hexid2 =
        VCLHexArray<TypeParam>::xy_to_hexid_with_remainder_f(xx2, yy2, true);
      ASSERT_EQ(hexid, hexid2[0]);
      ASSERT_NEAR(xx1,xx2[0],1e-6);
      ASSERT_NEAR(yy1,yy2[0],1e-6);
    }
}

TYPED_TEST(VCLHexArrayTest, UVToXY_Trans_Equals_Scalar) {
  float dx = 0.32426534147237f;
  float dy = 1.25432437427634f;
  float theta = 12.21752765/180.0*M_PI;
  float ctheta = std::cos(theta);
  float stheta = std::sin(theta);
  float scale = 0.9326575752f;
  for(unsigned v=-10;v<=10;v++)
    for(unsigned u=-10;u<=10;u++)
    {
      double x1,y1;
      typename TypeParam::float_vt x2,y2;
      uv_to_xy_trans(u, v, x1, y1, ctheta, stheta, scale, dx, dy);
      VCLHexArray<TypeParam>::uv_to_xy_trans_f(u, v, x2, y2, ctheta, stheta, scale, dx, dy);
      ASSERT_NEAR(float(x1),x2[0],abs(x1)*1e-6)
        << u << ' ' << v << ' ' << x1 << ' ' << y1 << ' '
        << x2 << ' ' << y2;
      ASSERT_NEAR(float(y1),y2[0],abs(y1)*1e-6)
        << u << ' ' << v << ' ' << x1 << ' ' << y1 << ' '
        << x2 << ' ' << y2;
    }
}

TYPED_TEST(VCLHexArrayTest, XYToUV_Trans_Equals_Scalar) {
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
      int u, v;
      xy_trans_to_uv_with_remainder(xx1, yy1, u, v, ctheta, stheta, scale, dx, dy);
      typename TypeParam::float_vt xx2 = x;
      typename TypeParam::float_vt yy2 = y;
      typename TypeParam::int32_vt u2;
      typename TypeParam::int32_vt v2;
      VCLHexArray<TypeParam>::xy_trans_to_uv_with_remainder_f(xx2, yy2, u2, v2, ctheta, stheta, scale, dx, dy);
      ASSERT_EQ(u, u2[0]);
      ASSERT_EQ(v, v2[0]);
      ASSERT_NEAR(xx1,xx2[0],1e-5);
      ASSERT_NEAR(yy1,yy2[0],1e-5);
    }
}

TYPED_TEST(VCLHexArrayTest, HexIDToXY_Trans_Equals_Scalar) {
  float dx = 0.32426534147237f;
  float dy = 1.25432437427634f;
  float theta = 12.21752765/180.0*M_PI;
  float ctheta = std::cos(theta);
  float stheta = std::sin(theta);
  float scale = 0.9326575752f;
  unsigned hexid = 1;
  for(unsigned iring=1;iring<500;iring++)
    for(unsigned ichan=0;ichan<6*iring;ichan++)
    {
      double x1,y1;
      typename TypeParam::float_vt x2,y2;
      hexid_to_xy_trans(hexid, x1, y1, false, ctheta, stheta, scale, dx, dy);
      VCLHexArray<TypeParam>::hexid_to_xy_trans_ccw_f(hexid, x2, y2, ctheta, stheta, scale, dx, dy);
      ASSERT_NEAR(float(x1),x2[0],(1+abs(x1))*1e-5)
        << hexid << ' ' << iring << ' ' << ichan << ' ' << x1 << ' ' << y1 << ' '
        << x2 << ' ' << y2;
      ASSERT_NEAR(float(y1),y2[0],(1+abs(y1))*1e-5)
        << hexid << ' ' << iring << ' ' << ichan << ' ' << x1 << ' ' << y1 << ' '
        << x2 << ' ' << y2;
      hexid++;
    }
}


TYPED_TEST(VCLHexArrayTest, XYToHexID_Trans_Equals_Scalar) {
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
      unsigned hexid = xy_trans_to_hexid_with_remainder(xx1, yy1, false, ctheta, stheta, scale, dx, dy);
      typename TypeParam::float_vt xx2 = x;
      typename TypeParam::float_vt yy2 = y;
      typename TypeParam::int32_vt hexid2 =
        VCLHexArray<TypeParam>::xy_trans_to_hexid_with_remainder_ccw_f(xx2, yy2, ctheta, stheta, scale, dx, dy);
      ASSERT_EQ(hexid, hexid2[0]);
      ASSERT_NEAR(xx1,xx2[0],1e-5);
      ASSERT_NEAR(yy1,yy2[0],1e-5);
    }
}

TYPED_TEST(VCLHexArrayTest, XYToHexID_CW_Trans_Equals_Scalar) {
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
      typename TypeParam::float_vt xx2 = x;
      typename TypeParam::float_vt yy2 = y;
      typename TypeParam::int32_vt hexid2 =
        VCLHexArray<TypeParam>::xy_trans_to_hexid_with_remainder_cw_f(xx2, yy2, ctheta, stheta, scale, dx, dy);
      ASSERT_EQ(hexid, hexid2[0]);
      ASSERT_NEAR(xx1,xx2[0],1e-5);
      ASSERT_NEAR(yy1,yy2[0],1e-5);
    }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
