#include <iostream>
#include <iomanip>
#include <gtest/gtest.h>
#include <vector>

#include "math/hex_array.hpp"

using namespace calin::math::hex_array;
using namespace calin::math::hex_array::vvv;

TEST(TestHexArray, HexIDToRingIDLoop_SpeedTest50Rings) {
  for(unsigned iloop = 0; iloop<3000; iloop++)
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
  for(unsigned iloop = 0; iloop<3000; iloop++)
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

TEST(TestHexArray, HexIDToRingIDRoot_2000Rings) {
  unsigned hexid = 1;
  for(unsigned iring=1;iring<2000;iring++)
    for(unsigned ichan=0;ichan<6*iring;ichan++)
    {
      EXPECT_EQ(iring, positive_hexid_to_ringid_root(hexid));
      hexid++;
    }
}

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

TEST(TestHexArray, SomeClusterElements) {
  EXPECT_EQ(cluster_hexid_to_member_hexid(0,1),
            std::vector<unsigned>({0,1,2,3,4,5,6}));
  EXPECT_EQ(cluster_hexid_to_member_hexid(1,1),
            std::vector<unsigned>({20,38,39,21,8,7,19}));
  EXPECT_EQ(cluster_hexid_to_member_hexid(2,1),
            std::vector<unsigned>({23,22,42,43,24,10,9}));
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
