#include <iostream>
#include <iomanip>
#include <gtest/gtest.h>
#include <vector>

#include "math/hex_array.hpp"

using namespace calin::math::hex_array;
using namespace calin::math::hex_array::vvv;

TEST(TestHexArray, HexIDToRingIDLoop_50Rings) {
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

TEST(TestHexArray, HexIDToRingIDRoot_50Rings) {
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

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
