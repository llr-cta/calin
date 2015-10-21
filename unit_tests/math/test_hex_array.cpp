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

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
