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

using namespace calin::util::vcl;

TEST(TestVCL, Print) {
  Vec8f v(1.234);
  std::cout << v << '\n';
  Vec8i vi(-101);
  std::cout << vi << '\n';
  std::cout << to_float(vi) << '\n';
  std::cout << to_float(vi) + 0.5 << '\n';
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
