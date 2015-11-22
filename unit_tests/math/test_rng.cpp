/* 

   calin/unit_tests/math/test_rng.cpp -- Stephen Fegan -- 2015-11-22

   Unit tests for RNG

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, Ecole polytechnique, CNRS/IN2P3, Universite Paris-Saclay

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
#include <tuple>

#include <math/rng.hpp>

using namespace calin::math::rng;

template<typename T> std::tuple<double, double, double>
calc_moments(const T& generator, unsigned N = 10000)
{
  double sum_x = 0;
  double sum_xx = 0;
  double sum_xxx = 0;
  for(unsigned i=0; i<N; i++) {
    double x;
    generator(x);
    sum_x += x;
    sum_xx += x*x;
    sum_xxx += x*x*x;
  }
  double m1 = sum_x/double(N);
  double m2 = sum_xx/double(N) - m1*m1;
  double m3 = sum_xxx/double(N) - 3*sum_xx/double(N)*m1 + 2*m1*m1*m1;
  return std::make_tuple(m1,m2,m3);
}

TEST(TestRNG, Uniform) {
  RNG rng;
  auto res = calc_moments([&rng](double& x){ x=rng.uniform(); });
  EXPECT_NEAR(std::get<0>(res), 0.5, 0.01);
  EXPECT_NEAR(std::get<1>(res), 1.0/12.0, 0.01);
  EXPECT_NEAR(std::get<2>(res), 0.0, 0.001);
  std::cout << std::get<0>(res) << ' '
            << std::get<1>(res) << ' '
            << std::get<2>(res) << '\n';
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
