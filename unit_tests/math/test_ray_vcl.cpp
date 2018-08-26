/*

   calin/unit_tests/math/test_hex_array.cpp -- Stephen Fegan -- 2018-08-25

   Unit tests for VCL ray class

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

#include <random>
#include <iostream>
#include <iomanip>
#include <gtest/gtest.h>
#include <vector>

#include "math/ray_vcl.hpp"

using namespace calin::math::ray;
using namespace calin::util::vcl;

template<typename VCLArchitecture> class VCLRayTest :
  public VCLArchitecture, public testing::Test
{
public:
};

using RealTypes = ::testing::Types<VCL128FloatReal, VCL256FloatReal, VCL512FloatReal>;

TYPED_TEST_CASE(VCLRayTest, RealTypes);

TYPED_TEST(VCLRayTest, SimplePropagateToYPlane) {
  typename TypeParam::vec3_vt pos(0.0,-1.0,0);
  typename TypeParam::vec3_vt dir(0.0,std::sqrt(0.5f),std::sqrt(0.5f));
  VCLRay<TypeParam> ray(pos,dir);
  ASSERT_TRUE(horizontal_and(ray.propagate_to_y_plane(0,true,1.5)));
  ASSERT_TRUE(horizontal_and(ray.ux() == dir.x()));
  ASSERT_TRUE(horizontal_and(ray.uy() == dir.y()));
  ASSERT_TRUE(horizontal_and(ray.uz() == dir.z()));
  ASSERT_TRUE(horizontal_and(ray.x() == 0.0f)) << ray.x();
  ASSERT_TRUE(horizontal_and(abs(ray.y()) < 1e-7f)) << ray.y();
  ASSERT_TRUE(horizontal_and(abs(ray.z() - 1.0f) < 1e-7))
    << ray.z() << ' ' << ray.z()-1.0f;
  ASSERT_TRUE(horizontal_and(abs(ray.ct() - std::sqrt(2.0f)*1.5f) < 1e-7))
    << ray.ct();
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
