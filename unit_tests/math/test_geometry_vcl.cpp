/*

   calin/unit_tests/math/test_geometry_vcl.cpp -- Stephen Fegan -- 2018-09-14

   Unit tests for VCL geometry class

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

#include "math/rng_vcl.hpp"
#include "math/geometry_vcl.hpp"

using namespace calin::math::geometry;
using namespace calin::math::rng;
using namespace calin::util::vcl;

template<typename VCLRealType> class VCLGeometryTest :
  public VCLRealType, public testing::Test
{
public:
};

using RealTypes = ::testing::Types<VCL128FloatReal, VCL256FloatReal, VCL512FloatReal, VCL128DoubleReal, VCL256DoubleReal, VCL512DoubleReal>;

TYPED_TEST_CASE(VCLGeometryTest, RealTypes);


TYPED_TEST(VCLGeometryTest, Rotation_Rzy) {
  VCLRNG<typename TypeParam::architecture> rng;
  for(unsigned i=0; i<10000; i++) {
    typename TypeParam::vec3_vt dir;
    rng.uniform_on_unit_sphere_vec_real(dir);

    typename TypeParam::mat3_vt m;
    calin::math::geometry::VCL<TypeParam>::rotation_z_to_xyz_Rzy(m, dir.x(), dir.y(), dir.z());

    typename TypeParam::vec3_vt dir2 = m * TypeParam::vec3_vt::UnitZ();

    ASSERT_TRUE(horizontal_and(dir.x() == dir2.x())) << dir.x() << dir2.x();
    ASSERT_TRUE(horizontal_and(dir.y() == dir2.y())) << dir.y() << dir2.y();
    ASSERT_TRUE(horizontal_and(dir.z() == dir2.z())) << dir.y() << dir2.z();
    ASSERT_TRUE(horizontal_and(abs(m.determinant()-1)<1e-6)) << m.determinant();
  }
}

TYPED_TEST(VCLGeometryTest, Rotation_Rzyz) {
  VCLRNG<typename TypeParam::architecture> rng;
  for(unsigned i=0; i<10000; i++) {
    typename TypeParam::vec3_vt dir;
    rng.uniform_on_unit_sphere_vec_real(dir);

    typename TypeParam::mat3_vt m;
    calin::math::geometry::VCL<TypeParam>::rotation_z_to_xyz_Rzyz(m, dir.x(), dir.y(), dir.z());

    typename TypeParam::vec3_vt dir2 = m * TypeParam::vec3_vt::UnitZ();

    // Test axis is eigen-vector
    typename TypeParam::vec3_vt ax = TypeParam::vec3_vt::UnitZ().cross(dir);
    typename TypeParam::vec3_vt ax2 = m * ax;

    ASSERT_TRUE(horizontal_and(dir.x() == dir2.x())) << dir.x() << dir2.x();
    ASSERT_TRUE(horizontal_and(dir.y() == dir2.y())) << dir.y() << dir2.y();
    ASSERT_TRUE(horizontal_and(dir.z() == dir2.z())) << dir.y() << dir2.z();
    ASSERT_TRUE(horizontal_and(abs(m.determinant()-1)<1e-6)) << m.determinant();
    ASSERT_TRUE(horizontal_and(abs(ax.x() - ax2.x())<1e-6)) << ax.x() << ax2.x();
    ASSERT_TRUE(horizontal_and(abs(ax.y() - ax2.y())<1e-6)) << ax.y() << ax2.y();
    ASSERT_TRUE(horizontal_and(abs(ax.z() - ax2.z())<1e-6)) << ax.y() << ax2.z();

  }
}

TYPED_TEST(VCLGeometryTest, Rotation_Ryxy) {
  VCLRNG<typename TypeParam::architecture> rng;
  for(unsigned i=0; i<10000; i++) {
    typename TypeParam::vec3_vt dir;
    rng.uniform_on_unit_sphere_vec_real(dir);

    typename TypeParam::mat3_vt m;
    calin::math::geometry::VCL<TypeParam>::rotation_y_to_xyz_Ryxy(m, dir.x(), dir.y(), dir.z());

    typename TypeParam::vec3_vt dir2 = m * TypeParam::vec3_vt::UnitY();

    // Test axis is eigen-vector
    typename TypeParam::vec3_vt ax = TypeParam::vec3_vt::UnitY().cross(dir);
    typename TypeParam::vec3_vt ax2 = m * ax;

    ASSERT_TRUE(horizontal_and(dir.x() == dir2.x())) << dir.x() << dir2.x();
    ASSERT_TRUE(horizontal_and(dir.y() == dir2.y())) << dir.y() << dir2.y();
    ASSERT_TRUE(horizontal_and(dir.z() == dir2.z())) << dir.y() << dir2.z();
    ASSERT_TRUE(horizontal_and(abs(m.determinant()-1)<1e-6)) << m.determinant();
    ASSERT_TRUE(horizontal_and(abs(ax.x() - ax2.x())<1e-6)) << ax.x() << ax2.x();
    ASSERT_TRUE(horizontal_and(abs(ax.y() - ax2.y())<1e-6)) << ax.y() << ax2.y();
    ASSERT_TRUE(horizontal_and(abs(ax.z() - ax2.z())<1e-6)) << ax.y() << ax2.z();

  }
}

TYPED_TEST(VCLGeometryTest, Rotation_Rxzx) {
  VCLRNG<typename TypeParam::architecture> rng;
  for(unsigned i=0; i<10000; i++) {
    typename TypeParam::vec3_vt dir;
    rng.uniform_on_unit_sphere_vec_real(dir);

    typename TypeParam::mat3_vt m;
    calin::math::geometry::VCL<TypeParam>::rotation_x_to_xyz_Rxzx(m, dir.x(), dir.y(), dir.z());

    typename TypeParam::vec3_vt dir2 = m * TypeParam::vec3_vt::UnitX();

    // Test axis is eigen-vector
    typename TypeParam::vec3_vt ax = TypeParam::vec3_vt::UnitX().cross(dir);
    typename TypeParam::vec3_vt ax2 = m * ax;

    ASSERT_TRUE(horizontal_and(dir.x() == dir2.x())) << dir.x() << dir2.x();
    ASSERT_TRUE(horizontal_and(dir.y() == dir2.y())) << dir.y() << dir2.y();
    ASSERT_TRUE(horizontal_and(dir.z() == dir2.z())) << dir.y() << dir2.z();
    ASSERT_TRUE(horizontal_and(abs(m.determinant()-1)<1e-6)) << m.determinant();
    ASSERT_TRUE(horizontal_and(abs(ax.x() - ax2.x())<1e-6)) << ax.x() << ax2.x();
    ASSERT_TRUE(horizontal_and(abs(ax.y() - ax2.y())<1e-6)) << ax.y() << ax2.y();
    ASSERT_TRUE(horizontal_and(abs(ax.z() - ax2.z())<1e-6)) << ax.y() << ax2.z();

  }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
