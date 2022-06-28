/*

   calin/unit_tests/math/test_geometry_vcl.cpp -- Stephen Fegan -- 2018-09-14

   Unit tests for VCL geometry class

   Copyright 2018, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <math/moments_calc.hpp>
#include <math/rng_vcl.hpp>
#include <math/geometry_vcl.hpp>
#include <math/vector3d_util.hpp>

using namespace calin::math::geometry;
using namespace calin::math::rng;
using namespace calin::util::vcl;
using namespace calin::math::moments_calc;

template<typename VCLReal> class VCLGeometryTest :
  public VCLReal, public testing::Test
{
public:
};

using RealTypes = ::testing::Types<VCL128FloatReal, VCL256FloatReal, VCL512FloatReal, VCL128DoubleReal, VCL256DoubleReal, VCL512DoubleReal>;

TYPED_TEST_CASE(VCLGeometryTest, RealTypes);


TYPED_TEST(VCLGeometryTest, RotationMatrix_Rzy) {
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

TYPED_TEST(VCLGeometryTest, RotationMatrix_Rzyz) {
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

TYPED_TEST(VCLGeometryTest, RotationMatrix_Ryxy) {
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

TYPED_TEST(VCLGeometryTest, RotationMatrix_Rxzx) {
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

TYPED_TEST(VCLGeometryTest, RotateInPlace_Rzy) {
  VCLRNG<typename TypeParam::architecture> rng;
  for(unsigned i=0; i<10000; i++) {
    typename TypeParam::vec3_vt u;
    rng.uniform_on_unit_sphere_vec_real(u);

    typename TypeParam::vec3_vt vx = TypeParam::vec3_vt::UnitX();
    typename TypeParam::vec3_vt vy = TypeParam::vec3_vt::UnitY();
    typename TypeParam::vec3_vt vz = TypeParam::vec3_vt::UnitZ();
    calin::math::geometry::VCL<TypeParam>::rotate_in_place_z_to_u_Rzy(vx,u);
    calin::math::geometry::VCL<TypeParam>::rotate_in_place_z_to_u_Rzy(vy,u);
    calin::math::geometry::VCL<TypeParam>::rotate_in_place_z_to_u_Rzy(vz,u);
    typename TypeParam::vec3_vt vx_x_vy = vx.cross(vy);
    typename TypeParam::vec3_vt vax = vz.cross(u);
    typename TypeParam::vec3_vt vax_rot = vax;
    calin::math::geometry::VCL<TypeParam>::rotate_in_place_z_to_u_Rzy(vax_rot,u);

    ASSERT_TRUE(horizontal_and(abs(vz.x() - u.x())<1e-6)) << vz.x() << '\n' << u.x();
    ASSERT_TRUE(horizontal_and(abs(vz.y() - u.y())<1e-6)) << vz.y() << '\n' << u.y();
    ASSERT_TRUE(horizontal_and(abs(vz.z() - u.z())<1e-6)) << vz.y() << '\n' << u.z();
    ASSERT_TRUE(horizontal_and(abs(vx.dot(vy))<1e-6)) << vx.dot(vy);
    ASSERT_TRUE(horizontal_and(abs(vx.dot(vz))<1e-6)) << vx.dot(vy);
    ASSERT_TRUE(horizontal_and(abs(vy.dot(vz))<1e-6)) << vx.dot(vy);
    ASSERT_TRUE(horizontal_and(abs(vx_x_vy.x() - vz.x())<1e-6)) << vx_x_vy.x() << '\n' << vz.x();
    ASSERT_TRUE(horizontal_and(abs(vx_x_vy.y() - vz.y())<1e-6)) << vx_x_vy.y() << '\n' << vz.y();
    ASSERT_TRUE(horizontal_and(abs(vx_x_vy.z() - vz.z())<1e-6)) << vx_x_vy.y() << '\n' << vz.z();

    ASSERT_TRUE(horizontal_and(abs(vax.x() - vax_rot.x())<1e-6)) << vax.x() << '\n' << vax_rot.x();
    ASSERT_TRUE(horizontal_and(abs(vax.y() - vax_rot.y())<1e-6)) << vax.y() << '\n' << vax_rot.y();
    ASSERT_TRUE(horizontal_and(abs(vax.z() - vax_rot.z())<1e-6)) << vax.y() << '\n' << vax_rot.z();

    calin::math::geometry::VCL<TypeParam>::derotate_in_place_z_to_u_Rzy(vx,u);
    calin::math::geometry::VCL<TypeParam>::derotate_in_place_z_to_u_Rzy(vy,u);
    calin::math::geometry::VCL<TypeParam>::derotate_in_place_z_to_u_Rzy(vz,u);

    ASSERT_TRUE(horizontal_and(abs(vx.x() - 1.0)<1e-6)) << vz.x() << '\n' << u.x();
    ASSERT_TRUE(horizontal_and(abs(vx.y())<1e-6)) << vz.y() << '\n' << u.y();
    ASSERT_TRUE(horizontal_and(abs(vx.z())<1e-6)) << vz.y() << '\n' << u.z();
    ASSERT_TRUE(horizontal_and(abs(vy.x())<1e-6)) << vz.x() << '\n' << u.x();
    ASSERT_TRUE(horizontal_and(abs(vy.y() - 1.0)<1e-6)) << vz.y() << '\n' << u.y();
    ASSERT_TRUE(horizontal_and(abs(vy.z())<1e-6)) << vz.y() << '\n' << u.z();
    ASSERT_TRUE(horizontal_and(abs(vz.x())<1e-6)) << vz.x() << '\n' << u.x();
    ASSERT_TRUE(horizontal_and(abs(vz.y())<1e-6)) << vz.y() << '\n' << u.y();
    ASSERT_TRUE(horizontal_and(abs(vz.z() - 1.0)<1e-6)) << vz.y() << '\n' << u.z();
  }
}

TYPED_TEST(VCLGeometryTest, RotateInPlace_Rzyz) {
  VCLRNG<typename TypeParam::architecture> rng;
  for(unsigned i=0; i<10000; i++) {
    typename TypeParam::vec3_vt u;
    rng.uniform_on_unit_sphere_vec_real(u);

    typename TypeParam::vec3_vt vx = TypeParam::vec3_vt::UnitX();
    typename TypeParam::vec3_vt vy = TypeParam::vec3_vt::UnitY();
    typename TypeParam::vec3_vt vz = TypeParam::vec3_vt::UnitZ();
    calin::math::geometry::VCL<TypeParam>::rotate_in_place_z_to_u_Rzyz(vx,u);
    calin::math::geometry::VCL<TypeParam>::rotate_in_place_z_to_u_Rzyz(vy,u);
    calin::math::geometry::VCL<TypeParam>::rotate_in_place_z_to_u_Rzyz(vz,u);
    typename TypeParam::vec3_vt vx_x_vy = vx.cross(vy);
    typename TypeParam::vec3_vt vax = vz.cross(u);
    typename TypeParam::vec3_vt vax_rot = vax;
    calin::math::geometry::VCL<TypeParam>::rotate_in_place_z_to_u_Rzy(vax_rot,u);

    ASSERT_TRUE(horizontal_and(abs(vz.x() - u.x())<1e-6)) << vz.x() << '\n' << u.x();
    ASSERT_TRUE(horizontal_and(abs(vz.y() - u.y())<1e-6)) << vz.y() << '\n' << u.y();
    ASSERT_TRUE(horizontal_and(abs(vz.z() - u.z())<1e-6)) << vz.y() << '\n' << u.z();
    ASSERT_TRUE(horizontal_and(abs(vx.dot(vy))<1e-6)) << vx.dot(vy);
    ASSERT_TRUE(horizontal_and(abs(vx.dot(vz))<1e-6)) << vx.dot(vy);
    ASSERT_TRUE(horizontal_and(abs(vy.dot(vz))<1e-6)) << vx.dot(vy);
    ASSERT_TRUE(horizontal_and(abs(vx_x_vy.x() - vz.x())<1e-6)) << vx_x_vy.x() << '\n' << vz.x();
    ASSERT_TRUE(horizontal_and(abs(vx_x_vy.y() - vz.y())<1e-6)) << vx_x_vy.y() << '\n' << vz.y();
    ASSERT_TRUE(horizontal_and(abs(vx_x_vy.z() - vz.z())<1e-6)) << vx_x_vy.y() << '\n' << vz.z();

    ASSERT_TRUE(horizontal_and(abs(vax.x() - vax_rot.x())<1e-6)) << vax.x() << '\n' << vax_rot.x();
    ASSERT_TRUE(horizontal_and(abs(vax.y() - vax_rot.y())<1e-6)) << vax.y() << '\n' << vax_rot.y();
    ASSERT_TRUE(horizontal_and(abs(vax.z() - vax_rot.z())<1e-6)) << vax.y() << '\n' << vax_rot.z();

    calin::math::geometry::VCL<TypeParam>::derotate_in_place_z_to_u_Rzyz(vx,u);
    calin::math::geometry::VCL<TypeParam>::derotate_in_place_z_to_u_Rzyz(vy,u);
    calin::math::geometry::VCL<TypeParam>::derotate_in_place_z_to_u_Rzyz(vz,u);

    ASSERT_TRUE(horizontal_and(abs(vx.x() - 1.0)<1e-6)) << vz.x() << '\n' << u.x();
    ASSERT_TRUE(horizontal_and(abs(vx.y())<1e-6)) << vz.y() << '\n' << u.y();
    ASSERT_TRUE(horizontal_and(abs(vx.z())<1e-6)) << vz.y() << '\n' << u.z();
    ASSERT_TRUE(horizontal_and(abs(vy.x())<1e-6)) << vz.x() << '\n' << u.x();
    ASSERT_TRUE(horizontal_and(abs(vy.y() - 1.0)<1e-6)) << vz.y() << '\n' << u.y();
    ASSERT_TRUE(horizontal_and(abs(vy.z())<1e-6)) << vz.y() << '\n' << u.z();
    ASSERT_TRUE(horizontal_and(abs(vz.x())<1e-6)) << vz.x() << '\n' << u.x();
    ASSERT_TRUE(horizontal_and(abs(vz.y())<1e-6)) << vz.y() << '\n' << u.y();
    ASSERT_TRUE(horizontal_and(abs(vz.z() - 1.0)<1e-6)) << vz.y() << '\n' << u.z();
  }
}

TYPED_TEST(VCLGeometryTest, RotateInPlace_Ryxy) {
  VCLRNG<typename TypeParam::architecture> rng;
  for(unsigned i=0; i<10000; i++) {
    typename TypeParam::vec3_vt u;
    rng.uniform_on_unit_sphere_vec_real(u);

    typename TypeParam::vec3_vt vx = TypeParam::vec3_vt::UnitX();
    typename TypeParam::vec3_vt vy = TypeParam::vec3_vt::UnitY();
    typename TypeParam::vec3_vt vz = TypeParam::vec3_vt::UnitZ();
    calin::math::geometry::VCL<TypeParam>::rotate_in_place_y_to_u_Ryxy(vx,u);
    calin::math::geometry::VCL<TypeParam>::rotate_in_place_y_to_u_Ryxy(vy,u);
    calin::math::geometry::VCL<TypeParam>::rotate_in_place_y_to_u_Ryxy(vz,u);
    typename TypeParam::vec3_vt vx_x_vy = vx.cross(vy);
    typename TypeParam::vec3_vt vax = vy.cross(u);
    typename TypeParam::vec3_vt vax_rot = vax;
    calin::math::geometry::VCL<TypeParam>::rotate_in_place_y_to_u_Ryxy(vax_rot,u);

    ASSERT_TRUE(horizontal_and(abs(vy.x() - u.x())<1e-6)) << vy.x() << '\n' << u.x();
    ASSERT_TRUE(horizontal_and(abs(vy.y() - u.y())<1e-6)) << vy.y() << '\n' << u.y();
    ASSERT_TRUE(horizontal_and(abs(vy.z() - u.z())<1e-6)) << vy.y() << '\n' << u.z();
    ASSERT_TRUE(horizontal_and(abs(vx.dot(vy))<1e-6)) << vx.dot(vy);
    ASSERT_TRUE(horizontal_and(abs(vx.dot(vz))<1e-6)) << vx.dot(vy);
    ASSERT_TRUE(horizontal_and(abs(vy.dot(vz))<1e-6)) << vx.dot(vy);
    ASSERT_TRUE(horizontal_and(abs(vx_x_vy.x() - vz.x())<1e-6)) << vx_x_vy.x() << '\n' << vz.x();
    ASSERT_TRUE(horizontal_and(abs(vx_x_vy.y() - vz.y())<1e-6)) << vx_x_vy.y() << '\n' << vz.y();
    ASSERT_TRUE(horizontal_and(abs(vx_x_vy.z() - vz.z())<1e-6)) << vx_x_vy.y() << '\n' << vz.z();

    ASSERT_TRUE(horizontal_and(abs(vax.x() - vax_rot.x())<1e-6)) << vax.x() << '\n' << vax_rot.x();
    ASSERT_TRUE(horizontal_and(abs(vax.y() - vax_rot.y())<1e-6)) << vax.y() << '\n' << vax_rot.y();
    ASSERT_TRUE(horizontal_and(abs(vax.z() - vax_rot.z())<1e-6)) << vax.y() << '\n' << vax_rot.z();

    calin::math::geometry::VCL<TypeParam>::derotate_in_place_y_to_u_Ryxy(vx,u);
    calin::math::geometry::VCL<TypeParam>::derotate_in_place_y_to_u_Ryxy(vy,u);
    calin::math::geometry::VCL<TypeParam>::derotate_in_place_y_to_u_Ryxy(vz,u);

    ASSERT_TRUE(horizontal_and(abs(vx.x() - 1.0)<1e-6)) << vz.x() << '\n' << u.x();
    ASSERT_TRUE(horizontal_and(abs(vx.y())<1e-6)) << vz.y() << '\n' << u.y();
    ASSERT_TRUE(horizontal_and(abs(vx.z())<1e-6)) << vz.y() << '\n' << u.z();
    ASSERT_TRUE(horizontal_and(abs(vy.x())<1e-6)) << vz.x() << '\n' << u.x();
    ASSERT_TRUE(horizontal_and(abs(vy.y() - 1.0)<1e-6)) << vz.y() << '\n' << u.y();
    ASSERT_TRUE(horizontal_and(abs(vy.z())<1e-6)) << vz.y() << '\n' << u.z();
    ASSERT_TRUE(horizontal_and(abs(vz.x())<1e-6)) << vz.x() << '\n' << u.x();
    ASSERT_TRUE(horizontal_and(abs(vz.y())<1e-6)) << vz.y() << '\n' << u.y();
    ASSERT_TRUE(horizontal_and(abs(vz.z() - 1.0)<1e-6)) << vz.y() << '\n' << u.z();
  }
}

TYPED_TEST(VCLGeometryTest, RotateInPlace_Rxzx) {
  VCLRNG<typename TypeParam::architecture> rng;
  for(unsigned i=0; i<10000; i++) {
    typename TypeParam::vec3_vt u;
    rng.uniform_on_unit_sphere_vec_real(u);

    typename TypeParam::vec3_vt vx = TypeParam::vec3_vt::UnitX();
    typename TypeParam::vec3_vt vy = TypeParam::vec3_vt::UnitY();
    typename TypeParam::vec3_vt vz = TypeParam::vec3_vt::UnitZ();
    calin::math::geometry::VCL<TypeParam>::rotate_in_place_x_to_u_Rxzx(vx,u);
    calin::math::geometry::VCL<TypeParam>::rotate_in_place_x_to_u_Rxzx(vy,u);
    calin::math::geometry::VCL<TypeParam>::rotate_in_place_x_to_u_Rxzx(vz,u);
    typename TypeParam::vec3_vt vx_x_vy = vx.cross(vy);
    typename TypeParam::vec3_vt vax = vx.cross(u);
    typename TypeParam::vec3_vt vax_rot = vax;
    calin::math::geometry::VCL<TypeParam>::rotate_in_place_z_to_u_Rzy(vax_rot,u);

    ASSERT_TRUE(horizontal_and(abs(vx.x() - u.x())<1e-6)) << vx.x() << '\n' << u.x();
    ASSERT_TRUE(horizontal_and(abs(vx.y() - u.y())<1e-6)) << vx.y() << '\n' << u.y();
    ASSERT_TRUE(horizontal_and(abs(vx.z() - u.z())<1e-6)) << vx.y() << '\n' << u.z();
    ASSERT_TRUE(horizontal_and(abs(vx.dot(vy))<1e-6)) << vx.dot(vy);
    ASSERT_TRUE(horizontal_and(abs(vx.dot(vz))<1e-6)) << vx.dot(vy);
    ASSERT_TRUE(horizontal_and(abs(vy.dot(vz))<1e-6)) << vx.dot(vy);
    ASSERT_TRUE(horizontal_and(abs(vx_x_vy.x() - vz.x())<1e-6)) << vx_x_vy.x() << '\n' << vz.x();
    ASSERT_TRUE(horizontal_and(abs(vx_x_vy.y() - vz.y())<1e-6)) << vx_x_vy.y() << '\n' << vz.y();
    ASSERT_TRUE(horizontal_and(abs(vx_x_vy.z() - vz.z())<1e-6)) << vx_x_vy.y() << '\n' << vz.z();

    ASSERT_TRUE(horizontal_and(abs(vax.x() - vax_rot.x())<1e-6)) << vax.x() << '\n' << vax_rot.x();
    ASSERT_TRUE(horizontal_and(abs(vax.y() - vax_rot.y())<1e-6)) << vax.y() << '\n' << vax_rot.y();
    ASSERT_TRUE(horizontal_and(abs(vax.z() - vax_rot.z())<1e-6)) << vax.y() << '\n' << vax_rot.z();

    calin::math::geometry::VCL<TypeParam>::derotate_in_place_x_to_u_Rxzx(vx,u);
    calin::math::geometry::VCL<TypeParam>::derotate_in_place_x_to_u_Rxzx(vy,u);
    calin::math::geometry::VCL<TypeParam>::derotate_in_place_x_to_u_Rxzx(vz,u);

    ASSERT_TRUE(horizontal_and(abs(vx.x() - 1.0)<1e-6)) << vz.x() << '\n' << u.x();
    ASSERT_TRUE(horizontal_and(abs(vx.y())<1e-6)) << vz.y() << '\n' << u.y();
    ASSERT_TRUE(horizontal_and(abs(vx.z())<1e-6)) << vz.y() << '\n' << u.z();
    ASSERT_TRUE(horizontal_and(abs(vy.x())<1e-6)) << vz.x() << '\n' << u.x();
    ASSERT_TRUE(horizontal_and(abs(vy.y() - 1.0)<1e-6)) << vz.y() << '\n' << u.y();
    ASSERT_TRUE(horizontal_and(abs(vy.z())<1e-6)) << vz.y() << '\n' << u.z();
    ASSERT_TRUE(horizontal_and(abs(vz.x())<1e-6)) << vz.x() << '\n' << u.x();
    ASSERT_TRUE(horizontal_and(abs(vz.y())<1e-6)) << vz.y() << '\n' << u.y();
    ASSERT_TRUE(horizontal_and(abs(vz.z() - 1.0)<1e-6)) << vz.y() << '\n' << u.z();
  }
}

TEST(VCLGeometryTest, ScatterDirectionEqualsScalar) {
  VCLRealRNG<VCL256DoubleReal> v_rng;
  RNG s_rng;

  SecondMomentsCalc2D v_moments;
  SecondMomentsCalc2D s_moments;

  for(unsigned i=0; i<1000000; i++) {
    VCL256DoubleReal::vec3_vt v_vec = VCL256DoubleReal::vec3_vt::UnitZ();
    Eigen::Vector3d s_vec = Eigen::Vector3d::UnitZ();

    calin::math::geometry::VCL<VCL256DoubleReal>::scatter_direction_in_place(v_vec, 0.01, v_rng);
    calin::math::vector3d_util::scatter_direction(s_vec, 0.01, s_rng);
    v_moments.accumulate(v_vec.x()[0], v_vec.y()[0]);
    s_moments.accumulate(s_vec.x(), s_vec.y());
  }
  EXPECT_NEAR(v_moments.mean_x(), v_moments.mean_x(), 0.0001);
  EXPECT_NEAR(v_moments.mean_y(), v_moments.mean_y(), 0.0001);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
