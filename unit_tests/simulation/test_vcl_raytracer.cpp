/*

   calin/unit_tests/io/test_vcl_raytracer.cpp -- Stephen Fegan -- 2018-09-17

   Unit tests for VCL raytracer

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

#include <iostream>
#include <iomanip>
#include <fstream>
#include <limits>
#include <gtest/gtest.h>

#include <math/rng.hpp>
#include <math/hex_array.hpp>
#include <math/moments_calc.hpp>
#include <simulation/vso_array.hpp>
#include <simulation/vso_raytracer.hpp>
#include <simulation/vs_optics.pb.h>

#include <util/vcl.hpp>
#include <math/ray_vcl.hpp>
#include <simulation/vcl_raytracer.hpp>

using namespace calin::util::vcl;
using namespace calin::simulation::vcl_raytracer;
using namespace calin::ix::simulation::vs_optics;

IsotropicDCArrayParameters* make_array_config(bool obscure_camera = false)
{
  IsotropicDCArrayParameters* mst = new IsotropicDCArrayParameters;
  mst->mutable_array_origin()->set_latitude(28.0 + 45.0/60.0 + 47.36/3600.0);
  mst->mutable_array_origin()->set_longitude(-(17.0 + 53.0/60.0 + 23.93/3600.0));
  mst->mutable_array_origin()->set_elevation(2147 * 100.0);

  auto* scope = mst->mutable_prescribed_array_layout()->add_scope_positions();
  scope->set_x(0);
  scope->set_y(0);
  scope->set_z(mst->array_origin().elevation());

  mst->mutable_reflector_frame()->set_optic_axis_rotation(-90.0);

  auto* dc = mst->mutable_reflector();
  dc->set_curvature_radius(1920);
  dc->set_aperture(1230);
  dc->set_facet_num_hex_rings(5);
  dc->mutable_psf_align()->set_object_plane(std::numeric_limits<double>::infinity());
  dc->set_alignment_image_plane(1600);
  dc->set_facet_spacing(122);
  dc->set_facet_size(120);
  dc->set_facet_focal_length(1607);
  dc->set_facet_focal_length_dispersion(1);
  dc->set_facet_spot_size_probability(0.8);
  dc->set_facet_spot_size(0.5 * 2.8); // Spot size of 28mm at 2F
  dc->set_facet_spot_size_dispersion(0.5 * 0.02);
  dc->set_facet_labeling_parity(true);
  dc->set_weathering_factor(1.0);
  for(auto id : { 67,72,82,87 })dc->add_facet_missing_list(id-1);
  mst->mutable_focal_plane()->set_camera_diameter(235);
  //mst->mutable_focal_plane()->mutable_translation()->set_y(1.0/(1.0/dc->alignment_image_plane()-1.0/(10 * 1e5)));
  mst->mutable_focal_plane()->mutable_translation()->set_y(dc->alignment_image_plane());
  mst->mutable_pixel()->set_spacing(5);
  mst->mutable_pixel()->set_cone_inner_diameter(5);
  mst->mutable_pixel()->set_cone_survival_prob(1);
  mst->mutable_pixel()->set_hex_module_size(1);
  mst->mutable_pixel()->set_hex_module_layout_use_b_configuration(true);
  mst->mutable_pixel()->set_module_num_hex_rings(9);

  int u1, v1;
  double x1,y1;
  calin::math::hex_array::cluster_hexid_to_center_uv(1,1,u1,v1,false);
  calin::math::hex_array::uv_to_xy(u1,v1,x1,y1);
  mst->mutable_pixel()->set_grid_rotation(atan2(-y1,x1)/M_PI*180.0 + 30.0);

  if(obscure_camera) {
    auto* obs_camera_box = mst->add_pre_reflection_obscuration();
    obs_camera_box->mutable_aligned_box()->mutable_max_corner()->set_x(2918.0/20);
    obs_camera_box->mutable_aligned_box()->mutable_max_corner()->set_z(2918.0/20);
    obs_camera_box->mutable_aligned_box()->mutable_max_corner()->set_y(mst->focal_plane().translation().y()+(1532.25-513.0)/10);
    obs_camera_box->mutable_aligned_box()->mutable_min_corner()->set_x(-2918.0/20);
    obs_camera_box->mutable_aligned_box()->mutable_min_corner()->set_z(-2918.0/20);
    obs_camera_box->mutable_aligned_box()->mutable_min_corner()->set_y(mst->focal_plane().translation().y()-513.0/10);

    auto* obs_outer_aperture = mst->add_post_reflection_obscuration();
    obs_outer_aperture->mutable_rectangular_aperture()->mutable_center_pos()->set_x(0);
    obs_outer_aperture->mutable_rectangular_aperture()->mutable_center_pos()->set_z(-50.0/20);
    obs_outer_aperture->mutable_rectangular_aperture()->mutable_center_pos()->set_y(mst->focal_plane().translation().y()-513.0/10);
    obs_outer_aperture->mutable_rectangular_aperture()->set_flat_to_flat_x(2714.0/10);
    obs_outer_aperture->mutable_rectangular_aperture()->set_flat_to_flat_z((2585.0-50.0)/10);

    auto* obs_inner_aperture = mst->add_post_reflection_obscuration();
    obs_inner_aperture->mutable_circular_aperture()->mutable_center_pos()->set_x(0);
    obs_inner_aperture->mutable_circular_aperture()->mutable_center_pos()->set_z(0);
    obs_inner_aperture->mutable_circular_aperture()->mutable_center_pos()->set_y(mst->focal_plane().translation().y()-222.5/10);
    obs_inner_aperture->mutable_circular_aperture()->set_diameter(2304.0/10);
  }

#if 0
  if include_window:
      win = mst->mutable_spherical_window()
      win.set_front_y_coord(mst->focal_plane().translation().y() - 427.5/10)
      win.set_outer_radius(3443.0/10)
      win.set_thickness(5.25/10.0)
      win.set_refractive_index(1.5)
#endif

  return mst;
}

calin::simulation::vs_optics::VSOArray* make_test_array(IsotropicDCArrayParameters* mst = nullptr)
{
  IsotropicDCArrayParameters* my_mst = nullptr;
  if(mst == nullptr)my_mst = mst = make_array_config();

  calin::simulation::vs_optics::VSOArray* array =
    new calin::simulation::vs_optics::VSOArray;
  calin::math::rng::RNG rng(__PRETTY_FUNCTION__, "VSOArray population");
  array->generateFromArrayParameters(*mst, rng);
  array->pointTelescopesAzEl(0,M_PI/2);

  delete my_mst;
  return array;
}

calin::simulation::vs_optics::VSOArray* make_test_array(bool obscure_camera)
{
  IsotropicDCArrayParameters* mst = make_array_config(obscure_camera);

  calin::simulation::vs_optics::VSOArray* array =
    new calin::simulation::vs_optics::VSOArray;
  calin::math::rng::RNG rng(__PRETTY_FUNCTION__, "VSOArray population");
  array->generateFromArrayParameters(*mst, rng);
  array->pointTelescopesAzEl(0,M_PI/2);

  delete mst;
  return array;
}


template<typename VCLReal> class TestVCLRaytracer :
  public VCLReal, public testing::Test
{
public:
};

using RealTypes = ::testing::Types<VCL128FloatReal, VCL256FloatReal, VCL512FloatReal, VCL128DoubleReal, VCL256DoubleReal, VCL512DoubleReal>;

TYPED_TEST_CASE(TestVCLRaytracer, RealTypes);

TYPED_TEST(TestVCLRaytracer, RayTrace) {
  calin::simulation::vs_optics::VSOArray* array = make_test_array();
  VCLScopeRayTracer<TypeParam> raytracer(array->telescope(0));

  for(const auto* mirror : array->telescope(0)->all_mirrors()) {
    calin::math::ray::VCLRay<TypeParam> ray;
    ray.mutable_direction() << 0.0 , -1.0, 0.0;
    ray.mutable_position() << mirror->pos().x() , 2000.0, mirror->pos().z();
    typename TypeParam::bool_vt mask = true;
    mask.insert(1, false);
    typename VCLScopeRayTracer<TypeParam>::TraceInfo trace_info;
    mask = raytracer.trace_reflector_frame(mask, ray, trace_info);
    ASSERT_TRUE(mask[0]);
    ASSERT_FALSE(mask[1]);
    ASSERT_EQ(trace_info.status[0], calin::simulation::vcl_raytracer::STS_TS_FOUND_PIXEL);
    ASSERT_EQ(trace_info.status[1], calin::simulation::vcl_raytracer::STS_MASKED_ON_ENTRY);
    ASSERT_EQ(trace_info.mirror_hexid[0], int(mirror->hexID()));
    ASSERT_EQ(trace_info.mirror_id[0], int(mirror->id()));
    ASSERT_EQ(trace_info.mirror_hexid[1], int(array->telescope(0)->numMirrorHexSites()));
    ASSERT_EQ(trace_info.mirror_id[1], int(array->telescope(0)->numMirrors()));
  }

  delete array;
}

TEST(TestVCLRaytracer, CompareToScalar) {
  calin::simulation::vs_optics::VSOArray* array = make_test_array();
  calin::simulation::vs_optics::VSOTelescope* scope = array->telescope(0);
  calin::math::rng::RNG s_rng;
  calin::simulation::vs_optics::VSORayTracer s_raytracer(array,&s_rng);

  VCLScopeRayTracer<VCL256FloatReal> raytracer(scope);
  calin::math::rng::VCLRealRNG<VCL256FloatReal> rng(__PRETTY_FUNCTION__,"Ray position generator");
  double theta = 0.0 / 180.0*M_PI;
  double cos_theta = std::cos(theta);
  double sin_theta = std::sin(theta);

  calin::math::moments_calc::SecondMomentsCalc2D v_moments;
  calin::math::moments_calc::SecondMomentsCalc2D s_moments;

  for(unsigned iray=0; iray<100000; iray++) {
    VCL256FloatReal::vec3_vt dir(0.0, -1.0, 0.0);
    VCL256FloatReal::vec3_vt pos(
      rng.uniform_zc(scope->reflectorIP()),
      2000.0,
      rng.uniform_zc(scope->reflectorIP()));
    calin::math::geometry::VCL<VCL256FloatReal>::rotate_in_place_Rz(
      pos, cos_theta, sin_theta);
    calin::math::geometry::VCL<VCL256FloatReal>::rotate_in_place_Rz(
      dir, cos_theta, sin_theta);

    calin::math::ray::VCLRay<VCL256FloatReal> ray(pos, dir);
    typename VCL256FloatReal::bool_vt mask = true;
    typename VCLScopeRayTracer<VCL256FloatReal>::TraceInfo trace_info;
    mask = raytracer.trace_reflector_frame(mask, ray, trace_info);
    for(unsigned i=0;i<8;i++) {
      if(trace_info.status[i]>=STS_OUTSIDE_FOCAL_PLANE_APERTURE)
        v_moments.accumulate(trace_info.fplane_x[i],trace_info.fplane_z[i]);
    }

    calin::math::geometry::VCL<VCL256FloatReal>::rotate_in_place_Rx(pos, 0.0, 1.0);
    calin::math::geometry::VCL<VCL256FloatReal>::rotate_in_place_Rx(dir, 0.0, 1.0);

    for(unsigned i=0;i<8;i++) {
      Eigen::Vector3d s_pos(pos.x()[i], pos.y()[i], pos.z()[i]);
      s_pos += scope->pos();
      Eigen::Vector3d s_dir(dir.x()[i], dir.y()[i], dir.z()[i]);
      calin::math::ray::Ray s_ray(s_pos, s_dir);
      calin::simulation::vs_optics::VSOTraceInfo s_trace_info;
      s_raytracer.trace(s_ray, s_trace_info, scope);
      if(s_trace_info.rayHitFocalPlane())
        s_moments.accumulate(s_trace_info.fplane_x,s_trace_info.fplane_z);
    }
  }

  std::cout
    << "S: " << s_moments.sum_w() << ' '
    << s_moments.mean_x() << ' ' << s_moments.mean_y() << ' '
    << s_moments.var_x() << ' ' << s_moments.var_y() << '\n';
  std::cout
    << "V: " << v_moments.sum_w() << ' '
    << v_moments.mean_x() << ' ' << v_moments.mean_y() << ' '
    << v_moments.var_x() << ' ' << v_moments.var_y() << '\n';

  EXPECT_NEAR(v_moments.sum_w(), s_moments.sum_w(), 2.5);
  EXPECT_NEAR(v_moments.mean_x(), s_moments.mean_x(), 4e-3);
  EXPECT_NEAR(v_moments.mean_y(), s_moments.mean_y(), 4e-3);
  EXPECT_NEAR(v_moments.var_x(), s_moments.var_x(), 3e-3);
  EXPECT_NEAR(v_moments.var_y(), s_moments.var_y(), 3e-3);
}

TEST(TestVCLRaytracer, CompareToScalar_Obscured) {
  calin::simulation::vs_optics::VSOArray* array = make_test_array(true);
  calin::simulation::vs_optics::VSOTelescope* scope = array->telescope(0);
  calin::math::rng::RNG s_rng;
  calin::simulation::vs_optics::VSORayTracer s_raytracer(array,&s_rng);

  VCLScopeRayTracer<VCL256FloatReal> raytracer(scope);
  calin::math::rng::VCLRealRNG<VCL256FloatReal> rng(__PRETTY_FUNCTION__,"Ray position generator");
  double theta = 0.0 / 180.0*M_PI;
  double cos_theta = std::cos(theta);
  double sin_theta = std::sin(theta);

  calin::math::moments_calc::SecondMomentsCalc2D v_moments;
  calin::math::moments_calc::SecondMomentsCalc2D s_moments;

  for(unsigned iray=0; iray<100000; iray++) {
    VCL256FloatReal::vec3_vt dir(0.0, -1.0, 0.0);
    VCL256FloatReal::vec3_vt pos(
      rng.uniform_zc(scope->reflectorIP()),
      2000.0,
      rng.uniform_zc(scope->reflectorIP()));
    calin::math::geometry::VCL<VCL256FloatReal>::rotate_in_place_Rz(
      pos, cos_theta, sin_theta);
    calin::math::geometry::VCL<VCL256FloatReal>::rotate_in_place_Rz(
      dir, cos_theta, sin_theta);

    calin::math::ray::VCLRay<VCL256FloatReal> ray(pos, dir);
    typename VCL256FloatReal::bool_vt mask = true;
    typename VCLScopeRayTracer<VCL256FloatReal>::TraceInfo trace_info;
    mask = raytracer.trace_reflector_frame(mask, ray, trace_info);
    for(unsigned i=0;i<8;i++) {
      if(trace_info.status[i]>=STS_OUTSIDE_FOCAL_PLANE_APERTURE)
        v_moments.accumulate(trace_info.fplane_x[i],trace_info.fplane_z[i]);
    }

    calin::math::geometry::VCL<VCL256FloatReal>::rotate_in_place_Rx(pos, 0.0, 1.0);
    calin::math::geometry::VCL<VCL256FloatReal>::rotate_in_place_Rx(dir, 0.0, 1.0);

    for(unsigned i=0;i<8;i++) {
      Eigen::Vector3d s_pos(pos.x()[i], pos.y()[i], pos.z()[i]);
      s_pos += scope->pos();
      Eigen::Vector3d s_dir(dir.x()[i], dir.y()[i], dir.z()[i]);
      calin::math::ray::Ray s_ray(s_pos, s_dir);
      calin::simulation::vs_optics::VSOTraceInfo s_trace_info;
      s_raytracer.trace(s_ray, s_trace_info, scope);
      if(s_trace_info.rayHitFocalPlane())
        s_moments.accumulate(s_trace_info.fplane_x,s_trace_info.fplane_z);
    }
  }

  std::cout
    << "S: " << s_moments.sum_w() << ' '
    << s_moments.mean_x() << ' ' << s_moments.mean_y() << ' '
    << s_moments.var_x() << ' ' << s_moments.var_y() << '\n';
  std::cout
    << "V: " << v_moments.sum_w() << ' '
    << v_moments.mean_x() << ' ' << v_moments.mean_y() << ' '
    << v_moments.var_x() << ' ' << v_moments.var_y() << '\n';

  EXPECT_NEAR(v_moments.sum_w(), s_moments.sum_w(), 2.5);
  EXPECT_NEAR(v_moments.mean_x(), s_moments.mean_x(), 4e-3);
  EXPECT_NEAR(v_moments.mean_y(), s_moments.mean_y(), 4e-3);
  EXPECT_NEAR(v_moments.var_x(), s_moments.var_x(), 3e-3);
  EXPECT_NEAR(v_moments.var_y(), s_moments.var_y(), 3e-3);
}

TEST(TestVCLRaytracer, LaserPSF) {
  calin::simulation::vs_optics::VSOArray* array = make_test_array();
  VCLScopeRayTracer<VCL256FloatReal> raytracer(array->telescope(0));
  calin::math::rng::VCLRealRNG<VCL256FloatReal> rng(__PRETTY_FUNCTION__,"Ray position generator");
  calin::math::moments_calc::SecondMomentsCalc2D v_moments;
  for(unsigned iray=0; iray<1000000; iray++) {
    calin::math::ray::VCLRay<VCL256FloatReal> ray;
    ray.mutable_direction() << 0.0, -1.0, 0.0;
    ray.mutable_position() << 0.0, 2000.0, 0.0;
    typename VCL256FloatReal::bool_vt mask = true;
    typename VCLScopeRayTracer<VCL256FloatReal>::TraceInfo trace_info;
    mask = raytracer.trace_reflector_frame(mask, ray, trace_info);
    ASSERT_TRUE(horizontal_and(mask)) << mask;
    for(unsigned i=0;i<8;i++)
      v_moments.accumulate(trace_info.fplane_x[i],trace_info.fplane_z[i]);
  }
  auto* scope = array->telescope(0);
  auto* mirror = scope->mirror(0);
  double disp_expected =
    mirror->spotSize()/mirror->focalLength()*scope->focalPlanePosition().y()*0.5;
  EXPECT_NEAR(std::sqrt(v_moments.var_x()), disp_expected, 0.001);
  EXPECT_NEAR(std::sqrt(v_moments.var_y()), disp_expected, 0.001);
  delete array;
}

TEST(TestVCLRaytracer, PSF) {
  calin::simulation::vs_optics::VSOArray* array = make_test_array(true);
  VCLScopeRayTracer<VCL256FloatReal> raytracer(array->telescope(0));
  calin::math::rng::VCLRealRNG<VCL256FloatReal> rng(__PRETTY_FUNCTION__,"Ray position generator");
  std::ofstream stream("test.dat");
  double theta = 2.0 / 180.0*M_PI;
  for(unsigned iray=0; iray<10000; iray++) {
    calin::math::ray::VCLRay<VCL256FloatReal> ray;
    ray.mutable_direction() << 0.0, -1.0, 0.0;
    ray.mutable_position() <<
      rng.uniform_zc(array->telescope(0)->reflectorIP()),
      2000.0,
      rng.uniform_zc(array->telescope(0)->reflectorIP());
    calin::math::geometry::VCL<VCL256FloatReal>::rotate_in_place_Rz(
      ray.mutable_direction(), cos(theta), sin(theta));
    calin::math::geometry::VCL<VCL256FloatReal>::rotate_in_place_Rz(
      ray.mutable_position(), cos(theta), sin(theta));
    typename VCL256FloatReal::bool_vt mask = true;
    typename VCLScopeRayTracer<VCL256FloatReal>::TraceInfo trace_info;
    auto pos = ray.position();
    mask = raytracer.trace_reflector_frame(mask, ray, trace_info);
    for(unsigned i=0;i<8;i++)stream
      << ' ' << mask[i] << ' ' << trace_info.status[i]
      << ' ' << pos.x()[i] << ' ' << pos.z()[i]
      << ' ' << trace_info.fplane_x[i] << ' ' << trace_info.fplane_z[i]
      << ' ' << trace_info.pixel_hexid[i] << ' ' << trace_info.pixel_id[i]
      << ' ' << trace_info.mirror_hexid[i] << ' ' << trace_info.mirror_id[i]
      << '\n';
  }

  delete array;
}

TEST(TestVCLAlignedCircularAperture, CompareToScalar) {
  using Real = VCL256DoubleReal;
  calin::math::rng::RNG rng(__PRETTY_FUNCTION__,"Ray generator");
  unsigned nray = 100000;
  unsigned nhit = 0;
  for(unsigned iray=0; iray<nray; ++iray) {
    double x0 = rng.uniform()*2.0 - 1.0;
    double y0 = rng.uniform()*2.0 - 1.0;
    double z0 = rng.uniform()*2.0 - 1.0;
    double d0 = rng.uniform();

    Eigen::Vector3d center(x0,y0,z0);
    calin::simulation::vs_optics::VSOAlignedCircularAperture sobs(center, d0);
    VCLAlignedCircularAperture<Real> vobs(center, d0);

    double x = rng.uniform()*2.0 - 1.0;
    double y = rng.uniform()*2.0 - 1.0;
    double z = rng.uniform()*2.0 - 1.0;
    double uy = rng.uniform()*2.0 - 1.0;
    double theta = rng.uniform()*2.0*M_PI;
    double ur = std::sqrt(1.0-uy*uy);
    double ux = ur*std::cos(theta);
    double uz = ur*std::sin(theta);
    Eigen::Vector3d r(x,y,z);
    Eigen::Vector3d u(ux,uy,uz);

    calin::math::ray::Ray sray_in(r,u);
    calin::math::ray::Ray sray_out;
    bool sintersects = sobs.doesObscure(sray_in, sray_out, 1.0);

    calin::math::ray::VCLRay<Real> vray_in(
      r.template cast<Real::real_vt>(), u.template cast<Real::real_vt>());
    calin::math::ray::VCLRay<Real> vray_out;
    Real::bool_vt vintersects = vobs.doesObscure(vray_in, vray_out, 1.0);

    ASSERT_EQ(sintersects, vintersects[0])
      << sintersects << ' ' << vintersects << ' ' << iray << '\n'
      << r.transpose() << '\n' << u.transpose();

    if(sintersects)++nhit;
  }
  ASSERT_GE(nhit*10, nray);
  ASSERT_LE(nhit*10, 9*nray);
}

TEST(TestAlignedBoxObscuration, CompareToScalar) {
  using Real = VCL256DoubleReal;
  calin::math::rng::RNG rng(__PRETTY_FUNCTION__,"Ray generator");
  unsigned nray = 100000;
  unsigned nhit = 0;
  for(unsigned iray=0; iray<nray; ++iray) {
    double min_x = -rng.uniform();
    double min_y = -rng.uniform();
    double min_z = -rng.uniform();
    double max_x = rng.uniform();
    double max_y = rng.uniform();
    double max_z = rng.uniform();

    Eigen::Vector3d min_xyz(min_x,min_y,min_z);
    Eigen::Vector3d max_xyz(max_x,max_y,max_z);

    calin::simulation::vs_optics::VSOAlignedBoxObscuration sobs(min_xyz, max_xyz);
    VCLAlignedBoxObscuration<Real> vobs(min_xyz, max_xyz);

    double x = rng.uniform()*2.0 - 1.0;
    double y = rng.uniform()*2.0 - 1.0;
    double z = rng.uniform()*2.0 - 1.0;
    double uy = rng.uniform()*2.0 - 1.0;
    double theta = rng.uniform()*2.0*M_PI;
    double ur = std::sqrt(1.0-uy*uy);
    double ux = ur*std::cos(theta);
    double uz = ur*std::sin(theta);
    Eigen::Vector3d r(x,y,z);
    Eigen::Vector3d u(ux,uy,uz);

    calin::math::ray::Ray sray_in(r,u);
    calin::math::ray::Ray sray_out;
    bool sintersects = sobs.doesObscure(sray_in, sray_out, 1.0);

    calin::math::ray::VCLRay<Real> vray_in(
      r.template cast<Real::real_vt>(), u.template cast<Real::real_vt>());
    calin::math::ray::VCLRay<Real> vray_out;
    Real::bool_vt vintersects = vobs.doesObscure(vray_in, vray_out, 1.0);

    ASSERT_EQ(sintersects, vintersects[0])
      << sintersects << ' ' << vintersects << ' ' << iray << '\n'
      << r.transpose() << '\n' << u.transpose();

    if(sintersects)++nhit;
  }
  ASSERT_GE(nhit*10, nray);
  ASSERT_LE(nhit*10, 9*nray);
}

TEST(TestVCLVCLAlignedRectangularAperture, CompareToScalar) {
  using Real = VCL256DoubleReal;
  calin::math::rng::RNG rng(__PRETTY_FUNCTION__,"Ray generator");
  unsigned nray = 100000;
  unsigned nhit = 0;
  for(unsigned iray=0; iray<nray; ++iray) {
    double x0 = rng.uniform()*2.0 - 1.0;
    double y0 = rng.uniform()*2.0 - 1.0;
    double z0 = rng.uniform()*2.0 - 1.0;
    double dx = rng.uniform();
    double dz = rng.uniform();

    Eigen::Vector3d center(x0,y0,z0);
    calin::simulation::vs_optics::VSOAlignedRectangularAperture sobs(center, dx, dz);
    VCLAlignedRectangularAperture<Real> vobs(center, dx, dz);

    double x = rng.uniform()*2.0 - 1.0;
    double y = rng.uniform()*2.0 - 1.0;
    double z = rng.uniform()*2.0 - 1.0;
    double uy = rng.uniform()*2.0 - 1.0;
    double theta = rng.uniform()*2.0*M_PI;
    double ur = std::sqrt(1.0-uy*uy);
    double ux = ur*std::cos(theta);
    double uz = ur*std::sin(theta);
    Eigen::Vector3d r(x,y,z);
    Eigen::Vector3d u(ux,uy,uz);

    calin::math::ray::Ray sray_in(r,u);
    calin::math::ray::Ray sray_out;
    bool sintersects = sobs.doesObscure(sray_in, sray_out, 1.0);

    calin::math::ray::VCLRay<Real> vray_in(
      r.template cast<Real::real_vt>(), u.template cast<Real::real_vt>());
    calin::math::ray::VCLRay<Real> vray_out;
    Real::bool_vt vintersects = vobs.doesObscure(vray_in, vray_out, 1.0);

    ASSERT_EQ(sintersects, vintersects[0])
      << sintersects << ' ' << vintersects << ' ' << iray << '\n'
      << r.transpose() << '\n' << u.transpose();

    if(sintersects)++nhit;
  }
  ASSERT_GE(nhit*10, nray);
  ASSERT_LE(nhit*10, 9*nray);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
