/*

   calin/unit_tests/io/test_vcl_raytracer.cpp -- Stephen Fegan -- 2018-09-17

   Unit tests for VCL raytracer

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
#include <fstream>
#include <limits>
#include <gtest/gtest.h>

#include <math/rng.hpp>
#include <math/hex_array.hpp>
#include <simulation/vso_array.hpp>
#include <simulation/vs_optics.pb.h>

#include <util/vcl.hpp>
#include <math/ray_vcl.hpp>
#include <simulation/vcl_raytracer.hpp>

using namespace calin::util::vcl;
using namespace calin::simulation::vcl_raytracer;

calin::simulation::vs_optics::VSOArray* make_test_array()
{
  calin::ix::simulation::vs_optics::IsotropicDCArrayParameters mst;
  mst.mutable_array_origin()->set_latitude(28.0 + 45.0/60.0 + 47.36/3600.0);
  mst.mutable_array_origin()->set_longitude(-(17.0 + 53.0/60.0 + 23.93/3600.0));
  mst.mutable_array_origin()->set_elevation(2147 * 100.0);

  auto* scope = mst.mutable_prescribed_array_layout()->add_scope_positions();
  scope->set_x(0);
  scope->set_y(0);
  scope->set_z(mst.array_origin().elevation());

  mst.mutable_reflector_frame()->set_optic_axis_rotation(-90.0);

  auto* dc = mst.mutable_reflector();
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
  for(auto id : { 1,67,72,82,87 })dc->add_facet_missing_list(id-1);
  mst.mutable_focal_plane()->set_camera_diameter(235);
  mst.mutable_focal_plane()->mutable_translation()->set_y(1.0/(1.0/dc->alignment_image_plane()-1.0/(10 * 1e5)));
  mst.mutable_pixel()->set_spacing(5);
  mst.mutable_pixel()->set_cone_inner_diameter(5);
  mst.mutable_pixel()->set_cone_survival_prob(1);
  mst.mutable_pixel()->set_hex_module_size(1);
  mst.mutable_pixel()->set_hex_module_layout_use_b_configuration(true);
  mst.mutable_pixel()->set_module_num_hex_rings(9);

  int u1, v1;
  double x1,y1;
  calin::math::hex_array::cluster_hexid_to_center_uv(1,1,u1,v1,false);
  calin::math::hex_array::uv_to_xy(u1,v1,x1,y1);
  mst.mutable_pixel()->set_grid_rotation(atan2(-y1,x1)/M_PI*180.0 + 30.0);

#if 0
  if(obscure_camera):
      obs_camera_box = mst.add_pre_reflection_obscuration()
      obs_camera_box.mutable_aligned_box().mutable_max_corner().set_x(2918.0/20)
      obs_camera_box.mutable_aligned_box().mutable_max_corner().set_z(2918.0/20)
      obs_camera_box.mutable_aligned_box().mutable_max_corner().set_y(mst.focal_plane().translation().y()+(1532.25-513.0)/10)
      obs_camera_box.mutable_aligned_box().mutable_min_corner().set_x(-2918.0/20)
      obs_camera_box.mutable_aligned_box().mutable_min_corner().set_z(-2918.0/20)
      obs_camera_box.mutable_aligned_box().mutable_min_corner().set_y(mst.focal_plane().translation().y()-513.0/10)

      obs_outer_aperture = mst.add_post_reflection_obscuration()
      obs_outer_aperture.mutable_rectangular_aperture().mutable_center_pos().set_x(0)
      obs_outer_aperture.mutable_rectangular_aperture().mutable_center_pos().set_z(-50.0/20)
      obs_outer_aperture.mutable_rectangular_aperture().mutable_center_pos().set_y(mst.focal_plane().translation().y()-513.0/10)
      obs_outer_aperture.mutable_rectangular_aperture().set_flat_to_flat_x(2714.0/10)
      obs_outer_aperture.mutable_rectangular_aperture().set_flat_to_flat_z((2585.0-50.0)/10)

      obs_inner_aperture = mst.add_post_reflection_obscuration()
      obs_inner_aperture.mutable_circular_aperture().mutable_center_pos().set_x(0)
      obs_inner_aperture.mutable_circular_aperture().mutable_center_pos().set_z(0)
      obs_inner_aperture.mutable_circular_aperture().mutable_center_pos().set_y(mst.focal_plane().translation().y()-222.5/10)
      obs_inner_aperture.mutable_circular_aperture().set_diameter(2304.0/10)

  if include_window:
      win = mst.mutable_spherical_window()
      win.set_front_y_coord(mst.focal_plane().translation().y() - 427.5/10)
      win.set_outer_radius(3443.0/10)
      win.set_thickness(5.25/10.0)
      win.set_refractive_index(1.5)
  return mst
  #endif

  calin::simulation::vs_optics::VSOArray* array =
    new calin::simulation::vs_optics::VSOArray;
  calin::math::rng::RNG rng(__PRETTY_FUNCTION__, "VSOArray population");
  array->generateFromArrayParameters(mst, rng);

  return array;
}

template<typename VCLRealType> class TestVCLRaytracer :
  public VCLRealType, public testing::Test
{
public:
};

using RealTypes = ::testing::Types<VCL128FloatReal, VCL256FloatReal, VCL512FloatReal, VCL128DoubleReal, VCL256DoubleReal, VCL512DoubleReal>;

TYPED_TEST_CASE(TestVCLRaytracer, RealTypes);

TYPED_TEST(TestVCLRaytracer, RayTrace) {
  calin::simulation::vs_optics::VSOArray* array = make_test_array();
  ScopeRayTracer<TypeParam> raytracer(array->telescope(0));

  for(const auto* mirror : array->telescope(0)->all_mirrors()) {
    calin::math::ray::VCLRay<TypeParam> ray;
    ray.mutable_direction() << 0.0 , -1.0, 0.0;
    ray.mutable_position() << mirror->pos().x() , 2000.0, mirror->pos().z();
    typename TypeParam::bool_vt mask = true;
    mask.insert(1, false);
    typename ScopeRayTracer<TypeParam>::TraceInfo trace_info;
    mask = raytracer.trace_reflector_frame(mask, ray, trace_info);
    ASSERT_TRUE(mask[0]);
    ASSERT_FALSE(mask[1]);
    ASSERT_EQ(trace_info.status[0], calin::simulation::vcl_raytracer::STS_TS_FOUND_PIXEL);
    ASSERT_EQ(trace_info.status[1], calin::simulation::vcl_raytracer::STS_MASKED_ON_ENTRY);
    ASSERT_EQ(trace_info.mirror_hexid[0], mirror->hexID());
    ASSERT_EQ(trace_info.mirror_id[0], mirror->id());
    ASSERT_EQ(trace_info.mirror_hexid[1], array->telescope(0)->numMirrorHexSites());
    ASSERT_EQ(trace_info.mirror_id[1], array->telescope(0)->numMirrors());
  }

  delete array;
}

TEST(TestVCLRaytracer, PSF) {
  calin::simulation::vs_optics::VSOArray* array = make_test_array();
  ScopeRayTracer<VCL256FloatReal> raytracer(array->telescope(0));
  calin::math::rng::VCLRealRNG<VCL256FloatReal> rng(__PRETTY_FUNCTION__,"Ray position generator");

  for(unsigned iray=0; iray<100000; iray++) {
    calin::math::ray::VCLRay<VCL256FloatReal> ray;
    ray.mutable_direction() << 0.0 , -1.0, 0.0;
    ray.mutable_position() <<
      rng.uniform_zc(array->telescope(0)->reflectorIP()),
      2000.0,
      rng.uniform_zc(array->telescope(0)->reflectorIP());
    typename VCL256FloatReal::bool_vt mask = true;
    typename ScopeRayTracer<VCL256FloatReal>::TraceInfo trace_info;
    auto pos = ray.position();
    mask = raytracer.trace_reflector_frame(mask, ray, trace_info);
    for(unsigned i=0;i<8;i++)std::cout << pos.x()[i] << ' ' << pos.z()[i] << ' ' << mask[i] << ' ' << trace_info.status[i] << '\n';
  }

  delete array;
}


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
