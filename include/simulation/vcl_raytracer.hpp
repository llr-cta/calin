/*

   calin/simulation/vcl_raytracer.hpp -- Stephen Fegan -- 2018-09-10

   Class for raytracing on a single VSOTelescope using VCL

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

#pragma once

#include <util/vcl.hpp>
#include <math/ray_vcl.hpp>
#include <math/rng_vcl.hpp>
#include <simulation/vso_telescope.hpp>

namespace calin { namespace simulation { namespace vcl_raytracer {

enum ScopeTraceStatus {
  STS_TRAVELLING_AWAY_REFLECTOR,
  STS_MISSED_REFLECTOR_SPHERE,
  STS_OUTSIDE_REFLECTOR_APERTURE,
  STS_NO_MIRROR,
  STS_MISSED_MIRROR_SPHERE,
  STS_MISSED_MIRROR_EDGE,
  STS_OBSCURED_BEFORE_MIRROR,
      TS_ABSORBED_AT_MIRROR,
      TS_MISSED_WINDOW,
      TS_OBSCURED_BEFORE_FOCAL_PLANE,
  STS_TRAVELLING_AWAY_FROM_FOCAL_PLANE,
      TS_NO_PIXEL,
      TS_ABSORBED_AT_CONCENTRATOR,
      TS_PE_GENERATED
};

template<typename VCLReal> class ScopeTraceInfo: public VCLReal
{
public:
  using typename VCLReal::real_t;
  using typename VCLReal::real_vt;
  using typename VCLReal::bool_vt;
  using typename VCLReal::int_vt;
  using typename VCLReal::uint_vt;
  using typename VCLReal::vec3_vt;
  using typename VCLReal::mat3_vt;
  using typename calin::math::ray::VCLRay<VCLReal> Ray;

  int_vt              status;
  real_vt             reflec_x;
  real_vt             reflec_z;
  int_vt              mirror_hexid;
  int_vt              mirror_id;
  vec3_vt             mirror_pos;
  mat3_vt             mirror_rot;
  real_vt             mirror_rad;

#if 0
  real_vt             reflec_dx;
  real_vt             reflec_dz;
  const VSOMirror*    mirror;
  double              mirror_normal_dispersion;
  Eigen::Vector3d     mirror_scattered;
  double              mirror_reflection_angle;
  double              fplane_x;
  double              fplane_z;
  double              fplane_dx;
  double              fplane_dz;
  double              fplane_t;
  double              fplane_uy;     // y-axis directional cosine at FP
  int                 pixel_hexid;
  const VSOPixel*     pixel;
  double              pixel_dist;
  bool                concentrator_hit;
  unsigned            obscuration_id;
  const VSOObscuration* obscuration;

  void reset();
  std::ostream& write(std::ostream& stream = std::cout,
                      bool convert_to_physical_units = true,
                      bool end_of_line = true) const;

  bool rayWasReflected() const
  { return (int)status >= (int)TS_TRAVELLING_AWAY_FROM_FOCAL_PLANE; }
  bool rayHitFocalPlane() const
  { return (int)status >= (int)TS_NO_PIXEL; }
#endif
};

template<typename VCLRealType> class ScopeRayTracer: public VCLRealType
{
 public:
   using typename VCLRealType::real_t;
   using typename VCLRealType::int_t;
   using typename VCLRealType::uint_t;
   using typename VCLRealType::real_vt;
   using typename VCLRealType::bool_vt;
   using typename VCLRealType::vec3_vt;
   using typename VCLRealType::mat3_vt;
   using typename calin::math::ray::VCLRay<VCLRealType> Ray;

  // RayTracer(const VSOArray* array, math::rng::RNG* rng):
  //   fArray(array), fRNG(rng) { /* nothing to see here */ }
  // RayTracer(math::rng::RNG* rng, const VSOArray* array = nullptr):
  //   fArray(array), fRNG(rng) { /* nothing to see here */ }
  // ~VSORayTracer();

  bool_vt trace_scope_centered_global_frame(Ray& ray, ScopeTraceInfo& info)
  {
    // *************************************************************************
    // *************** RAY STARTS IN SCOPE CENTERED GLOBAL FRAME ***************
    // *************************************************************************

    ray.rotate(global_to_reflector_rot_);
    bool_vt mask = trace_reflector_frame(ray, info);
    ray.derotate(global_to_reflector_rot_);
    return mask;
  }

  bool_vt trace_reflector_frame(Ray& ray, ScopeTraceInfo& info)
  {
    bool_vt mask = ray.uz() < 0;
    info.status = STS_TRAVELLING_AWAY_REFLECTOR,

    // *************************************************************************
    // ****************** RAY STARTS IN RELECTOR COORDINATES *******************
    // *************************************************************************

    // Test for obscuration - to come

    // Propagate to intersection with the reflector sphere
    info.status = select(mask, STS_MISSED_REFLECTOR_SPHERE, info.status);
    mask = ray.propagate_to_y_sphere_2nd_interaction_fwd_only_with_mask(mask,
      curvature_radius_, 0, ref_index_);

    info.reflec_x     = select(mask, ray.position().x(), 0);
    info.reflec_z     = select(mask, ray.position().z(), 0);

    // Test aperture
    info.status = select(mask, STS_OUTSIDE_REFLECTOR_APERTURE, info.status);
    real_vt reflec_r2 =
    mask &= (info.reflec_x*info.reflec_x + info.reflec_y*info.reflec_y) <= reflec_aperture2;

    // Assume mirrors on hexagonal grid - use hex_array routines to find which hit
    info.status = select(mask, STS_NO_MIRROR, info.status);
    info.mirror_hexid = calin::math::hex_array::VLCReal<VCLRealType>::
      xy_trans_to_hexid_scaleinv(info.reflec_x, info.reflec_z,
        reflec_crot_, reflec_srot_, reflec_scaleinv_, reflec_shift_x_, reflec_shift_z_,
        reflec_cw_);

    // Test we have a valid mirror hexid
    mask &= info.mirror_hexid < mirror_hexid_end_;

    // Find the mirror ID
    info.mirror_id =
      vcl::lookup<0x40000000>(select(mask_, info.mirror_hexid, mirror_hexid_end_)),
        mirror_id_lookup_);

    // Test we have a valid mirror id
    mask &= info.mirror_id < mirror_id_end_;

    info.mirror_pos
      << vcl::lookup<0x40000000>(info.mirror_id, mirror_x_lookup_)
      << vcl::lookup<0x40000000>(info.mirror_id, mirror_y_lookup_)
      << vcl::lookup<0x40000000>(info.mirror_id, mirror_z_lookup_);

    ray.translate_origin(info.mirror_pos)

    info.mirror_rot
      << vcl::lookup<0x40000000>(info.mirror_id, mirror_m00_lookup_)
      << vcl::lookup<0x40000000>(info.mirror_id, mirror_m01_lookup_)
      << vcl::lookup<0x40000000>(info.mirror_id, mirror_m02_lookup_)
      << vcl::lookup<0x40000000>(info.mirror_id, mirror_m10_lookup_)
      << vcl::lookup<0x40000000>(info.mirror_id, mirror_m11_lookup_)
      << vcl::lookup<0x40000000>(info.mirror_id, mirror_m12_lookup_)
      << vcl::lookup<0x40000000>(info.mirror_id, mirror_m20_lookup_)
      << vcl::lookup<0x40000000>(info.mirror_id, mirror_m21_lookup_)
      << vcl::lookup<0x40000000>(info.mirror_id, mirror_m22_lookup_);

    ray.derotate(info.mirror_rot);

    // *************************************************************************
    // ******************* RAY IS NOW IN MIRROR COORDINATES ********************
    // *************************************************************************

    info.mirror_rad = vcl::lookup<0x40000000>(info.mirror_id, mirror_r_lookup_);

    // Propagate to intersection with the reflector sphere
    info.status = select(mask, STS_MISSED_MIRROR_SPHERE, info.status);
    mask = ray.propagate_to_y_sphere_2nd_interaction_fwd_bwd_with_mask(mask,
      info.mirror_rad, 0, ref_index_);

    // Verify that ray impacts inside of hexagonal mirror surface
    const real_vt cos60 = 0.5;
    const real_vt sin60 = 0.5*CALIN_HEX_ARRAY_SQRT3;

    const real_vt x_cos60 = ray.x() * cos60;
    const real_vt z_sin60 = ray.z() * cos60;

    const real_vt dhex_pos60 = abs(x_cos60 - z_sin60);
    const real_vt dhex_neg60 = abs(x_cos60 + z_sin60);

    mask &= max(max(dhex_pos60, dhex_neg60), abs(ray.x())) < mirror_dhex_max_;

    // Calculate mirror normal at impact point
    info.mirror_normal =
      (vec3_vt(0,info.mirror_rad,0) - ray.position()).normalized();

    // Scatter the normal to account for the spot size ot the focal length of the
    // radius. The spot size is given as the DIAMETER at the focal distance.
    // Must divide by 2 (for reflection) and another 2 for diameter -> radius
    info.mirror_normal_dispersion =
      vcl::lookup<0x40000000>(info.mirror_id, mirror_normdisp_lookup_);

    info.mirror_scattered = info.mirror_normal;
    calin::math::geometry::VCL<VCLRealType>::scatter_direction(
      info.mirror_scattered, info.mirror_normal_dispersion, *rng_);

    // Reflect ray
    ray.reflect_from_surface_with_mask(mask, info.mirror_scattered);

    info.mirror->mirrorToReflector(ray);

    // Translate back to reflector frame
    ray.rotate(mirror_rot);
    ray.untranslate_origin(mirror_pos)

    // *************************************************************************
    // *************** RAY IS NOW BACK IN REFLECTOR COORDINATES ****************
    // *************************************************************************

    // Check obsucrations

    // Refract in window

    ray.translate(fp_pos_)
    if(fp_has_rot_)ray.rotate(fp_rot_);

    // *************************************************************************
    // ***************** RAY IS NOW IN FOCAL PLANE COORDINATES *****************
    // *************************************************************************

    // Propagate to focal plane
    info.status = select(mask, STS_TRAVELLING_AWAY_FROM_FOCAL_PLANE, info.status);
    mask = ray.propagate_to_y_plane_with_mask(mask, 0, false, ref_index);

    // Test for interaction with obscuration before focal plane was hit

    // We good, record position on focal plane etc
    info.fplane_x = select(mask, ray.x(), 0);
    info.fplane_z = select(mask, ray.z(), 0);
    info.fplane_t = select(mask, ray.ct(), 0) * math::constants::cgs_1_c;
    info.fplane_uy = select(mask, ray.uy(), 0);



    info.pixel_hexid =
      math::hex_array::VCLReal<VCLRealType>::xy_trans_to_hexid(
        info.fplane_x, info.fplane_z,
        pixel_crot_, pixel_srot_, pixel_scaleinv_, pixel_shift_x_, pixel_shift_z_,
        pixel_cw_);

    // Find pixel (if there is a real pixel at that site)
    info.pixel = info.scope->pixelByHexID(info.pixel_hexid);
    if(info.pixel==0)
    {
      info.status = TS_NO_PIXEL;
      info.scope->focalPlaneToGlobal(ray);
      return 0;
    }

    if(fp_has_rot_)ray.derotate(fp_rot_);
    ray.untranslate(fp_pos_)

    // *************************************************************************
    // ************ RAY IS NOW BACK IN REFLECTOR COORDINATES AGAIN *************
    // *************************************************************************



  }

 private:
   mat3_vt         global_to_reflector_rot_;
};

// std::ostream& operator <<(std::ostream& stream, const VSORayTracer::TraceInfo& o);

} } } // namespace calin::simulations::vs_optics
