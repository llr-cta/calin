/*

   calin/simulation/vcl_panoseti_raytracer.hpp -- Stephen Fegan -- 2026-07-23

   Class for raytracing on a single panoseti telescope using VCL

   Copyright 2026, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#pragma once

#include <algorithm>
#include <limits>

#include <util/memory.hpp>
#include <util/vcl.hpp>
#include <math/constants.hpp>
#include <math/special.hpp>
#include <math/ray_vcl.hpp>
#include <math/rng_vcl.hpp>
#include <math/geometry_vcl.hpp>
#include <util/log.hpp>
#include <simulation/panoseti.pb.h>

namespace calin { namespace simulation { namespace vcl_raytracer {

// These classes are too complex for SWIG to handle directly from the header
// files. We define a simplified subset directly in the SWIG interface file.

#ifndef SWIG

enum VCLPanosetiScopeTraceStatus {
  PSTS_MASKED_ON_ENTRY,                        // 0
  PSTS_TRAVELLING_AWAY_REFLECTOR,              // 1
  PSTS_MISSED_REFLECTOR_SPHERE,                // 2
  PSTS_OUTSIDE_REFLECTOR_APERTURE,             // 3
  PSTS_NO_MIRROR,                              // 4
  PSTS_MISSED_MIRROR_SPHERE,                   // 5
  PSTS_MISSED_MIRROR_EDGE,                     // 6
  PSTS_OBSCURED_BEFORE_MIRROR,                 // 7
  PSTS_OBSCURED_BEFORE_FOCAL_PLANE,            // 9
  PSTS_TRAVELLING_AWAY_FROM_FOCAL_PLANE,       // 10
  PSTS_OUTSIDE_FOCAL_PLANE_APERTURE,           // 11
  PSTS_TS_NO_PIXEL,                            // 12
  PSTS_TS_FOUND_PIXEL                          // 13
};

template<typename VCLReal> class alignas(VCLReal::vec_bytes) VCLPanosetiScopeTraceInfo: public VCLReal
{
public:
  using typename VCLReal::real_t;
  using typename VCLReal::real_vt;
  using typename VCLReal::real_bvt;
  using typename VCLReal::int_vt;
  using typename VCLReal::uint_vt;
  using typename VCLReal::vec3_vt;
  using typename VCLReal::mat3_vt;

  int_vt              status;          // Status of ray at end of tracing

  real_vt             reflec_x;        // Ray intersection point on reflector sphere
  real_vt             reflec_y;        // Ray intersection point on reflector sphere
  real_vt             reflec_z;        // Ray intersection point on reflector sphere

  int_vt              lens_groove      // Lens groove number for
  real_vt             lens_x_out;      // Ray exit point from lens
  real_vt             lens_y_out;      // Ray exit point from lens
  real_vt             lens_z_out;      // Ray exit point from lens
  real_vt             lens_nin_dot_u;  // Cosine if angle between ray and incoming lens normal
  real_vt             lens_nout_dot_u; // Cosine if angle between ray and utgoing lens normal

  real_vt             fplane_x;         // Ray intersection point on focal plane
  real_vt             fplane_z;         // Ray intersection point on focal plane
  real_vt             fplane_t;         // Ray intersection time on focal plane
  real_vt             fplane_ux;        // X directional cosine of ray at focal plane
  real_vt             fplane_uy;        // Cosine of angle between ray and focal plane (normal)
  real_vt             fplane_uz;        // Z directional cosine of ray at focal plane

  int_vt              pixel_hexid;      // Grid hex ID of pixel on focal plane
  int_vt              pixel_id;         // Sequential ID of pixel on focal plane (or -1)
};

template<typename VCLReal> class alignas(VCLReal::vec_bytes) VCLPanosetiThinLensScopeRayTracer: public VCLReal
{
public:
  using typename VCLReal::real_t;
  using typename VCLReal::int_t;
  using typename VCLReal::uint_t;
  using typename VCLReal::int_bvt;
  using typename VCLReal::uint_bvt;
  using typename VCLReal::mat3_t;
  using typename VCLReal::vec3_t;
  using typename VCLReal::real_vt;
  using typename VCLReal::real_bvt;
  using typename VCLReal::int_vt;
  using typename VCLReal::uint_vt;
  using typename VCLReal::vec3_vt;
  using typename VCLReal::mat3_vt;
  using Ray = calin::math::ray::VCLRay<VCLReal>;
  using TraceInfo = VCLScopeTraceInfo<VCLReal>;
  using RNG = calin::math::rng::VCLRealRNG<VCLReal>;

  VCLPanosetiThinLensScopeRayTracer(const& calin::ix::simulation::panoseti::ArrayParameters& array_params,
      unsigned scope_id, real_t air_refractive_index = 1.0,
      RNG* rng = nullptr, bool adopt_rng = false):
    VCLReal(), rng_(rng==nullptr ? new RNG(__PRETTY_FUNCTION__) : rng),
    adopt_rng_(rng==nullptr ? true : adopt_rng)
  {
    using calin::math::special::SQR;

    global_to_reflector_off_ = scope->translationGlobalToReflector().cast<real_t>();
    global_to_reflector_rot_ = scope->rotationGlobalToReflector().cast<real_t>();
    air_ref_index_           = refractive_index;

    reflec_curvature_radius_ = scope->curvatureRadius();
    reflec_aperture2_        = 0;
    reflec_crot_             = scope->cosReflectorRotation();
    reflec_srot_             = scope->sinReflectorRotation();

    fp_pos_                  = scope->focalPlanePosition().cast<real_t>();
    fp_has_rot_              = scope->hasFPRotation();
    fp_rot_                  = scope->rotationReflectorToFP().cast<real_t>();
    fp_aperture2_            = 0;

    pixel_crot_              = scope->cosPixelRotation();
    pixel_srot_              = scope->sinPixelRotation();
    pixel_scaleinv_          = 1.0/scope->pixelSpacing();
    pixel_shift_x_           = scope->pixelGridShiftX();
    pixel_shift_z_           = scope->pixelGridShiftZ();
    pixel_cw_                = scope->pixelParity();

    pixel_hexid_end_         = scope->numPixelHexSites();
    pixel_id_end_            = scope->numPixels();

    fp_aperture2_ = SQR(std::sqrt(fp_aperture2_) + scope->pixelSpacing());

#if 0
    std::cout << pixel_crot_ << ' ' << pixel_srot_ << ' ' << pixel_scaleinv_ << ' '
      << pixel_shift_x_ << ' ' << pixel_shift_z_ << ' ' << pixel_cw_ << ' '
      << pixel_hexid_end_ << ' ' << pixel_id_end_ << ' '
      << fp_aperture2_ << ' ' << std::sqrt(fp_aperture2_) << '\n';
#endif
  }

  ~VCLPanosetiThinLensScopeRayTracer()
  {
    if(adopt_rng_)delete rng_;
  }

  void point_telescope(const calin::simulation::vs_optics::VSOTelescope* scope) {
    global_to_reflector_off_ = scope->translationGlobalToReflector().cast<real_t>();
    global_to_reflector_rot_ = scope->rotationGlobalToReflector().cast<real_t>();
  }

  static void transform_to_scope_reflector_frame(Ray& ray,
      const calin::simulation::vs_optics::VSOTelescope* scope)
  {
    ray.translate_origin(scope->translationGlobalToReflector().cast<real_vt>());
    ray.rotate(scope->rotationGlobalToReflector().cast<real_vt>());
  }

  real_bvt trace_global_frame(real_bvt mask, Ray& ray, TraceInfo& info,
    bool do_derotation = true)
  {
    // *************************************************************************
    // ********************** RAY STARTS IN GLOBAL FRAME ***********************
    // *************************************************************************

    ray.translate_origin(global_to_reflector_off_.template cast<real_vt>());
    ray.rotate(global_to_reflector_rot_.template cast<real_vt>());
    mask = trace_reflector_frame(mask, ray, info);
    if(do_derotation) {
      ray.derotate(global_to_reflector_rot_.template cast<real_vt>());
      ray.untranslate_origin(global_to_reflector_off_.template cast<real_vt>());
    }
    return mask;
  }

  real_bvt trace_scope_centered_global_frame(real_bvt mask, Ray& ray, TraceInfo& info,
    bool do_derotation = true)
  {
    // *************************************************************************
    // *************** RAY STARTS IN SCOPE CENTERED GLOBAL FRAME ***************
    // *************************************************************************

    ray.rotate(global_to_reflector_rot_.template cast<real_vt>());
    mask = trace_reflector_frame(mask, ray, info);
    if(do_derotation) {
      ray.derotate(global_to_reflector_rot_.template cast<real_vt>());
    }
    return mask;
  }

//#define DEBUG_STATUS

  real_bvt trace_reflector_frame(real_bvt mask, Ray& ray, TraceInfo& info)
  {
    info.status = STS_MASKED_ON_ENTRY;
#ifdef DEBUG_STATUS
    std::cout << mask[0] << '/' << info.status[0];
#endif

    info.status = select(int_bvt(mask), STS_TRAVELLING_AWAY_REFLECTOR, info.status);
    mask &= ray.uy() < 0;
#ifdef DEBUG_STATUS
    std::cout << ' ' << mask[0] << '/' << info.status[0];
#endif

    // *************************************************************************
    // ****************** RAY STARTS IN RELECTOR COORDINATES *******************
    // *************************************************************************

    // Test for obscuration of incoming ray
    real_bvt was_obscured = false;
    real_vt ct_obscured = std::numeric_limits<real_t>::infinity();
    uint_vt hitmask = 1;
    info.pre_reflection_obs_hitmask = uint_vt(0);
    for(const auto* obs : pre_reflection_obscuration) {
      Ray ray_out;
      real_bvt was_obscured_here = obs->doesObscure(ray, ray_out, ref_index_);
      ct_obscured = vcl::select(was_obscured_here,
        vcl::min(ct_obscured, ray_out.ct()), ct_obscured);
      was_obscured |= was_obscured_here;
      info.pre_reflection_obs_hitmask |= vcl::select(uint_bvt(was_obscured_here), hitmask, 0);
      hitmask <<= 1;
    }

    // Remember initial ct to test reflection happens after emission
    real_vt ct0 = ray.ct();

    // Propagate to intersection with the reflector sphere (allow to go backwards a bit)
    info.status = select(int_bvt(mask), STS_MISSED_REFLECTOR_SPHERE, info.status);
    mask = ray.propagate_to_y_sphere_2nd_interaction_mostly_fwd_with_mask(mask,
      reflec_curvature_radius_, 0, (-2.0/CALIN_HEX_ARRAY_SQRT3)*mirror_dhex_max_, ref_index_);
#ifdef DEBUG_STATUS
    std::cout << ' ' << mask[0] << '/' << info.status[0];
#endif

    info.reflec_x     = select(mask, ray.position().x(), 0);
    info.reflec_y     = select(mask, ray.position().y(), 0);
    info.reflec_z     = select(mask, ray.position().z(), 0);

    // Test aperture
    info.status = select(int_bvt(mask), STS_OUTSIDE_REFLECTOR_APERTURE, info.status);
    mask &= (info.reflec_x*info.reflec_x + info.reflec_z*info.reflec_z) <= reflec_aperture2_;
#ifdef DEBUG_STATUS
    std::cout << ' ' << mask[0] << '/' << info.status[0];
#endif

    if(not horizontal_or(mask)) {
      // We outie ...
      return mask;
    }

    // Assume mirrors on hexagonal grid - use hex_array routines to find which hit
    info.status = select(int_bvt(mask), STS_NO_MIRROR, info.status);
    info.mirror_hexid = calin::math::hex_array::VCLReal<VCLReal>::
      xy_trans_to_hexid_scaleinv(info.reflec_x, info.reflec_z,
        reflec_crot_, reflec_srot_, reflec_scaleinv_, reflec_shift_x_, reflec_shift_z_,
        reflec_cw_);

    // Test we have a valid mirror hexid
    mask &= typename VCLReal::real_bvt(info.mirror_hexid < mirror_hexid_end_);
#ifdef DEBUG_STATUS
    std::cout << ' ' << mask[0] << '/' << info.status[0];
#endif

    if(not horizontal_or(mask)) {
      // We outie ...
      return mask;
    }

    info.mirror_hexid = select(int_bvt(mask), info.mirror_hexid, mirror_hexid_end_);

    // Find the mirror ID
    info.mirror_id = vcl::lookup<0x40000000>(info.mirror_hexid, mirror_id_lookup_);

    // Test we have a valid mirror id
    mask &= typename VCLReal::real_bvt(info.mirror_id < mirror_id_end_);
#ifdef DEBUG_STATUS
    std::cout << ' ' << mask[0] << '/' << info.status[0];
#endif

    vec3_vt mirror_dir;
    mirror_dir.x() = vcl::lookup<0x40000000>(info.mirror_id, mirror_nx_lookup_);
    mirror_dir.z() = vcl::lookup<0x40000000>(info.mirror_id, mirror_nz_lookup_);
#if 1
    // Is it faster to use lookup table than to compute ?
    mirror_dir.y() = vcl::lookup<0x40000000>(info.mirror_id, mirror_ny_lookup_);
#else
    mirror_dir.y() = sqrt(nmul_add(mirror_dir.z(), mirror_dir.z(),
      nmul_add( mirror_dir.x(), mirror_dir.x(),1.0)));
#endif

    real_vt mirror_r = vcl::lookup<0x40000000>(info.mirror_id, mirror_r_lookup_);

    vec3_vt mirror_pos;
    mirror_pos.x() = vcl::lookup<0x40000000>(info.mirror_id, mirror_x_lookup_);
    mirror_pos.y() = vcl::lookup<0x40000000>(info.mirror_id, mirror_y_lookup_);
    mirror_pos.z() = vcl::lookup<0x40000000>(info.mirror_id, mirror_z_lookup_);

    vec3_vt mirror_center = mirror_pos + mirror_dir * mirror_r;

    ray.translate_origin(mirror_center);

    // *************************************************************************
    // ******************* RAY IS NOW IN MIRROR COORDINATES ********************
    // *************************************************************************

    // Propagate to intersection with the mirror sphere
    info.status = select(int_bvt(mask), STS_MISSED_MIRROR_SPHERE, info.status);
    mask = ray.propagate_to_y_sphere_2nd_interaction_fwd_bwd_with_mask(mask,
      mirror_r, -mirror_r, ref_index_);
    mask &= ray.ct() >= ct0;
#ifdef DEBUG_STATUS
    std::cout << ' ' << mask[0] << '/' << info.status[0];
#endif

    // Impact point relative to facet attchment point
    vec3_vt ray_pos = ray.position() + mirror_center - mirror_pos;
    calin::math::geometry::VCL<VCLReal>::
      derotate_in_place_Ry(ray_pos, reflec_crot_, reflec_srot_);

    calin::math::geometry::VCL<VCLReal>::
      derotate_in_place_Ry(mirror_dir, reflec_crot_, reflec_srot_);

    calin::math::geometry::VCL<VCLReal>::
      derotate_in_place_y_to_u_Ryxy(ray_pos, mirror_dir);

    // Verify that ray impacts inside of hexagonal mirror surface
    const real_vt cos60 = 0.5;
    const real_vt sin60 = 0.5*CALIN_HEX_ARRAY_SQRT3;

    const real_vt x_cos60 = ray_pos.x() * cos60;
    const real_vt z_sin60 = ray_pos.z() * sin60;

    const real_vt dhex_pos60 = abs(x_cos60 - z_sin60);
    const real_vt dhex_neg60 = abs(x_cos60 + z_sin60);

    info.status = select(int_bvt(mask), STS_MISSED_MIRROR_EDGE, info.status);
    mask &= max(max(dhex_pos60, dhex_neg60), abs(ray_pos.x())) < mirror_dhex_max_;
#ifdef DEBUG_STATUS
    std::cout << ' ' << mask[0] << '/' << info.status[0];
#endif

    // Calculate mirror normal at impact point
    vec3_vt mirror_normal = -ray.position() * (1.0/mirror_r);

    // Scatter the normal to account for the spot size ot the focal length of the
    // radius. The spot size is given as the DIAMETER at the focal distance.
    // Must divide by 2 (for reflection) and another 2 for diameter -> radius
    real_vt mirror_normal_dispersion =
      vcl::lookup<0x40000000>(info.mirror_id, mirror_normdisp_lookup_);

    // Scatter the normal direction randomly
    calin::math::geometry::VCL<VCLReal>::scatter_direction_in_place(
      mirror_normal, mirror_normal_dispersion, *rng_);

    // Reflect ray
#if 1
    info.mirror_n_dot_u = ray.direction().dot(mirror_normal);
    ray.mutable_direction() -= mirror_normal * select(mask, 2.0*info.mirror_n_dot_u, 0);
#else
    // Do not use this function any longer as we wish to keep u dot n
    ray.reflect_from_surface_with_mask(mask, info.mirror_normal_scattered);
#endif

    // Translate back to reflector frame
    ray.untranslate_origin(mirror_center);

    info.mirror_x     = select(mask, ray.position().x(), 0);
    info.mirror_y     = select(mask, ray.position().y(), 0);
    info.mirror_z     = select(mask, ray.position().z(), 0);

    // *************************************************************************
    // *************** RAY IS NOW BACK IN REFLECTOR COORDINATES ****************
    // *************************************************************************

    // Finish checking obscuration before mirror hit
    info.status = select(int_bvt(mask), STS_OBSCURED_BEFORE_MIRROR, info.status);
    mask &= ~(was_obscured & (ct_obscured < ray.ct()));

    if(not horizontal_or(mask)) {
      // We outie ...
      return mask;
    }

    // Test for obscuration on way to focal plane - first with obscurations
    // that are given in reflector coordinates (telescope arms etc)
    was_obscured = false;
    ct_obscured = std::numeric_limits<real_t>::infinity();
    hitmask = uint_vt(1);
    info.post_reflection_obs_hitmask = uint_vt(0);
    for(const auto* obs : post_reflection_obscuration) {
      Ray ray_out;
      real_bvt was_obscured_here = obs->doesObscure(ray, ray_out, ref_index_);
      ct_obscured = vcl::select(was_obscured_here,
        vcl::min(ct_obscured, ray_out.ct()), ct_obscured);
      was_obscured |= was_obscured_here;
      info.post_reflection_obs_hitmask |= vcl::select(uint_bvt(was_obscured_here), hitmask, 0);
      hitmask <<= 1;
    }

    // Refract in window

    ray.translate_origin(fp_pos_.template cast<real_vt>());
    if(fp_has_rot_)ray.rotate(fp_rot_.template cast<real_vt>());

    // *************************************************************************
    // ***************** RAY IS NOW IN FOCAL PLANE COORDINATES *****************
    // *************************************************************************

    // Test for obscuration on way to focal plane - second with obscurations
    // that are given in focal plane coordinates
    hitmask = uint_vt(1);
    info.camera_obs_hitmask = uint_vt(0);
    for(const auto* obs : camera_obscuration) {
      Ray ray_out;
      real_bvt was_obscured_here = obs->doesObscure(ray, ray_out, ref_index_);
      ct_obscured = vcl::select(was_obscured_here,
        vcl::min(ct_obscured, ray_out.ct()), ct_obscured);
      was_obscured |= was_obscured_here;
      info.camera_obs_hitmask |= vcl::select(uint_bvt(was_obscured_here), hitmask, 0);
      hitmask <<= 1;
    }

    // Propagate to focal plane
    info.status = select(int_bvt(mask), STS_TRAVELLING_AWAY_FROM_FOCAL_PLANE, info.status);
    mask = ray.propagate_to_y_plane_with_mask(mask, 0, false, ref_index_);
#ifdef DEBUG_STATUS
    std::cout << ' ' << mask[0] << '/' << info.status[0];
#endif

    // Finish checking obscuration after mirror reflection
    info.status = select(int_bvt(mask), STS_OBSCURED_BEFORE_FOCAL_PLANE, info.status);
    mask &= ~(was_obscured & (ct_obscured < ray.ct()));

    // We good, record position on focal plane etc
    info.fplane_x = select(mask, ray.x(), 0);
    info.fplane_z = select(mask, ray.z(), 0);
    info.fplane_t = select(mask, ray.time(), 0);
    info.fplane_ux = select(mask, ray.ux(), 0);
    info.fplane_uy = select(mask, ray.uy(), 0);
    info.fplane_uz = select(mask, ray.uz(), 0);

    info.status = select(int_bvt(mask), STS_OUTSIDE_FOCAL_PLANE_APERTURE, info.status);
    mask &= (info.fplane_x*info.fplane_x + info.fplane_z*info.fplane_z) <= fp_aperture2_;
#ifdef DEBUG_STATUS
    std::cout << ' ' << mask[0] << '/' << info.status[0];
#endif

    info.pixel_hexid =
    calin::math::hex_array::VCLReal<VCLReal>::
      xy_trans_to_hexid_scaleinv(info.fplane_x, info.fplane_z,
        pixel_crot_, pixel_srot_, pixel_scaleinv_, pixel_shift_x_, pixel_shift_z_,
        pixel_cw_);

    // Test we have a valid pixel hexid
    info.status = select(int_bvt(mask), STS_TS_NO_PIXEL, info.status);
    mask &= typename VCLReal::real_bvt(info.pixel_hexid < pixel_hexid_end_);
#ifdef DEBUG_STATUS
    std::cout << ' ' << mask[0] << '/' << info.status[0];
#endif

    // Find the pixel ID
    info.pixel_id =
      vcl::lookup<0x40000000>(select(int_bvt(mask), info.pixel_hexid, pixel_hexid_end_),
        pixel_id_lookup_);

    mask &= typename VCLReal::real_bvt(info.pixel_id < pixel_id_end_);
#ifdef DEBUG_STATUS
    std::cout << ' ' << mask[0] << '/' << info.status[0];
#endif

    info.status = select(int_bvt(mask), STS_TS_FOUND_PIXEL, info.status);
#ifdef DEBUG_STATUS
    std::cout << ' ' << mask[0] << '/' << info.status[0];
#endif

    if(fp_has_rot_)ray.derotate(fp_rot_.template cast<real_vt>());
    ray.untranslate_origin(fp_pos_.template cast<real_vt>());

    // *************************************************************************
    // ************ RAY IS NOW BACK IN REFLECTOR COORDINATES AGAIN *************
    // *************************************************************************

#ifdef DEBUG_STATUS
    std::cout << '\n';
#endif
    return mask;
  }

  unsigned psf(Eigen::VectorXd& x_out, Eigen::VectorXd& y_out, Eigen::VectorXd& t_out, unsigned nray,
    double theta = 0, double phi = 0, double distance = std::numeric_limits<double>::infinity(), double radius = 0) 
  {
    theta *= M_PI/180.0;
    double sintheta = std::sin(theta);
    phi *= M_PI/180.0;
    mat3_vt rot;
    calin::math::geometry::VCL<VCLReal>::rotation_y_to_xyz_Ryxy(rot, 
      sintheta*std::cos(phi), std::cos(theta), sintheta*std::sin(phi));

    if(radius <= 0) {
      radius = std::sqrt(reflec_aperture2_);
    }
    if(distance <= 0) {
      throw std::runtime_error("Light emission distance must be positive or infinity");
    }

    double dcospolar = 1.0 - std::cos(radius / distance);
    double fppos_2y = 2.0 * fp_pos_.y();

    unsigned ntraced = 0;
    unsigned iray = 0;
    x_out.resize(nray);
    y_out.resize(nray);
    t_out.resize(nray);

    while(iray < nray) {
      vec3_vt x;
      vec3_vt u;

      if(distance < std::numeric_limits<double>::infinity()) {
        x.x() = 0; 
        x.y() = distance;
        x.z() = 0;
        real_vt cospolar = 1.0 - dcospolar * rng_->uniform();
        real_vt sinpolar = vcl::sqrt(vcl::nmul_add(cospolar, cospolar, 1.0));
        real_vt cosazimuth;
        real_vt sinazimuth;
        rng_->sincos(sinazimuth, cosazimuth);
        u.x() = sinpolar * cosazimuth;
        u.y() = -cospolar;
        u.z() = sinpolar * sinazimuth;
      } else {
        u.x() = 0;
        u.y() = -1.0;
        u.z() = 0;
        real_vt rho = radius * vcl::sqrt(rng_->uniform());
        real_vt cosphi;
        real_vt sinphi;
        rng_->sincos(sinphi, cosphi);
        x.x() = rho * cosphi;
        x.y() = fppos_2y;   // Launch from twice the FP distance
        x.z() = rho * sinphi;
      }

      Ray ray(x, u, -x.y() * calin::math::constants::g4_1_c, 0);
      ray.rotate(rot);

      TraceInfo info;
      trace_reflector_frame(true, ray, info);

      typename VCLReal::int_at status;
      typename VCLReal::real_at xfp;
      typename VCLReal::real_at yfp;
      typename VCLReal::real_at tfp;
      info.status.store_a(status);
      info.fplane_x.store_a(xfp);
      info.fplane_z.store_a(yfp);
      info.fplane_t.store_a(tfp);

      for(unsigned i=0; i< VCLReal::num_real; i++) {
        ntraced++;
        if(status[i] >= STS_OUTSIDE_FOCAL_PLANE_APERTURE) {
          x_out[iray] = xfp[i];
          y_out[iray] = yfp[i];
          t_out[iray] = tfp[i];
          iray++;
          if(iray >= nray)break;
        }
      }
    }
    return ntraced;
  }

private:

  void populate_obscuration(std::vector<VCLObscuration<VCLReal>*>& to,
    const std::vector<const calin::simulation::vs_optics::VSOObscuration*>& from,
    const std::string& type)
  {
    using namespace calin::simulation::vs_optics;
    for(const auto* obs : from) {
      if(const auto* dc_obs = dynamic_cast<const VSOAlignedBoxObscuration*>(obs)) {
        to.push_back(new VCLAlignedBoxObscuration<VCLReal>(*dc_obs));
      } else if(const auto* dc_obs = dynamic_cast<const VSOAlignedRectangularAperture*>(obs)) {
        to.push_back(new VCLAlignedRectangularAperture<VCLReal>(*dc_obs));
      } else if(const auto* dc_obs = dynamic_cast<const VSOAlignedCircularAperture*>(obs)) {
        to.push_back(new VCLAlignedCircularAperture<VCLReal>(*dc_obs));
      } else if(const auto* dc_obs = dynamic_cast<const VSOTubeObscuration*>(obs)) {
        to.push_back(new VCLTubeObscuration<VCLReal>(*dc_obs));
      } else {
        throw std::runtime_error("Unsupported " + type + " obscuration type");
      }
    }
  }

  vec3_t          global_to_reflector_off_;
  mat3_t          global_to_reflector_rot_;
  real_t          ref_index_;

  real_t          reflec_curvature_radius_;
  real_t          reflec_aperture2_;
  real_t          reflec_crot_;
  real_t          reflec_srot_;
  real_t          reflec_scaleinv_;
  real_t          reflec_shift_x_;
  real_t          reflec_shift_z_;
  bool            reflec_cw_;

  int_t           mirror_hexid_end_;
  int_t           mirror_id_end_;
  int_t*          mirror_id_lookup_ = nullptr;
  real_t*         mirror_nx_lookup_ = nullptr;
  real_t*         mirror_nz_lookup_ = nullptr;
  real_t*         mirror_ny_lookup_ = nullptr;
  real_t*         mirror_r_lookup_ = nullptr;
  real_t*         mirror_x_lookup_ = nullptr;
  real_t*         mirror_z_lookup_ = nullptr;
  real_t*         mirror_y_lookup_ = nullptr;
  real_t          mirror_dhex_max_;
  real_t*         mirror_normdisp_lookup_ = nullptr;

  vec3_t          fp_pos_;
  bool            fp_has_rot_;
  mat3_t          fp_rot_;
  real_t          fp_aperture2_;

  real_t          pixel_crot_;
  real_t          pixel_srot_;
  real_t          pixel_scaleinv_;
  real_t          pixel_shift_x_;
  real_t          pixel_shift_z_;
  bool            pixel_cw_;

  int_t           pixel_hexid_end_;
  int_t           pixel_id_end_;
  int_t*          pixel_id_lookup_ = nullptr;

  std::vector<VCLObscuration<VCLReal>*> pre_reflection_obscuration;
  std::vector<VCLObscuration<VCLReal>*> post_reflection_obscuration;
  std::vector<VCLObscuration<VCLReal>*> camera_obscuration;

  RNG* rng_ = nullptr;
  bool adopt_rng_ = false;
};

} } } // namespace calin::simulations::vcl_raytracer
