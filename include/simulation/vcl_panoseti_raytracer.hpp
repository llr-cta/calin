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
#include <math/spline_interpolation.hpp>
#include <util/log.hpp>
#include <simulation/panoseti_optics.pb.h>

namespace calin { namespace simulation { namespace vcl_raytracer {

// These classes are too complex for SWIG to handle directly from the header
// files. We define a simplified subset directly in the SWIG interface file.

#ifndef SWIG

enum VCLPanosetiScopeTraceStatus {
  PSTS_MASKED_ON_ENTRY,                        // 0
  PSTS_TRAVELLING_AWAY_LENS,                   // 1
  PSTS_OUTSIDE_LENS_APERTURE,                  // 2
  PSTS_NO_LENS_EXIT,                           // 3
  PSTS_TRAVELLING_AWAY_FROM_DETECTOR_PLANE,    // 4
  PSTS_TS_OUTSIDE_DETECTOR,                    // 5
  PSTS_TS_FOUND_PIXEL                          // 6
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

  int_vt              status;           // Status of ray at end of tracing

  real_vt             lens_in_x;        // Ray intersection with lens aperture
  real_vt             lens_in_y;        // Ray intersection with lens aperture
  real_vt             lens_in_z;        // Ray intersection with lens aperture
  real_vt             lens_in_n_dot_u;  // Cosine of angle between ray and incoming lens normal

  int_vt              lens_groove;      // Lens groove number for
  real_vt             lens_out_x;       // Ray exit point from lens
  real_vt             lens_out_y;       // Ray exit point from lens
  real_vt             lens_out_z;       // Ray exit point from lens
  real_vt             lens_out_n_dot_u; // Cosine of angle between ray and outgoing lens normal

  real_vt             detector_x;        // Ray intersection point on detector plane
  real_vt             detector_z;        // Ray intersection point on detector plane
  real_vt             detector_t;        // Ray intersection time on detector plane
  real_vt             detector_ux;       // X directional cosine of ray at detector plane
  real_vt             detector_uy;       // Y directional cosine of ray at detector plane
  real_vt             detector_uz;       // Z directional cosine of ray at detector plane

  int_vt              pixel_id;          // Sequential ID of pixel on focal plane (or -1)
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
  using typename VCLReal::vecX_t;
  using Ray = calin::math::ray::VCLRay<VCLReal>;
  using TraceInfo = VCLPanosetiScopeTraceInfo<VCLReal>;
  using RNG = calin::math::rng::VCLRealRNG<VCLReal>;

  VCLPanosetiThinLensScopeRayTracer(const calin::simulation::panoseti_optics::ArrayParameters& array_params,
      unsigned scope_id, const calin::math::spline_interpolation::CubicSpline& lens_refractive_index_spline,
      real_t air_refractive_index = 1.0,
      RNG* rng = nullptr, bool adopt_rng = false):
    VCLReal(), lens_refractive_index_spline_(lens_refractive_index_spline),
    rng_(rng==nullptr ? new RNG(__PRETTY_FUNCTION__) : rng),
    adopt_rng_(rng==nullptr ? true : adopt_rng)
  {
    using calin::math::special::SQR;

    if(scope_id>=array_params.scope_positions()) {
      throw std::runtime_error("VCLPanosetiThinLensScopeRayTracer: scope_id out of range");
    }

    scope_position_.x()      = array_params.scope_positions(scope_id).x();
    scope_position_.y()      = array_params.scope_positions(scope_id).y();
    scope_position_.z()      = array_params.scope_positions(scope_id).z();
    pointTelescopeAzElPhi(0.0, 0.0, 0.0);

    air_ref_index_           = air_refractive_index;

    Eigen::VectorXd lens_polynomial = calin::protobuf_to_eigenvec(array_params.fresnel_lens_polynomial());
    lens_derivative_polynomial_ = calin::math::least_squares::polyder(lens_polynomial).cast<real_t>();

    lens_aperture2_          = SQR(array_params.fresnel_lens_aperture());

    detector_distance_       = array_params.detector_separation();
    detector_origin_         = calin::protobuf_to_eigenvec(array_params.detector_shift()).cast<real_t>();
    detector_origin_.y()     -= detector_distance_;
    if(calin::math::geometry::euler_is_zero(array_params.detector_rotation())) {
      detector_has_rotation_ = false;
      detector_rotation_ = mat3_t::Identity(); // Unused in this case, but set to identity for safety
    } else {
      detector_has_rotation_ = true;
      detector_rotation_ = calin::math::geometry::euler_to_matrix(array_params.detector_rotation()).cast<real_t>();
    }

    real_t roughness = array_params.fresnel_lens_roughness();
    scattering_sigma_theta_ = (detector_distance_ > 0) ? (roughness / detector_distance_) : 0.0;

    pixel_spacing_           = array_params.pixel_pitch();
    pixel_spacing_inv_       = 1.0/pixel_spacing_;
    pixel_nside_             = array_params.num_pixels_per_axis();
    pixel_array_halfwidth_   = pixel_spacing_ * pixel_nside_ / 2.0;
  }

  ~VCLPanosetiThinLensScopeRayTracer()
  {
    if(adopt_rng_)delete rng_;
  }

  bool pointTelescope(const Eigen::Vector3d& v)
  {
    if(v.squaredNorm()==0)return false;
    return pointTelescopeAzElPhi(atan2(v.x(),v.y()), atan2(v.z(),sqrt(v.x()*v.x() + v.y()*v.y())), 0.0);
  }

  bool pointTelescopeAzEl(const double az_rad, const double el_rad)
  {
    return pointTelescopeAzElPhi(az_rad, el_rad, 0.0);
  }

  bool pointTelescopeAzElPhi(double az_rad, double el_rad, double phi_rad)
  {
    az_rad_ = az_rad;
    el_rad_ = el_rad;
    phi_rad_ = phi_rad;
    Eigen::Matrix3d rot_reflector_to_global =
      Eigen::AngleAxisd(-az_rad_,   Eigen::Vector3d::UnitZ()) *
      Eigen::AngleAxisd(el_rad_,  Eigen::Vector3d::UnitX()) *
      Eigen::AngleAxisd(phi_rad_,   Eigen::Vector3d::UnitZ());
    Eigen::Matrix3d rot_global_to_reflector = rot_reflector_to_global.transpose();
    Eigen::Vector3d off_global_to_reflector_ = scope_position_;

    global_to_reflector_off_ = off_global_to_reflector_.template cast<real_t>(); // Cast double to float if necessary
    global_to_reflector_rot_ = rot_global_to_reflector.template cast<real_t>();
    return true;
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
    info.status = PSTS_MASKED_ON_ENTRY;
#ifdef DEBUG_STATUS
    std::cout << mask[0] << '/' << info.status[0];
#endif

    info.status = select(int_bvt(mask), PSTS_TRAVELLING_AWAY_LENS, info.status);
    mask &= ray.uy() < 0;
#ifdef DEBUG_STATUS
    std::cout << ' ' << mask[0] << '/' << info.status[0];
#endif

    // *************************************************************************
    // ****************** RAY STARTS IN RELECTOR COORDINATES *******************
    // *************************************************************************

    // Remember initial ct to test reflection happens after emission
    real_vt ct0 = ray.ct();

    // Propagate to lens plane
    info.status = select(int_bvt(mask), PSTS_OUTSIDE_LENS_APERTURE, info.status);
    mask = ray.propagate_to_y_plane_with_mask(mask, 0.0, /* time_reversal_ok = */ false, air_ref_index_);
#ifdef DEBUG_STATUS
    std::cout << ' ' << mask[0] << '/' << info.status[0];
#endif

    info.lens_in_x    = select(mask, ray.position().x(), 0);
    info.lens_in_y    = select(mask, ray.position().y(), 0);
    info.lens_in_z    = select(mask, ray.position().z(), 0);

    // Test aperture
    mask &= (info.lens_in_x*info.lens_in_x + info.lens_in_z*info.lens_in_z) <= lens_aperture2_;
#ifdef DEBUG_STATUS
    std::cout << ' ' << mask[0] << '/' << info.status[0];
#endif

    if(not horizontal_or(mask)) {
      // We outie ...
      return mask;
    }

    // Calculate reftactive index of lens
    real_vt lens_ref_index = lens_refractive_index_spline_.value(ray.energy());

    // Refract into lens
    ray.refract_at_surface_in_with_mask(mask, vec3d::UnitY(), lens_ref_index);

    // In thin lens approximation, we assume the ray exits the lens at the same point it 
    // entered, but with a new direction. Calculate normal from polynomial.
    vec3_vt lens_norm = calin::math::geometry::vcl<VCLReal>::norm_of_common_derivative_polynomial_surface(
      ray.x(), ray.z(), lens_derivative_polynomial_.data(), lens_derivative_polynomial_.size(),
      /* convex= */ true);

    // Refract out of lens
    info.status = select(int_bvt(mask), PSTS_NO_LENS_EXIT, info.status);
    mask = ray.refract_at_surface_out_with_mask(mask, lens_norm, air_ref_index_);
#ifdef DEBUG_STATUS
    std::cout << ' ' << mask[0] << '/' << info.status[0];
#endif

    if(scattering_sigma_theta_ > 0) {
      ray.scatter_direction(select(mask, scattering_sigma_theta_, 0.0), *rng_);
    }

    ray.translate_origin(detector_origin_.template cast<real_vt>());
    if(detector_has_rotation_) {
      ray.rotate(detector_rotation_.template cast<real_vt>());
    }

    // *************************************************************************
    // **************** RAY IS NOW IN DETECTOR PLANE COORDINATES ***************
    // *************************************************************************

    if(detector_has_rotation_) {
      // Test ray is travelling to detector plane.. should always be case unless
      // plane is strongly tilted
      info.status = select(int_bvt(mask), PSTS_TRAVELLING_AWAY_FROM_DETECTOR_PLANE, info.status);
      mask &= ray.uy() < 0;
#ifdef DEBUG_STATUS
      std::cout << ' ' << mask[0] << '/' << info.status[0];
#endif
    }

    // Propagate to detector plane
    ray.propagate_to_y_plane_with_mask(mask, 0.0, /* time_reversal_ok = */ false, air_ref_index_);

    // We good, record position on detector plane etc
    info.detector_x = select(mask, ray.x(), 0);
    info.detector_z = select(mask, ray.z(), 0);
    info.detector_t = select(mask, ray.time(), 0);
    info.detector_ux = select(mask, ray.ux(), 0);
    info.detector_uy = select(mask, ray.uy(), 0);
    info.detector_uz = select(mask, ray.uz(), 0);

    // Test we hit the detector array
    info.status = select(int_bvt(mask), PSTS_TS_OUTSIDE_DETECTOR, info.status);
    mask &= vcl::max(vcl::abs(info.detector_x), vcl::abs(info.detector_z)) <= pixel_array_halfwidth_;
#ifdef DEBUG_STATUS
    std::cout << ' ' << mask[0] << '/' << info.status[0];
#endif

    // Find the pixel ID
    int_vt pixel_j = select(int_bvt(mask), vcl::floor((info.detector_z + pixel_array_halfwidth_) * pixel_spacing_inv_), 0);
    int_vt pixel_i = select(int_bvt(mask), vcl::floor((info.detector_x + pixel_array_halfwidth_) * pixel_spacing_inv_), 0);
    info.pixel_id = select(int_bvt(mask), pixel_j * pixel_nside_ + pixel_i, -1);

    info.status = select(int_bvt(mask), PSTS_TS_FOUND_PIXEL, info.status);
#ifdef DEBUG_STATUS
    std::cout << ' ' << mask[0] << '/' << info.status[0];
#endif

    if(detector_has_rotation_) {
      ray.derotate(detector_rotation_.template cast<real_vt>());
    }
    ray.untranslate_origin(detector_origin_.template cast<real_vt>());

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
      radius = std::sqrt(lens_aperture2_);
    }
    if(distance <= 0) {
      throw std::runtime_error("Light emission distance must be positive or infinity");
    }

    double dcospolar = 1.0 - std::cos(radius / distance);
    double fppos_2y = 2.0 * detector_origin_.y();

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
      info.detector_x.store_a(xfp);
      info.detector_z.store_a(yfp);
      info.detector_t.store_a(tfp);

      for(unsigned i=0; i< VCLReal::num_real; i++) {
        ntraced++;
        if(status[i] >= PSTS_TS_OUTSIDE_DETECTOR) {
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

  const calin::math::spline_interpolation::CubicSpline& lens_refractive_index_spline_;

  vec3_t          scope_position_;
  real_t          az_rad_;
  real_t          el_rad_;
  real_t          phi_rad_;
  vec3_t          global_to_reflector_off_;
  mat3_t          global_to_reflector_rot_;
  real_t          air_ref_index_;

  vecX_t          lens_derivative_polynomial_;

  real_t          lens_aperture2_;

  real_t          detector_distance_;
  vec3_t          detector_origin_;
  bool            detector_has_rotation_;
  mat3_t          detector_rotation_;
  
  double          pixel_spacing_;
  double          pixel_spacing_inv_;
  unsigned        pixel_nside_;
  double          pixel_array_halfwidth_;
  
  real_t          scattering_sigma_theta_;

  RNG* rng_ = nullptr;
  bool adopt_rng_ = false;
};

#endif // SWIG

} } } // namespace calin::simulation::vcl_raytracer
`