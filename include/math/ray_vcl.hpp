/*

   calin/math/ray_vcl.hpp -- Stephen Fegan -- 2018-08-22

   Class describing light ray : with position, time, direction and energy.
   A simplification of the old UCLA "Particle" class which hopefully will
   be faster. VCL version.

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

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <calin_global_definitions.hpp>
#include <math/constants.hpp>
#include <util/vcl.hpp>

namespace calin { namespace math { namespace ray {

#ifndef SWIG
using calin::math::constants::cgs_c;
#endif

template<typename VCLReal> class VCLRay: public VCLReal
{
public:
  using typename VCLReal::real_t;
  using typename VCLReal::real_vt;
  using typename VCLReal::real_bvt;
  using typename VCLReal::vec3_vt;
  using typename VCLReal::mat3_vt;

  VCLRay() { }
  VCLRay(const vec3_vt& pos, const vec3_vt& dir, const real_vt& time = 0,
      const real_vt& energy = 0):
    pos_(pos), dir_(dir), ct_(time*cgs_c), energy_(energy) {
    /* nothing to see here */
  }

  const vec3_vt& position() const { return pos_; }
  const vec3_vt& direction() const { return dir_; }
  const real_vt& ct() const { return ct_; }
  const real_vt& time() const { return ct_/cgs_c; }
  const real_vt& energy() const { return energy_; }
  const real_vt& x() const { return pos_.x(); }
  const real_vt& y() const { return pos_.y(); }
  const real_vt& z() const { return pos_.z(); }
  const real_vt& ux() const { return dir_.x(); }
  const real_vt& uy() const { return dir_.y(); }
  const real_vt& uz() const { return dir_.z(); }

  void set_position(const vec3_vt& pos) { pos_ = pos; }
  void set_direction(const vec3_vt& dir) { clear_dir_inv(); dir_ = dir; }
  void set_ct(const real_vt& ct) {  ct_ = ct; }
  void set_time(const real_vt& t) { ct_ = t*cgs_c; }
  void set_energy(const real_vt& e) { energy_ = e; }

  void translate_origin(const vec3_vt& origin) { pos_ -= origin; }
  void untranslate_origin(const vec3_vt& origin) { pos_ += origin; }

  void rotate(const mat3_vt& rot) {
    clear_dir_inv(); pos_ = rot * pos_; dir_ = rot * dir_; }
  void derotate(const mat3_vt& rot) {
    clear_dir_inv(); pos_ = rot.transpose() * pos_; dir_ = rot.transpose() * dir_; }

  void reflect_from_surface(const vec3_vt& surface_norm) {
    clear_dir_inv();
    dir_ -= surface_norm * (2.0*(dir_.dot(surface_norm)));
  }

  // Refract at incoming surface (where eta=1/n<1 and norm.dir<0)
  void refract_at_surface_in_eta(const vec3_vt& surface_norm, real_vt eta) {
    clear_dir_inv();
    real_vt cosi = -dir_.dot(surface_norm);
    real_vt c2 = 1.0-eta*eta*(1.0-cosi*cosi);
    dir_ = eta*dir_ + (eta*cosi - sqrt(c2))*surface_norm;
  };

  // Refract at incoming surface (where n>1 and norm.dir<0)
  void refract_at_surface_in(const vec3_vt& surface_norm, real_vt n) {
    refract_at_surface_in_eta(surface_norm, 1.0/n);
  };

  // Refract at outgoing surface (where n>1 and norm.dir>0)
  real_bvt refract_at_surface_out(const vec3_vt& surface_norm, real_vt n) {
    clear_dir_inv();
    const real_vt eta = n;
    const real_vt cosi = dir_.dot(surface_norm);
    const real_vt c2 = 1.0-eta*eta*(1.0-cosi*cosi);
    real_bvt mask = c2>0;
    dir_ = eta*dir_ - (eta*cosi - sqrt(select(mask,c2,0)))*surface_norm;
    return mask;
  };

  //! Propagate ray fixed distance
  void propagate_dist(const real_vt& dist, const real_vt& n = 1.0) {
    pos_ += dir_ * dist;
    ct_  += dist * n;
  }

  //! Propagate ray fixed distance with mask (true if we are to propagate)
  void propagate_dist_with_mask(const real_bvt& mask, const real_vt& dist, const real_vt& n = 1.0) {
    real_vt masked_dist = select(mask, dist, 0);
    propagate_dist(masked_dist, n);
  }

  //! Propagates free particle to plane with y-normal
  real_bvt propagate_to_y_plane(real_vt d, bool time_reversal_ok=true, real_vt n = 1.0)
  {
    calc_uy_inv();

    // Compute the distance between the plane and one parallel to it
    // which goes through the particle's current position
    real_vt plane_sep = -d - pos_.y();

    // Distance the particle must travel to reach the plane i.e.
    // n_hat * vhat = cos( theta ) = plane_dist / propagation_dist
    real_vt propagation_dist = plane_sep * uy_inv_;

    // Test whether particle is parallel to y plane
    real_bvt mask = is_finite(propagation_dist);

    // Test whether particles would need to travel backwards, if requested
    if(!time_reversal_ok)mask &= propagation_dist>=0;

    // Propagate particles that have mask=true
    propagate_dist_with_mask(mask, propagation_dist, n);

    return true;
  }

#if 0
  //! Propagates free particle to the given plane
  bool propagate_to_plane(const Eigen::Vector3d& normal, real_vt d,
    bool time_reversal_ok=true, real_vt n = 1.0);


  //! Propagates free particle to the closest approach with line
  bool propagate_to_point_closest_approach(const Eigen::Vector3d& r0,
    bool time_reversal_ok=true, double n = 1.0);

  //! Propagates free particle to the closest approach with line
  bool propagate_to_line_closest_approach(const Eigen::Vector3d& normal,
    const Eigen::Vector3d& r0, bool time_reversal_ok=true, double n = 1.0);

  //! Propagates particle to a given sphere
  enum IntersectionPoint { IP_CLOSEST, IP_FARTHEST,
                           IP_NEXT, IP_PREVIOUS,
                           IP_EARLIEST, IP_LATEST };

  enum IPOut { IPO_NONE, IPO_FIRST, IPO_SECOND };

  bool propagate_to_standard_sphere_1st_interaction_fwd_bwd(double radius,
    double n = 1.0);

  bool propagate_to_standard_sphere_1st_interaction_fwd_only(double radius,
    double n = 1.0);

  bool propagate_to_y_sphere_1st_interaction_fwd_bwd(double radius,
    double surface_y_min = 0, double n = 1.0);

  bool propagate_to_y_sphere_1st_interaction_fwd_only(double radius,
    double surface_y_min = 0, double n = 1.0);

  bool propagate_to_standard_sphere_2nd_interaction_fwd_bwd(double radius,
    double n = 1.0);

  bool propagate_to_standard_sphere_2nd_interaction_fwd_only(double radius,
    double n = 1.0);

  bool propagate_to_y_sphere_2nd_interaction_fwd_bwd(double radius,
    double surface_y_min = 0, double n = 1.0);

  bool propagate_to_y_sphere_2nd_interaction_fwd_only(double radius,
    double surface_y_min = 0, double n = 1.0);

  bool propagate_to_standard_sphere_2nd_interaction(double radius,
    bool time_reversal_ok = true, double n = 1.0)
  {
    if(time_reversal_ok)
      return propagate_to_standard_sphere_2nd_interaction_fwd_bwd(radius, n);
    else
      return propagate_to_standard_sphere_2nd_interaction_fwd_only(radius, n);
  }

  IPOut propagate_to_sphere(const Eigen::Vector3d& center, double radius,
    IntersectionPoint ip = IP_CLOSEST, bool time_reversal_ok = true,
    double n = 1.0);

  IPOut propagate_to_cylinder(const Eigen::Vector3d& center,
    const Eigen::Vector3d& normal, double radius,
    IntersectionPoint ip = IP_CLOSEST, bool time_reversal_ok = true, double n = 1.0);
#endif

  void calc_ux_inv() {
    if(!has_ux_inv_) { has_ux_inv_ = true; ux_inv_ = 1.0/dir_.x(); } }
  void calc_uy_inv() {
    if(!has_uy_inv_) { has_uy_inv_ = true; uy_inv_ = 1.0/dir_.y(); } }
  void calc_uz_inv() {
    if(!has_uz_inv_) { has_uz_inv_ = true; uz_inv_ = 1.0/dir_.z(); } }
  void calc_dir_inv() { calc_ux_inv(); calc_uy_inv(); calc_uz_inv(); }

private:
  void clear_dir_inv() { has_ux_inv_ = has_uy_inv_ = has_uz_inv_ = false; }

  vec3_vt pos_ = vec3_vt::Zero();
  vec3_vt dir_ = vec3_vt::UnitY();
  real_vt ct_ = 0;
  real_vt energy_ = 0;
  bool has_ux_inv_ = false;
  bool has_uy_inv_ = false;
  bool has_uz_inv_ = false;
  real_vt ux_inv_ = 0;
  real_vt uy_inv_ = 0;
  real_vt uz_inv_ = 0;
};

} } } // namespace calin::math::ray
