/*

   calin/math/ray_vcl.hpp -- Stephen Fegan -- 2018-08-22

   Class describing light ray : with position, time, direction and energy.
   A simplification of the old UCLA "Particle" class which hopefully will
   be faster. VCL version.

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

#pragma once

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <calin_global_definitions.hpp>
#include <math/constants.hpp>
#include <util/vcl.hpp>
#include <math/ray.hpp>
#include <math/geometry_vcl.hpp>

namespace calin { namespace math { namespace ray {

#ifndef SWIG
using calin::math::constants::cgs_c;
#endif

template<typename VCLReal> class VCLRay: public VCLReal
{
public:
  using typename VCLReal::real_t;
  using typename VCLReal::real_vt;
  using typename VCLReal::bool_vt;
  using typename VCLReal::vec3_vt;
  using typename VCLReal::mat3_vt;

  VCLRay() { }
  VCLRay(const vec3_vt& pos, const vec3_vt& dir, const real_vt& time = 0,
      const real_vt& energy = 0):
    pos_(pos), dir_(dir), ct_(time*cgs_c), energy_(energy) {
    /* nothing to see here */
  }
  VCLRay(const Ray& ray) {
    pos_ = ray.position().cast<real_vt>();
    dir_ = ray.direction().cast<real_vt>();
    ct_ = ray.ct();
    energy_ = ray.energy();
  }
  Ray extract(unsigned i) const {
    Ray ray;
    ray.mutable_x() = x().extract(i);
    ray.mutable_y() = y().extract(i);
    ray.mutable_z() = z().extract(i);
    ray.mutable_ux() = ux().extract(i);
    ray.mutable_uy() = uy().extract(i);
    ray.mutable_uz() = uz().extract(i);
    ray.mutable_ct() = ct().extract(i);
    ray.mutable_energy() = energy().extract(i);
    return ray;
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

  real_vt ux_inv() const { calc_ux_inv(); return ux_inv_; }
  real_vt uy_inv() const { calc_uy_inv(); return uy_inv_; }
  real_vt uz_inv() const { calc_uz_inv(); return uz_inv_; }

  vec3_vt& mutable_position() { return pos_; }
  vec3_vt& mutable_direction() { clear_dir_inv(); return dir_; }
  real_vt& mutable_ct() { return ct_; }
  real_vt& mutable_energy() { return energy_; }
  real_vt& mutable_x() { return pos_.x(); }
  real_vt& mutable_y() { return pos_.y(); }
  real_vt& mutable_z() { return pos_.z(); }
  real_vt& mutable_ux() { has_ux_inv_=false; return dir_.x(); }
  real_vt& mutable_uy() { has_uy_inv_=false; return dir_.y(); }
  real_vt& mutable_uz() { has_uz_inv_=false; return dir_.z(); }

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

  bool_vt box_has_future_intersection(real_vt& tmin, real_vt& tmax,
    const vec3_vt& min_corner, const vec3_vt& max_corner) const
  {
    calc_dir_inv();
    return calin::math::geometry::VCL<VCLReal>::
      box_has_future_intersection_dirinv(tmin, tmax, min_corner, max_corner,
        pos_, ux_inv_, uy_inv_, uz_inv_);
  }

  void reflect_from_surface_with_mask(const bool_vt& mask, const vec3_vt& surface_norm) {
    clear_dir_inv();
    dir_ -= surface_norm * select(mask, 2.0*(dir_.dot(surface_norm)), 0);
  }

  // Refract at incoming surface (where eta=1/n<1 and norm.dir<0)
  void refract_at_surface_in_eta_with_mask(const bool_vt& mask, const vec3_vt& surface_norm, real_vt eta) {
    clear_dir_inv();
    eta = select(mask, eta, 1.0);
    const real_vt cosi = -dir_.dot(surface_norm);
    const real_vt etacosi = eta*cosi;
    // const real_vt c2 = 1.0-eta*eta*(1.0-cosi*cosi);
    // const real_vt c2 = nmul_add(eta*eta,nmul_add(cosi,cosi,1.0),1.0);
    const real_vt c2 = nmul_add(etacosi, etacosi, 1.0 - eta*eta);
    dir_ = eta*dir_ + (etacosi - sqrt(c2))*surface_norm;
  };

  // Refract at incoming surface (where n>1 and norm.dir<0)
  void refract_at_surface_in_with_mask(const bool_vt& mask, const vec3_vt& surface_norm, real_vt n) {
    refract_at_surface_in_eta(mask, surface_norm, 1.0/n);
  };

  // Refract at outgoing surface (where n>1 and norm.dir>0)
  bool_vt refract_at_surface_out_with_mask(const bool_vt& mask_in, const vec3_vt& surface_norm, real_vt n) {
    clear_dir_inv();
    real_vt eta = n;
    const real_vt cosi = dir_.dot(surface_norm);
    real_vt etacosi = n*cosi;
    const real_vt c2 = nmul_add(etacosi, etacosi, 1.0 - eta*eta);
    bool_vt mask = mask_in & (c2>0);
    eta = select(mask, eta, 1.0);
    etacosi = select(mask, etacosi, cosi);
    c2 = select(mask, c2, cosi*cosi);
    dir_ = eta*dir_ - (etacosi - sqrt(c2))*surface_norm;
    return mask;
  };

  //! Propagate ray fixed distance
  void propagate_dist(const real_vt& dist, const real_vt& n = 1.0) {
    pos_ += dir_ * dist;
    ct_  += dist * n;
  }

  //! Propagate ray fixed ct interval
  void propagate_ct(const real_vt& ct, const real_vt& n_inv = 1.0) {
    pos_ += dir_ * ct * n_inv;
    ct_  += ct;
  }

  //! Propagate ray fixed distance with mask (true if we are to propagate)
  void propagate_dist_with_mask(const bool_vt& mask, const real_vt& dist, const real_vt& n = 1.0) {
    real_vt masked_dist = select(mask, dist, 0);
    propagate_dist(masked_dist, n);
  }

  //! Propagate ray fixed distance with mask (true if we are to propagate)
  void propagate_ct_with_mask(const bool_vt& mask, const real_vt& ct, const real_vt& n_inv = 1.0) {
    real_vt masked_ct = select(mask, ct, 0);
    propagate_ct(masked_ct, n_inv);
  }

  //! Propagates free particle to plane with y-normal : y + d = 0
  bool_vt propagate_to_y_plane_with_mask(bool_vt mask, const real_vt& d,
    bool time_reversal_ok=true, const real_vt& n = 1.0)
  {
    calc_uy_inv();

    // Compute the distance between the plane and one parallel to it
    // which goes through the particle's current position
    real_vt plane_sep = -d - pos_.y();

    // Distance the particle must travel to reach the plane i.e.
    // n_hat * vhat = cos( theta ) = plane_dist / propagation_dist
    real_vt propagation_dist = plane_sep * uy_inv_;

    // Test whether particle is parallel to y plane
    mask &= is_finite(propagation_dist);

    // Test whether particles would need to travel backwards, if requested
    if(!time_reversal_ok)mask &= propagation_dist>=0;

    // Propagate particles that have mask=true
    propagate_dist_with_mask(mask, propagation_dist, n);

    return mask;
  }

  bool_vt propagate_to_y_plane(const real_vt& d,
    bool time_reversal_ok=true, const real_vt& n = 1.0)
  {
    bool_vt mask = true;
    return propagate_to_y_plane_with_mask(mask, d, time_reversal_ok, n);
  }

  //! Propagates free particle to plane with z-normal : z + d = 0
  bool_vt propagate_to_z_plane_with_mask(bool_vt mask, const real_vt& d,
    bool time_reversal_ok=true, const real_vt& n = 1.0)
  {
    calc_uz_inv();

    // Compute the distance between the plane and one parallel to it
    // which goes through the particle's current position
    real_vt plane_sep = -d - pos_.z();

    // Distance the particle must travel to reach the plane i.e.
    // n_hat * vhat = cos( theta ) = plane_dist / propagation_dist
    real_vt propagation_dist = plane_sep * uz_inv_;

    // Test whether particle is parallel to z plane
    mask &= is_finite(propagation_dist);

    // Test whether particles would need to travel backwards, if requested
    if(!time_reversal_ok)mask &= propagation_dist>=0;

    // Propagate particles that have mask=true
    propagate_dist_with_mask(mask, propagation_dist, n);

    return mask;
  }

  bool_vt propagate_to_z_plane(const real_vt& d,
    bool time_reversal_ok=true, const real_vt& n = 1.0)
  {
    bool_vt mask = true;
    return propagate_to_z_plane_with_mask(mask, d, time_reversal_ok, n);
  }

  bool_vt propagate_to_y_sphere_1st_interaction_fwd_bwd_with_mask(bool_vt mask,
    const real_vt& radius, const real_vt& surface_y_min = 0, const real_vt& n = 1.0)
  {
    vec3_vt pos_rel(pos_.x(), pos_.y()-(radius+surface_y_min), pos_.z());
    const real_vt b_2 = pos_rel.dot(dir_);
    const real_vt c = nmul_add(radius, radius, pos_rel.squaredNorm());

    const real_vt disc_4 = mul_sub(b_2, b_2, c);

    mask &= disc_4 >= 0;

    real_vt dist = -sqrt(disc_4) - b_2;
    propagate_dist_with_mask(mask, dist, n);
    return mask;
  }

  bool_vt propagate_to_y_sphere_1st_interaction_fwd_only_with_mask(bool_vt mask,
    const real_vt& radius, const real_vt& surface_y_min = 0, const real_vt& n = 1.0)
  {
    vec3_vt pos_rel(pos_.x(), pos_.y()-(radius+surface_y_min), pos_.z());
    const real_vt b_2 = pos_rel.dot(dir_);
    const real_vt c = nmul_add(radius, radius, pos_rel.squaredNorm());

    const real_vt disc_4 = mul_sub(b_2, b_2, c);
    mask &= disc_4 >= 0;

    real_vt dist = -sqrt(disc_4) - b_2;
    mask &= dist >= 0;

    propagate_dist_with_mask(mask, dist, n);
    return mask;
  }

  bool_vt propagate_to_y_sphere_2nd_interaction_fwd_bwd_with_mask(bool_vt mask,
    const real_vt& radius, const real_vt& surface_y_min = 0, const real_vt& n = 1.0)
  {
    vec3_vt pos_rel(pos_.x(), pos_.y()-(radius+surface_y_min), pos_.z());
    const real_vt b_2 = pos_rel.dot(dir_);
    const real_vt c = nmul_add(radius, radius, pos_rel.squaredNorm());

    const real_vt disc_4 = mul_sub(b_2, b_2, c);
    mask &= disc_4 >= 0;

    real_vt dist = sqrt(disc_4) - b_2;
    propagate_dist_with_mask(mask, dist, n);
    return mask;
  }

  bool_vt propagate_to_y_sphere_2nd_interaction_fwd_only_with_mask(bool_vt mask,
    const real_vt& radius, const real_vt& surface_y_min = 0, const real_vt& n = 1.0)
  {
    vec3_vt pos_rel(pos_.x(), pos_.y()-(radius+surface_y_min), pos_.z());
    const real_vt b_2 = pos_rel.dot(dir_);
    const real_vt c = nmul_add(radius, radius, pos_rel.squaredNorm());

    const real_vt disc_4 = mul_sub(b_2, b_2, c);
    mask &= disc_4 >= 0;

    real_vt dist = sqrt(disc_4) - b_2;
    mask &= dist >= 0;

    propagate_dist_with_mask(mask, dist, n);
    return mask;
  }

  bool_vt propagate_to_y_sphere_2nd_interaction_mostly_fwd_with_mask(bool_vt mask,
    const real_vt& radius, const real_vt& surface_y_min = 0, const real_vt& fwd_margin = 0.0,
    const real_vt& n = 1.0)
  {
    vec3_vt pos_rel(pos_.x(), pos_.y()-(radius+surface_y_min), pos_.z());
    const real_vt b_2 = pos_rel.dot(dir_);
    const real_vt c = nmul_add(radius, radius, pos_rel.squaredNorm());

    const real_vt disc_4 = mul_sub(b_2, b_2, c);
    mask &= disc_4 >= 0;

    real_vt dist = sqrt(disc_4) - b_2;
    mask &= dist >= fwd_margin;

    propagate_dist_with_mask(mask, dist, n);
    return mask;
  }

  //! Propagates free particle to the closest approach with point
  bool_vt propagate_to_point_closest_approach_with_mask(bool_vt mask,
    const vec3_vt& r0, bool time_reversal_ok, const real_vt& n = 1.0)
  {
    // Distance the particle must travel to reach the closest approach
    const real_vt dist = (r0-pos_).dot(dir_);

    if(not time_reversal_ok) {
      mask &= dist >= 0;
    }

    propagate_dist_with_mask(mask, dist, n);
    return true;
  }

  //! Propagates free particle to the closest approach with point with an
  // additional offset in distance (positive means particle overruns the
  // closest point, negative meens it doesn't go far enough to reach it)
  bool_vt propagate_to_point_nearly_closest_approach_with_mask(bool_vt mask,
    const vec3_vt& r0, const real_vt& overrun, bool time_reversal_ok, const real_vt& n = 1.0)
  {
    // Distance the particle must travel to reach the closest approach + overrun
    const real_vt dist = (r0-pos_).dot(dir_) + overrun;

    if(not time_reversal_ok) {
      mask &= dist >= 0;
    }

    propagate_dist_with_mask(mask, dist, n);
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

#ifndef SWIG
  static void* operator new(size_t nbytes) {
    return VCLReal::aligned_malloc(nbytes);
  }
  static void* operator new(size_t nbytes, void* p) {
    return p;
  }
  static void operator delete(void *p) {
    VCLReal::aligned_free(p);
  }
#endif

private:
  void calc_ux_inv() const {
    if(!has_ux_inv_) { has_ux_inv_ = true; ux_inv_ = 1.0/dir_.x(); } }
  void calc_uy_inv() const {
    if(!has_uy_inv_) { has_uy_inv_ = true; uy_inv_ = 1.0/dir_.y(); } }
  void calc_uz_inv() const {
    if(!has_uz_inv_) { has_uz_inv_ = true; uz_inv_ = 1.0/dir_.z(); } }
  void calc_dir_inv() const { calc_ux_inv(); calc_uy_inv(); calc_uz_inv(); }

  void clear_dir_inv() const { has_ux_inv_ = has_uy_inv_ = has_uz_inv_ = false; }

  vec3_vt pos_ = vec3_vt::Zero();
  vec3_vt dir_ = vec3_vt::UnitY();
  real_vt ct_ = 0;
  real_vt energy_ = 0;
  mutable bool has_ux_inv_ = false;
  mutable bool has_uy_inv_ = false;
  mutable bool has_uz_inv_ = false;
  mutable real_vt ux_inv_ = 0;
  mutable real_vt uy_inv_ = 0;
  mutable real_vt uz_inv_ = 0;
};

template<typename VCLReal> class VCLRayArray: public VCLReal
{
public:
  using typename VCLReal::real_t;
  using typename VCLReal::real_vt;
  using typename VCLReal::real_at;
  using typename VCLReal::vec3_vt;

  VCLRayArray() { }

  VCLRayArray(const VCLRayArray<VCLReal>& ray) {
    set_rays(ray);
  }

  VCLRayArray(const vec3_vt& pos, const vec3_vt& dir, const real_vt& ct = 0,
      const real_vt& energy = 0) {
    set_positions(pos);
    set_directions(dir);
    set_cts(ct);
    set_energies(energy);
  }

  const real_at& x() const { return x_; }
  const real_at& y() const { return y_; }
  const real_at& z() const { return z_; }
  const real_at& ux() const { return ux_; }
  const real_at& uy() const { return uy_; }
  const real_at& uz() const { return uz_; }
  const real_at& ct() const { return ct_; }
  const real_at& energy() const { return energy_; }

  real_at& mutable_x() { return x_; }
  real_at& mutable_y() { return y_; }
  real_at& mutable_z() { return z_; }
  real_at& mutable_ux() { return ux_; }
  real_at& mutable_uy() { return uy_; }
  real_at& mutable_uz() { return uz_; }
  real_at& mutable_ct() { return ct_; }
  real_at& mutable_energy() { return energy_; }

  void set_rays(const VCLRayArray<VCLReal>& ray) {
    set_positions(ray.position());
    set_directions(ray.direction());
    set_cts(ray.ct());
    set_energies(ray.energy());
  }

  void get_rays(VCLRayArray<VCLReal>& ray) const {
    get_positions(ray.mutable_position());
    get_directions(ray.mutable_direction());
    get_cts(ray.mutable_ct());
    get_energies(ray.mutable_energy());
  }

  VCLRayArray<VCLReal> get_rays() const {
    VCLRayArray<VCLReal> ray;
    get_rays(ray);
    return ray;
  }

  void extract_one_ray(unsigned i, Ray& ray) const {
    ray.mutable_x() = x_[i];
    ray.mutable_y() = y_[i];
    ray.mutable_z() = z_[i];
    ray.mutable_ux() = ux_[i];
    ray.mutable_uy() = uy_[i];
    ray.mutable_uz() = uz_[i];
    ray.mutable_ct() = ct_[i];
    ray.mutable_energy() = energy_[i];
  }

  Ray extract_one_ray(unsigned i) const {
    VCLRayArray<VCLReal> ray;
    extract_ray(i, ray);
    return ray;
  }

  void insert_one_ray(unsigned i, const Ray& ray) {
    x_[i] = ray.x();
    y_[i] = ray.y();
    z_[i] = ray.z();
    ux_[i] = ray.ux();
    uy_[i] = ray.uy();
    uz_[i] = ray.uz();
    ct_[i] = ray.ct();
    energy_[i] = ray.energy();
  }

  void get_positions(vec3_vt& pos) const { pos.x().load(x_); pos.y().load(y_); pos.z().load(z_); }
  void get_directions(vec3_vt& dir) const { dir.x().load(ux_); dir.y().load(uy_); dir.z().load(uz_); }
  void get_cts(real_vt& ct) const { ct.load(ct_); }
  void get_energies(real_vt& energy) const { energy.load(energy_); }

  vec3_vt get_positions() const { vec3_vt pos; get_positions(pos); return pos; }
  vec3_vt get_directions() const { vec3_vt dir; get_directions(dir); return dir; }
  real_vt get_cts() const { real_vt ct; get_cts(ct); return ct; }
  real_vt get_energies() const { real_vt energy; get_energies(energy); return energy; }

  void set_positions(const vec3_vt& pos) { pos.x().store(x_); pos.y().store(y_); pos.z().store(z_); }
  void set_directions(const vec3_vt& dir) { dir.x().store(ux_); dir.y().store(uy_); dir.z().store(z_); }
  void set_cts(const real_vt& ct) { ct.store(ct_); }
  void set_energies(const real_vt& energy) { energy.store(energy_); }

#ifndef SWIG
  static void* operator new(size_t nbytes) {
    return VCLReal::aligned_malloc(nbytes);
  }
  static void* operator new(size_t nbytes, void* p) {
    return p;
  }
  static void operator delete(void *p) {
    VCLReal::aligned_free(p);
  }
#endif

private:
  real_at x_;
  real_at y_;
  real_at z_;
  real_at ux_;
  real_at uy_;
  real_at uz_;
  real_at ct_;
  real_at energy_;
};

} } } // namespace calin::math::ray
