/*

   calin/math/ray.hpp -- Stephen Fegan -- 2016-10-23

   Class describing light ray : with position, time, direction and energy.
   A simplification of the old UCLA "Particle" class which hopefully will
   be faster.

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, Ecole polytechnique, CNRS/IN2P3, Universite Paris-Saclay

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

namespace calin { namespace math { namespace ray {

using calin::math::constants::cgs_c;

class Ray
{
public:
  Ray() { /* nothing to see here */ }
  Ray(const Eigen::Vector3d& pos, const Eigen::Vector3d& dir, double time = 0,
      double energy = 0): pos_(pos), dir_(dir), ct_(time*cgs_c), energy_(energy) {
    /* nothing to see here */ }

  Eigen::Vector3d& position() { return pos_; }
  Eigen::Vector3d& direction() { return dir_; }
  double& ct() { return ct_; }
  double& energy() { return energy_; }

  const Eigen::Vector3d& position() const { return pos_; }
  const Eigen::Vector3d& direction() const { return dir_; }
  double ct() const { return ct_; }
  double time() const { return ct_/cgs_c; }
  double energy() const { return energy_; }

  void translate_origin(const Eigen::Vector3d& origin) { pos_ -= origin; }
  void untranslate_origin(const Eigen::Vector3d& origin) { pos_ += origin; }

  void rotate(const Eigen::Matrix3d& rot) {
    pos_ = rot * pos_; dir_ = rot * dir_; }
  void derotate(const Eigen::Matrix3d& rot) {
    pos_ = rot.transpose() * pos_; dir_ = rot.transpose() * dir_; }

  //! Propagate free particle fixed distance
  inline void propagate_dist(double dist, double n = 1.0) {
    pos_ += dir_ * dist;
    ct_  += dist * n;
  }

  //! Propagates free particle to the given plane
  bool propagate_to_plane(const Eigen::Vector3d& normal, double d,
    bool time_reversal_ok=true, double n = 1.0);

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

  IPOut propagate_to_sphere(const Eigen::Vector3d& center, double radius,
    IntersectionPoint ip = IP_CLOSEST, bool time_reversal_ok = true,
    double n = 1.0);

  IPOut propagate_to_cylinder(const Eigen::Vector3d& center,
    const Eigen::Vector3d& normal, double radius,
    IntersectionPoint ip = IP_CLOSEST, bool time_reversal_ok = true, double n = 1.0);

private:
  Eigen::Vector3d pos_ = Eigen::Vector3d::Zero();
  Eigen::Vector3d dir_ = Eigen::Vector3d::Zero();
  double ct_ = 0;
  double energy_ = 0;
};

} } } // namespace calin::math::ray
