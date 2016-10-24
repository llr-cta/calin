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

namespace calin { namespace math { namespace ray {

class Ray
{
public:
  Ray() { /* nothing to see here */ }
  Ray(const Eigen::Vector3d& pos, const Eigen::Vector3d& dir, double time = 0,
      double energy = 0): pos_(pos), dir_(dir), time_(time), energy_(energy) {
    /* nothing to see here */ }

  Eigen::Vector3d& position() { return pos_; }
  Eigen::Vector3d& direction() { return dir_; }
  double& time() { return time_; }
  double& energy() { return energy_; }

  const Eigen::Vector3d& position() const { return pos_; }
  const Eigen::Vector3d& direction() const { return dir_; }
  double time() const { return time_; }
  double energy() const { return energy_; }

  void translate_origin(const Eigen::Vector3d& origin) { pos_ -= origin; }
  void untranslate_origin(const Eigen::Vector3d& origin) { pos_ += origin; }

  void rotate(const Eigen::Matrix3d& rot) {
    pos_ = rot * pos_; dir_ = rot * dir_; }
  void derotate(const Eigen::Matrix3d& rot) {
    pos_ = rot.transpose() * pos_; dir_ = rot.transpose() * dir_; }

private:
  Eigen::Vector3d pos_ = Eigen::Vector3d::Zero();
  Eigen::Vector3d dir_ = Eigen::Vector3d::Zero();
  double time_ = 0;
  double energy_ = 0;
};

} } } // namespace calin::math::ray
