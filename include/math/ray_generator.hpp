/*

   calin/math/ray_generator.hpp -- Stephen Fegan -- 2017-01-20

   Geanerate rays in space using some algorithm.

   Copyright 2017, Stephen Fegan <sfegan@llr.in2p3.fr>
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
#include <math/ray.hpp>
#include <math/position_generator.hpp>
#include <math/direction_generator.hpp>

namespace calin { namespace math { namespace ray_generator {

class RayGenerator
{
public:
  virtual ~RayGenerator();
  virtual void reset() = 0;
  virtual bool next(calin::math::ray::Ray& ray, double& weight) = 0;
};

class PositionNestedByDirectionRayGenerator: public RayGenerator
{
public:
  PositionNestedByDirectionRayGenerator(const Eigen::Vector3d& shift_post_rot,
    const Eigen::Matrix3d& rot, const Eigen::Vector3d& shift_pre_rot,
    calin::math::direction_generator::DirectionGenerator* dir_gen,
    calin::math::position_generator::PositionGenerator* pos_gen,
    bool adopt_dir_gen = false, bool adopt_pos_gen = false);
  virtual ~PositionNestedByDirectionRayGenerator();
  void reset() override;
  bool next(calin::math::ray::Ray& ray, double& weight) override;
protected:
  bool compute_rot_and_shift();
  bool has_dir_ = false;
  Eigen::Matrix3d rot_ = Eigen::Matrix3d::Identity();
  Eigen::Vector3d shift_post_rot_ = Eigen::Vector3d::Zero();
  double dir_weight_ = 1;
  Eigen::Matrix3d base_rot_ = Eigen::Matrix3d::Identity();
  Eigen::Vector3d base_shift_post_rot_ = Eigen::Vector3d::Zero();
  Eigen::Vector3d base_shift_pre_rot_ = Eigen::Vector3d::Zero();
  calin::math::direction_generator::DirectionGenerator* dir_gen_ = nullptr;
  bool adopt_dir_gen_ = false;
  calin::math::position_generator::PositionGenerator* pos_gen_ = nullptr;
  bool adopt_pos_gen_ = false;
};

} } } // namespace calin::math::position_generator
