/*

   calin/math/ray_generator.cpp -- Stephen Fegan -- 2017-01-20

   Geanerate rays in space using some algorithm.

   Copyright 2017, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <math/geometry.hpp>
#include <math/ray_generator.hpp>

#include <util/log.hpp>
using namespace calin::util::log;

using namespace calin::math::ray_generator;

RayGenerator::~RayGenerator()
{
  // nothing to see here
}

PositionNestedByDirectionRayGenerator::
PositionNestedByDirectionRayGenerator(const Eigen::Vector3d& shift_post_rot,
    const Eigen::Matrix3d& rot, const Eigen::Vector3d& shift_pre_rot,
    calin::math::direction_generator::DirectionGenerator* dir_gen,
    calin::math::position_generator::PositionGenerator* pos_gen,
    bool adopt_dir_gen, bool adopt_pos_gen):
  RayGenerator(), base_rot_(rot), base_shift_post_rot_(shift_post_rot),
  base_shift_pre_rot_(shift_pre_rot),
  dir_gen_(dir_gen), adopt_dir_gen_(adopt_dir_gen),
  pos_gen_(pos_gen), adopt_pos_gen_(adopt_pos_gen)
{
  // nothing to see here
}

PositionNestedByDirectionRayGenerator::~PositionNestedByDirectionRayGenerator()
{
  if(adopt_pos_gen_)delete pos_gen_;
  if(adopt_dir_gen_)delete dir_gen_;
}

void PositionNestedByDirectionRayGenerator::reset()
{
  has_dir_ = false;
  pos_gen_->reset();
  dir_gen_->reset();
}

bool PositionNestedByDirectionRayGenerator::
next(calin::math::ray::Ray& ray, double& weight)
{
do_over:
  if(not has_dir_)has_dir_ = compute_rot_and_shift();
  if(not has_dir_)return false;
  if(not pos_gen_->next(ray.position(), weight))
  {
    has_dir_ = false;
    pos_gen_->reset();
    goto do_over;
  }
  weight *= dir_weight_;
  ray.position()  = rot_ * ray.position() + shift_post_rot_;
  ray.direction() = rot_.col(2);
  ray.ct() = 0;
  return true;
}

bool PositionNestedByDirectionRayGenerator::compute_rot_and_shift()
{
  Eigen::Vector3d dir;
  if(not dir_gen_->next_as_vector(dir, dir_weight_))return false;
  rot_ = base_rot_ * calin::math::geometry::rotation_z_to_vec(dir);
  //LOG(INFO) << rot_;
  shift_post_rot_ = base_shift_post_rot_ + rot_*base_shift_pre_rot_;
  return true;
}
