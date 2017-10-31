/*

   calin/math/position_generator.cpp -- Stephen Fegan -- 2017-01-18

   Geanerate positions in space using some algorithm.

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

#include <cmath>
#include <limits>
#include <math/position_generator.hpp>
#include <math/hex_array.hpp>
#include <math/special.hpp>

using namespace calin::math::position_generator;
using calin::math::special::SQR;

PositionGenerator::~PositionGenerator()
{
  // nothing to see here
}

MCPlanePositionGenerator::MCPlanePositionGenerator(double r_max, unsigned nray,
    calin::math::rng::RNG* rng, bool scale_weight_by_area, double base_weight,
    bool adopt_rng):
  PositionGenerator(),
  nray_(nray), r2_max_(SQR(r_max)), weight_(base_weight),
  rng_(rng), adopt_rng_(adopt_rng)
{
  if(scale_weight_by_area)weight_ *= M_PI*SQR(r_max)/double(nray);
}

MCPlanePositionGenerator::~MCPlanePositionGenerator()
{
  if(adopt_rng_)delete rng_;
}

void MCPlanePositionGenerator::reset()
{
  iray_ = 0;
}

bool MCPlanePositionGenerator::next(Eigen::Vector3d& pos, double& weight)
{
  if(iray_ == nray_)return false;
  double r = std::sqrt(rng_->uniform() * r2_max_);
  double phi = rng_->uniform() * 2.0*M_PI;
  pos.x() = r*cos(phi);
  pos.y() = r*sin(phi);
  pos.z() = 0;
  weight = weight_;
  ++iray_;
  return true;
}

HexGridPlanePositionGenerator::
HexGridPlanePositionGenerator(double r_max, double dx,
    bool scale_weight_by_area, double base_weight):
  PositionGenerator(),
  r2_max_(SQR(r_max)), dx_(dx), weight_(base_weight)
{
  if(scale_weight_by_area)weight_ *= calin::math::hex_array::cell_area(dx);
}

HexGridPlanePositionGenerator::~HexGridPlanePositionGenerator()
{
  // nothing to see here
}

void HexGridPlanePositionGenerator::reset()
{
  hexid_ = 0;
}

bool HexGridPlanePositionGenerator::next(Eigen::Vector3d& pos, double& weight)
{
  double r2_last = 0;
  bool increasing = true;
  unsigned hexid = hexid_;
  while(true)
  {
    double x;
    double y;
    calin::math::hex_array::hexid_to_xy(hexid++, x, y);
    x *= dx_;
    y *= dx_;
    double r2 = x*x+y*y;
    if(r2 <= r2_max_) {
      pos.x() = x;
      pos.y() = y;
      pos.z() = 0;
      weight = weight_;
      hexid_ = hexid;
      return true;
    } else if(r2 >= r2_last) {
      if(!increasing)return false;
    } else {
      increasing = false;
    }
    r2_last = r2;
    /* try again at next hex site! */
  }
  return false;
}

TransformedPositionGenerator::
TransformedPositionGenerator(const Eigen::Vector3d& shift_post_rot,
    const Eigen::Matrix3d& rot, const Eigen::Vector3d& shift_pre_rot,
    PositionGenerator* gen, bool adopt_gen):
  PositionGenerator(), gen_(gen), adopt_gen_(adopt_gen), rot_(rot),
  shift_post_rot_(rot*shift_pre_rot + shift_post_rot)
{
  // nothing to see here
}

TransformedPositionGenerator::
TransformedPositionGenerator(const Eigen::Vector3d& shift_post_rot,
    const Eigen::Matrix3d& rot, PositionGenerator* gen, bool adopt_gen):
  TransformedPositionGenerator(shift_post_rot, rot, Eigen::Vector3d::Zero(),
    gen, adopt_gen)
{
  // nothing to see here
}

TransformedPositionGenerator::~TransformedPositionGenerator()
{
  if(adopt_gen_)delete gen_;
}

void TransformedPositionGenerator::reset()
{
  gen_->reset();
}

bool TransformedPositionGenerator::next(Eigen::Vector3d& pos, double& weight)
{
  if(gen_->next(pos, weight)) {
    pos = rot_ * pos + shift_post_rot_;
    return true;
  } else {
    return false;
  }
}
