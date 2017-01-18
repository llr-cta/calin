/*

   calin/math/position_generator.cpp -- Stephen Fegan -- 2017-01-18

   Geanerate positions in space using some algorithm.

   Copyright 2017, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <cmath>
#include <math/position_generator.hpp>
#include <math/hex_array.hpp>
#include <math/special.hpp>

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
  if(scale_weight_by_area and nray>0)weight_ *= 4.0*M_PI/double(nray);
}

MCPlanePositionGenerator::~MCPlanePositionGenerator()
{
  // nothing to see here
}

void MCPlanePositionGenerator::reset()
{
  iray_ = 0;
}

bool MCPlanePositionGenerator::next(Eigen::Vector3d& x, double& weight)
{
  if(iray_ >= nray_ and nray > 0)return false;
  double r = std::sqrt(rng_.uniform() * r2_max_);
  double phi = rng_.uniform() * 2.0*M_PI;
  x.x() = r*cos(phi);
  x.y() = r*sin(phi);
  x.z() = 0;
  weight = weight_;
  ++iray_;
  return true;
}

HexGridPlanePositionGenerator::
HexGridPlanePositionGenerator(double r_max, double dx,
    bool scale_weight_by_area, double base_weight):
  PositionGenerator(),
  r2_max_(SQR(r_max)), weight_(base_weight), dx_(dx)
{
  if(scale_weight_by_area)weight_ *= calin::math::hex_array::cell_area(dx);
}

HexGridPlanePositionGenerator::~HexGridPlanePositionGenerator()
{
  // nothing to see here
}

void HexGridPlanePositionGenerator::reset()
{
  hex_id = 0;
}

bool HexGridPlanePositionGenerator::next(Eigen::Vector3d& x, double& weight)
{

  weight = weight_;
  ++iray_;
  return true;
}
