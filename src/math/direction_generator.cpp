/*

   calin/math/direction_generator.cpp -- Stephen Fegan -- 2017-01-19

   Geanerate directions in space using some algorithm.

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
#include <limits>
#include <math/direction_generator.hpp>
#include <math/healpix_array.hpp>
#include <math/special.hpp>

using namespace calin::math::direction_generator;
using calin::math::special::SQR;

DirectionGenerator::~DirectionGenerator()
{
  // nothing to see here
}

bool DirectionGenerator::next_as_vector(Eigen::Vector3d& dir, double& weight)
{
  double theta;
  double phi;
  if(this->next_as_theta_phi(theta,phi,weight))
  {
    double sin_theta = std::sin(theta);
    dir.x() = sin_theta * std::cos(phi);
    dir.y() = sin_theta * std::sin(phi);
    dir.z() = std::cos(theta);
    return true;
  }
  return false;
}

bool DirectionGenerator::next_as_matrix(Eigen::Matrix3d& trans_mat, double& weight)
{
  double theta;
  double phi;
  if(this->next_as_theta_phi(theta,phi,weight))
  {
    trans_mat =
      Eigen::AngleAxisd(theta, Eigen::Vector3d::UnitX()) *
      Eigen::AngleAxisd(phi,   Eigen::Vector3d::UnitZ());
    return true;
  }
  return false;
}

MCSphereDirectionGenerator::
MCSphereDirectionGenerator(double theta_max, unsigned nray,
    calin::math::rng::RNG* rng, bool scale_weight_by_area, double base_weight,
    bool adopt_rng):
  DirectionGenerator(), nray_(nray), cos_theta_max_(std::cos(theta_max)),
  weight_(base_weight), rng_(rng), adopt_rng_(adopt_rng)
{
  if(scale_weight_by_area)weight_ *= 2*M_PI*(1.0-cos_theta_max_)/double(nray);
}

MCSphereDirectionGenerator::~MCSphereDirectionGenerator()
{
  if(adopt_rng_)delete rng_;
}

void MCSphereDirectionGenerator::reset()
{
  iray_ = 0;
}

bool MCSphereDirectionGenerator::
next_as_theta_phi(double& theta, double& phi, double& weight)
{
  if(iray_ == nray_)return false;
  theta = std::acos(1.0 - rng_->uniform()*(1.0-cos_theta_max_));
  phi = rng_->uniform() * 2.0*M_PI;
  weight = weight_;
  ++iray_;
  return true;
}

#if 0
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
#endif
