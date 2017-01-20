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
#include <math/geometry.hpp>

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
      Eigen::AngleAxisd(phi,   Eigen::Vector3d::UnitZ()) *
      Eigen::AngleAxisd(theta, Eigen::Vector3d::UnitY());

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

HEALPixDirectionGenerator::
HEALPixDirectionGenerator(double theta_max, unsigned nside,
    bool scale_weight_by_area, double base_weight):
  DirectionGenerator(), cos_theta_max_(std::cos(theta_max)), nside_(nside),
  weight_(base_weight)
{
  if(scale_weight_by_area)
    weight_ *= calin::math::healpix_array::cell_area(nside);
}

HEALPixDirectionGenerator::~HEALPixDirectionGenerator()
{
  // nothing to see here
}

void HEALPixDirectionGenerator::reset()
{
  pixid_ = 0;
}

bool HEALPixDirectionGenerator::
next_as_theta_phi(double& theta, double& phi, double& weight)
{
  double x,y,z;
  calin::math::healpix_array::pixid_to_xyz(nside_, pixid_, x, y, z);
  if(z < cos_theta_max_)return false;
  ++pixid_;
  phi = std::atan2(y,x);
  theta = std::atan2(std::sqrt(x*x+y*y),z);
  weight = weight_;
  return true;
}

bool HEALPixDirectionGenerator::
next_as_vector(Eigen::Vector3d& dir, double& weight)
{
  double x,y,z;
  calin::math::healpix_array::pixid_to_xyz(nside_, pixid_, x, y, z);
  if(z < cos_theta_max_)return false;
  ++pixid_;
  dir << x, y, z;
  weight = weight_;
  return true;
}

bool HEALPixDirectionGenerator::
next_as_matrix(Eigen::Matrix3d& trans_mat, double& weight)
{
  double x,y,z;
  calin::math::healpix_array::pixid_to_xyz(nside_, pixid_, x, y, z);
  if(z < cos_theta_max_)return false;
  ++pixid_;
  calin::math::geometry::rotation_z_to_xyz(trans_mat, x, y, z);
  weight = weight_;
  return true;
}
