/*

   calin/math/direction_generator.cpp -- Stephen Fegan -- 2017-01-19

   Geanerate directions in space using some algorithm.

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

SingleDirectionGenerator::
SingleDirectionGenerator(double theta, double phi, double base_weight):
  DirectionGenerator(), theta_(theta), phi_(phi), weight_(base_weight)
{
  // nothing to see here
}

SingleDirectionGenerator::
SingleDirectionGenerator(const Eigen::Vector3d& dir, double base_weight):
  DirectionGenerator(),
  theta_(std::atan2(std::sqrt(SQR(dir.x())+SQR(dir.y())),dir.z())),
  phi_(std::atan2(dir.y(),dir.x())), weight_(base_weight)
{
  // nothing to see here
}

SingleDirectionGenerator::~SingleDirectionGenerator()
{
  // nothing to see here
}

void SingleDirectionGenerator::reset()
{
  direction_generated_ = false;
}

bool SingleDirectionGenerator::
next_as_theta_phi(double& theta, double& phi, double& weight)
{
  if(direction_generated_)return false;
  theta = theta_;
  phi = phi_;
  weight = weight_;
  direction_generated_ = true;
  return true;
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

TransformedDirectionGenerator::
TransformedDirectionGenerator(const Eigen::Matrix3d& rot,
    DirectionGenerator* gen, bool adopt_gen):
  DirectionGenerator(), gen_(gen), adopt_gen_(adopt_gen), rot_(rot)
{
  // nothing to see here
}

TransformedDirectionGenerator::~TransformedDirectionGenerator()
{
  if(adopt_gen_)delete gen_;
}

void TransformedDirectionGenerator::reset()
{
  gen_->reset();
}

bool TransformedDirectionGenerator::
next_as_theta_phi(double& theta, double& phi, double& weight)
{
  Eigen::Vector3d dir;
  if(gen_->next_as_vector(dir, weight)) {
    dir = rot_ * dir;
    phi = std::atan2(dir.y(),dir.x());
    theta = std::atan2(std::sqrt(SQR(dir.x())+SQR(dir.y())),dir.z());
    return true;
  } else {
    return false;
  }
}

bool TransformedDirectionGenerator::
next_as_vector(Eigen::Vector3d& dir, double& weight)
{
  if(gen_->next_as_vector(dir, weight)) {
    dir = rot_ * dir;
    return true;
  } else {
    return false;
  }
}
