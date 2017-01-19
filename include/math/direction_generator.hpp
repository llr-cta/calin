/*

   calin/math/direction_generator.hpp -- Stephen Fegan -- 2017-01-19

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

#pragma once

#include <Eigen/Dense>
#include <math/rng.hpp>

namespace calin { namespace math { namespace direction_generator {

class DirectionGenerator
{
public:
  virtual ~DirectionGenerator();
  virtual void reset() = 0;
  virtual bool next_as_theta_phi(double& theta, double& phi, double& weight) = 0;
  virtual bool next_as_vector(Eigen::Vector3d& dir, double& weight);
  virtual bool next_as_matrix(Eigen::Matrix3d& trans_mat, double& weight);
};

class MCSphereDirectionGenerator: public DirectionGenerator
{
public:
  MCSphereDirectionGenerator(double theta_max, unsigned nray,
    calin::math::rng::RNG* rng,
    bool scale_weight_by_area = true, double base_weight = 1.0,
    bool adopt_rng = false);
  virtual ~MCSphereDirectionGenerator();
  void reset() override;
  bool next_as_theta_phi(double& theta, double& phi, double& weight) override;
protected:
  unsigned iray_ = 0;
  unsigned nray_ = 0;
  double cos_theta_max_ = 0.0;
  double weight_ = 1.0;
  calin::math::rng::RNG* rng_ = nullptr;
  bool adopt_rng_ = false;
};

class HEALPixDirectionGenerator: public DirectionGenerator
{
public:
  HEALPixDirectionGenerator(double theta_max, unsigned nside,
    bool scale_weight_by_area = true, double base_weight = 1.0);
  virtual ~HEALPixDirectionGenerator();
  void reset() override;
  bool next_as_theta_phi(double& theta, double& phi, double& weight) override;
  bool next_as_vector(Eigen::Vector3d& dir, double& weight) override;
  bool next_as_matrix(Eigen::Matrix3d& trans_mat, double& weight) override;
protected:
  unsigned pixid_ = 0;
  double cos_theta_max_ = 1.0;
  unsigned nside_ = 1;
  double weight_ = 1.0;
};

} } } // namespace calin::math::position_generator
