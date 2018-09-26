/*

   calin/math/position_generator.hpp -- Stephen Fegan -- 2017-01-18

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

#pragma once

#include <Eigen/Dense>
#include <math/rng.hpp>

namespace calin { namespace math { namespace position_generator {

class PositionGenerator
{
public:
  virtual ~PositionGenerator();
  virtual void reset() = 0;
  virtual bool next(Eigen::Vector3d& pos, double& weight) = 0;
};

class SinglePositionGenerator: public PositionGenerator
{
public:
  SinglePositionGenerator(double x, double y, double z, double base_weight = 1.0);
  SinglePositionGenerator(const Eigen::Vector3d& pos, double base_weight = 1.0);
  virtual ~SinglePositionGenerator();
  void reset() override;
  bool next(Eigen::Vector3d& pos, double& weight) override;
  double weight() const { return weight_; }
private:
  const Eigen::Vector3d pos_;
  double weight_ = 1.0;
  bool position_generated_ = false;
};

class MCPlanePositionGenerator: public PositionGenerator
{
public:
  MCPlanePositionGenerator(double r_max, unsigned nray,
    calin::math::rng::RNG* rng,
    bool scale_weight_by_area = true, double base_weight = 1.0,
    bool adopt_rng = false);
  virtual ~MCPlanePositionGenerator();
  void reset() override;
  bool next(Eigen::Vector3d& pos, double& weight) override;
  double weight() const { return weight_; }
protected:
  unsigned iray_ = 0;
  unsigned nray_ = 0;
  double r2_max_ = 0.0;
  double weight_ = 1.0;
  calin::math::rng::RNG* rng_ = nullptr;
  bool adopt_rng_ = false;
};

class HexGridPlanePositionGenerator: public PositionGenerator
{
public:
  HexGridPlanePositionGenerator(double r_max, double dx = 1.0,
    bool scale_weight_by_area = true, double base_weight = 1.0);
  virtual ~HexGridPlanePositionGenerator();
  void reset() override;
  bool next(Eigen::Vector3d& pos, double& weight) override;
  double weight() const { return weight_; }
protected:
  unsigned hexid_ = 0;
  double r2_max_ = 0.0;
  double dx_ = 1.0;
  double weight_ = 1.0;
};

class TransformedPositionGenerator: public PositionGenerator
{
public:
  TransformedPositionGenerator(const Eigen::Vector3d& shift_post_rot,
    const Eigen::Matrix3d& rot, const Eigen::Vector3d& shift_pre_rot,
    PositionGenerator* gen, bool adopt_gen = false);
  TransformedPositionGenerator(const Eigen::Vector3d& shift_post_rot,
    const Eigen::Matrix3d& rot, PositionGenerator* gen, bool adopt_gen = false);
  virtual ~TransformedPositionGenerator();
  void reset() override;
  bool next(Eigen::Vector3d& pos, double& weight) override;
protected:
  PositionGenerator* gen_ = nullptr;
  bool adopt_gen_ = false;
  Eigen::Matrix3d rot_ = Eigen::Matrix3d::Identity();
  Eigen::Vector3d shift_post_rot_ = Eigen::Vector3d::Zero();
};

} } } // namespace calin::math::position_generator
