/*

   calin/math/polyfit.hpp -- Stephen Fegan -- 2019-03-07

   Implemantion of the clasic polynomial fit by least-squares

   Copyright 2019, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#pragma once

#include <Eigen/Dense>
#include <math/accumulator.hpp>
#include <util/log.hpp>

namespace calin { namespace math { namespace least_squares {

// The purpose of this class is to do linear regression on highly correlated
// integer values (such as clock tick values). In this case high (full)
// precision must be maintained to allow for cancellation of almost all the most
// significant bits. Accumulation is done in 128 bits and the calculation of the
// least-square distance in 256 bits. For the clock case this is not overkill !

class I64LinearRegressionAccumulator
{
public:
  I64LinearRegressionAccumulator() { /* nothing to see here */ }
  I64LinearRegressionAccumulator(int64_t x0, int64_t y0): zero_set_(true),
    x0_(x0), y0_(y0) { /* nothing to see here */ }

  void accumulate(int64_t x, int64_t y);
  void integrate_into(I64LinearRegressionAccumulator& other);

  void rebalance();
  void shift_origin(int64_t x0, int64_t y0);

  void fit_parameters(double& a, double& b) const;
  void fit_parameters_and_d2(double& a, double& b, double& D2) const;

  int64_t x0() const { return x0_; }
  int64_t y0() const { return y0_; }
  int64_t num_entries() const { return W_; }

  double mean_x() const;
  double mean_y() const;

  void moments(double& entries, double& mean_x, double& mean_y,
    double& sigma_xx, double& sigma_yy, double& sigma_xy) const;

  void dump_to_log() const;
private:
#ifndef SWIG
  bool zero_set_ = false;
  int64_t x0_ = 0;
  int64_t y0_ = 0;
  __int128_t W_ = 0;
  __int128_t X_ = 0;
  __int128_t Y_ = 0;
  __int128_t XX_ = 0;
  __int128_t XY_ = 0;
  __int128_t YY_ = 0;
#endif
};

class KahanLinearRegressionAccumulator
{
public:
  KahanLinearRegressionAccumulator() { /* nothing to see here */ }
  KahanLinearRegressionAccumulator(double x0, double y0): zero_set_(true),
    x0_(x0), y0_(y0) { /* nothing to see here */ }
  void accumulate(double x, double y);
  void integrate_into(KahanLinearRegressionAccumulator& other);
  void fit_parameters(double& a, double& b);
  void fit_parameters_and_d2(double& a, double& b, double& D2);

private:
  bool zero_set_ = false;
  double x0_ = 0;
  double y0_ = 0;
  calin::math::accumulator::SimpleAccumulator W_;
  calin::math::accumulator::SimpleAccumulator X_;
  calin::math::accumulator::SimpleAccumulator Y_;
  calin::math::accumulator::KahanAccumulator XX_;
  calin::math::accumulator::KahanAccumulator XY_;
  calin::math::accumulator::KahanAccumulator YY_;
};

Eigen::VectorXd polyfit(const Eigen::VectorXd& x, const Eigen::VectorXd& y, unsigned order);

} } } // namespace calin::math::least_squares
