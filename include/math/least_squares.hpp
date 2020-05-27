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

class I64LinearRegressionAccumulator
{
public:
  I64LinearRegressionAccumulator() { /* nothing to see here */ }
  I64LinearRegressionAccumulator(int64_t x0, int64_t y0): zero_set_(true),
    x0_(x0), y0_(y0) { /* nothing to see here */ }
  void accumulate(int64_t x, int64_t y);
  void integrate_into(I64LinearRegressionAccumulator& other);
  void fit_parameters(double& a, double& b) const;
  void fit_parameters_and_d2(double& a, double& b, double& D2) const;
private:
#ifndef SWIG
  bool zero_set_ = false;
  int64_t x0_ = 0;
  int64_t y0_ = 0;
  int64_t W_ = 0;
  int64_t X_ = 0;
  int64_t Y_ = 0;
  __int128_t XX_ = 0;
  __int128_t XY_ = 0;
  __int128_t YY_ = 0;
#endif
};

class KahanLinearRegressionAccumulator
{
public:
  KahanLinearRegressionAccumulator() { /* nothing to see here */ }
  KahanLinearRegressionAccumulator(double x0_, double y0_): zero_set(true),
    x0(x0_), y0(y0_) { /* nothing to see here */ }
  void accumulate(double x, double y);
  void integrate_into(KahanLinearRegressionAccumulator& other);
  void fit_parameters(double& a, double& b);
  void fit_parameters_and_d2(double& a, double& b, double& D2);

private:
  bool zero_set = false;
  double x0 = 0;
  double y0 = 0;
  calin::math::accumulator::SimpleAccumulator W;
  calin::math::accumulator::SimpleAccumulator X;
  calin::math::accumulator::SimpleAccumulator Y;
  calin::math::accumulator::KahanAccumulator XX;
  calin::math::accumulator::KahanAccumulator XY;
  calin::math::accumulator::KahanAccumulator YY;
};

Eigen::VectorXd polyfit(const Eigen::VectorXd& x, const Eigen::VectorXd& y, unsigned order);

} } } // namespace calin::math::least_squares
