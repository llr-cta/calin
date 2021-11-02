/*

   calin/math/m_estimate.hpp -- Stephen Fegan -- 2017-04-17

   Functions for impleminting M-Estimates

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

#pragma once

#include <string>

#include "function.hpp"

namespace calin { namespace math { namespace m_estimate {

class LikelihoodRhoFunction: public calin::math::function::SingleAxisFunction
{
public:
  virtual ~LikelihoodRhoFunction();
  calin::math::function::DomainAxis domain_axis() override;
  double value_1d(double x) override = 0;
  bool can_calculate_gradient() override; // default : true
  double value_and_gradient_1d(double x,  double& dfdx) override = 0;
  bool can_calculate_hessian() override; // default : true
  double value_gradient_and_hessian_1d(double x, double& dfdx, double& d2fdx2) override = 0;
  virtual double asymptotic_value() = 0;
};

class NullLikelihoodRhoFunction: public LikelihoodRhoFunction
{
public:
  NullLikelihoodRhoFunction();
  virtual ~NullLikelihoodRhoFunction();
  double value_1d(double x) override;
  double value_and_gradient_1d(double x,  double& dfdx) override;
  double value_gradient_and_hessian_1d(double x, double& dfdx,
                                       double& d2fdx2) override;
  double asymptotic_value() override;
};

class HyperbolicLikelihoodRhoFunction: public LikelihoodRhoFunction
{
public:
  HyperbolicLikelihoodRhoFunction(double asymptotic_value, double turnover_scale);
  virtual ~HyperbolicLikelihoodRhoFunction();
  double value_1d(double x) override;
  double value_and_gradient_1d(double x,  double& dfdx) override;
  double value_gradient_and_hessian_1d(double x, double& dfdx,
                                       double& d2fdx2) override;
  double asymptotic_value() override;
protected:
  double C_;
  double D2_;
};

class ModifiedHyperbolicLikelihoodRhoFunction: public LikelihoodRhoFunction
{
public:
  ModifiedHyperbolicLikelihoodRhoFunction(double asymptotic_value, double turnover_scale);
  virtual ~ModifiedHyperbolicLikelihoodRhoFunction();
  double value_1d(double x) override;
  double value_and_gradient_1d(double x,  double& dfdx) override;
  double value_gradient_and_hessian_1d(double x, double& dfdx,
                                       double& d2fdx2) override;
  double asymptotic_value() override;
protected:
  double D2_;
  double C_;
  double scale_;
  double offset_;
};

} } } // namespace calin::math::m_estimate
