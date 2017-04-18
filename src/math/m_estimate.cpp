/*

   calin/math/m_estimates.cpp -- Stephen Fegan -- 2017-04-17

   Functions for impleminting M-Estimates

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

#include "math/m_estimate.hpp"
#include "math/special.hpp"

using calin::math::special::SQR;

using namespace calin::math::m_estimate;

LikelihoodRhoFunction::~LikelihoodRhoFunction()
{
  // nothing to see here
}

calin::math::function::DomainAxis LikelihoodRhoFunction::domain_axis()
{
  return { "-log(probability)", "", false, 0, false, 0 };
}

bool LikelihoodRhoFunction::can_calculate_gradient()
{
  return true;
}

bool LikelihoodRhoFunction::can_calculate_hessian()
{
  return true;
}

HyperbolicLikelihoodRhoFunction::
HyperbolicLikelihoodRhoFunction(double asymptotic_value, double turnover_scale):
  LikelihoodRhoFunction(), C_(asymptotic_value), D2_(SQR(turnover_scale))
{
  // nothing to see here
}

HyperbolicLikelihoodRhoFunction::~HyperbolicLikelihoodRhoFunction()
{
  // nothing to see here
}

double HyperbolicLikelihoodRhoFunction::value_1d(double x)
{
  double R2 = SQR(x-C_) + D2_;
  double R = std::sqrt(R2);
  return 0.5*(x+C_-R);
}

double HyperbolicLikelihoodRhoFunction::value_and_gradient_1d(double x, double& dfdx)
{
  double R2 = SQR(x-C_) + D2_;
  double R = std::sqrt(R2);
  dfdx = 0.5*(1 - (x-C_)/R);
  return 0.5*(x+C_-R);
}

double HyperbolicLikelihoodRhoFunction::
value_gradient_and_hessian_1d(double x, double& dfdx, double& d2fdx2)
{
  double R2 = SQR(x-C_) + D2_;
  double R = std::sqrt(R2);
  dfdx = 0.5*(1 - (x-C_)/R);
  d2fdx2 = -0.5*D2_/(R2*R);
  return 0.5*(x+C_-R);
}

ModifiedHyperbolicLikelihoodRhoFunction::
ModifiedHyperbolicLikelihoodRhoFunction(double asymptotic_value, double turnover_scale):
  LikelihoodRhoFunction(),
  D2_(SQR(turnover_scale)), C_(std::sqrt(SQR(asymptotic_value) - D2_)),
  scale_(1.0/(1.0 + C_/std::sqrt(SQR(C_)+D2_))), offset_(std::sqrt(SQR(C_)+D2_))
{
  // nothing to see here
}

ModifiedHyperbolicLikelihoodRhoFunction::~ModifiedHyperbolicLikelihoodRhoFunction()
{
  // nothing to see here
}

double ModifiedHyperbolicLikelihoodRhoFunction::value_1d(double x)
{

  double R2 = SQR(x-C_) + D2_;
  double R = std::sqrt(R2);
  return scale_*(x-R+offset_);
}

double ModifiedHyperbolicLikelihoodRhoFunction::value_and_gradient_1d(double x, double& dfdx)
{
  double R2 = SQR(x-C_) + D2_;
  double R = std::sqrt(R2);
  dfdx = scale_*(1 - (x-C_)/R);
  return scale_*(x-R+offset_);
}

double ModifiedHyperbolicLikelihoodRhoFunction::
value_gradient_and_hessian_1d(double x, double& dfdx, double& d2fdx2)
{
  double R2 = SQR(x-C_) + D2_;
  double R = std::sqrt(R2);
  dfdx = scale_*(1 - (x-C_)/R);
  d2fdx2 = -scale_*D2_/(R2*R);
  return scale_*(x-R+offset_);
}
