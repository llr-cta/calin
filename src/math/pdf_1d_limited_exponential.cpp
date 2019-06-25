/*

   calin/math/pdf_1d_limited_exponential.cpp -- Stephen Fegan -- 2015-04-02

   Base classes for one-dimensional PDF functions. PDF functions are
   based on ParameterizableSingleAxisFunction, and are assumed to
   integrate out to 1.0. They also (optionally) provide analytic
   moments.

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <iostream>
#include <iomanip>
#include <stdexcept>

#include <math/special.hpp>
#include <math/pdf_1d.hpp>

using namespace calin::math;
using namespace calin::math::pdf_1d;
using calin::math::special::SQR;
using function::assign_parameters;

// *****************************************************************************
//
// LimitedExponentialPDF
//
// *****************************************************************************

LimitedExponentialPDF::~LimitedExponentialPDF()
{
  // nothing to see here
}

unsigned LimitedExponentialPDF::num_parameters()
{
  return 1;
}

std::vector<function::ParameterAxis> LimitedExponentialPDF::parameters()
{
  return { { "scale", xunits_, limit_a_lo_>-inf, limit_a_lo_, limit_a_hi_<inf, limit_a_hi_ } };
}

Eigen::VectorXd LimitedExponentialPDF::parameter_values()
{
  Eigen::VectorXd p(1);
  p << a_;
  return p;
}

void LimitedExponentialPDF::set_parameter_values(ConstVecRef values)
{
  verify_set_parameter_values(values, "LimitedExponentialPDF");
  assign_parameters(values, a_);
  set_cache();
}

function::DomainAxis LimitedExponentialPDF::domain_axis()
{
  return { "x-value", xunits_, xlo_>-inf, xlo_, xhi_<inf, xhi_ };
}

bool LimitedExponentialPDF::can_calculate_gradient()
{
  return true;
}

bool LimitedExponentialPDF::can_calculate_hessian()
{
  return true;
}

bool LimitedExponentialPDF::can_calculate_parameter_gradient()
{
  return true;
}

bool LimitedExponentialPDF::can_calculate_parameter_hessian()
{
  return true;
}

double LimitedExponentialPDF::value_1d(double x)
{
  if(x<xlo_ or x>=xhi_)return 0;
  const double xs = x/a_;
  return norm_ * std::exp(-xs);
}

double LimitedExponentialPDF::value_and_gradient_1d(double x,  double& dfdx)
{
  if(x<xlo_ or x>=xhi_)
  {
    dfdx = 0;
    return 0;
  }
  const double xs = x/a_;
  const double val = norm_*std::exp(-xs);
  dfdx = -val/a_;
  return val;
}

double LimitedExponentialPDF::
value_gradient_and_hessian_1d(double x, double& dfdx, double& d2fdx2)
{
  if(x<xlo_ or x>=xhi_)
  {
    dfdx = 0;
    d2fdx2 = 0;
    return 0;
  }
  const double xs = x/a_;
  const double val = norm_*std::exp(-xs);
  dfdx = -val/a_;
  d2fdx2 = val/SQR(a_);
  return val;
}

double LimitedExponentialPDF::
value_and_parameter_gradient_1d(double x,  VecRef gradient)
{
  gradient.resize(1);
  if(x<xlo_ or x>=xhi_)
  {
    gradient(0) = 0;
    return 0;
  }
  const double xs = x/a_;
  const double val = std::exp(-xs);
  gradient(0) = val*(norm_gradient_ + norm_*xs/a_);
  return norm_*val;
}

double LimitedExponentialPDF::
value_parameter_gradient_and_hessian_1d(double x, VecRef gradient,
                                        MatRef hessian)
{
  gradient.resize(1);
  hessian.resize(1,1);
  if(x<xlo_ or x>=xhi_)
  {
    gradient(0) = 0;
    hessian(0,0) = 0;
    return 0;
  }
  const double xs = x/a_;
  const double val = std::exp(-xs);
  gradient(0) = val*(norm_gradient_ + norm_*xs/a_);
  hessian(0,0) = val*(norm_hessian_ + 2.0*norm_gradient_*xs/a_
                      + norm_*xs*(xs - 2.0)/SQR(a_));
  return norm_*val;
}

void LimitedExponentialPDF::set_cache()
{
  double ehi         = 0.0;
  double dehi_da     = 0.0;
  double d2ehi_da2   = 0.0;

  if(xhi_!=inf)
  {
    double xs = xhi_/a_;
    double exs = exp(-xs);

    ehi         = exs*a_;
    dehi_da     = exs*(1.0 + xs);
    d2ehi_da2   = exs*SQR(xs)/a_;
  }
  else if(a_<=0)
  {
    throw std::out_of_range("LimitedExponentialPDF: scale must be strictly "
                            "positive when xhi=inf");
  }

  double elo         = 0.0;
  double delo_da     = 0.0;
  double d2elo_da2   = 0.0;

  if(xlo_!=-inf)
  {
    double xs = xlo_/a_;
    double exs = exp(-xs);

    elo         = exs*a_;
    delo_da     = exs*(1.0 + xs);
    d2elo_da2   = exs*SQR(xs)/a_;
  }
  else if(a_>=0)
  {
    throw std::out_of_range("LimitedExponentialPDF: scale must be strictly "
                            "negative when xlo=-inf");
  }

  norm_              = -1.0/(ehi - elo);
  const double norm2 = SQR(norm_);
  norm_gradient_     = norm2*(dehi_da - delo_da);
  norm_hessian_      = norm2*(d2ehi_da2 - d2elo_da2)
                       + 2.0*SQR(norm_gradient_)/norm_;

  if(dx_ != 0)
  {
    double dxs = 0.5*dx_/a_;
    double bf = -a_*(std::exp(-dxs) - exp(dxs))/dx_; // binning factor
    double dbfda = -(std::exp(-dxs)*(1.0+dxs) - exp(dxs)*(1.0-dxs))/dx_;
    double d2bfda2 = -0.5*(std::exp(-dxs) - exp(dxs))*dxs/SQR(a_);

    norm_hessian_ = norm_hessian_*bf + 2.0*norm_gradient_*dbfda + norm_*d2bfda2;
    norm_gradient_ = norm_gradient_*bf + norm_*dbfda;
    norm_ *= bf;
  }
}
