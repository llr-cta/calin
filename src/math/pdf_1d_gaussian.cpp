/*

   calin/math/pdf_1d_gaussian.cpp -- Stephen Fegan -- 2015-04-02

   Base classes for one-dimensional PDF functions. PDF functions are
   based on ParameterizableSingleAxisFunction, and are assumed to
   integrate out to 1.0. They also (optionally) provide analytic
   moments.

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <iostream>
#include <iomanip>
#include <stdexcept>


#include <math/special.hpp>
#include <math/pdf_1d.hpp>

using namespace calin::math;
using namespace calin::math::pdf_1d;
using calin::math::special::SQR;
using function::assign_parameters;

namespace {

constexpr double c_gauss_norm = 0.5*M_2_SQRTPI*M_SQRT1_2;

} // anonymous namespace

// *****************************************************************************
//
// GaussianPDF
//
// *****************************************************************************

GaussianPDF::~GaussianPDF()
{
  // nothing to see here
}

unsigned GaussianPDF::num_parameters()
{
  return 2;
}

std::vector<function::ParameterAxis> GaussianPDF::parameters()
{
  constexpr double inf = std::numeric_limits<double>::infinity();
  constexpr double tiny_val = std::numeric_limits<double>::min();
  return { { "mean", xunits_, false, -inf, false, inf },
      { "rms", xunits_, true, tiny_val, false, inf } };
}

Eigen::VectorXd GaussianPDF::parameter_values()
{
  Eigen::VectorXd pval(2);
  pval << x0_, s_;
  return pval;
}

void GaussianPDF::set_parameter_values(ConstVecRef values)
{
  assign_parameters(values, x0_, s_);
}

bool GaussianPDF::can_calculate_gradient()
{
  return true;
}

bool GaussianPDF::can_calculate_hessian()
{
  return true;
}

bool GaussianPDF::can_calculate_parameter_gradient()
{
  return true;
}

bool GaussianPDF::can_calculate_parameter_hessian()
{
  return true;
}

function::DomainAxis GaussianPDF::domain_axis()
{
  constexpr double inf = std::numeric_limits<double>::infinity();
  return { "x-value", xunits_, false, -inf, false, inf };
}

double GaussianPDF::value_1d(double x)
{
  const double xc = x-x0_;
  const double xs = xc/s_;
  const double xs2 = SQR(xs);
  double val = c_gauss_norm/s_*exp(-0.5*xs2);
  //std::cout << std::setprecision(17) << x << ' ' << val << '\n';
  return val;
}

double GaussianPDF::value_and_gradient_1d(double x,  double& dfdx)
{
  const double xc = x-x0_;
  const double xs = xc/s_;
  const double xs2 = SQR(xs);
  double val = c_gauss_norm/s_*exp(-0.5*xs2);
  dfdx = -val*xs/s_;
  return val;
}

double GaussianPDF::
value_gradient_and_hessian_1d(double x, double& dfdx, double& d2fdx2)
{
  const double xc = x-x0_;
  const double xs = xc/s_;
  const double xs2 = SQR(xs);
  double val = c_gauss_norm/s_*exp(-0.5*xs2);
  dfdx = -val*xs/s_;
  d2fdx2 = -(dfdx*xc + val)/SQR(s_);
  return val;
}

double GaussianPDF::
value_and_parameter_gradient_1d(double x,  VecRef gradient)
{
  const double xc = x-x0_;
  const double xs = xc/s_;
  const double xs2 = SQR(xs);
  double val = c_gauss_norm/s_*exp(-0.5*xs2);
  gradient.resize(2);
  gradient[0] = val*xs/s_; // df/dx0
  gradient[1] = val*(xs2 - 1.0)/s_;
  return val;
}

double GaussianPDF::
value_parameter_gradient_and_hessian_1d(double x, VecRef gradient,
                                        MatRef hessian)
{
  const double xc = x-x0_;
  const double xs = xc/s_;
  const double xs2 = SQR(xs);
  double val = c_gauss_norm/s_*exp(-0.5*xs2);
  gradient.resize(2);
  gradient[0] = val*xs/s_; // df/dx0
  gradient[1] = val*(xs2 - 1.0)/s_;
  hessian(0,0) = gradient[1]/s_; // val*(xs2 - 1.0)/SQR(s_);
  hessian(0,1) = val*(xs2 - 3.0)*xs/SQR(s_);
  hessian(1,0) = hessian(0,1);
  hessian(1,1) = val*(SQR(xs2) - 5.0*xs2 + 2.0)/SQR(s_);
  return val;
}

// *****************************************************************************
//
// LimitedGaussianPDF
//
// *****************************************************************************

LimitedGaussianPDF::~LimitedGaussianPDF()
{
  // nothing to see here
}

function::DomainAxis LimitedGaussianPDF::domain_axis()
{
  function::DomainAxis a = BinnedGaussianPDF::domain_axis();
  a.has_lo_bound = xlo_>-inf;
  a.lo_bound = xlo_;
  a.has_hi_bound = xhi_<inf;
  a.hi_bound = xhi_;
  return a;
}

void LimitedGaussianPDF::set_parameter_values(ConstVecRef values)
{
  BinnedGaussianPDF::set_parameter_values(values);
  set_cache();
}

double LimitedGaussianPDF::value_1d(double x)
{
  if(x<xlo_ or x>=xhi_)return 0;
  return norm_*BinnedGaussianPDF::value_1d(x);
}

double LimitedGaussianPDF::value_and_gradient_1d(double x,  double& dfdx)
{
  if(x<xlo_ or x>=xhi_)
  {
    dfdx = 0;
    return 0;
  }
  double val = norm_*BinnedGaussianPDF::value_and_gradient_1d(x, dfdx);
  dfdx *= norm_;
  return val;
}

double LimitedGaussianPDF::value_gradient_and_hessian_1d(double x, double& dfdx,
                                              double& d2fdx2)
{
  if(x<xlo_ or x>=xhi_)
  {
    d2fdx2 = 0;
    dfdx = 0;
    return 0;
  }
  double val = BinnedGaussianPDF::value_gradient_and_hessian_1d(x, dfdx, d2fdx2);
  dfdx *= norm_;
  d2fdx2 *= norm_;
  return norm_*val;
}

double LimitedGaussianPDF::
value_and_parameter_gradient_1d(double x,  VecRef gradient)
{
  gradient.resize(2);
  if(x<xlo_ or x>=xhi_)
  {
    gradient[0] = gradient[1] = 0;
    return 0;
  }
  double val = BinnedGaussianPDF::value_and_parameter_gradient_1d(x, gradient);
  gradient = norm_*gradient + val*norm_gradient_;
  return norm_*val;
}

double LimitedGaussianPDF::
value_parameter_gradient_and_hessian_1d(double x, VecRef gradient,
                                        MatRef hessian)
{
  gradient.resize(2);
  hessian.resize(2,2);
  if(x<xlo_ or x>=xhi_)
  {
    gradient[0] = gradient[1] = 0;
    hessian(0,0) = hessian(0,1) = hessian(1,0) = hessian(1,1) = 0;
    return 0;
  }
  double val =
     BinnedGaussianPDF::value_parameter_gradient_and_hessian_1d(x, gradient, hessian);
  hessian = norm_*hessian + val*norm_hessian_;
  hessian(0,0) += 2.0*norm_gradient_[0]*gradient[0];
  hessian(0,1) += norm_gradient_[0]*gradient[1] + norm_gradient_[1]*gradient[0];
  hessian(1,0) = hessian(0,1);
  hessian(1,1) += 2.0*norm_gradient_[1]*gradient[1];
  gradient = norm_*gradient + val*norm_gradient_;
  return norm_*val;
}

void LimitedGaussianPDF::set_cache()
{
  double ehi         = 1.0;
  double dehi_dx0    = 0.0;
  double dehi_ds     = 0.0;
  double d2ehi_dx02  = 0.0;
  double d2ehi_ds2   = 0.0;
  double d2ehi_dx0ds = 0.0;

  if(xhi_!=inf)
  {
    double xc = xhi_-x0_;
    double xs = xc/s_;
    double xs2 = SQR(xs);

    ehi         = 0.5*(1.0+std::erf(M_SQRT1_2*xs));

    dehi_dx0    = -c_gauss_norm/s_*exp(-0.5*xs2);
    dehi_ds     = xs*dehi_dx0;

    d2ehi_dx02  = dehi_ds/s_;
    d2ehi_dx0ds = dehi_dx0*(xs2-1)/s_;
    d2ehi_ds2   = xs*d2ehi_dx0ds - d2ehi_dx02;
  }

  double elo         = 0.0;
  double delo_dx0    = 0.0;
  double delo_ds     = 0.0;
  double d2elo_dx02  = 0.0;
  double d2elo_ds2   = 0.0;
  double d2elo_dx0ds = 0.0;

  if(xlo_!=-inf)
  {
    double xc = xlo_-x0_;
    double xs = xc/s_;
    double xs2 = SQR(xs);

    elo         = 0.5*(1.0+std::erf(M_SQRT1_2*xs));

    delo_dx0    = -c_gauss_norm/s_*exp(-0.5*xs2);
    delo_ds     = xs*delo_dx0;

    d2elo_dx02  = delo_ds/s_;
    d2elo_dx0ds = delo_dx0*(xs2-1)/s_;
    d2elo_ds2   = xs*d2elo_dx0ds - d2elo_dx02;
  }

  norm_              = 1.0/(ehi - elo);
  const double norm2 = SQR(norm_);
  norm_gradient_(0)  = -norm2*(dehi_dx0 - delo_dx0);
  norm_gradient_(1)  = -norm2*(dehi_ds - delo_ds);
  norm_hessian_(0,0) = -norm2*(d2ehi_dx02 - d2elo_dx02)
                       + 2.0*SQR(norm_gradient_(0))/norm_;
  norm_hessian_(0,1) = -norm2*(d2ehi_dx0ds - d2elo_dx0ds)
                       + 2.0*norm_gradient_(0)*norm_gradient_(1)/norm_;
  norm_hessian_(1,0) = norm_hessian_(0,1);
  norm_hessian_(1,1) = -norm2*(d2ehi_ds2 - d2elo_ds2)
                       + 2.0*SQR(norm_gradient_(1))/norm_;
}
