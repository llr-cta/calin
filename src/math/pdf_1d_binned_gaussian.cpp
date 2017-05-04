/*

   calin/math/pdf_1d_binned_gaussian.cpp -- Stephen Fegan -- 2016-07-04

   Binned Gaussian PDF

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
using calin::math::special::CUBE;
using function::assign_parameters;

namespace {

constexpr double C_SQRT1_2PI = 0.5*M_2_SQRTPI*M_SQRT1_2;

} // anonymous namespace

// *****************************************************************************
//
// BinnedGaussianPDF
//
// *****************************************************************************

BinnedGaussianPDF::~BinnedGaussianPDF()
{
  // nothing to see here
}

unsigned BinnedGaussianPDF::num_parameters()
{
  return 2;
}

std::vector<function::ParameterAxis> BinnedGaussianPDF::parameters()
{
  constexpr double inf = std::numeric_limits<double>::infinity();
  constexpr double tiny_val = std::numeric_limits<double>::min();
  return { { "mean", xunits_, false, -inf, false, inf },
      { "rms", xunits_, true, tiny_val, false, inf } };
}

Eigen::VectorXd BinnedGaussianPDF::parameter_values()
{
  Eigen::VectorXd pval(2);
  pval << x0_, s_;
  return pval;
}

void BinnedGaussianPDF::set_parameter_values(ConstVecRef values)
{
  assign_parameters(values, x0_, s_);
}

bool BinnedGaussianPDF::can_calculate_gradient()
{
  return true;
}

bool BinnedGaussianPDF::can_calculate_hessian()
{
  return true;
}

bool BinnedGaussianPDF::can_calculate_parameter_gradient()
{
  return true;
}

bool BinnedGaussianPDF::can_calculate_parameter_hessian()
{
  return true;
}

function::DomainAxis BinnedGaussianPDF::domain_axis()
{
  constexpr double inf = std::numeric_limits<double>::infinity();
  return { "x-value", xunits_, false, -inf, false, inf };
}

double BinnedGaussianPDF::value_1d(double x)
{
  const double xc = x-x0_;
  double val;
  if(dx_2_ == 0)
  {
    const double xs = xc/s_;
    const double xs2 = SQR(xs);
    val = C_SQRT1_2PI/s_*exp(-0.5*xs2);
  }
  else
  {
    const double xsr = (xc+dx_2_)/s_*M_SQRT1_2;
    const double xsl = (xc-dx_2_)/s_*M_SQRT1_2;
    if(xsl > 5)val = erfc(xsl) - erfc(xsr);
    else if(xsr<-5)val = erfc(-xsr) - erfc(-xsr);
    else val = erf(xsr) - erf(xsl);
    val *= 0.25/dx_2_;
  }
  return val;
}

double BinnedGaussianPDF::value_and_gradient_1d(double x,  double& dfdx)
{
  const double xc = x-x0_;
  double val;
  if(dx_2_ == 0)
  {
    const double xs = xc/s_;
    const double xs2 = SQR(xs);
    val = C_SQRT1_2PI/s_*exp(-0.5*xs2);
    dfdx = -val*xs/s_;
  }
  else
  {
    const double xsr = (xc+dx_2_)/s_*M_SQRT1_2;
    const double xsl = (xc-dx_2_)/s_*M_SQRT1_2;
    if(xsl > 5)val = erfc(xsl) - erfc(xsr);
    else if(xsr<-5)val = erfc(-xsr) - erfc(-xsr);
    else val = erf(xsr) - erf(xsl);
    val *= 0.25/dx_2_;
    dfdx = 0.5*C_SQRT1_2PI/dx_2_/s_*(exp(-SQR(xsr))-exp(-SQR(xsl)));
  }
  return val;
}

double BinnedGaussianPDF::
value_gradient_and_hessian_1d(double x, double& dfdx, double& d2fdx2)
{
  const double xc = x-x0_;
  double val;
  if(dx_2_ == 0)
  {
    const double xs = xc/s_;
    const double xs2 = SQR(xs);
    val = C_SQRT1_2PI/s_*exp(-0.5*xs2);
    dfdx = -val*xs/s_;
    d2fdx2 = -(dfdx*xc + val)/SQR(s_);
  }
  else
  {
    const double xsr = (xc+dx_2_)/s_*M_SQRT1_2;
    const double xsl = (xc-dx_2_)/s_*M_SQRT1_2;
    if(xsl > 5)val = erfc(xsl) - erfc(xsr);
    else if(xsr<-5)val = erfc(-xsr) - erfc(-xsr);
    else val = erf(xsr) - erf(xsl);
    val *= 0.25/dx_2_;
    const double exp_xsr2 = exp(-SQR(xsr));
    const double exp_xsl2 = exp(-SQR(xsl));
    dfdx = 0.5*C_SQRT1_2PI/dx_2_/s_*(exp_xsr2-exp_xsl2);
    d2fdx2 = -0.5*M_SQRT2*C_SQRT1_2PI/SQR(s_)/dx_2_*(xsr*exp_xsr2-xsl*exp_xsl2);
  }
  return val;
}

double BinnedGaussianPDF::
value_and_parameter_gradient_1d(double x,  VecRef gradient)
{
  const double xc = x-x0_;
  double val;
  if(dx_2_ == 0)
  {
    const double xs = xc/s_;
    const double xs2 = SQR(xs);
    val = C_SQRT1_2PI/s_*exp(-0.5*xs2);
    gradient.resize(2);
    gradient[0] = val*xs/s_; // df/dx0
    gradient[1] = val*(xs2 - 1.0)/s_;
  }
  else
  {
    const double xsr = (xc+dx_2_)/s_*M_SQRT1_2;
    const double xsl = (xc-dx_2_)/s_*M_SQRT1_2;
    if(xsl > 5)val = erfc(xsl) - erfc(xsr);
    else if(xsr<-5)val = erfc(-xsr) - erfc(-xsr);
    else val = erf(xsr) - erf(xsl);
    const double exp_xsr2 = exp(-SQR(xsr));
    const double exp_xsl2 = exp(-SQR(xsl));
    val *= 0.25/dx_2_;
    gradient[0] = -0.5*C_SQRT1_2PI/dx_2_/s_*(exp_xsr2-exp_xsl2);
    gradient[1] = -0.25*M_2_SQRTPI/dx_2_/s_*(xsr*exp_xsr2-xsl*exp_xsl2);
  }
  return val;
}

double BinnedGaussianPDF::
value_parameter_gradient_and_hessian_1d(double x, VecRef gradient,
                                        MatRef hessian)
{
  const double xc = x-x0_;
  double val;
  if(dx_2_ == 0)
  {
    const double xs = xc/s_;
    const double xs2 = SQR(xs);
    val = C_SQRT1_2PI/s_*exp(-0.5*xs2);
    gradient.resize(2);
    gradient[0] = val*xs/s_; // df/dx0
    gradient[1] = val*(xs2 - 1.0)/s_;
    hessian(0,0) = gradient[1]/s_; // val*(xs2 - 1.0)/SQR(s_);
    hessian(0,1) = val*(xs2 - 3.0)*xs/SQR(s_);
    hessian(1,0) = hessian(0,1);
    hessian(1,1) = val*(SQR(xs2) - 5.0*xs2 + 2.0)/SQR(s_);
  }
  else
  {
    const double xsr = (xc+dx_2_)/s_*M_SQRT1_2;
    const double xsl = (xc-dx_2_)/s_*M_SQRT1_2;
    if(xsl > 5)val = erfc(xsl) - erfc(xsr);
    else if(xsr<-5)val = erfc(-xsr) - erfc(-xsr);
    else val = erf(xsr) - erf(xsl);
    const double exp_xsr2 = exp(-SQR(xsr));
    const double exp_xsl2 = exp(-SQR(xsl));
    val *= 0.25/dx_2_;
    gradient[0] = -0.5*C_SQRT1_2PI/dx_2_/s_*(exp_xsr2-exp_xsl2);
    gradient[1] = -0.25*M_2_SQRTPI/dx_2_/s_*(xsr*exp_xsr2-xsl*exp_xsl2);
    hessian(0,0) = gradient[1]/s_;
    hessian(0,1) = -(gradient[0]/s_ +
      C_SQRT1_2PI/dx_2_/SQR(s_)*(SQR(xsr)*exp_xsr2 - SQR(xsl)*exp_xsl2));
    hessian(1,0) = hessian(0,1);
    hessian(1,1) = -2.0*gradient[1]/s_ -
      0.5*M_2_SQRTPI/dx_2_/SQR(s_)*(CUBE(xsr)*exp_xsr2 - CUBE(xsl)*exp_xsl2);
  }
  return val;
}
