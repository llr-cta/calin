/* 

   calin/math/pdf_1d.cpp -- Stephen Fegan -- 2015-04-02

   Base classes for one-dimensional PDF functions. PDF functions are
   based on ParameterizableSingleAxisFunction, and are assumed to
   integrate out to 1.0. They also (optionally) provide analytic
   moments.

   Copyright 2015, Stephen Fegan <sfegan@gmail.com>

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

#include "math/pdf_1d.hpp"

using namespace calin::math;
using namespace calin::math::pdf_1d;

namespace {

inline static double SQR(double x) { return x*x; }
constexpr double c_gauss_norm = 0.5*M_2_SQRTPI*M_SQRT1_2;

} // anonymous namespace

using function::assign_parameters;

Parameterizable1DPDF::~Parameterizable1DPDF()
{
  // nothing to see here
}

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

double GaussianPDF::error_up()
{
  return error_up_;
}

bool GaussianPDF::can_calculate_mean_and_variance()
{
  return true;
}

void GaussianPDF::mean_and_variance(double& mean, double& var)
{
  mean = x0_;
  var = SQR(s_);
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
  function::DomainAxis a = GaussianPDF::domain_axis();
  a.has_lo_bound = xlo_>-inf;
  a.lo_bound = xlo_;
  a.has_hi_bound = xhi_<inf;
  a.hi_bound = xhi_;
  return a;
}

void LimitedGaussianPDF::set_parameter_values(ConstVecRef values)
{
  GaussianPDF::set_parameter_values(values);
  set_cache();
}

double LimitedGaussianPDF::value_1d(double x)
{
  if(x<xlo_ or x>=xhi_)return 0;
  return norm_*GaussianPDF::value_1d(x);
}
    
double LimitedGaussianPDF::value_and_gradient_1d(double x,  double& dfdx)
{
  if(x<xlo_ or x>=xhi_)
  {
    dfdx = 0;
    return 0;
  }
  double val = norm_*GaussianPDF::value_and_gradient_1d(x, dfdx);
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
  double val = GaussianPDF::value_gradient_and_hessian_1d(x, dfdx, d2fdx2);
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
  double val = GaussianPDF::value_and_parameter_gradient_1d(x, gradient);
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
     GaussianPDF::value_parameter_gradient_and_hessian_1d(x, gradient, hessian);
  hessian = norm_*hessian + val*norm_hessian_;
  hessian(0,0) += 2.0*norm_gradient_[0]*gradient[0];
  hessian(0,1) += norm_gradient_[0]*gradient[1] + norm_gradient_[1]*gradient[0];
  hessian(1,0) = hessian(0,1);
  hessian(1,1) += 2.0*norm_gradient_[1]*gradient[1];
  gradient = norm_*gradient + val*norm_gradient_;
  return norm_*val;
}

bool LimitedGaussianPDF::can_calculate_mean_and_variance()
{
  return false;
}

void LimitedGaussianPDF::mean_and_variance(double& mean, double& var)
{
  assert(0);
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

double LimitedExponentialPDF::error_up()
{
  return error_up_;
}

bool LimitedExponentialPDF::can_calculate_mean_and_variance()
{
  return false;
}

void LimitedExponentialPDF::mean_and_variance(double& mean, double& var)
{
  assert(0);
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

// ============================================================================
//
// TwoComponentPDF
//
// ============================================================================

TwoComponentPDF::
TwoComponentPDF(Parameterizable1DPDF* pdf1, const std::string& cpt1_name,
                Parameterizable1DPDF* pdf2, const std::string& cpt2_name,
                bool adopt_pdf1, bool adopt_pdf2, double error_up):
    Parameterizable1DPDF(),
    pdf1_(pdf1), adopt_pdf1_(adopt_pdf1), cpt1_name_(cpt1_name), 
    pdf2_(pdf2), adopt_pdf2_(adopt_pdf2), cpt2_name_(cpt2_name),
    error_up_(error_up)
{
  // nothing to see here
}

TwoComponentPDF::~TwoComponentPDF()
{
  // nothing to see here
}

unsigned TwoComponentPDF::num_parameters()
{
  return 1 + pdf1_->num_parameters() + pdf2_->num_parameters();
}
 
std::vector<function::ParameterAxis> TwoComponentPDF::parameters()
{
  std::vector<function::ParameterAxis> pvec {
    { cpt1_name_+std::string("_probability"), "1", true, 0, true, 1 } };
  std::vector<function::ParameterAxis> pdf1_pvec = pdf1_->parameters();
  for(auto& p : pdf1_pvec)p.name = cpt1_name_ + std::string(".") + p.name;
  pvec.insert(pvec.end(), pdf1_pvec.begin(), pdf1_pvec.end());
  std::vector<function::ParameterAxis> pdf2_pvec = pdf2_->parameters();
  for(auto& p : pdf2_pvec)p.name = cpt2_name_ + std::string(".") + p.name;
  pvec.insert(pvec.end(), pdf2_pvec.begin(), pdf2_pvec.end());
  return pvec;
}

Eigen::VectorXd TwoComponentPDF::parameter_values()
{
  Eigen::VectorXd param(num_parameters());
  param[0] = prob_cpt1_;
  unsigned num_cpt1_params = pdf1_->num_parameters();
  param.segment(1,num_cpt1_params)
      = pdf1_->parameter_values();
  param.segment(num_cpt1_params+1,pdf2_->num_parameters())
      = pdf2_->parameter_values();
  return param;
}

void TwoComponentPDF::set_parameter_values(ConstVecRef values)
{
  assign_parameters(values, prob_cpt1_);
  unsigned num_cpt1_params = pdf1_->num_parameters();
  pdf1_->set_parameter_values(values.segment(1,num_cpt1_params));
  pdf2_->set_parameter_values(values.segment(num_cpt1_params+1,
                                             pdf2_->num_parameters()));
}

function::DomainAxis TwoComponentPDF::domain_axis()
{
  return pdf1_->domain_axis();
}

bool TwoComponentPDF::can_calculate_gradient()
{
  return pdf1_->can_calculate_gradient() and pdf2_->can_calculate_gradient();
}

bool TwoComponentPDF::can_calculate_hessian()
{
  return pdf1_->can_calculate_hessian() and pdf2_->can_calculate_hessian();
}

bool TwoComponentPDF::can_calculate_parameter_gradient()
{
  return pdf1_->can_calculate_parameter_gradient() and
      pdf2_->can_calculate_parameter_gradient();
}

bool TwoComponentPDF::can_calculate_parameter_hessian()
{
  return pdf1_->can_calculate_parameter_hessian() and
      pdf2_->can_calculate_parameter_hessian();
}

double TwoComponentPDF::value_1d(double x)
{
  //std::cout << pdf1_->value_1d(x) << ' ' << pdf2_->value_1d(x) << '\n';
  return prob_cpt1_*pdf1_->value_1d(x) + (1.0-prob_cpt1_)*pdf2_->value_1d(x);
}

double TwoComponentPDF::value_and_gradient_1d(double x,  double& dfdx)
{
  double dfdx1;
  double dfdx2;
  double val = prob_cpt1_*pdf1_->value_and_gradient_1d(x,dfdx1) +
               (1.0-prob_cpt1_)*pdf2_->value_and_gradient_1d(x,dfdx2);
  dfdx = prob_cpt1_*dfdx1 + (1.0-prob_cpt1_)*dfdx2;
  return val;
}

double TwoComponentPDF::value_gradient_and_hessian_1d(double x, double& dfdx,
                                           double& d2fdx2)
{
  double dfdx1;
  double dfdx2;
  double d2fdx21;
  double d2fdx22;
  double val = prob_cpt1_*pdf1_->value_gradient_and_hessian_1d(x,dfdx1,d2fdx21)+
        (1.0-prob_cpt1_)*pdf2_->value_gradient_and_hessian_1d(x,dfdx2,d2fdx22);
  dfdx = prob_cpt1_*dfdx1 + (1.0-prob_cpt1_)*dfdx2;
  d2fdx2 = prob_cpt1_*d2fdx21 + (1.0-prob_cpt1_)*d2fdx22;
  return val;
}

double
TwoComponentPDF::value_and_parameter_gradient_1d(double x,  VecRef gradient)
{
  const double omp = 1.0-prob_cpt1_;
  const unsigned npar1 = pdf1_->num_parameters();
  const unsigned npar2 = pdf2_->num_parameters();
  gradient.resize(1+npar1+npar2);

#ifndef CALIN_USE_EIGEN_REF
  Eigen::VectorXd grad1(npar1);
  Eigen::VectorXd grad2(npar2);
#else
  auto grad1 = gradient.segment(1,npar1);
  auto grad2 = gradient.segment(1+npar1,npar2);
#endif
  
  double val1 = pdf1_->value_and_parameter_gradient_1d(x, grad1);
  double val2 = pdf2_->value_and_parameter_gradient_1d(x, grad2);

  grad1 *= prob_cpt1_;
  grad2 *= omp;
  gradient[0] = val1 - val2;
#ifndef CALIN_USE_EIGEN_REF
  gradient.segment(1,npar1) = grad1;
  gradient.segment(1+npar1,npar2) = grad2;
#endif
  return prob_cpt1_*val1 + omp*val2;
}

double TwoComponentPDF::
value_parameter_gradient_and_hessian_1d(double x, VecRef gradient,
                                        MatRef hessian)
{
  const double omp = 1.0-prob_cpt1_;
  const unsigned npar1 = pdf1_->num_parameters();
  const unsigned npar2 = pdf2_->num_parameters();
  gradient.resize(1+npar1+npar2);
  hessian.resize(1+npar1+npar2,1+npar1+npar2);
  hessian.setZero();

#ifndef CALIN_USE_EIGEN_REF
  Eigen::VectorXd grad1(npar1);
  Eigen::VectorXd grad2(npar2);
  Eigen::MatrixXd hess1(npar1, npar1);
  Eigen::MatrixXd hess2(npar2, npar2);
#else
  auto grad1 = gradient.segment(1,npar1);
  auto grad2 = gradient.segment(1+npar1,npar2);
  auto hess1 = hessian.block(1,1,npar1,npar1);
  auto hess2 = hessian.block(1+npar1,1+npar1,npar2,npar2);
#endif

  double val1 = pdf1_->value_parameter_gradient_and_hessian_1d(x, grad1, hess1);
  double val2 = pdf2_->value_parameter_gradient_and_hessian_1d(x, grad2, hess2);

  hess1 *= prob_cpt1_;
  hess2 *= omp;
#ifndef CALIN_USE_EIGEN_REF
  hessian.block(1,1,npar1,npar1) = hess1;
  hessian.block(1+npar1,1+npar1,npar2,npar2) = hess2;
#endif
  hessian.block(1,0,npar1,1) = grad1;
  hessian.block(1+npar1,0,npar2,1) = -grad2;
  hessian.block(0,1,1,npar1) = grad1.transpose();
  hessian.block(0,1+npar1,1,npar2) = -grad2.transpose();

  grad1 *= prob_cpt1_;
  grad2 *= omp;
  gradient[0] = val1 - val2;
#ifndef CALIN_USE_EIGEN_REF
  gradient.segment(1,npar1) = grad1;
  gradient.segment(1+npar1,npar2) = grad2;
#endif
  
  return prob_cpt1_*val1 + omp*val2;
}

double TwoComponentPDF::error_up()
{
  return error_up_;
}

bool TwoComponentPDF::can_calculate_mean_and_variance()
{
  return false;
}

void TwoComponentPDF::mean_and_variance(double& mean, double& var)
{
  assert(0);
}
