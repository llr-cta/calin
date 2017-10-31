/*

   calin/math/data_modeling.cpp -- Stephen Fegan -- 2017-04-07

   Data modeling functions for various data types

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

#include "math/data_modeling.hpp"
#include "math/accumulator.hpp"
#include "math/special.hpp"

using namespace calin::math::data_modeling;
using calin::math::special::SQR;

IID1DDataLikelihoodFunction::~IID1DDataLikelihoodFunction()
{
  if(adopt_pdf_)delete pdf_;
}

unsigned IID1DDataLikelihoodFunction::num_domain_axes()
{
  return pdf_->num_parameters();
}

std::vector<calin::math::function::DomainAxis> IID1DDataLikelihoodFunction::domain_axes()
{
  return pdf_->parameters();
}

double IID1DDataLikelihoodFunction::value(ConstVecRef x)
{
  pdf_->set_parameter_values(x);
  math::accumulator::LikelihoodAccumulator acc;
  unsigned nx = x_.size();
  for(unsigned ix=0; ix<nx; ix++)
  {
    const double w = w_[ix];
    if(w>0)
    {
      double pdf = pdf_->value_1d(x_[ix]);
      if(pdf<=0)continue;
      acc.accumulate(std::log(pdf)*w);
    }
  }
  return -acc.total();
}

bool IID1DDataLikelihoodFunction::can_calculate_gradient()
{
  return pdf_->can_calculate_parameter_gradient();
}

double IID1DDataLikelihoodFunction::value_and_gradient(ConstVecRef x, VecRef gradient)
{
  pdf_->set_parameter_values(x);
  gradient.resize(npar_);
  math::accumulator::LikelihoodAccumulator acc;
  std::vector<math::accumulator::LikelihoodAccumulator> gradient_acc(npar_);
  unsigned nx = x_.size();
  for(unsigned ix=0; ix<nx; ix++)
  {
    const double w = w_[ix];
    if(w>0)
    {
      double pdf = pdf_->value_and_parameter_gradient_1d(x_[ix], gradient);
      if(pdf<=0)continue;
      acc.accumulate(std::log(pdf)*w);
      for(unsigned ipar=0;ipar<npar_;ipar++)
        gradient_acc[ipar].accumulate(gradient(ipar)/pdf*w);
    }
  }
  for(unsigned ipar=0;ipar<npar_;ipar++)
    gradient(ipar) = -gradient_acc[ipar].total();
  return -acc.total();
}

bool IID1DDataLikelihoodFunction::can_calculate_hessian()
{
  return pdf_->can_calculate_parameter_hessian();
}

double IID1DDataLikelihoodFunction::
value_gradient_and_hessian(ConstVecRef x, VecRef gradient, MatRef hessian)
{
  pdf_->set_parameter_values(x);
  gradient.resize(npar_);
  hessian.resize(npar_,npar_);
  math::accumulator::LikelihoodAccumulator acc;
  std::vector<math::accumulator::LikelihoodAccumulator> gradient_acc(npar_);
  std::vector<math::accumulator::LikelihoodAccumulator>
      hessian_acc(npar_*(npar_+1)/2);
  unsigned nx = x_.size();
  for(unsigned ix=0; ix<nx; ix++)
  {
    const double w = w_[ix];
    if(w>0)
    {
      double pdf = pdf_->value_parameter_gradient_and_hessian_1d(x_[ix],
        gradient, hessian);
      if(pdf<=0)continue;
      acc.accumulate(std::log(pdf)*w);
      for(unsigned ipar=0;ipar<npar_;ipar++)
        gradient_acc[ipar].accumulate(gradient(ipar)/pdf*w);
      unsigned itri = 0;
      for(unsigned icol=0;icol<npar_;icol++)
        for(unsigned irow=icol;irow<npar_;irow++)
        {
          double summand = (hessian(icol,irow)
                            - gradient[icol]*gradient[irow]/pdf)/pdf;
          hessian_acc[itri++].accumulate(summand*w);
        }
    }
  }

  for(unsigned ipar=0;ipar<npar_;ipar++)
    gradient(ipar) = -gradient_acc[ipar].total();

  unsigned iacc = 0;
  for(unsigned icol=0;icol<npar_;icol++)
  {
    for(unsigned irow=icol;irow<npar_;irow++)
      hessian(icol,irow) = -hessian_acc[iacc++].total();
    for(unsigned irow=0;irow<icol;irow++)
      hessian(icol,irow) = hessian(irow,icol);
  }

  return -acc.total();
}

// -----------------------------------------------------------------------------
//
// M-Estimate version of Likelihood function fitter
//
// -----------------------------------------------------------------------------

IID1DDataMEstimateLikelihoodFunction::~IID1DDataMEstimateLikelihoodFunction()
{
  if(adopt_pdf_)delete pdf_;
  if(adopt_rho_)delete rho_;
}

unsigned IID1DDataMEstimateLikelihoodFunction::num_domain_axes()
{
  return pdf_->num_parameters();
}

std::vector<calin::math::function::DomainAxis> IID1DDataMEstimateLikelihoodFunction::domain_axes()
{
  return pdf_->parameters();
}

double IID1DDataMEstimateLikelihoodFunction::value(ConstVecRef x)
{
  pdf_->set_parameter_values(x);
  math::accumulator::LikelihoodAccumulator acc;
  unsigned nx = x_.size();
  for(unsigned ix=0; ix<nx; ix++)
  {
    const double w = w_[ix];
    if(w>0)
    {
      double pdf = pdf_->value_1d(x_[ix]);
      if(pdf<=0)continue;
      acc.accumulate(rho_->value_1d(-std::log(pdf))*w);
    }
  }
  return acc.total();
}

bool IID1DDataMEstimateLikelihoodFunction::can_calculate_gradient()
{
  return pdf_->can_calculate_parameter_gradient() and
    rho_->can_calculate_gradient();
}

double IID1DDataMEstimateLikelihoodFunction::value_and_gradient(ConstVecRef x, VecRef gradient)
{
  pdf_->set_parameter_values(x);
  gradient.resize(npar_);
  math::accumulator::LikelihoodAccumulator acc;
  std::vector<math::accumulator::LikelihoodAccumulator> gradient_acc(npar_);
  unsigned nx = x_.size();
  for(unsigned ix=0; ix<nx; ix++)
  {
    const double w = w_[ix];
    if(w>0)
    {
      double pdf = pdf_->value_and_parameter_gradient_1d(x_[ix], gradient);
      if(pdf<=0)continue;
      double drho_dx;
      double rho = rho_->value_and_gradient_1d(-std::log(pdf), drho_dx);
      acc.accumulate(rho*w);
      for(unsigned ipar=0;ipar<npar_;ipar++)
        gradient_acc[ipar].accumulate(drho_dx*gradient(ipar)/pdf*w);
    }
  }
  for(unsigned ipar=0;ipar<npar_;ipar++)
    gradient(ipar) = -gradient_acc[ipar].total();
  return acc.total();
}

bool IID1DDataMEstimateLikelihoodFunction::can_calculate_hessian()
{
  return pdf_->can_calculate_parameter_hessian() and
    rho_->can_calculate_hessian();
}

double IID1DDataMEstimateLikelihoodFunction::
value_gradient_and_hessian(ConstVecRef x, VecRef gradient, MatRef hessian)
{
  pdf_->set_parameter_values(x);
  gradient.resize(npar_);
  hessian.resize(npar_,npar_);
  math::accumulator::LikelihoodAccumulator acc;
  std::vector<math::accumulator::LikelihoodAccumulator> gradient_acc(npar_);
  std::vector<math::accumulator::LikelihoodAccumulator>
      hessian_acc(npar_*(npar_+1)/2);
  unsigned nx = x_.size();
  for(unsigned ix=0; ix<nx; ix++)
  {
    const double w = w_[ix];
    if(w>0)
    {
      double pdf = pdf_->value_parameter_gradient_and_hessian_1d(x_[ix],
        gradient, hessian);
      if(pdf<=0)continue;
      double drho_dx;
      double d2rho_dx2;
      double rho = rho_->value_gradient_and_hessian_1d(-std::log(pdf), drho_dx, d2rho_dx2);
      acc.accumulate(rho*w);
      for(unsigned ipar=0;ipar<npar_;ipar++)
        gradient_acc[ipar].accumulate(drho_dx*gradient(ipar)/pdf*w);

      unsigned itri = 0;
      for(unsigned icol=0;icol<npar_;icol++)
        for(unsigned irow=icol;irow<npar_;irow++)
        {
          double summand = (drho_dx*hessian(icol,irow)
            - (d2rho_dx2 + drho_dx)*gradient[icol]*gradient[irow]/pdf)/pdf;
          hessian_acc[itri++].accumulate(summand*w);
        }
    }
  }

  for(unsigned ipar=0;ipar<npar_;ipar++)
    gradient(ipar) = -gradient_acc[ipar].total();

  unsigned iacc = 0;
  for(unsigned icol=0;icol<npar_;icol++)
  {
    for(unsigned irow=icol;irow<npar_;irow++)
      hessian(icol,irow) = -hessian_acc[iacc++].total();
    for(unsigned irow=0;irow<icol;irow++)
      hessian(icol,irow) = hessian(irow,icol);
  }

  return acc.total();
}

template<typename F> double integrate_binned_in_extendable_window(
  F f, calin::math::function::ParameterizableSingleAxisFunction* pdf,
  double x0, double x1, double dx, double norm_accuracy = 1e-6)
{
  calin::math::accumulator::LikelihoodAccumulator norm_acc;
  calin::math::accumulator::LikelihoodAccumulator integrand_acc;
  double p = std::numeric_limits<double>::infinity();
  double x;
  if(x0 > x1)std::swap(x0,x1);
  dx = std::fabs(dx);
  for(x=x0; x<=x1; x+=dx) {
    p = pdf->value_1d(x);
    if(p<=0)continue;
    norm_acc.accumulate(p);
    integrand_acc.accumulate(f(x,p) * p);
  }
  double p0 = std::numeric_limits<double>::infinity();
  double p1 = p;
  x1 = x;
  double norm_target = (1.0-norm_accuracy)/dx;
  unsigned nzero = 0;
  while(norm_acc.total()<norm_target and nzero<100) {
    double* pp = nullptr;
    if(p0 > p1 or (std::max(p0,p1)==0.0 and nzero%2==0)) {
      x0 -= dx;
      x = x0;
      pp = &p0;
    } else {
      x1 += dx;
      x = x1;
      pp = &p1;
    }
    *pp = std::max(pdf->value_1d(x),0.0);
#if 0
    calin::util::log::LOG(calin::util::log::INFO) << x << ' ' << *pp << ' '
      << nzero << ' ' << std::scientific << 1.0-norm_acc.total() << ' '
      << 1.0-norm_target;
#endif
    if(*pp==0){
      ++nzero;
      continue;
    }
    norm_acc.accumulate(p);
    integrand_acc.accumulate(f(x,p) * p);
    nzero = 0;
  }
  if(norm_acc.total() < norm_target)
    throw std::runtime_error("Integration did not converge, final normalization: 1.0-"
      + std::to_string(1.0 - norm_acc.total()*dx));
  return integrand_acc.total()*dx;
}

double IID1DDataMEstimateLikelihoodFunction::
expectation_value(ConstVecRef x, double norm_tolerance)
{
  if(dx_ == 0)
    throw std::runtime_error("Calculation of unbinned expectation_value not supported");
  pdf_->set_parameter_values(x);
  return integrate_binned_in_extendable_window([this](double x, double p) {
      return rho_->value_1d(-std::log(p)); }, pdf_,
    x_.front(), x_.back()+0.5*dx_, dx_, norm_tolerance);
}

double IID1DDataMEstimateLikelihoodFunction::
expectation_variance(ConstVecRef x, double norm_tolerance)
{
  if(dx_ == 0)
    throw std::runtime_error("Calculation of unbinned expectation_variance not supported");
  double exp_val = expectation_value(x, norm_tolerance);
  double exp_var =
    integrate_binned_in_extendable_window([this](double x, double p) {
        return SQR(rho_->value_1d(-std::log(p))); }, pdf_,
      x_.front(), x_.back()+0.5*dx_, dx_, norm_tolerance);
  exp_var -= SQR(exp_val);
  return exp_var;
}

// -----------------------------------------------------------------------------
//
// Chi-squared version of fitter
//
// -----------------------------------------------------------------------------

IID1DDataChi2Function::~IID1DDataChi2Function()
{
  if(adopt_pdf_)delete pdf_;
}

unsigned IID1DDataChi2Function::num_domain_axes()
{
  return pdf_->num_parameters()+1;
}

std::vector<calin::math::function::DomainAxis> IID1DDataChi2Function::domain_axes()
{
  std::vector<calin::math::function::DomainAxis> params = pdf_->parameters();
  params.insert(params.begin(), { "normalization", "events", true, 0, false, 0 });
  return params;
}

double IID1DDataChi2Function::value(ConstVecRef x)
{
  double norm;
  function::assign_parameters(x.data(), norm);
  pdf_->set_parameter_values(x.segment(1,npar_));
  math::accumulator::LikelihoodAccumulator acc;
  unsigned nx = x_.size();
  for(unsigned ix=0; ix<nx; ix++)
  {
    const double w = w_[ix];
    if(w>0)
    {
      double pdf = pdf_->value_1d(x_[ix]);
      acc.accumulate(SQR(norm*pdf-w)/w);
    }
  }
  return acc.total();
}

bool IID1DDataChi2Function::can_calculate_gradient()
{
  return pdf_->can_calculate_parameter_gradient();
}

double IID1DDataChi2Function::value_and_gradient(ConstVecRef x, VecRef gradient)
{
  double norm;
  function::assign_parameters(x.data(), norm);
  pdf_->set_parameter_values(x.segment(1,npar_));
  gradient.resize(npar_+1);
  math::accumulator::LikelihoodAccumulator acc;
  std::vector<math::accumulator::LikelihoodAccumulator> gradient_acc(npar_+1);
  Eigen::VectorXd pdf_gradient(npar_);
  unsigned nx = x_.size();
  for(unsigned ix=0; ix<nx; ix++)
  {
    const double w = w_[ix];
    if(w>0)
    {
      double pdf = pdf_->value_and_parameter_gradient_1d(x_[ix], pdf_gradient);
      double diff = norm*pdf-w;
      acc.accumulate(SQR(diff)/w);
      gradient_acc[0].accumulate(2.0*diff*pdf/w);
      for(unsigned ipar=0;ipar<npar_;ipar++)
        gradient_acc[ipar+1].accumulate(2*diff*norm*pdf_gradient(ipar)/w);
    }
  }
  for(unsigned ipar=0;ipar<npar_+1;ipar++)
    gradient(ipar) = gradient_acc[ipar].total();
  return acc.total();
}

bool IID1DDataChi2Function::can_calculate_hessian()
{
  return pdf_->can_calculate_parameter_hessian();
}

double IID1DDataChi2Function::
value_gradient_and_hessian(ConstVecRef x, VecRef gradient, MatRef hessian)
{
  double norm;
  function::assign_parameters(x.data(), norm);
  pdf_->set_parameter_values(x.segment(1,npar_));
  gradient.resize(npar_+1);
  hessian.resize(npar_+1, npar_+1);

  math::accumulator::LikelihoodAccumulator acc;
  std::vector<math::accumulator::LikelihoodAccumulator> gradient_acc(npar_+1);
  std::vector<math::accumulator::LikelihoodAccumulator>
    hessian_acc((npar_+1)*(npar_+2)/2);

  Eigen::VectorXd pdf_gradient(npar_);
  Eigen::MatrixXd pdf_hessian(npar_, npar_);
  unsigned nx = x_.size();
  for(unsigned ix=0; ix<nx; ix++)
  {
    const double w = w_[ix];
    if(w>0)
    {
      double pdf = pdf_->value_parameter_gradient_and_hessian_1d(x_[ix],
        pdf_gradient, pdf_hessian);
      double diff = norm*pdf-w;
      acc.accumulate(SQR(diff)/w);
      gradient_acc[0].accumulate(2.0*diff*pdf/w);
      for(unsigned ipar=0;ipar<npar_;ipar++)
        gradient_acc[ipar+1].accumulate(2.0*diff*norm*pdf_gradient(ipar)/w);
      hessian_acc[0].accumulate(2.0*pdf*pdf/w);
      for(unsigned ipar=0;ipar<npar_;ipar++)
        hessian_acc[ipar+1].accumulate(2.0*(pdf*norm+diff)*pdf_gradient(ipar)/w);
      unsigned itri = npar_+1;
      for(unsigned ipar=0;ipar<npar_;ipar++)
        for(unsigned jpar=ipar;jpar<npar_;jpar++)
        {
          double summand = 2*norm*(diff*pdf_hessian(ipar,jpar)
            + norm*pdf_gradient(ipar)*pdf_gradient(jpar))/w;
          hessian_acc[itri++].accumulate(summand);
        }
    }
  }

  for(unsigned ipar=0;ipar<npar_+1;ipar++)
    gradient(ipar) = gradient_acc[ipar].total();

  unsigned iacc = 0;
  for(unsigned icol=0;icol<npar_+1;icol++)
  {
    for(unsigned irow=icol;irow<npar_+1;irow++)
      hessian(icol,irow) = hessian_acc[iacc++].total();
    for(unsigned irow=0;irow<icol;irow++)
      hessian(icol,irow) = hessian(irow,icol);
  }

  return acc.total();
}
