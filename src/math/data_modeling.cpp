/*

   calin/math/data_modeling.cpp -- Stephen Fegan -- 2017-04-07

   Data modeling functions for various data types

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
  assert(x.size() == 1);
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
  assert(x.size() == 1);
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
  assert(x.size() == 1);

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
  assert(x.size() == 1);
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
  assert(x.size() == 1);
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
  assert(x.size() == 1);

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
            + (d2rho_dx2 - drho_dx)*gradient[icol]*gradient[irow]/pdf)/pdf;
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
  assert(x.size() == 1);
  function::assign_parameters(x.data(), norm_);
  pdf_->set_parameter_values(x.segment(1,npar_));
  math::accumulator::LikelihoodAccumulator acc;
  unsigned nx = x_.size();
  for(unsigned ix=0; ix<nx; ix++)
  {
    const double w = w_[ix];
    if(w>0)
    {
      double pdf = pdf_->value_1d(x_[ix]);
      acc.accumulate(SQR(norm_*pdf-w)/w);
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
  assert(x.size() == 1);
  function::assign_parameters(x.data(), norm_);
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
      double diff = norm_*pdf-w;
      acc.accumulate(SQR(diff)/w);
      gradient_acc[0].accumulate(2.0*diff*pdf/w);
      for(unsigned ipar=0;ipar<npar_;ipar++)
        gradient_acc[ipar+1].accumulate(2*diff*norm_*pdf_gradient(ipar)/w);
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
  assert(x.size() == 1);
  function::assign_parameters(x.data(), norm_);
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
      double diff = norm_*pdf-w;
      acc.accumulate(SQR(diff)/w);
      gradient_acc[0].accumulate(2.0*diff*pdf/w);
      for(unsigned ipar=0;ipar<npar_;ipar++)
        gradient_acc[ipar+1].accumulate(2.0*diff*norm_*pdf_gradient(ipar)/w);
      hessian_acc[0].accumulate(2.0*pdf/w);
      for(unsigned ipar=0;ipar<npar_;ipar++)
        hessian_acc[ipar+1].accumulate(2.0*(pdf*norm_+diff)*pdf_gradient(ipar)/w);
      unsigned itri = npar_+1;
      for(unsigned ipar=0;ipar<npar_;ipar++)
        for(unsigned jpar=ipar;jpar<npar_;jpar++)
        {
          double summand = 2*norm_*(diff*pdf_hessian(ipar,jpar)
            + norm_*pdf_gradient(ipar)*pdf_gradient(jpar))/w;
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
