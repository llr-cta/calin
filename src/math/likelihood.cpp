/*

   calin/math/likelihood.cpp -- Stephen Fegan -- 2017-04-07

   Likelihood functions for various data types

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

#include "math/likelihood.hpp"
#include "math/accumulator.hpp"

using namespace calin::math::likelihood;

IID1DValueLikelihoodFunction::~IID1DValueLikelihoodFunction()
{
  if(adopt_pdf_)delete pdf_;
}

unsigned IID1DValueLikelihoodFunction::num_domain_axes()
{
  return pdf_->num_parameters();
}

std::vector<calin::math::function::DomainAxis> IID1DValueLikelihoodFunction::domain_axes()
{
  return pdf_->parameters();
}

double IID1DValueLikelihoodFunction::value(ConstVecRef x)
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

bool IID1DValueLikelihoodFunction::can_calculate_gradient()
{
  return pdf_->can_calculate_parameter_gradient();
}

double IID1DValueLikelihoodFunction::value_and_gradient(ConstVecRef x, VecRef gradient)
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

bool IID1DValueLikelihoodFunction::can_calculate_hessian()
{
  return pdf_->can_calculate_parameter_hessian();
}

double IID1DValueLikelihoodFunction::
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
