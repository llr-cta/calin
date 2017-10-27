/*
   calin/calib/spe_fit.cpp -- Stephen Fegan -- 2015-03-01

   Functions to do fit to multi-electron spectrum in the "single PE"
   domain. Base class and Likelihood.

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <iostream>
#include <cmath>
#include <limits>
#include <cassert>
#include <memory>

#include <math/accumulator.hpp>
#include <math/special.hpp>
#include <calib/spe_fit.hpp>
#include <calib/pmt_model_pg.hpp>
#include <io/log.hpp>

using namespace calin::math;
using namespace calin::calib::spe_fit;
using namespace calin::calib::pmt_model_pg;
using namespace calin::io::log;

using calin::math::special::SQR;
using calin::math::function::assign_parameters;

#if 0
SingleElectronSpectrum::~SingleElectronSpectrum()
{
  // nothing to see here
}
#endif

// ============================================================================
//
// MultiElectronSpectrum base class
//
// ============================================================================

MultiElectronSpectrum::~MultiElectronSpectrum()
{
  // nothing to see here
}

// ============================================================================
//
// SPELikelihood -- calculate likelihood for optimizer
//
// ============================================================================

SPELikelihood::SPELikelihood(MultiElectronSpectrum& mes_model,
                             const math::histogram::SimpleHist& mes_data):
    function::MultiAxisFunction(), mes_model_(&mes_model),
    npar_(mes_model.num_parameters()),
    mes_data_(mes_data), has_ped_data_(false), ped_data_(1.0)
{
  // nothing to see here
}

SPELikelihood::SPELikelihood(MultiElectronSpectrum& mes_model,
                             const math::histogram::SimpleHist& mes_data,
                             const math::histogram::SimpleHist& ped_data):
    function::MultiAxisFunction(), mes_model_(&mes_model),
    npar_(mes_model.num_parameters()),
    mes_data_(mes_data), has_ped_data_(true), ped_data_(ped_data)
{
  // nothing to see here
}

SPELikelihood::~SPELikelihood()
{
  // nothing to see here
}

unsigned SPELikelihood::num_domain_axes()
{
  return mes_model_->num_parameters();
}

auto SPELikelihood::domain_axes() -> std::vector<math::function::DomainAxis>
{
  return mes_model_->parameters();
}

double SPELikelihood::value(ConstVecRef x)
{
  mes_model_->set_parameter_values(x);
  math::accumulator::LikelihoodAccumulator acc;
  for(auto& ibin : mes_data_)
    if(ibin.weight())
    {
      double pdf = mes_model_->pdf_mes(ibin.xval_center());
      if(pdf<=0)continue;
      acc.accumulate(std::log(pdf)*ibin.weight());
    }

  if(has_ped_data_)
    for(auto& ibin : ped_data_)
      if(ibin.weight())
      {
        double pdf = mes_model_->pdf_ped(ibin.xval_center());
        if(pdf<=0)continue;
        acc.accumulate(std::log(pdf)*ibin.weight());
      }
  return -acc.total();
}

bool SPELikelihood::can_calculate_gradient()
{
  return mes_model_->can_calculate_parameter_gradient();
}

bool SPELikelihood::can_calculate_hessian()
{
  return mes_model_->can_calculate_parameter_hessian();
}

double SPELikelihood::value_and_gradient(ConstVecRef x, VecRef gradient)
{
  mes_model_->set_parameter_values(x);
  gradient.resize(npar_);
  math::accumulator::LikelihoodAccumulator acc;
  std::vector<math::accumulator::LikelihoodAccumulator> gradient_acc(npar_);
  for(auto& ibin : mes_data_)
    if(ibin.weight())
    {
      double pdf = mes_model_->pdf_gradient_mes(ibin.xval_center(), gradient);
      if(pdf<=0)continue;
      acc.accumulate(std::log(pdf)*ibin.weight());
      for(unsigned ipar=0;ipar<npar_;ipar++)
        gradient_acc[ipar].accumulate(gradient(ipar)/pdf*ibin.weight());
    }

  if(has_ped_data_)
    for(auto& ibin : ped_data_)
      if(ibin.weight())
      {
        double pdf = mes_model_->pdf_gradient_ped(ibin.xval_center(), gradient);
        if(pdf<=0)continue;
        acc.accumulate(std::log(pdf)*ibin.weight());
        for(unsigned ipar=0;ipar<npar_;ipar++)
          gradient_acc[ipar].accumulate(gradient(ipar)/pdf*ibin.weight());
      }

  for(unsigned ipar=0;ipar<npar_;ipar++)
    gradient(ipar) = -gradient_acc[ipar].total();
  return -acc.total();
}

double SPELikelihood::
value_gradient_and_hessian(ConstVecRef x, VecRef gradient, MatRef hessian)
{
  mes_model_->set_parameter_values(x);
  gradient.resize(npar_);
  hessian.resize(npar_,npar_);
  math::accumulator::LikelihoodAccumulator acc;
  std::vector<math::accumulator::LikelihoodAccumulator> gradient_acc(npar_);
  std::vector<math::accumulator::LikelihoodAccumulator>
      hessian_acc(npar_*(npar_+1)/2);
  for(auto& ibin : mes_data_)
    if(ibin.weight())
    {
      double pdf =
          mes_model_->pdf_gradient_hessian_mes(ibin.xval_center(),
                                               gradient, hessian);
      if(pdf<=0)continue;
      acc.accumulate(std::log(pdf)*ibin.weight());
      for(unsigned ipar=0;ipar<npar_;ipar++)
        gradient_acc[ipar].accumulate(gradient(ipar)/pdf*ibin.weight());
      unsigned itri = 0;
      for(unsigned icol=0;icol<npar_;icol++)
        for(unsigned irow=icol;irow<npar_;irow++)
        {
          double summand = (hessian(icol,irow)
                            - gradient[icol]*gradient[irow]/pdf)/pdf;
          hessian_acc[itri++].accumulate(summand*ibin.weight());
        }
    }

  if(has_ped_data_)
    for(auto& ibin : ped_data_)
      if(ibin.weight())
      {
        double pdf =
            mes_model_->pdf_gradient_hessian_ped(ibin.xval_center(),
                                                 gradient, hessian);
        if(pdf<=0)continue;
        acc.accumulate(std::log(pdf)*ibin.weight());
        for(unsigned ipar=0;ipar<npar_;ipar++)
          gradient_acc[ipar].accumulate(gradient(ipar)/pdf*ibin.weight());
        unsigned iacc = 0;
        for(unsigned icol=0;icol<npar_;icol++)
          for(unsigned irow=icol;irow<npar_;irow++)
          {
            double summand = (hessian(icol,irow)
                            - gradient[icol]*gradient[irow]/pdf)/pdf;
            hessian_acc[iacc++].accumulate(summand*ibin.weight());
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
