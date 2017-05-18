/*
   calin/calib/spe_fit_fast_single_value_general_poisson.cpp
                                                -- Stephen Fegan -- 2017-04-04

   Functions to do fit to multi-electron spectrum in the "single PE"
   domain. Fast Poisson-Gaussian model.

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

#include <cmath>
#include <limits>
#include <cassert>

#include <math/accumulator.hpp>
#include <calib/spe_fit.hpp>

using namespace calin::math;
using namespace calin::calib::spe_fit;

using calin::math::function::assign_parameters;

// ============================================================================
//
// PoissonGaussMES -Fast PG model for single PE regime where maximum
// number of PEs is specified
//
// ============================================================================

FastSingleValueGeneralPoissonMES::
FastSingleValueGeneralPoissonMES(GeneralPoissonMES* mes, bool adopt_mes):
  MultiElectronSpectrum(), mes_(mes), adopt_mes_(adopt_mes),
  nes_pmf_(mes_->nsample(), mes_->num_electrons_in_model()+1),
  ped_pmf_(mes_->off_pedestal_spectrum()),
  nes_weight_(mes_->num_electrons_in_model()+1),
  nes_weight_deriv_(mes_->num_electrons_in_model()+1)
{
  update_from_general_mes();
}

void FastSingleValueGeneralPoissonMES::update_from_general_mes()
{
  Eigen::VectorXd param(1);
  param << mes_->intensity_pe();
  set_parameter_values(param);
  ped_pmf_ = mes_->off_pedestal_spectrum();
  for(unsigned ipe=0; ipe<nes_pmf_.cols(); ipe++)
    nes_pmf_.col(ipe) = mes_->n_electron_spectrum_with_pedestal(ipe);
}

FastSingleValueGeneralPoissonMES::~FastSingleValueGeneralPoissonMES()
{
  if(adopt_mes_)delete mes_;
}

unsigned FastSingleValueGeneralPoissonMES::num_parameters()
{
  return 1;
}

std::vector<calin::math::function::ParameterAxis>
FastSingleValueGeneralPoissonMES::parameters()
{
  constexpr double inf = std::numeric_limits<double>::infinity();
  return {{ "light_intensity", "PE", true, 0, false, inf }};
}

Eigen::VectorXd FastSingleValueGeneralPoissonMES::parameter_values()
{
  Eigen::VectorXd param(1);
  param << intensity_pe_;
  return param;
}

void FastSingleValueGeneralPoissonMES::
set_parameter_values(ConstVecRef values)
{
  assign_parameters(values.data(), intensity_pe_);
#if 0
  //double log_intensity = std::log(intensity_pe_);
  double weight = std::exp(-intensity_pe_);
  nes_weight_[0] = weight;
  nes_weight_deriv_[0] = -weight;
  for(unsigned ipe=1;ipe<nes_weight_.size();ipe++)
  {
    nes_weight_deriv_[ipe] = weight;
    weight *= intensity_pe_/double(ipe);
    nes_weight_[ipe] = weight;
    nes_weight_deriv_[ipe] -= weight;
  }
#else
  double log_intensity = std::log(intensity_pe_);
  for(unsigned ipe=0;ipe<nes_weight_.size();ipe++)
  {
    double dbl_n { double(ipe) };
    double weight =
      std::exp(dbl_n*log_intensity - intensity_pe_ - lgamma(dbl_n+1.0));
    nes_weight_[ipe] = weight;
    nes_weight_deriv_[ipe] = weight*(dbl_n/intensity_pe_ - 1);
  }
#endif
}

bool FastSingleValueGeneralPoissonMES::can_calculate_parameter_gradient()
{
  return true;
}

bool FastSingleValueGeneralPoissonMES::can_calculate_parameter_hessian()
{
  return false;
}

double FastSingleValueGeneralPoissonMES::pdf_ped(double x)
{
  return std::max(ped_pmf_[mes_->ibin(x)],0.0);
}

double FastSingleValueGeneralPoissonMES::pdf_gradient_ped(double x, VecRef gradient)
{
  gradient.resize(1);
  gradient << 0.0;
  return std::max(ped_pmf_[mes_->ibin(x)],0.0);
}

double FastSingleValueGeneralPoissonMES::
pdf_gradient_hessian_ped(double x, VecRef gradient, MatRef hessian)
{
  throw std::runtime_error(
    "Attempt to call FastSingleValueGeneralPoissonMES::pdf_gradient_hessian_ped");
  return 0;
}

double FastSingleValueGeneralPoissonMES::pdf_mes(double x)
{
  return std::max(nes_pmf_.row(mes_->ibin(x)).dot(nes_weight_), 0.0);
}

double FastSingleValueGeneralPoissonMES::pdf_gradient_mes(double x, VecRef gradient)
{
  gradient.resize(1);
  auto therow = nes_pmf_.row(mes_->ibin(x));
  gradient << therow.dot(nes_weight_deriv_);
  return std::max(therow.dot(nes_weight_), 0.0);
}

double FastSingleValueGeneralPoissonMES::
pdf_gradient_hessian_mes(double x, VecRef gradient, MatRef hessian)
{
  throw std::runtime_error(
    "Attempt to call FastSingleValueGeneralPoissonMES::pdf_gradient_hessian_mes");
  return 0;
}

double FastSingleValueGeneralPoissonMES::ped_rms_dc()
{
  throw std::runtime_error(
    "Attempt to call FastSingleValueGeneralPoissonMES::ped_rms_dc");
  return 0;
}

double FastSingleValueGeneralPoissonMES::ped_zero_dc()
{
  throw std::runtime_error(
    "Attempt to call FastSingleValueGeneralPoissonMES::ped_zero_dc");
  return 0;
}

double FastSingleValueGeneralPoissonMES::ses_mean_dc()
{
  throw std::runtime_error(
    "Attempt to call FastSingleValueGeneralPoissonMES::ses_mean_dc");
  return 0;
}

double FastSingleValueGeneralPoissonMES::ses_rms_pe()
{
  throw std::runtime_error(
    "Attempt to call FastSingleValueGeneralPoissonMES::ses_rms_pe");
  return 0;
}
