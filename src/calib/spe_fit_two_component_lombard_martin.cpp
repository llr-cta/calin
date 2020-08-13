/*
   calin/calib/spe_fit_two_component_lombard_martin.cpp -- Stephen Fegan -- 2020-08-13

   Functions to do fit to multi-electron spectrum in the "single PE"
   domain. Two-component Lombard and Martin model.

   Copyright 2020, Stephen Fegan <sfegan@llr.in2p3.fr>
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
#include <cmath>
#include <limits>
#include <cassert>
#include <memory>

#include <fftw3.h>

#include <math/fftw_util.hpp>
#include <math/accumulator.hpp>
#include <calib/spe_fit.hpp>
#include <simulation/pmt.hpp>
#include <util/log.hpp>

using namespace calin::calib::spe_fit;

TwoComponentLombardMartinMES::
TwoComponentLombardMartinMES(
    const calin::ix::calib::spe_fit::TwoComponentLombardMartinMESConfig& config):
  MultiElectronSpectrum(), config_(config)
{
  // nothing to see here
}

TwoComponentLombardMartinMES::~TwoComponentLombardMartinMES()
{
  // nothing to see here
}

unsigned TwoComponentLombardMartinMES::num_parameters()
{
  unsigned nparam = 1; // light_intensity
  return nparam;
}

std::vector<calin::math::function::ParameterAxis>
TwoComponentLombardMartinMES::parameters()
{
  constexpr double inf = std::numeric_limits<double>::infinity();
  std::vector<math::function::ParameterAxis> pvec;
  pvec.push_back({ "light_intensity", "PE", true, 0, false, inf });
  return pvec;
}

Eigen::VectorXd TwoComponentLombardMartinMES::parameter_values()
{
  Eigen::VectorXd param(num_parameters());
  param[0] = intensity_pe_;
  return param;
}

void TwoComponentLombardMartinMES::set_parameter_values(ConstVecRef values)
{

}

bool TwoComponentLombardMartinMES::can_calculate_parameter_gradient()
{
  return false;
}

bool TwoComponentLombardMartinMES::can_calculate_parameter_hessian()
{
  return false;
}

double TwoComponentLombardMartinMES::pdf_ped(double x)
{
  return 0;
}

double TwoComponentLombardMartinMES::pdf_gradient_ped(double x, VecRef gradient)
{
  throw std::runtime_error("TwoComponentLombardMartinMES: pdf_gradient_ped not implemented");
  return 0;
}

double TwoComponentLombardMartinMES::pdf_gradient_hessian_ped(double x, VecRef gradient, MatRef hessian)
{
  throw std::runtime_error("TwoComponentLombardMartinMES: pdf_gradient_hessian_ped not implemented");
  return 0;
}

double TwoComponentLombardMartinMES::pdf_mes(double x)
{
  return 0;
}

double TwoComponentLombardMartinMES::pdf_gradient_mes(double x, VecRef gradient)
{
  throw std::runtime_error("TwoComponentLombardMartinMES: pdf_gradient_mes not implemented");
  return 0;
}

double TwoComponentLombardMartinMES::pdf_gradient_hessian_mes(double x, VecRef gradient, MatRef hessian)
{
  throw std::runtime_error("TwoComponentLombardMartinMES: pdf_gradient_hessian_mes not implemented");
  return 0;
}

double TwoComponentLombardMartinMES::intensity_pe()
{
  return 0;
}

double TwoComponentLombardMartinMES::ped_rms_dc()
{
  return 0;
}

double TwoComponentLombardMartinMES::ped_zero_dc()
{
  return 0;
}

double TwoComponentLombardMartinMES::ses_mean_dc()
{
  return 0;
}

double TwoComponentLombardMartinMES::ses_rms_pe()
{
  return 0;
}

calin::ix::calib::spe_fit::TwoComponentLombardMartinMESConfig
TwoComponentLombardMartinMES::default_config()
{
  calin::ix::calib::spe_fit::TwoComponentLombardMartinMESConfig config;

  auto* pmt = config.mutable_pmt();
  pmt->set_num_stage(7);
  pmt->set_total_gain(40000.0);
  pmt->set_stage_0_hi_gain(12);
  pmt->set_stage_0_hi_gain_rms_frac(0.0);
  pmt->set_stage_0_lo_gain(3);
  pmt->set_stage_0_lo_gain_rms_frac(0.0);
  pmt->set_stage_0_lo_prob(0.15);
  pmt->set_stage_n_gain_rms_frac(0.0);

  config.set_sensitivity(58.0/40000.0);
  config.set_free_total_gain(true);
  config.set_free_stage_0_hi_gain(false);
  config.set_free_stage_0_hi_gain_rms_frac(false);
  config.set_free_stage_0_lo_prob(false);
  config.set_free_stage_0_lo_gain(false);
  config.set_free_stage_0_lo_gain_rms_frac(false);
  config.set_free_stage_n_gain_rms_frac(false);
  config.set_intensity_rms_frac(0.0);
  config.set_free_intensity_rms_frac(false);
  config.set_on_off_ped_shift(0.0);
  config.set_free_on_off_ped_shift(false);
  config.set_precision(1e-10);

  return config;
}
