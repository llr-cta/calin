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
#include <math/special.hpp>

using namespace calin::math::fftw_util;
using namespace calin::calib::spe_fit;
using calin::math::special::round_up_power_of_two;
using namespace calin::util::log;

TwoComponentLombardMartinMES::
TwoComponentLombardMartinMES(double x0, unsigned npoint,
    calin::math::function::ParameterizableSingleAxisFunction* ped_pdf,
    const calin::ix::calib::spe_fit::TwoComponentLombardMartinMESConfig& config,
    bool adopt_ped_pdf):
  MultiElectronSpectrum(), config_(config), x0_(x0), npoint_(npoint),
  ped_pdf_(ped_pdf), adopt_ped_pdf_(adopt_ped_pdf), mes_pmf_(npoint), off_pmf_(npoint)
{
  calculate_mes();
}

TwoComponentLombardMartinMES::~TwoComponentLombardMartinMES()
{
  if(adopt_ped_pdf_)delete ped_pdf_;
}

unsigned TwoComponentLombardMartinMES::num_parameters()
{
  unsigned nparam = 1; // light_intensity
  if(config_.free_intensity_rms_frac()) { ++nparam; }
  if(config_.free_total_gain()) { ++nparam; }
  if(config_.free_stage_0_lo_prob()) { ++nparam; }
  if(config_.free_stage_0_hi_gain()) { ++nparam; }
  if(config_.free_stage_0_lo_gain()) { ++nparam; }
  if(config_.free_stage_0_hi_gain_rms_frac()) { ++nparam; }
  if(config_.free_stage_0_lo_gain_rms_frac()) { ++nparam; }
  if(config_.free_stage_n_gain_rms_frac()) { ++nparam; }
  if(ped_pdf_) {
    if(config_.free_on_off_ped_shift()) { ++nparam; }
    nparam += ped_pdf_->num_parameters();
  }
  return nparam;
}

std::vector<calin::math::function::ParameterAxis>
TwoComponentLombardMartinMES::parameters()
{
  constexpr double inf = std::numeric_limits<double>::infinity();
  std::vector<math::function::ParameterAxis> pvec;
  pvec.push_back({ "light_intensity", "PE", true, 0, false, inf });
  if(config_.free_intensity_rms_frac()) {
    pvec.push_back({ "light_intensity_rms_frac", "1", true, 0, false, inf }); }
  if(config_.free_total_gain()) {
    pvec.push_back({ "total_gain", "DC", true, 0, false, inf }); }
  if(config_.free_stage_0_lo_prob()) {
    pvec.push_back({ "stage_0_lo_prob", "1", true, 0, true, 1 }); }
  if(config_.free_stage_0_hi_gain()) {
    pvec.push_back({ "stage_0_hi_gain", "1", true, 0, false, inf }); }
  if(config_.free_stage_0_lo_gain()) {
    pvec.push_back({ "stage_0_lo_gain", "1", true, 0, false, inf }); }
  if(config_.free_stage_0_hi_gain_rms_frac()) {
    pvec.push_back({ "stage_0_hi_gain_rms_frac", "1", true, 0, false, inf }); }
  if(config_.free_stage_0_lo_gain_rms_frac()) {
    pvec.push_back({ "stage_0_lo_gain_rms_frac", "1", true, 0, false, inf }); }
  if(config_.free_stage_n_gain_rms_frac()) {
    pvec.push_back({ "stage_n_gain_rms_frac", "1", true, 0, false, inf }); }
  if(ped_pdf_) {
    if(config_.free_on_off_ped_shift()) {
      pvec.push_back({ "on_off_ped_shift", "1", false, -inf, false, inf }); }
    std::vector<math::function::ParameterAxis> pped { ped_pdf_->parameters() };
    for(auto& ip : pped)ip.name = std::string("ped.") + ip.name;
    pvec.insert(pvec.end(), pped.begin(), pped.end());
  }
  return pvec;
}

Eigen::VectorXd TwoComponentLombardMartinMES::parameter_values()
{
  Eigen::VectorXd param(num_parameters());
  unsigned iparam = 0;
  param[iparam++] = intensity_pe_;
  if(config_.free_intensity_rms_frac()) {
    param[iparam++] = config_.intensity_rms_frac(); }
  if(config_.free_total_gain()) {
    param[iparam++] = config_.pmt().total_gain()*config_.sensitivity(); }
  if(config_.free_stage_0_lo_prob()) {
    param[iparam++] = config_.pmt().stage_0_lo_prob(); }
  if(config_.free_stage_0_hi_gain()) {
    param[iparam++] = config_.pmt().stage_0_hi_gain(); }
  if(config_.free_stage_0_lo_gain()) {
    param[iparam++] = config_.pmt().stage_0_lo_gain(); }
  if(config_.free_stage_0_hi_gain_rms_frac()) {
    param[iparam++] = config_.pmt().stage_0_hi_gain_rms_frac(); }
  if(config_.free_stage_0_lo_gain_rms_frac()) {
    param[iparam++] = config_.pmt().stage_0_lo_gain_rms_frac(); }
  if(config_.free_stage_n_gain_rms_frac()) {
    param[iparam++] = config_.pmt().stage_n_gain_rms_frac(); }
  if(ped_pdf_) {
    if(config_.free_on_off_ped_shift()) {
      param[iparam++] = config_.on_off_ped_shift(); }
    unsigned num_ped_params = ped_pdf_->num_parameters();
    param.segment(iparam,num_ped_params) = ped_pdf_->parameter_values();
    iparam += num_ped_params;
  }
  return param;
}

void TwoComponentLombardMartinMES::set_parameter_values(ConstVecRef values)
{
  verify_set_parameter_values(values, "TwoComponentLombardMartinMES");
  unsigned iparam = 0;
  intensity_pe_ = values[iparam++];
  if(config_.free_intensity_rms_frac()) {
    config_.set_intensity_rms_frac(values[iparam++]); }
  if(config_.free_total_gain()) {
    config_.mutable_pmt()->set_total_gain(values[iparam++]/config_.sensitivity()); }
  if(config_.free_stage_0_lo_prob()) {
    config_.mutable_pmt()->set_stage_0_lo_prob(values[iparam++]); }
  if(config_.free_stage_0_hi_gain()) {
    config_.mutable_pmt()->set_stage_0_hi_gain(values[iparam++]); }
  if(config_.free_stage_0_lo_gain()) {
    config_.mutable_pmt()->set_stage_0_lo_gain(values[iparam++]); }
  if(config_.free_stage_0_hi_gain_rms_frac()) {
    config_.mutable_pmt()->set_stage_0_hi_gain_rms_frac(values[iparam++]); }
  if(config_.free_stage_0_lo_gain_rms_frac()) {
    config_.mutable_pmt()->set_stage_0_lo_gain_rms_frac(values[iparam++]); }
  if(config_.free_stage_n_gain_rms_frac()) {
    config_.mutable_pmt()->set_stage_n_gain_rms_frac(values[iparam++]); }
  if(ped_pdf_) {
    if(config_.free_on_off_ped_shift()) {
      config_.set_on_off_ped_shift(values[iparam++]); }
    unsigned num_ped_params = ped_pdf_->num_parameters();
    ped_pdf_->set_parameter_values(values.segment(iparam,num_ped_params));
    iparam += num_ped_params;
  }
  calculate_mes();
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
  if(ped_pdf_ == nullptr) {
    throw std::runtime_error("TwoComponentLombardMartinMES::pdf_ped: no pedestal pdf supplied");
  }
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

  config.set_sensitivity(58.0/40000.0);
  config.set_dx(1.0);

  auto* pmt = config.mutable_pmt();
  pmt->set_num_stage(7);
  pmt->set_total_gain(40000.0);
  pmt->set_stage_0_hi_gain(12);
  pmt->set_stage_0_hi_gain_rms_frac(0.0);
  pmt->set_stage_0_lo_gain(3);
  pmt->set_stage_0_lo_gain_rms_frac(0.0);
  pmt->set_stage_0_lo_prob(0.15);
  pmt->set_stage_n_gain_rms_frac(0.0);

  config.set_free_total_gain(true);
  config.set_free_stage_0_hi_gain(true);
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

void TwoComponentLombardMartinMES::calculate_mes()
{
  unsigned mes_npoint = round_up_power_of_two(double(npoint_+1)/config_.sensitivity()+1);

  uptr_fftw_data ped_hc_dft { nullptr, fftw_free };
  if(ped_pdf_ != nullptr) {
    // FFT of pedestal
    ped_hc_dft.reset(fftw_alloc_real(mes_npoint));
    assert(ped_hc_dft);
    for(unsigned ipoint=0; ipoint!=mes_npoint; ++ipoint) {
      ped_hc_dft.get()[ipoint] = ped_pdf_->value_1d(mes_x(ipoint));
    }

    // If there is no shift in the on & off pedestal position then store the
    // pedestal spectrum we calculated for later before doing FFT
    if(config_.on_off_ped_shift() == 0) {
      rebin_spectrum(off_pmf_, ped_hc_dft.get(), mes_npoint);
    }

    // Do the forward DFT transforming fped
    int plan_flags = proto_planning_enum_to_fftw_flag(config_.fftw_planning());
    uptr_fftw_plan fwd_plan = {
      fftw_plan_r2r_1d(mes_npoint, ped_hc_dft.get(), ped_hc_dft.get(),
        FFTW_R2HC, plan_flags), fftw_destroy_plan };
    assert(fwd_plan);
    fftw_execute(fwd_plan.get());
  }

  calin::simulation::pmt::PMTSimTwoPopulation pmt(config_.pmt(), /* rng = */ nullptr,
    /* use_new_stage_n_algorithm = */ false, /* adopt_rng = */ false);

  calin::ix::simulation::pmt::PMTSimPMF mes;
  mes = pmt.calc_multi_electron_spectrum(intensity_pe_, config_.intensity_rms_frac(),
    ped_hc_dft.get(), mes_npoint, 0, config_.precision(), /* log_progress= */ false,
    config_.fftw_planning());

  rebin_spectrum(mes_pmf_, mes.pn().data(), mes_npoint);

  // If pedestal is defined and on/off shift is specified then calculate off spectrum
  if(ped_pdf_ != nullptr and config_.on_off_ped_shift() != 0) {
    for(unsigned ipoint=0; ipoint!=mes_npoint; ++ipoint) {
      ped_hc_dft.get()[ipoint] = ped_pdf_->value_1d(off_x(ipoint));
    }
    rebin_spectrum(off_pmf_, ped_hc_dft.get(), mes_npoint);
  }
}

Eigen::VectorXd TwoComponentLombardMartinMES::ses_pmf() const
{
  unsigned mes_npoint = round_up_power_of_two(double(npoint_+1)/config_.sensitivity()+1);
  calin::simulation::pmt::PMTSimTwoPopulation pmt(config_.pmt(), /* rng = */ nullptr,
    /* use_new_stage_n_algorithm = */ false, /* adopt_rng = */ false);

  calin::ix::simulation::pmt::PMTSimPMF ses;
  ses = pmt.calc_pmf_fft(mes_npoint, 0, config_.precision(),
    /* log_progress= */ false, /* skip_inverse_fft = */ false,
    config_.fftw_planning());

  Eigen::VectorXd ses_pmf(npoint_);
  rebin_spectrum(ses_pmf, ses.pn().data(), mes_npoint);
  return ses_pmf;
}

Eigen::VectorXd TwoComponentLombardMartinMES::ses_pmf_full_resolution() const
{
  unsigned mes_npoint = round_up_power_of_two(double(npoint_+1)/config_.sensitivity()+1);
  calin::simulation::pmt::PMTSimTwoPopulation pmt(config_.pmt(), /* rng = */ nullptr,
    /* use_new_stage_n_algorithm = */ false, /* adopt_rng = */ false);

  calin::ix::simulation::pmt::PMTSimPMF ses;
  ses = pmt.calc_pmf_fft(mes_npoint, 0, config_.precision(),
    /* log_progress= */ false, /* skip_inverse_fft = */ false,
    config_.fftw_planning());

  Eigen::VectorXd ses_pmf(mes_npoint);
  std::copy(ses.pn().data(), ses.pn().data()+mes_npoint, ses_pmf.data());
  return ses_pmf;
}

void TwoComponentLombardMartinMES::
rebin_spectrum(Eigen::VectorXd& pmf_out, const double* mes_in, unsigned nmes) const
{
  double sensitivity_over_dx = config_.sensitivity()/config_.dx();
  double dx_over_sensitivity = config_.dx()/config_.sensitivity();

  pmf_out.setZero();
  for(unsigned imes=0; imes<nmes; imes++) {
    double xpmf_l = (double(imes)-0.5)*sensitivity_over_dx;
    double xpmf_r = (double(imes)+0.5)*sensitivity_over_dx;
    int ipmf_l = floor(xpmf_l);
    int ipmf_r = floor(xpmf_r);
    if(ipmf_l >= int(npoint_)) {
      continue;
    }
    if(ipmf_l == ipmf_r) {
      pmf_out[ipmf_l] += mes_in[imes];
    } else {
      if(ipmf_l>=0) {
        pmf_out[ipmf_l] += (ipmf_r - xpmf_l)*dx_over_sensitivity*mes_in[imes];
      }
      if(ipmf_r<npoint_) {
        pmf_out[ipmf_r] += (xpmf_r - ipmf_r)*dx_over_sensitivity*mes_in[imes];
      }
    }
#if 0
    if(imes<10)LOG(INFO) << imes << ' ' << mes_in[imes] << ' '
      << xpmf_l << ' ' << ipmf_l << ' ' << (ipmf_l>=0 ? pmf_out[ipmf_l] : 0) << ' '
      << xpmf_r << ' ' << ipmf_r << ' ' << (ipmf_r>=0 ? pmf_out[ipmf_r] : 0);
#endif
  }
}
