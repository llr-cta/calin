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
#include <calib/pmt_ses_models.hpp>
#include <util/log.hpp>
#include <math/special.hpp>

using namespace calin::math::fftw_util;
using namespace calin::calib::spe_fit;
using calin::math::special::round_up_power_of_two;
using namespace calin::util::log;

LombardMartinPrescottMES::
LombardMartinPrescottMES(double x0, unsigned npoint,
    calin::math::function::ParameterizableSingleAxisFunction* ped_pdf,
    const calin::ix::calib::spe_fit::LombardMartinPrescottMESConfig& config,
    bool adopt_ped_pdf):
  MultiElectronSpectrum(), config_(config), x0_(x0), dx_inv_(1.0/config.dx()), npoint_(npoint),
  ped_pdf_(ped_pdf), adopt_ped_pdf_(adopt_ped_pdf), mes_pmf_(npoint), off_pmf_(npoint)
{
  calculate_mes();
}

LombardMartinPrescottMES::~LombardMartinPrescottMES()
{
  if(adopt_ped_pdf_)delete ped_pdf_;
}

unsigned LombardMartinPrescottMES::num_parameters()
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
LombardMartinPrescottMES::parameters()
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

Eigen::VectorXd LombardMartinPrescottMES::parameter_values()
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

void LombardMartinPrescottMES::set_parameter_values(ConstVecRef values)
{
  verify_set_parameter_values(values, "LombardMartinPrescottMES");
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

bool LombardMartinPrescottMES::can_calculate_parameter_gradient()
{
  return false;
}

bool LombardMartinPrescottMES::can_calculate_parameter_hessian()
{
  return false;
}

double LombardMartinPrescottMES::pdf_ped(double x)
{
  if(ped_pdf_ == nullptr) {
    throw std::runtime_error("LombardMartinPrescottMES::pdf_ped: no pedestal pdf supplied");
  }
  int i = ibin(x);
  if(i<0 or i>=npoint_)return 0.0;
  return std::max(off_pmf_[i], 0.0)*dx_inv_;
}

double LombardMartinPrescottMES::pdf_gradient_ped(double x, VecRef gradient)
{
  throw std::runtime_error("LombardMartinPrescottMES: pdf_gradient_ped not implemented");
  return 0;
}

double LombardMartinPrescottMES::pdf_gradient_hessian_ped(double x, VecRef gradient, MatRef hessian)
{
  throw std::runtime_error("LombardMartinPrescottMES: pdf_gradient_hessian_ped not implemented");
  return 0;
}

double LombardMartinPrescottMES::pdf_mes(double x)
{
  int i = ibin(x);
  if(i<0 or i>=npoint_)return 0.0;
  return std::max(mes_pmf_[i], 0.0)*dx_inv_;
}

double LombardMartinPrescottMES::pdf_gradient_mes(double x, VecRef gradient)
{
  throw std::runtime_error("LombardMartinPrescottMES: pdf_gradient_mes not implemented");
  return 0;
}

double LombardMartinPrescottMES::pdf_gradient_hessian_mes(double x, VecRef gradient, MatRef hessian)
{
  throw std::runtime_error("LombardMartinPrescottMES: pdf_gradient_hessian_mes not implemented");
  return 0;
}

double LombardMartinPrescottMES::intensity_pe()
{
  return intensity_pe_;
}

double LombardMartinPrescottMES::ped_rms_dc()
{
  double m = ped_sum_px_/ped_sum_p_;
  return std::sqrt(ped_sum_pxx_/ped_sum_p_ - m*m);
}

double LombardMartinPrescottMES::ped_zero_dc()
{
  return ped_sum_px_/ped_sum_p_;
}

double LombardMartinPrescottMES::ses_mean_dc()
{
  return pmt_total_gain_*config_.sensitivity();
}

double LombardMartinPrescottMES::ses_rms_pe()
{
  return pmt_resolution_;
}

calin::ix::calib::spe_fit::LombardMartinPrescottMESConfig
LombardMartinPrescottMES::default_config()
{
  calin::ix::calib::spe_fit::LombardMartinPrescottMESConfig config;

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

void LombardMartinPrescottMES::calculate_mes()
{
  unsigned mes_npoint = round_up_power_of_two(double(npoint_+1)/config_.sensitivity()+1);

  uptr_fftw_data ped { nullptr, fftw_free };
  if(ped_pdf_ != nullptr) {
    ped_sum_p_ = 0;
    ped_sum_px_ = 0;
    ped_sum_pxx_ = 0;

    // FFT of pedestal
    ped.reset(fftw_alloc_real(mes_npoint));
    assert(ped);
    for(unsigned ipoint=0; ipoint!=mes_npoint; ++ipoint) {
      double x = mes_x(ipoint);
      double p = ped_pdf_->value_1d(x)*config_.sensitivity();
      ped.get()[ipoint] = p;
      ped_sum_p_ += p;
      ped_sum_px_ += p*x;
      ped_sum_pxx_ += p*x*x;
    }

    // If there is no shift in the on & off pedestal position then store the
    // pedestal spectrum we calculated for later before doing FFT
    if(config_.on_off_ped_shift() == 0) {
      rebin_spectrum(off_pmf_, ped.get(), mes_npoint);
    }
  }

  calin::calib::pmt_ses_models::LombardMartinPrescottPMTModel pmt(config_.pmt(),
    config_.precision(), config_.fftw_planning());

  pmt_total_gain_ = pmt.total_gain();
  pmt_resolution_ = pmt.resolution();

  Eigen::VectorXd mes;
  mes = pmt.calc_mes(intensity_pe_, config_.intensity_rms_frac(), mes_npoint, ped.get());

  rebin_spectrum(mes_pmf_, mes.data(), mes_npoint);

  // If pedestal is defined and on/off shift is specified then calculate off spectrum
  if(ped_pdf_ != nullptr and config_.on_off_ped_shift() != 0) {
    for(unsigned ipoint=0; ipoint!=mes_npoint; ++ipoint) {
      ped.get()[ipoint] = ped_pdf_->value_1d(off_x(ipoint))*config_.sensitivity();
    }
    rebin_spectrum(off_pmf_, ped.get(), mes_npoint);
  }
}

Eigen::VectorXd LombardMartinPrescottMES::ses_pmf() const
{
  unsigned mes_npoint = round_up_power_of_two(double(npoint_+1)/config_.sensitivity()+1);

  calin::calib::pmt_ses_models::LombardMartinPrescottPMTModel pmt(config_.pmt(),
    config_.precision(), config_.fftw_planning());

  Eigen::VectorXd ses;
  ses = pmt.calc_ses(mes_npoint);

  Eigen::VectorXd ses_pmf(npoint_);
  rebin_spectrum(ses_pmf, ses.data(), mes_npoint);
  return ses_pmf;
}

Eigen::VectorXd LombardMartinPrescottMES::ses_pmf_full_resolution() const
{
  unsigned mes_npoint = round_up_power_of_two(double(npoint_+1)/config_.sensitivity()+1);

  calin::calib::pmt_ses_models::LombardMartinPrescottPMTModel pmt(config_.pmt(),
    config_.precision(), config_.fftw_planning());

  Eigen::VectorXd ses;
  ses = pmt.calc_ses(mes_npoint);

  return ses;
}

void LombardMartinPrescottMES::
rebin_spectrum(Eigen::VectorXd& pmf_out, const double* mes_in, unsigned nmes) const
{
  double sensitivity_over_dx = config_.sensitivity()*dx_inv_;
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
