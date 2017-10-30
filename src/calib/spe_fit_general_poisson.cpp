/*
   calin/calib/spe_fit_general_poisson.cpp -- Stephen Fegan -- 2015-03-01

   Functions to do fit to multi-electron spectrum in the "single PE"
   domain. General Poisson MES.

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

#include <fftw3.h>

#include <math/fftw_util.hpp>
#include <math/accumulator.hpp>
#include <math/special.hpp>
#include <calib/spe_fit.hpp>
#include <calib/pmt_model_pg.hpp>
#include <util/log.hpp>

using namespace calin::math;
using namespace calin::math::fftw_util;
using namespace calin::calib::spe_fit;
using namespace calin::calib::pmt_model_pg;
using namespace calin::util::log;

using calin::math::special::SQR;
using calin::math::function::assign_parameters;

using uptr_fftw_plan = std::unique_ptr<fftw_plan_s,void(*)(fftw_plan_s*)>;
using uptr_fftw_data = std::unique_ptr<double,void(*)(void*)>;

// ============================================================================
//
// GeneralPoissonMES - Poisson model using generic SES and pedestal
// distribution. Uses FFTs to do convolution
//
// ============================================================================

GeneralPoissonMES::
GeneralPoissonMES(double x0, double dx, unsigned npoint,
                  calin::math::function::ParameterizableSingleAxisFunction* ses,
                  calin::math::function::ParameterizableSingleAxisFunction* ped,
                  calin::ix::calib::spe_fit::GeneralPoissonMESConfig config,
                  bool adopt_ses, bool adopt_ped):
    MultiElectronSpectrum(),
    ses_pdf_(ses), ped_pdf_(ped),
    adopt_ses_pdf_(adopt_ses), adopt_ped_pdf_(adopt_ped),
    config_(config), x0_(x0), dx_(dx), nsample_(npoint),
    nes_fft_(config.num_pe_convolutions())
{
  int plan_flags = proto_planning_enum_to_fftw_flag(config_.fftw_planning());

  ped_spec_ = fftw_alloc_real(nsample_);
  if(config_.include_on_off_ped_shift())
  {
    off_spec_ = fftw_alloc_real(nsample_);
    off_dfdx_ = fftw_alloc_real(nsample_);
  }
  ped_fft_ = fftw_alloc_real(nsample_);
  for(unsigned ines=0; ines<config_.num_pe_convolutions();ines++)
    nes_fft_[ines] = fftw_alloc_real(nsample_);
  mes_spec_ = fftw_alloc_real(nsample_);

  ped_plan_fwd_ =
      fftw_plan_r2r_1d(nsample_, ped_fft_, ped_fft_, FFTW_R2HC, plan_flags);
  ses_plan_fwd_ =
      fftw_plan_r2r_1d(nsample_, nes_fft_[0], nes_fft_[0], FFTW_R2HC, plan_flags);
  mes_plan_rev_ =
      fftw_plan_r2r_1d(nsample_, mes_spec_, mes_spec_, FFTW_HC2R, plan_flags);

  if(can_calculate_parameter_gradient())
  {
    unsigned ped_npar = ped_pdf_->num_parameters();

    ped_grad_.resize(ped_npar);
    ped_grad_fft_.resize(ped_npar);
    ped_grad_plan_fwd_.resize(ped_npar);
    for(unsigned ipar=0;ipar<ped_npar;ipar++)
    {
      ped_grad_[ipar] = fftw_alloc_real(nsample_);
      ped_grad_fft_[ipar] = fftw_alloc_real(nsample_);
      ped_grad_plan_fwd_[ipar] =
          fftw_plan_r2r_1d(nsample_, ped_grad_[ipar], ped_grad_fft_[ipar],
                           FFTW_R2HC, plan_flags);
    }

    if(config_.include_on_off_ped_shift())
    {
      off_grad_.resize(ped_npar);
      for(unsigned ipar=0;ipar<ped_npar;ipar++)
        off_grad_[ipar] = fftw_alloc_real(nsample_);
    }

    unsigned ses_npar = ses_pdf_->num_parameters();
    ses_grad_fft_.resize(ses_npar);
    ses_grad_plan_fwd_.resize(ses_npar);
    for(unsigned ipar=0;ipar<ses_npar;ipar++)
    {
      ses_grad_fft_[ipar] = fftw_alloc_real(nsample_);
      ses_grad_plan_fwd_[ipar] =
          fftw_plan_r2r_1d(nsample_, ses_grad_fft_[ipar], ses_grad_fft_[ipar],
                           FFTW_R2HC, plan_flags);
    }

    unsigned mes_npar = 1+ses_npar+ped_npar;
    mes_grad_.resize(mes_npar);
    mes_grad_plan_rev_.resize(mes_npar);
    for(unsigned ipar=0;ipar<mes_npar;ipar++)
    {
      mes_grad_[ipar] = fftw_alloc_real(nsample_);
      mes_grad_plan_rev_[ipar] =
          fftw_plan_r2r_1d(nsample_, mes_grad_[ipar], mes_grad_[ipar],
                           FFTW_HC2R, plan_flags);
    }
  }

  set_cache(/* force = */ true);
}

GeneralPoissonMES::~GeneralPoissonMES()
{
  for(auto x : mes_grad_plan_rev_)fftw_destroy_plan(x);
  for(auto x : ses_grad_plan_fwd_)fftw_destroy_plan(x);
  for(auto x : ped_grad_plan_fwd_)fftw_destroy_plan(x);
  for(auto x : mes_grad_)fftw_free(x);
  for(auto x : ped_grad_)fftw_free(x);
  for(auto x : off_grad_)fftw_free(x);
  for(auto x : ses_grad_fft_)fftw_free(x);
  for(auto x : ped_grad_fft_)fftw_free(x);

  fftw_destroy_plan(mes_plan_rev_);
  fftw_destroy_plan(ses_plan_fwd_);
  fftw_destroy_plan(ped_plan_fwd_);

  fftw_free(ped_fft_);
  fftw_free(ped_spec_);
  fftw_free(off_spec_);
  fftw_free(off_dfdx_);
  for(auto x : nes_fft_)fftw_free(x);
  fftw_free(mes_spec_);

  if(adopt_ses_pdf_)delete ses_pdf_;
  if(adopt_ped_pdf_)delete ped_pdf_;
}

unsigned GeneralPoissonMES::num_parameters()
{
  return num_intrinsic_parameters()
    + ses_pdf_->num_parameters() + ped_pdf_->num_parameters();
}

unsigned GeneralPoissonMES::num_intrinsic_parameters()
{
  unsigned npar = 1;
  if(config_.include_on_off_ped_shift())++npar;
  return npar;
}

auto GeneralPoissonMES::parameters() ->
    std::vector<math::function::ParameterAxis>
{
  //constexpr double tiny_val = std::numeric_limits<double>::min();
  constexpr double inf = std::numeric_limits<double>::infinity();
  std::vector<math::function::ParameterAxis> pvec;
  pvec.push_back({ "light_intensity", "PE", true, 0, false, inf });
  if(config_.include_on_off_ped_shift())
    pvec.push_back({ "off_ped_shift", "DC", false, -inf, false, inf });
  std::vector<math::function::ParameterAxis> pped { ped_pdf_->parameters() };
  for(auto& ip : pped)ip.name = std::string("ped.") + ip.name;
  pvec.insert(pvec.end(), pped.begin(), pped.end());
  std::vector<math::function::ParameterAxis> pses { ses_pdf_->parameters() };
  for(auto& ip : pses)ip.name = std::string("ses.") + ip.name;
  pvec.insert(pvec.end(), pses.begin(), pses.end());
  return pvec;
}

Eigen::VectorXd GeneralPoissonMES::parameter_values()
{
  Eigen::VectorXd param(num_parameters());
  param[iparam_light_intensity()] = intensity_pe_;
  if(config_.include_on_off_ped_shift())
    param[iparam_off_ped_shift()] = off_ped_shift_dc_;
  unsigned ip = iparam_ped();
  unsigned num_ped_params = ped_pdf_->num_parameters();
  param.segment(ip,num_ped_params) = ped_pdf_->parameter_values();
  ip += num_ped_params;
  param.segment(ip,ses_pdf_->num_parameters()) = ses_pdf_->parameter_values();
  return param;
}

void GeneralPoissonMES::set_parameter_values(ConstVecRef values)
{
  verify_set_parameter_values(values, "GeneralPoissonMES");
  assign_parameters(values.data()+iparam_light_intensity(), intensity_pe_);
  if(config_.include_on_off_ped_shift())
    assign_parameters(values.data()+iparam_off_ped_shift(), off_ped_shift_dc_);
  unsigned ip = iparam_ped();
  unsigned num_ped_params = ped_pdf_->num_parameters();
  ped_pdf_->set_parameter_values(values.segment(ip,num_ped_params));
  ip += num_ped_params;
  ses_pdf_->set_parameter_values(values.segment(ip,ses_pdf_->num_parameters()));
  set_cache();
}

bool GeneralPoissonMES::can_calculate_parameter_gradient()
{
  //return false;
  return ped_pdf_->can_calculate_parameter_gradient() and
    ses_pdf_->can_calculate_parameter_gradient() and
    (not config_.include_on_off_ped_shift() or
      ped_pdf_->can_calculate_gradient());
}

bool GeneralPoissonMES::can_calculate_parameter_hessian()
{
  return false; // for the moment we are lazy
}

double GeneralPoissonMES::pdf_ped(double x)
{
  if(config_.include_on_off_ped_shift())
    return std::max(off_spec_[ibin(x)],0.0);
  else
    return std::max(ped_spec_[ibin(x)],0.0);
}

double GeneralPoissonMES::pdf_gradient_ped(double x, VecRef gradient)
{
  assert(can_calculate_parameter_gradient());
  int thebin = ibin(x);
  //std::cout << "HELLO : " << x << ' ' << thebin << '\n';
  unsigned npar = ped_pdf_->num_parameters();
  gradient.resize(num_parameters());
  gradient.setZero();
  if(config_.include_on_off_ped_shift()) {
    gradient[iparam_off_ped_shift()] = off_dfdx_[thebin];
    unsigned ipedpar = iparam_ped();
    for(unsigned ipar=0;ipar<npar;ipar++)
      gradient[ipedpar + ipar] = off_grad_[ipar][thebin];
    return std::max(off_spec_[ibin(x)],0.0);
  } else {
    unsigned ipedpar = iparam_ped();
    for(unsigned ipar=0;ipar<npar;ipar++)
      gradient[ipedpar + ipar] = ped_grad_[ipar][thebin];
    return std::max(ped_spec_[thebin],0.0);
  }
}

double GeneralPoissonMES::pdf_gradient_hessian_ped(double x, VecRef gradient,
                                                   MatRef hessian)
{
  throw std::logic_error("GeneralPoissonMES: cannot calculate ped hessian");
}

double GeneralPoissonMES::pdf_mes(double x)
{
  //std::cout << x << ' ' << ibin(x) << ' ' << mes_spec_[ibin(x)] << '\n';
  return std::max(mes_spec_[ibin(x)],0.0);
}

double GeneralPoissonMES::pdf_gradient_mes(double x, VecRef gradient)
{
  assert(can_calculate_parameter_gradient());
  int thebin = ibin(x);
  unsigned npar = num_parameters();
  unsigned ipedpar = iparam_ped();
  gradient.resize(npar);
  gradient.setZero();
  gradient[iparam_light_intensity()] = mes_grad_[0][thebin];
  for(unsigned ipar=ipedpar;ipar<npar;ipar++)
    gradient[ipar] = mes_grad_[ipar-ipedpar+1][thebin];
  return std::max(mes_spec_[thebin],0.0);
}

double GeneralPoissonMES::pdf_gradient_hessian_mes(double x, VecRef gradient,
                                                   MatRef hessian)
{
  throw std::logic_error("GeneralPoissonMES: cannot calculate ses hessian");
}

double GeneralPoissonMES::ped_rms_dc()
{
  assert(false);
  return 0;
}

double GeneralPoissonMES::ped_zero_dc()
{
  assert(false);
  return 0;
}

double GeneralPoissonMES::ses_mean_dc()
{
  assert(false);
  return 0;
}

double GeneralPoissonMES::ses_rms_pe()
{
  assert(false);
  return 0;
}

Eigen::VectorXd GeneralPoissonMES::multi_electron_spectrum() const
{
  Eigen::VectorXd spec(nsample_);
  std::copy(mes_spec_, mes_spec_+nsample_, spec.data());
  return spec;
}

Eigen::VectorXd GeneralPoissonMES::pedestal_spectrum() const
{
  Eigen::VectorXd spec(nsample_);
  std::copy(ped_spec_, ped_spec_+nsample_, spec.data());
  return spec;
}

Eigen::VectorXd GeneralPoissonMES::off_pedestal_spectrum() const
{
  if(not config_.include_on_off_ped_shift())return pedestal_spectrum();
  Eigen::VectorXd spec(nsample_);
  std::copy(off_spec_, off_spec_+nsample_, spec.data());
  return spec;
}

Eigen::VectorXd GeneralPoissonMES::
multi_electron_spectrum_gradient(unsigned iparam) const
{
  if(iparam >= const_cast<GeneralPoissonMES*>(this)->num_parameters())
    throw std::out_of_range(
        "GeneralPoissonMES::multi_electron_spectrum_gradient: "
        "iparam out of range");

  unsigned ipedpar = iparam_ped();
  Eigen::VectorXd grad(nsample_);
  if(iparam == iparam_light_intensity())
    std::copy(mes_grad_[0], mes_grad_[0]+nsample_, grad.data());
  else if(iparam>=ipedpar)
    std::copy(mes_grad_[iparam-ipedpar+1], mes_grad_[iparam-ipedpar+1]+nsample_, grad.data());
  else
    grad.setZero();
  return grad;
}

Eigen::VectorXd GeneralPoissonMES::
pedestal_spectrum_gradient(unsigned iparam) const
{
  if(iparam >= const_cast<GeneralPoissonMES*>(this)->num_parameters())
    throw std::out_of_range("GeneralPoissonMES::pedestal_spectrum_gradient: "
                            "iparam out of range");

  Eigen::VectorXd grad(nsample_);
  if(iparam==0 || iparam>ped_pdf_->num_parameters())
    grad.setZero();
  else
    std::copy(ped_grad_[iparam-1], ped_grad_[iparam-1]+nsample_, grad.data());
  return grad;
}

Eigen::VectorXd GeneralPoissonMES::
off_pedestal_spectrum_gradient(unsigned iparam) const
{
  if(not config_.include_on_off_ped_shift())
    return off_pedestal_spectrum_gradient(iparam);

  if(iparam >= const_cast<GeneralPoissonMES*>(this)->num_parameters())
    throw std::out_of_range("GeneralPoissonMES::off_pedestal_spectrum_gradient: "
                            "iparam out of range");

  Eigen::VectorXd grad(nsample_);
  unsigned ipedpar = iparam_ped();
  if(iparam == iparam_off_ped_shift())
    std::copy(off_dfdx_, off_dfdx_+nsample_, grad.data());
  else if(iparam >= ipedpar and iparam < iparam_ses())
    std::copy(off_grad_[iparam-ipedpar], off_grad_[iparam-ipedpar]+nsample_, grad.data());
  else
    grad.setZero();
  return grad;
}

Eigen::VectorXd GeneralPoissonMES::
single_electron_spectrum_gradient(unsigned iparam) const
{
  int plan_flags = proto_planning_enum_to_fftw_flag(config_.fftw_planning());

  if(iparam >= const_cast<GeneralPoissonMES*>(this)->num_parameters())
    throw std::out_of_range("GeneralPoissonMES::pedestal_spectrum_gradient: "
                            "iparam out of range");

  unsigned isespar = iparam_ses();
  if(iparam<=isespar)return Eigen::VectorXd::Zero(nsample_);
  iparam -= isespar;

  uptr_fftw_data spec_buffer { fftw_alloc_real(nsample_), fftw_free };
  assert(spec_buffer);

  uptr_fftw_plan spec_plan = {
    fftw_plan_r2r_1d(nsample_, spec_buffer.get(), spec_buffer.get(),
                     FFTW_HC2R, plan_flags), fftw_destroy_plan };
  assert(spec_plan);

  std::copy(ses_grad_fft_[iparam], ses_grad_fft_[iparam]+nsample_,
            spec_buffer.get());
  fftw_execute(spec_plan.get());

  Eigen::VectorXd spec_gradient(nsample_);
  double norm { 1.0/double(nsample_) };
  std::transform(spec_buffer.get(), spec_buffer.get()+nsample_,
                 spec_gradient.data(),
                 [norm](double x){return x*norm;});
  return spec_gradient;
}

Eigen::VectorXd GeneralPoissonMES::n_electron_spectrum(unsigned n) const
{
  int plan_flags = proto_planning_enum_to_fftw_flag(config_.fftw_planning());

  if(n==0 or n>config_.num_pe_convolutions())
    throw std::out_of_range("GeneralPoissonMES::n_electron_spectrum: "
                            "number of PEs out of range");

  uptr_fftw_data spec_buffer { fftw_alloc_real(nsample_), fftw_free };
  assert(spec_buffer);

  uptr_fftw_plan spec_plan = {
    fftw_plan_r2r_1d(nsample_, spec_buffer.get(), spec_buffer.get(),
                     FFTW_HC2R, plan_flags), fftw_destroy_plan };
  assert(spec_plan);

  std::copy(nes_fft_[n-1], nes_fft_[n-1]+nsample_, spec_buffer.get());
  fftw_execute(spec_plan.get());

  Eigen::VectorXd spec(nsample_);
  double norm { 1.0/double(nsample_) };
  std::transform(spec_buffer.get(), spec_buffer.get()+nsample_, spec.data(),
                 [norm](double x){return x*norm;});
  return spec;
}

Eigen::VectorXd GeneralPoissonMES::n_electron_spectrum_with_pedestal(unsigned n) const
{
  int plan_flags = proto_planning_enum_to_fftw_flag(config_.fftw_planning());

  if(n==0)return pedestal_spectrum();

  if(n>config_.num_pe_convolutions())
    throw std::out_of_range("GeneralPoissonMES::n_electron_spectrum_with_pedestal: "
                            "number of PEs out of range");

  uptr_fftw_data spec_buffer { fftw_alloc_real(nsample_), fftw_free };
  assert(spec_buffer);

  uptr_fftw_plan spec_plan = {
    fftw_plan_r2r_1d(nsample_, spec_buffer.get(), spec_buffer.get(),
                     FFTW_HC2R, plan_flags), fftw_destroy_plan };
  assert(spec_plan);

  hcvec_scale_and_multiply(spec_buffer.get(), nes_fft_[n-1], ped_fft_,
                           nsample_, 1.0/double(nsample_));

  fftw_execute(spec_plan.get());

  Eigen::VectorXd spec(nsample_);
  std::copy(spec_buffer.get(), spec_buffer.get()+nsample_, spec.data());
  return spec;
}

Eigen::VectorXd GeneralPoissonMES::mes_n_electron_cpt(unsigned n) const
{
  int plan_flags = proto_planning_enum_to_fftw_flag(config_.fftw_planning());

  if(n>config_.num_pe_convolutions())
    throw std::out_of_range("GeneralPoissonMES::mes_n_electron_cpt: "
                            "number of PEs out of range");

  Eigen::VectorXd spec(nsample_);
  double log_nsample = std::log(double(nsample_));

  if(n==0)
  {
    double scale = std::exp(-intensity_pe_);
    std::transform(ped_spec_, ped_spec_+nsample_, spec.data(),
                   [scale](const double& x){ return x*scale; });
    return spec;
  }

  uptr_fftw_data spec_buffer { fftw_alloc_real(nsample_), fftw_free };
  assert(spec_buffer);

  uptr_fftw_plan spec_plan = {
    fftw_plan_r2r_1d(nsample_, spec_buffer.get(), spec_buffer.get(),
                     FFTW_HC2R, plan_flags), fftw_destroy_plan };
  assert(spec_plan);

  double log_intensity = std::log(intensity_pe_);
  double dbl_n { double(n) };
  double poisson_factor {
    std::exp(dbl_n*log_intensity - intensity_pe_ - lgamma(dbl_n+1.0) -
             log_nsample) };

  hcvec_scale_and_multiply(spec_buffer.get(), nes_fft_[n-1], ped_fft_,
                           nsample_, poisson_factor*dx_);

  fftw_execute(spec_plan.get());
  std::copy(spec_buffer.get(), spec_buffer.get()+nsample_, spec.data());
  return spec;
}

Eigen::VectorXd GeneralPoissonMES::
extract_ped_gradient_values(ConstVecRef gradient)
{
  return gradient.segment(iparam_ped(), ped_pdf_->num_parameters());
}

Eigen::VectorXd GeneralPoissonMES::
extract_ses_gradient_values(ConstVecRef gradient)
{
  return gradient.segment(iparam_ses(), ses_pdf_->num_parameters());
}

Eigen::MatrixXd GeneralPoissonMES::
extract_ped_hessian_values(ConstMatRef hessian)
{
  unsigned ped_npar = ped_pdf_->num_parameters();
  unsigned ipedpar = iparam_ped();
  return hessian.block(ipedpar, ipedpar, ped_npar, ped_npar);
}

Eigen::MatrixXd GeneralPoissonMES::
extract_ses_hessian_values(ConstMatRef hessian)
{
  unsigned ses_npar = ses_pdf_->num_parameters();
  unsigned isespar = iparam_ses();
  return hessian.block(isespar, isespar, ses_npar, ses_npar);
}

int GeneralPoissonMES::ibin(double x) const
{
  int thebin = std::floor((x-x0_)/dx_);
  if(thebin<0 or thebin>=(int)nsample_)
  {
    std::ostringstream str;
    str << "GeneralPoissonMES::ibin: x=" << x
        << " out of range (bin=" << thebin << ", nsample=" << nsample_ << ")";
    throw std::out_of_range(str.str());
  }
  return thebin;
}

void GeneralPoissonMES::set_cache(bool force)
{
  // THIS FUNCTION IS TOO LONG

  bool calc_gradient = can_calculate_parameter_gradient();
  unsigned ses_npar = ses_pdf_->num_parameters();

  // If the SES has no paramerters then we assume it doesn't change, so no
  // need to do the FFTs each time
  if(ses_npar or force)
  {
    Eigen::VectorXd ses_gradient(ses_npar);
    ses_gradient.setZero();
    calin::math::accumulator::LikelihoodAccumulator ses_acc;
    for(unsigned isample = 0;isample<nsample_;isample++)
    {
      // The SES describes charge delivered by the PMT and is assumed be
      // positive. Hence the function is sampled from zero here and not
      // from x0. We could allow negative charges by using the domain
      // axis limits to specify the minimum bound. We would need to
      // adjust x0 to compensate for where the pedestal distribution
      // would end up.
      const double x = ses_x(isample); // (0.5+double(isample))*dx_;
      double val;
      if(calc_gradient)
      {
        val = ses_pdf_->value_and_parameter_gradient_1d(x, ses_gradient);
        if(!std::isfinite(val)){ val = 0; ses_gradient.setZero(); }
        for(unsigned ipar=0;ipar<ses_npar;ipar++)
          ses_grad_fft_[ipar][isample] = ses_gradient[ipar];
      }
      else
      {
        val = ses_pdf_->value_1d(x);
        if(!std::isfinite(val))val = 0;
      }
      nes_fft_[0][isample] = val;
      ses_acc.accumulate(val);
    }

    if(config_.ses_norm_warning_threshold()>0 and
      (config_.max_ses_norm_warning()==0 or
        n_ses_norm_warning_<config_.max_ses_norm_warning()) and
      std::abs(ses_acc.total() * dx_ - 1.0) >
        config_.ses_norm_warning_threshold()/double(nsample_))
    {
      // The owls are not what they seem
      LOG(WARNING) << "SES normalization is significantly different from 1.0: "
                   << ses_acc.total() * dx_ << '\n'
                   << "SES parameter values : "
                   << ses_pdf_->parameter_values().transpose();
      n_ses_norm_warning_++;
      if(n_ses_norm_warning_ == config_.ses_norm_warning_threshold())
        LOG(INFO) << "Further SES normalization warnings will be suppressed!";
    }

    fftw_execute(ses_plan_fwd_);
    for(unsigned ines=1;ines<config_.num_pe_convolutions();ines++)
      hcvec_scale_and_multiply(nes_fft_[ines], nes_fft_[ines-1], nes_fft_[0],
                               nsample_, dx_);
  } // if(ses_npar or force)

  // Same for PED, if it has no paramerters then we assume it doesn't change,
  // so, again, no need to do the FFTs each time
  unsigned ped_npar = ped_pdf_->num_parameters();
  if(ped_npar or force)
  {
    Eigen::VectorXd ped_gradient(ped_npar);
    ped_gradient.setZero();
    for(unsigned isample = 0;isample<nsample_;isample++)
    {
      const double x = ped_x(isample); // x0_ + (0.5+double(isample))*dx_;
      double val;
      if(calc_gradient)
      {
        val = ped_pdf_->value_and_parameter_gradient_1d(x, ped_gradient);
        if(!std::isfinite(val)){ val = 0; ped_gradient.setZero(); }
        for(unsigned ipar=0;ipar<ped_npar;ipar++)
          ped_grad_[ipar][isample] = ped_gradient[ipar];
      }
      else
      {
        val = ped_pdf_->value_1d(x);
        if(!std::isfinite(val))val = 0;
      }
      ped_spec_[isample] = ped_fft_[isample] = val;
    }
    fftw_execute(ped_plan_fwd_);
  } // if(ped_npar or force)

  if(config_.include_on_off_ped_shift())
  {
    Eigen::VectorXd ped_gradient(ped_npar);
    ped_gradient.setZero();
    for(unsigned isample = 0;isample<nsample_;isample++)
    {
      const double x = ped_x(isample) + off_ped_shift_dc_;
      double val;
      if(calc_gradient)
      {
        double dfdx;
        val = ped_pdf_->value_and_parameter_gradient_1d(x, ped_gradient);
        ped_pdf_->value_and_gradient_1d(x, dfdx);
        if(!std::isfinite(val)){ val = dfdx = 0; ped_gradient.setZero(); }
        for(unsigned ipar=0;ipar<ped_npar;ipar++)
          off_grad_[ipar][isample] = ped_gradient[ipar];
        off_dfdx_[isample] = dfdx;
      }
      else
      {
        val = ped_pdf_->value_1d(x);
        if(!std::isfinite(val))val = 0;
      }
      off_spec_[isample] = val;
    }
  }

  unsigned mes_npar = 1+ses_npar+ped_npar;
  if(calc_gradient)
  {
    for(unsigned ipar=0;ipar<ses_npar;ipar++)
      fftw_execute(ses_grad_plan_fwd_[ipar]);
    for(unsigned ipar=0;ipar<ped_npar;ipar++)
      fftw_execute(ped_grad_plan_fwd_[ipar]);
    for(unsigned ipar=0;ipar<mes_npar;ipar++)
      std::fill(mes_grad_[ipar], mes_grad_[ipar]+nsample_, 0.0);
  }

  std::fill(mes_spec_, mes_spec_+nsample_, 0.0);

  double log_intensity = std::log(intensity_pe_);
  double log_nsample = std::log(double(nsample_));
  for(unsigned ines = 0;ines<config_.num_pe_convolutions();ines++)
  {
    double dbl_n { double(ines+1) };
    double poisson_factor {
      std::exp(dbl_n*log_intensity - intensity_pe_ - lgamma(dbl_n+1.0) -
               log_nsample) };
    hcvec_scale_and_add(mes_spec_, nes_fft_[ines], nsample_, poisson_factor);

    if(calc_gradient)
    {
      hcvec_scale_and_add(mes_grad_[0],
        nes_fft_[ines], nsample_, poisson_factor*(dbl_n/intensity_pe_ - 1.0));
      if(ines>0)
        for(unsigned ipar=0;ipar<ses_npar;ipar++)
          hcvec_scale_and_add(mes_grad_[1+ped_npar+ipar], nes_fft_[ines-1],
                              nsample_, dbl_n*poisson_factor);
    }
  }

  if(calc_gradient)
  {
    hcvec_scale_and_multiply(mes_grad_[0],
      mes_grad_[0], ped_fft_, nsample_, dx_);
    hcvec_scale_and_add(mes_grad_[0],
      ped_fft_, nsample_, -std::exp(-intensity_pe_-log_nsample));
    for(unsigned ipar=0;ipar<ped_npar;ipar++)
    {
      hcvec_scale_and_multiply(mes_grad_[1+ipar],
        mes_spec_, ped_grad_fft_[ipar], nsample_, dx_);
      hcvec_scale_and_add(mes_grad_[1+ipar],
        ped_grad_fft_[ipar], nsample_, std::exp(-intensity_pe_-log_nsample));
    }
    for(unsigned ipar=0;ipar<ses_npar;ipar++)
    {
      hcvec_scale_and_multiply(mes_grad_[1+ped_npar+ipar],
        mes_grad_[1+ped_npar+ipar], ses_grad_fft_[ipar], nsample_, dx_);
      hcvec_scale_and_add(mes_grad_[1+ped_npar+ipar],
        ses_grad_fft_[ipar], nsample_,
        std::exp(log_intensity -intensity_pe_ - log_nsample));
      hcvec_scale_and_multiply(mes_grad_[1+ped_npar+ipar],
        mes_grad_[1+ped_npar+ipar], ped_fft_, nsample_, dx_);
    }
  }

  hcvec_scale_and_multiply(mes_spec_, mes_spec_, ped_fft_, nsample_, dx_);
  hcvec_scale_and_add(mes_spec_, ped_fft_, nsample_,
                      std::exp(-intensity_pe_-log_nsample));
  fftw_execute(mes_plan_rev_);

  if(calc_gradient)
    for(unsigned ipar=0;ipar<mes_npar;ipar++)
      fftw_execute(mes_grad_plan_rev_[ipar]);
}

calin::ix::calib::spe_fit::GeneralPoissonMESConfig
GeneralPoissonMES::default_config()
{
  calin::ix::calib::spe_fit::GeneralPoissonMESConfig config;
  config.set_num_pe_convolutions(10);
  config.set_max_ses_norm_warning(10);
  config.set_ses_norm_warning_threshold(1.0);
  //config.set_fftw_planning(calin::ix::math::fftw_util::PATIENT);
  return config;
}
