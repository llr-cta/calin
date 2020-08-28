/*
   calin/calib/pmt_ses_models.cpp -- Stephen Fegan -- 2017-04-24

   PMT single-electron spectrum models

   Copyright 2017, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <limits>
#include <cmath>
#include <algorithm>
#include <fftw3.h>

#include <Eigen/Dense>

#include <calib/pmt_ses_models.hpp>
#include <util/log.hpp>
#include <math/fftw_util.hpp>
#include <math/special.hpp>
#include <math/accumulator.hpp>
#include <math/least_squares.hpp>

using namespace calin::util::log;
using namespace calin::calib::pmt_ses_models;
using namespace calin::math::fftw_util;
using calin::math::special::SQR;

namespace {
  void validate_config(const calin::ix::calib::pmt_ses_models::LombardMartinPrescottPMTModelConfig& config) {
    if(config.num_stage()<2)
      throw std::runtime_error("LombardMartinPrescottPMTModel: number of stages must be greater than 1.");
    if(config.total_gain()<=0)
      throw std::runtime_error("LombardMartinPrescottPMTModel: total gain must be positive.");
    if(config.stage_0_lo_prob()<0 or config.stage_0_lo_prob()>1)
      throw std::runtime_error("LombardMartinPrescottPMTModel: stage 0 low-charge probability must be in the range [0,1].");
    if(config.stage_0_lo_prob()<1) {
      if(config.stage_0_hi_gain()<=0)
        throw std::runtime_error("LombardMartinPrescottPMTModel: stage 0 high-charge gain must be positive.");
      if(config.stage_0_hi_gain_rms_frac()<0)
        throw std::runtime_error("LombardMartinPrescottPMTModel: stage 0 high-charge excess RMS must be zero or positive.");
    }
    if(config.stage_0_lo_prob()>0) {
      if(config.stage_0_lo_gain()<=0)
        throw std::runtime_error("LombardMartinPrescottPMTModel: stage 0 low-charge gain must be positive.");
      if(config.stage_0_lo_gain_rms_frac()<0)
        throw std::runtime_error("LombardMartinPrescottPMTModel: stage 0 low-charge excess RMS must be zero or positive.");
    }
    if(config.stage_n_gain_rms_frac()<0)
      throw std::runtime_error("PMTSimTwoPopulation: stage 1+ excess RMS must be zero or positive.");
  }

  inline double polya_Exx(double mean, double rms_frac) {
    return mean * (1.0 + mean * (1.0 + SQR(rms_frac)));
  }
} // anonymous namespace

LombardMartinPrescottPMTModel::
LombardMartinPrescottPMTModel(
  const calin::ix::calib::pmt_ses_models::LombardMartinPrescottPMTModelConfig& config,
  double precision, calin::ix::math::fftw_util::FFTWPlanningRigor fftw_rigor):
    config_(config), precision_(precision), fftw_rigor_(fftw_rigor)
{
  validate_config(config);

  if(config.stage_0_lo_prob()<1) {
    double Phi = (1 - config.stage_0_lo_prob());
    stage_0_gain_ = Phi * config.stage_0_hi_gain();
    stage_0_Exx_ = Phi * polya_Exx(config.stage_0_hi_gain(), config.stage_0_hi_gain_rms_frac());
    stage_0_pmf_ = polya_pmf(config.stage_0_hi_gain(), config.stage_0_hi_gain_rms_frac(), precision);
    if(Phi < 1.0) {
      std::transform(stage_0_pmf_.begin(), stage_0_pmf_.end(), stage_0_pmf_.begin(),
        [Phi](double x) { return Phi*x; });
    }
  }

  if(config.stage_0_lo_prob()>0) {
    double Plo = config.stage_0_lo_prob();
    if(config.stage_0_lo_pmf_size() == 0) {
      stage_0_gain_ += Plo * config.stage_0_lo_gain();
      stage_0_Exx_ += Plo * polya_Exx(config.stage_0_lo_gain(), config.stage_0_lo_gain_rms_frac());
      std::vector<double> pmf = polya_pmf(config.stage_0_lo_gain(), config.stage_0_lo_gain_rms_frac(), precision);
      if(pmf.size() > stage_0_pmf_.size()) {
        stage_0_pmf_.resize(pmf.size(), 0.0);
      }
      std::transform(pmf.begin(), pmf.end(), stage_0_pmf_.begin(), stage_0_pmf_.begin(),
        [Plo](double x, double y) { return Plo*x + y; });
    } else {
      if(config.stage_0_lo_pmf_size() > stage_0_pmf_.size()) {
        stage_0_pmf_.resize(config.stage_0_lo_pmf_size(), 0.0);
      }
      calin::math::accumulator::KahanAccumulator psum;
      for(int i=0; i<config.stage_0_lo_pmf_size(); ++i) {
        double pi = config.stage_0_lo_pmf(i);
        psum.accumulate(pi);
        stage_0_pmf_[i] += Plo * pi;
        stage_0_gain_ += Plo * double(i) * pi;
        stage_0_Exx_ += Plo * SQR(double(i)) * pi;
      }
      psum.accumulate(-1.0);
      if(std::abs(psum.total()) > precision) {
        LOG(WARNING) << "LombardMartinPrescottPMTModel: stage 0 low-gain PMF normalization error exceeds specified precision.";
      }
    }
  }

  // Compute the stage-n gain required to achieve desired total gain, iterating
  // for zero-suppression if required
  set_stage_n_gain(
    std::pow(config.total_gain()/stage_0_gain_, 1/double(config_.num_stage()-1)));
  unsigned iter = 10;
  while(std::abs(total_gain_/config.total_gain() - 1)>precision) {
    if(iter == 0) {
      throw std::runtime_error("LombardMartinPrescottPMTModel: iteration limit reached in calculating stage N gain");
    }
    --iter;
    set_stage_n_gain(stage_n_gain_ *
      std::pow(config.total_gain()/total_gain_, 1/double(config_.num_stage()-1)));
  }

  // Compute "zero suppression adjusted" stage 0 PMF
  if(config_.suppress_zero()) {
    double scale = 1.0 / (1.0-p0_);
    double shift = -p0_*scale;
    stage_0_pmf_zsa_.resize(stage_0_pmf_.size());
    std::transform(stage_0_pmf_.begin(), stage_0_pmf_.end(), stage_0_pmf_zsa_.begin(),
      [scale](double x) { return x*scale; });
    stage_0_pmf_zsa_[0] += shift;
  } else {
    stage_0_pmf_zsa_ = stage_0_pmf_;
  }
}

void LombardMartinPrescottPMTModel::set_stage_n_gain(double stage_n_gain)
{
  stage_n_gain_ = stage_n_gain;
  stage_n_pmf_ = polya_pmf(stage_n_gain_, config_.stage_n_gain_rms_frac(), precision_);
  stage_n_Exx_ = polya_Exx(stage_n_gain_, config_.stage_n_gain_rms_frac());

  total_gain_ = stage_n_gain_;
  p0_ = stage_n_pmf_[0];
  double Exx = stage_n_Exx_;

  for(int istage=0; istage<config_.num_stage()-2; ++istage) {
    Exx = (stage_n_Exx_-stage_n_gain_)*SQR(total_gain_) + Exx*stage_n_gain_;
    total_gain_ *= stage_n_gain_;
    p0_ = calin::math::least_squares::polyval(stage_n_pmf_, p0_);
  }

  Exx = (stage_0_Exx_-stage_0_gain_)*SQR(total_gain_) + Exx*stage_0_gain_;
  total_gain_ *= stage_0_gain_;
  p0_ = calin::math::least_squares::polyval(stage_0_pmf_, p0_);

  if(config_.suppress_zero()) {
    double renorm = 1.0/(1.0 - p0_);
    total_gain_ *= renorm;
    Exx *= renorm;
  }

  resolution_ = std::sqrt(Exx/SQR(total_gain_) - 1.0);
}

std::vector<double> LombardMartinPrescottPMTModel::
polya_pmf(double mean, double rms_frac, double precision)
{
  const double mu      = mean;
  const double b       = SQR(rms_frac);

  const double C1      = mu*b/(1.0 + mu*b);
  const double C2      = (mu - mu*b)/(1.0 + mu*b);

  const double pcutoff = 1.0 - precision;

  calin::math::accumulator::KahanAccumulator psum;

  std::vector<double> pk;
  pk.reserve(5*int(ceil(mu)));
  pk.resize(1);

  // Do one iteration of Prescott's algorithm (1965) for case of pk-1 = { 0, 1 }

  double pkx;
  if(b == 0) {
    pkx = std::exp(-mu);
  } else {
    pkx = std::pow(1.0+mu*b,-1.0/b);
  }
  pk[0] = pkx;
  psum.accumulate(pkx);

  for(unsigned ix=1;psum.total()<pcutoff;ix++)
  {
    pkx *= C1 + C2/double(ix);
    pk.push_back(pkx);
    psum.accumulate(pkx);
    if(psum.total() > 0.99 and pkx < precision)break;
  }
  return pk;
}

Eigen::VectorXd LombardMartinPrescottPMTModel::calc_spectrum(unsigned npoint,
  const std::vector<double>* pe_spec)
{
  int plan_flags = proto_planning_enum_to_fftw_flag(fftw_rigor_);

  // Assign a buffer to hold DFT of SES / MES
  uptr_fftw_data dft { fftw_alloc_real(npoint), fftw_free };
  assert(dft);

  // Assign a buffer to hold temporary data
  uptr_fftw_data buffer { fftw_alloc_real(npoint), fftw_free };
  assert(buffer);

  // Temporarily make zero-padded copy of single n-stage PMF in the buffer
  std::copy(stage_n_pmf_.begin(), stage_n_pmf_.end(), buffer.get());
  std::fill(buffer.get()+stage_n_pmf_.size(), buffer.get()+npoint, 0.0);

  // Do the out-of-place forward DFT transforming pmf to Fourier space in "dft"
  uptr_fftw_plan fwd_plan = {
    fftw_plan_r2r_1d(npoint, buffer.get(), dft.get(), FFTW_R2HC, plan_flags), fftw_destroy_plan };
  assert(fwd_plan);
  fftw_execute(fwd_plan.get());

  std::vector<const std::vector<double>*> all_stages;
  all_stages.reserve(config_.num_stage());

  // Do num_stage()-2 convolutions with the "n-stage" PMF
  std::vector<double> stage_0_pmf_adjusted;
  for(unsigned istage=0; istage<config_.num_stage()-2; istage++) {
    all_stages.push_back(&stage_n_pmf_);
  }

  // Then convolve with stage 0 PMF (zero suppress adjusted if necessary)
  all_stages.push_back(&stage_0_pmf_zsa_);

  if(pe_spec != nullptr) {
    all_stages.push_back(pe_spec);
  }

  // Do the convolutions
  hcvec_multi_stage_polynomial(dft.get(), dft.get(), all_stages, npoint);

  // Prepare the inverse DFT
  uptr_fftw_plan bwd_plan = {
    fftw_plan_r2r_1d(npoint, dft.get(), buffer.get(), FFTW_HC2R, plan_flags),
    fftw_destroy_plan };
  assert(bwd_plan);
  fftw_execute(bwd_plan.get());

  // Copy the spectrum to the output vector, normalizing for FFTW comventions
  Eigen::VectorXd spec(npoint);
  double norm = 1.0/double(npoint);
  std::transform(buffer.get(), buffer.get()+npoint, spec.data(),
    [norm](double x) { return x*norm; });

  return spec;
}

Eigen::VectorXd LombardMartinPrescottPMTModel::calc_ses(unsigned npoint)
{
  if(npoint==0) {
    npoint = calin::math::special::round_up_power_of_two(total_gain_ * 4.0);
  }
  return calc_spectrum(npoint);
}
