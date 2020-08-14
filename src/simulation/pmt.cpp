/*

   calin/simulation/pmt.cpp -- Stephen Fegan -- 2016-03-21

   Class to genereate random deviates from a PMT spectrum

   Copyright 2014, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <algorithm>
#include <cmath>
#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cctype>
#include <fftw3.h>

#include <calin_global_definitions.hpp>
#include <simulation/pmt.hpp>
#include <math/special.hpp>
#include <math/accumulator.hpp>
#include <util/log.hpp>
#include <math/fftw_util.hpp>
#include <provenance/chronicle.hpp>
#include <math/fftw_util.hpp>

using namespace calin::simulation::pmt;
using calin::math::special::SQR;
using calin::math::rng::RNG;
using namespace calin::util::log;
using namespace calin::math::fftw_util;

using uptr_fftw_plan = std::unique_ptr<fftw_plan_s,void(*)(fftw_plan_s*)>;
using uptr_fftw_data = std::unique_ptr<double,void(*)(void*)>;

SignalSource::~SignalSource()
{
  // nothing to see here
}

Eigen::VectorXd SignalSource::rvs(unsigned size)
{
  Eigen::VectorXd x(size);
  for(unsigned irv=0;irv<size;irv++)x[irv] = rv();
  return x;
}

calin::ix::simulation::pmt::PMTSimAbbreviatedConfig PMTSimPolya::cta_model_1()
{
  calin::ix::simulation::pmt::PMTSimAbbreviatedConfig config;
  config.set_num_stage(7);
  config.set_total_gain(4.0E4);
  config.set_stage_0_gain(10.0);
  config.set_suppress_zero(true);
  return config;
}

calin::ix::simulation::pmt::PMTSimAbbreviatedConfig PMTSimPolya::cta_model_2()
{
  auto config = cta_model_1();
  config.set_stage_0_gain_rms_frac(0.30);
  return config;
}

calin::ix::simulation::pmt::PMTSimAbbreviatedConfig PMTSimPolya::cta_model_3()
{
  auto config = cta_model_2();
  config.set_stage_0_prob_skip(0.1);
  return config;
}

PMTSimPolya::
PMTSimPolya(const calin::ix::simulation::pmt::PMTSimAbbreviatedConfig& config,
    math::rng::RNG* rng):
  SignalSource(), rng_(rng), my_rng_(0), napprox_limit_(50)
{
  if(config.num_stage()==0)
    throw std::runtime_error("PMTSimPolya: number of stages must be positive.");

  if(rng==nullptr)rng_ = my_rng_ = new RNG(__PRETTY_FUNCTION__);
  rng_->save_to_proto(config_.mutable_rng_config());

  config_.set_suppress_zero(config.suppress_zero());
  config_.set_signal_in_pe(config.signal_in_pe());
  auto* stage_0 = config_.add_stage();
  stage_0->set_gain_mean(config.stage_0_gain());
#if 0
  // Must decide what definition of s1_gain should be - either gain of
  // the dynode, or the average of the number of electrons entering
  // second stage (including skips)
  stage_0->set_gain_mean((s1_gain-s1_pskip)/(1.0-s1_pskip));
#endif
  stage_0->set_prob_skip(config.stage_0_prob_skip());

  if(config.stage_0_gain_rms_frac()>0)
  {
    stage_0->set_gain_rms_frac(config.stage_0_gain_rms_frac());
    gauss_a_.emplace_back(1.0/SQR(stage_0->gain_rms_frac()));
    gauss_b_.emplace_back(gauss_a_.back()/stage_0->gain_mean());
  }

  double s2_mean_n0 =
    stage_0->prob_skip() + stage_0->gain_mean()*(1.0-stage_0->prob_skip());
  double sn_gain = std::pow(config.total_gain()/s2_mean_n0,
    1.0/double(config.num_stage()-1));
  sn_gain = (sn_gain - config.stage_n_prob_skip())/(1.0 - config.stage_n_prob_skip());
  for(unsigned istage=1; istage<config.num_stage(); istage++)
  {
    auto* stage_n = config_.add_stage();
    stage_n->set_gain_mean(sn_gain);
    stage_n->set_prob_skip(config.stage_n_prob_skip());
    if(config.stage_n_gain_rms_frac()>0)
    {
      stage_n->set_gain_rms_frac(config.stage_n_gain_rms_frac());
      gauss_a_.emplace_back(1.0/SQR(stage_n->gain_rms_frac()));
      gauss_b_.emplace_back(gauss_a_.back()/stage_n->gain_mean());
    }
  }

  total_gain_ = 1;
  p0_         = 0;
  for(int istage=0; istage<config_.stage_size(); istage++)
  {
    total_gain_ *= stage_gain(istage);

    // This little mess of code implements the calculation of p0
    // from Prescott (1965) to account for the suppressed zeros in
    // the calculation of the total gain
    const auto& stage = config_.stage(config_.stage_size() - istage - 1);
    double p0_last = p0_;
    double mu      = stage.gain_mean();
    double b       = SQR(stage.gain_rms_frac());
    double bmu     = b*mu;
    double C1      = 1.0 + bmu*(1.0-p0_last);
    if(b == 0)
      p0_ = std::exp(mu*(p0_last-1.0));
    else
      p0_ = pow(C1,-1.0/b);
    p0_ = p0_*(1.0-stage.prob_skip()) + p0_last*stage.prob_skip();
  }
  if(config_.suppress_zero())
    total_gain_ /= 1.0-p0_;
}

double PMTSimPolya::stage_gain(unsigned istage) const
{
  const auto& stage = config_.stage(istage);
  double mu      = stage.gain_mean();
  return mu*(1.0-stage.prob_skip()) + stage.prob_skip();
}

PMTSimPolya::~PMTSimPolya()
{
  delete my_rng_;
}

double PMTSimPolya::rv()
{
  unsigned n;
do_over:
  n = 1;
  for(int istage=0; istage<config_.stage_size(); istage++)
  {
    const auto& stage = config_.stage(istage);
    unsigned n_in = n;

    if(stage.prob_skip() > 0)
  	{
  	  n = rng_->binomial(stage.prob_skip(), n_in);
  	  n_in -= n;
  	  if(n_in == 0)continue;
  	}
    else n = 0;

    double nmean = 0;
    if(stage.gain_rms_frac()>0)
  	{
  	  if(n_in > napprox_limit_)
  	    nmean = rng_->gamma_by_alpha_and_beta(double(n_in)*gauss_a_[istage],
  				 gauss_b_[istage]);
  	  else
  	    for(unsigned in_in = 0; in_in<n_in;in_in++)
  	      nmean +=
            rng_->gamma_by_alpha_and_beta(gauss_a_[istage],gauss_b_[istage]);
  	}
    else nmean = double(n_in)*stage.gain_mean();
    n += rng_->poisson(nmean);

    if(n==0)
  	{
  	  if(config_.suppress_zero())goto do_over;
  	  else return 0.0;
  	}
  }
  if(config_.signal_in_pe())return double(n)/total_gain_;
  else return double(n);
}

// Implements the calculation of the PDF for the SES from Prescott
// (1965). The calculation is O(n^2) in the number of total entries in
// the final SES, which is itself some factor (4-5) times the total
// gain. This calculation could also be done with FFTs which may be
// faster. Interestingly the FFT was "rediscovered" in ~1965, so
// Prescott may not have known about it.

calin::ix::simulation::pmt::PMTSimPMF
PMTSimPolya::calc_pmf(double precision, bool log_progress) const
{
  std::vector<double> pk(2);
  pk[0] = 0;
  pk[1] = 1;
  for(int ik=0;ik<config_.stage_size();ik++)
  {
    const auto& stage(config_.stage(config_.stage_size()-ik-1));

    std::vector<double> pkmo(pk);

    double mu      = stage.gain_mean();
    double b       = SQR(stage.gain_rms_frac());
    double bmu     = b*mu;
    double bmo     = b-1.0;
    double C1      = 1.0 + bmu*(1.0-pkmo[0]);
    double pcutoff = 1.0 - precision*((ik==config_.stage_size()-1)?0.1:1.0);
    double psum;

    pk.resize(1);
    pk.reserve(pkmo.size()*5*int(ceil(mu)));
    if(b == 0)
      pk[0] = std::exp(mu*(pkmo[0]-1.0));
    else
      pk[0] = pow(C1,-1.0/b);
    psum = pk[0];

    for(unsigned ix=1;psum<pcutoff;ix++)
  	{
  	  //std::cout << "- " << ix << ' ' << psum << '\n';
  	  double pkx = 0;
  	  for(unsigned ii=ix-std::min(ix,unsigned(pkmo.size()-1));ii<ix;ii++)
  	    {
  	      //	      std::cout << "-- " <<  ii << '\n';
  	      //	      if(ix-ii < pkmo.size())
      		pkx += pk[ii]*pkmo[ix-ii]*(double(ix)+double(ii)*bmo);
  	    }
  	  pkx *= mu/double(ix)/C1;
  	  psum += pkx;
  	  pk.push_back(pkx);
  	  if(pkx/psum < 1e-10)break;
  	}

    if(stage.prob_skip()>0)
  	{
  	  for(unsigned ix=0;ix<pk.size();ix++)
  	    pk[ix] *= (1.0-stage.prob_skip())/psum;
  	  for(unsigned ix=0;ix<std::min(pkmo.size(),pk.size());ix++)
  	    pk[ix] += stage.prob_skip()*pkmo[ix];
  	}
    else
  	{
  	  for(unsigned ix=0;ix<pk.size();ix++)
  	    pk[ix] /= psum;
  	}

    if(log_progress)
    {
      double pn = 0;
      double pnn = 0;
      unsigned ix0 = config_.suppress_zero() ? 1 : 0;
      double psum0 = config_.suppress_zero() ? (psum-pk[0]) : psum;
      for(unsigned ix=ix0;ix<pk.size();ix++)
        pn += double(ix)*pk[ix], pnn += SQR(double(ix))*pk[ix];
      LOG(INFO) << "Stage " << config_.stage_size()-ik << ", p0=" << pk[0]
        << ", <x>=" << pn/psum0
        << " res(x)=" << sqrt(pnn*psum0/SQR(pn)-1);
    }

  }

  calin::ix::simulation::pmt::PMTSimPMF OUTPUT;
  OUTPUT.set_suppress_zero(config_.suppress_zero());
  OUTPUT.set_signal_in_pe(config_.signal_in_pe());
  OUTPUT.mutable_pn()->Resize(pk.size(),0);
  std::copy(pk.begin(), pk.end(), OUTPUT.mutable_pn()->begin());
  return OUTPUT;
}

PMTSimTwoPopulation::
PMTSimTwoPopulation(const calin::ix::simulation::pmt::PMTSimTwoPopulationConfig& config,
    math::rng::RNG* rng, bool use_new_stage_n_algorithm, bool adopt_rng):
  config_(config), rng_(rng), adopt_rng_(adopt_rng),
  use_new_stage_n_algorithm_(use_new_stage_n_algorithm)
{
  if(config.num_stage()==0)
    throw std::runtime_error("PMTSimTwoPopulation: number of stages must be positive.");
  if(config.total_gain()<=0)
    throw std::runtime_error("PMTSimTwoPopulation: total gain must be positive.");
  if(config.stage_0_hi_gain()<=0)
    throw std::runtime_error("PMTSimTwoPopulation: stage 0 high-charge gain must be positive.");
  if(config.stage_0_hi_gain_rms_frac()<0)
    throw std::runtime_error("PMTSimTwoPopulation: stage 0 high-charge excess RMS must be zero or positive.");
  if(config.stage_0_lo_prob()<0 or config.stage_0_lo_prob()>=1)
    throw std::runtime_error("PMTSimTwoPopulation: stage 0 low-charge probability must be in the range [0,1).");
  if(config.stage_0_lo_prob()>0) {
    if(config.stage_0_lo_gain()<=0)
      throw std::runtime_error("PMTSimTwoPopulation: stage 0 low-charge gain must be positive.");
    if(config.stage_0_lo_gain_rms_frac()<0)
      throw std::runtime_error("PMTSimTwoPopulation: stage 0 low-charge excess RMS must be zero or positive.");
  }
  if(config.stage_n_gain_rms_frac()<0)
    throw std::runtime_error("PMTSimTwoPopulation: stage 1+ excess RMS must be zero or positive.");

  if(config_.stage_0_hi_gain_rms_frac() > 0) {
    gamma_a_0_hi_ = 1.0/SQR(config_.stage_0_hi_gain_rms_frac());
    gamma_b_0_hi_ = gamma_a_0_hi_/config_.stage_0_hi_gain();
  }
  if(config_.stage_0_lo_gain_rms_frac() > 0) {
    gamma_a_0_lo_ = 1.0/SQR(config_.stage_0_lo_gain_rms_frac());
    gamma_b_0_lo_ = gamma_a_0_lo_/config_.stage_0_lo_gain();
  }

  stage_n_gain_ = std::pow(config.total_gain()/stage_0_gain(), 1/double(config_.num_stage()-1));
  recalc_total_gain_and_p0();
  unsigned iter = 10;
  while(std::abs(total_gain_/config.total_gain() - 1)>1e-8) {
    if(iter == 0) {
      throw std::runtime_error("PMTSimTwoPopulation: iteration limit reached in calculating stage gain");
    }
    --iter;
    stage_n_gain_ *= std::pow(config.total_gain()/total_gain_, 1/double(config_.num_stage()-1));
    recalc_total_gain_and_p0();
  }

  if(config_.stage_n_gain_rms_frac() > 0) {
    gamma_a_n_    = 1.0/SQR(config_.stage_n_gain_rms_frac());
    gamma_b_n_    = gamma_a_n_/stage_n_gain_;

    if(use_new_stage_n_algorithm_) {
      stage_n_x1_cdf_ = stage_pmf({0.0, 1.0}, stage_n_gain_, config_.stage_n_gain_rms_frac(), 1e-8);
      stage_n_x2_cdf_ = stage_pmf(stage_n_x1_cdf_, stage_n_gain_, config_.stage_n_gain_rms_frac(), 1e-8);
      stage_n_x3_cdf_ = stage_pmf(stage_n_x2_cdf_, stage_n_gain_, config_.stage_n_gain_rms_frac(), 1e-8);

      renorm_pmf(stage_n_x1_cdf_);
      renorm_pmf(stage_n_x2_cdf_);
      renorm_pmf(stage_n_x3_cdf_);

      double psum = 0;
      for(auto & ibin : stage_n_x1_cdf_) { psum += ibin; ibin = psum; }
      if(psum<1.0) { stage_n_x1_cdf_.push_back(1.0) ; }

      psum = 0;
      for(auto & ibin : stage_n_x2_cdf_) { psum += ibin; ibin = psum; }
      if(psum<1.0) { stage_n_x2_cdf_.push_back(1.0) ; }

      psum = 0;
      for(auto & ibin : stage_n_x3_cdf_) { psum += ibin; ibin = psum; }
      if(psum<1.0) { stage_n_x3_cdf_.push_back(1.0) ; }
    }
  }
}

PMTSimTwoPopulation::~PMTSimTwoPopulation()
{
  if(adopt_rng_)delete rng_;
}

double PMTSimTwoPopulation::
stage_p0(double p0_last, double gain, double gain_rms_frac)
{
  // This little mess of code implements one iteration of the calculation of p0
  // from Prescott (1965) to account for the suppressed zeros in the calculation
  // of the total gain

  double mu      = gain;
  double b       = SQR(gain_rms_frac);
  double bmu     = b*mu;
  double C1      = 1.0 + bmu*(1.0-p0_last);
  double p0;
  if(b == 0) {
    p0 = std::exp(mu*(p0_last-1.0));
  } else {
    p0 = pow(C1,-1.0/b);
  }
  return p0;
}

void PMTSimTwoPopulation::recalc_total_gain_and_p0()
{
  total_gain_ = 1;
  p0_         = 0;

  for(unsigned istage=1; istage<config_.num_stage(); istage++)
  {
    total_gain_ *= stage_n_gain_;
    p0_ = stage_p0(p0_, stage_n_gain_, config_.stage_n_gain_rms_frac());
  }

  total_gain_ *= stage_0_gain();
  double p0_hi = stage_p0(p0_, config_.stage_0_hi_gain(), config_.stage_0_hi_gain_rms_frac());
  if(config_.stage_0_lo_prob() > 0) {
    double p0_lo = stage_p0(p0_, config_.stage_0_lo_gain(), config_.stage_0_lo_gain_rms_frac());
    p0_ = config_.stage_0_lo_prob()*p0_lo + (1-config_.stage_0_lo_prob())*p0_hi;
  } else {
    p0_ = p0_hi;
  }

  if(config_.suppress_zero())
    total_gain_ /= 1.0-p0_;
}

double PMTSimTwoPopulation::rv()
{
do_over:
  double nmean;
  if(config_.stage_0_lo_prob()>0 and rng_->uniform()<config_.stage_0_lo_prob()) {
    if(config_.stage_0_lo_gain_rms_frac() > 0) {
      nmean = rng_->gamma_by_alpha_and_beta(gamma_a_0_lo_, gamma_b_0_lo_);
    } else {
      nmean = config_.stage_0_lo_gain();
    }
  } else {
    if(config_.stage_0_hi_gain_rms_frac() > 0) {
      nmean = rng_->gamma_by_alpha_and_beta(gamma_a_0_hi_, gamma_b_0_hi_);
    } else {
      nmean = config_.stage_0_hi_gain();
    }
  }

  unsigned n = rng_->poisson(nmean);
  if(n==0)
  {
    if(config_.suppress_zero()) {
      goto do_over;
    } else {
      return 0.0;
    }
  }

  if(config_.stage_n_gain_rms_frac() == 0) {
    n = stage_n_poisson(n);
  } else if(use_new_stage_n_algorithm_) {
    n = stage_n_new(n);
  } else {
    n = stage_n_old(n);
  }

  if(n==0 and config_.suppress_zero()) {
    goto do_over;
  }

  if(config_.signal_in_pe()) {
    return double(n)/total_gain_;
  } else {
    return double(n);
  }
}

unsigned PMTSimTwoPopulation::stage_n_new(unsigned n) const
{
  unsigned istage = 1;

  if((config_.num_stage()-1) % 2 == 1) {
    unsigned n_in = n;
    n = 0;
    for(unsigned i = 0; i<n_in; i++) {
      n += std::upper_bound(stage_n_x1_cdf_.begin(), stage_n_x1_cdf_.end(),
        rng_->uniform()) - stage_n_x1_cdf_.begin();
    }
    istage += 1;
  } else if((config_.num_stage()-1) % 2 == 2) {
    unsigned n_in = n;
    n = 0;
    for(unsigned i = 0; i<n_in; i++) {
      n += std::upper_bound(stage_n_x2_cdf_.begin(), stage_n_x2_cdf_.end(),
        rng_->uniform()) - stage_n_x2_cdf_.begin();
    }
    istage += 2;
  }

  if(n==0)return 0;

  for(; istage<config_.num_stage(); istage+=3)
  {
    unsigned n_in = n;
    n = 0;
    for(unsigned i = 0; i<n_in; i++) {
      n += std::upper_bound(stage_n_x3_cdf_.begin(), stage_n_x3_cdf_.end(),
        rng_->uniform()) - stage_n_x3_cdf_.begin();
    }
    if(n==0)return 0;
  }
  return n;
}

unsigned PMTSimTwoPopulation::stage_n_old(unsigned n) const
{
  for(unsigned istage=1; istage<config_.num_stage(); istage++)
  {
    double nmean = 0;
    for(unsigned i = 0; i<n;i++) {
      nmean += rng_->gamma_by_alpha_and_beta(gamma_a_n_,gamma_b_n_);
    }
    n = rng_->poisson(nmean);
    if(n==0)break;
  }
  return n;
}

unsigned PMTSimTwoPopulation::stage_n_poisson(unsigned n) const
{
  for(unsigned istage=1; istage<config_.num_stage(); istage++)
  {
    n = rng_->poisson(double(n)*stage_n_gain_);
    if(n==0)break;
  }
  return n;
}

std::vector<double> PMTSimTwoPopulation::stage_pmf(
  const std::vector<double>& pkmo, double gain, double gain_rms_frac, double precision) const
{
  // Do one iteration of Prescott's algorithm (1965).

  double mu      = gain;
  double b       = SQR(gain_rms_frac);
  double bmu     = b*mu;
  double bmo     = b-1.0;
  double C1      = 1.0 + bmu*(1.0-pkmo[0]);
  double pcutoff = 1.0 - precision;
  calin::math::accumulator::KahanAccumulator psum;

  std::vector<double> pk;
  pk.reserve(pkmo.size()*5*int(ceil(mu)));
  pk.resize(1);
  if(b == 0)
    pk[0] = std::exp(mu*(pkmo[0]-1.0));
  else
    pk[0] = std::pow(C1,-1.0/b);
  psum.accumulate(pk[0]);

  for(unsigned ix=1;psum.total()<pcutoff;ix++)
  {
    double ix_dbl = double(ix);
    calin::math::accumulator::KahanAccumulator pkx_acc;
    unsigned ii0 = ix-std::min(ix,unsigned(pkmo.size()-1));
    double ii_dbl = ix_dbl + double(ii0)*bmo;
    for(unsigned ii=ii0; ii<ix; ++ii, ii_dbl+=bmo) {
      nflop_ += 2;
      // LOG(INFO) << ix << ' ' << ii << ' ' << ix-ii;
      pkx_acc.accumulate(pk[ii]*pkmo[ix-ii]*ii_dbl);
    }
    double pkx = pkx_acc.total()*mu/ix_dbl/C1;
    psum.accumulate(pkx);
    pk.push_back(pkx);
    if(psum.total() > 0.99 and pkx < precision)break;
  }

  return pk;
}

void PMTSimTwoPopulation::renorm_pmf(std::vector<double>& pk)
{
  calin::math::accumulator::KahanAccumulator psum = 0;
  for(unsigned ix=0;ix<pk.size();ix++) {
    psum.accumulate(pk[ix]);
  }
  double norm = 1.0/psum.total();
  for(unsigned ix=0;ix<pk.size();ix++) {
    pk[ix] *= norm;
  }
}

std::vector<double> PMTSimTwoPopulation::stage_0_pmf(double precision) const
{
  std::vector<double> p0 = stage_0_hi_pmf();
  if(config_.stage_0_lo_prob() > 0) {
    std::vector<double> p0_lo = stage_0_lo_pmf();
    for(unsigned ipe=0; ipe<p0.size(); ipe++) {
      p0[ipe] *= (1-config_.stage_0_lo_prob());
    }
    if(p0_lo.size() > p0.size()) {
      p0.resize(p0_lo.size());
    }
    for(unsigned ipe=0; ipe<p0_lo.size(); ipe++) {
      p0[ipe] += config_.stage_0_lo_prob()*p0_lo[ipe];
    }
  }
  return p0;
}

std::string PMTSimTwoPopulation::stage_summary(unsigned istage, const std::vector<double>& pk) const
{
  calin::math::accumulator::KahanAccumulator p1;
  calin::math::accumulator::KahanAccumulator pn;
  calin::math::accumulator::KahanAccumulator pnn;
  unsigned ix0 = config_.suppress_zero() ? 1 : 0;
  for(unsigned ix=ix0;ix<pk.size();ix++) {
    p1.accumulate(pk[ix]);
    pn.accumulate(double(ix)*pk[ix]);
    pnn.accumulate(SQR(double(ix))*pk[ix]);
  }
  double ptot = (config_.suppress_zero() ? pk[0] : 0.0) + p1.total();
  std::ostringstream stream;
  stream << "Stage " << istage
    << ", p(0)=" << pk[0]
    << ", <x>=" << pn.total()/p1.total()
    << " res(x)=" << sqrt(pnn.total()*p1.total()/SQR(pn.total())-1)
    << ", 1-norm=" << 1-ptot
    << ", N=" << pk.size()
    ;
  return stream.str();
}

std::string PMTSimTwoPopulation::
fft_log_progress(unsigned istage, double* fk, unsigned npoint, int plan_flags) const
{
  // FFT of multi-stage PMF including the k lowest stages
  uptr_fftw_data fk_copy { fftw_alloc_real(npoint), fftw_free };
  assert(fk_copy);

  uptr_fftw_data pk { fftw_alloc_real(npoint), fftw_free };
  assert(pk);

  // Prepare the backward DFT
  uptr_fftw_plan bwd_plan = {
    fftw_plan_r2r_1d(npoint, fk_copy.get(), pk.get(),
                    FFTW_HC2R , plan_flags), fftw_destroy_plan };
  assert(bwd_plan);

  std::copy(fk, fk+npoint, fk_copy.get());

  fftw_execute(bwd_plan.get());

  double p1 = 0;
  double pi = 0;
  double pii = 0;
  for(unsigned i=(config_.suppress_zero()?1:0); i<npoint; i++) {
    p1 += pk.get()[i];
    pi += double(i)*pk.get()[i];
    pii += SQR(double(i))*pk.get()[i];
  }
  std::ostringstream stream;
  stream << "Stage " << istage
    << ", p(0)=" << *pk/npoint << ", <x>=" << pi/p1
    //<< ", EVF(x)=" << pii*p1/(pi*pi)
    << ", res(x)=" << sqrt(pii*p1/(pi*pi) - 1)
    << ", 1-norm=" << 1.0-(p1 + (config_.suppress_zero()?(*pk):0))/double(npoint)
    << ", nflop=" << nflop_;
  return stream.str();
}

// Fast function to calculate PMF using FFTs
calin::ix::simulation::pmt::PMTSimPMF
PMTSimTwoPopulation::calc_pmf_fft(unsigned npoint, unsigned nstage, double precision,
  bool log_progress, bool skip_inverse_fft,
  calin::ix::math::fftw_util::FFTWPlanningRigor fftw_rigor) const
{
  int plan_flags = proto_planning_enum_to_fftw_flag(fftw_rigor);
  nflop_ = 0;

  std::string stage_summary;

  if(nstage == 0) {
    nstage = config_.num_stage();
  }

  if(nstage == 1) {
    return calc_pmf_prescott(nstage, precision);
  }

  if(npoint==0) {
    npoint = int(total_gain_ * 4.0);
    unsigned log_npoint = 0;
    while(npoint) {
      npoint >>= 1;
      ++log_npoint;
    }
    npoint = 1 << log_npoint;
  }

  // FFT of multi-stage PMF including the k lowest stages
  uptr_fftw_data fk { fftw_alloc_real(npoint), fftw_free };
  assert(fk);

  // FFT of multi-stage PMF including the (k - 1) lowest stages
  uptr_fftw_data fkmo { fftw_alloc_real(npoint), fftw_free };
  assert(fkmo);

  // Get the stage-n PMF
  std::vector<double> pn =
    stage_pmf({0.0, 1.0}, stage_n_gain_, config_.stage_n_gain_rms_frac(), precision);
  renorm_pmf(pn);
  assert(pn.size() < npoint);

  // Temporarily make zero-padded copy of single n-stage PMF in fk
  std::fill(fk.get(), fk.get()+npoint, 0.0);
  std::copy(pn.begin(), pn.end(), fk.get());

  // Do the forward DFT transforming pn into fkmo
  uptr_fftw_plan fwd_plan = {
    fftw_plan_r2r_1d(npoint, fk.get(), fkmo.get(),
                    FFTW_R2HC, plan_flags), fftw_destroy_plan };
  assert(fwd_plan);
  fftw_execute(fwd_plan.get());

  double nadd;
  double nmul;
  double nfma;
  fftw_flops(fwd_plan.get(), &nadd, &nmul, &nfma);
  nflop_ += uint64_t(nmul) + uint64_t(nfma);

  if(log_progress) {
    std::string s = fft_log_progress(nstage-1, fkmo.get(), npoint, plan_flags);
    LOG(INFO) << s;
    stage_summary += s;
  }

  // Do the convolutions nstage-2 times from fkmo to fk, swapping back each time
  for(unsigned ik=1; ik<(nstage-1); ik++)
  {
    hcvec_polynomial(fk.get(), fkmo.get(), pn, npoint);
    nflop_ += pn.size() * npoint * 2;
    std::swap(fk, fkmo);

    if(log_progress) {
      std::string s = fft_log_progress(nstage-1-ik, fkmo.get(), npoint, plan_flags);
      LOG(INFO) << s;
      stage_summary += s;
    }
  }

  // Calculate the stage-0 PMF
  std::vector<double> p0 = stage_0_pmf();
  renorm_pmf(p0);

  // Do the final stage convolution
  hcvec_polynomial(fk.get(), fkmo.get(), p0, npoint);
  nflop_ += p0.size() * npoint * 2;
  std::swap(fk, fkmo);

  if(log_progress) {
    std::string s = fft_log_progress(0, fkmo.get(), npoint, plan_flags);
    LOG(INFO) << s;
    stage_summary += s;
  }

  calin::ix::simulation::pmt::PMTSimPMF OUTPUT;
  OUTPUT.set_suppress_zero(config_.suppress_zero());
  OUTPUT.set_signal_in_pe(config_.signal_in_pe());
  OUTPUT.mutable_pn()->Resize(npoint,0);
  OUTPUT.set_stage_statistics_summary(stage_summary);

  if(skip_inverse_fft) {
    std::copy(fkmo.get(), fkmo.get()+npoint, OUTPUT.mutable_pn()->begin());
  } else {
    // Prepare the backward DFT
    uptr_fftw_plan bwd_plan = {
      fftw_plan_r2r_1d(npoint, fkmo.get(), fk.get(),
                      FFTW_HC2R , plan_flags), fftw_destroy_plan };
    assert(bwd_plan);
    fftw_execute(bwd_plan.get());

    fftw_flops(bwd_plan.get(), &nadd, &nmul, &nfma);
    nflop_ += uint64_t(nmul) + uint64_t(nfma);

    const double norm = 1/double(npoint);
    std::transform(fk.get(), fk.get()+npoint, OUTPUT.mutable_pn()->begin(),
      [norm](double x) { return x*norm; });
    nflop_ += npoint;
  }

  return OUTPUT;
}

// Slow function to calculate PMF using Prescott (1965).
calin::ix::simulation::pmt::PMTSimPMF
PMTSimTwoPopulation::calc_pmf_prescott(unsigned nstage, double precision, bool log_progress) const
{
  nflop_ = 0;
  std::string stats_summary;

  if(nstage == 0)nstage = config_.num_stage();

  std::vector<double> pk { 0.0, 1.0 };
  for(unsigned ik=1;ik<nstage;ik++)
  {
    pk = stage_pmf(pk, stage_n_gain_, config_.stage_n_gain_rms_frac(), precision);

    std::string this_stage_summary = stage_summary(nstage-ik, pk);
    stats_summary += this_stage_summary + '\n';
    if(log_progress) {
      LOG(INFO) << this_stage_summary;
    }
  }

  std::vector<double> pk_hi =
    stage_pmf(pk, config_.stage_0_hi_gain(), config_.stage_0_hi_gain_rms_frac(), precision*0.1);
  if(config_.stage_0_lo_prob() > 0) {
    std::vector<double> pk_lo =
      stage_pmf(pk, config_.stage_0_lo_gain(), config_.stage_0_lo_gain_rms_frac(), precision*0.1);

    double p_lo = config_.stage_0_lo_prob();
    double p_hi = 1-p_lo;

    if(pk_hi.size() < pk_lo.size()) {
      std::swap(pk_hi, pk_lo);
      std::swap(p_hi, p_lo);
    }

    for(unsigned i=0;i<pk_hi.size();i++)pk_hi[i] *= p_hi;
    for(unsigned i=0;i<pk_lo.size();i++)pk_hi[i] += pk_lo[i] * p_lo;
  };

  std::swap(pk, pk_hi);

  std::string this_stage_summary = stage_summary(0, pk);
  stats_summary += this_stage_summary + '\n';
  if(log_progress) {
    LOG(INFO) << this_stage_summary;
  }

  renorm_pmf(pk);

  calin::ix::simulation::pmt::PMTSimPMF OUTPUT;
  OUTPUT.set_suppress_zero(config_.suppress_zero());
  OUTPUT.set_signal_in_pe(config_.signal_in_pe());
  OUTPUT.mutable_pn()->Resize(pk.size(),0);
  std::copy(pk.begin(), pk.end(), OUTPUT.mutable_pn()->begin());
  OUTPUT.set_stage_statistics_summary(stats_summary);
  return OUTPUT;
}

calin::ix::simulation::pmt::PMTSimPMF PMTSimTwoPopulation::calc_multi_electron_spectrum(
  double intensity_mean, double intensity_rms_frac,
  const Eigen::VectorXd& ped_hc_dft,
  unsigned npoint, unsigned nstage,
  double precision, bool log_progress,
  calin::ix::math::fftw_util::FFTWPlanningRigor fftw_rigor) const
{
  int plan_flags = proto_planning_enum_to_fftw_flag(fftw_rigor);
  double nadd;
  double nmul;
  double nfma;

  if(npoint == 0) {
    if(ped_hc_dft.size()==0) {
      npoint = int(std::min(1.0, intensity_mean) * total_gain_ * 4.0);
      unsigned log_npoint = 0;
      while(npoint) {
        npoint >>= 1;
        ++log_npoint;
      }
      npoint = 1 << log_npoint;
    } else {
      npoint = ped_hc_dft.size();
    }
  }

  // Calculate the single electron spectrum, leaving it in Fourier space
  calin::ix::simulation::pmt::PMTSimPMF OUTPUT = calc_pmf_fft(npoint, nstage,
    precision, log_progress, /* skip_inverse_fft = */ true, fftw_rigor);

  if(OUTPUT.suppress_zero()) {
    double p0 = hcvec_avg_real(OUTPUT.pn().begin(), npoint);
    hcvec_scale_and_add_real(OUTPUT.mutable_pn()->begin(), 1/(1-p0), -p0/(1-p0), npoint);
  }

  // FFT of multi-electron spectrum
  uptr_fftw_data fmes { fftw_alloc_real(npoint), fftw_free };
  assert(fmes);

  // Get the photo-electron PMF
  std::vector<double> ppe =
    stage_pmf({0.0, 1.0}, intensity_mean, intensity_rms_frac, precision);
  renorm_pmf(ppe);
  assert(ppe.size() < npoint);

  // Do the PE convolution
  hcvec_polynomial(fmes.get(), OUTPUT.pn().begin(), ppe, npoint);
  nflop_ += ppe.size() * npoint * 2;

  // Add the arbitrary pedestal noise if requested
  if(ped_hc_dft.size() > 0) {
    if(ped_hc_dft.size() != npoint) {
      throw std::runtime_error("Pedestal PDF must be same size as npoint.");
    }
    hcvec_scale_and_multiply(fmes.get(), fmes.get(), ped_hc_dft.data(), npoint);
    nflop_ += npoint * 2;
  }

  // Prepare the backward DFT
  uptr_fftw_plan bwd_plan = {
    fftw_plan_r2r_1d(npoint, fmes.get(), fmes.get(), FFTW_HC2R , plan_flags),
    fftw_destroy_plan };
  assert(bwd_plan);
  fftw_execute(bwd_plan.get());

  fftw_flops(bwd_plan.get(), &nadd, &nmul, &nfma);
  nflop_ += uint64_t(nmul) + uint64_t(nfma);

  const double norm = 1/double(npoint);
  std::transform(fmes.get(), fmes.get()+npoint, OUTPUT.mutable_pn()->begin(),
    [norm](double x) { return x*norm; });
  nflop_ += npoint;

  return OUTPUT;
}

calin::ix::simulation::pmt::PMTSimTwoPopulationConfig PMTSimTwoPopulation::cta_model_4()
{
  calin::ix::simulation::pmt::PMTSimTwoPopulationConfig cfg;
  cfg.set_num_stage(7);
  cfg.set_total_gain(40000);
  cfg.set_stage_0_hi_gain(12);
  cfg.set_stage_0_hi_gain_rms_frac(0);
  cfg.set_stage_0_lo_gain(3);
  cfg.set_stage_0_lo_gain_rms_frac(0);
  cfg.set_stage_0_lo_prob(0.15);
  cfg.set_stage_n_gain_rms_frac(0.0);
  cfg.set_suppress_zero(true);
  return cfg;
}

// =============================================================================
//
// PMTSimInvCDF
//
// =============================================================================

PMTSimInvCDF::PMTSimInvCDF(const calin::ix::simulation::pmt::PMTSimPMF& cmf,
    unsigned npoint, math::rng::RNG* rng):
  SignalSource(), rng_(rng), my_rng_(nullptr)
{
  if(rng==nullptr)rng_ = my_rng_ = new math::rng::RNG(__PRETTY_FUNCTION__);

  int iy0 = 0;
  if(cmf.suppress_zero())iy0 = 1;
  inv_cdf_.resize(cmf.pn_size() + 1 - iy0);
  double ysum = 0;
  double ynsum = 0;
  for(int iy=iy0;iy<cmf.pn_size();iy++)
    ysum += cmf.pn(iy), ynsum += double(iy)*cmf.pn(iy);
  double gain = cmf.signal_in_pe() ? (ysum/ynsum) : 1.0;
  double ycumsum = 0;
  inv_cdf_[0] = std::make_pair<double,double>((double(iy0)-0.5)*gain,0.0);
  for(int iy=iy0;iy<cmf.pn_size();iy++)
  {
    ycumsum += cmf.pn(iy);
    inv_cdf_[iy + 1 - iy0] =
      std::make_pair<double,double>((double(iy)+0.5)*gain,ycumsum/ysum);
  }
#if 0
  for(unsigned iy=0;iy<10;iy++)
    std::cout << m_inv_cdf[iy].first << ' ' << m_inv_cdf[iy].second << '\n';
#endif
  rng_->generate_inverse_cdf(inv_cdf_, npoint);
}

PMTSimInvCDF::PMTSimInvCDF(const std::string& filename, math::rng::RNG* rng):
  SignalSource(), rng_(rng), my_rng_(nullptr)
{
  if(rng==nullptr)rng_ = my_rng_ = new math::rng::RNG(__PRETTY_FUNCTION__);
  std::ifstream s(filename.c_str());
  if(!s)throw std::runtime_error(
    std::string("PMTSimInvCDF: Could not open file: ")+filename);
  auto* file_record = calin::provenance::chronicle::register_file_open(filename,
    calin::ix::provenance::chronicle::AT_READ, __PRETTY_FUNCTION__);
  std::string line;
  std::string comment;
  std::getline(s,line);
  while(s)
  {
    unsigned ic=0;
    while(ic<line.size() && isspace(line[ic]))ic++;
    if(ic==line.size())goto next_line;
    if(line[ic] == '#') {
      comment += line.substr(ic);
      comment += '\n';
      goto next_line;
    }
    double x,y;
    if(std::istringstream(line) >> x >> y)
    inv_cdf_.emplace_back(x,y);
  next_line:
    std::getline(s,line);
  }
  file_record->set_comment(comment);
  calin::provenance::chronicle::register_file_close(file_record);
}

PMTSimInvCDF::~PMTSimInvCDF()
{
  delete my_rng_;
}

double PMTSimInvCDF::rv()
{
  return rng_->inverse_cdf(inv_cdf_);
}

void PMTSimInvCDF::
save_inv_cdf_to_file(const std::string& filename,
  const std::string& comment) const
{
  std::ofstream s(filename.c_str());
  if(!s)throw std::runtime_error(
    std::string("PMTSimInvCDF: Could not open file: ")+filename);
  auto* file_record = calin::provenance::chronicle::register_file_open(filename,
    calin::ix::provenance::chronicle::AT_WRITE, __PRETTY_FUNCTION__, comment);
  time_t rawtime;
  time (&rawtime);
  s << "# Inverse CDF written " << ctime(&rawtime);
  if(!comment.empty())s << "# " << comment << '\n';
  for(unsigned iy=0;iy<inv_cdf_.size();iy++)
    s << std::to_string(inv_cdf_[iy].first) << ' '
      << std::to_string(inv_cdf_[iy].second) <<  '\n';
  file_record->set_file_size(s.tellp());
  calin::provenance::chronicle::register_file_close(file_record);
}

PMTSimGammaDistribution::PMTSimGammaDistribution(double alpha, double beta,
    calin::math::rng::RNG* rng, bool adopt_rng):
  SignalSource(), alpha_(alpha), beta_(beta),
  rng_((rng==nullptr)?new math::rng::RNG(__PRETTY_FUNCTION__) : rng),
  adopt_rng_((rng==nullptr)?true:adopt_rng)
{
  // nothing to see here
}

PMTSimGammaDistribution::~PMTSimGammaDistribution()
{
  if(adopt_rng_)delete rng_;
}

double PMTSimGammaDistribution::rv()
{
  return rng_->gamma_by_alpha_and_beta(alpha_, beta_);
}

// =============================================================================
//
// MultiPESpectrum
//
// =============================================================================

MultiPESpectrum::
MultiPESpectrum(SignalSource* pmt,
    const calin::ix::simulation::pmt::MultiPESpectrumConfig& config,
    math::rng::RNG* rng):
  SignalSource(),  pmt_(pmt), rng_(rng), my_rng_(0), config_(config),
  signal_gamma_a_(), signal_gamma_b_()
{
  if(rng==nullptr)rng_ = my_rng_ = new math::rng::RNG(__PRETTY_FUNCTION__);
  if(config_.signal_rms_frac()>0)
    {
      signal_gamma_a_ = 1.0/SQR(config_.signal_rms_frac());
      signal_gamma_b_ = signal_gamma_a_/config_.signal_mean();
    }
}

MultiPESpectrum::~MultiPESpectrum()
{
  delete my_rng_;
}

double MultiPESpectrum::rv()
{
  double x;
  double npe_mean;
  int npe;
  rv_xnm(x,npe,npe_mean);
  return x;
}

void MultiPESpectrum::rv_xnm(double& x, int& n, double& m)
{
  x = rng_->normal()*config_.pedestal_rms() + config_.pedestal_mean();
  if(config_.signal_rms_frac()>0)
    m = rng_->gamma_by_alpha_and_beta(signal_gamma_a_,signal_gamma_b_);
  else
    m = config_.signal_mean();
  n = rng_->poisson(m);
  for(int ipe=0;ipe<n;ipe++)x += pmt_->rv();
}

void MultiPESpectrum::
rvs_xnm(VecRef x, IntVecRef n, VecRef m, unsigned size)
{
  x.resize(size);
  n.resize(size);
  m.resize(size);
  for(unsigned irv=0;irv<size;irv++)rv_xnm(x[irv],n[irv],m[irv]);
}

double MultiPESpectrum::rv_ped()
{
  return rng_->normal()*config_.pedestal_rms() + config_.pedestal_mean();
}

Eigen::VectorXd MultiPESpectrum::rvs_ped(unsigned size)
{
  Eigen::VectorXd x(size);
  for(unsigned irv=0;irv<size;irv++)x[irv]=rv_ped();
  return x;
}

PoissonSignalSim::PoissonSignalSim(SignalSource* pmt,
    const calin::ix::simulation::pmt::PoissonSignalSimConfig& config,
    math::rng::RNG* rng):
  pmt_(pmt), rng_(rng), my_rng_(0), config_(config)
{
  if(rng==nullptr)rng_ = my_rng_ = new math::rng::RNG(__PRETTY_FUNCTION__);
}

PoissonSignalSim::~PoissonSignalSim()
{
  delete my_rng_;
}

double PoissonSignalSim::rv(double lambda)
{
  double x = 0;
  if(config_.pedestal_rms())x = rng_->normal()*config_.pedestal_rms();
  x += config_.pedestal_mean();
  int npe = rng_->poisson(lambda);
  for(int ipe=0;ipe<npe;ipe++)x += pmt_->rv();
  return x;
}

Eigen::VectorXd PoissonSignalSim::rvs(double lambda, unsigned size)
{
  Eigen::VectorXd x(size);
  for(unsigned irv=0;irv<size;irv++)x[irv]=rv(lambda);
  return x;
}

Eigen::VectorXd PoissonSignalSim::rvs(const Eigen::VectorXd& lambdas)
{
  unsigned size = lambdas.size();
  Eigen::VectorXd x(size);
  for(unsigned irv=0;irv<size;irv++)x[irv]=rv(lambdas[irv]);
  return x;
}

ExponentialTraceSim::
ExponentialTraceSim(SignalSource* pmt, const Eigen::VectorXd& pmt_pulse, double rate,
    SignalSource* pmt_ap, double rate_ap,
    math::rng::RNG* rng, bool adopt_pmt, bool adopt_pmt_ap, bool adopt_rng):
  pmt_(pmt), tmean_(1.0/rate), pmt_ap_(pmt_ap), tmean_ap_(rate_ap>0 ? (1.0/rate_ap) : 0.0),
  nsample_(pmt_pulse.size()),
  pmt_pulse_fft_(fftw_alloc_real(nsample_)), trace_(fftw_alloc_real(nsample_)),
  rng_((rng==nullptr)?new math::rng::RNG(__PRETTY_FUNCTION__) : rng),
  adopt_pmt_(adopt_pmt), adopt_pmt_ap_(adopt_pmt_ap), adopt_rng_((rng==nullptr)?true:adopt_rng)
{
  int plan_flags = FFTW_ESTIMATE; // proto_planning_enum_to_fftw_flag(config_.fftw_planning());
  trace_plan_fwd_ =
    fftw_plan_r2r_1d(nsample_, trace_, trace_, FFTW_R2HC, plan_flags);
  trace_plan_rev_ =
    fftw_plan_r2r_1d(nsample_, trace_, trace_, FFTW_HC2R, plan_flags);

  std::copy(pmt_pulse.data(), pmt_pulse.data()+nsample_, trace_);
  fftw_execute(trace_plan_fwd_);
  std::copy(trace_, trace_+nsample_, pmt_pulse_fft_);
}

ExponentialTraceSim::~ExponentialTraceSim()
{
  if(adopt_pmt_)delete pmt_;
  if(adopt_pmt_ap_)delete pmt_ap_;
  if(adopt_rng_)delete rng_;
  fftw_destroy_plan(trace_plan_fwd_);
  fftw_destroy_plan(trace_plan_rev_);
  fftw_free(pmt_pulse_fft_);
  fftw_free(trace_);
}

Eigen::VectorXd ExponentialTraceSim::trace()
{
  double t = rng_->exponential(tmean_);
  unsigned it = 0;
  while(it < nsample_) {
    double& sample = trace_[it];
    sample = 0;
    ++it;
    while(it > t) {
      sample += (pmt_ == nullptr) ? 1.0 : pmt_->rv();
      t += rng_->exponential(tmean_);
    }
  }

  if(pmt_ap_ and tmean_ap_>0)
  {
    double t = rng_->exponential(tmean_ap_);
    unsigned it = std::ceil(t);
    while(it < nsample_) {
      trace_[it] += pmt_ap_->rv();
      t += rng_->exponential(tmean_ap_);
      it = std::ceil(t);
    }
  }

  fftw_execute(trace_plan_fwd_);

    hcvec_scale_and_multiply(trace_, trace_, pmt_pulse_fft_, nsample_, 1.0/double(nsample_));
  fftw_execute(trace_plan_rev_);

  Eigen::VectorXd trace(nsample_);
  std::copy(trace_, trace_+nsample_, trace.data());
  return trace;
}
