/*

   calin/simulation/pmt.cpp -- Stephen Fegan -- 2016-03-21

   Class to genereate random deviates from a PMT spectrum

   Copyright 2014, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <algorithm>
#include <cmath>
#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cctype>

#include <calin_global_definitions.hpp>
#include <simulation/pmt.hpp>
#include <math/special.hpp>
#include <io/log.hpp>

using namespace calin::simulation::pmt;
using calin::math::special::SQR;
using calin::math::rng::RNG;
using namespace calin::io::log;

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

  if(rng==nullptr)rng_ = my_rng_ = new RNG();
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

// =============================================================================
//
// PMTSimInvCDF
//
// =============================================================================

PMTSimInvCDF::PMTSimInvCDF(const calin::ix::simulation::pmt::PMTSimPMF& cmf,
    unsigned npoint, math::rng::RNG* rng):
  SignalSource(), rng_(rng), my_rng_(nullptr)
{
  if(rng==nullptr)rng_ = my_rng_ = new math::rng::RNG();

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
  if(rng==nullptr)rng_ = my_rng_ = new math::rng::RNG();
  std::ifstream s(filename.c_str());
  if(!s)throw std::runtime_error(
    std::string("PMTSimInvCDF: Could not open file: ")+filename);
  std::string line;
  std::getline(s,line);
  while(s)
  {
    unsigned ic=0;
    while(ic<line.size() && isspace(line[ic]))ic++;
    if(ic==line.size())goto next_line;
    if(line[ic] == '#')goto next_line;
    double x,y;
    if(std::istringstream(line) >> x >> y)
    inv_cdf_.emplace_back(x,y);
  next_line:
    std::getline(s,line);
  }
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
  time_t rawtime;
  time (&rawtime);
  s << "# Inverse CDF written " << ctime(&rawtime);
  if(!comment.empty())s << "# " << comment << '\n';
  for(unsigned iy=0;iy<inv_cdf_.size();iy++)
    s << std::to_string(inv_cdf_[iy].first) << ' '
      << std::to_string(inv_cdf_[iy].second) <<  '\n';
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
  if(rng==nullptr)rng_ = my_rng_ = new math::rng::RNG();
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
