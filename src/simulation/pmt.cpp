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

void SignalSource::rvs(std::vector<double>& n, unsigned size)
{
  if(size==0)size = n.size();
  else n.resize(size);
  for(unsigned irv=0;irv<size;irv++)n[irv] = rv();
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
  config.set_stage_0_gain(9.1);
  config.set_stage_0_prob_skip(0.1);
  return config;
}


PMTSimPolya::
PMTSimPolya(const calin::ix::simulation::pmt::PMTSimAbbreviatedConfig& config,
    math::rng::RNG* rng):
  SignalSource(),
  rng_(rng), my_rng_(0), napprox_limit_(50)
{
  if(config.num_stage()==0)
    throw std::runtime_error("PMTSimPolya: number of stages must be positive");

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
  for(unsigned istage=0; istage<config_.stage_size(); istage++)
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
  for(unsigned istage=0; istage<config_.stage_size(); istage++)
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

Eigen::VectorXd PMTSimPolya::
calc_pmf(bool log_progress, double precision) const
{
  std::vector<double> pk(2);
  pk[0] = 0;
  pk[1] = 1;
  for(unsigned ik=0;ik<config_.stage_size();ik++)
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

    if(log_progress)LOG(INFO) << ik << ' ' << psum;
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

    if(log_progress)
    {
      auto log_info = LOG(INFO);
      log_info << "-- " << psum;
      for(unsigned ix=0;ix<std::min(unsigned(config_.stage_size()+1),
          unsigned(pk.size()));ix++)
        log_info << ' ' << pk[ix];
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
  }

  return std_to_eigenvec(pk);
}

#if 0
// =============================================================================
//
// PMTSimInvCDF
//
// =============================================================================

PMTSimInvCDF::PMTSimInvCDF(RandomNumbers* rng):
  SignalSource(), m_rng(rng), m_my_rng(0)
{
  if(rng==0)
    m_rng = m_my_rng = new RandomNumbers(RandomNumbers::defaultFilename());
}

PMTSimInvCDF::~PMTSimInvCDF()
{
  delete m_my_rng;
}

double PMTSimInvCDF::rv()
{
  return m_rng->InverseCDF(m_inv_cdf);
}

void PMTSimInvCDF::
setInvCDFFromPMTPMF(const std::vector<double>& y, double gain,
		    bool suppress_zero, unsigned npoint)
{
  double iy0 = 0;
  if(suppress_zero)iy0 = 1;
  m_inv_cdf.resize(y.size() + 1 - iy0);
  double ysum = 0;
  for(unsigned iy=iy0;iy<y.size();iy++)ysum += y[iy];
  double ycumsum = 0;
  m_inv_cdf[0] = std::make_pair<double,double>((double(iy0)-0.5)/gain,0.0);
  for(unsigned iy=iy0;iy<y.size();iy++)
    {
      ycumsum += y[iy];
      m_inv_cdf[iy + 1 - iy0] =
	std::make_pair<double,double>((double(iy)+0.5)/gain,ycumsum/ysum);
    }
#if 0
  for(unsigned iy=0;iy<10;iy++)
    std::cout << m_inv_cdf[iy].first << ' ' << m_inv_cdf[iy].second << '\n';
#endif
  m_rng->GenerateInverseCDF(m_inv_cdf, npoint);
}

void PMTSimInvCDF::
setInvCDFFromFile(const std::string& filename)
{
  std::ifstream s(filename.c_str());
  if(!s)throw std::string("PMTSimInvCDF::setInvCDFFromFile: Could not open file: ")+filename;
  m_inv_cdf.clear();
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
	m_inv_cdf.push_back(std::make_pair<double,double>(x,y));
    next_line:
      std::getline(s,line);
    }
}

void PMTSimInvCDF::
saveInvCDFToFile(const std::string& filename,
		 const std::string& comment) const
{
  std::ofstream s(filename.c_str());
  if(!s)throw std::string("PMTSimInvCDF::saveInvCDFToFile: Could not open file: ")+filename;
  time_t rawtime;
  time (&rawtime);
  s << "# Inverse CDF written " << ctime(&rawtime);
  if(!comment.empty())s << "# " << comment << '\n';
  for(unsigned iy=0;iy<m_inv_cdf.size();iy++)
    s << VERITAS::VSDataConverter::toString(m_inv_cdf[iy].first) << ' '
      << VERITAS::VSDataConverter::toString(m_inv_cdf[iy].second) <<  '\n';
}

// =============================================================================
//
// MultiPESpectrum
//
// =============================================================================

MultiPESpectrum::
MultiPESpectrum(SignalSource* pmt,
		double signal_mean, double pedestal_rms,
		double signal_rms_frac,  double pedestal_mean,
		double quantization, RandomNumbers* rng):
  SignalSource(),  m_pmt(pmt), m_rng(rng), m_my_rng(0),
  m_signal_mean(signal_mean), m_signal_rms_frac(signal_rms_frac),
  m_signal_gamma_a(), m_signal_gamma_b(),
  m_pedestal_rms(pedestal_rms), m_pedestal_mean(pedestal_mean),
  m_quantization(quantization)
{
  if(rng==0)
    m_rng = m_my_rng = new RandomNumbers(RandomNumbers::defaultFilename());
  if(signal_rms_frac>0)
    {
      m_signal_gamma_a = 1.0/SQR(signal_rms_frac);
      m_signal_gamma_b = m_signal_gamma_a/signal_mean;
    }
}

MultiPESpectrum::~MultiPESpectrum()
{
  delete m_my_rng;
}

double MultiPESpectrum::rv()
{
  double x;
  double npe_mean;
  unsigned npe;
  rv_xnm(x,npe,npe_mean);
  return x;
}

void MultiPESpectrum::rv_xnm(double& x, unsigned& n, double& m)
{
  x = m_rng->Normal()*m_pedestal_rms + m_pedestal_mean;
  if(m_signal_rms_frac>0)
    m = m_rng->Gamma(m_signal_gamma_a,m_signal_gamma_b);
  else
    m = m_signal_mean;
  n = m_rng->Poisson(m);
  for(unsigned ipe=0;ipe<n;ipe++)x += m_pmt->rv();
}

void MultiPESpectrum::
rvs_xnm(std::vector<double>& x, std::vector<unsigned>& n,
	std::vector<double>& m, unsigned size)
{
  if(size==0)size = x.size();
  else x.resize(size);
  n.resize(size);
  m.resize(size);
  for(unsigned irv=0;irv<size;irv++)rv_xnm(x[irv],n[irv],m[irv]);
}

double MultiPESpectrum::rv_ped()
{
  return m_rng->Normal()*m_pedestal_rms + m_pedestal_mean;
}

void MultiPESpectrum::rvs_ped(std::vector<double>& x, unsigned size)
{
  if(size==0)size = x.size();
  else x.resize(size);
  for(unsigned irv=0;irv<size;irv++)x[irv]=rv_ped();
}
#endif
