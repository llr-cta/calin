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

#include <cmath>
#include <cassert>
#include <simulation/pmt.hpp>
#include <math/special.hpp>

using namespace calin::simulation::pmt;
using calin::math::special::SQR

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

ix::simulation::pmt::PMTSimAbbreviatedConfig PMTSimPolya::default_config()
{
  ix::simulation::pmt::PMTSimAbbreviatedConfig config;
  config.set_num_stages(7);
  config.set_total_gain(4E4);
  config.set_stage_1_gain(8);
  config.set_suppress_zero(true);
  return config;
}

PMTSimPolya::
PMTSimPolya(unsigned nstage, double total_gain,
	    double s1_gain, double s1_gain_rms_frac, double s1_pskip,
	    double sn_gain_rms_frac, double sn_pskip,
	    bool suppress_zero, bool signal_in_pe,
	    RandomNumbers* rng):
  SignalSource(),
  m_stages(nstage), m_rng(rng), m_my_rng(0), m_napprox_limit(50),
  m_suppress_zero(suppress_zero), m_signal_in_pe(signal_in_pe)
{
  if(nstage==0)throw std::string("PMTSimPolya: nstages must be positive");

  if(rng==0)
    m_rng = m_my_rng = new RandomNumbers(RandomNumbers::defaultFilename());

  m_stages[0].gain_mean                  = s1_gain;
#if 0
  // Must decide what definition of s1_gain should be - either gain of
  // the dynode, or the average of the number of electrons entering
  // second stage (including skips)
  m_stages[0].gain_mean                  = (s1_gain-s1_pskip)/(1.0-s1_pskip);
#endif
  m_stages[0].pskip                      = s1_pskip;
  if(s1_gain_rms_frac>0)
    {
      m_stages[0].gain_rms_frac          = s1_gain_rms_frac;
      m_stages[0].gauss_a                = 1.0/SQR(s1_gain_rms_frac);
      m_stages[0].gauss_b                = m_stages[0].gauss_a/m_stages[0].gain_mean;
    }

  double s2_mean_n0 = s1_pskip + m_stages[0].gain_mean*(1.0-s1_pskip);
  double sn_gain = std::pow(total_gain/s2_mean_n0,1.0/double(nstage-1));
  sn_gain = (sn_gain - sn_pskip)/(1.0 - sn_pskip);
  for(unsigned istage=1; istage<nstage; istage++)
    {
      m_stages[istage].gain_mean         = sn_gain;
      m_stages[istage].pskip             = sn_pskip;
      if(sn_gain_rms_frac>0)
	{
	  m_stages[istage].gain_rms_frac = sn_gain_rms_frac;
	  m_stages[istage].gauss_a       = 1.0/SQR(sn_gain_rms_frac);
	  m_stages[istage].gauss_b       = m_stages[istage].gauss_a/sn_gain;
	}
    }

  m_total_gain = 1;
  m_p0         = 0;
  for(unsigned istage=0; istage<m_stages.size(); istage++)
    {
      m_total_gain *= stageGain(istage);

      // This little mess of code implements the calculation of p0
      // from Prescott (1965) to account for the suppressed zeros in
      // the calculation of the total gain
      const PMTSimStage& stage(m_stages[m_stages.size()-istage-1]);
      double p0_last = m_p0;
      double mu      = stage.gain_mean;
      double b       = SQR(stage.gain_rms_frac);
      double bmu     = b*mu;
      double C1      = 1.0 + bmu*(1.0-p0_last);
      if(b == 0)
	m_p0 = std::exp(mu*(p0_last-1.0));
      else
	m_p0 = pow(C1,-1.0/b);
      m_p0 = m_p0*(1.0-stage.pskip) + p0_last*stage.pskip;
    }
  if(m_suppress_zero)
    m_total_gain /= 1.0-m_p0;
}

double PMTSimPolya::stageGain(unsigned istage) const
{
  const PMTSimStage& stage(m_stages[istage]);
  double mu      = stage.gain_mean;
  return mu*(1.0-stage.pskip) + stage.pskip;
}

PMTSimPolya::~PMTSimPolya()
{
  delete m_my_rng;
}

double PMTSimPolya::rv()
{
  unsigned n;
 do_over:
  n = 1;
  for(unsigned istage=0; istage<m_stages.size(); istage++)
    {
      const PMTSimStage& stage(m_stages[istage]);
      unsigned n_in = n;

      if(stage.pskip>0)
	{
	  n = m_rng->Binomial(stage.pskip, n_in);
	  n_in -= n;
	  if(n_in == 0)continue;
	}
      else n = 0;

      double nmean = 0;
      if(stage.gain_rms_frac>0)
	{
	  if(n_in > m_napprox_limit)
	    nmean = m_rng->Gamma(double(n_in)*stage.gauss_a,
				 stage.gauss_b);
	  else
	    for(unsigned in_in = 0; in_in<n_in;in_in++)
	      nmean += m_rng->Gamma(stage.gauss_a,stage.gauss_b);
	}
      else nmean = double(n_in)*stage.gain_mean;
      n += m_rng->Poisson(nmean);

      if(n==0)
	{
	  if(m_suppress_zero)goto do_over;
	  else return 0.0;
	}
    }
  if(m_signal_in_pe)return double(n)/m_total_gain;
  else return double(n);
}

// Implements the calculation of the PDF for the SES from Prescott
// (1965). The calculation is O(n^2) in the number of total entries in
// the final SES, which is itself some factor (4-5) times the total
// gain. This calculation could also be done with FFTs which may be
// faster. Interestingly the FFT was "rediscovered" in ~1965, so
// Prescott may not have known about it.

void PMTSimPolya::
pmf(std::vector<double>& pk, double precision) const
{
  pk.resize(2);
  pk[0] = 0;
  pk[1] = 1;
  for(unsigned ik=0;ik<m_stages.size();ik++)
    {
      const PMTSimStage& stage(m_stages[m_stages.size()-ik-1]);

      std::vector<double> pkmo(pk);

      double mu      = stage.gain_mean;
      double b       = SQR(stage.gain_rms_frac);
      double bmu     = b*mu;
      double bmo     = b-1.0;
      double C1      = 1.0 + bmu*(1.0-pkmo[0]);
      double pcutoff = 1.0 - precision*((ik==m_stages.size()-1)?0.1:1.0);
      double psum;

      pk.resize(1);
      pk.reserve(pkmo.size()*5*int(ceil(mu)));
      if(b == 0)
	pk[0] = std::exp(mu*(pkmo[0]-1.0));
      else
	pk[0] = pow(C1,-1.0/b);
      psum = pk[0];

      std::cout << ik << ' ' << psum << '\n';
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
      std::cout << "-- " << psum;
      for(unsigned ix=0;ix<std::min(m_stages.size()+1,pk.size());ix++)
	std::cout << ' ' << pk[ix];
      std::cout << '\n';

      if(stage.pskip>0)
	{
	  for(unsigned ix=0;ix<pk.size();ix++)
	    pk[ix] *= (1.0-stage.pskip)/psum;
	  for(unsigned ix=0;ix<std::min(pkmo.size(),pk.size());ix++)
	    pk[ix] += stage.pskip*pkmo[ix];
	}
      else
	{
	  for(unsigned ix=0;ix<pk.size();ix++)
	    pk[ix] /= psum;
	}
    }
}

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
