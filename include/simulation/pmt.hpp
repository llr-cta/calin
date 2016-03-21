/*

   calin/simulation/pmt.hpp -- Stephen Fegan -- 2016-03-21

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

// Original header:
// PMTSim.hpp - Stephen Fegan - sfegan@llr.in2p3.fr - May 2014

#pragma once

#include <math/rng.hpp>
#include <simulation/pmt.pb.h>

namespace calin { namespace simulation { namespace pmt {

class SignalSource
{
public:
  virtual ~SignalSource();
  virtual double rv() = 0;
  virtual void rvs(std::vector<double>& n, unsigned size = 1);
};

class PMTSimPolya: public SignalSource
{
public:
  PMTSimPolya(const ix::simulation::pmt::PMTSimAbbreviatedConfig& config,
	  math::rng::RNG* rng = nullptr);
  virtual ~PMTSimPolya();
  virtual double rv();
  //virtual void rvs(std::vector<double>& n, unsigned size = 1);
  void calc_pmf(std::vector<double>& p, double precision = 0.0001) const;
  double stage_gain(unsigned istage) const;
  double total_gain() const { return m_total_gain; }
  double p0() const { return m_p0; }
  RandomNumbers* rng() { return m_rng; }
  void set_rng(RandomNumbers* rng) { delete m_my_rng; m_my_rng=0; m_rng=rng; }

  ix::simulation::pmt::PMTSimAbbreviatedConfig default_config();
private:
  std::vector<ix::simulation::pmt::PMTSimConfig> config_;
  RandomNumbers*                                 rng_;
  RandomNumbers*                                 my_rng_;
  unsigned                                       napprox_limit_;
  double                                         p0_;
  double                                         total_gain_;
};

#if 0
class PMTSimInvCDF: public SignalSource
{
public:
  PMTSimInvCDF(RandomNumbers* rng = 0);
  virtual ~PMTSimInvCDF();
  virtual double rv();
  void setInvCDFFromPMTPMF(const std::vector<double>& y, double gain = 1.0,
			   bool suppress_zero = true, unsigned npoint = 0);
  void setInvCDFFromFile(const std::string& filename);
  void saveInvCDFToFile(const std::string& filename,
			const std::string& comment = "") const;
  RandomNumbers* rng() { return m_rng; }
  void setRNG(RandomNumbers* rng) { delete m_my_rng; m_my_rng=0; m_rng=rng; }
private:
  RandomNumbers*           m_rng;
  RandomNumbers*           m_my_rng;
  std::vector<std::pair<double,double> > m_inv_cdf;
};

class MultiPESpectrum: public SignalSource
{
public:
  MultiPESpectrum(SignalSource* pmt,
		  double signal_mean, double pedestal_rms,
		  double signal_rms_frac = 0.0, double pedestal_mean = 0.0,
		  double quantization = 0.0, RandomNumbers* rng = 0);
  virtual ~MultiPESpectrum();
  virtual double rv();
  void rv_xnm(double& x, unsigned& n, double& m);
  void rvs_xnm(std::vector<double>& x, std::vector<unsigned>& n,
	       std::vector<double>& m, unsigned size=1);
  double rv_ped();
  void rvs_ped(std::vector<double>& x, unsigned size = 1);
  RandomNumbers* rng() { return m_rng; }
  void setRNG(RandomNumbers* rng) { delete m_my_rng; m_my_rng=0; m_rng=rng; }
private:
  SignalSource*            m_pmt;
  RandomNumbers*           m_rng;
  RandomNumbers*           m_my_rng;
  double                   m_signal_mean;
  double                   m_signal_rms_frac;
  double                   m_signal_gamma_a;
  double                   m_signal_gamma_b;
  double                   m_pedestal_rms;
  double                   m_pedestal_mean;
  double                   m_quantization;
};
#endif
} } } // namespace calin::simulation::pmt
