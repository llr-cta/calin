/*

   calin/simulation/pmt.hpp -- Stephen Fegan -- 2016-03-21

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

// Original header:
// PMTSim.hpp - Stephen Fegan - sfegan@llr.in2p3.fr - May 2014

#pragma once

#include <fftw3.h>

#include <Eigen/Core>
#include <math/rng.hpp>
#include <simulation/pmt.pb.h>

namespace calin { namespace simulation { namespace pmt {

class SignalSource
{
public:
  virtual ~SignalSource();
  virtual double rv() = 0;
  virtual Eigen::VectorXd rvs(unsigned size = 1);
};

// Original Polya PMT simulation with Polya distributon for all stages and an
// option for electron to bypass a stage to generate a low-gain population
class PMTSimPolya: public SignalSource
{
public:
  PMTSimPolya(const calin::ix::simulation::pmt::PMTSimAbbreviatedConfig& config,
	  math::rng::RNG* rng = nullptr);
  virtual ~PMTSimPolya();
  virtual double rv();

  // Slow function to calculate PMF using Prescott (1965).
  calin::ix::simulation::pmt::PMTSimPMF calc_pmf(double precision = 0.0001,
    bool log_progress = false) const;

  double stage_gain(unsigned istage) const;
  double total_gain() const { return total_gain_; }
  double p0() const { return p0_; }

  math::rng::RNG* rng() { return rng_; }
  void set_rng(math::rng::RNG* rng) { delete my_rng_; my_rng_=0; rng_=rng; }

  static calin::ix::simulation::pmt::PMTSimAbbreviatedConfig cta_model_1();
  static calin::ix::simulation::pmt::PMTSimAbbreviatedConfig cta_model_2();
  static calin::ix::simulation::pmt::PMTSimAbbreviatedConfig cta_model_3();

protected:
  calin::ix::simulation::pmt::PMTSimConfig       config_;
  std::vector<double>                            gauss_a_;
  std::vector<double>                            gauss_b_;
  math::rng::RNG*                                rng_;
  math::rng::RNG*                                my_rng_;
  unsigned                                       napprox_limit_;
  double                                         p0_;
  double                                         total_gain_;
};

// Newer Polya PMT simulation with option for multiple populations in the first
// stage but with simple Poisson statistics for all subsequent stages

class PMTSimTwoPopulation: public SignalSource
{
public:
  PMTSimTwoPopulation(const calin::ix::simulation::pmt::PMTSimTwoPopulationConfig& config,
    math::rng::RNG* rng = nullptr, bool use_new_stage_n_algorithm = true);
  virtual ~PMTSimTwoPopulation();
  virtual double rv() override;

  // Slow function to calculate PMF using Prescott (1965).
  calin::ix::simulation::pmt::PMTSimPMF calc_pmf(unsigned nstage = 0,
    double precision = 0.0001, bool log_progress = false) const;

  double total_gain() const { return total_gain_; }
  double p0() const { return p0_; }
  double stage_0_gain() const {
    return config_.stage_0_lo_prob()*config_.stage_0_lo_gain()
      + (1-config_.stage_0_lo_prob())*config_.stage_0_hi_gain();
  }
  double stage_n_gain() const { return stage_n_gain_; }
  std::vector<double> stage1_n_cdf() const { return stage1_n_cdf_; }
  std::vector<double> stage2_n_cdf() const { return stage2_n_cdf_; }
  std::vector<double> stage3_n_cdf() const { return stage3_n_cdf_; }

  math::rng::RNG* rng() { return rng_; }
  void set_rng(math::rng::RNG* rng) { delete my_rng_; my_rng_=0; rng_=rng; }

  // static calin::ix::simulation::pmt::PMTSimAbbreviatedConfig cta_model_4();

protected:
  unsigned stage_n_poisson(unsigned n_in) const;
  unsigned stage_n_new(unsigned n_in) const;
  unsigned stage_n_old(unsigned n_in) const;

  static double stage_p0(double p0_last, double gain, double gain_rms_frac);
  static std::vector<double> stage_pmf(const std::vector<double>& pmf_last,
    double gain, double gain_rms_frac, double precision);
  void recalc_total_gain_and_p0();
  std::string stage_summary(unsigned istage, const std::vector<double>& pk) const;

  calin::ix::simulation::pmt::PMTSimTwoPopulationConfig config_;
  double          p0_                       = 0.0;
  double          total_gain_               = 1.0;
  double          stage_n_gain_             = 1.0;
  double          gamma_a_0_lo_             = 0.0;
  double          gamma_b_0_lo_             = 0.0;
  double          gamma_a_0_hi_             = 0.0;
  double          gamma_b_0_hi_             = 0.0;
  double          gamma_a_n_                = 0.0;
  double          gamma_b_n_                = 0.0;
  math::rng::RNG* rng_                      = nullptr;
  math::rng::RNG* my_rng_                   = nullptr;
  std::vector<double> stage1_n_cdf_;
  std::vector<double> stage2_n_cdf_;
  std::vector<double> stage3_n_cdf_;
  bool            use_new_stage_n_algorithm = false;
};

class PMTSimInvCDF: public SignalSource
{
public:
  PMTSimInvCDF(const calin::ix::simulation::pmt::PMTSimPMF& cmf,
    unsigned npoint = 0, math::rng::RNG* rng = nullptr);
  PMTSimInvCDF(const std::string& filename, math::rng::RNG* rng = nullptr);

  virtual ~PMTSimInvCDF();
  virtual double rv();

  void save_inv_cdf_to_file(const std::string& filename,
			const std::string& comment = "") const;

  math::rng::RNG* rng() { return rng_; }
  void set_rng(math::rng::RNG* rng) { delete my_rng_; my_rng_=0; rng_=rng; }

protected:
  math::rng::RNG*                                rng_;
  math::rng::RNG*                                my_rng_;
  std::vector<std::pair<double,double>>          inv_cdf_;
};

class PMTSimGammaDistribution: public SignalSource
{
public:
  PMTSimGammaDistribution(double alpha, double beta,
    calin::math::rng::RNG* rng = nullptr, bool adopt_rng = false);

  virtual ~PMTSimGammaDistribution();
  virtual double rv();

  calin::math::rng::RNG* rng() { return rng_; }

protected:
  double alpha_;
  double beta_;
  math::rng::RNG* rng_ = nullptr;
  bool adopt_rng_ = false;
};

class MultiPESpectrum: public SignalSource
{
public:
  MultiPESpectrum(SignalSource* pmt,
    const calin::ix::simulation::pmt::MultiPESpectrumConfig& config,
    math::rng::RNG* rng = nullptr);
  virtual ~MultiPESpectrum();

  virtual double rv();
  void rv_xnm(double& x, int& n, double& m);
  void rvs_xnm(VecRef x, IntVecRef n, VecRef m, unsigned size=1);

  double rv_ped();
  Eigen::VectorXd rvs_ped(unsigned size = 1);

  math::rng::RNG* rng() { return rng_; }
  void set_rng(math::rng::RNG* rng) { delete my_rng_; my_rng_=0; rng_=rng; }

private:
  SignalSource*                                     pmt_;
  math::rng::RNG*                                   rng_;
  math::rng::RNG*                                   my_rng_;
  calin::ix::simulation::pmt::MultiPESpectrumConfig config_;
  double                                            signal_gamma_a_ = 0;
  double                                            signal_gamma_b_ = 0;
};

class PoissonSignalSim
{
public:
  PoissonSignalSim(SignalSource* pmt,
    const calin::ix::simulation::pmt::PoissonSignalSimConfig& config,
    math::rng::RNG* rng = nullptr);
  virtual ~PoissonSignalSim();

  double rv(double lambda);
  Eigen::VectorXd rvs(double lambda, unsigned size);
  Eigen::VectorXd rvs(const Eigen::VectorXd& lambdas);

  math::rng::RNG* rng() { return rng_; }
  void set_rng(math::rng::RNG* rng) { delete my_rng_; my_rng_=0; rng_=rng; }

  double pedestal_mean() const { return config_.pedestal_mean(); }
  double pedestal_rms() const { return config_.pedestal_rms(); }

private:
  SignalSource*                                      pmt_;
  math::rng::RNG*                                    rng_;
  math::rng::RNG*                                    my_rng_;
  calin::ix::simulation::pmt::PoissonSignalSimConfig config_;
};

class ExponentialTraceSim
{
public:
  ExponentialTraceSim(const Eigen::VectorXd& pmt_pulse, double rate,
      math::rng::RNG* rng = nullptr, bool adopt_rng = false) :
    ExponentialTraceSim(nullptr, pmt_pulse, rate, nullptr, 0, rng, false, false, adopt_rng) { }

  ExponentialTraceSim(SignalSource* pmt, const Eigen::VectorXd& pmt_pulse, double rate,
      math::rng::RNG* rng = nullptr, bool adopt_pmt = false, bool adopt_rng = false):
    ExponentialTraceSim(pmt, pmt_pulse, rate, nullptr, 0, rng, adopt_pmt, false, adopt_rng) { }

  ExponentialTraceSim(SignalSource* pmt, const Eigen::VectorXd& pmt_pulse, double rate,
    SignalSource* pmt_ap, double rate_ap, math::rng::RNG* rng = nullptr,
    bool adopt_pmt = false, bool adopt_pmt_ap = false, bool adopt_rng = false);

  ~ExponentialTraceSim();

  Eigen::VectorXd trace();

  double rate() const { return 1.0/tmean_; }
  void set_rate(double rate) { tmean_ = 1.0/rate; }
private:
  SignalSource* pmt_ = nullptr;
  double tmean_ = 0;
  SignalSource* pmt_ap_;
  double tmean_ap_ = 0;
  unsigned nsample_;
  double* pmt_pulse_fft_ = nullptr;
  double* trace_ = nullptr;
  fftw_plan trace_plan_fwd_;
  fftw_plan trace_plan_rev_;
  math::rng::RNG* rng_ = nullptr;
  bool adopt_pmt_ = false;
  bool adopt_pmt_ap_ = false;
  bool adopt_rng_ = false;
};

} } } // namespace calin::simulation::pmt
