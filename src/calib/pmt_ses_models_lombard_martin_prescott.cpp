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
    if(config.apply_downsampling()) {
      if(config.downsampling_num_stage() >= config.num_stage())
        throw std::runtime_error("PMTSimTwoPopulation: number of stages for down-sampling must be smaller than total.");
      if(config.downsampling_factor() <= 1)
        throw std::runtime_error("PMTSimTwoPopulation: down-sampling factor must be greater than one.");
    }
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
    std::vector<double> pmf;
    if(config.stage_0_lo_pmf_size()) {
      pmf.resize(config.stage_0_lo_pmf_size());
      std::copy(config.stage_0_lo_pmf().begin(), config.stage_0_lo_pmf().end(), pmf.begin());
    } else if(config.stage_0_lo_half_gaussian_pmf()) {
      pmf = half_gaussian_pmf(config.stage_0_lo_gain(), precision);
    } else {
      pmf = polya_pmf(config.stage_0_lo_gain(), config.stage_0_lo_gain_rms_frac(), precision);
    }

    if(pmf.size() > stage_0_pmf_.size()) {
      stage_0_pmf_.resize(pmf.size(), 0.0);
    }

    calin::math::accumulator::KahanAccumulator psum;
    for(unsigned i=0; i<pmf.size(); ++i) {
      double pi = pmf[i];
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

  for(int istage=0; istage<int(config_.num_stage())-2; ++istage) {
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

// Implements the calculation of the PDF for the multi-stage Polya distribution
// using the algoirithm from Prescott (1965). The calculation is O(n^2) in the
// number of total entries in the final PDF, which is itself some factor (4-5)
// times the total gain.

std::vector<double> LombardMartinPrescottPMTModel::
multi_stage_polya_pmf(unsigned nstage, double mean, double rms_frac, unsigned rebinning, double precision)
{
  std::vector<double> pk;
  std::vector<double> pkmo;
  pk.reserve(5*int(pow(mean,nstage)));
  pkmo.reserve(5*int(pow(mean,nstage)));

  pk.push_back(0.0);
  pk.push_back(1.0);

  const double mu      = mean;
  const double b       = SQR(rms_frac);
  const double bmu     = b*mu;
  const double bmo     = b-1.0;

  for(unsigned ik=0;ik<nstage;ik++)
  {
    std::swap(pk, pkmo);

    double C1      = 1.0 + bmu*(1.0-pkmo[0]);
    double pcutoff = 1.0 - precision*((ik==nstage-1)?0.1:1.0);
    double psum;

    pk.resize(1);
    if(b == 0)
      pk[0] = std::exp(mu*(pkmo[0]-1.0));
    else
      pk[0] = pow(C1,-1.0/b);
    psum = pk[0];

    for(unsigned ix=1;psum<pcutoff;ix++)
  	{
  	  double pkx = 0;
  	  for(unsigned ii=ix-std::min(ix,unsigned(pkmo.size()-1));ii<ix;ii++) {
    		pkx += pk[ii]*pkmo[ix-ii]*(double(ix)+double(ii)*bmo);
	    }
  	  pkx *= mu/double(ix)/C1;
  	  psum += pkx;
  	  pk.push_back(pkx);
  	  if(pkx/psum < 1e-10)break;
  	}
  }

  rebin_pmf(pk, rebinning);

  return pk;
}

std::vector<double> LombardMartinPrescottPMTModel::
multi_stage_pmf(Tableau& tableau, unsigned nstage, const std::vector<double>& pmf,
  unsigned rebinning, double precision, bool suppress_wraparound_warning,
  unsigned* wraparound_warning_count)
{
  // Set up the convolutions
  std::vector<const std::vector<double>*> all_stages(nstage, &pmf);

  // Do the convolutions
  hcvec_multi_stage_polynomial(tableau.dft.get(), tableau.basis_dft.get(),
    all_stages, tableau.npoint);

  // Execute the inverse DFT
  fftw_execute(tableau.dft_to_pmf_bwd_plan.get());

  // Find the positive part of the PMF
  precision *= nstage*tableau.npoint;
  const double pcutoff1 = tableau.npoint - precision;
  const double pcutoff2 = 0.99 * tableau.npoint;

  calin::math::accumulator::KahanAccumulator psum;
  unsigned ipoint;
  for(ipoint=0;ipoint<tableau.npoint && psum.total()<pcutoff1;ipoint++) {
    double p = tableau.pmf.get()[ipoint];
    psum.accumulate(p);
    if(psum.total() > pcutoff2 and p < precision)break;
  }

  if(!suppress_wraparound_warning) {
    if(ipoint*100 > tableau.npoint*95) {
      LOG(WARNING) << "multi_stage_pmf: calculated PMF length close to buffer size; possible wrap around.";
      if(wraparound_warning_count) {
        ++(*wraparound_warning_count);
      }
    }
  }

  // Rebin the PMF is necessary
  unsigned npmf = rebin_pmf(tableau.pmf.get(), ipoint, rebinning);

  // Copy the spectrum to the output vector, normalizing for FFTW conventions
  const double norm = 1.0/double(tableau.npoint);
  std::vector<double> pk;
  pk.reserve(npmf);
  for(unsigned ipmf=0; ipmf<npmf; ipmf++) {
    pk.push_back(norm * tableau.pmf.get()[ipmf]);
  }

  return pk;
}

std::vector<double> LombardMartinPrescottPMTModel::
multi_stage_pmf(unsigned npoint, unsigned nstage, const std::vector<double>& pmf,
  unsigned rebinning, double precision, calin::ix::math::fftw_util::FFTWPlanningRigor fftw_rigor)
{
  if(npoint==0) {
    throw std::runtime_error("LombardMartinPrescottPMTModel::multi_stage_pmf: npoint must be positive");
  }
  Tableau tableau(npoint, fftw_rigor);
  return multi_stage_pmf(tableau, nstage, pmf, rebinning, precision);
}

std::vector<double> LombardMartinPrescottPMTModel::
half_gaussian_pmf(double mean, double precision)
{
  const double C1      = 0.5*M_2_SQRTPI/(std::abs(mean));

  const double pcutoff = 1.0 - precision;

  std::vector<double> pk;
  pk.reserve(5*int(ceil(mean)));

  double p_last = 0;
  for(unsigned ix=0;p_last<pcutoff;ix++)
  {
    double p = std::erf((double(ix)+0.5)*C1);
    pk.push_back(p-p_last);
    p_last = p;
  }
  return pk;
}

unsigned LombardMartinPrescottPMTModel::
rebin_pmf(double* pmf, unsigned npmf, unsigned binning)
{
  if(npmf==0 or binning <= 1) {
    return npmf;
  } else if(binning % 2 == 1) {
    // odd binning - all PMF entries go to one bin only
    const unsigned newbin = (binning + 1)/2;
    unsigned obin = 0;
    for(unsigned ibin=1; ibin<npmf; ++ibin) {
      if(ibin % binning == newbin) {
        pmf[++obin] = pmf[ibin];
      } else {
        pmf[obin] += pmf[ibin];
      }
    }
    return obin+1;
  } else {
    const unsigned newbin = binning/2;
    unsigned obin = 0;
    for(unsigned ibin=1; ibin<npmf; ++ibin) {
      if(ibin % binning == newbin) {
        double p = 0.5*pmf[ibin];
        pmf[obin] += p;
        pmf[++obin] = p;
      } else {
        pmf[obin] += pmf[ibin];
      }
    }
    return obin+1;
  }
}

void LombardMartinPrescottPMTModel::
rebin_pmf(std::vector<double>& pmf, unsigned binning)
{
  unsigned new_size = rebin_pmf(pmf.data(), pmf.size(), binning);
  pmf.resize(new_size);
  return;
}

std::vector<double> LombardMartinPrescottPMTModel::
stage_n_pmf_downsampled(Tableau& tableau) const
{
  if(not config_.apply_downsampling())return {};

  if(config_.downsampling_num_stage()<4 or
      (config_.downsampling_num_stage()==4 and stage_n_gain_<4.5) or
      (config_.downsampling_num_stage()==5 and stage_n_gain_<3.0)) {
    return multi_stage_polya_pmf(config_.downsampling_num_stage(),
      stage_n_gain_, config_.stage_n_gain_rms_frac(),
      config_.downsampling_factor(), precision_);
  } else {
    return multi_stage_pmf(tableau, config_.downsampling_num_stage(), stage_n_pmf_,
      config_.downsampling_factor(), precision_,
      suppress_wraparound_warning_,  &num_wraparound_warnings_sent_);
  }
}

void LombardMartinPrescottPMTModel::calc_spectrum(
  Tableau& tableau, const std::vector<double>* pe_spec, const double* ped, bool ped_is_fft)
{
  int plan_flags = proto_planning_enum_to_fftw_flag(fftw_rigor_);

  // Set up the convolutions
  std::vector<const std::vector<double>*> all_stages;
  all_stages.reserve(config_.num_stage()+1);

  std::vector<double> _stage_n_pmf_downsampled;
  if(config_.apply_downsampling()) {
    // Do one convolution with downsampled multi-stage PMF followed by ...
    _stage_n_pmf_downsampled = stage_n_pmf_downsampled(tableau);
    all_stages.push_back(&_stage_n_pmf_downsampled);
    // num_stage()-downsampling_num_stage()-1 convolutions with the "n-stage" PMF
    for(unsigned istage=config_.downsampling_num_stage(); istage<config_.num_stage()-1; istage++) {
      all_stages.push_back(&stage_n_pmf_);
    }
  } else {
    // Do num_stage()-1 convolutions with the "n-stage" PMF
    for(unsigned istage=0; istage<config_.num_stage()-1; istage++) {
      all_stages.push_back(&stage_n_pmf_);
    }
  }

  // Then convolve with stage 0 PMF (zero suppress adjusted if necessary)
  all_stages.push_back(&stage_0_pmf_zsa_);

  // Finally convolve with the PE spectrum PMF if given
  if(pe_spec != nullptr) {
    all_stages.push_back(pe_spec);
  }

  // Do the convolutions
  hcvec_multi_stage_polynomial(tableau.dft.get(), tableau.basis_dft.get(),
    all_stages, tableau.npoint);

  // If we have a pedestal then apply it
  if(ped != nullptr) {
    if(ped_is_fft) {
      hcvec_scale_and_multiply(tableau.dft.get(), tableau.dft.get(),
        ped, tableau.npoint);
    } else {
      uptr_fftw_plan fwd_plan = {
        fftw_plan_r2r_1d(tableau.npoint, (double*)ped, tableau.pmf.get(),
          FFTW_R2HC, plan_flags), fftw_destroy_plan };
      assert(fwd_plan);
      fftw_execute(fwd_plan.get());
      hcvec_scale_and_multiply(tableau.dft.get(), tableau.dft.get(),
        tableau.pmf.get(), tableau.npoint);
    }
  }

  // Prepare the inverse DFT
  fftw_execute(tableau.dft_to_pmf_bwd_plan.get());

  // Copy the spectrum to the output vector, normalizing for FFTW conventions
  double norm = 1.0/double(tableau.npoint);
  std::transform(tableau.pmf.get(), tableau.pmf.get()+tableau.npoint, tableau.pmf.get(),
    [norm](double x) { return x*norm; });
}

void LombardMartinPrescottPMTModel::calc_ses(Tableau& tableau)
{
  calc_spectrum(tableau);
}

Eigen::VectorXd LombardMartinPrescottPMTModel::calc_ses(unsigned npoint)
{
  if(npoint==0) {
    npoint = calin::math::special::round_up_power_of_two(total_gain_ * 4);
    if(config_.apply_downsampling())npoint /= config_.downsampling_factor();
  }
  Tableau tableau(npoint, fftw_rigor_);
  calc_spectrum(tableau);
  return tableau.pmf_as_vec();
}

void LombardMartinPrescottPMTModel::calc_mes(
  Tableau& tableau, double mean, double rms_frac,
  const double* ped, bool ped_is_fft)
{
  std::vector<double> pe_pmf = polya_pmf(mean, rms_frac, precision_);
  calc_spectrum(tableau, &pe_pmf, ped, ped_is_fft);
}

void LombardMartinPrescottPMTModel::calc_mes(
  Tableau& tableau, const std::vector<double>& pe_pmf,
  const double* ped, bool ped_is_fft)
{
  calc_spectrum(tableau, &pe_pmf, ped, ped_is_fft);
}

Eigen::VectorXd LombardMartinPrescottPMTModel::calc_mes(
  double mean, double rms_frac, unsigned npoint, const double* ped, bool ped_is_fft)
{
  if(npoint==0) {
    throw std::runtime_error("LombardMartinPrescottPMTModel::calc_mes: npoint must be positive");
  }
  std::vector<double> pe_pmf = polya_pmf(mean, rms_frac, precision_);
  Tableau tableau(npoint, fftw_rigor_);
  calc_spectrum(tableau, &pe_pmf, ped, ped_is_fft);
  return tableau.pmf_as_vec();
}

Eigen::VectorXd LombardMartinPrescottPMTModel::calc_mes(
  const std::vector<double>& pe_pmf, unsigned npoint, const double* ped, bool ped_is_fft)
{
  if(npoint==0) {
    throw std::runtime_error("LombardMartinPrescottPMTModel::calc_mes: npoint must be positive");
  }
  Tableau tableau(npoint, fftw_rigor_);
  calc_spectrum(tableau, &pe_pmf, ped, ped_is_fft);
  return tableau.pmf_as_vec();
}

Eigen::VectorXd LombardMartinPrescottPMTModel::calc_mes(
  double mean, double rms_frac, unsigned npoint)
{
  if(npoint==0) {
    npoint = calin::math::special::round_up_power_of_two(total_gain_ * 4 * mean);
    if(config_.apply_downsampling())npoint /= config_.downsampling_factor();
  }
  return calc_mes(mean, rms_frac, npoint, nullptr, false);
}

Eigen::VectorXd LombardMartinPrescottPMTModel::calc_mes(
  const std::vector<double>& pe_pmf, unsigned npoint)
{
  if(npoint==0) {
    double mean = 0;
    for(unsigned ipe=0; ipe<pe_pmf.size(); ipe++) {
      mean += pe_pmf[ipe]*ipe;
    }
    npoint = calin::math::special::round_up_power_of_two(total_gain_ * 4 * mean);
    if(config_.apply_downsampling())npoint /= config_.downsampling_factor();
  }
  return calc_mes(pe_pmf, npoint, nullptr, false);
}

Eigen::VectorXd LombardMartinPrescottPMTModel::calc_mes(
  double mean, double rms_frac, const Eigen::VectorXd& ped, bool ped_is_fft)
{
  // Assign a buffer to hold (properly aligned) pedestal data
  uptr_fftw_data ped_copy { fftw_alloc_real(ped.size()), fftw_free };
  assert(ped_copy);
  std::copy(ped.data(), ped.data()+ped.size(), ped_copy.get());
  return calc_mes(mean, rms_frac, ped.size(), ped_copy.get(), ped_is_fft);
}

Eigen::VectorXd LombardMartinPrescottPMTModel::calc_mes(
  const std::vector<double>& pe_pmf, const Eigen::VectorXd& ped, bool ped_is_fft)
{
  // Assign a buffer to hold (properly aligned) pedestal data
  uptr_fftw_data ped_copy { fftw_alloc_real(ped.size()), fftw_free };
  assert(ped_copy);
  std::copy(ped.data(), ped.data()+ped.size(), ped_copy.get());
  return calc_mes(pe_pmf, ped.size(), ped_copy.get(), ped_is_fft);
}

calin::ix::calib::pmt_ses_models::LombardMartinPrescottPMTModelConfig
LombardMartinPrescottPMTModel::nectarcam_config()
{
  calin::ix::calib::pmt_ses_models::LombardMartinPrescottPMTModelConfig cfg;
  cfg.set_num_stage(7);
  cfg.set_total_gain(40000.0);
  cfg.set_stage_0_hi_gain(12.0);
  cfg.set_stage_0_lo_prob(0.15);
  cfg.set_stage_0_lo_gain(3.0);
  cfg.set_suppress_zero(true);
  cfg.set_apply_downsampling(false);
  cfg.set_downsampling_num_stage(4);
  cfg.set_downsampling_factor(16);
  return cfg;
}
