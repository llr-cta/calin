/*

   calin/simulation/pe_processor.cpp -- Stephen Fegan -- 2017-01-16

   Multi-purpose PE (weight, scope, pixel & time) processor.

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

#include <cmath>
#include <stdexcept>
#include <algorithm>

#include <fftw3.h>

#include <simulation/pe_processor.hpp>
#include <math/fftw_util.hpp>

using namespace calin::simulation::pe_processor;

#include <util/log.hpp>
using namespace calin::util::log;

PEProcessor::~PEProcessor()
{
  // nothing to see here
}

void PEProcessor::start_processing()
{
  // nothing to see here
}

void PEProcessor::process_focal_plane_hit(unsigned scope_id, int pixel_id,
    double x, double y, double ux, double uy, double t, double pe_weight)
{
  // nothing to see here
}

void PEProcessor::finish_processing()
{
  // nothing to see here
}

RecordingPEProcessor::RecordingPEProcessor(unsigned nmax, bool auto_clear):
  PEProcessor(), nmax_(nmax), auto_clear_(auto_clear)
{
  // nothing to see here
}

RecordingPEProcessor::~RecordingPEProcessor()
{
  // nothing to see here
}

void RecordingPEProcessor::start_processing()
{
  if(auto_clear_) {
    clear_all_pes();
  }
}

void RecordingPEProcessor::process_focal_plane_hit(unsigned scope_id, int pixel_id,
    double x, double y, double ux, double uy, double t, double pe_weight)
{
  if(nmax_==0 or x_.size()<nmax_) {
    sid_.push_back(scope_id);
    pid_.push_back(pixel_id);
    x_.push_back(x);
    y_.push_back(y);
    ux_.push_back(ux);
    uy_.push_back(uy);
    t_.push_back(t);
    w_.push_back(pe_weight);
  }
}

void RecordingPEProcessor::clear_all_pes()
{
  sid_.clear();
  pid_.clear();
  x_.clear();
  y_.clear();
  ux_.clear();
  uy_.clear();
  t_.clear();
  w_.clear();
}

SimpleImagePEProcessor::
SimpleImagePEProcessor(unsigned nscope, unsigned npix, bool auto_clear):
  PEProcessor(), auto_clear_(auto_clear),
  images_(nscope, std::vector<Accumulator>(npix))
{
  // nothing to see here
}

SimpleImagePEProcessor::
SimpleImagePEProcessor(const std::vector<unsigned> npix, bool auto_clear):
  PEProcessor(), auto_clear_(auto_clear), images_()
{
  for(auto inpix : npix)images_.emplace_back(inpix);
}

SimpleImagePEProcessor::~SimpleImagePEProcessor()
{
  // nothing to see here
}

void SimpleImagePEProcessor::start_processing()
{
  if(auto_clear_)clear_all_images();
}

void SimpleImagePEProcessor::
process_focal_plane_hit(unsigned scope_id, int pixel_id,
  double x, double y, double ux, double uy, double t, double pe_weight)
{
#if 0
  static unsigned counter = 0;
  if(counter++ < 10)
    LOG(INFO) << scope_id << ' ' << pixel_id << ' ' << x << ' ' << y << ' '
              << t << ' ' << pe_weight;
#endif

  if(pixel_id<0)return;
  if(scope_id >= images_.size())
    throw std::out_of_range("SimpleImagePEProcessor::process_focal_plane_hit: scope_id out "
      "of range, " + std::to_string(scope_id));
  if(unsigned(pixel_id) >= images_[scope_id].size())
    throw std::out_of_range("SimpleImagePEProcessor::process_focal_plane_hit: pixel_id out "
      "of range, " + std::to_string(pixel_id));
  images_[scope_id][pixel_id].accumulate(pe_weight);
}

const std::vector<double>
SimpleImagePEProcessor::scope_image(unsigned iscope) const
{
  if(iscope >= images_.size())
    throw std::out_of_range("SimpleImagePEProcessor::scope_image: iscope out "
      "of range");
  std::vector<double> image(images_[iscope].size());
  std::transform(images_[iscope].begin(), images_[iscope].end(), image.begin(),
    [](const Accumulator& acc){ return acc.total(); });
  return image;
}

void SimpleImagePEProcessor::clear_all_images()
{
  for(auto& image : images_)
    std::for_each(image.begin(), image.end(),
      [](Accumulator& acc){ acc.reset(); });
}

TelescopePSFCalcPEProcessor::
TelescopePSFCalcPEProcessor(unsigned iscope, bool auto_clear): PEProcessor(),
  auto_clear_(auto_clear), iscope_(iscope), mom_()
{
  // nothing to see here
}

TelescopePSFCalcPEProcessor::~TelescopePSFCalcPEProcessor()
{
  // nothing to see here
}

void TelescopePSFCalcPEProcessor::start_processing()
{
  if(auto_clear_)clear();
}

void TelescopePSFCalcPEProcessor::process_focal_plane_hit(unsigned scope_id, int pixel_id,
  double x, double y, double ux, double uy, double t, double pe_weight)
{
  if(scope_id == iscope_)mom_.accumulate(x, y, pe_weight);
}

TelescopePSFCalcThirdMomentPEProcessor::
TelescopePSFCalcThirdMomentPEProcessor(unsigned iscope, bool auto_clear):
  PEProcessor(), auto_clear_(auto_clear), iscope_(iscope), mom_()
{
  // nothing to see here
}

TelescopePSFCalcThirdMomentPEProcessor::
~TelescopePSFCalcThirdMomentPEProcessor()
{
  // nothing to see here
}

void TelescopePSFCalcThirdMomentPEProcessor::start_processing()
{
  if(auto_clear_)clear();
}

void TelescopePSFCalcThirdMomentPEProcessor::
process_focal_plane_hit(unsigned scope_id, int pixel_id,
  double x, double y, double ux, double uy, double t, double pe_weight)
{
  if(scope_id == iscope_)mom_.accumulate(x, y, pe_weight);
}

WaveformPEProcessor::WaveformPEProcessor(unsigned nscope, unsigned npix,
    unsigned nsamp, double delta_t, calin::math::rng::RNG* rng, bool auto_clear):
  PEProcessor(), nsamp_(nsamp), npix_(npix), delta_t_inv_(1.0/delta_t),
  traces_(nscope, { npix, nsamp }), t0_(nscope),
  nmin_(nscope), nmax_(nscope), overflow_(nscope, npix_),
  auto_clear_(auto_clear), rng_(rng)
{
  if(nsamp & (nsamp-1)) {
    throw std::runtime_error("nsamp must be power of two : "+std::to_string(nsamp));
  }
  clear_all_traces();
}

WaveformPEProcessor::~WaveformPEProcessor()
{
  // nothing to see here
}

void WaveformPEProcessor::start_processing()
{
  if(auto_clear_) {
    clear_all_traces();
  }
  warning_sent_ = false;
}

void WaveformPEProcessor::process_focal_plane_hit(unsigned scope_id, int pixel_id,
  double x, double y, double ux, double uy, double t, double pe_weight)
{
  if(pixel_id < 0) {
    return;
  } else if (pixel_id >= npix_) {
    throw std::out_of_range("WaveformPEProcessor::process_focal_plane_hit : pixel_id out of range : "
      + std::to_string(pixel_id) + " >= " + std::to_string(npix_));
  }

  double trel;
  if(std::isfinite(t0_[scope_id])) {
    trel = t - t0_[scope_id];
  } else {
    t0_[scope_id] = t;
    trel = 0;
  }

  int n = std::round(trel * delta_t_inv_);
  double nmin = std::min(nmin_[scope_id], n);
  double nmax = std::max(nmax_[scope_id], n);
  if(nmax - nmin < nsamp_) {
    n &= (nsamp_-1);
    traces_[scope_id](pixel_id, n) += pe_weight;
    nmin_[scope_id] = nmin;
    nmax_[scope_id] = nmax;
  } else {
    overflow_(scope_id, pixel_id) += pe_weight;
    if(not warning_sent_) {
      LOG(INFO)
        << "WaveformPEProcessor::process_focal_plane_hit : Circular trace buffer overflow\n"
        << "scope_id=" << scope_id << ", pixel_id=" << pixel_id << ", t=" << t
        << ", t0=" << t0_[scope_id] << ", n=" << n << ", nmin=" << nmin_[scope_id]
        << ", nmax=" << nmax_[scope_id];
      warning_sent_ = true;
    }
  }
}

void WaveformPEProcessor::clear_all_traces()
{
  for(auto& traces : traces_) {
    traces.setZero();
  }
  t0_.setConstant(std::numeric_limits<double>::quiet_NaN());
  nmin_.setZero();
  nmax_.setZero();
  overflow_.setZero();
}

void WaveformPEProcessor::add_nsb(double rate_ghz,
  calin::simulation::detector_efficiency::PEAmplitudeGenerator* pegen)
{
  if(rng_ == nullptr) {
    rng_ = new calin::math::rng::RNG(__PRETTY_FUNCTION__, "NSB photon time generator");
  }
  double delta_t_nsb = delta_t_inv_/rate_ghz;
  for(auto& traces : traces_) {
    double t_max = traces.cols();
    for(unsigned ipix=0;ipix<npix_; ++ipix) {
      double t = delta_t_nsb*rng_->exponential();
      while(t < t_max) {
        double amp = pegen==nullptr? 1.0 : pegen->generate_amplitude();
        traces(ipix, unsigned(floor(t))) += amp;
        t += delta_t_nsb*rng_->exponential();
      }
    }
  }
}

void WaveformPEProcessor::add_nsb(const Eigen::VectorXd rate_per_pixel_ghz,
  calin::simulation::detector_efficiency::PEAmplitudeGenerator* pegen)
{
  if(rate_per_pixel_ghz.size() != npix_) {
    throw std::runtime_error("WaveformPEProcessor::add_nsb : rate_per_pixel_ghz array must have "
      + std::to_string(npix_) + " entries");
  }
  if(rng_ == nullptr) {
    rng_ = new calin::math::rng::RNG(__PRETTY_FUNCTION__, "NSB photon time generator");
  }
  for(auto& traces : traces_) {
    double t_max = traces.cols();
    for(unsigned ipix=0;ipix<traces.rows(); ++ipix) {
      double delta_t_nsb = delta_t_inv_/rate_per_pixel_ghz[ipix];
      double t = delta_t_nsb*rng_->exponential();
      while(t < t_max) {
        double amp = pegen==nullptr? 1.0 : pegen->generate_amplitude();
        traces(ipix, unsigned(floor(t))) += amp;
        t += delta_t_nsb*rng_->exponential();
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
//
// UnbinnedWaveformPEProcessor
//
//
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

UnbinnedWaveformPEProcessor::
UnbinnedWaveformPEProcessor(unsigned nscope, unsigned npix,
  double scope_trace_delta_t, unsigned scope_trace_nsamp,
  bool auto_clear): PEProcessor(), nscope_(nscope), npix_(npix),
  scope_trace_nsamp_(scope_trace_nsamp),
  scope_trace_delta_t_inv_(1.0/scope_trace_delta_t),
  scope_trace_(nscope, scope_trace_nsamp), t0_(nscope), nmin_(nscope), nmax_(nscope),
  scope_trace_overflow_(nscope), auto_clear_(auto_clear)
{
  if(scope_trace_nsamp_ & (scope_trace_nsamp_-1)) {
    throw std::runtime_error("UnbinnedWaveformPEProcessor : scope_trace_nsamp_ must be power of two : "+std::to_string(scope_trace_nsamp_));
  }
  clear_all_traces();
}

UnbinnedWaveformPEProcessor::~UnbinnedWaveformPEProcessor()
{
  // nothing to see here
}

void UnbinnedWaveformPEProcessor::start_processing()
{
  if(auto_clear_) {
    clear_all_traces();
  }
}

void UnbinnedWaveformPEProcessor::process_focal_plane_hit(unsigned scope_id, int pixel_id,
  double x, double y, double ux, double uy, double t, double pe_weight)
{
  if(pixel_id < 0) {
    return;
  } else if (pixel_id >= npix_) {
    throw std::out_of_range("UnbinnedWaveformPEProcessor::process_focal_plane_hit : pixel_id out of range : "
      + std::to_string(pixel_id) + " >= " + std::to_string(npix_));
  }

  pe_iscope_.push_back(scope_id);
  pe_ipix_.push_back(pixel_id);
  pe_t_.push_back(t);
  pe_q_.push_back(pe_weight);

  if(not std::isfinite(t0_[scope_id])) {
    t0_[scope_id] = std::round(t * scope_trace_delta_t_inv_);
  }

  double trel = std::round(t * scope_trace_delta_t_inv_) - t0_[scope_id];

  int n = trel;
  double nmin = std::min(nmin_[scope_id], n);
  double nmax = std::max(nmax_[scope_id], n);
  if(nmax - nmin < scope_trace_nsamp_) {
    n &= (scope_trace_nsamp_-1);
    scope_trace_(scope_id, n) += pe_weight;
    nmin_[scope_id] = nmin;
    nmax_[scope_id] = nmax;
  } else {
    scope_trace_overflow_(scope_id) += pe_weight;
    if(not warning_sent_) {
      LOG(INFO)
        << "UnbinnedWaveformPEProcessor::process_focal_plane_hit : Circular trace buffer overflow\n"
        << "scope_id=" << scope_id << ", pixel_id=" << pixel_id << ", t=" << t
        << ", t0=" << t0_(scope_id) << ", n=" << n << ", nmin=" << nmin_(scope_id)
        << ", nmax=" << nmax_(scope_id);
      warning_sent_ = true;
    }
  }

}

Eigen::MatrixXd UnbinnedWaveformPEProcessor::pixel_traces(
  double& trace_t0, Eigen::VectorXd& trace_overflow,
  unsigned iscope,
  double trace_delta_t, unsigned trace_nsamp, double trace_advance_time,
  double nsb_rate_ghz, calin::math::rng::RNG* rng,
  calin::simulation::detector_efficiency::PEAmplitudeGenerator* nsb_pegen,
  bool ac_couple) const
{
  check_iscope(iscope);

  Eigen::MatrixXd pix_traces(trace_nsamp, npix_);
  Eigen::VectorXd pix_overflow(npix_);
  if(ac_couple and nsb_rate_ghz>0 and rng!=nullptr) {
    double mean_amplitude = nsb_pegen==nullptr ? 1.0 : nsb_pegen->mean_amplitude();
    pix_traces.setConstant(-trace_delta_t*nsb_rate_ghz*mean_amplitude);
  } else {
    pix_traces.setZero();
  }
  pix_overflow.setZero();

  double trace_delta_t_inv = 1.0/trace_delta_t;

  if(std::isfinite(t0_[iscope])) {
    Eigen::VectorXd pe_integral(nmax_(iscope) - nmin_(iscope) + 2);
    double I = 0;
    for(int n=nmin_(iscope);n<=nmax_(iscope);++n) {
      pe_integral[n - nmin_(iscope)] = I;
      I += scope_trace_(iscope, n & (scope_trace_nsamp_-1));
    }
    pe_integral[nmax_(iscope) - nmin_(iscope) + 1] = I;

    auto nmedian = std::upper_bound(pe_integral.begin(), pe_integral.end(), I*0.5);
    double tfrac = (I*0.5 - *(nmedian-1))/(*nmedian - *(nmedian-1));
    double tmedian = (t0_[iscope] + double(nmedian-pe_integral.begin()+nmin_(iscope)-1) + tfrac - 0.5)/scope_trace_delta_t_inv_;

    double tstart = std::round((tmedian - trace_advance_time)*trace_delta_t_inv);

#if 0
    LOG(INFO) << "t0=" << t0_[iscope] << " nmin=" << nmin_(iscope)
      << " nmedian=" << nmedian-pe_integral.begin() << " tfrac=" << tfrac
      << " tmedian=" << tmedian << " tstart=" << tstart;
#endif

    for(unsigned ipe=0;ipe<pe_iscope_.size();++ipe) {
      if(pe_iscope_[ipe] == iscope) {
        int n = std::round(pe_t_[ipe] * trace_delta_t_inv) - tstart;
        if(n>=0 and n<=trace_nsamp) {
          pix_traces(n, pe_ipix_[ipe]) += pe_q_[ipe];
        } else {
          pix_overflow(pe_ipix_[ipe]) += pe_q_[ipe];
        }
      }
    }
    trace_t0 = tstart*trace_delta_t;
    trace_overflow = pix_overflow;
  } else {
    trace_t0 = t0_[iscope];
  }

  if(nsb_rate_ghz>0 and rng!=nullptr) {
    double dx = trace_delta_t_inv/nsb_rate_ghz;
    double xmax = pix_traces.size();
    double x = dx*rng->exponential();
    while(x < xmax) {
      double amp = nsb_pegen==nullptr? 1.0 : nsb_pegen->generate_amplitude();
      pix_traces.data()[unsigned(floor(x))] += amp;
      x += dx*rng->exponential();
    }
  }

  return pix_traces;
}

void UnbinnedWaveformPEProcessor::clear_all_traces()
{
  pe_iscope_.clear();
  pe_ipix_.clear();
  pe_t_.clear();
  pe_q_.clear();
  scope_trace_.setZero();
  t0_.setConstant(std::numeric_limits<double>::quiet_NaN());
  nmin_.setZero();
  nmax_.setZero();
  scope_trace_overflow_.setZero();
  warning_sent_ = false;
}

Eigen::MatrixXd UnbinnedWaveformPEProcessor::convolve_instrument_response(
  const Eigen::MatrixXd& traces, const Eigen::VectorXd& impulse_response_dft,
  double pedestal)
{
  if(traces.rows() != impulse_response_dft.size()) {
    throw std::length_error("convolve_instrument_response: number of columns in traces does not match impulse response");
  }

  double* traces_a = fftw_alloc_real(traces.size());
  double* traces_b = fftw_alloc_real(traces.size());

  int rank = 1;
  int n[] = { static_cast<int>(traces.rows()) };
  int howmany = traces.cols();
  double* in = traces_a;
  int* inembed = n;
  int istride = 1;
  int idist = traces.rows();
  double* out = traces_b;
  int* onembed = n;
  int ostride = 1;
  int odist = traces.rows();
  fftw_r2r_kind kind[] = { FFTW_R2HC };

  auto fwd_plan = fftw_plan_many_r2r(rank, n, howmany,
    in, inembed, istride, idist, out, onembed, ostride, odist, kind, 0);

  kind[0] = FFTW_HC2R;

  auto rev_plan = fftw_plan_many_r2r(rank, n, howmany,
    in, inembed, istride, idist, out, onembed, ostride, odist, kind, 0);

  std::copy(traces.data(), traces.data()+traces.size(), traces_a);

  fftw_execute(fwd_plan);

  for(unsigned icol=0; icol<traces.cols(); ++icol) {
    calin::math::fftw_util::hcvec_scale_and_multiply(traces_a + icol*traces.rows(),
      traces_b + icol*traces.rows(), impulse_response_dft.data(),
      traces.rows(), 1.0/traces.rows());
    *(traces_a + icol*traces.rows()) += pedestal;
  }

  fftw_execute(rev_plan);

  Eigen::MatrixXd convolved_traces(traces.rows(), traces.cols());
  std::copy(traces_b, traces_b+traces.size(), convolved_traces.data());

  fftw_destroy_plan(fwd_plan);
  fftw_destroy_plan(rev_plan);
  fftw_free(traces_a);
  fftw_free(traces_b);

  return convolved_traces;
}
