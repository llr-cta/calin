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
#include <simulation/pe_processor.hpp>

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
    double x, double y, double ux, double uy, double t0, double pe_weight)
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
    double x, double y, double ux, double uy, double t0, double pe_weight)
{
  if(nmax_==0 or x_.size()<nmax_) {
    sid_.push_back(scope_id);
    pid_.push_back(pixel_id);
    x_.push_back(x);
    y_.push_back(y);
    ux_.push_back(ux);
    uy_.push_back(uy);
    t_.push_back(t0);
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
  double x, double y, double ux, double uy, double t0, double pe_weight)
{
#if 0
  static unsigned counter = 0;
  if(counter++ < 10)
    LOG(INFO) << scope_id << ' ' << pixel_id << ' ' << x << ' ' << y << ' '
              << t0 << ' ' << pe_weight;
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
  double x, double y, double ux, double uy, double t0, double pe_weight)
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
  double x, double y, double ux, double uy, double t0, double pe_weight)
{
  if(scope_id == iscope_)mom_.accumulate(x, y, pe_weight);
}

WaveformPEProcessor::WaveformPEProcessor(unsigned nscope, unsigned npix,
    unsigned nsamp, double delta_t, calin::math::rng::RNG* rng, bool auto_clear):
  PEProcessor(), nsamp_(nsamp), npix_(npix), delta_t_inv_(1.0/delta_t),
  traces_(nscope, { npix, nsamp }), t0_(nscope),
  nmin_(nscope), nmax_(nscope), auto_clear_(auto_clear), rng_(rng)
{
  if(nsamp & (nsamp-1)) {
    throw std::runtime_error("nsamp must be power of two : "+std::to_string(nsamp));
  }
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
  double x, double y, double ux, double uy, double t0, double pe_weight)
{
  if(pixel_id < 0) {
    return;
  } else if (pixel_id >= npix_) {
    throw std::out_of_range("WaveformPEProcessor::process_focal_plane_hit : pixel_id out of range : "
      + std::to_string(pixel_id) + " >= " + std::to_string(npix_));
  }

  double trel;
  if(std::isfinite(t0_[scope_id])) {
    trel = t0 - t0_[scope_id];
  } else {
    t0_[scope_id] = t0;
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
  } else if(not warning_sent_) {
    LOG(INFO)
      << "WaveformPEProcessor::process_focal_plane_hit : Circular trace buffer overflow\n"
      << "scope_id=" << scope_id << ", pixel_id=" << pixel_id << ", t=" << t0
      << ", t0=" << t0_[scope_id] << ", n=" << n << ", nmin=" << nmin_[scope_id]
      << ", nmax=" << nmax_[scope_id];
    warning_sent_ = true;
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
}

void WaveformPEProcessor::add_nsb(double rate_ghz)
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
        traces(ipix, unsigned(floor(t))) += 1.0;
        t += delta_t_nsb*rng_->exponential();
      }
    }
  }
}

void WaveformPEProcessor::add_nsb(const Eigen::VectorXd rate_per_pixel_ghz)
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
        traces(ipix, unsigned(floor(t))) += 1.0;
        t += delta_t_nsb*rng_->exponential();
      }
    }
  }
}
