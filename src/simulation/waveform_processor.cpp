/*

   calin/simulation/waveform_processor.hpp -- Stephen Fegan -- 2022-09-16

   Class to process waveforms

   Copyright 2022, Stephen Fegan <sfegan@llr.in2p3.fr>
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
#include <alloca.h>

#include <fftw3.h>

#include <simulation/waveform_processor.hpp>
#include <math/fftw_util.hpp>

#include <util/log.hpp>
using namespace calin::util::log;

using namespace calin::simulation::waveform_processor;

WaveformProcessor::
WaveformProcessor(unsigned npixels, double trace_sampling_ns, unsigned trace_nsamples,
    double trace_advance_time, calin::math::rng::RNG* rng, unsigned fftw_flags, bool adopt_rng):
  npixels_(npixels), trace_sampling_ns_(trace_sampling_ns),
  trace_sampling_inv_(1.0/trace_sampling_ns),
  trace_nsamples_(trace_nsamples), trace_advance_time_(trace_advance_time),
  rng_(rng==nullptr ? new calin::math::rng::RNG(__PRETTY_FUNCTION__, "NSB and electronics noise generation") : rng),
  adopt_rng_(rng==nullptr ? true : adopt_rng),
  pe_waveform_(fftw_alloc_real(npixels * trace_nsamples)),
  pe_waveform_dft_(fftw_alloc_real(npixels * trace_nsamples)),
  el_waveform_dft_(fftw_alloc_real(npixels * trace_nsamples)),
  el_waveform_(fftw_alloc_real(npixels * trace_nsamples))
{
  int rank = 1;
  int n[] = { static_cast<int>(trace_nsamples_) };
  int howmany = npixels_;
  double* in = pe_waveform_;
  int* inembed = n;
  int istride = 1;
  int idist = trace_nsamples_;
  double* out = pe_waveform_dft_;
  int* onembed = n;
  int ostride = 1;
  int odist = trace_nsamples_;
  fftw_r2r_kind kind[] = { FFTW_R2HC };

  fwd_plan_ = fftw_plan_many_r2r(rank, n, howmany,
    in, inembed, istride, idist, out, onembed, ostride, odist, kind, fftw_flags);

  in = el_waveform_dft_;
  out = el_waveform_;
  kind[0] = FFTW_HC2R;

  rev_plan_ = fftw_plan_many_r2r(rank, n, howmany,
    in, inembed, istride, idist, out, onembed, ostride, odist, kind, fftw_flags);
}

WaveformProcessor::~WaveformProcessor()
{
  fftw_destroy_plan(fwd_plan_);
  fftw_destroy_plan(rev_plan_);
  fftw_free(pe_waveform_);
  fftw_free(pe_waveform_dft_);
  fftw_free(el_waveform_dft_);
  fftw_free(el_waveform_);
  if(adopt_rng_) {
    delete rng_;
  }
}

void WaveformProcessor::load_pes_from_processor(
  calin::simulation::pe_processor::UnbinnedWaveformPEProcessor* pe_processor,
  unsigned iscope)
{
  clear_pes();
  pe_processor->pixel_traces_into_buffer(pe_waveform_, nullptr, wavewform_t0_,
    iscope, trace_sampling_ns_, trace_nsamples_, trace_advance_time_,
    /* pe_waveform_buffer_stride= */ trace_nsamples_);
}

void WaveformProcessor::add_nsb(double nsb_rate_ghz,
  calin::simulation::detector_efficiency::PEAmplitudeGenerator* nsb_pegen,
  bool ac_couple)
{
  double dx = trace_sampling_inv_/nsb_rate_ghz;
  double xmax = npixels_*trace_nsamples_;
  double x = dx*rng_->exponential();
  while(x < xmax) {
    double amp = nsb_pegen==nullptr? 1.0 : nsb_pegen->generate_amplitude();
    pe_waveform_[unsigned(floor(x))] += amp;
    x += dx*rng_->exponential();
  }
  if(ac_couple) {
    double mean_amp = nsb_pegen==nullptr? 1.0 : nsb_pegen->mean_amplitude();
    ac_coupling_constant_ += nsb_rate_ghz*trace_sampling_ns_*mean_amp;
  }
  pe_waveform_dft_valid_ = false;
}

void WaveformProcessor::convolve_unity_impulse_response(double pedestal)
{
  compute_pe_waveform_dft();
  pedestal -= ac_coupling_constant_;
  double scale = 1.0/trace_nsamples_;
  for(unsigned ipixel=0; ipixel<npixels_; ++ipixel) {
    std::transform(pe_waveform_dft_ + ipixel*trace_nsamples_,
      pe_waveform_dft_ + ipixel*trace_nsamples_ + trace_nsamples_,
      el_waveform_dft_ + ipixel*trace_nsamples_,
      [scale](double x) { return x*scale; });
    *(el_waveform_dft_ + ipixel*trace_nsamples_) += pedestal;
  }
  el_waveform_dft_valid_ = true;
  el_waveform_valid_ = false;
}

void WaveformProcessor::convolve_impulse_response(const double* impulse_response_dft, double pedestal)
{
  compute_pe_waveform_dft();
  pedestal -= ac_coupling_constant_*impulse_response_dft[0];
  double scale = 1.0/trace_nsamples_;
  for(unsigned ipixel=0; ipixel<npixels_; ++ipixel) {
    calin::math::fftw_util::hcvec_scale_and_multiply(
      el_waveform_dft_ + ipixel*trace_nsamples_,
      pe_waveform_dft_ + ipixel*trace_nsamples_,
      impulse_response_dft,
      trace_nsamples_, scale);
    *(el_waveform_dft_ + ipixel*trace_nsamples_) += pedestal;
  }
  el_waveform_dft_valid_ = true;
  el_waveform_valid_ = false;
}

void WaveformProcessor::convolve_impulse_response(
  const Eigen::VectorXd& impulse_response_dft, double pedestal)
{
  if(impulse_response_dft.size() != trace_nsamples_) {
    throw std::runtime_error("convolve_impulse_response : impulse response DFT have " + std::to_string(trace_nsamples_) + " points");
  }
  convolve_impulse_response(impulse_response_dft.data(), pedestal);
}


void WaveformProcessor::add_electronics_noise(const double* noise_spectrum_amplitude)
{
  if(not el_waveform_dft_valid_) {
    throw std::runtime_error("add_electronics_noise : impulse response must be applied before noise added");
  }
  for(unsigned ipixel=0; ipixel<npixels_; ++ipixel) {
    for(unsigned isample=0; isample<trace_nsamples_; ++isample) {
      *(el_waveform_dft_ + ipixel*trace_nsamples_ + isample) +=
        *(noise_spectrum_amplitude + isample) * rng_->normal();
    }
  }
  el_waveform_valid_ = false;
}

void WaveformProcessor::add_electronics_noise(const Eigen::VectorXd& noise_spectrum_amplitude)
{
  if(noise_spectrum_amplitude.size() != trace_nsamples_) {
    throw std::runtime_error("add_electronics_noise : noise spectrum amplitude must have " + std::to_string(trace_nsamples_) + " points");
  }
  add_electronics_noise(noise_spectrum_amplitude.data());
}

Eigen::MatrixXd WaveformProcessor::pe_waveform() const
{
  Eigen::MatrixXd pe_traces(trace_nsamples_, npixels_);
  std::copy(pe_waveform_, pe_waveform_+trace_nsamples_*npixels_, pe_traces.data());
  return pe_traces;
}

Eigen::MatrixXd WaveformProcessor::el_waveform()
{
  compute_el_waveform();
  Eigen::MatrixXd el_traces(trace_nsamples_, npixels_);
  std::copy(el_waveform_, el_waveform_+trace_nsamples_*npixels_, el_traces.data());
  return el_traces;
}

void WaveformProcessor::compute_pe_waveform_dft()
{
  if(not pe_waveform_dft_valid_) {
    fftw_execute(fwd_plan_);
    pe_waveform_dft_valid_ = true;
  }
}

void WaveformProcessor::compute_el_waveform()
{
  if(not el_waveform_valid_) {
    if(not el_waveform_dft_valid_) {
      throw std::runtime_error("compute_el_waveform : impulse response must be waveform can be computed");
    }
    fftw_execute(rev_plan_);
    el_waveform_valid_ = true;
  }
}

void WaveformProcessor::clear_pes()
{
  pe_waveform_dft_valid_ = false;
  el_waveform_dft_valid_ = false;
  el_waveform_valid_ = false;
  ac_coupling_constant_ = 0;
  wavewform_t0_ = std::numeric_limits<double>::quiet_NaN();
  std::fill(pe_waveform_, pe_waveform_ + npixels_*trace_nsamples_, 0);
}

int WaveformProcessor::digital_multipicity_trigger(double threshold,
  unsigned time_over_threshold_samples, unsigned coherence_time_samples,
  unsigned multiplicity_threshold)
{
  compute_el_waveform();
  unsigned* l0_eop = static_cast<unsigned*>(alloca(npixels_ * sizeof(unsigned)));
  unsigned* l0_tot = static_cast<unsigned*>(alloca(npixels_ * sizeof(unsigned)));
  std::fill(l0_eop, l0_eop+npixels_, 0);
  std::fill(l0_tot, l0_eop+npixels_, 0);
  for(unsigned isamp=0; isamp<trace_nsamples_; ++isamp) {
    unsigned multiplicity = 0;
    for(unsigned ipixel=0; ipixel<npixels_; ++ipixel) {
      if(el_waveform_[ipixel*trace_nsamples_ + isamp] > threshold) {
        ++l0_tot[ipixel];
      } else {
        l0_tot[ipixel] = 0;
      }
      if(l0_tot[ipixel] >= time_over_threshold_samples) {
        l0_eop[ipixel] = isamp+coherence_time_samples;
      }
      if(l0_eop[ipixel] > isamp) {
        ++multiplicity;
        if(multiplicity >= multiplicity_threshold) {
          return isamp;
        }
      }
    }
  }
  return -1;
}
