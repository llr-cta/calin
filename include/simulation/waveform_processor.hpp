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

#pragma once

// Units:
// Time, t:      ns
// Height, z:    cm
// Energy, e:    MeV

#include <vector>
#include <ostream>
#include <map>

#include <fftw3.h>

#include <Eigen/Core>
#include <simulation/pe_processor.hpp>
#include <iact_data/instrument_layout.pb.h>
#include <math/fftw_util.hpp>

namespace calin { namespace simulation { namespace waveform_processor {

class WaveformProcessor
{
public:
  WaveformProcessor(unsigned npixels, double trace_sampling_ns, unsigned trace_nsamples,
    double trace_advance_time, calin::math::rng::RNG* rng = nullptr,
    unsigned fftw_flags = 0, bool adopt_rng = false);
  WaveformProcessor(const calin::ix::iact_data::instrument_layout::CameraLayout* camera,
    double trace_sampling_ns, unsigned trace_nsamples,
    double trace_advance_time, calin::math::rng::RNG* rng = nullptr,
    unsigned fftw_flags = 0, bool adopt_rng = false);
  ~WaveformProcessor();
  void load_pes_from_processor(
    calin::simulation::pe_processor::UnbinnedWaveformPEProcessor* pe_processor,
    unsigned iscope);
  void add_nsb(double nsb_rate_ghz,
    calin::simulation::detector_efficiency::PEAmplitudeGenerator* nsb_pegen = nullptr,
    bool ac_couple=true);
  void convolve_unity_impulse_response(double pedestal = 0);
#ifndef SWIG
  void convolve_impulse_response(const double* impulse_response_dft, double pedestal = 0);
  void add_electronics_noise(const double* noise_spectrum_amplitude);
#endif
  void convolve_impulse_response(const Eigen::VectorXd& impulse_response_dft, double pedestal = 0);
  void add_electronics_noise(const Eigen::VectorXd& noise_spectrum_amplitude);
  void clear_pes();

  Eigen::MatrixXd pe_waveform() const;
  Eigen::MatrixXd el_waveform();

  int digital_multipicity_trigger(double threshold,
    unsigned time_over_threshold_samples, unsigned coherence_time_samples,
    unsigned multiplicity_threshold);
  int digital_nn_trigger(double threshold,
    unsigned time_over_threshold_samples, unsigned coherence_time_samples,
    unsigned nn_threshold);

  double wavewform_t0() { return wavewform_t0_; }
  double ac_coupling_constant() { return ac_coupling_constant_; }

  template<typename VCLArchitecture> void vcl_add_nsb(
    calin::math::rng::VCLRNG<VCLArchitecture>& vcl_rng, double nsb_rate_ghz,
    calin::simulation::detector_efficiency::SplinePEAmplitudeGenerator* nsb_pegen = nullptr,
    bool ac_couple=true)
  {
    if(npixels_*trace_nsamples_ % VCLArchitecture::num_double != 0) {
      throw std::logic_error("vcl_add_nsb : vector size must divide evenly into array size");
    }
    int64_t base_size = npixels_*trace_nsamples_/VCLArchitecture::num_double;
    typename VCLArchitecture::int64_vt pes_base = VCLArchitecture::int64_iota() * base_size;
    double dx = trace_sampling_inv_/nsb_rate_ghz;
    double xmax = base_size;
    typename VCLArchitecture::double_vt x = dx * vcl_rng.exponential_double();
    typename VCLArchitecture::double_at pes_array;
    typename VCLArchitecture::uint64_at pes_index_array;
    while(vcl::horizontal_or(x<xmax)) {
      typename VCLArchitecture::int64_vt pes_index = pes_base +
        vcl::min(base_size, vcl::truncate_to_int64_limited(vcl::floor(x)));
      typename VCLArchitecture::double_vt pes;
      pes = vcl::lookup<0x40000000>(pes_index, pe_waveform_);
      typename VCLArchitecture::double_vt amp =
        nsb_pegen==nullptr? 1.0 : nsb_pegen->vcl_generate_amplitude(vcl_rng);
      pes += amp;
      pes.store(pes_array);
      pes_index.store(pes_index_array);
      if(vcl::horizontal_and(x<xmax)) {
        for(unsigned i=0;i<VCLArchitecture::num_double;++i) {
          pe_waveform_[pes_index_array[i]] = pes_array[i];
        }
      } else {
        unsigned xgood = vcl::to_bits(x<xmax);
        for(unsigned i=0;i<VCLArchitecture::num_double;++i) {
          if(xgood & 0x1<<i) {
            pe_waveform_[pes_index_array[i]] = pes_array[i];
          }
        }
      }
      x += dx * vcl_rng.exponential_double();
    }
    if(ac_couple) {
      double mean_amp = nsb_pegen==nullptr? 1.0 : nsb_pegen->mean_amplitude();
      ac_coupling_constant_ += nsb_rate_ghz*trace_sampling_ns_*mean_amp;
    }
    pe_waveform_dft_valid_ = false;
  }

#if MAX_VECTOR_SIZE >= 256 && INSTRSET >= 7
  void vcl256_add_nsb(calin::math::rng::VCLRNG<calin::util::vcl::VCL256Architecture>& vcl_rng, double nsb_rate_ghz,
      calin::simulation::detector_efficiency::SplinePEAmplitudeGenerator* nsb_pegen = nullptr,
      bool ac_couple=true) {
    vcl_add_nsb<calin::util::vcl::VCL256Architecture>(vcl_rng, nsb_rate_ghz, nsb_pegen, ac_couple);
  }
#endif

private:
  void compute_pe_waveform_dft();
  void compute_el_waveform();
  unsigned npixels_;
  double trace_sampling_ns_;
  double trace_sampling_inv_;
  unsigned trace_nsamples_;
  double trace_advance_time_;
  calin::math::rng::RNG* rng_ = nullptr;
  bool adopt_rng_;

  double wavewform_t0_ = std::numeric_limits<double>::quiet_NaN();
  double* pe_waveform_ = nullptr;
  double* pe_waveform_dft_ = nullptr;
  double* el_waveform_dft_ = nullptr;
  double* el_waveform_ = nullptr;
  fftw_plan fwd_plan_;
  fftw_plan rev_plan_;
  bool pe_waveform_dft_valid_ = false;
  bool el_waveform_dft_valid_ = false;
  bool el_waveform_valid_ = false;
  double ac_coupling_constant_ = 0;

  int* neighbour_map_ = nullptr;
  int max_num_neighbours_ = 0;
};

} } } // namespace calin::simulation::waveform_processor
