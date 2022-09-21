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
#include <math/rng_vcl.hpp>

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

#ifndef SWIG
  template<typename VCLArchitecture> void vcl_add_nsb(
    calin::math::rng::VCLRNG<VCLArchitecture>& vcl_rng, double nsb_rate_ghz,
    calin::simulation::detector_efficiency::SplinePEAmplitudeGenerator* nsb_pegen = nullptr,
    bool ac_couple=true)
  {
    const double dx = trace_sampling_inv_/nsb_rate_ghz;
    const double xmax = npixels_*trace_nsamples_;
    typename VCLArchitecture::double_vt vx = dx * vcl_rng.exponential_double();
    typename VCLArchitecture::double_at ax;
    vx.store(ax);

    double x = ax[0];
    while(x < xmax) {
      typename VCLArchitecture::double_vt vamp =
        nsb_pegen==nullptr? 1.0 : nsb_pegen->vcl_generate_amplitude(vcl_rng);
      typename VCLArchitecture::double_at aamp;
      vamp.store(aamp);
      pe_waveform_[unsigned(floor(x))] += aamp[0];
      for(unsigned i=1;i<VCLArchitecture::num_double;++i) {
        x += ax[i];
        if(x<xmax) {
          pe_waveform_[unsigned(floor(x))] += aamp[i];
        } else {
          goto break_to_outer_loop;
        }
      }
      vx = dx * vcl_rng.exponential_double();
      vx.store(ax);
      x += ax[0];
    }
break_to_outer_loop:
    if(ac_couple) {
      double mean_amp = nsb_pegen==nullptr? 1.0 : nsb_pegen->mean_amplitude();
      ac_coupling_constant_ += nsb_rate_ghz*trace_sampling_ns_*mean_amp;
    }
    pe_waveform_dft_valid_ = false;
  }

  template<typename VCLArchitecture> void vcl_add_nsb_unroll(
    calin::math::rng::VCLRNG<VCLArchitecture>& vcl_rng_a,
    calin::math::rng::VCLRNG<VCLArchitecture>& vcl_rng_b,
    double nsb_rate_ghz,
    calin::simulation::detector_efficiency::SplinePEAmplitudeGenerator* nsb_pegen = nullptr,
    bool ac_couple=true)
  {
    const double dx = trace_sampling_inv_/nsb_rate_ghz;
    const double xmax = npixels_*trace_nsamples_;

    typename VCLArchitecture::double_vt vx_a;
    typename VCLArchitecture::double_at ax_a;
    typename VCLArchitecture::int32_at axi_a;

    typename VCLArchitecture::double_vt vx_b;
    typename VCLArchitecture::double_at ax_b;
    typename VCLArchitecture::int32_at axi_b;

    vx_a = dx * vcl_rng_a.exponential_double();
    vx_a.store(ax_a);

    axi_a[0] = ax_a[0];
    __builtin_prefetch(pe_waveform_ + axi_a[0]*sizeof(double));
    for(unsigned i=1;i<VCLArchitecture::num_double;++i) {
      ax_a[i] += ax_a[i-1];
      axi_a[i] = ax_a[i];
      __builtin_prefetch(pe_waveform_ + axi_a[i]*sizeof(double));
    }

    vx_b = dx * vcl_rng_b.exponential_double();
    vx_b.store(ax_b);

    ax_b[0] += ax_a[VCLArchitecture::num_double-1];
    axi_b[0] = ax_b[0];
    __builtin_prefetch(pe_waveform_ + axi_b[0]*sizeof(double));
    for(unsigned i=1;i<VCLArchitecture::num_double;++i) {
      ax_b[i] += ax_b[i-1];
      axi_b[i] = ax_b[i];
      __builtin_prefetch(pe_waveform_ + axi_b[i]*sizeof(double));
    }

    typename VCLArchitecture::double_vt vamp_a;
    typename VCLArchitecture::double_at aamp_a;
    typename VCLArchitecture::double_vt vamp_b;
    typename VCLArchitecture::double_at aamp_b;

    while(ax_a[0] < xmax) {
      vamp_a = nsb_pegen==nullptr? 1.0 : nsb_pegen->vcl_generate_amplitude(vcl_rng_a);
      vamp_a.store(aamp_a);

      vamp_b = nsb_pegen==nullptr? 1.0 : nsb_pegen->vcl_generate_amplitude(vcl_rng_b);
      vamp_b.store(aamp_b);

      pe_waveform_[axi_a[0]] += aamp_a[0];
      for(unsigned i=1;i<VCLArchitecture::num_double;++i) {
        if(ax_a[i]<xmax) {
          pe_waveform_[axi_a[i]] += aamp_a[i];
        } else {
          goto break_to_outer_loop;
        }
      }

      pe_waveform_[axi_b[0]] += aamp_b[0];
      for(unsigned i=1;i<VCLArchitecture::num_double;++i) {
        if(ax_b[i]<xmax) {
          pe_waveform_[axi_b[i]] += aamp_b[i];
        } else {
          goto break_to_outer_loop;
        }
      }

      vx_a = dx * vcl_rng_a.exponential_double();
      vx_a.store(ax_a);

      ax_a[0] += ax_b[VCLArchitecture::num_double-1];
      axi_a[0] = ax_a[0];
      __builtin_prefetch(pe_waveform_ + axi_a[0]*sizeof(double));
      for(unsigned i=1;i<VCLArchitecture::num_double;++i) {
        ax_a[i] += ax_a[i-1];
        axi_a[i] = ax_a[i];
        __builtin_prefetch(pe_waveform_ + axi_a[i]*sizeof(double));
      }

      vx_b = dx * vcl_rng_b.exponential_double();
      vx_b.store(ax_b);

      ax_b[0] += ax_a[VCLArchitecture::num_double-1];
      axi_b[0] = ax_b[0];
      __builtin_prefetch(pe_waveform_ + axi_b[0]*sizeof(double));
      for(unsigned i=1;i<VCLArchitecture::num_double;++i) {
        ax_b[i] += ax_b[i-1];
        axi_b[i] = ax_b[i];
        __builtin_prefetch(pe_waveform_ + axi_b[i]*sizeof(double));
      }
    }
break_to_outer_loop:
    if(ac_couple) {
      double mean_amp = nsb_pegen==nullptr? 1.0 : nsb_pegen->mean_amplitude();
      ac_coupling_constant_ += nsb_rate_ghz*trace_sampling_ns_*mean_amp;
    }
    pe_waveform_dft_valid_ = false;
  }

#endif

  void vcl128_add_nsb(calin::math::rng::VCLRNG<calin::util::vcl::VCL128Architecture>& vcl_rng, double nsb_rate_ghz,
      calin::simulation::detector_efficiency::SplinePEAmplitudeGenerator* nsb_pegen = nullptr,
      bool ac_couple=true);
  void vcl256_add_nsb(calin::math::rng::VCLRNG<calin::util::vcl::VCL256Architecture>& vcl_rng, double nsb_rate_ghz,
      calin::simulation::detector_efficiency::SplinePEAmplitudeGenerator* nsb_pegen = nullptr,
      bool ac_couple=true);
  void vcl512_add_nsb(calin::math::rng::VCLRNG<calin::util::vcl::VCL512Architecture>& vcl_rng, double nsb_rate_ghz,
      calin::simulation::detector_efficiency::SplinePEAmplitudeGenerator* nsb_pegen = nullptr,
      bool ac_couple=true);

  void vcl256_add_nsb_unroll(calin::math::rng::VCLRNG<calin::util::vcl::VCL256Architecture>& vcl_rng_a,
      calin::math::rng::VCLRNG<calin::util::vcl::VCL256Architecture>& vcl_rng_b, double nsb_rate_ghz,
      calin::simulation::detector_efficiency::SplinePEAmplitudeGenerator* nsb_pegen = nullptr,
      bool ac_couple=true);

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
