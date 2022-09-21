/*

   calin/simulation/waveform_processor_vcl_function_implemntations.hpp -- Stephen Fegan -- 2022-09-21

   Class to process waveforms - VCL function implementations

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

#ifndef SWIG
template<typename VCLArchitecture> void WaveformProcessor::vcl_add_nsb(
  calin::math::rng::VCLRNG<VCLArchitecture>& vcl_rng, double nsb_rate_ghz,
  calin::simulation::detector_efficiency::SplinePEAmplitudeGenerator* nsb_pegen,
  bool ac_couple)
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

template<typename VCLArchitecture> void WaveformProcessor::vcl_add_nsb(
  calin::math::rng::VCLRNG<VCLArchitecture>& vcl_rng_a,
  calin::math::rng::VCLRNG<VCLArchitecture>& vcl_rng_b,
  double nsb_rate_ghz,
  calin::simulation::detector_efficiency::SplinePEAmplitudeGenerator* nsb_pegen,
  bool ac_couple)
{
  const double dx = trace_sampling_inv_/nsb_rate_ghz;
  const double xmax = npixels_*trace_nsamples_;

  typename VCLArchitecture::double_vt vx_a;
  typename VCLArchitecture::double_at ax_a;
  int32_t axi_a[VCLArchitecture::num_double];

  typename VCLArchitecture::double_vt vx_b;
  typename VCLArchitecture::double_at ax_b;
  int32_t axi_b[VCLArchitecture::num_double];


  typename VCLArchitecture::double_vt vamp_a;
  typename VCLArchitecture::double_at aamp_a;
  typename VCLArchitecture::double_vt vamp_b;
  typename VCLArchitecture::double_at aamp_b;

  // Initialze position within buffer
  ax_b[VCLArchitecture::num_double-1] = 0;

  while(true) {
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

    vamp_a = nsb_pegen==nullptr? 1.0 : nsb_pegen->vcl_generate_amplitude(vcl_rng_a);
    vamp_a.store(aamp_a);

    vamp_b = nsb_pegen==nullptr? 1.0 : nsb_pegen->vcl_generate_amplitude(vcl_rng_b);
    vamp_b.store(aamp_b);

    for(unsigned i=0;i<VCLArchitecture::num_double;++i) {
      if(ax_a[i]<xmax) {
        pe_waveform_[axi_a[i]] += aamp_a[i];
      } else {
        goto break_to_outer_loop;
      }
    }

    for(unsigned i=0;i<VCLArchitecture::num_double;++i) {
      if(ax_b[i]<xmax) {
        pe_waveform_[axi_b[i]] += aamp_b[i];
      } else {
        goto break_to_outer_loop;
      }
    }
  }
break_to_outer_loop:
  if(ac_couple) {
    double mean_amp = nsb_pegen==nullptr? 1.0 : nsb_pegen->mean_amplitude();
    ac_coupling_constant_ += nsb_rate_ghz*trace_sampling_ns_*mean_amp;
  }
  pe_waveform_dft_valid_ = false;
}

template<typename VCLArchitecture> void WaveformProcessor::vcl_add_nsb(
  calin::math::rng::VCLRNG<VCLArchitecture>& vcl_rng_a,
  calin::math::rng::VCLRNG<VCLArchitecture>& vcl_rng_b,
  calin::math::rng::VCLRNG<VCLArchitecture>& vcl_rng_c,
  double nsb_rate_ghz,
  calin::simulation::detector_efficiency::SplinePEAmplitudeGenerator* nsb_pegen,
  bool ac_couple)
{
  const double dx = trace_sampling_inv_/nsb_rate_ghz;
  const double xmax = npixels_*trace_nsamples_;

  typename VCLArchitecture::double_vt vx_a;
  typename VCLArchitecture::double_at ax_a;
  int32_t axi_a[VCLArchitecture::num_double];

  typename VCLArchitecture::double_vt vx_b;
  typename VCLArchitecture::double_at ax_b;
  int32_t axi_b[VCLArchitecture::num_double];

  typename VCLArchitecture::double_vt vx_c;
  typename VCLArchitecture::double_at ax_c;
  int32_t axi_c[VCLArchitecture::num_double];

  typename VCLArchitecture::double_vt vamp_a;
  typename VCLArchitecture::double_at aamp_a;
  typename VCLArchitecture::double_vt vamp_b;
  typename VCLArchitecture::double_at aamp_b;
  typename VCLArchitecture::double_vt vamp_c;
  typename VCLArchitecture::double_at aamp_c;

  // Initialze position within buffer
  ax_c[VCLArchitecture::num_double-1] = 0;

  while(true) {
    vx_a = dx * vcl_rng_a.exponential_double();
    vx_a.store(ax_a);

    ax_a[0] += ax_c[VCLArchitecture::num_double-1];
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

    vx_c = dx * vcl_rng_c.exponential_double();
    vx_c.store(ax_c);

    ax_c[0] += ax_b[VCLArchitecture::num_double-1];
    axi_c[0] = ax_c[0];
    __builtin_prefetch(pe_waveform_ + axi_c[0]*sizeof(double));
    for(unsigned i=1;i<VCLArchitecture::num_double;++i) {
      ax_c[i] += ax_c[i-1];
      axi_c[i] = ax_c[i];
      __builtin_prefetch(pe_waveform_ + axi_c[i]*sizeof(double));
    }

    vamp_a = nsb_pegen==nullptr? 1.0 : nsb_pegen->vcl_generate_amplitude(vcl_rng_a);
    vamp_a.store(aamp_a);

    vamp_b = nsb_pegen==nullptr? 1.0 : nsb_pegen->vcl_generate_amplitude(vcl_rng_b);
    vamp_b.store(aamp_b);

    vamp_c = nsb_pegen==nullptr? 1.0 : nsb_pegen->vcl_generate_amplitude(vcl_rng_c);
    vamp_c.store(aamp_c);

    for(unsigned i=0;i<VCLArchitecture::num_double;++i) {
      if(ax_a[i]<xmax) {
        pe_waveform_[axi_a[i]] += aamp_a[i];
      } else {
        goto break_to_outer_loop;
      }
    }

    for(unsigned i=0;i<VCLArchitecture::num_double;++i) {
      if(ax_b[i]<xmax) {
        pe_waveform_[axi_b[i]] += aamp_b[i];
      } else {
        goto break_to_outer_loop;
      }
    }

    for(unsigned i=0;i<VCLArchitecture::num_double;++i) {
      if(ax_c[i]<xmax) {
        pe_waveform_[axi_c[i]] += aamp_c[i];
      } else {
        goto break_to_outer_loop;
      }
    }
  }
break_to_outer_loop:
  if(ac_couple) {
    double mean_amp = nsb_pegen==nullptr? 1.0 : nsb_pegen->mean_amplitude();
    ac_coupling_constant_ += nsb_rate_ghz*trace_sampling_ns_*mean_amp;
  }
  pe_waveform_dft_valid_ = false;
}

template<typename VCLArchitecture> void WaveformProcessor::vcl_add_nsb(
  calin::math::rng::VCLRNG<VCLArchitecture>& vcl_rng_a,
  calin::math::rng::VCLRNG<VCLArchitecture>& vcl_rng_b,
  calin::math::rng::VCLRNG<VCLArchitecture>& vcl_rng_c,
  calin::math::rng::VCLRNG<VCLArchitecture>& vcl_rng_d,
  double nsb_rate_ghz,
  calin::simulation::detector_efficiency::SplinePEAmplitudeGenerator* nsb_pegen,
  bool ac_couple)
{
  const double dx = trace_sampling_inv_/nsb_rate_ghz;
  const double xmax = npixels_*trace_nsamples_;

  typename VCLArchitecture::double_vt vx_a;
  typename VCLArchitecture::double_at ax_a;
  int32_t axi_a[VCLArchitecture::num_double];

  typename VCLArchitecture::double_vt vx_b;
  typename VCLArchitecture::double_at ax_b;
  int32_t axi_b[VCLArchitecture::num_double];

  typename VCLArchitecture::double_vt vx_c;
  typename VCLArchitecture::double_at ax_c;
  int32_t axi_c[VCLArchitecture::num_double];

  typename VCLArchitecture::double_vt vx_d;
  typename VCLArchitecture::double_at ax_d;
  int32_t axi_d[VCLArchitecture::num_double];

  typename VCLArchitecture::double_vt vamp_a;
  typename VCLArchitecture::double_at aamp_a;
  typename VCLArchitecture::double_vt vamp_b;
  typename VCLArchitecture::double_at aamp_b;
  typename VCLArchitecture::double_vt vamp_c;
  typename VCLArchitecture::double_at aamp_c;
  typename VCLArchitecture::double_vt vamp_d;
  typename VCLArchitecture::double_at aamp_d;

  // Initialze position within buffer
  ax_d[VCLArchitecture::num_double-1] = 0;

  while(true) {
    vx_a = dx * vcl_rng_a.exponential_double();
    vx_a.store(ax_a);

    ax_a[0] += ax_d[VCLArchitecture::num_double-1];
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

    vx_c = dx * vcl_rng_c.exponential_double();
    vx_c.store(ax_c);

    ax_c[0] += ax_b[VCLArchitecture::num_double-1];
    axi_c[0] = ax_c[0];
    __builtin_prefetch(pe_waveform_ + axi_c[0]*sizeof(double));
    for(unsigned i=1;i<VCLArchitecture::num_double;++i) {
      ax_c[i] += ax_c[i-1];
      axi_c[i] = ax_c[i];
      __builtin_prefetch(pe_waveform_ + axi_c[i]*sizeof(double));
    }

    vx_d = dx * vcl_rng_d.exponential_double();
    vx_d.store(ax_d);

    ax_d[0] += ax_c[VCLArchitecture::num_double-1];
    axi_d[0] = ax_d[0];
    __builtin_prefetch(pe_waveform_ + axi_d[0]*sizeof(double));
    for(unsigned i=1;i<VCLArchitecture::num_double;++i) {
      ax_d[i] += ax_d[i-1];
      axi_d[i] = ax_d[i];
      __builtin_prefetch(pe_waveform_ + axi_d[i]*sizeof(double));
    }

    vamp_a = nsb_pegen==nullptr? 1.0 : nsb_pegen->vcl_generate_amplitude(vcl_rng_a);
    vamp_a.store(aamp_a);

    vamp_b = nsb_pegen==nullptr? 1.0 : nsb_pegen->vcl_generate_amplitude(vcl_rng_b);
    vamp_b.store(aamp_b);

    vamp_c = nsb_pegen==nullptr? 1.0 : nsb_pegen->vcl_generate_amplitude(vcl_rng_c);
    vamp_c.store(aamp_c);

    vamp_d = nsb_pegen==nullptr? 1.0 : nsb_pegen->vcl_generate_amplitude(vcl_rng_d);
    vamp_d.store(aamp_d);

    for(unsigned i=0;i<VCLArchitecture::num_double;++i) {
      if(ax_a[i]<xmax) {
        pe_waveform_[axi_a[i]] += aamp_a[i];
      } else {
        goto break_to_outer_loop;
      }
    }

    for(unsigned i=0;i<VCLArchitecture::num_double;++i) {
      if(ax_b[i]<xmax) {
        pe_waveform_[axi_b[i]] += aamp_b[i];
      } else {
        goto break_to_outer_loop;
      }
    }

    for(unsigned i=0;i<VCLArchitecture::num_double;++i) {
      if(ax_c[i]<xmax) {
        pe_waveform_[axi_c[i]] += aamp_c[i];
      } else {
        goto break_to_outer_loop;
      }
    }

    for(unsigned i=0;i<VCLArchitecture::num_double;++i) {
      if(ax_d[i]<xmax) {
        pe_waveform_[axi_d[i]] += aamp_d[i];
      } else {
        goto break_to_outer_loop;
      }
    }
  }
break_to_outer_loop:
  if(ac_couple) {
    double mean_amp = nsb_pegen==nullptr? 1.0 : nsb_pegen->mean_amplitude();
    ac_coupling_constant_ += nsb_rate_ghz*trace_sampling_ns_*mean_amp;
  }
  pe_waveform_dft_valid_ = false;
}

template<typename VCLArchitecture> void WaveformProcessor::vcl_add_electronics_noise(
  calin::math::rng::VCLRNG<VCLArchitecture>& vcl_rng,
  const double* noise_spectrum_amplitude, double scale)
{
  if(not el_waveform_dft_valid_) {
    throw std::runtime_error("add_electronics_noise : impulse response must be applied before noise added");
  }
  if(trace_nsamples_%(2*VCLArchitecture::num_double) != 0) {
    throw std::runtime_error("vcl_add_electronics_noise : number of samples must be multiple of "
      + std::to_string(2*VCLArchitecture::num_double));
  }
  double*__restrict__ buffer = el_waveform_dft_;
  scale *= std::sqrt(1.0/trace_nsamples_);
  for(unsigned ipixel=0; ipixel<npixels_; ++ipixel) {
    for(unsigned isample=0; isample<trace_nsamples_; isample+=2*VCLArchitecture::num_double) {
      typename VCLArchitecture::double_vt norm_1;
      typename VCLArchitecture::double_vt norm_2;
      vcl_rng.normal_two_double_bm(norm_1, norm_2);
      typename VCLArchitecture::double_vt x_1;
      typename VCLArchitecture::double_vt x_2;
      x_1.load(buffer + ipixel*trace_nsamples_ + isample);
      x_2.load(buffer + ipixel*trace_nsamples_ + isample + VCLArchitecture::num_double);
      typename VCLArchitecture::double_vt amp_1;
      typename VCLArchitecture::double_vt amp_2;
      amp_1.load(noise_spectrum_amplitude + isample);
      amp_2.load(noise_spectrum_amplitude + isample + VCLArchitecture::num_double);
      x_1 += scale * amp_1 * norm_1;
      x_2 += scale * amp_2 * norm_2;
      x_1.store(buffer + ipixel*trace_nsamples_ + isample);
      x_2.store(buffer + ipixel*trace_nsamples_ + isample + VCLArchitecture::num_double);
    }
  }
  el_waveform_valid_ = false;
}

template<typename VCLArchitecture> void WaveformProcessor::vcl_add_electronics_noise(
  calin::math::rng::VCLRNG<VCLArchitecture>& vcl_rng_a,
  calin::math::rng::VCLRNG<VCLArchitecture>& vcl_rng_b,
  const double* noise_spectrum_amplitude, double scale)
{
  if(not el_waveform_dft_valid_) {
    throw std::runtime_error("add_electronics_noise : impulse response must be applied before noise added");
  }
  if(trace_nsamples_%(4*VCLArchitecture::num_double) != 0) {
    throw std::runtime_error("vcl_add_electronics_noise : number of samples must be multiple of "
      + std::to_string(4*VCLArchitecture::num_double));
  }
  double*__restrict__ buffer = el_waveform_dft_;
  scale *= std::sqrt(1.0/trace_nsamples_);
  for(unsigned ipixel=0; ipixel<npixels_; ++ipixel) {
    for(unsigned isample=0; isample<trace_nsamples_; isample+=4*VCLArchitecture::num_double) {
      typename VCLArchitecture::double_vt norm_a1;
      typename VCLArchitecture::double_vt norm_a2;
      typename VCLArchitecture::double_vt norm_b1;
      typename VCLArchitecture::double_vt norm_b2;
      vcl_rng_a.normal_two_double_bm(norm_a1, norm_a2);
      vcl_rng_b.normal_two_double_bm(norm_b1, norm_b2);

      typename VCLArchitecture::double_vt x_a1;
      typename VCLArchitecture::double_vt x_a2;
      typename VCLArchitecture::double_vt x_b1;
      typename VCLArchitecture::double_vt x_b2;
      x_a1.load(buffer + ipixel*trace_nsamples_ + isample);
      x_a2.load(buffer + ipixel*trace_nsamples_ + isample + VCLArchitecture::num_double);
      x_b1.load(buffer + ipixel*trace_nsamples_ + isample + 2*VCLArchitecture::num_double);
      x_b2.load(buffer + ipixel*trace_nsamples_ + isample + 3*VCLArchitecture::num_double);

      typename VCLArchitecture::double_vt amp_a1;
      typename VCLArchitecture::double_vt amp_a2;
      typename VCLArchitecture::double_vt amp_b1;
      typename VCLArchitecture::double_vt amp_b2;
      amp_a1.load(noise_spectrum_amplitude + isample);
      amp_a2.load(noise_spectrum_amplitude + isample + VCLArchitecture::num_double);
      amp_b1.load(noise_spectrum_amplitude + isample + 2*VCLArchitecture::num_double);
      amp_b2.load(noise_spectrum_amplitude + isample + 3*VCLArchitecture::num_double);

      x_a1 += scale * amp_a1 * norm_a1;
      x_a2 += scale * amp_a2 * norm_a2;
      x_b1 += scale * amp_b1 * norm_b1;
      x_b2 += scale * amp_b2 * norm_b2;

      x_a1.store(buffer + ipixel*trace_nsamples_ + isample);
      x_a2.store(buffer + ipixel*trace_nsamples_ + isample + VCLArchitecture::num_double);
      x_b1.store(buffer + ipixel*trace_nsamples_ + isample + 2*VCLArchitecture::num_double);
      x_b2.store(buffer + ipixel*trace_nsamples_ + isample + 3*VCLArchitecture::num_double);
    }
  }
  el_waveform_valid_ = false;
}

#endif
