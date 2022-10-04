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
  if(trace_nsamples_%VCLArchitecture::num_double != 0) {
    throw std::runtime_error("vcl_add_electronics_noise : number of samples must be multiple of "
      + std::to_string(VCLArchitecture::num_double));
  }
  double*__restrict__ buffer = el_waveform_dft_;
  scale *= std::sqrt(1.0/trace_nsamples_);
  for(unsigned ipixel=0; ipixel<npixels_; ++ipixel) {
    for(unsigned isample=0; isample<trace_nsamples_; isample+=VCLArchitecture::num_double) {
      typename VCLArchitecture::double_vt norm;
      norm = vcl_rng.normal_double();
      typename VCLArchitecture::double_vt x;
      x.load(buffer + ipixel*trace_nsamples_ + isample);
      typename VCLArchitecture::double_vt amp;
      amp.load(noise_spectrum_amplitude + isample);
      x += scale * amp * norm;
      x.store(buffer + ipixel*trace_nsamples_ + isample);
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
  if(trace_nsamples_%(2*VCLArchitecture::num_double) != 0) {
    throw std::runtime_error("vcl_add_electronics_noise : number of samples must be multiple of "
      + std::to_string(2*VCLArchitecture::num_double));
  }
  double*__restrict__ buffer = el_waveform_dft_;
  scale *= std::sqrt(1.0/trace_nsamples_);
  unsigned mpixels = npixels_>>1;
  for(unsigned ipixel=0; ipixel<mpixels; ++ipixel) {
    for(unsigned isample=0; isample<trace_nsamples_; isample+=VCLArchitecture::num_double) {
      typename VCLArchitecture::double_vt norm_a;
      typename VCLArchitecture::double_vt norm_b;
      norm_a = vcl_rng_a.normal_double();
      norm_b = vcl_rng_b.normal_double();

      typename VCLArchitecture::double_vt x_a;
      typename VCLArchitecture::double_vt x_b;
      x_a.load(buffer + ipixel*trace_nsamples_ + isample);
      x_b.load(buffer + (ipixel+mpixels)*trace_nsamples_ + isample);

      typename VCLArchitecture::double_vt amp;
      amp.load(noise_spectrum_amplitude + isample);
      amp *= scale;

      x_a += amp * norm_a;
      x_b += amp * norm_b;

      x_a.store(buffer + ipixel*trace_nsamples_ + isample);
      x_b.store(buffer + (ipixel+mpixels)*trace_nsamples_ + isample);
    }
  }
  for(unsigned ipixel=mpixels<<1; ipixel<npixels_; ++ipixel) {
    for(unsigned isample=0; isample<trace_nsamples_; isample+=2*VCLArchitecture::num_double) {
      typename VCLArchitecture::double_vt norm_a;
      typename VCLArchitecture::double_vt norm_b;
      norm_a = vcl_rng_a.normal_double();
      norm_b = vcl_rng_b.normal_double();

      typename VCLArchitecture::double_vt x_a;
      typename VCLArchitecture::double_vt x_b;
      x_a.load(buffer + ipixel*trace_nsamples_ + isample);
      x_b.load(buffer + ipixel*trace_nsamples_ + isample + VCLArchitecture::num_double);

      typename VCLArchitecture::double_vt amp_a;
      typename VCLArchitecture::double_vt amp_b;
      amp_a.load(noise_spectrum_amplitude + isample);
      amp_b.load(noise_spectrum_amplitude + isample + VCLArchitecture::num_double);

      x_a += scale * amp_a * norm_a;
      x_b += scale * amp_b * norm_b;

      x_a.store(buffer + ipixel*trace_nsamples_ + isample);
      x_b.store(buffer + ipixel*trace_nsamples_ + isample + VCLArchitecture::num_double);
    }
  }
  el_waveform_valid_ = false;
}

template<typename VCLArchitecture> int WaveformProcessor::vcl_digital_multiplicity_trigger_alt(
  double threshold,
  unsigned time_over_threshold_samples, unsigned coherence_time_samples,
  unsigned multiplicity_threshold, bool loud)
{
  compute_el_waveform();
  uint16_t* multiplicity = static_cast<uint16_t*>(alloca(trace_nsamples_ * sizeof(uint16_t)));
  std::fill(multiplicity, multiplicity+trace_nsamples_, 0);
  for(unsigned ipixel=0; ipixel<npixels_; ++ipixel) {
    double*__restrict__ pixel_waveform = el_waveform_ + ipixel*trace_nsamples_;
    unsigned l0_tot = 0;
    unsigned l0_eop = 0;
    for(unsigned isamp=0; isamp<trace_nsamples_; isamp += 4*VCLArchitecture::num_double) {
      typename VCLArchitecture::double_vt samples_a;
      samples_a.load(pixel_waveform + isamp);

      typename VCLArchitecture::double_vt samples_b;
      samples_b.load(pixel_waveform + isamp + VCLArchitecture::num_double);

      typename VCLArchitecture::double_vt samples_c;
      samples_c.load(pixel_waveform + isamp + 2*VCLArchitecture::num_double);

      typename VCLArchitecture::double_vt samples_d;
      samples_d.load(pixel_waveform + isamp + 3*VCLArchitecture::num_double);

      uint32_t above_threshold =
        static_cast<uint32_t>(vcl::to_bits(samples_a > threshold)) |
        (static_cast<uint32_t>(vcl::to_bits(samples_b > threshold)) << VCLArchitecture::num_double) |
        (static_cast<uint32_t>(vcl::to_bits(samples_c > threshold)) << (2*VCLArchitecture::num_double)) |
        (static_cast<uint32_t>(vcl::to_bits(samples_d > threshold)) << (3*VCLArchitecture::num_double));
      uint32_t triggered = 0;
      if(l0_eop > isamp) {
        triggered |= (1<<std::min(l0_eop-isamp, 4*VCLArchitecture::num_double))-1;
      }

      if(above_threshold == 0) {
        l0_tot = 0;
      } else {
        unsigned jsamp = 0;
        uint32_t value = above_threshold & 0x1;
        while(jsamp < 4*VCLArchitecture::num_double) {
          if(value == 0) {
            l0_tot = 0;
            uint32_t ksamp = ffs(above_threshold);
            ksamp = std::min(ksamp-1, 4*VCLArchitecture::num_double-jsamp);
            above_threshold >>= ksamp;
            jsamp += ksamp;
            value = 0x01;
          } else { /* value == 1 */
            uint32_t ksamp = ffs(~above_threshold);
            ksamp = std::min(ksamp-1, 4*VCLArchitecture::num_double-jsamp);
            l0_tot += ksamp;
            if(l0_tot >= time_over_threshold_samples) {
              unsigned new_l0_eop = isamp+jsamp+ksamp+coherence_time_samples-1;
              if(l0_eop <= isamp+4*VCLArchitecture::num_double) {
                unsigned ksop = time_over_threshold_samples - (l0_tot - ksamp) - 1;
                triggered |= ((1<<(std::min(new_l0_eop-isamp, 4*VCLArchitecture::num_double)-jsamp-ksop))-1)<<(jsamp+ksop);
              }
              l0_eop = new_l0_eop;
            }
            above_threshold >>= ksamp;
            jsamp += ksamp;
            value = 0x00;
          }
        }
      }
      if(triggered) {
        typename VCLArchitecture::uint16_vt mult;
        mult.load(multiplicity + isamp);
        typename VCLArchitecture::uint16_bvt mult_add_mask;
        mult_add_mask.load_bits(triggered);
        mult = vcl::if_add(mult_add_mask, mult, 1);
        mult.store(multiplicity + isamp);
      }
    }
  }
  if(loud) {
    for(unsigned isamp=0; isamp<trace_nsamples_; ++isamp) {
      calin::util::log::LOG(calin::util::log::INFO) << isamp << ' ' << multiplicity[isamp];
    }
  }
  for(unsigned isamp=0; isamp<trace_nsamples_; ++isamp) {
    if(multiplicity[isamp] >= multiplicity_threshold) {
      return isamp;
    }
  }
  return -1;
}

template<typename VCLArchitecture> int WaveformProcessor::vcl_digital_nn_trigger_alt(
  double threshold,
  unsigned time_over_threshold_samples, unsigned coherence_time_samples,
  unsigned multiplicity_threshold, WaveformProcessorTriggerMemoryBuffers* buffer)
{
  compute_el_waveform();

  std::unique_ptr<WaveformProcessorTriggerMemoryBuffers> my_buffer;
  if(buffer == nullptr) {
    buffer = new WaveformProcessorTriggerMemoryBuffers(npixels_, trace_nsamples_);
    my_buffer.reset(buffer);
  }
  buffer->clear();

  for(unsigned ipixel=0; ipixel<npixels_; ++ipixel) {
    double*__restrict__ pixel_waveform = el_waveform_ + ipixel*trace_nsamples_;
    unsigned l0_tot = 0;
    unsigned l0_eop = 0;
    uint32_t triggered_32;
    uint32_t newly_triggered;
    for(unsigned isamp=0; isamp<trace_nsamples_; isamp += 4*VCLArchitecture::num_double) {
      if(isamp % 32 == 0) {
        newly_triggered = buffer->newly_triggered_bitmask[isamp/32];
        triggered_32 = 0;
      }
      typename VCLArchitecture::double_vt samples_a;
      samples_a.load(pixel_waveform + isamp);

      typename VCLArchitecture::double_vt samples_b;
      samples_b.load(pixel_waveform + isamp + VCLArchitecture::num_double);

      typename VCLArchitecture::double_vt samples_c;
      samples_c.load(pixel_waveform + isamp + 2*VCLArchitecture::num_double);

      typename VCLArchitecture::double_vt samples_d;
      samples_d.load(pixel_waveform + isamp + 3*VCLArchitecture::num_double);

      uint32_t above_threshold =
        static_cast<uint32_t>(vcl::to_bits(samples_a > threshold)) |
        (static_cast<uint32_t>(vcl::to_bits(samples_b > threshold)) << VCLArchitecture::num_double) |
        (static_cast<uint32_t>(vcl::to_bits(samples_c > threshold)) << (2*VCLArchitecture::num_double)) |
        (static_cast<uint32_t>(vcl::to_bits(samples_d > threshold)) << (3*VCLArchitecture::num_double));
      uint32_t triggered = 0;
      if(l0_eop > isamp) {
        triggered |= (1<<std::min(l0_eop-isamp, 4*VCLArchitecture::num_double))-1;
      }

      if(above_threshold == 0) {
        l0_tot = 0;
      } else {
        unsigned jsamp = 0;
        uint32_t value = above_threshold & 0x1;
        while(jsamp < 4*VCLArchitecture::num_double) {
          if(value == 0) {
            l0_tot = 0;
            uint32_t ksamp = ffs(above_threshold);
            ksamp = std::min(ksamp-1, 4*VCLArchitecture::num_double-jsamp);
            above_threshold >>= ksamp;
            jsamp += ksamp;
            value = 0x01;
          } else { /* value == 1 */
            uint32_t ksamp = ffs(~above_threshold);
            ksamp = std::min(ksamp-1, 4*VCLArchitecture::num_double-jsamp);
            l0_tot += ksamp;
            if(l0_tot >= time_over_threshold_samples) {
              unsigned new_l0_eop = isamp+jsamp+ksamp+coherence_time_samples-1;
              if(l0_eop <= isamp+4*VCLArchitecture::num_double) {
                unsigned ksop = time_over_threshold_samples - (l0_tot - ksamp) - 1;
                triggered |= ((1<<(std::min(new_l0_eop-isamp, 4*VCLArchitecture::num_double)-jsamp-ksop))-1)<<(jsamp+ksop);
                newly_triggered |= 0x1<<((isamp+jsamp+ksop)%32);
              }
              l0_eop = new_l0_eop;
            }
            above_threshold >>= ksamp;
            jsamp += ksamp;
            value = 0x00;
          }
        }
      }
      if(triggered) {
        typename VCLArchitecture::uint16_vt mult;
        mult.load(buffer->multiplicity + isamp);
        typename VCLArchitecture::uint16_bvt mult_add_mask;
        mult_add_mask.load_bits(triggered);
        mult = vcl::if_add(mult_add_mask, mult, 1);
        mult.store(buffer->multiplicity + isamp);
      }

      triggered_32 |= triggered << (isamp%32);
      if((isamp+4*VCLArchitecture::num_double) % 32 == 0) {
        buffer->newly_triggered_bitmask[isamp/32] = newly_triggered;
        buffer->triggered_bitmask[(ipixel*trace_nsamples_ + isamp)/32] = triggered_32;
      }
    }
  }

  for(unsigned isamp=0; isamp<trace_nsamples_; ++isamp) {
    bool found_new_l0_triggers = buffer->is_newly_triggered(isamp);
    if(buffer->multiplicity[isamp] >= multiplicity_threshold and found_new_l0_triggers) {
      // The simple multiplicity threshold has been met, and we have some newly
      // triggered channels - test neighbours
      switch(multiplicity_threshold) {
      case 0:
      case 1:
        return isamp;
      case 2:
        for(unsigned ipixel=0; ipixel<npixels_; ++ipixel) {
          if(buffer->is_triggered(ipixel, isamp)) {
            for(unsigned ineighbour=0; ineighbour<max_num_neighbours_; ++ineighbour) {
              int jpixel = neighbour_map_[ipixel*max_num_neighbours_ + ineighbour];
              if(jpixel>0 and buffer->is_triggered(jpixel,isamp)) {
                return isamp;
              }
            }
          }
        }
        break;
      case 3:
        for(unsigned ipixel=0; ipixel<npixels_; ++ipixel) {
          unsigned nneighbour = 1;
          if(buffer->is_triggered(ipixel, isamp)) {
            for(unsigned ineighbour=0; ineighbour<max_num_neighbours_; ++ineighbour) {
              int jpixel = neighbour_map_[ipixel*max_num_neighbours_ + ineighbour];
              if(jpixel>0 and buffer->is_triggered(jpixel,isamp)) {
                ++nneighbour;
              }
            }
            if(nneighbour >= multiplicity_threshold) {
              return isamp;
            }
          }
        }
        break;
      default:
        throw std::runtime_error("digital_nn_trigger : multiplicity "
          + std::to_string(multiplicity_threshold) + " unsupported");
      }
    }
  }

  return -1;
}

#endif
