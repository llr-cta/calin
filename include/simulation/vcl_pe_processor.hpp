/*

   calin/simulation/vcl_pe_processor.hpp -- Stephen Fegan -- 2025-10-08

   VCL class for processing PEs into waveforms: convolvolution with impulse
   response, trigger algorithms, and others

   Copyright 2025, Stephen Fegan <sfegan@llr.in2p3.fr>
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

// ===============================
// NO USER SERVICABLE PARTS INSIDE
// ===============================

#pragma once

#include <algorithm>
#include <limits>
#include <bit>
#include <fftw3.h>

#include <util/log.hpp>
#include <util/vcl.hpp>
#include <util/memory.hpp>
#include <util/string.hpp>
#include <math/special.hpp>
#include <math/ray_vcl.hpp>
#include <math/rng_vcl.hpp>
#include <math/geometry_vcl.hpp>
#include <math/fftw_util.hpp>
#include <math/interpolation_1d.hpp>
#include <math/brent.hpp>
#include <simulation/pe_processor.hpp>
#include <simulation/vcl_detector_efficiency.hpp>
#include <iact_data/instrument_layout.pb.h>
#include <iact_data/instrument_layout.hpp>

namespace calin { namespace simulation { namespace vcl_pe_processor {

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
// VCLWaveformPEProcessor
//
// A PE processor that uses VCL types and operations to process PEs into a
// binned waveform, convolved with a unit response function, and with noise
//
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template<typename VCLReal> class alignas(VCLReal::vec_bytes) VCLWaveformPEProcessor:
  public calin::simulation::pe_processor::SimpleListPEProcessor
{
public:
  CALIN_TYPEALIAS(real_t,    typename VCLReal::real_t);
  CALIN_TYPEALIAS(real_vt,   typename VCLReal::real_vt);
  CALIN_TYPEALIAS(real_bvt,  typename VCLReal::real_bvt);
  CALIN_TYPEALIAS(real_at,   typename VCLReal::real_at);
  CALIN_TYPEALIAS(int_t,     typename VCLReal::int_t);
  CALIN_TYPEALIAS(int_vt,    typename VCLReal::int_vt);
  CALIN_TYPEALIAS(int_at,    typename VCLReal::int_at);
  CALIN_TYPEALIAS(int_bvt,   typename VCLReal::int_bvt);
  CALIN_TYPEALIAS(uint_t,    typename VCLReal::uint_t);
  CALIN_TYPEALIAS(uint_vt,   typename VCLReal::uint_vt);
  CALIN_TYPEALIAS(uint_bvt,  typename VCLReal::uint_bvt);
  CALIN_TYPEALIAS(vecX_t,    typename VCLReal::vecX_t);
  CALIN_TYPEALIAS(matX_t,    typename VCLReal::matX_t);

  CALIN_TYPEALIAS(VCLArchitecture, typename VCLReal::architecture);

  CALIN_TYPEALIAS(RNG, calin::math::rng::VCLRNG<VCLArchitecture>);

  VCLWaveformPEProcessor(unsigned nscope, unsigned npix, unsigned nsample, 
      double time_resolution_ns, unsigned nsample_advance,
      RNG* rng = nullptr, bool auto_clear = true, bool adopt_rng = false):
    calin::simulation::pe_processor::SimpleListPEProcessor(nscope, npix, auto_clear),
    nsample_(round_nreal_to_vector(nsample)),  nadvance_(nsample_advance),
    time_resolution_ns_(time_resolution_ns), sampling_freq_ghz_(1.0/time_resolution_ns_),
    time_advance_(double(nadvance_)*time_resolution_ns_),
    pe_waveform_(nsample_, round_nreal_to_vector(npix)), 
    v_waveform_(nsample_, round_nreal_to_vector(npix)), 
    rng_(rng), adopt_rng_(adopt_rng)
  {
    if(nsample_ % VCLReal::num_real != 0) {
      throw std::domain_error("Number of samples " + std::to_string(nsample_) 
        + " must be a multiple of vector size " 
        + std::to_string(VCLReal::num_real));
    }

    if((uintptr_t)(const void *)pe_waveform_.data() % VCLReal::vec_bytes != 0) {
      calin::util::log::LOG(calin::util::log::WARNING) 
        << "pe_waveform_ not aligned on " << VCLReal::vec_bytes << " boundary.";
    }

    if((uintptr_t)(const void *)v_waveform_.data() % VCLReal::vec_bytes != 0) {
      calin::util::log::LOG(calin::util::log::WARNING) 
        << "v_waveform_ not aligned on " << VCLReal::vec_bytes << " boundary.";
    }

    if(rng_ == nullptr) {
      rng_ = new RNG(__PRETTY_FUNCTION__, "VCLWaveformPEProcessor RNG");
      adopt_rng_ = true;
    }

    clear_waveforms();
  }

  virtual ~VCLWaveformPEProcessor() 
  {
    if(adopt_rng_) {
      delete rng_;
    }
  }

  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  //
  //    db       .d88b.   .d8b.  d8888b. 
  //    88      .8P  Y8. d8' `8b 88  `8D 
  //    88      88    88 88ooo88 88   88 
  //    88      88    88 88~~~88 88   88 
  //    88booo. `8b  d8' 88   88 88  .8D 
  //    Y88888P  `Y88P'  YP   YP Y8888D' 
  //
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  void load_from_simulated_event_with_pmt_noise(
    const calin::ix::simulation::simulated_event::DetectorGroupEvent& detector_group_event,
    calin::simulation::detector_efficiency::VCLPEAmplitudeGenerator<VCLArchitecture>* pegen,
    double time_spread = 0)
  {
    this->start_processing();

    constexpr unsigned num_double = VCLArchitecture::num_double;
    typename VCLArchitecture::double_at pe_dt;
    typename VCLArchitecture::double_at pe_w;
    std::fill(pe_w, pe_w+num_double, 1.0);
    std::fill(pe_dt, pe_dt+num_double, 0.0);
    unsigned idev = 0;

    for(int idetector=0; idetector<detector_group_event.detector_size(); ++idetector) {
      const auto& detector_event = detector_group_event.detector(idetector);
      unsigned iscope = detector_event.detector_id();
      double ref_time = detector_event.reference_time();
      for(int jpixel=0; jpixel<detector_event.pixel_size(); ++jpixel) {
        const auto& pixel_event = detector_event.pixel(jpixel);
        int ipix = pixel_event.pixel_id();
        for(int kpe=0; kpe<pixel_event.weight_size(); ++kpe) {
          if(idev == 0) {
            pegen->vcl_generate_amplitude(*rng_).store_a(pe_w);
            if(time_spread > 0) {
              (rng_->normal_double() * time_spread).store_a(pe_dt);
            } else if (time_spread < 0) {
              rng_->uniform_double_zc(time_spread).store_a(pe_dt);
            }
          }
          double t = pixel_event.time(kpe) + ref_time + pe_dt[idev];
          double w = pixel_event.weight(kpe) * pe_w[idev];
          this->process_focal_plane_hit(iscope, ipix, 0.0, 0.0, 0.0, 0.0, t, w);
          idev = (idev+1) % num_double;
        }
        for(int kpe=pixel_event.weight_size(); kpe<pixel_event.time_size(); ++kpe) {
          if(idev == 0) {
            pegen->vcl_generate_amplitude(*rng_).store_a(pe_w);
            if(time_spread > 0) {
              (rng_->normal_double() * time_spread).store_a(pe_dt);
            } else if (time_spread < 0) {
              rng_->uniform_double_zc(time_spread).store_a(pe_dt);
            }
          }
          double t = pixel_event.time(kpe) + ref_time + pe_dt[idev];
          this->process_focal_plane_hit(iscope, ipix, 0.0, 0.0, 0.0, 0.0, t, pe_w[idev]);
          idev = (idev+1) % num_double;
        }
        for(int kpe=0; kpe<pixel_event.integer_time_size(); ++kpe) {
          if(idev == 0) {
            pegen->vcl_generate_amplitude(*rng_).store_a(pe_w);
            if(time_spread > 0) {
              (rng_->normal_double() * time_spread).store_a(pe_dt);
            } else if (time_spread < 0) {
              rng_->uniform_double_zc(time_spread).store_a(pe_dt);
            }
          }
          double t = pixel_event.integer_time(kpe)*detector_event.integer_time_resolution() + ref_time + pe_dt[idev];
          this->process_focal_plane_hit(iscope, ipix, 0.0, 0.0, 0.0, 0.0, t, pe_w[idev]);
          idev = (idev+1) % num_double;
        }
      }
    }
    this->finish_processing();
  }

  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  //
  //    d8888b. d88888b .d8888. 
  //    88  `8D 88'     88'  YP 
  //    88oodD' 88ooooo `8bo.   
  //    88~~~   88~~~~~   `Y8b. 
  //    88      88.     db   8D 
  //    88      Y88888P `8888Y' 
  //
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  // This function transfers PEs from the list into the waveform array. It can be used
  // with any instance of SimpleListPEProcessor, not just this one. This allows
  // the electronics simulation to be separated from the PE generation.
  
  void transfer_scope_pes_to_waveform(const calin::simulation::pe_processor::SimpleListPEProcessor& pe_list, unsigned iscope, 
    const Eigen::VectorXd& channel_time_offset_ns = Eigen::VectorXd())
  {
    // Non-vectorized version. Takes advantage of "approximate" ordering of PEs
    // in time within each pixel to minimize memory accesses. Vectorized version
    // was more complex and was tested to be not much faster.
    
    if(channel_time_offset_ns.size()!=0 and channel_time_offset_ns.size()!=npix_) {
      throw std::domain_error("Time offset vector length is not equal to number of pixels " 
        + std::to_string(channel_time_offset_ns.size()) + " != " + std::to_string(npix_));
    }

    // Find the PE time to use as a reference for binning. Use the minimum PE time if
    // the time range is small (<8ns), othewise use the time of the 8ns bin in which
    // most PEs arrive. Later a full trigger simulation will calculate a better time.

    t0_ = pe_list.tmin(iscope);
    constexpr double max_jitter = 8.0;
    constexpr double inv_jitter = 1.0 / max_jitter;
    if(pe_list.tmax(iscope) - t0_ > max_jitter) {
      unsigned nbin = std::min(unsigned((pe_list.tmax(iscope) - t0_)*inv_jitter) + 1u, 128u);
      double* bins = static_cast<double*>(alloca(nbin*sizeof(double)));
      std::fill_n(bins, nbin, 0);
      for(unsigned ipix=0; ipix<npix_; ++ipix) {
        const double* t_ptr = nullptr;
        const double* w_ptr = nullptr;
        unsigned npe = pe_list.pe_ptrs(iscope, ipix, t_ptr, w_ptr);
        for(unsigned ipe=0; ipe<npe; ++ipe) {
          int itime = int(floor((t_ptr[ipe] - t0_)*inv_jitter));
          if(itime >= 0 and itime < int(nbin)) {
            bins[itime] += w_ptr[ipe];
          }
        }
      }
      auto itmax = std::max_element(bins, bins+nbin)-bins;
      auto it0 = itmax;
      while(it0 > 0 and bins[it0-1] > std::min(100.0,0.5*bins[itmax])) {
        --it0;
      }
      t0_ += it0*max_jitter;
    }
    t0_ -= time_advance_;

    for(unsigned ipix=0; ipix<npix_; ++ipix) {
      double t0ipix = t0_;
      if(channel_time_offset_ns.size()) {
        t0ipix -= channel_time_offset_ns[ipix];
      }
      real_t *__restrict__ pe_waveform_ptr_ = &pe_waveform_(0, ipix);
      const double* t_ptr = nullptr;
      const double* w_ptr = nullptr;
      unsigned npe = pe_list.pe_ptrs(iscope, ipix, t_ptr, w_ptr);
      if(t_ptr==nullptr) {
        continue;
      }
      int it = int(floor((t_ptr[0] - t0ipix) * sampling_freq_ghz_));
      double wt = w_ptr[0];
      for(unsigned ipe=1; ipe<npe; ++ipe) {
        int jt = int(floor((t_ptr[ipe] - t0ipix) * sampling_freq_ghz_));
        if(it==jt) {
          wt += w_ptr[ipe];
        } else {
          if(it>=0 and it<int(nsample_)) {
            pe_waveform_ptr_[it] += wt;
          }
          it = jt;
          wt = w_ptr[ipe];
        }
      }
      if(it>=0 and it<int(nsample_)) {
        pe_waveform_ptr_[it] += wt;
      }
    }
  }

  void transfer_scope_pes_to_waveform(unsigned iscope, const Eigen::VectorXd& channel_time_offset_ns = Eigen::VectorXd())
  {
    return transfer_scope_pes_to_waveform(*this, iscope, channel_time_offset_ns);
  }

  vecX_t estimate_pes_in_window(unsigned window_size, unsigned window_start=0)
  {
    if(window_size == 0) {
      throw std::out_of_range("Window size must be non-zero");
    }
    if(window_size > nsample_) {
      throw std::out_of_range("Window size longer than samples array: " + 
        std::to_string(window_start) + ">" + std::to_string(nsample_));
    }
    if(window_size+window_start > nsample_) {
      throw std::out_of_range("End of window exceeds samples array: " + 
        std::to_string(window_size+window_start) + ">" + std::to_string(nsample_));
    }
    if(nsample_ % VCLReal::num_real != 0) {
      throw std::domain_error("Number of samples " + std::to_string(nsample_) 
        + " is not a multiple of vector size " 
        + std::to_string(VCLReal::num_real));
    }

    vecX_t pes_per_channel_in_window(npix_);

    unsigned isample0 = window_start - window_start%VCLReal::num_real;
    real_vt* a_vec = calin::util::memory::aligned_calloc<real_vt>(nsample_);
    for(unsigned ipix=0; ipix<npix_; ipix+=VCLReal::num_real) {
      // Load data from "pe_waveform_" into "a_vec", block by block, transposing as we go along
      for(unsigned isample=isample0; isample<nsample_; isample+=VCLReal::num_real) {
        real_vt block[VCLReal::num_real]; // square matrix of doubles
        for(unsigned jpix=0, mpix=std::min(npix_-ipix,VCLReal::num_real); jpix<mpix; jpix++) {
          block[jpix].load_a(pe_waveform_.data() + (ipix+jpix)*nsample_ + isample);
        }
        calin::util::vcl::transpose(block);
        for(unsigned jsample = 0; jsample<VCLReal::num_real; jsample++) {
          a_vec[isample+jsample] = block[jsample];
        }
      }

      real_vt max_integral = 0;
      real_vt integral = 0;

      real_vt *__restrict__ iptr = a_vec + window_start;
      real_vt *__restrict__ jptr = iptr;
      real_vt *__restrict__ zptr = iptr + window_size;
      while(iptr < zptr) {
        integral += *(iptr++);
      }
      zptr += nsample_ - window_size - window_start;
      max_integral = integral;
      while(iptr < zptr) {
        integral += *(iptr++);
        integral -= *(jptr++);
        max_integral = vcl::max(max_integral, integral);
      }

      real_at max_integral_array;
      max_integral.store(max_integral_array);
      for(unsigned jpix=0, mpix=std::min(npix_-ipix,VCLReal::num_real); jpix<mpix; jpix++) {
        pes_per_channel_in_window(ipix+jpix) = max_integral_array[jpix];
      }
    }
    ::free(a_vec);
    return pes_per_channel_in_window;
  }

  double t0() const {
    return t0_;
  } 

  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  //
  //    d8b   db .d8888. d8888b. 
  //    888o  88 88'  YP 88  `8D 
  //    88V8o 88 `8bo.   88oooY' 
  //    88 V8o88   `Y8b. 88~~~b. 
  //    88  V888 db   8D 88   8D 
  //    VP   V8P `8888Y' Y8888P' 
  //
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  void add_nsb_noise_to_waveform(const Eigen::VectorXd& nsb_freq_per_pixel_ghz,
    calin::simulation::detector_efficiency::VCLPEAmplitudeGenerator<VCLArchitecture>* pegen = nullptr,
    double t0_samples = 0)
  {
    if(nsb_freq_per_pixel_ghz.size() != int(npix_)) {
      throw std::domain_error("NSB frequency vector length is not equal to number of pixels " 
        + std::to_string(nsb_freq_per_pixel_ghz.size()) + " != " + std::to_string(npix_));
    }

    t0_samples = std::max(t0_samples, 0.0);

    // Non-vectorized version - a vectorized version would have to use scatter/gather
    // possibly not worth it
    double t_max = double(nsample_);
    for(unsigned ipix=0; ipix<npix_; ++ipix) {
      double rate_samples = sampling_freq_ghz_/nsb_freq_per_pixel_ghz[ipix];
      if(rate_samples<=0) {
        continue;
      }

      real_t *__restrict__ pe_waveform_ptr_ = &pe_waveform_(0, ipix);
      double t_samples = t0_samples;
      while(t_samples < t_max) {
        typename VCLArchitecture::double_vt dt_samples = rng_->exponential_double() * rate_samples;
        typename VCLArchitecture::double_vt charge = 1.0;
        if(pegen) {
          charge = pegen->vcl_generate_amplitude(*rng_);
        }
        typename VCLArchitecture::double_at dt_samples_a;
        dt_samples.store(dt_samples_a);
        typename VCLArchitecture::double_at charge_a;
        charge.store(charge_a);
        for(unsigned insb=0; insb<VCLArchitecture::num_double; ++insb) {
          t_samples += dt_samples_a[insb];
          if(t_samples >= t_max) {
            goto next_pixel;
          }
          int it = int(floor(t_samples));
          pe_waveform_ptr_[it] += charge_a[insb];
        }
      }
      next_pixel:
      ;
    }
  }

  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  //
  //    d888888b .88b  d88. d8888b. db    db db      .d8888. d88888b 
  //      `88'   88'YbdP`88 88  `8D 88    88 88      88'  YP 88'     
  //       88    88  88  88 88oodD' 88    88 88      `8bo.   88ooooo 
  //       88    88  88  88 88~~~   88    88 88        `Y8b. 88~~~~~ 
  //      .88.   88  88  88 88      88b  d88 88booo. db   8D 88.     
  //    Y888888P YP  YP  YP 88      ~Y8888P' Y88888P `8888Y' Y88888P 
  //
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  unsigned register_impulse_response(const vecX_t& impulse_response, 
    const std::string& units, double window_fraction=0.6)
  {
    ImpulseResponse ir;

    ir.response.resize(nsample_);
    ir.transform.resize(nsample_);

    auto plan = calin::math::fftw_util::plan_r2hc(nsample_, &ir.response[0], &ir.transform[0], FFTW_MEASURE);

    // 1. Store impulse response, truncating or extending if necessary
    ir.response_size = impulse_response.size();    
    if(ir.response_size>nsample_) {
      calin::util::log::LOG(calin::util::log::WARNING) 
        << "Truncating impulse response of size "
        << impulse_response.size() << " to " << nsample_;
      ir.response_size = nsample_;
      ir.response = impulse_response.head(nsample_);
    } else {
      ir.response.setZero();
      ir.response.head(ir.response_size) = impulse_response;
    }

    // 2. Calculate peak, integral and percentiles
    std::vector<double> x(ir.response_size);
    std::vector<double> y(ir.response_size);
    for(unsigned i=0; i<ir.response_size; ++i) {
      double ri = impulse_response(i);
      ir.response_total_integral += ri;
      x[i] = i;
      y[i] = ir.response_total_integral;
      if(ri > ir.response_peak_value) {
        ir.response_peak_value = ri;
        ir.response_peak_index = i;
      }
      if(ir.response_total_integral > ir.response_integral_peak_value) {
        ir.response_integral_peak_value = ir.response_total_integral;
        ir.response_integral_peak_index = i;
      }
    }

    calin::math::interpolation_1d::InterpLinear1D f_of_x(x,y);
    ir.response_window_frac = window_fraction;
    double target = ir.response_integral_peak_value * (1.0-window_fraction)/2.0;
    ir.response_window_lo = 
      calin::math::brent::brent_zero(0, ir.response_integral_peak_index, 
        [&f_of_x,target](double x) { return f_of_x(x)-target; });
    target = ir.response_integral_peak_value * (1.0+window_fraction)/2.0;
    ir.response_window_hi = 
      calin::math::brent::brent_zero(ir.response_window_lo, ir.response_integral_peak_index, 
        [&f_of_x,target](double x) { return f_of_x(x)-target; });

    // 3. Calculate FFT of response function
    ir.transform.resize(nsample_);
    plan.execute();

    // 4. Store data in IR array
    ir.units = units;

    impulse_responses_.emplace_back(ir);
    return impulse_responses_.size()-1;
  }

  unsigned copy_and_shift_impulse_response(unsigned impulse_response_id, 
    unsigned ishift, double scale = 1.0)
  {
    validate_impulse_response_id(impulse_response_id);
    const ImpulseResponse& iro = impulse_responses_[impulse_response_id];

    vecX_t shifted_impulse_response(std::min(iro.response_size+ishift, nsample_));
    shifted_impulse_response.setZero();
    shifted_impulse_response.tail(std::min(iro.response_size, nsample_-ishift)) =
      iro.response.head(std::min(iro.response_size, nsample_-ishift)) * scale;

    return register_impulse_response(shifted_impulse_response, iro.units, iro.response_window_frac);
  }
  
  vecX_t impulse_response(unsigned impulse_response_id) const
  {
    validate_impulse_response_id(impulse_response_id);
    const ImpulseResponse& ir = impulse_responses_[impulse_response_id];
    return ir.response;
  }

  vecX_t impulse_response_fft(unsigned impulse_response_id) const
  {
    validate_impulse_response_id(impulse_response_id);
    const ImpulseResponse& ir = impulse_responses_[impulse_response_id];
    return ir.transform;
  }

  std::string impulse_response_summary(unsigned impulse_response_id) const
  {
    validate_impulse_response_id(impulse_response_id);
    const ImpulseResponse& ir = impulse_responses_[impulse_response_id];
    std::string s;
    s += "Dur=" + calin::util::string::to_string_with_commas(ir.response_size*time_resolution_ns_,1) + "ns";
    s += " Peak=" + calin::util::string::to_string_with_commas(ir.response_peak_value,1) + ir.units;
    s += "@" + calin::util::string::to_string_with_commas(ir.response_peak_index*time_resolution_ns_,1) + "ns";
    s += " IntMax=" + calin::util::string::to_string_with_commas(ir.response_integral_peak_value*time_resolution_ns_,1) + ir.units + ".ns";
    s += "@" + calin::util::string::to_string_with_commas(ir.response_integral_peak_index*time_resolution_ns_,1) + "ns";
    s += " IntWin(" + calin::util::string::to_string_with_commas(ir.response_window_frac*100.0,0);
    s += "%)=" + calin::util::string::to_string_with_commas(ir.response_window_lo*time_resolution_ns_,1);
    s += "->" + calin::util::string::to_string_with_commas(ir.response_window_hi*time_resolution_ns_,1) + "ns";
    return s;
  }

  double impulse_response_peak_value(unsigned impulse_response_id) const
  {
    validate_impulse_response_id(impulse_response_id);
    const ImpulseResponse& ir = impulse_responses_[impulse_response_id];
    return ir.response_peak_value;
  }

  double impulse_response_peak_time_ns(unsigned impulse_response_id) const
  {
    validate_impulse_response_id(impulse_response_id);
    const ImpulseResponse& ir = impulse_responses_[impulse_response_id];
    return ir.response_peak_index * time_resolution_ns_;
  }

  double impulse_response_integral_peak_value(unsigned impulse_response_id) const
  {
    validate_impulse_response_id(impulse_response_id);
    const ImpulseResponse& ir = impulse_responses_[impulse_response_id];
    return ir.response_integral_peak_value * time_resolution_ns_;
  }

  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  //
  //    d8888b. d888888b d8888b. d88888b  .o88b. d888888b 
  //    88  `8D   `88'   88  `8D 88'     d8P  Y8 `~~88~~' 
  //    88   88    88    88oobY' 88ooooo 8P         88    
  //    88   88    88    88`8b   88~~~~~ 8b         88    
  //    88  .8D   .88.   88 `88. 88.     Y8b  d8    88    
  //    Y8888D' Y888888P 88   YD Y88888P  `Y88P'    YP    
  //
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  void convolve_impulse_response_direct(unsigned impulse_response_id, unsigned first_sample_of_interest,
    const vecX_t& pedestal = vecX_t(), const vecX_t& relative_gain = vecX_t(), double white_noise_rms = 0.0)
  {
    if(pedestal.size()!=0 and pedestal.size()!=npix_) {
      throw std::domain_error("Pedestal vector length is not equal to number of pixels " 
        + std::to_string(pedestal.size()) + " != " + std::to_string(npix_));
    }

    if(relative_gain.size()!=0 and relative_gain.size()!=npix_) {
      throw std::domain_error("Relative gain vector length is not equal to number of pixels " 
        + std::to_string(relative_gain.size()) + " != " + std::to_string(npix_));
    }

    validate_impulse_response_id(impulse_response_id);
    const ImpulseResponse& ir = impulse_responses_[impulse_response_id];
  
    first_sample_of_interest -= first_sample_of_interest%VCLReal::num_real;

    for(unsigned ipix=0; ipix<npix_; ++ipix) {
      real_t *__restrict__ v_waveform_ptr = &v_waveform_(0, ipix);
      const real_t *__restrict__ pe_waveform_ptr = &pe_waveform_(0, ipix);
      for(unsigned isample=first_sample_of_interest;isample<nsample_;isample+=VCLReal::num_real) {
        real_vt x;
        if(pedestal.size()) {
          x = pedestal[ipix];
        } else {
          x = 0;
        }
        if(white_noise_rms != 0.0) {
          x += white_noise_rms * rng_->template normal_real<VCLReal>();
        }
        x.store_a(v_waveform_ptr + isample);
      }
      
      real_t gain = 1.0;
      if(relative_gain.size()) {
        gain = relative_gain[ipix];
      }
      unsigned isample0 = std::max(first_sample_of_interest, ir.response_size) - ir.response_size;
      for(unsigned isample=isample0; isample<nsample_; ++isample) {
        real_t wt = gain * pe_waveform_ptr[isample];
        if(wt==0) {
          continue;
        }

        unsigned jmin = std::max(first_sample_of_interest, isample) - isample;
        unsigned jmax = std::min(unsigned(ir.response_size), nsample_ - isample);
        unsigned jcount = std::max(jmax, jmin) - jmin;
        real_t *__restrict__ w_waveform_ptr = v_waveform_ptr + isample + jmin;
        const real_t *__restrict__ impulse_response_ptr = ir.response.data() + jmin;
        for(unsigned jsample=0; jsample<jcount; ++jsample) {
          *(w_waveform_ptr++) += wt * (*(impulse_response_ptr++));
        }
      }
    }
    v_waveform_is_packed_ = false;
  }

  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  //
  //    d88888b d88888b d888888b db   d8b   db 
  //    88'     88'     `~~88~~' 88   I8I   88 
  //    88ooo   88ooo      88    88   I8I   88 
  //    88~~~   88~~~      88    Y8   I8I   88 
  //    88      88         88    `8b d8'8b d8' 
  //    YP      YP         YP     `8b8' `8d8'  
  //
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  void convolve_impulse_response_fft(unsigned impulse_response_id,
    const vecX_t& pedestal = vecX_t(),
    const vecX_t& relative_gain = vecX_t(),
    const vecX_t& noise_spectrum = vecX_t())
  {
    if(pedestal.size()!=0 and pedestal.size()!=npix_) {
      throw std::domain_error("Pedestal vector length is not equal to number of pixels " 
        + std::to_string(pedestal.size()) + " != " + std::to_string(npix_));
    }

    if(relative_gain.size()!=0 and relative_gain.size()!=npix_) {
      throw std::domain_error("Relative gain vector length is not equal to number of pixels " 
        + std::to_string(relative_gain.size()) + " != " + std::to_string(npix_));
    }

    if(noise_spectrum.size()!=0 and noise_spectrum.size()!=nsample_) {
      throw std::domain_error("Noise spectrum vector length is not equal to number of samples " 
        + std::to_string(noise_spectrum.size()) + " != " + std::to_string(nsample_));
    }

    validate_impulse_response_id(impulse_response_id);
    const ImpulseResponse& ir = impulse_responses_[impulse_response_id];

    real_t* a_vec = calin::util::memory::aligned_calloc<real_t>(nsample_);
    real_t* b_vec = calin::util::memory::aligned_calloc<real_t>(nsample_);

    auto fwd = calin::math::fftw_util::plan_r2hc(nsample_, a_vec, b_vec, FFTW_MEASURE);
    auto bwd = calin::math::fftw_util::plan_hc2r(nsample_, a_vec, b_vec, FFTW_MEASURE);
    const real_t scale = 1.0/nsample_;

    for(unsigned ipix=0; ipix<npix_; ++ipix) {
      std::copy(&pe_waveform_(0,ipix), &pe_waveform_(0,ipix)+nsample_, a_vec);

      fwd.execute();

      real_t gain = 1.0;
      if(relative_gain.size()) {
        gain = relative_gain[ipix];
      }

      calin::math::fftw_util::hcvec_scale_and_multiply(a_vec, b_vec, ir.transform.data(),
        nsample_, scale*gain);

      if(pedestal.size() > ipix) {
        a_vec[0] += pedestal(ipix);
      }

      if(noise_spectrum.size()) {
        for(unsigned isample=0;isample<nsample_;isample+=VCLReal::num_real) {
          real_vt x;
          x.load(&noise_spectrum[isample]);
          if(to_bits(x != 0.0)) {
            x *= rng_->template normal_real<VCLReal>();
            real_vt a;
            a.load_a(&a_vec[isample]);
            a += x;
            a.store_a(&a_vec[isample]);
          }
        }
      }

      bwd.execute();
      std::copy(b_vec, b_vec+nsample_, &v_waveform_(0,ipix));
    }
    v_waveform_is_packed_ = false;
    free(a_vec);
    free(b_vec);
  }

  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  //
  //     .o88b.  .d88b.  d8888b. d88888b db      d88888b d888888b 
  //    d8P  Y8 .8P  Y8. 88  `8D 88'     88      88'     `~~88~~' 
  //    8P      88    88 88   88 88ooooo 88      88ooooo    88    
  //    8b      88    88 88   88 88~~~~~ 88      88~~~~~    88    
  //    Y8b  d8 `8b  d8' 88  .8D 88.     88booo. 88.        88    
  //     `Y88P'  `Y88P'  Y8888D' Y88888P Y88888P Y88888P    YP    
  //
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  void convolve_impulse_response_fftw_codelet(unsigned impulse_response_id,
    const vecX_t& pedestal = vecX_t(),
    const vecX_t& relative_gain = vecX_t(),
    const vecX_t& noise_spectrum = vecX_t())
  {
    calin::math::fftw_util::FFTWCodelet<VCLReal> fft;
    if(!fft.has_codelet(nsample_)) {
      throw std::domain_error("FFTW codelet for size " + std::to_string(nsample_) 
        + " not available.");
    }

    if(pedestal.size()!=0 and pedestal.size()!=npix_) {
      throw std::domain_error("Pedestal vector length is not equal to number of pixels " 
        + std::to_string(pedestal.size()) + " != " + std::to_string(npix_));
    }

    if(relative_gain.size()!=0 and relative_gain.size()!=npix_) {
      throw std::domain_error("Relative gain vector length is not equal to number of pixels " 
        + std::to_string(relative_gain.size()) + " != " + std::to_string(npix_));
    }

    if(noise_spectrum.size()!=0 and noise_spectrum.size()!=nsample_) {
      throw std::domain_error("Noise spectrum vector length is not equal to number of samples " 
        + std::to_string(noise_spectrum.size()) + " != " + std::to_string(nsample_));
    }

    validate_impulse_response_id(impulse_response_id);
    const ImpulseResponse& ir = impulse_responses_[impulse_response_id];

    real_vt* a_vec = calin::util::memory::aligned_calloc<real_vt>(nsample_);
    real_vt* b_vec = calin::util::memory::aligned_calloc<real_vt>(nsample_);

    real_t scale = 1.0/nsample_;

    for(unsigned ipix=0; ipix<npix_; ipix+=VCLReal::num_real) {
      // Load data from "pe_waveform_" into "a_vec", block by block, transposing as we go along
      for(unsigned isample=0; isample<nsample_; isample += VCLReal::num_real) {
        real_vt block[VCLReal::num_real]; // square matrix of doubles
        for(unsigned jpix=0, mpix=std::min(npix_-ipix,VCLReal::num_real); jpix<mpix; jpix++) {
          block[jpix].load_a(pe_waveform_.data() + (ipix+jpix)*nsample_ + isample);
        }
        calin::util::vcl::transpose(block);
        for(unsigned jsample = 0; jsample<VCLReal::num_real; jsample++) {
          a_vec[isample+jsample] = block[jsample];
        }
      }

      // Do FFT of "a_vec" into "b_vec"
      fft.r2hc(nsample_, a_vec, b_vec);

      // Extract relative gain from array
      real_vt gain = 1.0;
      if(relative_gain.size()) {
        if(ipix + VCLReal::num_real <= npix_) {
          gain.load(relative_gain.data() + ipix);
        } else {
          real_at gain_array;
          for(unsigned i=0; i+ipix<npix_; i++) {
            gain_array[i] = relative_gain(ipix+i);
          }
          gain.load(gain_array);
        }
      }

      // Multiply FFT in "b_vec" by impulse response and gain, tranferring back into "a_vec"
      calin::math::fftw_util::hcvec_scale_and_multiply_non_overlapping(a_vec, b_vec, ir.transform.data(), 
        nsample_, gain*scale);

      // Add pedestal in frequency domain as DC offset before inverse transform
      if(pedestal.size()) {
        real_vt ped;
        if(ipix + VCLReal::num_real <= npix_) {
          ped.load(pedestal.data() + ipix);
        } else {
          real_at ped_array;
          for(unsigned i=0; i+ipix<npix_; i++) {
            ped_array[i] = pedestal(ipix+i);
          }
          ped.load(ped_array);
        }
        *a_vec += ped;
      }

      // Add Gaussian noise with given spectrum
      if(noise_spectrum.size()) {
        const real_t *__restrict__ noise_spectrum_ptr = noise_spectrum.data();
        int isample = 0;
        for(int zsample = int(nsample_)-3; isample<zsample; isample+=4) {
          // Unroll here and use normal pair function which is (might be) faster
          real_vt x0, x1, x2, x3;
          rng_->normal_pair_real(x0, x1);
          rng_->normal_pair_real(x2, x3);
          a_vec[isample]   += x0 * noise_spectrum_ptr[isample];
          a_vec[isample+1] += x1 * noise_spectrum_ptr[isample+1];
          a_vec[isample+2] += x2 * noise_spectrum_ptr[isample+2];
          a_vec[isample+3] += x3 * noise_spectrum_ptr[isample+3];
        }
        for(int zsample = int(nsample_); isample<zsample; isample++) {
          real_vt x = rng_->template normal_real<VCLReal>();
          a_vec[isample] += x * noise_spectrum_ptr[isample];
        }        
      }

      // Do inverse-FFT of "a_vec" into "b_vec"
      fft.hc2r(nsample_, a_vec, b_vec);

      for(unsigned isample=0;isample<nsample_;isample++) {
        b_vec[isample].store_a(v_waveform_.data() + ipix*nsample_ + isample*VCLReal::num_real);
      }
    }

    v_waveform_is_packed_ = true;

    free(a_vec);
    free(b_vec);
  }

  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  //
  //    .88b  d88. db    db db      d888888b d888888b d8888b. db      d88888b 
  //    88'YbdP`88 88    88 88      `~~88~~'   `88'   88  `8D 88      88'     
  //    88  88  88 88    88 88         88       88    88oodD' 88      88ooooo 
  //    88  88  88 88    88 88         88       88    88~~~   88      88~~~~~ 
  //    88  88  88 88b  d88 88booo.    88      .88.   88      88booo. 88.     
  //    YP  YP  YP ~Y8888P' Y88888P    YP    Y888888P 88      Y88888P Y88888P 
  //
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  void convolve_multiple_impulse_response_fftw_codelet(
    const Eigen::VectorXi& impulse_response_id,
    const vecX_t& pedestal = vecX_t(),
    const vecX_t& relative_gain = vecX_t(),
    const matX_t& noise_spectrum = matX_t())
  {
    calin::math::fftw_util::FFTWCodelet<VCLReal> fft;
    if(!fft.has_codelet(nsample_)) {
      throw std::domain_error("FFTW codelet for size " + std::to_string(nsample_) 
        + " not available.");
    }

    if(impulse_response_id.size() != npix_) {
      throw std::domain_error("Impulse response id vector length is not equal to number of pixels " 
        + std::to_string(impulse_response_id.size()) + " != " + std::to_string(npix_));
    }

    if(pedestal.size()!=0 and pedestal.size()!=npix_) {
      throw std::domain_error("Pedestal vector length is not equal to number of pixels " 
        + std::to_string(pedestal.size()) + " != " + std::to_string(npix_));
    }
    
    if(relative_gain.size()!=0 and relative_gain.size()!=npix_) {
      throw std::domain_error("Relative gain vector length is not equal to number of pixels " 
        + std::to_string(relative_gain.size()) + " != " + std::to_string(npix_));
    }

    if(noise_spectrum.size()!=0 and (noise_spectrum.rows()!=nsample_
        or (noise_spectrum.cols()!=1 and noise_spectrum.cols()!=npix_))) {
      throw std::domain_error("Noise spectrum matrix has incoorect size: [" 
        + std::to_string(noise_spectrum.rows()) + "," + std::to_string(noise_spectrum.cols()) + "] != ["
        + std::to_string(nsample_) + ",1] or ["
        + std::to_string(nsample_) + "," + std::to_string(npix_) + "]");
    }

    int iri_max = std::numeric_limits<int>::min();
    int iri_min = std::numeric_limits<int>::max();
    for(unsigned ipix=0; ipix<npix_; ipix++) {
      int iri = impulse_response_id(ipix);
      iri_max = std::max(iri_max, iri);
      iri_min = std::min(iri_min, iri);
    }

    validate_impulse_response_id(iri_max);
    if(iri_min<0) {
      throw std::out_of_range("Invalid impulse response id cannot be negative: " 
        + std::to_string(iri_min));
    }

    real_vt* a_vec = calin::util::memory::aligned_calloc<real_vt>(nsample_);
    real_vt* b_vec = calin::util::memory::aligned_calloc<real_vt>(nsample_);
    real_vt* c_vec = calin::util::memory::aligned_calloc<real_vt>(nsample_);
  
    real_t scale = 1.0/nsample_;

    for(unsigned ipix=0; ipix<npix_; ipix+=VCLReal::num_real) {
      // Load data from "pe_waveform_" into "a_vec", block by block, transposing as we go along
      for(unsigned isample=0; isample<nsample_; isample += VCLReal::num_real) {
        real_vt block[VCLReal::num_real]; // square matrix of doubles
        for(unsigned jpix=0, mpix=std::min(npix_-ipix,VCLReal::num_real); jpix<mpix; jpix++) {
          block[jpix].load_a(pe_waveform_.data() + (ipix+jpix)*nsample_ + isample);
        }
        calin::util::vcl::transpose(block);
        for(unsigned jsample = 0; jsample<VCLReal::num_real; jsample++) {
          a_vec[isample+jsample] = block[jsample];
        }
      }
    
      // Do FFT of "a_vec" into "b_vec"
      fft.r2hc(nsample_, a_vec, b_vec);

      // Load data from "impulse_respose" into "c_vec", block by block, transposing as we go along
      for(unsigned isample=0; isample<nsample_; isample += VCLReal::num_real) {
        real_vt block[VCLReal::num_real]; // square matrix of doubles
        for(unsigned jpix=0, mpix=std::min(npix_-ipix,VCLReal::num_real); jpix<mpix; jpix++) {
          int irid = impulse_response_id[ipix + jpix];
          block[jpix].load_a(impulse_responses_[irid].transform.data() + isample);
        }
        calin::util::vcl::transpose(block);
        for(unsigned jsample = 0; jsample<VCLReal::num_real; jsample++) {
          c_vec[isample+jsample] = block[jsample];
        }
      }

      // Extract relative gain from array
      real_vt gain = 1.0;
      if(relative_gain.size()) {
        if(ipix + VCLReal::num_real <= npix_) {
          gain.load(relative_gain.data() + ipix);
        } else {
          real_at gain_array;
          for(unsigned i=0; i+ipix<npix_; i++) {
            gain_array[i] = relative_gain(ipix+i);
          }
          gain.load(gain_array);
        }
      }

      // Multiply FFT in "b_vec" by impulse response in "c_vec" and gain, then transfer back into "a_vec"
      calin::math::fftw_util::hcvec_scale_and_multiply_non_overlapping(a_vec, b_vec, c_vec, 
        nsample_, gain*scale);

      // Add pedestal in frequency domain as DC offset before inverse transform
      if(pedestal.size()) {
        real_vt ped;
        if(ipix + VCLReal::num_real <= npix_) {
          ped.load(pedestal.data() + ipix);
        } else {
          real_at ped_array;
          for(unsigned i=0; i+ipix<npix_; i++) {
            ped_array[i] = pedestal(ipix+i);
          }
          ped.load(ped_array);
        }
        *a_vec += ped;
      }

      // Add Gaussian noise with given spectrum
      if(noise_spectrum.cols() == 1) {
        // Unique noise spectrum for all channels
        int isample = 0;
        const real_t *__restrict__ noise_spectrum_ptr = noise_spectrum.data();
        for(int zsample = int(nsample_)-3; isample<zsample; isample+=4) {
          // Unroll by two here and use Gaussian pair function to produce 2 deviates
          real_vt x0, x1, x2, x3;
          rng_->normal_pair_real(x0, x1);
          rng_->normal_pair_real(x2, x3);
          a_vec[isample]   += x0 * noise_spectrum_ptr[isample];
          a_vec[isample+1] += x1 * noise_spectrum_ptr[isample+1];
          a_vec[isample+2] += x2 * noise_spectrum_ptr[isample+2];
          a_vec[isample+3] += x3 * noise_spectrum_ptr[isample+3];
        }
        for(int zsample = int(nsample_); isample<zsample; isample++) {
          real_vt x = rng_->template normal_real<VCLReal>();
          a_vec[isample] += x * noise_spectrum_ptr[isample];
        }        
      } else if (noise_spectrum.size()) {
        // Load data from "noise_spectrum", block by block, transposing as we go along
        const real_t *__restrict__ noise_spectrum_ptr = noise_spectrum.data();
        for(unsigned isample=0; isample<nsample_; isample += VCLReal::num_real) {
          real_vt block[VCLReal::num_real]; // square matrix of doubles
          for(unsigned jpix=0, mpix=std::min(npix_-ipix,VCLReal::num_real); jpix<mpix; jpix++) {
            block[jpix].load_a(noise_spectrum_ptr + (ipix+jpix)*nsample_ + isample);
          }
          calin::util::vcl::transpose(block);
          // One of these three possible cases which is chosen at compile time
          if(VCLReal::num_real % 4 == 0) {        
            // If there is multiple of 4 reals in the vector then do an unrolled 
            // loop generating two pairs of Gaussian deviates each iteration.
            // This branch is chose for all AVX/AVX-2 and AVX-512 float and double 
            // architectures, and in the SSE (128 bit) float architecture.
            for(unsigned jsample = 0; jsample<VCLReal::num_real; jsample+=4) {
              real_vt x0, x1, x2, x3;
              rng_->normal_pair_real(x0, x1);
              rng_->normal_pair_real(x2, x3);
              a_vec[isample]   += x0 * block[jsample];
              a_vec[isample+1] += x1 * block[jsample+1];
              a_vec[isample+2] += x2 * block[jsample+2];
              a_vec[isample+3] += x3 * block[jsample+3];
            }
          } else if (VCLReal::num_real % 2 == 0) {
            // If there is multiple of 2 reals in the vector then do loop 
            // generating one pairs of Gaussian deviates each iteration. This 
            // branch is chosen is for the 128-bit SSE double architecture.
            for(unsigned jsample = 0; jsample<VCLReal::num_real; jsample+=2) {
              real_vt x0, x1;
              rng_->normal_pair_real(x0, x1);
              a_vec[isample]   += x0 * block[jsample];
              a_vec[isample+1] += x1 * block[jsample+1];
            }
          } else {
            // Otherwise just generate one Gaussian per iteration. There is
            // no case in VCL where this branch is chosen.
            for(unsigned jsample = 0; jsample<VCLReal::num_real; jsample++) {
              real_vt x = rng_->template normal_real<VCLReal>();
              a_vec[isample] += x * block[jsample];
            }
          }
        }
      }

      // Do inverse-FFT of "a_vec" into "b_vec"
      fft.hc2r(nsample_, a_vec, b_vec);

      for(unsigned isample=0;isample<nsample_;isample++) {
        b_vec[isample].store_a(v_waveform_.data() + ipix*nsample_ + isample*VCLReal::num_real);
      }
    }
    v_waveform_is_packed_ = true;
    free(a_vec);
    free(b_vec);
    free(c_vec);
  }

  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  //
  //    .d8888. d888888b d8b   db  d888b  db      d88888b .d8888. 
  //    88'  YP   `88'   888o  88 88' Y8b 88      88'     88'  YP 
  //    `8bo.      88    88V8o 88 88      88      88ooooo `8bo.   
  //      `Y8b.    88    88 V8o88 88  ooo 88      88~~~~~   `Y8b. 
  //    db   8D   .88.   88  V888 88. ~8~ 88booo. 88.     db   8D 
  //    `8888Y' Y888888P VP   V8P  Y888P  Y88888P Y88888P `8888Y' 
  //                                                              
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  unsigned trigger_singles(Eigen::VectorXi& pixel_trigger_time, const vecX_t& threshold,
    unsigned first_sample_of_interest=0)
  {
    if(threshold.size()!=1 and threshold.size()!=npix_) {
      throw std::domain_error("Trigger threshold vector length is not equal to one or number of pixels " 
        + std::to_string(threshold.size()) + " != 1 or " + std::to_string(npix_));
    }
    pixel_trigger_time.resize(npix_);
    int *__restrict__ ptt_ptr = pixel_trigger_time.data();

    unsigned ntriggered = 0;
    for(unsigned ipix=0; ipix<npix_; ipix+=VCLReal::num_real) { 
      real_vt pix_threshold = 0;
      if(threshold.size() == 1) {
        pix_threshold = threshold[0];
      } else if(ipix + VCLReal::num_real <= npix_) {
        pix_threshold.load_a(threshold.data() + ipix);
      } else {
        real_at threshold_array;
        pix_threshold.store_a(threshold_array);
        for(unsigned i=0; i+ipix<npix_; i++) {
          threshold_array[i] = threshold(ipix+i);
        }
        pix_threshold.load_a(threshold_array);
      }
      pix_threshold = select(VCLReal::int_iota()+ipix < npix_, pix_threshold, std::numeric_limits<real_t>::infinity());        

      int_vt ptt = -1;
      if(v_waveform_is_packed_) {
        unsigned isample0 = first_sample_of_interest;
        const real_t *__restrict__ v_waveform_ptr = v_waveform_.data() + ipix*nsample_ + isample0*VCLReal::num_real;
        for(unsigned isample = isample0; isample<nsample_; ++isample) {
          real_vt v;
          v.load_a(v_waveform_ptr);
          v_waveform_ptr += VCLReal::num_real;
          ptt = vcl::select(int_bvt(v > pix_threshold) && (ptt < 0), int_vt(isample), ptt);
        }
      } else {
        unsigned isample0 = first_sample_of_interest;
        isample0 = (isample0/VCLReal::num_real)*VCLReal::num_real;
        for(unsigned isample = isample0; isample<nsample_; isample+=VCLReal::num_real) {
          const real_t *__restrict__ v_waveform_ptr = v_waveform_.data() + ipix*nsample_ + isample;
          real_vt block[VCLReal::num_real]; // square matrix of reals
          for(unsigned jpix=0;jpix<VCLReal::num_real;jpix++) {
            block[jpix].load_a(v_waveform_ptr);
            v_waveform_ptr += nsample_;
          }
          calin::util::vcl::transpose(block);
          for(unsigned jsample=std::max(first_sample_of_interest, isample)-isample;jsample<VCLReal::num_real;++jsample) {
            unsigned ksample = isample+jsample;
            ptt = vcl::select(int_bvt(block[jsample] > pix_threshold) && (ptt < 0), int_vt(ksample), ptt);
          }
        }
      }
      int_at ptt_array;
      ptt.store_a(ptt_array);
      for(unsigned jpix=0; jpix+ipix<npix_; jpix++) {
        ptt_ptr[ipix+jpix] = ptt_array[jpix];
      }
      ntriggered += vcl::horizontal_count(ptt >= 0);
    }
    return ntriggered;
  }

  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  //
  //    .88b  d88. db      d888888b d8888b. db       .o88b. d888888b db    db 
  //    88'YbdP`88 88      `~~88~~' 88  `8D 88      d8P  Y8 `~~88~~' `8b  d8' 
  //    88  88  88 88         88    88oodD' 88      8P         88     `8bd8'  
  //    88  88  88 88         88    88~~~   88      8b         88       88    
  //    88  88  88 88booo.    88    88      88booo. Y8b  d8    88       88    
  //    YP  YP  YP Y88888P    YP    88      Y88888P  `Y88P'    YP       YP    
  //
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  int trigger_multiplicity(const vecX_t& threshold, unsigned multiplicity, unsigned coincidence_window, unsigned first_sample_of_interest=0)
  {
    if(threshold.size()!=1 and threshold.size()!=npix_) {
      throw std::domain_error("Trigger threshold vector length is not equal to one or number of pixels " 
        + std::to_string(threshold.size()) + " != 1 or " + std::to_string(npix_));
    }

    coincidence_window = std::max(coincidence_window, 1U); // Zero makes no sense
    
    int_t *__restrict__ cwin_tend = calin::util::memory::aligned_calloc<int_t>(v_waveform_.cols());
    std::fill(cwin_tend, cwin_tend+v_waveform_.cols(), -1);

    unsigned isample0 = std::max(first_sample_of_interest,coincidence_window)-coincidence_window;
    isample0 = (isample0/VCLReal::num_real)*VCLReal::num_real;
    for(unsigned isample = isample0; isample<nsample_; isample+=VCLReal::num_real) {
      unsigned ntriggered[VCLReal::num_real];
      std::fill(ntriggered, ntriggered+VCLReal::num_real, 0);
      real_vt block[VCLReal::num_real]; // square matrix of reals
      for(unsigned ipix=0; ipix<npix_; ipix+=VCLReal::num_real) {
        real_vt pix_threshold = 0;
        if(threshold.size() == 1) {
          pix_threshold = threshold[0];
        } else if(ipix + VCLReal::num_real <= npix_) {
          pix_threshold.load_a(threshold.data() + ipix);
        } else {
          real_at threshold_array;
          pix_threshold.store_a(threshold_array);
          for(unsigned i=0; i+ipix<npix_; i++) {
            threshold_array[i] = threshold(ipix+i);
          }
          pix_threshold.load_a(threshold_array);
        }
        pix_threshold = select(VCLReal::int_iota()+ipix < npix_, pix_threshold, std::numeric_limits<real_t>::infinity());

        if(v_waveform_is_packed_) {
          const real_t *__restrict__ v_waveform_ptr = v_waveform_.data() + ipix*nsample_ + isample*VCLReal::num_real;
          for(unsigned jsample=0;jsample<VCLReal::num_real;jsample++) {
            block[jsample].load_a(v_waveform_ptr);
            v_waveform_ptr += VCLReal::num_real;
          }
        } else {
          const real_t *__restrict__ v_waveform_ptr = v_waveform_.data() + ipix*nsample_ + isample;
          for(unsigned jpix=0;jpix<VCLReal::num_real;jpix++) {
            block[jpix].load_a(v_waveform_ptr);
            v_waveform_ptr += nsample_;
          }
          calin::util::vcl::transpose(block);
        }

        int_vt cwin;
        cwin.load_a(cwin_tend + ipix);
        for(unsigned jsample=0;jsample<VCLReal::num_real;jsample++) {
          int ksample = isample+jsample;
          cwin = vcl::select(int_bvt(block[jsample] > pix_threshold), ksample+coincidence_window, cwin);
          ntriggered[jsample] += vcl::horizontal_count(ksample < cwin);
        }
        cwin.store_a(cwin_tend + ipix);
      }

      for(unsigned jsample=0;jsample<VCLReal::num_real;jsample++) {
        unsigned ksample = isample+jsample;
        if(ksample>=first_sample_of_interest and ntriggered[jsample] >= multiplicity) {
          ::free(cwin_tend);
          return ksample;
        }
      }
    }
    ::free(cwin_tend);
    return -1;
  }

  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  //
  //    d8888b.  .d8b.  d8b   db  .d88b.  .d8888. d88888b d888888b d888888b 
  //    88  `8D d8' `8b 888o  88 .8P  Y8. 88'  YP 88'     `~~88~~'   `88'   
  //    88oodD' 88ooo88 88V8o 88 88    88 `8bo.   88ooooo    88       88    
  //    88~~~   88~~~88 88 V8o88 88    88   `Y8b. 88~~~~~    88       88    
  //    88      88   88 88  V888 `8b  d8' db   8D 88.        88      .88.   
  //    88      YP   YP VP   V8P  `Y88P'  `8888Y' Y88888P    YP    Y888888P 
  //
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  int trigger_panoseti_multiplicity(const vecX_t& threshold, unsigned multiplicity,
      unsigned coincidence_window, unsigned first_sample_of_interest=0)
  {
    // Multiplicity trigger algorithm with triggering of the four 16x16
    // channel quabos handled independently. Camera-level trigger is
    // satisfied when any quabo triggers.

    // Algorithm assumes all values in SIMD vector belong to same Quabo.
    // This is always true in current Intel SIMD vectors
    static_assert(16 % VCLReal::num_real == 0,
      "SIMD lane count must divide the 16-channel PANOSETI module width");

    if(npix_ != 1024) {
      throw std::domain_error("PANOSETI multiplicity trigger requires exactly 1024 channels; got "
        + std::to_string(npix_));
    }
    if(threshold.size()!=1 and threshold.size()!=npix_) {
      throw std::domain_error("Trigger threshold vector length is not equal to one or number of pixels "
        + std::to_string(threshold.size()) + " != 1 or " + std::to_string(npix_));
    }

    coincidence_window = std::max(coincidence_window, 1U); // Zero makes no sense

    int_t *__restrict__ cwin_tend = calin::util::memory::aligned_calloc<int_t>(v_waveform_.cols());
    std::fill(cwin_tend, cwin_tend+v_waveform_.cols(), -1);

    unsigned isample0 = std::max(first_sample_of_interest,coincidence_window)-coincidence_window;
    isample0 = (isample0/VCLReal::num_real)*VCLReal::num_real;
    for(unsigned isample = isample0; isample<nsample_; isample+=VCLReal::num_real) {
      unsigned ntriggered[4][VCLReal::num_real];
      for(unsigned iquabo=0; iquabo<4; iquabo++) {
        std::fill(ntriggered[iquabo], ntriggered[iquabo]+VCLReal::num_real, 0);
      }
      real_vt block[VCLReal::num_real];
      for(unsigned ipix=0; ipix<npix_; ipix+=VCLReal::num_real) {
        real_vt pix_threshold = 0;
        if(threshold.size() == 1) {
          pix_threshold = threshold[0];
        } else {
          pix_threshold.load_a(threshold.data() + ipix);
        }

        if(v_waveform_is_packed_) {
          const real_t *__restrict__ v_waveform_ptr = v_waveform_.data() + ipix*nsample_ + isample*VCLReal::num_real;
          for(unsigned jsample=0;jsample<VCLReal::num_real;jsample++) {
            block[jsample].load_a(v_waveform_ptr);
            v_waveform_ptr += VCLReal::num_real;
          }
        } else {
          const real_t *__restrict__ v_waveform_ptr = v_waveform_.data() + ipix*nsample_ + isample;
          for(unsigned jpix=0;jpix<VCLReal::num_real;jpix++) {
            block[jpix].load_a(v_waveform_ptr);
            v_waveform_ptr += nsample_;
          }
          calin::util::vcl::transpose(block);
        }

        const unsigned row = ipix / 32;
        const unsigned col = ipix % 32;
        const unsigned iquabo = (row / 16) * 2 + (col / 16);

        int_vt cwin;
        cwin.load_a(cwin_tend + ipix);
        for(unsigned jsample=0;jsample<VCLReal::num_real;jsample++) {
          int ksample = isample+jsample;
          cwin = vcl::select(int_bvt(block[jsample] > pix_threshold),
            ksample+coincidence_window, cwin);
          ntriggered[iquabo][jsample] += vcl::horizontal_count(ksample < cwin);
        }
        cwin.store_a(cwin_tend + ipix);
      }

      for(unsigned jsample=0;jsample<VCLReal::num_real;jsample++) {
        unsigned ksample = isample+jsample;
        if(ksample>=first_sample_of_interest) {
          for(unsigned iquabo=0; iquabo<4; iquabo++) {
            if(ntriggered[iquabo][jsample] >= multiplicity) {
              ::free(cwin_tend);
              return ksample;
            }
          }
        }
      }
    }
    ::free(cwin_tend);
    return -1;
  }

  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  //
  //    d8888b. d8b   db d8b   db 
  //    VP  `8D 888o  88 888o  88 
  //      oooY' 88V8o 88 88V8o 88 
  //      ~~~b. 88 V8o88 88 V8o88 
  //    db   8D 88  V888 88  V888 
  //    Y8888P' VP   V8P VP   V8P 
  //
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  int trigger_3nn_old(const vecX_t& threshold, const Eigen::MatrixXi& neighbors, unsigned coincidence_window, unsigned first_sample_of_interest=0)
  {
    using mask_t = uint32_t;
    constexpr unsigned mask_t_size_bits = sizeof(mask_t)*8;
    static_assert(VCLReal::num_real <= mask_t_size_bits);
    static_assert(mask_t_size_bits % VCLReal::num_real == 0);

    if(threshold.size()!=1 and threshold.size()!=npix_) {
      throw std::domain_error("Trigger threshold vector length is not equal to one or number of pixels " 
        + std::to_string(threshold.size()) + " != 1 or " + std::to_string(npix_));
    }

    if(neighbors.cols()!=npix_) {
      throw std::domain_error("Number of columns in neighbors matrix must equal number of pixels " 
        + std::to_string(neighbors.cols()) + " != " + std::to_string(npix_));
    }

    coincidence_window = std::max(coincidence_window, 1U); // Zero makes no sense
    
    int_t *__restrict__ cwin_tend = calin::util::memory::aligned_calloc<int_t>(v_waveform_.cols());
    std::fill(cwin_tend, cwin_tend+v_waveform_.cols(), -1);

    unsigned ntriggered[VCLReal::num_real];
    bool new_triggers[VCLReal::num_real];

    unsigned trigger_hit_size = (v_waveform_.cols() + mask_t_size_bits - 1)/mask_t_size_bits;
    mask_t *__restrict__ trigger_hit[VCLReal::num_real];
    for(unsigned i=0; i<VCLReal::num_real; i++) {
      trigger_hit[i] = calin::util::memory::aligned_calloc<mask_t>(trigger_hit_size);
    }

    unsigned isample0 = std::max(first_sample_of_interest,coincidence_window)-coincidence_window;
    isample0 = (isample0/VCLReal::num_real)*VCLReal::num_real;
    for(unsigned isample = isample0; isample<nsample_; isample+=VCLReal::num_real) {
      std::fill(ntriggered, ntriggered+VCLReal::num_real, 0);
      std::fill(new_triggers, new_triggers+VCLReal::num_real, 0);
      for(unsigned i=0; i<VCLReal::num_real; i++) {
        std::fill(trigger_hit[i], trigger_hit[i]+trigger_hit_size, 0);
      }

      real_vt block[VCLReal::num_real]; // square matrix of reals
      for(unsigned ipix=0; ipix<npix_; ipix+=VCLReal::num_real) {
        const unsigned imask = ipix/mask_t_size_bits;
        const unsigned ishift = ((ipix - imask*mask_t_size_bits)/VCLReal::num_real)*VCLReal::num_real;
        real_vt pix_threshold = 0;
        if(threshold.size() == 1) {
          pix_threshold = threshold[0];
        } else if(ipix + VCLReal::num_real <= npix_) {
          pix_threshold.load_a(threshold.data() + ipix);
        } else {
          real_at threshold_array;
          pix_threshold.store_a(threshold_array);
          for(unsigned i=0; i+ipix<npix_; i++) {
            threshold_array[i] = threshold(ipix+i);
          }
          pix_threshold.load_a(threshold_array);
        }
        pix_threshold = select(VCLReal::int_iota()+ipix < npix_, pix_threshold, std::numeric_limits<real_t>::infinity());

        if(v_waveform_is_packed_) {
          const real_t *__restrict__ v_waveform_ptr = v_waveform_.data() + ipix*nsample_ + isample*VCLReal::num_real;
          for(unsigned jsample=0;jsample<VCLReal::num_real;jsample++) {
            block[jsample].load_a(v_waveform_ptr);
            v_waveform_ptr += VCLReal::num_real;
          }
        } else {
          const real_t *__restrict__ v_waveform_ptr = v_waveform_.data() + ipix*nsample_ + isample;
          for(unsigned jpix=0;jpix<VCLReal::num_real;jpix++) {
            block[jpix].load_a(v_waveform_ptr);
            v_waveform_ptr += nsample_;
          }
          calin::util::vcl::transpose(block);
        }

        int_vt cwin;
        cwin.load_a(cwin_tend + ipix);
        for(unsigned jsample=0;jsample<VCLReal::num_real;jsample++) {
          int ksample = isample+jsample;
          int_vt new_cwin = vcl::select(int_bvt(block[jsample] > pix_threshold), ksample+coincidence_window, cwin);
          int_bvt triggered_pix = ksample < new_cwin;
          mask_t triggered_pix_mask = vcl::to_bits(triggered_pix);
          new_triggers[jsample] |= vcl::horizontal_or(triggered_pix && (ksample-1 >= cwin));
          ntriggered[jsample] += std::popcount(triggered_pix_mask);
          trigger_hit[jsample][imask] |= triggered_pix_mask << ishift;
          cwin = new_cwin;
        }
        cwin.store_a(cwin_tend + ipix);
      }

      for(unsigned jsample=0;jsample<VCLReal::num_real;jsample++) {
        unsigned ksample = isample+jsample;
        if(ksample>=first_sample_of_interest and ntriggered[jsample] >= 3 and new_triggers[jsample]) {          
          for(unsigned ith=0; ith<trigger_hit_size; ith++) {
            mask_t th = trigger_hit[jsample][ith];
            const unsigned ipix0 = ith*mask_t_size_bits;
            while(th!=0) {
              const unsigned ihit = std::countr_zero(th);
              const unsigned ipix = ipix0 + ihit;
              th &= th - 1;
              unsigned nneighbor_hit = 0;
              // unsigned ipix_n1 = -1;
              for(auto jpix : neighbors.col(ipix)) {
                if(jpix>=0 && jpix<int(npix_)) {
                  const unsigned jmask = jpix/mask_t_size_bits;
                  const unsigned jshift = jpix - jmask*mask_t_size_bits;
                  if(trigger_hit[jsample][jmask] & (mask_t(1) << jshift)) {
                    ++nneighbor_hit;
                    if(nneighbor_hit >= 2) {
                      // Found 3 neighbor channels that triggered
                      // calin::util::log::LOG(calin::util::log::INFO) << ipix << ' ' << ipix_n1 << ' ' <<  jpix << " ( " << neighbors.col(ipix).transpose() << ')';
                      ::free(cwin_tend);
                      for(unsigned i=0; i<VCLReal::num_real; i++) {
                        ::free(trigger_hit[i]);
                      }
                      return isample+jsample;
                    }
                    // ipix_n1 = jpix;
                  }
                }
              }
            }
          }
        }
      }
    }
    ::free(cwin_tend);
    for(unsigned i=0; i<VCLReal::num_real; i++) {
      ::free(trigger_hit[i]);
    }
    return -1;
  }

  int trigger_3nn(const vecX_t& threshold, const Eigen::MatrixXi& neighbors, unsigned coincidence_window, unsigned first_sample_of_interest=0)
  {
    // Vectorized, unrolled, 3NN trigger algorithm that supports coincidence window
    // 1) For each sample: find which pixels are above threshold and set end of their coincidence window
    // 2) Record whether any pixel is newly triggered (was outside window in last sample, now inside)
    // 3) Count number of pixels triggered (i.e. within their coincidence windows)
    // 4) If there are more than 3 triggered pixels, and at least one is newly triggered then run 3NN algorithm:
    // 5) Run through triggered bitmask finding triggered channels, then
    // 6) Test whether neighbors are also triggered, counting how many
    // 7) If 2 neighbors are found then return trigger time

    using mask_t = typename VCLReal::uint_t;
    constexpr unsigned mask_t_size_bits = sizeof(mask_t)*8;
    static_assert(VCLReal::num_real <= mask_t_size_bits);
    static_assert(mask_t_size_bits % VCLReal::num_real == 0);

    if(threshold.size()!=1 and threshold.size()!=npix_) {
      throw std::domain_error("Trigger threshold vector length is not equal to one or number of pixels " 
        + std::to_string(threshold.size()) + " != 1 or " + std::to_string(npix_));
    }

    if(neighbors.cols()!=npix_) {
      throw std::domain_error("Number of columns in neighbors matrix must equal number of pixels " 
        + std::to_string(neighbors.cols()) + " != " + std::to_string(npix_));
    }

    coincidence_window = std::max(coincidence_window, 1U); // Zero makes no sense
    
    int_t *__restrict__ cwin_tend = calin::util::memory::aligned_calloc<int_t>(v_waveform_.cols());
    std::fill(cwin_tend, cwin_tend+v_waveform_.cols(), -1);

    unsigned ntriggered[VCLReal::num_real];
    bool new_triggers[VCLReal::num_real];

    unsigned trigger_hit_size = (v_waveform_.cols() + mask_t_size_bits - 1)/mask_t_size_bits;
    mask_t *__restrict__ trigger_hit_array = 
      calin::util::memory::aligned_calloc<mask_t>(trigger_hit_size * VCLReal::num_real);

    unsigned isample0 = std::max(first_sample_of_interest,coincidence_window)-coincidence_window;
    isample0 = (isample0/VCLReal::num_real)*VCLReal::num_real;
    for(unsigned isample = isample0; isample<nsample_; isample+=VCLReal::num_real) {
      std::fill(ntriggered, ntriggered+VCLReal::num_real, 0);
      std::fill(new_triggers, new_triggers+VCLReal::num_real, 0);

      uint_vt thv = 0;
      for(unsigned ipix=0; ipix<npix_; ipix+=VCLReal::num_real) {
        const unsigned imask = ipix/mask_t_size_bits;
        const unsigned ishift = ((ipix - imask*mask_t_size_bits)/VCLReal::num_real)*VCLReal::num_real;

        real_vt pix_threshold = 0;
        if(threshold.size() == 1) {
          pix_threshold = threshold[0];
        } else if(ipix + VCLReal::num_real <= npix_) {
          pix_threshold.load_a(threshold.data() + ipix);
        } else {
          real_at threshold_array;
          pix_threshold.store_a(threshold_array);
          for(unsigned i=0; i+ipix<npix_; i++) {
            threshold_array[i] = threshold(ipix+i);
          }
          pix_threshold.load_a(threshold_array);
        }
        pix_threshold = select(VCLReal::int_iota()+ipix < npix_, pix_threshold, std::numeric_limits<real_t>::infinity());

        int_vt cwin;
        cwin.load_a(cwin_tend + ipix);
        if(v_waveform_is_packed_) {
          const real_t *__restrict__ v_waveform_ptr = v_waveform_.data() + ipix*nsample_ + isample*VCLReal::num_real;
          for(unsigned jsample=0;jsample<VCLReal::num_real;jsample++) {
            real_vt v;
            v.load_a(v_waveform_ptr);
            v_waveform_ptr += VCLReal::num_real;
            // ****************************************************************
            // NOTE THAT THIS DUPICATES CODE IN UNPACKED SECTION BELOW
            // ****************************************************************
            int ksample = isample+jsample;
            int_vt new_cwin = vcl::select(int_bvt(v > pix_threshold), ksample+coincidence_window, cwin);
            int_bvt triggered_pix = ksample < new_cwin;
            mask_t triggered_pix_mask = vcl::to_bits(triggered_pix);
            new_triggers[jsample] |= vcl::horizontal_or(triggered_pix && (ksample-1 >= cwin));
            ntriggered[jsample] += std::popcount(triggered_pix_mask);
            thv |= vcl::select(VCLReal::uint_iota() == jsample, uint_vt(triggered_pix_mask << ishift), uint_vt(0U));
            cwin = new_cwin;
          }
        } else {
          const real_t *__restrict__ v_waveform_ptr = v_waveform_.data() + ipix*nsample_ + isample;
          real_vt block[VCLReal::num_real]; // square matrix of reals
          for(unsigned jpix=0;jpix<VCLReal::num_real;jpix++) {
            block[jpix].load_a(v_waveform_ptr);
            v_waveform_ptr += nsample_;
          }
          calin::util::vcl::transpose(block);
          for(unsigned jsample=0;jsample<VCLReal::num_real;jsample++) {
            // ****************************************************************
            // NOTE THAT THIS DUPICATES CODE IN PACKED SECTION ABOVE
            // ****************************************************************
            int ksample = isample+jsample;
            int_vt new_cwin = vcl::select(int_bvt(block[jsample] > pix_threshold), ksample+coincidence_window, cwin);
            int_bvt triggered_pix = ksample < new_cwin;
            mask_t triggered_pix_mask = vcl::to_bits(triggered_pix);
            new_triggers[jsample] |= vcl::horizontal_or(triggered_pix && (ksample-1 >= cwin));
            ntriggered[jsample] += std::popcount(triggered_pix_mask);
            thv |= vcl::select(VCLReal::uint_iota() == jsample, uint_vt(triggered_pix_mask << ishift), uint_vt(0U));
            cwin = new_cwin;
          }
        }
        cwin.store_a(cwin_tend + ipix);

        if(ishift == mask_t_size_bits-VCLReal::num_real or ipix+VCLReal::num_real>=npix_) {
          thv.store_a(trigger_hit_array + imask*VCLReal::num_real);
          thv = uint_vt(0U);
        }
      }

      for(unsigned jsample=0;jsample<VCLReal::num_real;jsample++) {
        unsigned ksample = isample+jsample;
        if(ksample>=first_sample_of_interest and ntriggered[jsample]>=3 and new_triggers[jsample]) {          
          for(unsigned ith=0; ith<trigger_hit_size; ith++) {
            mask_t th = trigger_hit_array[ith*VCLReal::num_real + jsample];
            const unsigned ipix0 = ith*mask_t_size_bits;
            while(th!=0) {
              const int ihit = std::countr_zero(th);
              const int ipix = ipix0 + ihit;
              th &= th - 1;
              unsigned nneighbor_hit = 0;
              for(auto jpix : neighbors.col(ipix)) {
                if(jpix>=0 && jpix<int(npix_)) {
                  const unsigned jmask = jpix/mask_t_size_bits;
                  const unsigned jshift = jpix - jmask*mask_t_size_bits;
                  if(trigger_hit_array[jmask*VCLReal::num_real + jsample] & (mask_t(1) << jshift)) {
                    ++nneighbor_hit;
                    if(nneighbor_hit >= 2) {
                      // Found 3 neighbor channels that triggered
                      // calin::util::log::LOG(calin::util::log::INFO) << ipix << ' ' << jpix << " ( " << neighbors.col(ipix).transpose() << ')';
                      ::free(cwin_tend);
                      ::free(trigger_hit_array);
                      return isample+jsample;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    ::free(cwin_tend);
    ::free(trigger_hit_array);
    return -1;
  }

  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  //
  //      j88D  d8b   db d8b   db 
  //     j8~88  888o  88 888o  88 
  //    j8' 88  88V8o 88 88V8o 88 
  //    V88888D 88 V8o88 88 V8o88 
  //        88  88  V888 88  V888 
  //        VP  VP   V8P VP   V8P 
  //
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  
  int trigger_4nn(const vecX_t& threshold, const Eigen::MatrixXi& neighbors, unsigned coincidence_window, unsigned first_sample_of_interest=0)
  {
    // Vectorized, unrolled, 4NN trigger algorithm that supports coincidence window
    // 1) For each sample: find which pixels are above threshold and set end of their coincidence window
    // 2) Record whether any pixel is newly triggered (was outside window in last sample, now inside)
    // 3) Count number of pixels triggered (i.e. within their coincidence windows)
    // 4) If there are more than 4 triggered pixels, and at least one is newly triggered then run 4NN algorithm:
    // 5) Run through triggered bitmask finding triggered channels, then
    // 6) Test whether neighbors are also triggered, counting how many, and recording first two triggered neighbors
    // 7) If 3 neighbors were found then return trigger time
    // 8) Otherwise, if 2 triggered neighbors were found then test the neighbors of these that haven't already been tested
    // 9) If one is found to be triggered then return trigger time

    using mask_t = typename VCLReal::uint_t;
    constexpr unsigned mask_t_size_bits = sizeof(mask_t)*8;
    static_assert(VCLReal::num_real <= mask_t_size_bits);
    static_assert(mask_t_size_bits % VCLReal::num_real == 0);

    if(threshold.size()!=1 and threshold.size()!=npix_) {
      throw std::domain_error("Trigger threshold vector length is not equal to one or number of pixels " 
        + std::to_string(threshold.size()) + " != 1 or " + std::to_string(npix_));
    }

    if(neighbors.cols()!=npix_) {
      throw std::domain_error("Number of columns in neighbors matrix must equal number of pixels " 
        + std::to_string(neighbors.cols()) + " != " + std::to_string(npix_));
    }

    coincidence_window = std::max(coincidence_window, 1U); // Zero makes no sense
    
    int_t *__restrict__ cwin_tend = calin::util::memory::aligned_calloc<int_t>(v_waveform_.cols());
    std::fill(cwin_tend, cwin_tend+v_waveform_.cols(), -1);

    unsigned ntriggered[VCLReal::num_real];
    bool new_triggers[VCLReal::num_real];

    unsigned trigger_hit_size = (v_waveform_.cols() + mask_t_size_bits - 1)/mask_t_size_bits;
    mask_t *__restrict__ trigger_hit_array = 
      calin::util::memory::aligned_calloc<mask_t>(trigger_hit_size * VCLReal::num_real);

    unsigned isample0 = std::max(first_sample_of_interest,coincidence_window)-coincidence_window;
    isample0 = (isample0/VCLReal::num_real)*VCLReal::num_real;
    for(unsigned isample = isample0; isample<nsample_; isample+=VCLReal::num_real) {
      std::fill(ntriggered, ntriggered+VCLReal::num_real, 0);
      std::fill(new_triggers, new_triggers+VCLReal::num_real, 0);

      uint_vt thv = 0;
      for(unsigned ipix=0; ipix<npix_; ipix+=VCLReal::num_real) {
        const unsigned imask = ipix/mask_t_size_bits;
        const unsigned ishift = ((ipix - imask*mask_t_size_bits)/VCLReal::num_real)*VCLReal::num_real;

        real_vt pix_threshold = 0;
        if(threshold.size() == 1) {
          pix_threshold = threshold[0];
        } else if(ipix + VCLReal::num_real <= npix_) {
          pix_threshold.load_a(threshold.data() + ipix);
        } else {
          real_at threshold_array;
          pix_threshold.store_a(threshold_array);
          for(unsigned i=0; i+ipix<npix_; i++) {
            threshold_array[i] = threshold(ipix+i);
          }
          pix_threshold.load_a(threshold_array);
        }
        pix_threshold = select(VCLReal::int_iota()+ipix < npix_, pix_threshold, std::numeric_limits<real_t>::infinity());

        int_vt cwin;
        cwin.load_a(cwin_tend + ipix);
        if(v_waveform_is_packed_) {
          const real_t *__restrict__ v_waveform_ptr = v_waveform_.data() + ipix*nsample_ + isample*VCLReal::num_real;
          for(unsigned jsample=0;jsample<VCLReal::num_real;jsample++) {
            real_vt v;
            v.load_a(v_waveform_ptr);
            v_waveform_ptr += VCLReal::num_real;
            // ****************************************************************
            // NOTE THAT THIS DUPICATES CODE IN UNPACKED SECTION BELOW
            // ****************************************************************
            int ksample = isample+jsample;
            int_vt new_cwin = vcl::select(int_bvt(v > pix_threshold), ksample+coincidence_window, cwin);
            int_bvt triggered_pix = ksample < new_cwin;
            mask_t triggered_pix_mask = vcl::to_bits(triggered_pix);
            new_triggers[jsample] |= vcl::horizontal_or(triggered_pix && (ksample-1 >= cwin));
            ntriggered[jsample] += std::popcount(triggered_pix_mask);
            thv |= vcl::select(VCLReal::uint_iota() == jsample, uint_vt(triggered_pix_mask << ishift), uint_vt(0U));
            cwin = new_cwin;
          }
        } else {
          const real_t *__restrict__ v_waveform_ptr = v_waveform_.data() + ipix*nsample_ + isample;
          real_vt block[VCLReal::num_real]; // square matrix of reals
          for(unsigned jpix=0;jpix<VCLReal::num_real;jpix++) {
            block[jpix].load_a(v_waveform_ptr);
            v_waveform_ptr += nsample_;
          }
          calin::util::vcl::transpose(block);
          for(unsigned jsample=0;jsample<VCLReal::num_real;jsample++) {
            // ****************************************************************
            // NOTE THAT THIS DUPICATES CODE IN PACKED SECTION ABOVE
            // ****************************************************************
            int ksample = isample+jsample;
            int_vt new_cwin = vcl::select(int_bvt(block[jsample] > pix_threshold), ksample+coincidence_window, cwin);
            int_bvt triggered_pix = ksample < new_cwin;
            mask_t triggered_pix_mask = vcl::to_bits(triggered_pix);
            new_triggers[jsample] |= vcl::horizontal_or(triggered_pix && (ksample-1 >= cwin));
            ntriggered[jsample] += std::popcount(triggered_pix_mask);
            thv |= vcl::select(VCLReal::uint_iota() == jsample, uint_vt(triggered_pix_mask << ishift), uint_vt(0U));
            cwin = new_cwin;
          }
        }
        cwin.store_a(cwin_tend + ipix);

        if(ishift == mask_t_size_bits-VCLReal::num_real or ipix+VCLReal::num_real>=npix_) {
          thv.store_a(trigger_hit_array + imask*VCLReal::num_real);
          thv = uint_vt(0U);
        }
      }

      for(unsigned jsample=0;jsample<VCLReal::num_real;jsample++) {
        unsigned ksample = isample+jsample;
        if(ksample>=first_sample_of_interest and ntriggered[jsample]>=4 and new_triggers[jsample]) {          
          for(unsigned ith=0; ith<trigger_hit_size; ith++) {
            mask_t th = trigger_hit_array[ith*VCLReal::num_real + jsample];
            const unsigned ipix0 = ith*mask_t_size_bits;
            while(th!=0) {
              const int ihit = std::countr_zero(th);
              const int ipix = ipix0 + ihit;
              th &= th - 1;
              unsigned nneighbor_hit = 0;
              int ipix_n1 = -1;
              int ipix_n2 = -1;
              for(auto jpix : neighbors.col(ipix)) {
                if(jpix>=0 && jpix<int(npix_)) {
                  const unsigned jmask = jpix/mask_t_size_bits;
                  const unsigned jshift = jpix - jmask*mask_t_size_bits;
                  if(trigger_hit_array[jmask*VCLReal::num_real + jsample] & (mask_t(1) << jshift)) {
                    ++nneighbor_hit;
                    if(nneighbor_hit >= 3) {
                      // Found 4 neighbor channels that triggered
                      // calin::util::log::LOG(calin::util::log::INFO) << ipix << ' ' << jpix << " ( " << neighbors.col(ipix).transpose() << ')';
                      ::free(cwin_tend);
                      ::free(trigger_hit_array);
                      return isample+jsample;
                    } else if(nneighbor_hit == 2) {
                      ipix_n2 = jpix;
                    } else {
                      ipix_n1 = jpix;
                    }                    
                  }
                }
              }
              if(nneighbor_hit == 2) {
                // Found 3 neighbor channels that triggered, try to find 4th
                for(auto jpix : neighbors.col(ipix_n1)) {
                  if(jpix>=0 and jpix<int(npix_) and jpix!=ipix and jpix!=ipix_n2) {
                    const unsigned jmask = jpix/mask_t_size_bits;
                    const unsigned jshift = jpix - jmask*mask_t_size_bits;
                    if(trigger_hit_array[jmask*VCLReal::num_real + jsample] & (mask_t(1) << jshift)) {
                      // Found a fourth triggered channel
                      ::free(cwin_tend);
                      ::free(trigger_hit_array);
                      return isample+jsample;
                    }
                  }
                }
                for(auto jpix : neighbors.col(ipix_n2)) {
                  if(jpix>=0 and jpix<int(npix_) and jpix!=ipix and jpix!=ipix_n1) {
                    const unsigned jmask = jpix/mask_t_size_bits;
                    const unsigned jshift = jpix - jmask*mask_t_size_bits;
                    if(trigger_hit_array[jmask*VCLReal::num_real + jsample] & (mask_t(1) << jshift)) {
                      // Found a fourth triggered channel
                      ::free(cwin_tend);
                      ::free(trigger_hit_array);
                      return isample+jsample;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    ::free(cwin_tend);
    ::free(trigger_hit_array);
    return -1;
  }

  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  //
  //    d888888b d8b   db    d88b d88888b  .o88b. d888888b 
  //      `88'   888o  88    `8P' 88'     d8P  Y8 `~~88~~' 
  //       88    88V8o 88     88  88ooooo 8P         88    
  //       88    88 V8o88     88  88~~~~~ 8b         88    
  //      .88.   88  V888 db. 88  88.     Y8b  d8    88    
  //    Y888888P VP   V8P Y8888P  Y88888P  `Y88P'    YP    
  //
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  void inject_pe(unsigned pixel_id, unsigned t_sample, double amplitude = 1.0) 
  {
    if(pixel_id>=npix_) {
      throw std::out_of_range("Pixel index out of range in inject_pe: " 
        + std::to_string(pixel_id) + " >= " + std::to_string(npix_));
    }
    if(t_sample>=nsample_) {
      throw std::out_of_range("Sample time index out of range in inject_pe: " 
        + std::to_string(t_sample) + " >= " + std::to_string(nsample_));
    }
    pe_waveform_(t_sample, pixel_id) += amplitude;
  }

  void inject_n_pes(unsigned pixel_id, unsigned npe, double t0_samples,
    calin::simulation::detector_efficiency::VCLPEAmplitudeGenerator<VCLArchitecture>* pegen = nullptr,
    double time_spread_ns = 0.0)
  {
    double time_spread_samples = time_spread_ns*sampling_freq_ghz_;
    real_t *__restrict__ pe_waveform_ptr = &pe_waveform_(0,pixel_id);
    while(npe) {
      typename VCLArchitecture::double_vt t = t0_samples;
      if(time_spread_ns > 0) {
        t += rng_->normal_double() * time_spread_samples;
      } else if (time_spread_ns < 0) {
        t -= rng_->uniform_double_zc(time_spread_samples);
      }        
      typename VCLArchitecture::double_vt q;
      if(pegen) {
        q = pegen->vcl_generate_amplitude(*rng_);
      } else {
        q = 1.0;
      }

      typename VCLArchitecture::double_at q_a;
      q.store_a(q_a);

      typename VCLArchitecture::int64_vt it = round_to_int64_limited(t);
      typename VCLArchitecture::int64_at it_a;
      it.store_a(it_a);
      for(unsigned jvec=0;jvec<VCLArchitecture::num_double && npe; jvec++,npe--) {
        unsigned isample = it_a[jvec];
        if(isample>=0 and isample<nsample_) {
          pe_waveform_ptr[isample] += q_a[jvec];
        }
      }
    }
  }

  void inject_poisson_pes(unsigned pixel_id, double lambda, double t0_samples,
    calin::simulation::detector_efficiency::VCLPEAmplitudeGenerator<VCLArchitecture>* pegen = nullptr,
    double time_spread_ns = 0.0)
  {
    calin::math::rng::VCLToScalarRNGCore scalar_core(rng_->core());
    calin::math::rng::RNG scalar_rng(&scalar_core);
    unsigned npe = scalar_rng.poisson(lambda);
    inject_n_pes(pixel_id, npe, t0_samples, pegen, time_spread_ns);
  }

  void inject_poisson_pes(const Eigen::VectorXd& lambda, const Eigen::VectorXd& t0_samples, 
    calin::simulation::detector_efficiency::VCLPEAmplitudeGenerator<VCLArchitecture>* pegen = nullptr,
    const Eigen::VectorXd& time_spread_ns = Eigen::VectorXd())
  {
    if(lambda.size() != npix_) {
      throw std::domain_error("Lambda vector length is not equal to number of pixels " 
        + std::to_string(lambda.size()) + " != " + std::to_string(npix_));
    }
    if(t0_samples.size() != 1 and t0_samples.size() != npix_) {
      throw std::domain_error("T0 vector length must be 1 or equal to number of pixels " 
        + std::to_string(t0_samples.size()) + " != 1 or " + std::to_string(npix_));
    }
    if(time_spread_ns.size() != 0 and time_spread_ns.size() != 1 and time_spread_ns.size() != npix_) {
      throw std::domain_error("Time spread vector length must be 0, 1 or equal to number of pixels " 
        + std::to_string(time_spread_ns.size()) + " != 0, 1 or " + std::to_string(npix_));
    }

    calin::math::rng::RNG& scalar_rng(rng_->scalar_rng());
    for(unsigned ipix=0; ipix<npix_; ipix++) {
      unsigned npe = scalar_rng.poisson(lambda[ipix]);
      double t0 = t0_samples[0];
      double time_spread = 0;
      if(t0_samples.size() > 1) {
        t0 = t0_samples[ipix];
      }
      if(time_spread_ns.size() > 1) {
        time_spread = time_spread_ns[ipix];
      } else if (time_spread_ns.size() == 1) {
        time_spread = time_spread_ns[0];
      }
      inject_n_pes(ipix, npe, t0, pegen, time_spread);
    }
  }

  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  //
  //     .o88b.  .d8b.  .88b  d88.      d8888b. d88888b .d8888. d8888b. 
  //    d8P  Y8 d8' `8b 88'YbdP`88      88  `8D 88'     88'  YP 88  `8D 
  //    8P      88ooo88 88  88  88      88oobY' 88ooooo `8bo.   88oodD' 
  //    8b      88~~~88 88  88  88      88`8b   88~~~~~   `Y8b. 88~~~   
  //    Y8b  d8 88   88 88  88  88      88 `88. 88.     db   8D 88      
  //     `Y88P' YP   YP YP  YP  YP      88   YD Y88888P `8888Y' 88      
  //
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  unsigned add_camera_response(const Eigen::VectorXi& impulse_response_id, 
    bool add_ac_coupling_offset = false) 
  {
    if(impulse_response_id.size() != 1 and impulse_response_id.size() != npix_) {
      throw std::domain_error("Impulse response id must be scalar or vector of length equal to number of pixels " 
        + std::to_string(impulse_response_id.size()) + " != 1 or " + std::to_string(npix_));
    }

    for(auto irid : impulse_response_id) {
      validate_impulse_response_id(irid);
    }

    CameraResponse cr;
    cr.impulse_response_id       = impulse_response_id;
    cr.add_ac_coupling_offset    = add_ac_coupling_offset;

    camera_responses_.emplace_back(std::move(cr));
    return camera_responses_.size() - 1;
  }

  void set_cr_nsb_rate(unsigned camera_response_id, 
    const Eigen::VectorXd& nsb_freq_per_pixel_ghz,
    calin::simulation::detector_efficiency::VCLPEAmplitudeGenerator<VCLArchitecture>* pegen = nullptr,
    bool adopt_pegen = false)
  {
    validate_camera_response_id(camera_response_id);
    auto& cr = camera_responses_[camera_response_id];

    if(nsb_freq_per_pixel_ghz.size() != int(npix_)) {
      throw std::domain_error("NSB frequency vector length is not equal to number of pixels " 
        + std::to_string(nsb_freq_per_pixel_ghz.size()) + " != " + std::to_string(npix_));
    }

    cr.nsb_freq_per_pixel_ghz    = nsb_freq_per_pixel_ghz;
    
    if(pegen!=nullptr) {
      if(cr.nsb_pegen!=nullptr and cr.adopt_nsb_pegen) {
        delete cr.nsb_pegen;
      }
      cr.nsb_pegen                   = pegen;
      cr.adopt_nsb_pegen             = adopt_pegen;
    }

    cr.recalc_pedestal(this);
  }

  void set_cr_relative_gain(unsigned camera_response_id, const vecX_t& relative_gain)
  {
    validate_camera_response_id(camera_response_id);
    auto& cr = camera_responses_[camera_response_id];

    if(relative_gain.size()!=0 and relative_gain.size()!=npix_) {
      throw std::domain_error("Relative gain vector length is not equal to number of pixels " 
        + std::to_string(relative_gain.size()) + " != " + std::to_string(npix_));
    }

    cr.relative_gain             = relative_gain;
    cr.recalc_pedestal(this);
  }

  void set_cr_pedestal(unsigned camera_response_id, const vecX_t& pedestal)
  {
    validate_camera_response_id(camera_response_id);
    auto& cr = camera_responses_[camera_response_id];

    if(pedestal.size()!=0 and pedestal.size()!=npix_) {
      throw std::domain_error("Pedestal vector length is not equal to number of pixels " 
        + std::to_string(pedestal.size()) + " != " + std::to_string(npix_));
    }

    cr.demand_pedestal           = pedestal;
    cr.recalc_pedestal(this);

  }

  void set_cr_noise_spectrum(unsigned camera_response_id, const matX_t& noise_spectrum)
  {
    validate_camera_response_id(camera_response_id);
    auto& cr = camera_responses_[camera_response_id];

    if(noise_spectrum.size()!=0 and (noise_spectrum.rows()!=nsample_
        or (noise_spectrum.cols()!=1 and noise_spectrum.cols()!=npix_))) {
      throw std::domain_error("Noise spectrum matrix has incoorect size: [" 
        + std::to_string(noise_spectrum.rows()) + "," + std::to_string(noise_spectrum.cols()) + "] != ["
        + std::to_string(nsample_) + ",1] or ["
        + std::to_string(nsample_) + "," + std::to_string(npix_) + "]");
    }

    cr.noise_spectrum            = noise_spectrum;
  }

  void set_cr_neighbors(unsigned camera_response_id, const Eigen::MatrixXi& neighbors)
  {
    validate_camera_response_id(camera_response_id);
    auto& cr = camera_responses_[camera_response_id];

    if(neighbors.cols()!=npix_) {
      throw std::domain_error("Number of columns in neighbors matrix must equal number of pixels " 
        + std::to_string(neighbors.cols()) + " != " + std::to_string(npix_));
    }

    cr.neighbors                 = neighbors;
  }

  void set_cr_neighbors(unsigned camera_response_id, 
    const calin::ix::iact_data::instrument_layout::CameraLayout& camera_layout)
  {
    set_cr_neighbors(camera_response_id, calin::iact_data::instrument_layout::neighbor_matrix(camera_layout));
  }

  void set_cr_multiplicity(unsigned camera_response_id, unsigned multiplicity)
  {
    validate_camera_response_id(camera_response_id);
    auto& cr = camera_responses_[camera_response_id];

    cr.multiplicity              = multiplicity;
  }

  void set_cr_threshold(unsigned camera_response_id, const vecX_t& threshold, unsigned conincidence_window = 0)
  {
    validate_camera_response_id(camera_response_id);
    auto& cr = camera_responses_[camera_response_id];

    if(threshold.size()!=1 and threshold.size()!=npix_) {
      throw std::domain_error("Trigger threshold vector length is not equal to one or number of pixels " 
        + std::to_string(threshold.size()) + " != 1 or " + std::to_string(npix_));
    }

    cr.threshold                 = threshold;
    if(conincidence_window > 0) {
      cr.coincidence_window      = conincidence_window;
    }
  }

  Eigen::VectorXd get_cr_nsb_rate(unsigned camera_response_id)
  {
    validate_camera_response_id(camera_response_id);
    auto& cr = camera_responses_[camera_response_id];
    return cr.nsb_freq_per_pixel_ghz;
  }

  calin::simulation::detector_efficiency::VCLPEAmplitudeGenerator<VCLArchitecture>* get_cr_nsb_pegen(unsigned camera_response_id)
  {
    validate_camera_response_id(camera_response_id);
    auto& cr = camera_responses_[camera_response_id];
    return cr.nsb_pegen;
  }

  Eigen::VectorXi get_cr_impulse_response_id(unsigned camera_response_id)
  {
    validate_camera_response_id(camera_response_id);
    auto& cr = camera_responses_[camera_response_id];
    return cr.impulse_response_id;
  }

  vecX_t get_cr_pedestal(unsigned camera_response_id)
  {
    validate_camera_response_id(camera_response_id);
    auto& cr = camera_responses_[camera_response_id];
    return cr.demand_pedestal;
  }
  
  vecX_t get_cr_offset(unsigned camera_response_id)
  {
    validate_camera_response_id(camera_response_id);
    auto& cr = camera_responses_[camera_response_id];
    return cr.pedestal;
  }

  vecX_t get_cr_relative_gain(unsigned camera_response_id)
  {
    validate_camera_response_id(camera_response_id);
    auto& cr = camera_responses_[camera_response_id];
    return cr.relative_gain;
  }

  matX_t get_cr_noise_spectrum(unsigned camera_response_id)
  {
    validate_camera_response_id(camera_response_id);
    auto& cr = camera_responses_[camera_response_id];
    return cr.noise_spectrum;
  }

  vecX_t get_cr_threshold(unsigned camera_response_id)
  {
    validate_camera_response_id(camera_response_id);
    auto& cr = camera_responses_[camera_response_id];
    return cr.threshold;
  }

  Eigen::MatrixXi get_cr_neighbors(unsigned camera_response_id)
  {
    validate_camera_response_id(camera_response_id);
    auto& cr = camera_responses_[camera_response_id];
    return cr.neighbors;
  }

  unsigned get_cr_multiplicity(unsigned camera_response_id)
  {
    validate_camera_response_id(camera_response_id);
    auto& cr = camera_responses_[camera_response_id];
    return cr.multiplicity;
  }

  unsigned get_cr_coincidence_window(unsigned camera_response_id)
  {
    validate_camera_response_id(camera_response_id);
    auto& cr = camera_responses_[camera_response_id];
    return cr.coincidence_window;
  }

  void add_nsb_noise_to_waveform_cr(unsigned camera_response_id, double t0_samples = 0)
  {
    validate_camera_response_id(camera_response_id);
    auto& cr = camera_responses_[camera_response_id];

    add_nsb_noise_to_waveform(cr.nsb_freq_per_pixel_ghz, cr.nsb_pegen, t0_samples);
  }

  void convolve_impulse_response_fftw_codelet_cr(unsigned camera_response_id)
  {
    validate_camera_response_id(camera_response_id);
    auto& cr = camera_responses_[camera_response_id];

    if(cr.impulse_response_id.size() == 1) {
      convolve_impulse_response_fftw_codelet(
        cr.impulse_response_id[0], cr.pedestal, cr.relative_gain, cr.noise_spectrum);
    } else {
      convolve_multiple_impulse_response_fftw_codelet(
        cr.impulse_response_id, cr.pedestal, cr.relative_gain, cr.noise_spectrum);
    }
  }

  unsigned trigger_singles_cr(Eigen::VectorXi& pixel_trigger_time, unsigned camera_response_id, unsigned first_sample_of_interest=0){
    validate_camera_response_id(camera_response_id);
    auto& cr = camera_responses_[camera_response_id];

    return trigger_singles(pixel_trigger_time, cr.threshold, first_sample_of_interest);
  }

  int trigger_multiplicity_cr(unsigned camera_response_id, unsigned first_sample_of_interest=0) 
  {
    validate_camera_response_id(camera_response_id);
    auto& cr = camera_responses_[camera_response_id];

    return trigger_multiplicity(cr.threshold, cr.multiplicity, cr.coincidence_window, 
      first_sample_of_interest);
  }

  int trigger_panoseti_multiplicity_cr(unsigned camera_response_id, unsigned first_sample_of_interest=0)
  {
    validate_camera_response_id(camera_response_id);
    auto& cr = camera_responses_[camera_response_id];

    return trigger_panoseti_multiplicity(cr.threshold, cr.multiplicity, cr.coincidence_window, first_sample_of_interest);
  }

  int trigger_3nn_cr(unsigned camera_response_id, unsigned first_sample_of_interest=0)
  {
    validate_camera_response_id(camera_response_id);
    auto& cr = camera_responses_[camera_response_id];

    return trigger_3nn(cr.threshold, cr.neighbors, cr.coincidence_window, first_sample_of_interest);
  }

  int trigger_4nn_cr(unsigned camera_response_id, unsigned first_sample_of_interest=0)
  {
    validate_camera_response_id(camera_response_id);
    auto& cr = camera_responses_[camera_response_id];

    return trigger_4nn(cr.threshold, cr.neighbors, cr.coincidence_window, first_sample_of_interest);
  }

  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  //
  //    .88b  d88. d888888b .d8888.  .o88b. 
  //    88'YbdP`88   `88'   88'  YP d8P  Y8 
  //    88  88  88    88    `8bo.   8P      
  //    88  88  88    88      `Y8b. 8b      
  //    88  88  88   .88.   db   8D Y8b  d8 
  //    YP  YP  YP Y888888P `8888Y'  `Y88P' 
  //
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  vecX_t ac_coupling_offset(const Eigen::VectorXi& impulse_response_id, 
    const Eigen::VectorXd& nsb_freq_per_pixel_ghz,
    calin::simulation::detector_efficiency::VCLPEAmplitudeGenerator<VCLArchitecture>* pegen = nullptr,
    const vecX_t& relative_gain = vecX_t())
  {
    int iri_max = std::numeric_limits<int>::min();
    int iri_min = std::numeric_limits<int>::max();
    for(unsigned ipix=0; ipix<impulse_response_id.size(); ipix++) {
      int iri = impulse_response_id(ipix);
      iri_max = std::max(iri_max, iri);
      iri_min = std::min(iri_min, iri);
    }

    validate_impulse_response_id(iri_max);
    if(iri_min<0) {
      throw std::out_of_range("Invalid impulse response id cannot be negative: " 
        + std::to_string(iri_min));
    }

    if(nsb_freq_per_pixel_ghz.size() != impulse_response_id.size()) {
      throw std::domain_error("NSB frequency and impulse response vectors must have same size " 
        + std::to_string(nsb_freq_per_pixel_ghz.size()) + " != " + std::to_string(impulse_response_id.size()));
    }

    if(relative_gain.size()!=0 and relative_gain.size()!=impulse_response_id.size()) {
      throw std::domain_error("Relative gain and impulse response vectors must have same lengths " 
        + std::to_string(relative_gain.size()) + " != " + std::to_string(impulse_response_id.size()));
    }

    double pe_amp = 1.0;
    if(pegen) {
      pe_amp = pegen->mean_amplitude();
    }

    vecX_t offset = (nsb_freq_per_pixel_ghz * time_resolution_ns_ * pe_amp).cast<real_t>();

    for(unsigned iir=0; iir<impulse_response_id.size(); iir++) {
      offset[iir] *= impulse_responses_[impulse_response_id[iir]].response.sum();
    }
    
    if(relative_gain.size()) {
      offset *= relative_gain;
    }

    return offset;
  }

  vecX_t ac_coupling_offset(int impulse_response_id, 
    const Eigen::VectorXd& nsb_freq_per_pixel_ghz,
    calin::simulation::detector_efficiency::VCLPEAmplitudeGenerator<VCLArchitecture>* pegen = nullptr,
    const vecX_t& relative_gain = vecX_t())
  {
    Eigen::VectorXi all_impulse_response_id(nsb_freq_per_pixel_ghz.size());
    all_impulse_response_id.setConstant(impulse_response_id);
    return ac_coupling_offset(all_impulse_response_id, nsb_freq_per_pixel_ghz, pegen,relative_gain);
  }

  void clear_waveforms()
  {
    pe_waveform_.setZero();
    v_waveform_.setZero();
    v_waveform_is_packed_ = false;
  }

  void shift_and_clear_waveforms(unsigned nsample)
  {
    if(nsample == 0) {
      clear_waveforms();
      return;
    }
    pe_waveform_.topRows(nsample) = pe_waveform_.bottomRows(nsample);
    pe_waveform_.bottomRows(nsample_-nsample).setZero();
    v_waveform_.setZero();
    v_waveform_is_packed_ = false;
  }

  void clear() override
  {
    SimpleListPEProcessor::clear();
    clear_waveforms();
  }

  std::string fftw_plans_summary() const
  {
    double* a_vec = fftw_alloc_real(nsample_);
    double* b_vec = fftw_alloc_real(nsample_);

    auto fwd = calin::math::fftw_util::plan_r2hc(nsample_, a_vec, b_vec, FFTW_MEASURE);
    auto bwd = calin::math::fftw_util::plan_hc2r(nsample_, a_vec, b_vec, FFTW_MEASURE);

    std::string os;
    os += "Photo-electron forward transformation:\n";
    os += "--------------------------------------\n";
    os += fwd.print();

    os += "\nAmplitude backward transformation:\n";
    os += "----------------------------------\n";
    os += bwd.print();

    fftw_free(a_vec);
    fftw_free(b_vec);

    return os;
  }
  
  real_t noise_spectrum_var(const vecX_t& noise_spectrum) const 
  {
    using calin::math::special::SQR;
    if(noise_spectrum.size()!=nsample_) {
      throw std::domain_error("Noise spectrum vector length is not equal to number of samples " 
        + std::to_string(noise_spectrum.size()) + " != " + std::to_string(nsample_));
    }

    const real_t* r = noise_spectrum.data();
    const real_t* c = noise_spectrum.data() + nsample_;
    real_t var = SQR(*r++);
    c--;
    while(r<c) {
      var += 2.0*(SQR(*r++) + SQR(*c--));
    }
    if(r==c) {
      var += SQR(*r++);
    }
    return var;
  }

  Eigen::VectorXd spectral_frequencies_ghz(bool imaginary_negative = false) const 
  {
    Eigen::VectorXd freq(nsample_);
    calin::math::fftw_util::hcvec_fftfreq(freq.data(), nsample_, time_resolution_ns_, imaginary_negative);
    return freq;
  }

  matX_t pe_waveform() const
  { 
    return pe_waveform_.block(0, 0, nsample_, npix_);
  }
  
  matX_t v_waveform() const 
  { 
    if(!v_waveform_is_packed_) {
      return v_waveform_.block(0, 0, nsample_, npix_);
    }

    matX_t wf(nsample_, npix_);
    const real_t *__restrict__ v_waveform_ptr = v_waveform_.data();
    for(unsigned ipix=0; ipix<npix_; ipix+=VCLReal::num_real) {
      real_t *__restrict__ wf_ptr = wf.data() + ipix*nsample_;
      for(unsigned isample = 0; isample<nsample_; isample+=VCLReal::num_real) {
        real_vt block[VCLReal::num_real]; // square matrix of reals
        for(unsigned jsample=0;jsample<VCLReal::num_real;jsample++) {
          block[jsample].load_a(v_waveform_ptr);
          v_waveform_ptr += VCLReal::num_real;
        }
        calin::util::vcl::transpose(block);
        for(unsigned jpix=0,zpix=std::min(npix_-ipix,VCLReal::num_real);jpix<zpix;jpix++) {
          block[jpix].store_a(wf_ptr + jpix*nsample_);
        }
        wf_ptr += VCLReal::num_real;
      }
    }

    return wf; 
  }


private:
#ifndef SWIG
  static inline unsigned round_nreal_to_vector(unsigned n) {
    return VCLReal::template round_nreal_up_to_vector_size<real_t>(n);
  } 

  inline void validate_impulse_response_id(unsigned impulse_response_id) const {
    if(impulse_response_id >= impulse_responses_.size()) {
      throw std::out_of_range("Invalid impulse response id: " 
        + std::to_string(impulse_response_id) 
        + " >= " + std::to_string(impulse_responses_.size()));
    }
  }

  inline void validate_camera_response_id(unsigned camera_response_id) const {
    if(camera_response_id >= camera_responses_.size()) {
      throw std::out_of_range("Camera impulse response id: " 
        + std::to_string(camera_response_id) 
        + " >= " + std::to_string(camera_responses_.size()));
    }
  }

  struct ImpulseResponse {
    unsigned response_size;
    vecX_t response;
    vecX_t transform;
    int response_peak_index = 0;
    double response_peak_value = -std::numeric_limits<double>::infinity();
    int response_integral_peak_index = 0;
    double response_integral_peak_value = -std::numeric_limits<double>::infinity();
    double response_window_frac;
    double response_window_lo;
    double response_window_hi;
    double response_total_integral = 0;
    std::string units = "";
  };

  struct CameraResponse {
    CameraResponse() = default;
    CameraResponse(CameraResponse&& o): 
      nsb_freq_per_pixel_ghz(o.nsb_freq_per_pixel_ghz),
      nsb_pegen(o.nsb_pegen), 
      add_ac_coupling_offset(o.add_ac_coupling_offset),
      impulse_response_id(o.impulse_response_id),
      pedestal(std::move(o.pedestal)),
      relative_gain(std::move(o.relative_gain)),
      noise_spectrum(std::move(o.noise_spectrum)),
      threshold(std::move(o.threshold)),
      neighbors(o.neighbors),
      coincidence_window(o.coincidence_window),
      adopt_nsb_pegen(o.adopt_nsb_pegen)
    {
      o.nsb_pegen = nullptr;
    }
    CameraResponse& operator=(CameraResponse&& o) {
      if(this != &o) {
        nsb_freq_per_pixel_ghz = o.nsb_freq_per_pixel_ghz;
        if(adopt_nsb_pegen) {
          delete nsb_pegen;
        }
        nsb_pegen = o.nsb_pegen;
        add_ac_coupling_offset = o.add_ac_coupling_offset;
        impulse_response_id = o.impulse_response_id;
        pedestal = std::move(o.pedestal);
        relative_gain = std::move(o.relative_gain);
        noise_spectrum = std::move(o.noise_spectrum);
        threshold = std::move(o.threshold);
        neighbors = o.neighbors;
        coincidence_window = o.coincidence_window;
        adopt_nsb_pegen = o.adopt_nsb_pegen;
        o.nsb_pegen = nullptr;
      }
      return *this;
    }
    ~CameraResponse() {
      if(adopt_nsb_pegen) {
        delete nsb_pegen;
      }
    }
    void recalc_pedestal(VCLWaveformPEProcessor* base) {
      pedestal = demand_pedestal;      
      if(add_ac_coupling_offset) {
        vecX_t offset;
        if(impulse_response_id.size() == 1) {
          offset = base->ac_coupling_offset(impulse_response_id[0], nsb_freq_per_pixel_ghz, nsb_pegen, relative_gain);
        } else {
          offset = base->ac_coupling_offset(impulse_response_id, nsb_freq_per_pixel_ghz, nsb_pegen, relative_gain);
        }
        if(pedestal.size() == 0) {
          pedestal = -offset;
        } else {
          pedestal -= offset;
        }
      }
    }
    Eigen::VectorXd nsb_freq_per_pixel_ghz;
    calin::simulation::detector_efficiency::VCLPEAmplitudeGenerator<VCLArchitecture>* nsb_pegen = nullptr;
    bool add_ac_coupling_offset = true;
    Eigen::VectorXi impulse_response_id;
    vecX_t demand_pedestal;
    vecX_t pedestal;
    vecX_t relative_gain;
    matX_t noise_spectrum;
    vecX_t threshold;
    Eigen::MatrixXi neighbors;
    unsigned multiplicity = 3;
    unsigned coincidence_window = 0;
    bool adopt_nsb_pegen = false;
  };

  unsigned nsample_;
  unsigned nadvance_;
  double time_resolution_ns_;
  double sampling_freq_ghz_;
  double time_advance_;
  double t0_ = 0;
  matX_t pe_waveform_;  // shape (nsample_, npix_): access as (it, ipix)
  matX_t v_waveform_;   // shape (nsample_, npix_): access as (it, ipix)
  bool v_waveform_is_packed_ = false;
  std::vector<ImpulseResponse> impulse_responses_;
  std::vector<CameraResponse> camera_responses_;
  RNG* rng_ = nullptr;
  bool adopt_rng_ = false;
#endif // ifndef SWIG
};


} } } // namespace calin::simulations::vcl_pe_processor
