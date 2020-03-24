/*

   calin/iact_data/waveform_treatment_event_visitor.hpp -- Stephen Fegan -- 2018-01-11

   Waveform treatment event data visitor - process waveforms

   Copyright 2018, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <string>

#include <iact_data/event_visitor.hpp>
#include <util/vcl.hpp>
#include <util/log.hpp>
#include <iact_data/waveform_treatment_event_visitor.pb.h>
#include <provenance/system_info.hpp>
#include <util/memory.hpp>

namespace calin { namespace iact_data { namespace waveform_treatment_event_visitor {

class OptimalWindowSumWaveformTreatmentEventVisitor:
  public calin::iact_data::event_visitor::ParallelEventVisitor
{
public:
  enum GainChannel { HIGH_GAIN, LOW_GAIN, SINGLE_OR_MIXED_GAIN };

  OptimalWindowSumWaveformTreatmentEventVisitor(
    calin::ix::iact_data::waveform_treatment_event_visitor::
      OptimalWindowSumWaveformTreatmentEventVisitorConfig config = default_config(),
      GainChannel gain_channel_to_treat = HIGH_GAIN);

  OptimalWindowSumWaveformTreatmentEventVisitor(
    GainChannel gain_channel_to_treat,
    calin::ix::iact_data::waveform_treatment_event_visitor::
      OptimalWindowSumWaveformTreatmentEventVisitorConfig config = default_config()):
    OptimalWindowSumWaveformTreatmentEventVisitor(config, gain_channel_to_treat)
  { /* nothing to see here */ }

  virtual ~OptimalWindowSumWaveformTreatmentEventVisitor();

  OptimalWindowSumWaveformTreatmentEventVisitor* new_sub_visitor(
    std::map<calin::iact_data::event_visitor::ParallelEventVisitor*,
        calin::iact_data::event_visitor::ParallelEventVisitor*>
      antecedent_visitors = { }) override;

  bool visit_telescope_run(
    const calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration* run_config,
    calin::iact_data::event_visitor::EventLifetimeManager* event_lifetime_manager) override;

  bool visit_telescope_event(uint64_t seq_index,
    calin::ix::iact_data::telescope_event::TelescopeEvent* event) override;

  static calin::ix::iact_data::waveform_treatment_event_visitor::
    OptimalWindowSumWaveformTreatmentEventVisitorConfig default_config()
  {
    calin::ix::iact_data::waveform_treatment_event_visitor::
      OptimalWindowSumWaveformTreatmentEventVisitorConfig cfg;
    cfg.set_integration_n(16);
    return cfg;
  }

  bool has_event() const { return has_event_; }

  unsigned nchan() const { return nchan_; }
  unsigned nsamp() const { return nsamp_; }

  std::vector<int> chan_max() const { return make_vec(chan_max_); }
  std::vector<int> chan_max_index() const { return make_vec(chan_max_index_); }
  std::vector<int> chan_bkg_win_sum() const { return make_vec(chan_bkg_win_sum_); }
  std::vector<int> chan_sig_win_sum() const { return make_vec(chan_sig_win_sum_); }
  std::vector<int> chan_opt_win_sum() const { return make_vec(chan_opt_win_sum_); }
  std::vector<int> chan_opt_win_sum_qt() const { return make_vec(chan_opt_win_sum_qt_); }
  std::vector<int> chan_opt_win_index() const { return make_vec(chan_opt_win_index_); }
  std::vector<int> chan_all_sum() const { return make_vec(chan_all_sum_); }
  std::vector<int> chan_signal_type() const { return chan_signal_type_; };

#ifndef SWIG
  const int*__restrict__ array_chan_max() const { return chan_max_; }
  const int*__restrict__ array_chan_max_index() const { return chan_max_index_; }
  const int*__restrict__ array_chan_bkg_win_sum() const { return chan_bkg_win_sum_; }
  const int*__restrict__ array_chan_sig_win_sum() const { return chan_sig_win_sum_; }
  const int*__restrict__ array_chan_opt_win_sum() const { return chan_opt_win_sum_; }
  const int*__restrict__ array_chan_opt_win_sum_qt() const { return chan_opt_win_sum_qt_; }
  const int*__restrict__ array_chan_opt_win_index() const { return chan_opt_win_index_; }
  const int*__restrict__ array_chan_all_sum_q() const { return chan_all_sum_; }
  const int*__restrict__ array_chan_signal_type() const { return chan_signal_type_.data(); };
#endif

  unsigned window_n() const { return window_n_; };
  int bkg_window_0() const { return bkg_window_0_; }
  std::vector<int> sig_window_0() const { return make_vec(sig_window_0_); }

protected:
#ifndef SWIG
  void reconfigure(unsigned nchan, unsigned nsamp);
  void scalar_analyze_waveforms(const uint16_t* __restrict__ data);

  template<typename T> std::vector<T> make_vec(const T* ptr) const {
    return std::vector<T>(ptr, ptr+nchan_);
  }

  calin::ix::iact_data::waveform_treatment_event_visitor::
    OptimalWindowSumWaveformTreatmentEventVisitorConfig config_;
  GainChannel gain_channel_to_treat_ = HIGH_GAIN;
  unsigned nchan_ = 0;
  unsigned nsamp_ = 0;
  unsigned window_n_;
  int bkg_window_0_;
  int*__restrict__ sig_window_0_ = nullptr;

  bool has_event_ = false;
  std::vector<int> chan_signal_type_;
  int*__restrict__ chan_max_ = nullptr;
  int*__restrict__ chan_max_index_ = nullptr;
  int*__restrict__ chan_bkg_win_sum_ = nullptr;
  int*__restrict__ chan_sig_win_sum_ = nullptr;
  int*__restrict__ chan_opt_win_sum_ = nullptr;
  int*__restrict__ chan_opt_win_sum_qt_ = nullptr;
  int*__restrict__ chan_opt_win_index_ = nullptr;
  int*__restrict__ chan_all_sum_ = nullptr;
#endif
};

template<typename VCLArchitecture>
class VCL_OptimalWindowSumWaveformTreatmentEventVisitor:
  public OptimalWindowSumWaveformTreatmentEventVisitor
{
public:
#ifndef SWIG
  CALIN_TYPEALIAS(int16_vt, typename VCLArchitecture::int16_vt);
  CALIN_TYPEALIAS(int16_bvt, typename VCLArchitecture::int16_bvt);
  CALIN_TYPEALIAS(uint16_vt, typename VCLArchitecture::uint16_vt);

  CALIN_TYPEALIAS(int32_vt, typename VCLArchitecture::int32_vt);
  CALIN_TYPEALIAS(int32_bvt, typename VCLArchitecture::int32_bvt);
  CALIN_TYPEALIAS(uint32_vt, typename VCLArchitecture::uint32_vt);

  constexpr static unsigned num_int16 = VCLArchitecture::num_int16;
  constexpr static unsigned num_int32 = VCLArchitecture::num_int32;
#endif

  VCL_OptimalWindowSumWaveformTreatmentEventVisitor(
      calin::ix::iact_data::waveform_treatment_event_visitor::
        OptimalWindowSumWaveformTreatmentEventVisitorConfig config = default_config(),
      GainChannel gain_channel_to_treat = HIGH_GAIN):
    OptimalWindowSumWaveformTreatmentEventVisitor(config, gain_channel_to_treat)
  {
    /* nothing to see here */
  }

  VCL_OptimalWindowSumWaveformTreatmentEventVisitor(
    GainChannel gain_channel_to_treat,
    calin::ix::iact_data::waveform_treatment_event_visitor::
      OptimalWindowSumWaveformTreatmentEventVisitorConfig config = default_config()):
    VCL_OptimalWindowSumWaveformTreatmentEventVisitor(config, gain_channel_to_treat)
  {
    /* nothing to see here */
  }

  virtual ~VCL_OptimalWindowSumWaveformTreatmentEventVisitor()
  {
    free(samples_);
  }

  VCL_OptimalWindowSumWaveformTreatmentEventVisitor<VCLArchitecture>* new_sub_visitor(
    std::map<calin::iact_data::event_visitor::ParallelEventVisitor*,
        calin::iact_data::event_visitor::ParallelEventVisitor*>
      antecedent_visitors = { }) override
  {
    return new VCL_OptimalWindowSumWaveformTreatmentEventVisitor<VCLArchitecture>(config_);
  }

  bool visit_telescope_run(
    const calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration* run_config,
    calin::iact_data::event_visitor::EventLifetimeManager* event_lifetime_manager) override
  {
    bool old_nsamp = nsamp_;
    bool good = OptimalWindowSumWaveformTreatmentEventVisitor::visit_telescope_run(run_config, event_lifetime_manager);
    if(nsamp_!=old_nsamp) {
      auto* host_info = calin::provenance::system_info::the_host_info();
      const unsigned nv_samp = (nsamp_+31)/32;
      const unsigned nv_block = nv_samp*32;
      calin::util::memory::safe_aligned_recalloc(samples_, nv_block, host_info->log2_simd_vec_size());
    }
    return good;
  }

  bool visit_telescope_event(uint64_t seq_index,
    calin::ix::iact_data::telescope_event::TelescopeEvent* event) override
  {
    const calin::ix::iact_data::telescope_event::Waveforms* wf = nullptr;
    has_event_ = false;
    switch(gain_channel_to_treat_) {
    case HIGH_GAIN:
      if(event->has_high_gain_image() and
          event->high_gain_image().has_camera_waveforms()) {
        wf = &event->high_gain_image().camera_waveforms();
      }
      break;
    case LOW_GAIN:
      if(event->has_low_gain_image() and
          event->low_gain_image().has_camera_waveforms()) {
        wf = &event->low_gain_image().camera_waveforms();
      }
      break;
    case SINGLE_OR_MIXED_GAIN:
      if(event->has_image() and
          event->image().has_camera_waveforms()) {
        wf = &event->image().camera_waveforms();
      }
      break;
    }
    if(wf == nullptr)return true;
    has_event_ = true;
    chan_signal_type_.assign(wf->channel_signal_type().begin(), wf->channel_signal_type().end());
    const uint16_t* data = reinterpret_cast<const uint16_t*>(
      wf->raw_samples_array().data());
    vcl_analyze_waveforms(data);
    return true;
  }

#ifndef SWIG
private:

  void vcl_analyze_waveforms(const uint16_t* __restrict__ data)
  {
    unsigned nchan_block = (nchan_+num_int16-1)/num_int16;
    unsigned nsamp_block = (nsamp_+num_int16-1)/num_int16;
    unsigned nchan_left = nchan_;

    for(unsigned ichan_block=0;ichan_block<nchan_block;++ichan_block)
    {
      unsigned ioffset_l = ichan_block*2*num_int32;
      unsigned ioffset_h = ichan_block*2*num_int32 + num_int32;

      unsigned nchan_proc = std::min(num_int16, nchan_left);
      nchan_left -= num_int16;

      for(unsigned isamp_block=0;isamp_block<nsamp_block;++isamp_block)
      {
        for(unsigned ichan_proc=0;ichan_proc<nchan_proc;++ichan_proc) {
          samples_[isamp_block*num_int16 + ichan_proc].load(data
            + (ichan_block*num_int16 + ichan_proc)*nsamp_
            + isamp_block*num_int16);
        }
        calin::util::vcl::transpose(samples_ + isamp_block*num_int16);
      }

      int16_vt samp_max = samples_[0];
      uint16_vt samp_max_index = 0;

      int32_vt q_l = calin::util::vcl::extend_16_to_32_low(samples_[0]);
      int32_vt q_h = calin::util::vcl::extend_16_to_32_high(samples_[0]);

      int32_vt qt_l = 0;
      int32_vt qt_h = 0;

      // Fill the integration window sum up from first samples
      unsigned isamp_win_end;
      for(isamp_win_end=1;isamp_win_end<window_n_;++isamp_win_end) {
        int16_vt sample = samples_[isamp_win_end];

        int16_bvt samp_bigger_than_max = sample > samp_max;
        samp_max = vcl::select(samp_bigger_than_max, sample, samp_max);
        samp_max_index = vcl::select(samp_bigger_than_max, uint16_t(isamp_win_end), samp_max_index);

        int32_vt sample32;

        sample32 = calin::util::vcl::extend_16_to_32_low(sample);
        q_l += sample32;
        qt_l += sample32 * isamp_win_end;

        sample32 = calin::util::vcl::extend_16_to_32_high(sample);
        q_h += sample32;
        qt_h += sample32 * isamp_win_end;
      }

      if(bkg_window_0_ == 0) {
        q_l.store(chan_bkg_win_sum_ + ichan_block*2*num_int32);
        q_h.store(chan_bkg_win_sum_ + ichan_block*2*num_int32 + num_int32);
      }

      int32_vt all_q_l = q_l;
      int32_vt all_q_h = q_h;

      int32_vt sig_q_l = q_l;
      int32_vt sig_q_h = q_h;

      int32_vt opt_q_l = q_l;
      int32_vt opt_q_h = q_h;

      int32_vt opt_qt_l = qt_l;
      int32_vt opt_qt_h = qt_h;

      uint32_vt opt_index_l = 0;
      uint32_vt opt_index_h = 0;

      uint32_vt sig_win_index_l;
      uint32_vt sig_win_index_h;
      sig_win_index_l.load(sig_window_0_ + ioffset_l);
      sig_win_index_h.load(sig_window_0_ + ioffset_h);

      // Now move integration window along waveform
      unsigned isamp_win_begin;
      for(isamp_win_begin=0; isamp_win_end<nsamp_; ++isamp_win_end) {
        int16_vt sample = samples_[isamp_win_end];

        int16_bvt samp_bigger_than_max = sample > samp_max;
        samp_max = vcl::select(samp_bigger_than_max, sample, samp_max);
        samp_max_index = vcl::select(samp_bigger_than_max, uint16_t(isamp_win_end), samp_max_index);

        int32_vt sample32;
        sample32 = calin::util::vcl::extend_16_to_32_low(sample);
        q_l += sample32;
        qt_l += sample32 * isamp_win_end;
        all_q_l += sample32;

        sample32 = calin::util::vcl::extend_16_to_32_high(sample);
        q_h += sample32;
        qt_h += sample32 * isamp_win_end;
        all_q_h += sample32;

        sample = samples_[isamp_win_begin];

        sample32 = calin::util::vcl::extend_16_to_32_low(sample);
        q_l -= sample32;
        qt_l -= sample32 * isamp_win_begin;

        sample32 = calin::util::vcl::extend_16_to_32_high(sample);
        q_h -= sample32;
        qt_h -= sample32 * isamp_win_begin;

        ++isamp_win_begin;

        sig_q_l = vcl::select(sig_win_index_l == isamp_win_begin, q_l, sig_q_l);
        sig_q_h = vcl::select(sig_win_index_h == isamp_win_begin, q_h, sig_q_h);

        int32_bvt q_bigger_than_max;
        q_bigger_than_max = q_l > opt_q_l;
        opt_q_l = vcl::select(q_bigger_than_max, q_l, opt_q_l);
        opt_qt_l = vcl::select(q_bigger_than_max, qt_l, opt_qt_l);
        opt_index_l = vcl::select(q_bigger_than_max, isamp_win_begin, opt_index_l);

        q_bigger_than_max = q_h > opt_q_h;
        opt_q_h = vcl::select(q_bigger_than_max, q_h, opt_q_h);
        opt_qt_h = vcl::select(q_bigger_than_max, qt_h, opt_qt_h);
        opt_index_h = vcl::select(q_bigger_than_max, isamp_win_begin, opt_index_h);

        if(bkg_window_0_ == int(isamp_win_begin)) {
          q_l.store(chan_bkg_win_sum_ + ichan_block*2*num_int32);
          q_h.store(chan_bkg_win_sum_ + ichan_block*2*num_int32 + num_int32);
        }
      }

      calin::util::vcl::extend_16_to_32_low(samp_max).store(chan_max_ + ioffset_l);
      calin::util::vcl::extend_16_to_32_high(samp_max).store(chan_max_ + ioffset_h);
      calin::util::vcl::extend_16_to_32_low(samp_max_index).store(chan_max_index_ + ioffset_l);
      calin::util::vcl::extend_16_to_32_high(samp_max_index).store(chan_max_index_ + ioffset_h);
      sig_q_l.store(chan_sig_win_sum_ + ioffset_l);
      sig_q_h.store(chan_sig_win_sum_ + ioffset_h);
      opt_q_l.store(chan_opt_win_sum_ + ioffset_l);
      opt_q_h.store(chan_opt_win_sum_ + ioffset_h);
      opt_qt_l.store(chan_opt_win_sum_qt_ + ioffset_l);
      opt_qt_h.store(chan_opt_win_sum_qt_ + ioffset_h);
      opt_index_l.store(chan_opt_win_index_ + ioffset_l);
      opt_index_h.store(chan_opt_win_index_ + ioffset_h);
      all_q_l.store(chan_all_sum_ + ioffset_l);
      all_q_h.store(chan_all_sum_ + ioffset_h);
    }
  }

  int16_vt*__restrict__ samples_ = nullptr;
#endif
};






// *****************************************************************************
// *****************************************************************************
//
// POSSIBLY OBSOLETE VERSIONS
//
// *****************************************************************************
// *****************************************************************************

class SingleGainDualWindowWaveformTreatmentEventVisitor:
  public calin::iact_data::event_visitor::ParallelEventVisitor
{
public:
  SingleGainDualWindowWaveformTreatmentEventVisitor(
    calin::ix::iact_data::waveform_treatment_event_visitor::
      SingleGainDualWindowWaveformTreatmentEventVisitorConfig config = default_config(),
    bool treat_high_gain = true);
  virtual ~SingleGainDualWindowWaveformTreatmentEventVisitor();

  SingleGainDualWindowWaveformTreatmentEventVisitor* new_sub_visitor(
    std::map<calin::iact_data::event_visitor::ParallelEventVisitor*,
        calin::iact_data::event_visitor::ParallelEventVisitor*>
      antecedent_visitors = { }) override;

  bool visit_telescope_run(
    const calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration* run_config,
    calin::iact_data::event_visitor::EventLifetimeManager* event_lifetime_manager) override;

  bool visit_telescope_event(uint64_t seq_index,
    calin::ix::iact_data::telescope_event::TelescopeEvent* event) override;

  static calin::ix::iact_data::waveform_treatment_event_visitor::
    SingleGainDualWindowWaveformTreatmentEventVisitorConfig default_config()
  {
    calin::ix::iact_data::waveform_treatment_event_visitor::
      SingleGainDualWindowWaveformTreatmentEventVisitorConfig cfg;
    cfg.set_integration_n(16);
    return cfg;
  }

#ifndef SWIG
  // Exposed as public to facilitate testing
  void reconfigure(unsigned nchan, unsigned nsamp);
  void scalar_analyze_waveforms(const uint16_t* __restrict__ data);
#endif

  std::vector<float> chan_ped() const { return make_vec(chan_ped_est_); };

  std::vector<int> chan_max_index() const { return make_vec(chan_max_index_); }
  std::vector<int> chan_max() const { return make_vec(chan_max_); }
  std::vector<int> chan_bkg_win_sum() const { return make_vec(chan_bkg_win_sum_); }
  std::vector<int> chan_sig_win_sum() const { return make_vec(chan_sig_win_sum_); }
  std::vector<int> chan_sig_max_sum() const { return make_vec(chan_sig_max_sum_); }
  std::vector<int> chan_sig_max_sum_index() const { return make_vec(chan_sig_max_sum_index_); }
  std::vector<int> chan_all_sum_q() const { return make_vec(chan_all_sum_q_); }
  std::vector<int> chan_all_sum_qt() const { return make_vec(chan_all_sum_qt_); }

  std::vector<float> chan_sig() const { return make_vec(chan_sig_); }
  std::vector<float> chan_mean_t() const { return make_vec(chan_mean_t_); }

protected:
#ifndef SWIG

  template<typename T> std::vector<T> make_vec(const T* ptr) const {
    return std::vector<T>(ptr, ptr+nchan_);
  }

  calin::ix::iact_data::waveform_treatment_event_visitor::
    SingleGainDualWindowWaveformTreatmentEventVisitorConfig config_;
  bool treat_high_gain_ = true;
  unsigned nchan_ = 0;
  unsigned nsamp_ = 0;
  unsigned window_n_;
  int bkg_window_0_;
  int*__restrict__ sig_window_0_ = nullptr;

  float*__restrict__ chan_ped_est_ = nullptr;
  float ped_iir_old_;
  float ped_iir_new_;

  int*__restrict__ chan_max_ = nullptr;
  int*__restrict__ chan_max_index_ = nullptr;
  int*__restrict__ chan_bkg_win_sum_ = nullptr;
  int*__restrict__ chan_sig_win_sum_ = nullptr;
  int*__restrict__ chan_sig_max_sum_ = nullptr;
  int*__restrict__ chan_sig_max_sum_index_ = nullptr;
  int*__restrict__ chan_all_sum_q_ = nullptr;
  int*__restrict__ chan_all_sum_qt_ = nullptr;

  float*__restrict__ chan_sig_ = nullptr;
  float*__restrict__ chan_mean_t_ = nullptr;
#endif
};



template<typename VCLArchitecture>
class VCL_SingleGainDualWindowWaveformTreatmentEventVisitor:
  public SingleGainDualWindowWaveformTreatmentEventVisitor
{
public:
#ifndef SWIG
  CALIN_TYPEALIAS(int16_vt, typename VCLArchitecture::int16_vt);
  CALIN_TYPEALIAS(int16_bvt, typename VCLArchitecture::int16_bvt);
  CALIN_TYPEALIAS(uint16_vt, typename VCLArchitecture::uint16_vt);

  CALIN_TYPEALIAS(int32_vt, typename VCLArchitecture::int32_vt);
  CALIN_TYPEALIAS(int32_bvt, typename VCLArchitecture::int32_bvt);
  CALIN_TYPEALIAS(uint32_vt, typename VCLArchitecture::uint32_vt);

  CALIN_TYPEALIAS(float_vt, typename VCLArchitecture::float_vt);

  constexpr static unsigned num_int16 = VCLArchitecture::num_int16;
  constexpr static unsigned num_int32 = VCLArchitecture::num_int32;
  constexpr static unsigned num_float = VCLArchitecture::num_float;
#endif

  VCL_SingleGainDualWindowWaveformTreatmentEventVisitor(
      calin::ix::iact_data::waveform_treatment_event_visitor::
        SingleGainDualWindowWaveformTreatmentEventVisitorConfig config = default_config(),
      bool treat_high_gain = true):
    SingleGainDualWindowWaveformTreatmentEventVisitor(config, treat_high_gain)
  {
    /* nothing to see here */
  }

  virtual ~VCL_SingleGainDualWindowWaveformTreatmentEventVisitor()
  {
    free(samples_);
  }

  VCL_SingleGainDualWindowWaveformTreatmentEventVisitor<VCLArchitecture>* new_sub_visitor(
    std::map<calin::iact_data::event_visitor::ParallelEventVisitor*,
        calin::iact_data::event_visitor::ParallelEventVisitor*>
      antecedent_visitors = { }) override
  {
    return new VCL_SingleGainDualWindowWaveformTreatmentEventVisitor<VCLArchitecture>(config_);
  }

  bool visit_telescope_run(
    const calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration* run_config,
    calin::iact_data::event_visitor::EventLifetimeManager* event_lifetime_manager) override
  {
    bool old_nsamp = nsamp_;
    bool good = SingleGainDualWindowWaveformTreatmentEventVisitor::visit_telescope_run(run_config, event_lifetime_manager);
    if(nsamp_!=old_nsamp) {
      auto* host_info = calin::provenance::system_info::the_host_info();
      const unsigned nv_samp = (nsamp_+15)/16;
      const unsigned nv_block = nv_samp*16;
      calin::util::memory::safe_aligned_recalloc(samples_, nv_block, host_info->log2_simd_vec_size());
    }
    return good;
  }

  bool visit_telescope_event(uint64_t seq_index,
    calin::ix::iact_data::telescope_event::TelescopeEvent* event) override
  {
    const calin::ix::iact_data::telescope_event::Waveforms* wf = nullptr;
    if(treat_high_gain_) {
      if(event->has_high_gain_image() and
          event->high_gain_image().has_camera_waveforms()) {
        wf = &event->high_gain_image().camera_waveforms();
      }
    } else if(event->has_low_gain_image() and
        event->low_gain_image().has_camera_waveforms()) {
      wf = &event->low_gain_image().camera_waveforms();
    }
    if(wf == nullptr)return true;
    const uint16_t* data = reinterpret_cast<const uint16_t*>(
      wf->raw_samples_array().data());
    vcl_analyze_waveforms(data);
    return true;
  }

#ifndef SWIG
  void vcl_analyze_waveforms(const uint16_t* __restrict__ data)
  {
    unsigned nchan_block = (nchan_+num_int16-1)/num_int16;
    unsigned nsamp_block = (nsamp_+num_int16-1)/num_int16;
    unsigned nchan_left = nchan_;

    for(unsigned ichan_block=0;ichan_block<nchan_block;++ichan_block)
    {
      unsigned ioffset_l = ichan_block*2*num_int32;
      unsigned ioffset_h = ichan_block*2*num_int32 + num_int32;

      unsigned nchan_proc = std::min(num_int16, nchan_left);
      nchan_left -= num_int16;

      for(unsigned isamp_block=0;isamp_block<nsamp_block;++isamp_block)
      {
        for(unsigned ichan_proc=0;ichan_proc<nchan_proc;++ichan_proc) {
          samples_[isamp_block*num_int16 + ichan_proc].load(data
            + (ichan_block*num_int16 + ichan_proc)*nsamp_
            + isamp_block*num_int16);
        }
        calin::util::vcl::transpose(samples_ + isamp_block*num_int16);
      }

      int16_vt samp_max = samples_[0];
      uint16_vt samp_max_index = 0;

      int32_vt q_l = calin::util::vcl::extend_16_to_32_low(samples_[0]);
      int32_vt q_h = calin::util::vcl::extend_16_to_32_high(samples_[0]);

      int32_vt all_qt_l = 0;
      int32_vt all_qt_h = 0;

      // Fill the integration window sum up from first samples
      unsigned isamp_win_end;
      for(isamp_win_end=1;isamp_win_end<window_n_;++isamp_win_end) {
        int16_vt sample = samples_[isamp_win_end];

        int16_bvt samp_bigger_than_max = sample > samp_max;
        samp_max = vcl::select(samp_bigger_than_max, sample, samp_max);
        samp_max_index = vcl::select(samp_bigger_than_max, uint16_t(isamp_win_end), samp_max_index);

        int32_vt sample32;

        sample32 = calin::util::vcl::extend_16_to_32_low(sample);
        q_l += sample32;
        all_qt_l += sample32 * isamp_win_end;

        sample32 = calin::util::vcl::extend_16_to_32_high(sample);
        q_h += sample32;
        all_qt_h += sample32 * isamp_win_end;
      }

      if(bkg_window_0_ == 0) {
        q_l.store(chan_bkg_win_sum_ + ichan_block*2*num_int32);
        q_h.store(chan_bkg_win_sum_ + ichan_block*2*num_int32 + num_int32);
      }

      int32_vt all_q_l = q_l;
      int32_vt all_q_h = q_h;

      int32_vt sig_q_l = q_l;
      int32_vt sig_q_h = q_h;

      int32_vt max_q_l = q_l;
      int32_vt max_q_h = q_h;

      uint32_vt max_q_l_index = 0;
      uint32_vt max_q_h_index = 0;

      uint32_vt sig_win_index_l;
      uint32_vt sig_win_index_h;
      sig_win_index_l.load(sig_window_0_ + ioffset_l);
      sig_win_index_h.load(sig_window_0_ + ioffset_h);

      // Now move integration window along waveform
      unsigned isamp_win_begin;
      for(isamp_win_begin=0; isamp_win_end<nsamp_; ++isamp_win_end) {
        int16_vt sample = samples_[isamp_win_end];

        int16_bvt samp_bigger_than_max = sample > samp_max;
        samp_max = vcl::select(samp_bigger_than_max, sample, samp_max);
        samp_max_index = vcl::select(samp_bigger_than_max, uint16_t(isamp_win_end), samp_max_index);

        int32_vt sample32;
        sample32 = calin::util::vcl::extend_16_to_32_low(sample);
        q_l += sample32;
        all_q_l += sample32;
        all_qt_l += sample32 * isamp_win_end;

        sample32 = calin::util::vcl::extend_16_to_32_high(sample);
        q_h += sample32;
        all_q_h += sample32;
        all_qt_h += sample32 * isamp_win_end;

        sample = samples_[isamp_win_begin];

        q_l -= calin::util::vcl::extend_16_to_32_low(sample);
        q_h -= calin::util::vcl::extend_16_to_32_high(sample);

        ++isamp_win_begin;

        sig_q_l = vcl::select(sig_win_index_l == isamp_win_begin, q_l, sig_q_l);
        sig_q_h = vcl::select(sig_win_index_h == isamp_win_begin, q_h, sig_q_h);

        int32_bvt q_bigger_than_max;
        q_bigger_than_max = q_l > max_q_l;
        max_q_l = vcl::select(q_bigger_than_max, q_l, max_q_l);
        max_q_l_index = vcl::select(q_bigger_than_max, isamp_win_begin, max_q_l_index);

        q_bigger_than_max = q_h > max_q_h;
        max_q_h = vcl::select(q_bigger_than_max, q_h, max_q_h);
        max_q_h_index = vcl::select(q_bigger_than_max, isamp_win_begin, max_q_h_index);

        if(bkg_window_0_ == int(isamp_win_begin)) {
          q_l.store(chan_bkg_win_sum_ + ichan_block*2*num_int32);
          q_h.store(chan_bkg_win_sum_ + ichan_block*2*num_int32 + num_int32);
        }
      }

      calin::util::vcl::extend_16_to_32_low(samp_max).store(chan_max_ + ioffset_l);
      calin::util::vcl::extend_16_to_32_high(samp_max).store(chan_max_ + ioffset_h);
      calin::util::vcl::extend_16_to_32_low(samp_max_index).store(chan_max_index_ + ioffset_l);
      calin::util::vcl::extend_16_to_32_high(samp_max_index).store(chan_max_index_ + ioffset_h);
      sig_q_l.store(chan_sig_win_sum_ + ioffset_l);
      sig_q_h.store(chan_sig_win_sum_ + ioffset_h);
      max_q_l.store(chan_sig_max_sum_ + ioffset_l);
      max_q_h.store(chan_sig_max_sum_ + ioffset_h);
      max_q_l_index.store(chan_sig_max_sum_index_ + ioffset_l);
      max_q_h_index.store(chan_sig_max_sum_index_ + ioffset_h);
      all_q_l.store(chan_all_sum_q_ + ioffset_l);
      all_q_h.store(chan_all_sum_q_ + ioffset_h);
      all_qt_l.store(chan_all_sum_qt_ + ioffset_l);
      all_qt_h.store(chan_all_sum_qt_ + ioffset_h);

      float_vt ped;
      int32_vt bkg_i32;
      float_vt bkg;
      float_vt sig;
      float_vt mean_t;

      ped.load(chan_ped_est_ + ioffset_l);
      bkg_i32.load(chan_bkg_win_sum_ + ioffset_l);
      bkg = vcl::to_float(bkg_i32);
      ped = vcl::select(ped<0, bkg, ped_iir_old_*ped + ped_iir_new_*bkg);
      sig = vcl::to_float(sig_q_l) - ped;
      mean_t = (vcl::to_float(window_n_*all_qt_l) - ped*float(nsamp_*(nsamp_-1)/2))/
        (vcl::to_float(window_n_*all_q_l) - ped*float(nsamp_));

      ped.store(chan_ped_est_ + ioffset_l);
      sig.store(chan_sig_ + ioffset_l);
      mean_t.store(chan_mean_t_ + ioffset_l);

      ped.load(chan_ped_est_ + ioffset_h);
      bkg_i32.load(chan_bkg_win_sum_ + ioffset_h);
      bkg = vcl::to_float(bkg_i32);
      ped = vcl::select(ped<0, bkg, ped_iir_old_*ped + ped_iir_new_*bkg);
      sig = vcl::to_float(sig_q_h) - ped;
      mean_t = (vcl::to_float(window_n_*all_qt_h) - ped*float(nsamp_*(nsamp_-1)/2))/
        (vcl::to_float(window_n_*all_q_h) - ped*float(nsamp_));

      ped.store(chan_ped_est_ + ioffset_h);
      sig.store(chan_sig_ + ioffset_h);
      mean_t.store(chan_mean_t_ + ioffset_h);
    }
  }

private:
  int16_vt*__restrict__ samples_ = nullptr;
#endif
};



} } } // namespace calin::iact_data::waveform_treatment_event_visitor
