/*

   calin/iact_data/waveform_treatment_event_visitor.hpp -- Stephen Fegan -- 2018-01-11

   Waveform treatment event data visitor - process waveforms

   Copyright 2018, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, Ecole Polytechnique, CNRS/IN2P3

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
#include <math/simd.hpp>
#include <util/vcl.hpp>
#include <iact_data/waveform_treatment_event_visitor.pb.h>

namespace calin { namespace iact_data { namespace waveform_treatment_event_visitor {

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
    const std::map<calin::iact_data::event_visitor::ParallelEventVisitor*,
        calin::iact_data::event_visitor::ParallelEventVisitor*>&
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

class AVX2_SingleGainDualWindowWaveformTreatmentEventVisitor:
  public SingleGainDualWindowWaveformTreatmentEventVisitor
{
public:
  AVX2_SingleGainDualWindowWaveformTreatmentEventVisitor(
      calin::ix::iact_data::waveform_treatment_event_visitor::
        SingleGainDualWindowWaveformTreatmentEventVisitorConfig config = default_config(),
      bool treat_high_gain = true):
    SingleGainDualWindowWaveformTreatmentEventVisitor(config, treat_high_gain)
  {
    /* nothing to see here */
  }

  virtual ~AVX2_SingleGainDualWindowWaveformTreatmentEventVisitor();

  AVX2_SingleGainDualWindowWaveformTreatmentEventVisitor* new_sub_visitor(
    const std::map<calin::iact_data::event_visitor::ParallelEventVisitor*,
        calin::iact_data::event_visitor::ParallelEventVisitor*>&
      antecedent_visitors = { }) override;

  bool visit_telescope_run(
    const calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration* run_config,
    calin::iact_data::event_visitor::EventLifetimeManager* event_lifetime_manager) override;

  bool visit_telescope_event(uint64_t seq_index,
    calin::ix::iact_data::telescope_event::TelescopeEvent* event) override;

#ifndef SWIG
  void avx2_analyze_waveforms(const uint16_t* __restrict__ data);
  void avx2_analyze_waveforms_v2(const uint16_t* __restrict__ data);
  void avx2_analyze_waveforms_v3(const uint16_t* __restrict__ data);
#endif

protected:
#if defined(__AVX2__) and defined(__FMA__)
#ifndef SWIG
  __m256i*__restrict__ samples_ = nullptr;
  __m256i*__restrict__ q_l_ = nullptr;
  __m256i*__restrict__ q_u_ = nullptr;
  __m256i*__restrict__ qt_l_ = nullptr;
  __m256i*__restrict__ qt_u_ = nullptr;
#endif
#endif
};

template<typename VCLArchitecture>
class VCL_SingleGainDualWindowWaveformTreatmentEventVisitor:
  public SingleGainDualWindowWaveformTreatmentEventVisitor,
  public VCLArchitecture
{
public:
#ifndef SWIG
  using typename VCLArchitecture::int16_vt;
  using typename VCLArchitecture::int16_bvt;
  using typename VCLArchitecture::uint16_vt;

  using typename VCLArchitecture::int32_vt;
  using typename VCLArchitecture::int32_bvt;
  using typename VCLArchitecture::uint32_vt;

  using typename VCLArchitecture::float_vt;

  using VCLArchitecture::num_int16;
  using VCLArchitecture::num_int32;
  using VCLArchitecture::num_float;
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
    // nothing to see here
  }

  VCL_SingleGainDualWindowWaveformTreatmentEventVisitor<VCLArchitecture>* new_sub_visitor(
    const std::map<calin::iact_data::event_visitor::ParallelEventVisitor*,
        calin::iact_data::event_visitor::ParallelEventVisitor*>&
      antecedent_visitors = { }) override
  {
    return new VCL_SingleGainDualWindowWaveformTreatmentEventVisitor<VCLArchitecture>(config_);
  }

  bool visit_telescope_run(
    const calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration* run_config,
    calin::iact_data::event_visitor::EventLifetimeManager* event_lifetime_manager) override
  {
    bool good = SingleGainDualWindowWaveformTreatmentEventVisitor::visit_telescope_run(run_config, event_lifetime_manager);
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
    int16_vt samples[num_int16];
    unsigned nchan_block = (nchan_+num_int16-1)/num_int16;
    unsigned nsamp_block = (nsamp_+num_int16-1)/num_int16;
    unsigned nchan_left = nchan_;

    for(unsigned ichan_block=0;ichan_block<nchan_block;++ichan_block)
    {
      unsigned nchan_proc = std::min(num_int16, nchan_left);
      nchan_left -= num_int16;

      unsigned isamp = 0;
      unsigned nsamp_left = nsamp_;

      int16_vt samp_max = -32768;
      uint16_vt isamp_max = 0;

      for(unsigned isamp_block=0;isamp_block<nsamp_block;++isamp_block)
      {
        for(unsigned ichan_proc=0;ichan_proc<nchan_proc;++ichan_proc) {
          samples[ichan_proc].load(data
            + (ichan_block*num_int16 + ichan_proc)*nsamp_
            + isamp_block*num_int16);
        }
        calin::util::vcl::transpose(samples);

        unsigned nsamp_proc = std::min(num_int16, nsamp_left);
        for(unsigned isamp_proc = 0; isamp_proc < nsamp_proc; ++isamp_proc, ++isamp)
        {
          int16_bvt samp_bigger_than_max = samples[isamp_proc] > samp_max;
          samp_max = select(samp_bigger_than_max, samples[isamp_proc], samp_max);
          isamp_max = select(samp_bigger_than_max, uint16_t(isamp), isamp_max);
        }
        nsamp_left -= num_int16;
      }

      vcl::extend_low(samp_max).store(chan_max_ + ichan_block*16);
      vcl::extend_high(samp_max).store(chan_max_ + ichan_block*16 + 8);
      vcl::extend_low(isamp_max).store(chan_max_index_ + ichan_block*16);
      vcl::extend_high(isamp_max).store(chan_max_index_ + ichan_block*16 + 8);
    }
  }
#endif
};



} } } // namespace calin::iact_data::waveform_treatment_event_visitor
