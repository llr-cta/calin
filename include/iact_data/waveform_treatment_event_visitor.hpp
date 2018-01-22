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
#include <iact_data/waveform_treatment_event_visitor.pb.h>

namespace calin { namespace iact_data { namespace waveform_treatment_event_visitor {

class SingleGainDualWindowWaveformTreatmentEventVisitor:
  public calin::iact_data::event_visitor::TelescopeEventVisitor
{
public:
  SingleGainDualWindowWaveformTreatmentEventVisitor(
    calin::ix::iact_data::waveform_treatment_event_visitor::
      SingleGainDualWindowWaveformTreatmentEventVisitorConfig config = default_config(),
    bool treat_high_gain = true);
  virtual ~SingleGainDualWindowWaveformTreatmentEventVisitor();

  bool demand_waveforms() override;
  bool is_parallelizable() override;
  SingleGainDualWindowWaveformTreatmentEventVisitor* new_sub_visitor(
    const std::map<calin::iact_data::event_visitor::TelescopeEventVisitor*,
        calin::iact_data::event_visitor::TelescopeEventVisitor*>&
      antecedent_visitors = { }) override;

  bool visit_telescope_run(
    const calin::ix::iact_data::telescope_run_configuration::
      TelescopeRunConfiguration* run_config) override;

  bool visit_telescope_event(uint64_t seq_index,
    calin::ix::iact_data::telescope_event::TelescopeEvent* event) override;

#ifndef SWIG
  // This function allows testing of the code
  inline static void analyze_waveforms(
    const uint16_t* __restrict__ data, unsigned nchan, int nsamp,
    int window_n, int bkg_window_0, const int* sig_window_0,
    float*__restrict__ ped, float ped_iir_old, float ped_iir_new,
    int*__restrict__ chan_max_index, int*__restrict__ chan_max,
    int*__restrict__ chan_bkg_win_sum, int* chan_sig_win_sum,
    int*__restrict__ chan_sig_max_sum, int*__restrict__ chan_sig_max_sum_index,
    int*__restrict__ chan_all_sum_q, int*__restrict__ chan_all_sum_qt,
    float*__restrict__ chan_sig, float*__restrict__ chan_mean_t);
#endif

  static calin::ix::iact_data::waveform_treatment_event_visitor::
    SingleGainDualWindowWaveformTreatmentEventVisitorConfig default_config()
  {
    calin::ix::iact_data::waveform_treatment_event_visitor::
      SingleGainDualWindowWaveformTreatmentEventVisitorConfig cfg;
    cfg.set_integration_n(16);
    return cfg;
  }

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
  int* sig_window_0_ = nullptr;

  float* chan_ped_est_ = nullptr;
  float ped_iir_old_;
  float ped_iir_new_;

  int* chan_max_ = nullptr;
  int* chan_max_index_ = nullptr;
  int* chan_bkg_win_sum_ = nullptr;
  int* chan_sig_win_sum_ = nullptr;
  int* chan_sig_max_sum_ = nullptr;
  int* chan_sig_max_sum_index_ = nullptr;
  int* chan_all_sum_q_ = nullptr;
  int* chan_all_sum_qt_ = nullptr;

  float* chan_sig_ = nullptr;
  float* chan_mean_t_ = nullptr;
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
    const std::map<calin::iact_data::event_visitor::TelescopeEventVisitor*,
        calin::iact_data::event_visitor::TelescopeEventVisitor*>&
      antecedent_visitors = { }) override;

  bool visit_telescope_run(
    const calin::ix::iact_data::telescope_run_configuration::
      TelescopeRunConfiguration* run_config) override;

  bool visit_telescope_event(uint64_t seq_index,
    calin::ix::iact_data::telescope_event::TelescopeEvent* event) override;

#ifndef SWIG
  // This function allows testing of the code
  inline static void avx2_analyze_waveforms(
    const uint16_t* __restrict__ data, unsigned nchan, int nsamp,
    int window_n, int bkg_window_0, const int* sig_window_0,
    float*__restrict__ ped, float ped_iir_old, float ped_iir_new,
    int*__restrict__ chan_max_index, int*__restrict__ chan_max,
    int*__restrict__ chan_bkg_win_sum, int* chan_sig_win_sum,
    int*__restrict__ chan_sig_max_sum, int*__restrict__ chan_sig_max_sum_index,
    int*__restrict__ chan_all_sum_q, int*__restrict__ chan_all_sum_qt,
    float*__restrict__ chan_sig, float*__restrict__ chan_mean_t);
#endif
};

} } } // namespace calin::iact_data::waveform_treatment_event_visitor
