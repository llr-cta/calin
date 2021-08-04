/*

   calin/diagnostics/simple_charge_stats.hpp -- Stephen Fegan -- 2020-03-20

   Various information computed from channel data

   Copyright 2020, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <math/histogram.hpp>
#include <iact_data/event_visitor.hpp>
#include <iact_data/waveform_treatment_event_visitor.hpp>
#include <diagnostics/simple_charge_stats.pb.h>

namespace calin { namespace diagnostics { namespace simple_charge_stats {

class SimpleChargeStatsParallelEventVisitor:
  public calin::iact_data::event_visitor::ParallelEventVisitor
{
public:
  SimpleChargeStatsParallelEventVisitor(
    calin::iact_data::waveform_treatment_event_visitor::OptimalWindowSumWaveformTreatmentParallelEventVisitor* high_gain_visitor,
    calin::iact_data::waveform_treatment_event_visitor::OptimalWindowSumWaveformTreatmentParallelEventVisitor* low_gain_visitor,
    const calin::ix::diagnostics::simple_charge_stats::SimpleChargeStatsConfig& config = default_config());

  SimpleChargeStatsParallelEventVisitor(
      calin::iact_data::waveform_treatment_event_visitor::OptimalWindowSumWaveformTreatmentParallelEventVisitor* mixed_or_unique_gain_visitor,
      const calin::ix::diagnostics::simple_charge_stats::SimpleChargeStatsConfig& config = default_config()):
    SimpleChargeStatsParallelEventVisitor(mixed_or_unique_gain_visitor, nullptr, config)
  { /* nothing to see here */ }

  virtual ~SimpleChargeStatsParallelEventVisitor();

  SimpleChargeStatsParallelEventVisitor* new_sub_visitor(
    std::map<calin::iact_data::event_visitor::ParallelEventVisitor*,
        calin::iact_data::event_visitor::ParallelEventVisitor*>
      antecedent_visitors = { }) override;

  bool visit_telescope_run(
    const calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration* run_config,
    calin::iact_data::event_visitor::EventLifetimeManager* event_lifetime_manager) override;

  bool leave_telescope_run() override;

  bool visit_telescope_event(uint64_t seq_index,
    calin::ix::iact_data::telescope_event::TelescopeEvent* event) override;

  bool merge_results() override;

#ifndef SWIG
  calin::ix::diagnostics::simple_charge_stats::SimpleChargeStats* simple_charge_stats(
    calin::ix::diagnostics::simple_charge_stats::SimpleChargeStats* stats = nullptr) const;
#else
  calin::ix::diagnostics::simple_charge_stats::SimpleChargeStats* simple_charge_stats() const;
  void simple_charge_stats(calin::ix::diagnostics::simple_charge_stats::SimpleChargeStats* stats) const;
#endif

  static calin::ix::diagnostics::simple_charge_stats::SimpleChargeStatsConfig default_config();

  // const calin::ix::diagnostics::simple_charge_stats::PartialSimpleChargeStats& partials() const { return partials_; }

private:
  struct SingleGainChannelHists {
    SingleGainChannelHists(double time_resolution, double time_max = 3600.0):
      all_pedwin_1_sum_vs_time(new calin::math::histogram::Histogram1D(time_resolution, -60.0, time_max, 0.0)),
      all_pedwin_q_sum_vs_time(new calin::math::histogram::Histogram1D(time_resolution, -60.0, time_max, 0.0)),
      all_pedwin_q2_sum_vs_time(new calin::math::histogram::Histogram1D(time_resolution, -60.0, time_max, 0.0)),
      ped_wf_1_sum_vs_time(new calin::math::histogram::Histogram1D(time_resolution, -60.0, time_max, 0.0)),
      ped_wf_q_sum_vs_time(new calin::math::histogram::Histogram1D(time_resolution, -60.0, time_max, 0.0)),
      ped_wf_q2_sum_vs_time(new calin::math::histogram::Histogram1D(time_resolution, -60.0, time_max, 0.0))
    { /* nothing to see here */ }

    ~SingleGainChannelHists() {
      delete all_pedwin_1_sum_vs_time;
      delete all_pedwin_q_sum_vs_time;
      delete all_pedwin_q2_sum_vs_time;
      delete ped_wf_1_sum_vs_time;
      delete ped_wf_q_sum_vs_time;
      delete ped_wf_q2_sum_vs_time;
    }

    calin::math::histogram::Histogram1D* all_pedwin_1_sum_vs_time;
    calin::math::histogram::Histogram1D* all_pedwin_q_sum_vs_time;
    calin::math::histogram::Histogram1D* all_pedwin_q2_sum_vs_time;

    calin::math::histogram::Histogram1D* ped_wf_1_sum_vs_time;
    calin::math::histogram::Histogram1D* ped_wf_q_sum_vs_time;
    calin::math::histogram::Histogram1D* ped_wf_q2_sum_vs_time;
  };

  struct ChannelHists {
    ChannelHists(bool has_dual_gain_, double time_resolution, double time_max = 3600.0):
      high_gain(new SingleGainChannelHists(time_resolution, time_max)),
      low_gain(has_dual_gain_ ? new SingleGainChannelHists(time_resolution, time_max) : nullptr)
    { /* nothing to see here */ }

    ~ChannelHists() {
      delete high_gain;
      delete low_gain;
    }

    SingleGainChannelHists* high_gain = nullptr;
    SingleGainChannelHists* low_gain = nullptr;
  };

  struct CameraHists {
    CameraHists(bool has_dual_gain_, double time_resolution, double time_max = 3600.0):
      high_gain(new SingleGainChannelHists(time_resolution, time_max)),
      low_gain(has_dual_gain_ ? new SingleGainChannelHists(time_resolution, time_max) : nullptr),
      num_channel_triggered_hist(new calin::math::histogram::Histogram1D(1.0)),
      num_contiguous_channel_triggered_hist(new calin::math::histogram::Histogram1D(1.0)),
      phys_trig_num_channel_triggered_hist(new calin::math::histogram::Histogram1D(1.0)),
      phys_trig_num_contiguous_channel_triggered_hist(new calin::math::histogram::Histogram1D(1.0))
    { /* nothing to see here */ }

    ~CameraHists() {
      delete num_channel_triggered_hist;
      delete num_contiguous_channel_triggered_hist;
      delete phys_trig_num_channel_triggered_hist;
      delete phys_trig_num_contiguous_channel_triggered_hist;
      delete high_gain;
      delete low_gain;
    }

    SingleGainChannelHists* high_gain = nullptr;
    SingleGainChannelHists* low_gain = nullptr;

    calin::math::histogram::Histogram1D* num_channel_triggered_hist = nullptr;
    calin::math::histogram::Histogram1D* num_contiguous_channel_triggered_hist = nullptr;
    calin::math::histogram::Histogram1D* phys_trig_num_channel_triggered_hist = nullptr;
    calin::math::histogram::Histogram1D* phys_trig_num_contiguous_channel_triggered_hist = nullptr;
  };

  void integrate_one_gain_partials(
    calin::ix::diagnostics::simple_charge_stats::OneGainSimpleChargeStats* results_g,
    const calin::ix::diagnostics::simple_charge_stats::PartialOneGainChannelSimpleChargeStats& partials_gc);

  void integrate_one_gain_camera_partials(
    calin::ix::diagnostics::simple_charge_stats::OneGainSimpleChargeStats* results_g,
    const calin::ix::diagnostics::simple_charge_stats::PartialOneGainCameraSimpleChargeStats& partials_gc);

  void dump_single_gain_channel_hists_to_partials(
    const SingleGainChannelHists& hists,
    calin::ix::diagnostics::simple_charge_stats::PartialOneGainChannelSimpleChargeStats* partials);

  void dump_single_gain_camera_hists_to_partials(
    const SingleGainChannelHists& hists,
    calin::ix::diagnostics::simple_charge_stats::PartialOneGainCameraSimpleChargeStats* partials);

  void record_one_gain_channel_data(const calin::ix::iact_data::telescope_event::TelescopeEvent* event,
    const calin::iact_data::waveform_treatment_event_visitor::OptimalWindowSumWaveformTreatmentParallelEventVisitor* sum_visitor,
    unsigned ichan, double elapsed_event_time,
    calin::ix::diagnostics::simple_charge_stats::PartialOneGainChannelSimpleChargeStats* one_gain_stats,
    SingleGainChannelHists* one_gain_hists,
    unsigned& nsum, int64_t& opt_sum, int64_t& sig_sum, int64_t& bkg_sum, int64_t& wf_sum,
    int wf_clipping_value);

  void record_one_visitor_data(uint64_t seq_index, const calin::ix::iact_data::telescope_event::TelescopeEvent* event,
    const calin::iact_data::waveform_treatment_event_visitor::OptimalWindowSumWaveformTreatmentParallelEventVisitor* sum_visitor,
    calin::ix::diagnostics::simple_charge_stats::PartialSimpleChargeStats* partials);

  SimpleChargeStatsParallelEventVisitor* parent_ = nullptr;
  calin::ix::diagnostics::simple_charge_stats::SimpleChargeStatsConfig config_;
  calin::iact_data::waveform_treatment_event_visitor::OptimalWindowSumWaveformTreatmentParallelEventVisitor* high_gain_visitor_ = nullptr;
  calin::iact_data::waveform_treatment_event_visitor::OptimalWindowSumWaveformTreatmentParallelEventVisitor* low_gain_visitor_ = nullptr;

  bool has_dual_gain_ = false;
  calin::ix::diagnostics::simple_charge_stats::SimpleChargeStats results_;
  calin::ix::diagnostics::simple_charge_stats::PartialSimpleChargeStats partials_;

  std::vector<ChannelHists*> chan_hists_;
  CameraHists* camera_hists_ = nullptr;
  calin::ix::iact_data::instrument_layout::CameraLayout* data_order_camera_ = nullptr;
  std::vector<int> channel_island_id_;
  std::vector<int> channel_island_count_;
};

} } } // namespace calin::diagnostics::simple_charge_stats
