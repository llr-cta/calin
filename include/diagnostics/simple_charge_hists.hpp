/*

   calin/diagnostics/simple_charge_hists.hpp -- Stephen Fegan -- 2020-04-17

   Various histograms computed from channel data

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
#include <diagnostics/simple_charge_hists.pb.h>

namespace calin { namespace diagnostics { namespace simple_charge_hists {

class SimpleChargeHistsParallelEventVisitor:
  public calin::iact_data::event_visitor::ParallelEventVisitor
{
public:
  SimpleChargeHistsParallelEventVisitor(
    calin::iact_data::waveform_treatment_event_visitor::OptimalWindowSumWaveformTreatmentParallelEventVisitor* high_gain_visitor,
    calin::iact_data::waveform_treatment_event_visitor::OptimalWindowSumWaveformTreatmentParallelEventVisitor* low_gain_visitor,
    const calin::ix::diagnostics::simple_charge_hists::SimpleChargeHistsConfig& config = ped_trig_default_config());

  SimpleChargeHistsParallelEventVisitor(
      calin::iact_data::waveform_treatment_event_visitor::OptimalWindowSumWaveformTreatmentParallelEventVisitor* mixed_or_unique_gain_visitor,
      const calin::ix::diagnostics::simple_charge_hists::SimpleChargeHistsConfig& config = ped_trig_default_config()):
    SimpleChargeHistsParallelEventVisitor(mixed_or_unique_gain_visitor, nullptr, config)
  { /* nothing to see here */ }

  virtual ~SimpleChargeHistsParallelEventVisitor();

  SimpleChargeHistsParallelEventVisitor* new_sub_visitor(
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
  calin::ix::diagnostics::simple_charge_hists::SimpleChargeHists* simple_charge_hists(
    calin::ix::diagnostics::simple_charge_hists::SimpleChargeHists* stats = nullptr) const;
#else
  calin::ix::diagnostics::simple_charge_hists::SimpleChargeHists* simple_charge_hists() const;
  void simple_charge_hists(calin::ix::diagnostics::simple_charge_hists::SimpleChargeHists* stats) const;
#endif

  static calin::ix::diagnostics::simple_charge_hists::SimpleChargeHistsConfig ped_trig_default_config();
  static calin::ix::diagnostics::simple_charge_hists::SimpleChargeHistsConfig phy_trig_default_config();
  static calin::ix::diagnostics::simple_charge_hists::SimpleChargeHistsConfig ext_trig_default_config();
  static calin::ix::diagnostics::simple_charge_hists::SimpleChargeHistsConfig int_trig_default_config();

  // const calin::ix::diagnostics::simple_charge_stats::PartialSimpleChargeStats& partials() const { return partials_; }

private:
#if 0
  struct SingleGainChannelHists {
    SingleGainChannelHists():
      full_wf_qsum(new_hist());
      full_wf_max();
      opt_win_qsum();
      opt_win_index();
      ped_win_qsum();
      sig_win_qsum();

      ped_wf_q_sum(new calin::math::histogram::Histogram1D(dark_resolution)),
      ext_opt_win_q_sum(new calin::math::histogram::Histogram1D(bright_resolution)),
      ext_opt_win_index(new calin::math::histogram::Histogram1D(1.0)),
      ext_opt_win_max(new calin::math::histogram::Histogram1D(bright_resolution)),
      phys_opt_win_q_sum(new calin::math::histogram::Histogram1D(bright_resolution)),
      phys_opt_win_index(new calin::math::histogram::Histogram1D(1.0))
    { /* nothing to see here */ }

    ~SingleGainChannelHists() {
      delete ped_wf_q_sum;
      delete ped_wf_1_sum_vs_time;
      delete ped_wf_q_sum_vs_time;
      delete ped_wf_q2_sum_vs_time;
      delete ext_opt_win_q_sum;
      delete ext_opt_win_index;
      delete ext_opt_win_max;
      delete phys_opt_win_q_sum;
      delete phys_opt_win_index;
    }

    calin::math::histogram::Histogram1D* ped_wf_q_sum;

    calin::math::histogram::Histogram1D* ped_wf_1_sum_vs_time;
    calin::math::histogram::Histogram1D* ped_wf_q_sum_vs_time;
    calin::math::histogram::Histogram1D* ped_wf_q2_sum_vs_time;

    calin::math::histogram::Histogram1D* ext_opt_win_q_sum;
    calin::math::histogram::Histogram1D* ext_opt_win_index;
    calin::math::histogram::Histogram1D* ext_opt_win_max;

    calin::math::histogram::Histogram1D* phys_opt_win_q_sum;
    calin::math::histogram::Histogram1D* phys_opt_win_index;
  };

  struct DualGainChannelHists {
    DualGainChannelHists(double lg_resolution):
      ext_opt_win_1_sum_vs_qlg(new calin::math::histogram::Histogram1D(lg_resolution)),
      ext_opt_win_qhg_sum_vs_qlg(new calin::math::histogram::Histogram1D(lg_resolution)),
      ext_opt_win_q2hg_sum_vs_qlg(new calin::math::histogram::Histogram1D(lg_resolution)),
      phy_opt_1_sum_vs_qlg(new calin::math::histogram::Histogram1D(lg_resolution)),
      phy_opt_qhg_sum_vs_qlg(new calin::math::histogram::Histogram1D(lg_resolution)),
      phy_opt_q2hg_sum_vs_qlg(new calin::math::histogram::Histogram1D(lg_resolution))
    { /* nothing to see here */ }

    ~DualGainChannelHists() {
      delete ext_opt_win_1_sum_vs_qlg;
      delete ext_opt_win_qhg_sum_vs_qlg;
      delete ext_opt_win_q2hg_sum_vs_qlg;
      delete phy_opt_1_sum_vs_qlg;
      delete phy_opt_qhg_sum_vs_qlg;
      delete phy_opt_q2hg_sum_vs_qlg;
    }

    calin::math::histogram::Histogram1D* ext_opt_win_1_sum_vs_qlg;
    calin::math::histogram::Histogram1D* ext_opt_win_qhg_sum_vs_qlg;
    calin::math::histogram::Histogram1D* ext_opt_win_q2hg_sum_vs_qlg;
    calin::math::histogram::Histogram1D* phy_opt_1_sum_vs_qlg;
    calin::math::histogram::Histogram1D* phy_opt_qhg_sum_vs_qlg;
    calin::math::histogram::Histogram1D* phy_opt_q2hg_sum_vs_qlg;
  };

  struct ChannelHists {
    ChannelHists(bool has_dual_gain_, double time_resolution,
        double hg_dark_resolution, double hg_bright_resolution,
        double lg_dark_resolution, double lg_bright_resolution):
      high_gain(new SingleGainChannelHists(time_resolution, hg_dark_resolution, hg_bright_resolution)),
      low_gain(has_dual_gain_ ? new SingleGainChannelHists(time_resolution, lg_dark_resolution, lg_bright_resolution) : nullptr),
      dual_gain(has_dual_gain_ ? new DualGainChannelHists(lg_bright_resolution) : nullptr)
    { /* nothing to see here */ }

    ~ChannelHists() {
      delete high_gain;
      delete low_gain;
      delete dual_gain;
    }

    SingleGainChannelHists* high_gain = nullptr;
    SingleGainChannelHists* low_gain = nullptr;
    DualGainChannelHists* dual_gain = nullptr;
  };

  void integrate_one_gain_partials(
    calin::ix::diagnostics::simple_charge_stats::OneGainSimpleChargeStats* results_g,
    const calin::ix::diagnostics::simple_charge_stats::PartialOneGainChannelSimpleChargeStats& partials_gc);

  void dump_single_gain_channel_hists_to_partials(
    const SingleGainChannelHists& hists,
    calin::ix::diagnostics::simple_charge_stats::PartialOneGainChannelSimpleChargeStats* partials);

  void record_one_gain_channel_data(const calin::ix::iact_data::telescope_event::TelescopeEvent* event,
    const calin::iact_data::waveform_treatment_event_visitor::OptimalWindowSumWaveformTreatmentParallelEventVisitor* sum_visitor,
    unsigned ichan, double elapsed_event_time,
    calin::ix::diagnostics::simple_charge_stats::PartialOneGainChannelSimpleChargeStats* one_gain_stats,
    SingleGainChannelHists* one_gain_hists);

  void record_one_visitor_data(uint64_t seq_index, const calin::ix::iact_data::telescope_event::TelescopeEvent* event,
    const calin::iact_data::waveform_treatment_event_visitor::OptimalWindowSumWaveformTreatmentParallelEventVisitor* sum_visitor,
    calin::ix::diagnostics::simple_charge_stats::PartialSimpleChargeStats* partials);
#endif

  SimpleChargeHistsParallelEventVisitor* parent_ = nullptr;
  calin::ix::diagnostics::simple_charge_hists::SimpleChargeHistsConfig config_;
  calin::iact_data::waveform_treatment_event_visitor::OptimalWindowSumWaveformTreatmentParallelEventVisitor* high_gain_visitor_ = nullptr;
  calin::iact_data::waveform_treatment_event_visitor::OptimalWindowSumWaveformTreatmentParallelEventVisitor* low_gain_visitor_ = nullptr;

  bool has_dual_gain_ = false;
  // calin::ix::diagnostics::simple_charge_stats::SimpleChargeStats results_;
  // calin::ix::diagnostics::simple_charge_stats::PartialSimpleChargeStats partials_;

  // std::vector<ChannelHists*> chan_hists_;
  //ChannelHists* camera_hists_;
};

} } } // namespace calin::diagnostics::simple_charge_stats