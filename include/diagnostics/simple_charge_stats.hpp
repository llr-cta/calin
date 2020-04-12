/*

   calin/diagnostics/channel_info.hpp -- Stephen Fegan -- 2020-03-20

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
  void integrate_one_gain_partials(
    calin::ix::diagnostics::simple_charge_stats::OneGainSimpleChargeStats* results_g,
    const calin::ix::diagnostics::simple_charge_stats::PartialOneGainChannelSimpleChargeStats& partials_gc);

  SimpleChargeStatsParallelEventVisitor* parent_ = nullptr;
  calin::ix::diagnostics::simple_charge_stats::SimpleChargeStatsConfig config_;
  calin::iact_data::waveform_treatment_event_visitor::OptimalWindowSumWaveformTreatmentParallelEventVisitor* high_gain_visitor_ = nullptr;
  calin::iact_data::waveform_treatment_event_visitor::OptimalWindowSumWaveformTreatmentParallelEventVisitor* low_gain_visitor_ = nullptr;

  bool has_dual_gain_ = false;
  calin::ix::diagnostics::simple_charge_stats::SimpleChargeStats results_;
  calin::ix::diagnostics::simple_charge_stats::PartialSimpleChargeStats partials_;

  std::vector<calin::math::histogram::Histogram1D*> ped_hist_;
};

} } } // namespace calin::diagnostics::simple_charge_stats
