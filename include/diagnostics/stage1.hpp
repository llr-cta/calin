/*

   calin/diagnostics/stage1.hpp -- Stephen Fegan -- 2020-03-28

   Stage 1 analysis - calculate distributions and perform some low-level
                      diagnostics on events

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

#include <diagnostics/stage1.pb.h>
#include <iact_data/event_visitor.hpp>
#include <iact_data/waveform_treatment_event_visitor.hpp>
#include <diagnostics/run_info.hpp>
#include <diagnostics/simple_charge_stats.hpp>
#include <diagnostics/simple_charge_hists.hpp>
#include <diagnostics/clock_regression.hpp>
#include <diagnostics/waveform.hpp>
#include <iact_data/nectarcam_ancillary_data.hpp>

namespace calin { namespace diagnostics { namespace stage1 {

class Stage1ParallelEventVisitor:
  public calin::iact_data::event_visitor::FilteredDelegatingParallelEventVisitor
{
public:
  Stage1ParallelEventVisitor(const calin::ix::diagnostics::stage1::Stage1Config& config = default_config());

  virtual ~Stage1ParallelEventVisitor();

  bool visit_telescope_run(
    const calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration* run_config,
    calin::iact_data::event_visitor::EventLifetimeManager* event_lifetime_manager) override;

  bool leave_telescope_run() override;

#ifndef SWIG
  calin::ix::diagnostics::stage1::Stage1* stage1_results(
    calin::ix::diagnostics::stage1::Stage1* stage1 = nullptr) const;
#else
  calin::ix::diagnostics::stage1::Stage1* stage1_results() const;
  void stage1_results(calin::ix::diagnostics::stage1::Stage1* stage1) const;
#endif

  static calin::ix::diagnostics::stage1::Stage1Config default_config();

private:
  calin::ix::diagnostics::stage1::Stage1Config config_;
  const calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration* run_config_ = nullptr;

  calin::iact_data::waveform_treatment_event_visitor::OptimalWindowSumWaveformTreatmentParallelEventVisitor* hg_sum_pev_ = nullptr;
  calin::iact_data::waveform_treatment_event_visitor::OptimalWindowSumWaveformTreatmentParallelEventVisitor* lg_sum_pev_ = nullptr;

  calin::diagnostics::run_info::RunInfoDiagnosticsParallelEventVisitor* run_info_pev_ = nullptr;

  calin::diagnostics::simple_charge_stats::SimpleChargeStatsParallelEventVisitor* charge_stats_pev_ = nullptr;

  calin::diagnostics::simple_charge_hists::SimpleChargeHistsParallelEventVisitor* charge_hists_phy_pev_ = nullptr;
  calin::diagnostics::simple_charge_hists::SimpleChargeHistsParallelEventVisitor* charge_hists_ped_pev_ = nullptr;
  calin::diagnostics::simple_charge_hists::SimpleChargeHistsParallelEventVisitor* charge_hists_ext_pev_ = nullptr;
  calin::diagnostics::simple_charge_hists::SimpleChargeHistsParallelEventVisitor* charge_hists_int_pev_ = nullptr;

  calin::diagnostics::simple_charge_hists::SimpleChargeHistsParallelEventVisitor* charge_hists_trig_bit_set_pev_ = nullptr;
  calin::diagnostics::simple_charge_hists::SimpleChargeHistsParallelEventVisitor* charge_hists_trig_bit_clr_pev_ = nullptr;

  calin::diagnostics::waveform::WaveformSumParallelEventVisitor* wf_mean_phy_pev_ = nullptr;
  calin::diagnostics::waveform::WaveformSumParallelEventVisitor* wf_mean_ped_pev_ = nullptr;
  calin::diagnostics::waveform::WaveformSumParallelEventVisitor* wf_mean_ext_pev_ = nullptr;
  calin::diagnostics::waveform::WaveformSumParallelEventVisitor* wf_mean_int_pev_ = nullptr;

  calin::ix::iact_data::nectarcam_ancillary_data::NectarCAMAncillaryData* nectarcam_ancillary_data_ = nullptr;

  calin::diagnostics::clock_regression::ClockRegressionParallelEventVisitor* clock_regression_pev_ = nullptr;
};

} } } // namespace calin::diagnostics::stage1
