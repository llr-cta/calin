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

#include <util/log.hpp>
#include <diagnostics/stage1.hpp>
#include <provenance/anthology.hpp>

using namespace calin::util::log;
using namespace calin::diagnostics::stage1;

namespace {
  char LOGO[] =
    "  ____                ___       ____  \n"
    " / / /    _________ _/ (_)___   \\ \\ \\ \n" \
    "/ / /    / ___/ __ `/ / / __ \\   \\ \\ \\   __              \n" \
    "\\ \\ \\   / /__/ /_/ / / / / / /   / / /  (__|_ _. _  _  /|\n" \
    " \\_\\_\\  \\___/\\__,_/_/_/_/ /_/   /_/_/   __)|_(_|(_|(/_  |\n" \
    "                                                ._|      \n";
}

Stage1ParallelEventVisitor::Stage1ParallelEventVisitor():
  FilteredDelegatingParallelEventVisitor()
{
  hg_sum_pev_ = calin::iact_data::waveform_treatment_event_visitor::
    OptimalWindowSumWaveformTreatmentParallelEventVisitor::New();
  lg_sum_pev_ = calin::iact_data::waveform_treatment_event_visitor::
    OptimalWindowSumWaveformTreatmentParallelEventVisitor::New(
      calin::iact_data::waveform_treatment_event_visitor::
      OptimalWindowSumWaveformTreatmentParallelEventVisitor::LOW_GAIN);

  run_info_pev_ = new calin::diagnostics::run_info::RunInfoDiagnosticsParallelEventVisitor();
  charge_stats_pev_ = new calin::diagnostics::simple_charge_stats::
    SimpleChargeStatsParallelEventVisitor(hg_sum_pev_, lg_sum_pev_);

  wf_mean_phy_pev_ = new calin::diagnostics::waveform::WaveformSumParallelEventVisitor();
  wf_mean_ped_pev_ = new calin::diagnostics::waveform::WaveformSumParallelEventVisitor();
  wf_mean_ext_pev_ = new calin::diagnostics::waveform::WaveformSumParallelEventVisitor();
  wf_mean_int_pev_ = new calin::diagnostics::waveform::WaveformSumParallelEventVisitor();

  this->add_visitor(hg_sum_pev_);
  this->add_visitor(lg_sum_pev_);
  this->add_visitor(run_info_pev_);
  this->add_visitor(charge_stats_pev_);
  this->add_physics_trigger_visitor(wf_mean_phy_pev_);
  this->add_pedestal_trigger_visitor(wf_mean_ped_pev_);
  this->add_external_flasher_trigger_visitor(wf_mean_ext_pev_);
  this->add_internal_flasher_trigger_visitor(wf_mean_int_pev_);
}

Stage1ParallelEventVisitor::~Stage1ParallelEventVisitor()
{
  delete wf_mean_int_pev_;
  delete wf_mean_ext_pev_;
  delete wf_mean_ped_pev_;
  delete wf_mean_phy_pev_;
  delete charge_stats_pev_;
  delete run_info_pev_;
  delete lg_sum_pev_;
  delete hg_sum_pev_;
}

bool Stage1ParallelEventVisitor::visit_telescope_run(
  const calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration* run_config,
  calin::iact_data::event_visitor::EventLifetimeManager* event_lifetime_manager)
{
  LOG(INFO) << LOGO;
  return FilteredDelegatingParallelEventVisitor::visit_telescope_run(
    run_config, event_lifetime_manager);
}

bool Stage1ParallelEventVisitor::leave_telescope_run()
{
  bool good = FilteredDelegatingParallelEventVisitor::leave_telescope_run();
  return good;
}

calin::ix::diagnostics::stage1::Stage1* Stage1ParallelEventVisitor::stage1_results(
  calin::ix::diagnostics::stage1::Stage1* stage1) const
{
  if(stage1 == nullptr)stage1 = new calin::ix::diagnostics::stage1::Stage1;

  run_info_pev_->run_config(stage1->mutable_run_config());
  run_info_pev_->run_info(stage1->mutable_run_info());
  charge_stats_pev_->simple_charge_stats(stage1->mutable_charge_stats());
  wf_mean_phy_pev_->mean_waveforms(stage1->mutable_mean_wf_physics());
  wf_mean_ped_pev_->mean_waveforms(stage1->mutable_mean_wf_pedestal());
  wf_mean_ext_pev_->mean_waveforms(stage1->mutable_mean_wf_external_flasher());
  wf_mean_int_pev_->mean_waveforms(stage1->mutable_mean_wf_internal_flasher());

  stage1->set_run_number(stage1->run_config().run_number());
  stage1->set_run_start_time(stage1->run_config().run_start_time().time_ns());
  stage1->set_run_start_time_string(
    calin::util::timestamp::Timestamp(stage1->run_config().run_start_time().time_ns()).as_string());
  stage1->set_telescope_id(stage1->run_config().telescope_id());
  stage1->set_filename(stage1->run_config().filename());

  calin::provenance::anthology::get_current_anthology(stage1->mutable_provenance_anthology());

  return stage1;
}
