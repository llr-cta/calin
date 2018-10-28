/*

   calin/diagnostics/run_info.hpp -- Stephen Fegan -- 2018-10-26

   Run info diagnostics

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

#include <util/log.hpp>
#include <diagnostics/run_info.hpp>
#include <diagnostics/range.hpp>

using namespace calin::util::log;
using namespace calin::diagnostics::run_info;

RunInfoDiagnosticsVisitor::RunInfoDiagnosticsVisitor():
  calin::iact_data::event_visitor::ParallelEventVisitor(),
  arena_(new google::protobuf::Arena),
  results_(google::protobuf::Arena::CreateMessage<calin::ix::diagnostics::run_info::RunInfo>(arena_)),
  partials_(google::protobuf::Arena::CreateMessage<calin::ix::diagnostics::run_info::PartialRunInfo>(arena_))
{
  // nothing to see here
}

RunInfoDiagnosticsVisitor::~RunInfoDiagnosticsVisitor()
{
  delete arena_;
}

RunInfoDiagnosticsVisitor* RunInfoDiagnosticsVisitor::new_sub_visitor(
  const std::map<calin::iact_data::event_visitor::ParallelEventVisitor*,
      calin::iact_data::event_visitor::ParallelEventVisitor*>&
    antecedent_visitors)
{
  RunInfoDiagnosticsVisitor* child = new RunInfoDiagnosticsVisitor();
  child->parent_ = this;
  return child;
}

bool RunInfoDiagnosticsVisitor::visit_telescope_run(
  const calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration* run_config,
  calin::iact_data::event_visitor::EventLifetimeManager* event_lifetime_manager)
{
  results_->Clear();
  partials_->Clear();
  event_number_hist_.clear();
  elapsed_time_hist_.clear();
  return true;
}

bool RunInfoDiagnosticsVisitor::leave_telescope_run()
{
  return true;
}

bool RunInfoDiagnosticsVisitor::visit_telescope_event(uint64_t seq_index,
  calin::ix::iact_data::telescope_event::TelescopeEvent* event)
{
  results_->set_num_events_found(results_->num_events_found() + 1);
  event_number_hist_.insert(event->local_event_number());
  elapsed_time_hist_.insert(event->elapsed_event_time().time_ns()*1e-9);

  partials_->add_event_number_sequence(event->local_event_number());

  // UCTS
  if(!event->has_cdts_data()) {
    results_->set_num_events_missing_cdts(results_->num_events_missing_cdts() + 1);
  }
  calin::diagnostics::range::encode_value(
    partials_->mutable_cdts_presence(), event->has_cdts_data());

  // TIB
  if(!event->has_tib_data()) {
    results_->set_num_events_missing_tib(results_->num_events_missing_tib() + 1);
  }
  calin::diagnostics::range::encode_value(
    partials_->mutable_tib_presence(), event->has_tib_data());

  // SWAT
  if(!event->has_swat_data()) {
    results_->set_num_events_missing_swat(results_->num_events_missing_swat() + 1);
  }
  calin::diagnostics::range::encode_value(
    partials_->mutable_swat_presence(), event->has_swat_data());

  return true;
}

const calin::ix::diagnostics::run_info::RunInfo& RunInfoDiagnosticsVisitor::run_info()
{
  integrate_data();
  return *results_;
}

const calin::ix::diagnostics::run_info::PartialRunInfo& RunInfoDiagnosticsVisitor::partial_run_info()
{
  return *partials_;
}

bool RunInfoDiagnosticsVisitor::merge_results()
{
  integrate_data();
  parent_->results_->IntegrateFrom(*results_);
  parent_->partials_->IntegrateFrom(*partials_);
  return true;
}

void RunInfoDiagnosticsVisitor::integrate_data()
{
  auto* event_number_hist_data = event_number_hist_.dump_as_proto();
  results_->mutable_event_number_histogram()->IntegrateFrom(*event_number_hist_data);
  delete event_number_hist_data;
  event_number_hist_.clear();

  auto* elapsed_time_hist_data = elapsed_time_hist_.dump_as_proto();
  results_->mutable_elapsed_time_histogram()->IntegrateFrom(*elapsed_time_hist_data);
  delete elapsed_time_hist_data;
  elapsed_time_hist_.clear();
}
