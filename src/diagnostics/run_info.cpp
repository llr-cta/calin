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

#include <diagnostics/run_info.hpp>

using namespace calin::iact_data::diagnostics;

RunInfoDiagnosticsVisitor::RunInfoDiagnosticsVisitor()
{

}

RunInfoDiagnosticsVisitor::~RunInfoDiagnosticsVisitor()
{

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
  return true;
}

bool RunInfoDiagnosticsVisitor::leave_telescope_run()
{
  return true;
}

bool RunInfoDiagnosticsVisitor::visit_telescope_event(uint64_t seq_index,
  calin::ix::iact_data::telescope_event::TelescopeEvent* event)
{
  results_.set_num_events_found(results_.num_events_found() + 1);
  event_number_hist_.insert(event->local_event_number());
  event_number_hist_.insert(event->elapsed_event_time().time_ns()*1e-9);
  return true;
}

bool RunInfoDiagnosticsVisitor::merge_results()
{
  parent_->results_.IntegrateFrom(results_);
  return true;
}
