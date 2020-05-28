/*

   calin/diagnostics/clock_regression.cpp -- Stephen Fegan -- 2020-05-25

   Clock regression visitor

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

#include <Eigen/Dense>
#include <diagnostics/clock_regression.hpp>

using namespace calin::diagnostics::clock_regression;

ClockRegressionParallelEventVisitor::
ClockRegressionParallelEventVisitor(
    const calin::ix::diagnostics::clock_regression::ClockRegressionConfig& config):
  calin::iact_data::event_visitor::ParallelEventVisitor(), config_(config)
{
  // nothing to see here
}

ClockRegressionParallelEventVisitor::~ClockRegressionParallelEventVisitor()
{
  // nothing to see here
}

ClockRegressionParallelEventVisitor*
ClockRegressionParallelEventVisitor::new_sub_visitor(
  std::map<calin::iact_data::event_visitor::ParallelEventVisitor*,
    calin::iact_data::event_visitor::ParallelEventVisitor*> antecedent_visitors)
{
  auto* sub_visior = new ClockRegressionParallelEventVisitor(config_);
  sub_visior->parent_ = this;
  return sub_visior;
}

bool ClockRegressionParallelEventVisitor::visit_telescope_run(
  const calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration* run_config,
  calin::iact_data::event_visitor::EventLifetimeManager* event_lifetime_manager)
{
  camera_tests_.clear();

  const google::protobuf::RepeatedPtrField<calin::ix::diagnostics::
    clock_regression::SingleClockRegressionConfig>* vec_tests = nullptr;
  if(config_.camera_clocks_size() > 0) {
    vec_tests = &config_.camera_clocks();
  }

  // repeated SingleClockRegressionConfig camera_clocks       = 10;
  // repeated SingleClockRegressionConfig default_netarcam_camera_clocks = 100;
  // repeated SingleClockRegressionConfig default_lstcam_camera_clocks = 110;
  //
  // repeated SingleClockRegressionConfig module_clocks       = 11;
  // repeated SingleClockRegressionConfig default_netarcam_module_clocks = 101;
  // repeated SingleClockRegressionConfig default_lstcam_module_clocks = 111;
  //
  // std::vector<ModuleClockTest> module_tests_;


  return true;
}

bool ClockRegressionParallelEventVisitor::leave_telescope_run()
{
  return true;
}

bool ClockRegressionParallelEventVisitor::visit_telescope_event(uint64_t seq_index,
  calin::ix::iact_data::telescope_event::TelescopeEvent* event)
{
  return true;
}

bool ClockRegressionParallelEventVisitor::merge_results()
{
  return true;
}
