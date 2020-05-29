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
  module_tests_.clear();

  const google::protobuf::RepeatedPtrField<calin::ix::diagnostics::
    clock_regression::SingleClockRegressionConfig>* cam_tests = nullptr;
  const google::protobuf::RepeatedPtrField<calin::ix::diagnostics::
    clock_regression::SingleClockRegressionConfig>* mod_tests = nullptr;
  if(config_.camera_clocks_size() > 0) {
    cam_tests = &config_.camera_clocks();
  }
  if(config_.module_clocks_size() > 0) {
    mod_tests = &config_.module_clocks();
  }

  switch(run_config->camera_layout().camera_type()) {
  case calin::ix::iact_data::instrument_layout::CameraLayout::NECTARCAM:
  case calin::ix::iact_data::instrument_layout::CameraLayout::NECTARCAM_TESTBENCH_19CHANNEL:
  case calin::ix::iact_data::instrument_layout::CameraLayout::NECTARCAM_TESTBENCH_61CHANNEL:
    if(cam_tests == nullptr and config_.default_nectarcam_camera_clocks_size() > 0) {
      cam_tests = &config_.default_nectarcam_camera_clocks();
    }
    if(mod_tests == nullptr and config_.default_nectarcam_module_clocks_size() > 0) {
      mod_tests = &config_.default_nectarcam_module_clocks();
    }
    break;
  case calin::ix::iact_data::instrument_layout::CameraLayout::LSTCAM:
    if(cam_tests == nullptr and config_.default_lstcam_camera_clocks_size() > 0) {
      cam_tests = &config_.default_lstcam_camera_clocks();
    }
    if(mod_tests == nullptr and config_.default_lstcam_module_clocks_size() > 0) {
      mod_tests = &config_.default_lstcam_module_clocks();
    }
    break;
  default:
    break;
  }

  if(cam_tests) {
    for(auto& cam_test : (*cam_tests)) {
      ClockTest ct;
      ct.config = &cam_test;
      camera_tests_.push_back(ct);
    }
  }

  if(mod_tests) {
    for(auto& mod_test : (*mod_tests)) {
      ModuleClockTest mct;
      mct.config = &mod_test;
      mct.modules.resize(run_config->configured_module_id_size());
      module_tests_.push_back(mct);
    }
  }

  return true;
}

bool ClockRegressionParallelEventVisitor::leave_telescope_run()
{
  return true;
}

namespace {
  const calin::ix::iact_data::telescope_event::Clock* find_clock(unsigned id,
    const google::protobuf::RepeatedPtrField<calin::ix::iact_data::telescope_event::Clock>& clocks)
  {
    if(id < clocks.size() and clocks[id].clock_id() == id) {
      return &clocks[id];
    }
    for(auto& iclock : clocks) {
      if(iclock.clock_id() == id) {
        return &iclock;
      }
    }
    return nullptr;
  }

  void accumulate_clock(int64_t master_time, uint64_t local_event_number,
    const calin::ix::iact_data::telescope_event::Clock& clock,
    const calin::ix::diagnostics::clock_regression::SingleClockRegressionConfig& config,
    std::map<int, calin::math::least_squares::I64LinearRegressionAccumulator*>& bins)
  {
    int ibin;
    switch(config.partition_mode()) {
    case calin::ix::diagnostics::clock_regression::PARTITION_BY_CLOCK_SEQUENCE_ID:
      ibin = clock.time_sequence_id();
      break;
    case calin::ix::diagnostics::clock_regression::PARTITION_BY_LOCAL_EVENT_NUMBER:
      ibin = local_event_number/config.partition_bin_size();
      break;
    case calin::ix::diagnostics::clock_regression::PARTITION_BY_MASTER_CLOCK:
      ibin = master_time/config.partition_bin_size();
      break;
    case calin::ix::diagnostics::clock_regression::SINGLE_PARTITION:
    default:
      ibin = 0;
      break;
    }
    auto* accumulator = bins[ibin];
    if(accumulator == nullptr) {
      bins[ibin] = accumulator =
        new calin::math::least_squares::I64LinearRegressionAccumulator();
    }
    accumulator->accumulate(master_time, clock.time_value());
  }
}

bool ClockRegressionParallelEventVisitor::visit_telescope_event(uint64_t seq_index,
  calin::ix::iact_data::telescope_event::TelescopeEvent* event)
{
  const auto* master_clock = find_clock(config_.master_clock_id(),event->camera_clock());
  if(master_clock == nullptr) {
    return true;
  }

  for(auto& ct : camera_tests_)
  {
    const auto* test_clock = find_clock(ct.config->clock_id(),event->camera_clock());
    if(test_clock) {
      int64_t master_time = master_clock->time_value();
      if(ct.config->master_clock_divisor() > 1) {
        master_time /= ct.config->master_clock_divisor();
      }
      accumulate_clock(master_time, event->local_event_number(), *test_clock, *ct.config, ct.bins);
    }
  }

  for(auto& mt : module_tests_)
  {
    int64_t master_time = master_clock->time_value();
    if(mt.config->master_clock_divisor() > 1) {
      master_time /= mt.config->master_clock_divisor();
    }

    for(auto& imod : event->module_clock()) {
      const auto* test_clock = find_clock(mt.config->clock_id(),imod.clock());
      if(test_clock) {
        accumulate_clock(master_time, event->local_event_number(), *test_clock,
          *mt.config, mt.modules[imod.module_id()].bins);
      }
    }
  }

  return true;
}

void ClockRegressionParallelEventVisitor::merge_into(ClockTest* to, const ClockTest& from)
{
  for(auto& ibin : from.bins) {
    if(to->bins[ibin.first] == nullptr) {
      to->bins[ibin.first] = new calin::math::least_squares::I64LinearRegressionAccumulator();
    }
    ibin.second->integrate_into(*to->bins[ibin.first]);
  }
}

bool ClockRegressionParallelEventVisitor::merge_results()
{
  if(parent_) {
    for(unsigned ict=0; ict<camera_tests_.size(); ict++) {
      merge_into(&parent_->camera_tests_[ict], camera_tests_[ict]);
    }
    for(unsigned imt=0; imt<module_tests_.size(); imt++) {
      for(unsigned imod=0; imod<module_tests_[imt].modules.size(); imod++) {
        merge_into(&parent_->module_tests_[imt].modules[imod],
          module_tests_[imt].modules[imod]);
      }
    }
  }
  return true;
}
