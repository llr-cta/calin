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

#pragma once

#include <math/histogram.hpp>
#include <iact_data/event_visitor.hpp>
#include <diagnostics/run_info.pb.h>

namespace calin { namespace diagnostics { namespace run_info {

#ifndef SWIG
class ModuleCounterProcessor {
public:
  virtual ~ModuleCounterProcessor();
  virtual int64_t processCounters(std::vector<int64_t>& values, int64_t event_number) = 0;
};

class DirectModuleCounterProcessor: public ModuleCounterProcessor
{
public:
  virtual ~DirectModuleCounterProcessor();
  int64_t processCounters(std::vector<int64_t>& values, int64_t event_number) override;
};

class RelativeToEventNumberModuleCounterProcessor: public ModuleCounterProcessor
{
public:
  virtual ~RelativeToEventNumberModuleCounterProcessor();
  int64_t processCounters(std::vector<int64_t>& values, int64_t event_number) override;
};

class RelativeToMedianModuleCounterProcessor: public ModuleCounterProcessor
{
public:
  virtual ~RelativeToMedianModuleCounterProcessor();
  int64_t processCounters(std::vector<int64_t>& values, int64_t event_number) override;
private:
  std::map<int64_t, unsigned> median_find_;
};

#endif // not defined SWIG

class RunInfoDiagnosticsVisitor:
  public calin::iact_data::event_visitor::ParallelEventVisitor
{
public:
  RunInfoDiagnosticsVisitor(const calin::ix::diagnostics::run_info::RunInfoConfig& config = default_config());

  virtual ~RunInfoDiagnosticsVisitor();

  RunInfoDiagnosticsVisitor* new_sub_visitor(
    const std::map<calin::iact_data::event_visitor::ParallelEventVisitor*,
        calin::iact_data::event_visitor::ParallelEventVisitor*>&
      antecedent_visitors = { }) override;

  bool visit_telescope_run(
    const calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration* run_config,
    calin::iact_data::event_visitor::EventLifetimeManager* event_lifetime_manager) override;

  bool leave_telescope_run() override;

  bool visit_telescope_event(uint64_t seq_index,
    calin::ix::iact_data::telescope_event::TelescopeEvent* event) override;

  bool merge_results() override;

  const calin::ix::diagnostics::run_info::RunInfo& run_info();
  const calin::ix::diagnostics::run_info::PartialRunInfo& partial_run_info();

  static calin::ix::diagnostics::run_info::RunInfoConfig default_config() {
    calin::ix::diagnostics::run_info::RunInfoConfig config;
    return config;
  }

private:


  void integrate_histograms();
  void integrate_partials ();

  RunInfoDiagnosticsVisitor* parent_ = nullptr;
  calin::ix::diagnostics::run_info::RunInfoConfig config_ = default_config();

  std::vector<int64_t> mod_counter_values_;
  std::vector<unsigned> mod_counter_id_;
  std::vector<unsigned> mod_counter_mode_;
  std::vector<ModuleCounterProcessor*> mod_counter_processor_;

  calin::math::histogram::Histogram1D event_number_hist_ { 1.0e4, 0.0, 1.0e9, 0.0 };
  calin::math::histogram::Histogram1D elapsed_time_hist_ { 1.0, -60.0, 7200.0, 0.0 };
  google::protobuf::Arena* arena_ = nullptr;
  calin::ix::diagnostics::run_info::RunInfo* results_ = nullptr;
  calin::ix::diagnostics::run_info::PartialRunInfo* partials_ = nullptr;
  calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration* run_config_ = nullptr;
};

} } } // namespace calin::diagnostics::run_info
