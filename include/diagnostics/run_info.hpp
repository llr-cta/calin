/*

   calin/diagnostics/run_info.hpp -- Stephen Fegan -- 2018-10-26

   Run info diagnostics

   Copyright 2018, Stephen Fegan <sfegan@llr.in2p3.fr>
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

class RunInfoDiagnosticsParallelEventVisitor:
  public calin::iact_data::event_visitor::ParallelEventVisitor
{
public:
  RunInfoDiagnosticsParallelEventVisitor(const calin::ix::diagnostics::run_info::RunInfoConfig& config = default_config());

  virtual ~RunInfoDiagnosticsParallelEventVisitor();

  RunInfoDiagnosticsParallelEventVisitor* new_sub_visitor(
    std::map<calin::iact_data::event_visitor::ParallelEventVisitor*,
        calin::iact_data::event_visitor::ParallelEventVisitor*>
      antecedent_visitors = { }) override;

  bool visit_telescope_run(
    const calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration* run_config,
    calin::iact_data::event_visitor::EventLifetimeManager* event_lifetime_manager,
    calin::ix::provenance::chronicle::ProcessingRecord* processing_record = nullptr) override;

  bool leave_telescope_run(
    calin::ix::provenance::chronicle::ProcessingRecord* processing_record = nullptr) override;

  bool visit_telescope_event(uint64_t seq_index,
    calin::ix::iact_data::telescope_event::TelescopeEvent* event) override;

  bool merge_results() override;

  int64_t min_event_time() const { return results_->min_event_time(); }
  int64_t max_event_time() const { return results_->max_event_time(); }

#ifndef SWIG
  calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration* run_config(
    calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration* rc = nullptr) const;
  calin::ix::diagnostics::run_info::RunInfo* run_info(calin::ix::diagnostics::run_info::RunInfo* ri = nullptr) const;
#else
  calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration* run_config() const;
  void run_config(calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration* rc) const;
  calin::ix::diagnostics::run_info::RunInfo* run_info() const;
  void run_info(calin::ix::diagnostics::run_info::RunInfo* ri) const;
#endif

  const calin::ix::diagnostics::run_info::PartialRunInfo& partial_run_info() const;

  static calin::ix::diagnostics::run_info::RunInfoConfig default_config();

private:
  void integrate_partials ();

  RunInfoDiagnosticsParallelEventVisitor* parent_ = nullptr;
  calin::ix::diagnostics::run_info::RunInfoConfig config_ = default_config();

  std::vector<int64_t> mod_counter_values_;
  std::vector<unsigned> mod_counter_id_;
  std::vector<unsigned> mod_counter_mode_;
  std::vector<ModuleCounterProcessor*> mod_counter_processor_;
  std::vector<std::vector<calin::math::histogram::SimpleHist*> > mod_counter_hist_;

  google::protobuf::Arena* arena_ = nullptr;
  calin::ix::diagnostics::run_info::RunInfo* results_ = nullptr;
  calin::ix::diagnostics::run_info::PartialRunInfo* partials_ = nullptr;
  calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration* run_config_ = nullptr;
  calin::math::histogram::SimpleHist trigger_type_code_hist_ = { 1.0, -1.0, 256.0 };
  calin::math::histogram::SimpleHist trigger_type_code_diff_hist_ = { 1.0, -256.0, 256.0 };
};

} } } // namespace calin::diagnostics::run_info
