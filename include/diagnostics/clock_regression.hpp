/*

   calin/diagnostics/clock_regression.hpp -- Stephen Fegan -- 2020-05-25

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

#pragma once

// #include <diagnostics/clock_regression.pb.h>
#include <math/histogram.hpp>
#include <math/least_squares.hpp>
#include <iact_data/event_visitor.hpp>
#include <diagnostics/clock_regression.pb.h>

namespace calin { namespace diagnostics { namespace clock_regression {

class ClockRegressionParallelEventVisitor:
  public calin::iact_data::event_visitor::ParallelEventVisitor
{
public:
  using RegressionAccumulator = calin::math::least_squares::I64LinearRegressionAccumulatorIgnoringFirstDatum;

  ClockRegressionParallelEventVisitor(
    const calin::ix::diagnostics::clock_regression::ClockRegressionConfig& config = default_config());

  virtual ~ClockRegressionParallelEventVisitor();

  ClockRegressionParallelEventVisitor* new_sub_visitor(
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

  static calin::ix::diagnostics::clock_regression::ClockRegressionConfig default_config();

#ifndef SWIG
  calin::ix::diagnostics::clock_regression::ClockRegressionResults* clock_regression(
    calin::ix::diagnostics::clock_regression::ClockRegressionResults* results = nullptr) const;
#else
  calin::ix::diagnostics::clock_regression::ClockRegressionResults* clock_regression() const;
  void clock_regression(calin::ix::diagnostics::clock_regression::ClockRegressionResults* results) const;
#endif

private:

  struct ClockTest {
    ~ClockTest() { for(auto& ibin : bins) { delete ibin.second; } }
    calin::ix::diagnostics::clock_regression::SingleClockRegressionConfig config;
    std::map<int, RegressionAccumulator*> bins;
  };

  struct ModuleClockTest {
    calin::ix::diagnostics::clock_regression::SingleClockRegressionConfig config;
    std::vector<ClockTest> modules;
  };

  void merge_into(ClockTest* to, const ClockTest& from);

  void transfer_clock_results(calin::ix::diagnostics::clock_regression::SingleClockRegressionResults* res,
    const ClockTest& ct) const;

  ClockRegressionParallelEventVisitor* parent_ = nullptr;

  calin::ix::diagnostics::clock_regression::ClockRegressionConfig config_;
  int rebalance_ = 0;
  std::string principal_clock_name_ = {};
  int principal_clock_id_ = 0;
  std::vector<ClockTest> camera_tests_;
  std::vector<ModuleClockTest> module_tests_;
};

} } } // namespace calin::diagnostics::clock_regression
