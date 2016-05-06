  /*

   calin/diagnostics/functinal.hpp -- Stephen Fegan -- 2016-02-10

   Functional diagnostics visitor

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, Ecole polytechnique, CNRS/IN2P3, Universite Paris-Saclay

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

#include <iact_data/event_visitor.hpp>
#include <iact_data/functional_event_visitor.hpp>
#include <diagnostics/functional.pb.h>
#include <math/histogram.hpp>

namespace calin { namespace diagnostics { namespace functional {

class FunctionalStatsVisitor:
  public iact_data::event_visitor::TelescopeEventVisitor
{
public:
  FunctionalStatsVisitor(
    calin::iact_data::functional_event_visitor::
      DualValueInt32FunctionalTelescopeEventVisitor* value_supplier,
    const calin::ix::diagnostics::functional::
      FunctionalStatsVisitorConfig& config = default_config());

  virtual ~FunctionalStatsVisitor();

  bool demand_waveforms() override;
  bool is_parallelizable() override;
  FunctionalStatsVisitor* new_sub_visitor(
    const std::map<calin::iact_data::event_visitor::TelescopeEventVisitor*,
        calin::iact_data::event_visitor::TelescopeEventVisitor*>&
      antecedent_visitors) override;

  bool visit_telescope_run(
    const calin::ix::iact_data::telescope_run_configuration::
      TelescopeRunConfiguration* run_config) override;
  bool leave_telescope_run() override;

  bool visit_telescope_event(uint64_t seq_index,
    calin::ix::iact_data::telescope_event::TelescopeEvent* event) override;
  bool leave_telescope_event() override;

  bool visit_waveform(unsigned ichan,
    calin::ix::iact_data::telescope_event::ChannelWaveform* high_gain,
    calin::ix::iact_data::telescope_event::ChannelWaveform* low_gain) override;

  bool merge_results() override;

  calin::ix::diagnostics::functional::CameraIntFunctionalRawStats results()
  {
    return results_;
  }

  static Eigen::VectorXd channel_mean(
    const ix::diagnostics::functional::OneGainIntFunctionalRawStats* stat);

  static Eigen::VectorXd channel_var(
    const ix::diagnostics::functional::OneGainIntFunctionalRawStats* stat);

  static Eigen::MatrixXd channel_cov(
    const ix::diagnostics::functional::OneGainIntFunctionalRawStats* stat);

  static calin::ix::diagnostics::functional::FunctionalStatsVisitorConfig
  default_config()
  {
    calin::ix::diagnostics::functional::FunctionalStatsVisitorConfig cfg;
    cfg.set_calculate_covariance(true);
    return cfg;
  }

protected:
  void visit_one_waveform(
    const calin::ix::iact_data::telescope_event::ChannelWaveform* wf,
    unsigned index, std::vector<int>& mask, std::vector<int32_t>& signal);

  void process_one_gain(const std::vector<int>& mask,
      const std::vector<int32_t>& signal,
    std::vector<calin::math::histogram::SimpleHist>& hist,
    calin::ix::diagnostics::functional::OneGainIntFunctionalRawStats* stats);

  void merge_one_gain(
    const ix::diagnostics::functional::OneGainIntFunctionalRawStats* from,
    ix::diagnostics::functional::OneGainIntFunctionalRawStats* to);

  calin::iact_data::functional_event_visitor::
    DualValueInt32FunctionalTelescopeEventVisitor* value_supplier_;
  calin::ix::diagnostics::functional::
    FunctionalStatsVisitorConfig config_;

  std::vector<int> high_gain_mask_;
  std::vector<int32_t> high_gain_signal_;
  std::vector<int> low_gain_mask_;
  std::vector<int32_t> low_gain_signal_;

  std::vector<calin::math::histogram::SimpleHist> high_gain_hist_;
  std::vector<calin::math::histogram::SimpleHist> low_gain_hist_;

  FunctionalStatsVisitor* parent_ = nullptr;
  calin::ix::diagnostics::functional::CameraIntFunctionalRawStats results_;
  const ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration*
    run_config_ = nullptr;
  bool calculate_covariance_ = false;
};

#if 0

class FunctionalCaptureVisitor:
  public iact_data::event_visitor::TelescopeEventVisitor
{
public:
  FunctionalCaptureVisitor(
    calin::iact_data::functional_event_visitor::
      DualValueInt32FunctionalTelescopeEventVisitor* value_supplier,
    const calin::ix::diagnostics::functional::
      FunctionalCaptureVisitorConfig& config = default_config());

  virtual ~FunctionalCaptureVisitor();

  bool demand_waveforms() override;
  bool is_parallelizable() override;
  FunctionalCaptureVisitor* new_sub_visitor(
    const std::map<calin::iact_data::event_visitor::TelescopeEventVisitor*,
        calin::iact_data::event_visitor::TelescopeEventVisitor*>&
      antecedent_visitors) override;

  bool visit_telescope_run(
    const calin::ix::iact_data::telescope_run_configuration::
      TelescopeRunConfiguration* run_config) override;
  bool leave_telescope_run() override;

  bool visit_telescope_event(
    calin::ix::iact_data::telescope_event::TelescopeEvent* event) override;
  bool leave_telescope_event() override;

  bool visit_waveform(unsigned ichan,
    calin::ix::iact_data::telescope_event::ChannelWaveform* high_gain,
    calin::ix::iact_data::telescope_event::ChannelWaveform* low_gain) override;

  bool merge_results() override;

  calin::ix::diagnostics::functional::CameraIntFunctionalRawStats results()
  {
    return results_;
  }

  static calin::ix::diagnostics::functional::FunctionalCaptureVisitorConfig
  default_config()
  {
    calin::ix::diagnostics::functional::FunctionalCaptureVisitorConfig cfg;
    return cfg;
  }

protected:
  void visit_one_waveform(
    const calin::ix::iact_data::telescope_event::ChannelWaveform* wf,
    unsigned index, std::vector<int>& mask, std::vector<int32_t>& signal);

  void process_one_gain(const std::vector<int>& mask,
      const std::vector<int32_t>& signal,
    std::vector<calin::math::histogram::SimpleHist>& hist,
    calin::ix::diagnostics::functional::OneGainIntFunctionalRawStats* stats);

  void merge_one_gain(
    const ix::diagnostics::functional::OneGainIntFunctionalRawStats* from,
    ix::diagnostics::functional::OneGainIntFunctionalRawStats* to);

  calin::iact_data::functional_event_visitor::
    DualValueInt32FunctionalTelescopeEventVisitor* value_supplier_;
  calin::ix::diagnostics::functional::
    FunctionalStatsVisitorConfig config_;

  std::vector<int> high_gain_mask_;
  std::vector<int32_t> high_gain_signal_;
  std::vector<int> low_gain_mask_;
  std::vector<int32_t> low_gain_signal_;

  std::vector<calin::math::histogram::SimpleHist> high_gain_hist_;
  std::vector<calin::math::histogram::SimpleHist> low_gain_hist_;

  FunctionalStatsVisitor* parent_ = nullptr;
  calin::ix::diagnostics::functional::CameraIntFunctionalRawStats results_;
  const ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration*
    run_config_ = nullptr;
  bool calculate_covariance_ = false;
};

#endif

} } } // namespace calin::diagnostics::functional_diagnostics
