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
#include <diagnostics/value_capture.hpp>
#include <diagnostics/functional.pb.h>
#include <math/histogram.hpp>

namespace calin { namespace diagnostics { namespace functional {

class SingleFunctionalValueSupplierVisitor:
  public calin::diagnostics::value_capture::ValueSupplierVisitor<int32_t>
{
public:
  SingleFunctionalValueSupplierVisitor(
    calin::iact_data::functional_event_visitor::
      DualGainInt32FunctionalTelescopeEventVisitor* value_supplier,
    unsigned chan, bool low_gain = false);
  virtual ~SingleFunctionalValueSupplierVisitor();
  bool visit_telescope_event(uint64_t seq_index,
    calin::ix::iact_data::telescope_event::TelescopeEvent* event) override;
  bool visit_waveform(unsigned ichan,
    calin::ix::iact_data::telescope_event::ChannelWaveform* high_gain,
    calin::ix::iact_data::telescope_event::ChannelWaveform* low_gain) override;
  bool get_value(int32_t& value) override;
private:
  calin::iact_data::functional_event_visitor::
    DualGainInt32FunctionalTelescopeEventVisitor* value_supplier_;
  unsigned chan_;
  bool low_gain_;
  int32_t value_;
  bool has_value_ = false;
};

template<typename DualGainFunctionalVisitor, typename Results>
class FunctionalStatsVisitor:
  public iact_data::event_visitor::TelescopeEventVisitor
{
public:
  CALIN_TYPEALIAS(functional_value_type,
    typename DualGainFunctionalVisitor::value_type);

  FunctionalStatsVisitor(
    DualGainFunctionalVisitor* value_supplier,
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

  Results results()
  {
    return results_;
  }

  static calin::ix::diagnostics::functional::FunctionalStatsVisitorConfig
  default_config()
  {
    calin::ix::diagnostics::functional::FunctionalStatsVisitorConfig cfg;
    cfg.set_calculate_covariance(true);
    cfg.mutable_hist_config()->set_dxval(1.0);
    cfg.mutable_hist_config()->set_xval_align(0.5);
    return cfg;
  }

protected:
  template<typename OneGainRawStats>
  void process_one_gain(const std::vector<int>& mask,
    const std::vector<functional_value_type>& signal,
    std::vector<calin::math::histogram::SimpleHist>& hist,
    calin::math::histogram::SimpleHist& mean_hist,
    OneGainRawStats* stats);

  DualGainFunctionalVisitor* value_supplier_;
  calin::ix::diagnostics::functional::
    FunctionalStatsVisitorConfig config_;

  std::vector<int> high_gain_mask_;
  std::vector<functional_value_type> high_gain_signal_;
  std::vector<int> low_gain_mask_;
  std::vector<functional_value_type> low_gain_signal_;

  std::vector<calin::math::histogram::SimpleHist> high_gain_hist_;
  std::vector<calin::math::histogram::SimpleHist> low_gain_hist_;

  calin::math::histogram::SimpleHist high_gain_mean_hist_;
  calin::math::histogram::SimpleHist low_gain_mean_hist_;

  FunctionalStatsVisitor* parent_ = nullptr;
  Results results_;
  const ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration*
    run_config_ = nullptr;
  bool calculate_covariance_ = false;
};

template<typename OneGainRawStats>
Eigen::VectorXd channel_mean(const OneGainRawStats* stat);

template<typename OneGainRawStats>
Eigen::VectorXd channel_var(const OneGainRawStats* stat);

template<typename OneGainRawStats>
Eigen::MatrixXd channel_cov(const OneGainRawStats* stat);

template<typename OneGainRawStats>
Eigen::MatrixXd channel_cov_frac(const OneGainRawStats* stat);

template<typename OneGainRawStats>
double mean_of_mean_over_channels(const OneGainRawStats* stat);

template<typename OneGainRawStats>
double var_of_mean_over_channels(const OneGainRawStats* stat);

} } } // namespace calin::diagnostics::functional_diagnostics

#include <diagnostics/functional_stats_impl.hpp>

#ifndef SWIG
#ifndef CALIN_DIAGNOSTICS_FUNCTIONAL_NO_EXTERN

extern template
class calin::diagnostics::functional::FunctionalStatsVisitor<
  calin::iact_data::functional_event_visitor::
    DualGainInt32FunctionalTelescopeEventVisitor,
  calin::ix::diagnostics::functional::CameraIntFunctionalRawStats>;

extern template
class calin::diagnostics::functional::FunctionalStatsVisitor<
  calin::iact_data::functional_event_visitor::
    DualGainDoubleFunctionalTelescopeEventVisitor,
  calin::ix::diagnostics::functional::CameraDoubleFunctionalRawStats>;

#endif // #ifdef CALIN_VALUE_CAPTURE_NO_EXTERN
#endif // #ifdef SWIG
