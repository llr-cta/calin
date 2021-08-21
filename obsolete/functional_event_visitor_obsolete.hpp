/*

   calin/iact_data/functional_event_visitor.hpp -- Stephen Fegan -- 2016-04-13

   Vistor to event data that calculates some quantity from the waveforms

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <cstdlib>
#include <numeric>

#include <iact_data/event_visitor.hpp>
#include <iact_data/functional_event_visitor.pb.h>

namespace calin { namespace iact_data { namespace functional_event_visitor {

template<typename T> class BasicDualGainFunctionalTelescopeEventVisitor:
  public calin::iact_data::event_visitor::TelescopeEventVisitor
{
public:
  CALIN_TYPEALIAS(value_type, T);
  virtual ~BasicDualGainFunctionalTelescopeEventVisitor() {
    // nothing to see here
  }
  virtual T low_gain_value() = 0;
  virtual T high_gain_value() = 0;
};

CALIN_TYPEALIAS(DualGainInt32FunctionalTelescopeEventVisitor,
  BasicDualGainFunctionalTelescopeEventVisitor<int32_t>);
CALIN_TYPEALIAS(DualGainDoubleFunctionalTelescopeEventVisitor,
  BasicDualGainFunctionalTelescopeEventVisitor<double>);

#ifdef SWIG
} } } // namespace calin::iact_data::functional_event_visitor
%template(DualGainInt32FunctionalTelescopeEventVisitor)
  calin::iact_data::functional_event_visitor::BasicDualGainFunctionalTelescopeEventVisitor<int32_t>;
%template(DualGainDoubleFunctionalTelescopeEventVisitor)
  calin::iact_data::functional_event_visitor::BasicDualGainFunctionalTelescopeEventVisitor<double>;
namespace calin { namespace iact_data { namespace functional_event_visitor {
#endif

class FixedWindowSumFunctionalTelescopeEventVisitor:
  public DualGainInt32FunctionalTelescopeEventVisitor
{
public:
  FixedWindowSumFunctionalTelescopeEventVisitor(
    calin::ix::iact_data::functional_event_visitor::
      FixedWindowSumFunctionalTelescopeEventVisitorConfig config =
        default_config());
  virtual ~FixedWindowSumFunctionalTelescopeEventVisitor();

  bool demand_waveforms() override;
  bool is_parallelizable() override;

  FixedWindowSumFunctionalTelescopeEventVisitor* new_sub_visitor(
      const std::map<calin::iact_data::event_visitor::TelescopeEventVisitor*,
          calin::iact_data::event_visitor::TelescopeEventVisitor*>&
        antecedent_visitors) override;

  bool visit_telescope_run(
    const calin::ix::iact_data::telescope_run_configuration::
      TelescopeRunConfiguration* run_config) override;

  bool visit_waveform(unsigned ichan,
    calin::ix::iact_data::telescope_event::ChannelWaveform* high_gain,
    calin::ix::iact_data::telescope_event::ChannelWaveform* low_gain) override;

  int32_t low_gain_value() override;
  int32_t high_gain_value() override;

  static calin::ix::iact_data::functional_event_visitor::
    FixedWindowSumFunctionalTelescopeEventVisitorConfig default_config()
  {
    calin::ix::iact_data::functional_event_visitor::
      FixedWindowSumFunctionalTelescopeEventVisitorConfig cfg;
    cfg.set_integration_0(0);
    cfg.set_integration_n(0);
    return cfg;
  }

private:

  inline int32_t process_one_gain(
    const calin::ix::iact_data::telescope_event::ChannelWaveform* ch_wf) const
  {
    const auto* samples = ch_wf->samples().data() + window_0_;
    int32_t sum = 0;
    for(unsigned i=0; i<window_n_; i++)sum += *(samples++);
    return sum;
  }

  calin::ix::iact_data::functional_event_visitor::
    FixedWindowSumFunctionalTelescopeEventVisitorConfig config_;
  unsigned window_0_ = 0;
  unsigned window_n_ = 0;
  int32_t high_gain_value_ = 0;
  int32_t low_gain_value_ = 0;
};

class SlidingWindowSumFunctionalTelescopeEventVisitor:
  public DualGainInt32FunctionalTelescopeEventVisitor
{
public:
  SlidingWindowSumFunctionalTelescopeEventVisitor(
    calin::ix::iact_data::functional_event_visitor::
      SlidingWindowSumFunctionalTelescopeEventVisitorConfig config =
        default_config());
  virtual ~SlidingWindowSumFunctionalTelescopeEventVisitor();

  bool demand_waveforms() override;
  bool is_parallelizable() override;

  SlidingWindowSumFunctionalTelescopeEventVisitor* new_sub_visitor(
      const std::map<calin::iact_data::event_visitor::TelescopeEventVisitor*,
          calin::iact_data::event_visitor::TelescopeEventVisitor*>&
        antecedent_visitors) override;

  bool visit_telescope_run(
    const calin::ix::iact_data::telescope_run_configuration::
      TelescopeRunConfiguration* run_config) override;

  bool visit_waveform(unsigned ichan,
    calin::ix::iact_data::telescope_event::ChannelWaveform* high_gain,
    calin::ix::iact_data::telescope_event::ChannelWaveform* low_gain) override;

  int32_t low_gain_value() override;
  int32_t high_gain_value() override;

  static calin::ix::iact_data::functional_event_visitor::
    SlidingWindowSumFunctionalTelescopeEventVisitorConfig default_config()
  {
    calin::ix::iact_data::functional_event_visitor::
      SlidingWindowSumFunctionalTelescopeEventVisitorConfig cfg;
    cfg.set_integration_n(0);
    return cfg;
  }

private:

  inline int32_t process_one_gain(
    const calin::ix::iact_data::telescope_event::ChannelWaveform* ch_wf) const
  {
    const unsigned nsample = ch_wf->samples_size();
    const auto*const samples = ch_wf->samples().data();
    const auto*const samples_end = samples + nsample;

    const auto* window_end = samples;
    int32_t sum = *(window_end++);
    for(unsigned i=1; i<window_n_; i++)sum += *(window_end++);

    int32_t max_sum = sum;
    const auto* window_start = samples;
    while(window_end < samples_end) {
      sum += *(window_end++);
      sum -= *(window_start++);
      max_sum = std::max(max_sum, sum);
    }

    return max_sum;
  }

  calin::ix::iact_data::functional_event_visitor::
    SlidingWindowSumFunctionalTelescopeEventVisitorConfig config_;
  unsigned window_n_ = 0;
  int32_t high_gain_value_ = 0;
  int32_t low_gain_value_ = 0;
};

class DifferencingFunctionalTelescopeEventVisitor:
  public DualGainInt32FunctionalTelescopeEventVisitor
{
public:
  DifferencingFunctionalTelescopeEventVisitor(
    DualGainInt32FunctionalTelescopeEventVisitor* sig_value_supplier,
    DualGainInt32FunctionalTelescopeEventVisitor* bkg_value_supplier);
  virtual ~DifferencingFunctionalTelescopeEventVisitor();

  bool demand_waveforms() override;
  bool is_parallelizable() override;

  DifferencingFunctionalTelescopeEventVisitor* new_sub_visitor(
    const std::map<calin::iact_data::event_visitor::TelescopeEventVisitor*,
        calin::iact_data::event_visitor::TelescopeEventVisitor*>&
      antecedent_visitors) override;

  int32_t low_gain_value() override;
  int32_t high_gain_value() override;

private:
  DualGainInt32FunctionalTelescopeEventVisitor* sig_value_supplier_;
  DualGainInt32FunctionalTelescopeEventVisitor* bkg_value_supplier_;
};

class RisetimeTimingFunctionalTelescopeEventVisitor:
  public DualGainDoubleFunctionalTelescopeEventVisitor
{
public:
  RisetimeTimingFunctionalTelescopeEventVisitor();
  virtual ~RisetimeTimingFunctionalTelescopeEventVisitor();

  bool demand_waveforms() override;
  bool is_parallelizable() override;

  RisetimeTimingFunctionalTelescopeEventVisitor* new_sub_visitor(
      const std::map<calin::iact_data::event_visitor::TelescopeEventVisitor*,
          calin::iact_data::event_visitor::TelescopeEventVisitor*>&
        antecedent_visitors) override;

  bool visit_waveform(unsigned ichan,
    calin::ix::iact_data::telescope_event::ChannelWaveform* high_gain,
    calin::ix::iact_data::telescope_event::ChannelWaveform* low_gain) override;

  double low_gain_value() override;
  double high_gain_value() override;

protected:

  inline double process_one_gain(
    const calin::ix::iact_data::telescope_event::ChannelWaveform* ch_wf) const
  {
    const auto*const samples_start = ch_wf->samples().data();
    const auto* samples_end = samples_start + ch_wf->samples_size();
    const auto* samples_max = samples_start;
    for(const auto* samples = samples_start + 1;
        samples < samples_end; ++samples)
      if(*samples > *samples_max)samples_max = samples;
    const auto* samples_min = samples_start;
    samples_end = samples_max;
    for(const auto* samples = samples_start + 1;
        samples < samples_end; ++samples)
      if(*samples < *samples_min)samples_min = samples;
    if(samples_min == samples_max)return -1;
    double samples_mid = 0.5*(double(*samples_min)+double(*samples_max));
    while(double(*samples_max) > samples_mid)--samples_max;
    return (samples_max-samples_start) +
      (samples_mid - double(*samples_max))/double(*(samples_max+1)-*samples_max);
  }

  double high_gain_value_ = 0.0;
  double low_gain_value_ = 0.0;
};

class MeantimeTimingFunctionalTelescopeEventVisitor:
  public DualGainDoubleFunctionalTelescopeEventVisitor
{
public:
  MeantimeTimingFunctionalTelescopeEventVisitor(
    calin::ix::iact_data::functional_event_visitor::
      MeantimeTimingFunctionalTelescopeEventVisitorConfig config =
        default_config());
  virtual ~MeantimeTimingFunctionalTelescopeEventVisitor();

  bool demand_waveforms() override;
  bool is_parallelizable() override;

  MeantimeTimingFunctionalTelescopeEventVisitor* new_sub_visitor(
      const std::map<calin::iact_data::event_visitor::TelescopeEventVisitor*,
          calin::iact_data::event_visitor::TelescopeEventVisitor*>&
        antecedent_visitors) override;

  bool visit_telescope_run(
    const calin::ix::iact_data::telescope_run_configuration::
      TelescopeRunConfiguration* run_config) override;

  bool visit_waveform(unsigned ichan,
    calin::ix::iact_data::telescope_event::ChannelWaveform* high_gain,
    calin::ix::iact_data::telescope_event::ChannelWaveform* low_gain) override;

  double low_gain_value() override;
  double high_gain_value() override;

  static calin::ix::iact_data::functional_event_visitor::
    MeantimeTimingFunctionalTelescopeEventVisitorConfig default_config()
  {
    calin::ix::iact_data::functional_event_visitor::
      MeantimeTimingFunctionalTelescopeEventVisitorConfig cfg;
    cfg.set_pedestal_0(0);
    cfg.set_pedestal_n(0);
    cfg.set_signal_0(0);
    cfg.set_signal_n(0);
    cfg.set_pedestal_decay_constant(0.0);
    return cfg;
  }

protected:

  inline double process_one_gain(double& pedestal,
    const calin::ix::iact_data::telescope_event::ChannelWaveform* ch_wf) const
  {
    const auto*const samples_start = ch_wf->samples().data();

    if(ped_decay_factor_ > 0) {
      int32_t sum = 0;
      const auto* samples = samples_start + ped_window_0_;
      for(unsigned i=0; i<ped_window_n_; i++)sum += *(samples++);
      pedestal = pedestal * /* (1.0-ped_decay_factor_) */ ped_factor_A_
        + sum * /* ped_decay_factor_/double(ped_window_n_) */ ped_factor_B_ ;
    } else if(pedestal < 0) {
      int32_t sum = 0;
      const auto* samples = samples_start + ped_window_0_;
      for(unsigned i=0; i<ped_window_n_; i++)sum += *(samples++);
      pedestal = sum/double(ped_window_n_);
    }

    int32_t sum_tx = 0;
    int32_t sum_x = 0;
    const auto* samples = samples_start + sig_window_0_;
    for(unsigned i=0; i<sig_window_n_; i++) {
      sum_x += *samples;
      sum_tx += *(samples++) * i;
    }

    double numerator = sum_tx - pedestal * /* ped_window_n_*(ped_window_n_-1)/2 */ sum_t_;
    double denominator = sum_x - pedestal * sig_window_n_;

    if(denominator < 0)return 0;
    return numerator/denominator + sig_window_0_;
  }

  calin::ix::iact_data::functional_event_visitor::
    MeantimeTimingFunctionalTelescopeEventVisitorConfig config_;

  unsigned ped_window_0_ = 0;
  unsigned ped_window_n_ = 0;
  unsigned sig_window_0_ = 0;
  unsigned sig_window_n_ = 0;
  double ped_decay_factor_ = 0;

  double ped_factor_A_;
  double ped_factor_B_;
  double sum_t_;

  double high_gain_value_ = 0.0;
  double low_gain_value_ = 0.0;
  std::vector<double> high_gain_pedestal_;
  std::vector<double> low_gain_pedestal_;
};


} } } // namespace calin::iact_data::functional_event_visitor
