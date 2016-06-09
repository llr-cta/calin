/*

   calin/iact_data/functional_event_visitor.hpp -- Stephen Fegan -- 2016-04-13

   Vistor to event data that calculates some quantity from the waveforms

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
#include <iact_data/functional_event_visitor.pb.h>

namespace calin { namespace iact_data { namespace functional_event_visitor {

class DualGainInt32FunctionalTelescopeEventVisitor:
  public calin::iact_data::event_visitor::TelescopeEventVisitor
{
public:
  CALIN_TYPEALIAS(value_type, int32_t);
  virtual ~DualGainInt32FunctionalTelescopeEventVisitor();
  virtual int32_t low_gain_value() = 0;
  virtual int32_t high_gain_value() = 0;
};

class DualGainDoubleFunctionalTelescopeEventVisitor:
  public calin::iact_data::event_visitor::TelescopeEventVisitor
{
public:
  CALIN_TYPEALIAS(value_type, double);
  virtual ~DualGainDoubleFunctionalTelescopeEventVisitor();
  virtual double low_gain_value() = 0;
  virtual double high_gain_value() = 0;
};

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

class NoPedestalTimingFunctionalTelescopeEventVisitor:
  public DualGainDoubleFunctionalTelescopeEventVisitor
{
public:
  NoPedestalTimingFunctionalTelescopeEventVisitor();
  virtual ~NoPedestalTimingFunctionalTelescopeEventVisitor();

  bool demand_waveforms() override;
  bool is_parallelizable() override;

  NoPedestalTimingFunctionalTelescopeEventVisitor* new_sub_visitor(
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
    const auto* samples_end = ch_wf->samples().data() + ch_wf->samples_size();
    const auto* samples_max = ch_wf->samples().data();
    for(const auto* samples = ch_wf->samples().data() + 1;
        samples < samples_end; samples++)
      if(*samples > *samples_max)samples_max = samples;
    const auto* samples_min = ch_wf->samples().data();
    samples_end = samples_max;
    for(const auto* samples = ch_wf->samples().data() + 1;
        samples < samples_end; samples++)
      if(*samples < *samples_min)samples_min = samples;
    return 0;
  }

  double high_gain_value_ = 0.0;
  double low_gain_value_ = 0.0;
};


} } } // namespace calin::iact_data::functional_event_visitor
