/*

   calin/iact_data/functional_event_visitor.cpp -- Stephen Fegan -- 2016-04-13

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

#include <iact_data/functional_event_visitor.hpp>

using namespace calin::iact_data::functional_event_visitor;

DualValueInt32FunctionalTelescopeEventVisitor::
~DualValueInt32FunctionalTelescopeEventVisitor()
{
  // nothing to see here
}

FixedWindowSumFunctionalTelescopeEventVisitor::
FixedWindowSumFunctionalTelescopeEventVisitor(
  calin::ix::iact_data::functional_event_visitor::
    FixedWindowSumFunctionalTelescopeEventVisitorConfig config):
  DualValueInt32FunctionalTelescopeEventVisitor(), config_(config)
{
  // nothing to see here
}

FixedWindowSumFunctionalTelescopeEventVisitor::
~FixedWindowSumFunctionalTelescopeEventVisitor()
{
  // nothing to see here
}

bool FixedWindowSumFunctionalTelescopeEventVisitor::demand_waveforms()
{
  return true;
}

bool FixedWindowSumFunctionalTelescopeEventVisitor::is_parallelizable()
{
  return true;
}

FixedWindowSumFunctionalTelescopeEventVisitor*
FixedWindowSumFunctionalTelescopeEventVisitor::new_sub_visitor(
  const std::map<TelescopeEventVisitor*,TelescopeEventVisitor*>&
    antecedent_visitors)
{
  return new FixedWindowSumFunctionalTelescopeEventVisitor(config_);
}

bool FixedWindowSumFunctionalTelescopeEventVisitor::visit_telescope_run(
  const calin::ix::iact_data::telescope_run_configuration::
    TelescopeRunConfiguration* run_config)
{
  if(config_.integration_n()==0)window_n_ = run_config->num_samples();
  else window_n_ = config_.integration_n();
  if(window_n_ > run_config->num_samples())
    throw std::out_of_range("FixedWindowSumFunctionalTelescopeEventVisitor: "
      "requested window larger than number of samples: "
      + std::to_string(window_n_) + ">"
      + std::to_string(run_config->num_samples()));

  if(config_.integration_0() < -int(run_config->num_samples()-window_n_)
      or config_.integration_0() > run_config->num_samples()-window_n_)
    throw std::out_of_range("FixedWindowSumFunctionalTelescopeEventVisitor: "
      "requested window start "
      + std::to_string(config_.integration_0())
      + " out of allowed range: ["
      + std::to_string(-int(run_config->num_samples()-window_n_)) + ", "
      + std::to_string(run_config->num_samples()-window_n_) + "]");

  if(config_.integration_0()<0)
    window_0_ = config_.integration_0() + run_config->num_samples()
      - window_n_;
  else
    window_0_ = config_.integration_0();

  return true;
}

bool FixedWindowSumFunctionalTelescopeEventVisitor::
visit_waveform(unsigned ichan,
  calin::ix::iact_data::telescope_event::ChannelWaveform* high_gain,
  calin::ix::iact_data::telescope_event::ChannelWaveform* low_gain)
{
  if(high_gain)
    high_gain_value_ = process_one_gain(high_gain);
  if(low_gain)
    low_gain_value_ = process_one_gain(low_gain);
  return true;
}

int32_t FixedWindowSumFunctionalTelescopeEventVisitor::low_gain_value()
{
  return low_gain_value_;
}

int32_t FixedWindowSumFunctionalTelescopeEventVisitor::high_gain_value()
{
  return high_gain_value_;
}


DifferencingFunctionalTelescopeEventVisitor::
DifferencingFunctionalTelescopeEventVisitor(
    DualValueInt32FunctionalTelescopeEventVisitor* sig_value_supplier,
    DualValueInt32FunctionalTelescopeEventVisitor* bkg_value_supplier):
  DualValueInt32FunctionalTelescopeEventVisitor(),
  sig_value_supplier_(sig_value_supplier), bkg_value_supplier_(bkg_value_supplier)
{
  // nothing to see here
}

DifferencingFunctionalTelescopeEventVisitor::
~DifferencingFunctionalTelescopeEventVisitor()
{
  // nothing to see here
}

bool DifferencingFunctionalTelescopeEventVisitor::demand_waveforms()
{
  return false;
}

bool DifferencingFunctionalTelescopeEventVisitor::is_parallelizable()
{
  return true;
}

DifferencingFunctionalTelescopeEventVisitor*
DifferencingFunctionalTelescopeEventVisitor::new_sub_visitor(
  const std::map<calin::iact_data::event_visitor::TelescopeEventVisitor*,
    calin::iact_data::event_visitor::TelescopeEventVisitor*>&
  antecedent_visitors)
{
  auto i_sig_value_supplier = antecedent_visitors.find(sig_value_supplier_);
  if(i_sig_value_supplier == antecedent_visitors.end())
    throw std::logic_error("No thread specific value suppler for signal!");
  auto* sig_value_supplier =
    dynamic_cast<DualValueInt32FunctionalTelescopeEventVisitor*>(
      i_sig_value_supplier->second);
  auto i_bkg_value_supplier = antecedent_visitors.find(bkg_value_supplier_);
  if(i_bkg_value_supplier == antecedent_visitors.end())
    throw std::logic_error("No thread specific value suppler for background!");
  auto* bkg_value_supplier =
    dynamic_cast<DualValueInt32FunctionalTelescopeEventVisitor*>(
      i_bkg_value_supplier->second);
  auto* sub_visitor =
    new DifferencingFunctionalTelescopeEventVisitor(sig_value_supplier,
      bkg_value_supplier);
  return sub_visitor;
}

int32_t DifferencingFunctionalTelescopeEventVisitor::low_gain_value()
{
  return sig_value_supplier_->low_gain_value() -
    bkg_value_supplier_->low_gain_value();
}

int32_t DifferencingFunctionalTelescopeEventVisitor::high_gain_value()
{
  return sig_value_supplier_->high_gain_value() -
    bkg_value_supplier_->high_gain_value();
}
