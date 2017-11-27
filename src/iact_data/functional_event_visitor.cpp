/*

   calin/iact_data/functional_event_visitor.cpp -- Stephen Fegan -- 2016-04-13

   Vistor to event data that calculates some quantity from the waveforms

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <iact_data/functional_event_visitor.hpp>
#include <util/log.hpp>

using namespace calin::util::log;
using namespace calin::iact_data::functional_event_visitor;

FixedWindowSumFunctionalTelescopeEventVisitor::
FixedWindowSumFunctionalTelescopeEventVisitor(
  calin::ix::iact_data::functional_event_visitor::
    FixedWindowSumFunctionalTelescopeEventVisitorConfig config):
  DualGainInt32FunctionalTelescopeEventVisitor(), config_(config)
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
      or config_.integration_0() > int(run_config->num_samples()-window_n_))
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
    DualGainInt32FunctionalTelescopeEventVisitor* sig_value_supplier,
    DualGainInt32FunctionalTelescopeEventVisitor* bkg_value_supplier):
  DualGainInt32FunctionalTelescopeEventVisitor(),
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
    dynamic_cast<DualGainInt32FunctionalTelescopeEventVisitor*>(
      i_sig_value_supplier->second);
  auto i_bkg_value_supplier = antecedent_visitors.find(bkg_value_supplier_);
  if(i_bkg_value_supplier == antecedent_visitors.end())
    throw std::logic_error("No thread specific value suppler for background!");
  auto* bkg_value_supplier =
    dynamic_cast<DualGainInt32FunctionalTelescopeEventVisitor*>(
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

RisetimeTimingFunctionalTelescopeEventVisitor::
RisetimeTimingFunctionalTelescopeEventVisitor():
  DualGainDoubleFunctionalTelescopeEventVisitor()
{
  // nothing to see here
}

RisetimeTimingFunctionalTelescopeEventVisitor::
~RisetimeTimingFunctionalTelescopeEventVisitor()
{
  // nothing to see here
}

bool RisetimeTimingFunctionalTelescopeEventVisitor::demand_waveforms()
{
  return true;
}

bool RisetimeTimingFunctionalTelescopeEventVisitor::is_parallelizable()
{
  return true;
}

RisetimeTimingFunctionalTelescopeEventVisitor*
RisetimeTimingFunctionalTelescopeEventVisitor::new_sub_visitor(
  const std::map<calin::iact_data::event_visitor::TelescopeEventVisitor*,
      calin::iact_data::event_visitor::TelescopeEventVisitor*>&
    antecedent_visitors)
{
  return new RisetimeTimingFunctionalTelescopeEventVisitor;
}

bool RisetimeTimingFunctionalTelescopeEventVisitor::visit_waveform(
  unsigned ichan,
  calin::ix::iact_data::telescope_event::ChannelWaveform* high_gain,
  calin::ix::iact_data::telescope_event::ChannelWaveform* low_gain)
{
  if(high_gain)
    high_gain_value_ = process_one_gain(high_gain);
  if(low_gain)
    low_gain_value_ = process_one_gain(low_gain);
  return true;
}

double RisetimeTimingFunctionalTelescopeEventVisitor::low_gain_value()
{
  return low_gain_value_;
}

double RisetimeTimingFunctionalTelescopeEventVisitor::high_gain_value()
{
  return high_gain_value_;
}

MeantimeTimingFunctionalTelescopeEventVisitor::
MeantimeTimingFunctionalTelescopeEventVisitor(
  calin::ix::iact_data::functional_event_visitor::
    MeantimeTimingFunctionalTelescopeEventVisitorConfig config):
  DualGainDoubleFunctionalTelescopeEventVisitor(), config_(config)
{
  // nothing to see here
}

MeantimeTimingFunctionalTelescopeEventVisitor::
~MeantimeTimingFunctionalTelescopeEventVisitor()
{
  // nothing to see here
}

bool MeantimeTimingFunctionalTelescopeEventVisitor::demand_waveforms()
{
  return true;
}

bool MeantimeTimingFunctionalTelescopeEventVisitor::is_parallelizable()
{
  return true;
}

MeantimeTimingFunctionalTelescopeEventVisitor*
MeantimeTimingFunctionalTelescopeEventVisitor::new_sub_visitor(
  const std::map<calin::iact_data::event_visitor::TelescopeEventVisitor*,
      calin::iact_data::event_visitor::TelescopeEventVisitor*>&
    antecedent_visitors)
{
  return new MeantimeTimingFunctionalTelescopeEventVisitor(config_);
}

bool MeantimeTimingFunctionalTelescopeEventVisitor::visit_telescope_run(
  const calin::ix::iact_data::telescope_run_configuration::
    TelescopeRunConfiguration* run_config)
{
  if(config_.pedestal_n()==0)ped_window_n_ = run_config->num_samples();
  else ped_window_n_ = config_.pedestal_n();
  if(ped_window_n_ > run_config->num_samples())
    throw std::out_of_range("MeantimeTimingFunctionalTelescopeEventVisitor: "
      "requested pedestal window larger than number of samples: "
      + std::to_string(ped_window_n_) + ">"
      + std::to_string(run_config->num_samples()));

  if(config_.pedestal_0() < -int(run_config->num_samples()-ped_window_n_)
      or config_.pedestal_0() > int(run_config->num_samples()-ped_window_n_))
    throw std::out_of_range("MeantimeTimingFunctionalTelescopeEventVisitor: "
      "requested pedestal window start "
      + std::to_string(config_.pedestal_0())
      + " out of allowed range: ["
      + std::to_string(-int(run_config->num_samples()-ped_window_n_)) + ", "
      + std::to_string(run_config->num_samples()-ped_window_n_) + "]");

  if(config_.pedestal_0()<0)
    ped_window_0_ = config_.pedestal_0() + run_config->num_samples()
      - ped_window_n_;
  else
    ped_window_0_ = config_.pedestal_0();

  if(config_.signal_n()==0)sig_window_n_ = run_config->num_samples();
  else sig_window_n_ = config_.signal_n();
  if(sig_window_n_ > run_config->num_samples())
    throw std::out_of_range("MeantimeTimingFunctionalTelescopeEventVisitor: "
      "requested signal window larger than number of samples: "
      + std::to_string(sig_window_n_) + ">"
      + std::to_string(run_config->num_samples()));

  if(config_.signal_0() < -int(run_config->num_samples()-sig_window_n_)
      or config_.signal_0() > int(run_config->num_samples()-sig_window_n_))
    throw std::out_of_range("MeantimeTimingFunctionalTelescopeEventVisitor: "
      "requested signal window start "
      + std::to_string(config_.signal_0())
      + " out of allowed range: ["
      + std::to_string(-int(run_config->num_samples()-sig_window_n_)) + ", "
      + std::to_string(run_config->num_samples()-sig_window_n_) + "]");

  if(config_.signal_0()<0)
    sig_window_0_ = config_.signal_0() + run_config->num_samples()
      - sig_window_n_;
  else
    sig_window_0_ = config_.signal_0();

  ped_decay_factor_ = config_.pedestal_decay_constant();

  ped_factor_A_ = 1.0-ped_decay_factor_;
  ped_factor_B_ = ped_decay_factor_/double(ped_window_n_);
  sum_t_ = sig_window_n_*(sig_window_n_-1)/2;

  high_gain_pedestal_.clear();
  for(const auto& ped : config_.high_gain_pedestal())
    if(high_gain_pedestal_.size() < run_config->configured_channel_id_size())
      high_gain_pedestal_.push_back(ped);
  while(high_gain_pedestal_.size() < run_config->configured_channel_id_size())
    high_gain_pedestal_.push_back(-1);

  low_gain_pedestal_.clear();
  for(const auto& ped : config_.low_gain_pedestal())
    if(low_gain_pedestal_.size() < run_config->configured_channel_id_size())
      low_gain_pedestal_.push_back(ped);
  while(low_gain_pedestal_.size() < run_config->configured_channel_id_size())
    low_gain_pedestal_.push_back(-1);

  return true;
}

bool MeantimeTimingFunctionalTelescopeEventVisitor::visit_waveform(
  unsigned ichan,
  calin::ix::iact_data::telescope_event::ChannelWaveform* high_gain,
  calin::ix::iact_data::telescope_event::ChannelWaveform* low_gain)
{
  if(high_gain)
    high_gain_value_ = process_one_gain(high_gain_pedestal_[ichan], high_gain);
  static unsigned nprint = 10;
  if(ichan == 0 and nprint) {
    LOG(INFO) << ichan << ' '
      << high_gain_pedestal_[ichan] << ' ' << high_gain_value_ << ' '
      << ped_window_0_ << ' ' << ped_window_n_ << ' '
      << sig_window_0_ << ' ' << sig_window_n_ << ' '
      << ped_decay_factor_ << ' '
      << ped_factor_A_ << ' ' << ped_factor_B_ << ' ' << sum_t_;

    const auto*const samples_start = high_gain->samples().data();
    auto xlog = LOG(INFO);
    xlog << "->";
    for(unsigned isamp=0;isamp<60;isamp++)
      xlog << " " << samples_start[isamp];
    --nprint;
  }
  if(low_gain)
    low_gain_value_ = process_one_gain(low_gain_pedestal_[ichan], low_gain);
  return true;
}

double MeantimeTimingFunctionalTelescopeEventVisitor::low_gain_value()
{
  return low_gain_value_;
}

double MeantimeTimingFunctionalTelescopeEventVisitor::high_gain_value()
{
  return high_gain_value_;
}
