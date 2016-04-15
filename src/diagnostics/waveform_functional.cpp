/*

   calin/diagnostics/waveform_integrated.cpp -- Stephen Fegan -- 2016-05-10

   Waveform diagnostics visitor - integrated stats with histo and
   cross-channel covariance

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

#include <io/log.hpp>
#include <math/special.hpp>
#include <diagnostics/waveform.hpp>

using calin::math::special::SQR;
using namespace calin::diagnostics::waveform;
using calin::iact_data::functional_event_visitor::
  DualValueInt32FunctionalTelescopeEventVisitor;
using calin::ix::diagnostics::waveform::
  FunctionalWaveformStatsVisitorConfig;

FunctionalWaveformStatsVisitor::
FunctionalWaveformStatsVisitor(
    DualValueInt32FunctionalTelescopeEventVisitor* value_supplier,
    const FunctionalWaveformStatsVisitorConfig& config):
  TelescopeEventVisitor(), value_supplier_(value_supplier), config_(config)
{
  // nothing to see here
}

FunctionalWaveformStatsVisitor::~FunctionalWaveformStatsVisitor()
{
  // nothing to see here
}

bool FunctionalWaveformStatsVisitor::demand_waveforms()
{
  return true;
}

bool FunctionalWaveformStatsVisitor::is_parallelizable()
{
  return true;
}

FunctionalWaveformStatsVisitor*
FunctionalWaveformStatsVisitor::new_sub_visitor(
  const std::map<TelescopeEventVisitor*,TelescopeEventVisitor*>&
    antecedent_visitors)
{
  auto i_sub_value_supplier = antecedent_visitors.find(value_supplier_);
  if(i_sub_value_supplier == antecedent_visitors.end())
    throw std::logic_error("No thread specific value suppler!");
  auto* sub_value_supplier =
    dynamic_cast<DualValueInt32FunctionalTelescopeEventVisitor*>(
      i_sub_value_supplier->second);
  auto* sub_visitor =
    new FunctionalWaveformStatsVisitor(sub_value_supplier, config_);
  sub_visitor->parent_ = this;
  return sub_visitor;
}

bool FunctionalWaveformStatsVisitor::visit_telescope_run(
  const calin::ix::iact_data::telescope_run_configuration::
    TelescopeRunConfiguration* run_config)
{
  run_config_ = run_config;

  results_.Clear();

  int nchan = run_config_->configured_channel_id_size();
  for(auto chan_id : run_config_->configured_channel_id())
    results_.add_channel_id(chan_id);

  high_gain_mask_.resize(nchan);
  high_gain_signal_.resize(nchan);
  auto* hg_results = results_.mutable_high_gain();
  hg_results->mutable_num_sum_entries()->Resize(nchan,0);
  hg_results->mutable_num_sum_product_entries()->Resize(nchan*(nchan-1)/2,0);
  hg_results->mutable_sum()->Resize(nchan,0);
  hg_results->mutable_sum_squared()->Resize(nchan,0);
  hg_results->mutable_sum_product()->Resize(nchan*(nchan-1)/2,0);

  low_gain_mask_.resize(nchan);
  low_gain_signal_.resize(nchan);
  auto* lg_results = results_.mutable_low_gain();
  lg_results->mutable_num_sum_entries()->Resize(nchan,0);
  lg_results->mutable_num_sum_product_entries()->Resize(nchan*(nchan-1)/2,0);
  lg_results->mutable_sum()->Resize(nchan,0);
  lg_results->mutable_sum_squared()->Resize(nchan,0);
  lg_results->mutable_sum_product()->Resize(nchan*(nchan-1)/2,0);

  return true;
}

bool FunctionalWaveformStatsVisitor::leave_telescope_run()
{
  run_config_ = nullptr;
  return true;
}

bool FunctionalWaveformStatsVisitor::visit_telescope_event(
  calin::ix::iact_data::telescope_event::TelescopeEvent* event)
{
  for(auto& ichan : high_gain_mask_)ichan = 0;
  for(auto& ichan : low_gain_mask_)ichan = 0;
  return true;
}

bool FunctionalWaveformStatsVisitor::leave_telescope_event()
{
  if(results_.has_high_gain())
    process_one_gain(high_gain_mask_, high_gain_signal_,
      results_.mutable_high_gain());
  if(results_.has_low_gain())
    process_one_gain(low_gain_mask_, low_gain_signal_,
      results_.mutable_low_gain());
  return true;
}

bool FunctionalWaveformStatsVisitor::visit_waveform(unsigned ichan,
  calin::ix::iact_data::telescope_event::ChannelWaveform* high_gain,
  calin::ix::iact_data::telescope_event::ChannelWaveform* low_gain)
{
  const int index = run_config_->configured_channel_index(ichan);
  if(high_gain) {
    high_gain_mask_[index] = 1;
    high_gain_signal_[index] = value_supplier_->high_gain_value();
  }
  if(low_gain) {
    low_gain_mask_[index] = 1;
    low_gain_signal_[index] = value_supplier_->low_gain_value();
  }
  return true;
}

void FunctionalWaveformStatsVisitor::visit_one_waveform(
  const calin::ix::iact_data::telescope_event::ChannelWaveform* wf,
  unsigned index, std::vector<int>& mask, std::vector<int32_t>& signal)
{
  mask[index] = 1;
}

void FunctionalWaveformStatsVisitor::process_one_gain(
  const std::vector<int>& mask,
  const std::vector<int32_t>& signal,
  calin::ix::diagnostics::waveform::OneGainIntFunctionalWaveformRawStats* stats)
{
  auto* num_sum = stats->mutable_num_sum_entries()->mutable_data();
  auto* sum = stats->mutable_sum()->mutable_data();
  auto* sum_squared = stats->mutable_sum_squared()->mutable_data();
  const unsigned nchan = mask.size();
  for(unsigned ichan=0; ichan<nchan; ichan++)
    if(mask[ichan])
    {
      auto si = signal[ichan];
      num_sum[ichan]++;
      sum[ichan] += si;
      sum_squared[ichan] += SQR(si);
    }

  if(config_.calculate_covariance())
  {
    auto* num_sum_product =
      stats->mutable_num_sum_product_entries()->mutable_data();
    auto* sum_product = stats->mutable_sum_product()->mutable_data();
    for(unsigned ichan=0; ichan<nchan; ichan++)
      if(mask[ichan])
      {
        auto si = signal[ichan];
        for(unsigned jchan=ichan+1; jchan<nchan; jchan++)
          if(mask[jchan])
          {
            num_sum_product[ichan]++;
            sum_product[ichan] += si*signal[jchan];
          }
      }
    }
}

bool FunctionalWaveformStatsVisitor::merge_results()
{
  if(parent_)
  {
#if 0
    assert(results_.high_gain_size() == parent_->results_.high_gain_size());
    assert(results_.low_gain_size() == parent_->results_.low_gain_size());
    for(int ichan=0; ichan<results_.high_gain_size(); ichan++)
      merge_one_gain(&results_.high_gain(ichan),
        parent_->results_.mutable_high_gain(ichan));
    for(int ichan=0; ichan<results_.low_gain_size(); ichan++)
      merge_one_gain(&results_.low_gain(ichan),
        parent_->results_.mutable_low_gain(ichan));
#endif
  }
  return true;
}

void FunctionalWaveformStatsVisitor::merge_one_gain(
  const ix::diagnostics::waveform::OneGainIntFunctionalWaveformRawStats* from,
  ix::diagnostics::waveform::OneGainIntFunctionalWaveformRawStats* to)
{
#if 0
  assert(to->psd_sum_size() == from->psd_sum_size());
  assert(to->psd_sum_squared_size() == from->psd_sum_squared_size());
  to->set_num_entries(to->num_entries() + from->num_entries());
  for(int i=0; i<from->psd_sum_size(); i++)
    to->set_psd_sum(i,to->psd_sum(i) + from->psd_sum(i));
  for(int i=0; i<from->psd_sum_squared_size(); i++)
    to->set_psd_sum_squared(i,to->psd_sum_squared(i) + from->psd_sum_squared(i));
  for(int i=0; i<from->corr_sum_size(); i++)
    to->set_corr_sum(i,to->corr_sum(i) + from->corr_sum(i));
  for(int i=0; i<from->corr_sum_squared_size(); i++)
    to->set_corr_sum_squared(i,to->corr_sum_squared(i) + from->corr_sum_squared(i));
#endif
}

Eigen::MatrixXd FunctionalWaveformStatsVisitor::camera_cov(
  const ix::diagnostics::waveform::OneGainIntFunctionalWaveformRawStats* stat)
{

}
