/*

   calin/diagnostics/waveform.cpp -- Stephen Fegan -- 2016-04-13

   Waveform diagnostics visitor

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
#include <diagnostics/waveform.hpp>

using namespace calin::io::log;
using namespace calin::diagnostics::waveform;

WaveformStatsVisitor::WaveformStatsVisitor(bool calculate_covariance):
  TelescopeEventVisitor(), calculate_covariance_(calculate_covariance)
{
  // nothing to see here
}

WaveformStatsVisitor::~WaveformStatsVisitor()
{
  // nothing to see here
}

bool WaveformStatsVisitor::demand_waveforms()
{
  return true;
}

bool WaveformStatsVisitor::is_parallelizable()
{
  return true;
}

WaveformStatsVisitor* WaveformStatsVisitor::new_sub_visitor()
{
  auto* sub_visitor = new WaveformStatsVisitor;
  sub_visitor->parent_ = this;
  return sub_visitor;
}

bool WaveformStatsVisitor::visit_telescope_run(
  const calin::ix::iact_data::telescope_run_configuration::
    TelescopeRunConfiguration* run_config)
{
  run_config_ = run_config;
  results_.Clear();
  unsigned N = run_config->num_samples();
  for(int ichan = 0; ichan<run_config->configured_channel_id_size(); ichan++)
  {
    auto* hg_wf = results_.add_high_gain();
    hg_wf->mutable_sum()->Resize(N,0);
    hg_wf->mutable_sum_squared()->Resize(N,0);
    if(calculate_covariance_)
      hg_wf->mutable_sum_product()->Resize(N*(N-1)/2,0);
    auto* lg_wf = results_.add_low_gain();
    lg_wf->mutable_sum()->Resize(N,0);
    lg_wf->mutable_sum_squared()->Resize(N,0);
    if(calculate_covariance_)
      lg_wf->mutable_sum_product()->Resize(N*(N-1)/2,0);
  }
  return true;
}

bool WaveformStatsVisitor::leave_telescope_run()
{
  run_config_ = nullptr;
  return true;
}

bool WaveformStatsVisitor::visit_telescope_event(
  calin::ix::iact_data::telescope_event::TelescopeEvent* event)
{
  // nothing to see here
  return true;
}

bool WaveformStatsVisitor::visit_waveform(unsigned ichan,
  calin::ix::iact_data::telescope_event::ChannelWaveform* high_gain,
  calin::ix::iact_data::telescope_event::ChannelWaveform* low_gain)
{
  const int index = run_config_->configured_channel_index(ichan);
  if(high_gain)
    process_one_waveform(high_gain, results_.mutable_high_gain(index));
  if(low_gain)
    process_one_waveform(low_gain, results_.mutable_low_gain(index));
  return true;
}

void WaveformStatsVisitor::
process_one_waveform(
  const calin::ix::iact_data::telescope_event::ChannelWaveform* wf,
  ix::diagnostics::waveform::WaveformRawStats* stat)
{
  const unsigned nsample = run_config_->num_samples();
  assert(wf->samples_size() == int(nsample));
  const auto* sample = wf->samples().data();
  stat->set_num_entries(stat->num_entries()+1);
  auto* sum = stat->mutable_sum()->mutable_data();
  for(unsigned isample=0; isample<nsample; isample++)
    sum[isample] += sample[isample];
  auto* sum_squared = stat->mutable_sum_squared()->mutable_data();
  for(unsigned isample=0; isample<nsample; isample++)
    sum_squared[isample] += sample[isample] * sample[isample];
  if(calculate_covariance_)
  {
    auto* sum_product = stat->mutable_sum_product()->mutable_data();
    const auto* sample_j = sample;
    unsigned msample = nsample;
    for(unsigned isample=0; isample<nsample; isample++)
    {
      uint32_t sample_i = *sample_j;
      ++sample_j;
      --msample;
      for(unsigned jsample=0; jsample<msample; jsample++)
        sum_product[jsample] += sample_i * sample_j[jsample];
      sum_product += msample;
    }
  }
}

bool WaveformStatsVisitor::merge_results()
{
  if(parent_)
  {
    assert(results_.high_gain_size() == parent_->results_.high_gain_size());
    assert(results_.low_gain_size() == parent_->results_.low_gain_size());
    for(int ichan=0; ichan<results_.high_gain_size(); ichan++)
      merge_one_gain(&results_.high_gain(ichan),
        parent_->results_.mutable_high_gain(ichan));
    for(int ichan=0; ichan<results_.low_gain_size(); ichan++)
      merge_one_gain(&results_.low_gain(ichan),
        parent_->results_.mutable_low_gain(ichan));
  }
  return true;
}

void WaveformStatsVisitor::merge_one_gain(
  const ix::diagnostics::waveform::WaveformRawStats* from,
  ix::diagnostics::waveform::WaveformRawStats* to)
{
  assert(to->sum_size() == from->sum_size());
  assert(to->sum_squared_size() == from->sum_squared_size());
  assert(to->sum_product_size() == from->sum_product_size());
  to->set_num_entries(to->num_entries() + from->num_entries());
  for(int i=0; i<from->sum_size(); i++)
    to->set_sum(i,to->sum(i) + from->sum(i));
  for(int i=0; i<from->sum_squared_size(); i++)
    to->set_sum_squared(i,to->sum_squared(i) + from->sum_squared(i));
  for(int i=0; i<from->sum_product_size(); i++)
    to->set_sum_product(i,to->sum_product(i) + from->sum_product(i));
}
