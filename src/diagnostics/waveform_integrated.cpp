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

IntegratedWaveformStatsVisitor::
IntegratedWaveformStatsVisitor(bool calculate_covariance):
  TelescopeEventVisitor(), calculate_covariance_(calculate_covariance)
{
  // nothing to see here
}

IntegratedWaveformStatsVisitor::~IntegratedWaveformStatsVisitor()
{
  // nothing to see here
}

bool IntegratedWaveformStatsVisitor::demand_waveforms()
{
  return true;
}

bool IntegratedWaveformStatsVisitor::is_parallelizable()
{
  return true;
}

IntegratedWaveformStatsVisitor* IntegratedWaveformStatsVisitor::new_sub_visitor()
{
  auto* sub_visitor = new IntegratedWaveformStatsVisitor;
  sub_visitor->parent_ = this;
  return sub_visitor;
}

bool IntegratedWaveformStatsVisitor::visit_telescope_run(
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

bool IntegratedWaveformStatsVisitor::leave_telescope_run()
{
  run_config_ = nullptr;
  return true;
}

bool IntegratedWaveformStatsVisitor::visit_telescope_event(
  calin::ix::iact_data::telescope_event::TelescopeEvent* event)
{
  for(auto& ichan : high_gain_mask_)ichan = 0;
  for(auto& ichan : low_gain_mask_)ichan = 0;
  return true;
}

bool IntegratedWaveformStatsVisitor::visit_waveform(unsigned ichan,
  calin::ix::iact_data::telescope_event::ChannelWaveform* high_gain,
  calin::ix::iact_data::telescope_event::ChannelWaveform* low_gain)
{
  const int index = run_config_->configured_channel_index(ichan);
#if 0
  if(high_gain)
    process_one_waveform(high_gain, results_.mutable_high_gain(index));
  if(low_gain)
    process_one_waveform(low_gain, results_.mutable_low_gain(index));
#endif
  return true;
}

#if 0
void IntegratedWaveformStatsVisitor::
process_one_waveform(
  const calin::ix::iact_data::telescope_event::ChannelWaveform* wf,
  ix::diagnostics::waveform::WaveformRawPSD* psd)
{
}
#endif

bool IntegratedWaveformStatsVisitor::merge_results()
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

#if 0
void IntegratedWaveformStatsVisitor::merge_one_gain(
  const ix::diagnostics::waveform::WaveformRawPSD* from,
  ix::diagnostics::waveform::WaveformRawPSD* to)
{
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
}
#endif
