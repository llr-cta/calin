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
  for(int ichan = 0; ichan < run_config->configured_channel_id_size(); ichan++)
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
}

bool WaveformStatsVisitor::visit_waveform(unsigned ichan,
  calin::ix::iact_data::telescope_event::ChannelWaveform* high_gain,
  calin::ix::iact_data::telescope_event::ChannelWaveform* low_gain)
{

}

bool WaveformStatsVisitor::merge_results()
{

}
