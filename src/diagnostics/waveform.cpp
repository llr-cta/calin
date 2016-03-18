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

WaveformStatsVisitor::WaveformStatsVisitor()
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

TelescopeEventVisitor* new_sub_visitor()
{
  auto* sub_visitor = new TelescopeEventVisitor;
  sub_visitor_->parent_ = this;
  sub_visitor_->wf_results_.CopyFrom(wf_results_);
  return sub_visitor_;
}

bool visit_telescope_run(
  const calin::ix::iact_data::telescope_run_configuration::
    TelescopeRunConfiguration* run_config)
{
  wf_results_.Clear();
  //wf_results_.
  return true;
}

bool leave_telescope_run()
{
  return true;
}

bool WaveformStatsVisitor::visit_telescope_event(
  calin::ix::iact_data::telescope_event::TelescopeEvent* event)
{

}

bool WaveformStatsVisitor::visit_waveform(unsigned ichan,
  calin::ix::iact_data::telescope_event::ChannelWaveform* high_gain,
  calin::ix::iact_data::telescope_event::ChannelWaveform* low_gain)
{

}

bool WaveformStatsVisitor::merge_results()
{

}
