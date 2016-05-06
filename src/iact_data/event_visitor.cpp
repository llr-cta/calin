/*

   calin/iact_data/event_visitor.cpp -- Stephen Fegan -- 2016-02-10

   A visitor of run and event data

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

#include <stdexcept>

#include <iact_data/event_visitor.hpp>

using namespace calin::iact_data::event_visitor;
using namespace calin::ix::iact_data::telescope_event;
using namespace calin::ix::iact_data::telescope_run_configuration;

TelescopeEventVisitor::~TelescopeEventVisitor()
{
  // nothing to see here
}

bool TelescopeEventVisitor::demand_waveforms()
{
  return true;
}

bool TelescopeEventVisitor::is_parallelizable()
{
  return false;
}

TelescopeEventVisitor* TelescopeEventVisitor::new_sub_visitor(
  const std::map<TelescopeEventVisitor*,TelescopeEventVisitor*>&
    antecedent_visitors)
{
  throw std::runtime_error("TelescopeEventVisitor: parallel sub workers not "
    "implemented.");
  return nullptr;
}

bool TelescopeEventVisitor::
visit_telescope_run(const TelescopeRunConfiguration* run_config)
{
  return true;
}

bool TelescopeEventVisitor::leave_telescope_run()
{
  return true;
}

bool TelescopeEventVisitor::
visit_telescope_event(uint64_t seq_index, TelescopeEvent* event)
{
  return true;
}

bool TelescopeEventVisitor::leave_telescope_event()
{
  return true;
}

bool TelescopeEventVisitor::visit_waveform(unsigned ichan,
    ChannelWaveform* high_gain, ChannelWaveform* low_gain)
{
  return true;
}

bool TelescopeEventVisitor::merge_results()
{
  return true;
}
