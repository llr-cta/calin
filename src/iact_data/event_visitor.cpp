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

#include <iact_data/event_visitor.hpp>

using namespace calin::iact_data::event_visitor;
using namespace calin::ix::iact_data::telescope_event;

TelescopeEventVisitor::~TelescopeEventVisitor()
{
  // nothing to see here
}

bool TelescopeEventVisitor::visit_telescope_event(TelescopeEvent* event)
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
