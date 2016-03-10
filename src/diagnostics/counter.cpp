/*

   calin/diagnostics/counter.cpp -- Stephen Fegan -- 2016-03-04

   Counter diagnostics visitor

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
#include <diagnostics/counter.hpp>

using namespace calin::io::log;
using namespace calin::diagnostics::counter;

CounterDiagnostics::CounterDiagnostics(): TelescopeEventVisitor()
{
  // nothing to see here
}

CounterDiagnostics::~CounterDiagnostics()
{
  // nothing to see here
}

bool CounterDiagnostics::visit_telescope_event(
  calin::ix::iact_data::telescope_event::TelescopeEvent* event)
{
  event_data_.set_num_events(event_data_.num_events()+1);
  while(next_event_number_ < event->local_event_number())
    event_data_.add_missing_event_number(next_event_number_++);
  next_event_number_++;
  return true;
}
