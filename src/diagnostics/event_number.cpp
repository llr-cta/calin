/*

   calin/diagnostics/event_number.cpp -- Stephen Fegan -- 2016-03-09

   Event number diagnostics visitor

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
#include <diagnostics/event_number.hpp>

using namespace calin::io::log;
using namespace calin::diagnostics::event_number;
using namespace calin::ix::diagnostics::event_number;

SequentialEventNumberGlitchDetector::
SequentialEventNumberGlitchDetector(bool test_local_event_number):
  TelescopeEventVisitor(), test_local_event_number_(test_local_event_number)
{
  // nothing to see here
}

SequentialEventNumberGlitchDetector::~SequentialEventNumberGlitchDetector()
{
  // nothing to see here
}

bool SequentialEventNumberGlitchDetector::visit_telescope_event(
  calin::ix::iact_data::telescope_event::TelescopeEvent* event)
{
  int64_t event_number = test_local_event_number_ ?
    event->local_event_number() : event->array_event_number();
  if(event_number != last_event_number_+1)
  {
    auto* glitch = glitch_data_.add_glitch();
    glitch->set_source_event_index(event->source_event_index());
    glitch->set_event_number(event_number);
    glitch->set_last_event_number(last_event_number_);
  }
  last_event_number_ = event_number;
  return true;
}
