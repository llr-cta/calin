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

CountersEventNumberGlitchDetector::
CountersEventNumberGlitchDetector(unsigned counter_index):
  TelescopeEventVisitor(), counter_index_(counter_index)
{
  // nothing to see here
}

CountersEventNumberGlitchDetector::~CountersEventNumberGlitchDetector()
{
  // nothing to see here
}

bool CountersEventNumberGlitchDetector::visit_telescope_event(
    calin:: ix::iact_data::telescope_event::TelescopeEvent* event)
{
  if(event->module_counter_size() > counters_event_num_diff_.size())
    counters_event_num_diff_.resize(event->module_counter_size());
  int64_t index = event->source_event_index();
  bool found_glitch = false;
  if(event->local_event_number() != index+local_event_num_diff_)
  {
    local_event_num_diff_ = event->local_event_number()-index;
    found_glitch = true;
  }
  for(unsigned imodctr=0; imodctr!=event->module_counter_size(); imodctr++)
  {
    const auto& mod = event->module_counter(imodctr);
    unsigned imod = mod.module_id();
    if(mod.counter_size()<counter_index_)continue;
    if(mod.counter(counter_index_).value() !=
      index+counters_event_num_diff_[imod])
    {
      counters_event_num_diff_[imod] =
        mod.counter(counter_index_).value()-index;
      found_glitch = true;
    }
  }
  if(found_glitch)
  {
    auto* glitch = glitch_data_.add_glitch();
    glitch->set_source_event_index(index);
    glitch->set_delta_event_number(local_event_num_diff_);
    for(auto diff : counters_event_num_diff_)
      glitch->add_delta_counters_event_number(diff);
    for(auto mod_index : event->module_index())
      glitch->add_source_event_has_module(mod_index != -1);
  }
  return true;
}
