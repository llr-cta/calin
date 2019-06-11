/*

   calin/diagnostics/delta_t.cpp -- Stephen Fegan -- 2016-06-03

   Delta-T visitor

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, CNRS, Ecole Polytechnique, Institut Polytechnique de Paris

   This file is part of "calin"

   "calin" is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 or
   later, as published by the Free Software Foundation.

   "calin" is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

*/

#include "diagnostics/delta_t.hpp"

using namespace calin::diagnostics::delta_t;

bool OneModuleCounterDeltaTVisitor::
visit_telescope_event(uint64_t seq_index,
  calin::ix::iact_data::telescope_event::TelescopeEvent* event)
{
  int64_t last_counter_event_number = last_counter_event_number_;
  last_counter_event_number_ = -1;
  has_delta_t_val_ = false;
  if(module_id_ >= event->module_index_size())return true;
  auto module_index = event->module_index(module_id_);
  if(module_index == -1)return true;
  const auto& counters = event->module_counter(module_index);
  if(counter_index_ >= counters.counter_value_size())return true;
  int64_t counter_val = counters.counter_value(counter_index_);
  if(last_counter_event_number != -1 and
    event->local_event_number() == last_counter_event_number+1)
  {
    has_delta_t_val_ = true;
    delta_t_val_ = (counter_val - last_counter_val_)*conversion_to_ns_;
  }
  last_counter_event_number_ = event->local_event_number();
  last_counter_val_ = counter_val;
  return true;
}

bool OneModuleCounterDeltaTVisitor::get_value(double& value)
{
  if(!has_delta_t_val_)return false;
  value = delta_t_val_;
  return true;
}
