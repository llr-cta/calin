/*

   calin/diagnostics/module.cpp -- Stephen Fegan -- 2016-03-04

   Module diagnostics visitor

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
#include <diagnostics/module.hpp>

using namespace calin::io::log;
using namespace calin::diagnostics::module;

ModulePresentVisitor::ModulePresentVisitor(): TelescopeEventVisitor()
{
  // nothing to see here
}

ModulePresentVisitor::~ModulePresentVisitor()
{
  // nothing to see here
}

bool ModulePresentVisitor::demand_waveforms()
{
  return false;
}

bool ModulePresentVisitor::visit_telescope_run(
  const calin::ix::iact_data::telescope_run_configuration::
    TelescopeRunConfiguration* run_config)
{
  return true;
}

bool ModulePresentVisitor::visit_telescope_event(uint64_t seq_index,
    calin:: ix::iact_data::telescope_event::TelescopeEvent* event)
{
  auto index = event->source_event_index();
  module_data_.set_num_events(module_data_.num_events()+1);
  if(event->all_modules_present())return true;
  while(module_data_.module_missing_num_events_size() <
      event->module_index_size())
  {
    module_data_.add_module_missing_num_events(0);
    module_data_.add_module_missing_history();
  }
  for(int imod=0;imod<event->module_index_size();imod++)
    if(event->module_index(imod) == -1)
    {
      auto* hist = module_data_.mutable_module_missing_history(imod);
      if(hist->range_size()==0 or
        hist->range(hist->range_size()-1).last_source_index() != index)
      {
        auto* range = hist->add_range();
        range->set_first_source_index(index);
        range->set_last_source_index(index+1);
      }
      else
      {
        hist->mutable_range(hist->range_size()-1)->
          set_last_source_index(index+1);
      }
      module_data_.set_module_missing_num_events(imod,
        module_data_.module_missing_num_events(imod)+1);
    }
  return true;
}
