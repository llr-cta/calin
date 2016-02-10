/*

   calin/iact_data/event_dispatcher.cpp -- Stephen Fegan -- 2016-02-09

   A dispatcher of run and event data

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

#include <iact_data/event_dispatcher.hpp>

using namespace calin::iact_data::event_dispatcher;
using namespace calin::iact_data::telescope_data_source;
using namespace calin::ix::iact_data::telescope_event;

TelescopeEventDispatcher::TelescopeEventDispatcher()
{
  // nothing to see here
}

TelescopeEventDispatcher::~TelescopeEventDispatcher()
{

}

void TelescopeEventDispatcher::accept(TelescopeEvent* event)
{

}

void TelescopeEventDispatcher::
accept_from_src(TelescopeDataSource* src, unsigned num_event_max)
{
  while(TelescopeEvent* event = src->get_next())accept(event);
}
