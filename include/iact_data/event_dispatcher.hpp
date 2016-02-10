/*

   calin/iact_data/event_dispatcher.hpp -- Stephen Fegan -- 2016-01-30

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

#pragma once

#include <string>
#include <vector>

#include <calin_global_definitions.hpp>
#include <calin_global_config.hpp>
#include <iact_data/telescope_data_source.hpp>
#include <iact_data/event_visitor.hpp>

namespace calin { namespace iact_data { namespace event_dispatcher {

class TelescopeEventDispatcher
{
public:
  TelescopeEventDispatcher();
  ~TelescopeEventDispatcher();
  void add_visitor(
    calin::iact_data::event_visitor::TelescopeEventVisitor* visitor,
    bool adopt_visitor = false);
  void accept(
    calin::ix::iact_data::telescope_event::TelescopeEvent* event);
  void accept_from_src(
    calin::iact_data::telescope_data_source::TelescopeDataSource* src,
    unsigned num_event_max = 0);
private:
  std::vector<calin::iact_data::event_visitor::TelescopeEventVisitor*>
    adopted_visitors_;
  std::vector<calin::iact_data::event_visitor::TelescopeEventVisitor*>
    visitors_;
};

} } } // namespace calin::iact_data::event_dispatcher
