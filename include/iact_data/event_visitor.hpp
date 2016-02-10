/*

   calin/iact_data/event_visitor.hpp -- Stephen Fegan -- 2016-01-30

   Vistor to event data

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

#include <calin_global_definitions.hpp>
#include <calin_global_config.hpp>
#include <iact_data/telescope_event.pb.h>

namespace calin { namespace iact_data { namespace event_visitor {

class TelescopeEventVisitor
{
public:
  virtual ~TelescopeEventVisitor();

  virtual bool visit_telescope_event(
    calin::ix::iact_data::telescope_event::TelescopeEvent* event);
  virtual bool leave_telescope_event();

  virtual bool visit_waveform(unsigned ichan,
    calin::ix::iact_data::telescope_event::ChannelWaveform* high_gain,
    calin::ix::iact_data::telescope_event::ChannelWaveform* low_gain);
};

} } } // namespace calin::iact_data::event_visitor
