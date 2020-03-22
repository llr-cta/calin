/*

   calin/iact_data/event_visitor.hpp -- Stephen Fegan -- 2016-01-30

   Vistor to event data

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
   Laboratoire Leprince-Ringuet, CNRS/IN2P3, Ecole Polytechnique, Institut Polytechnique de Paris

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
#include <iact_data/telescope_run_configuration.pb.h>
#include <iact_data/telescope_event.pb.h>

namespace calin { namespace iact_data { namespace event_visitor {

class EventLifetimeManager
{
public:
  virtual ~EventLifetimeManager();
  virtual void keep_event(
    calin::ix::iact_data::telescope_event::TelescopeEvent* event) = 0;
  virtual void release_event(
    calin::ix::iact_data::telescope_event::TelescopeEvent* event) = 0;
};

class TelescopeEventVisitor
{
public:
  virtual ~TelescopeEventVisitor();

  virtual bool demand_waveforms();
  virtual bool is_parallelizable();
  virtual TelescopeEventVisitor* new_sub_visitor(
    const std::map<calin::iact_data::event_visitor::TelescopeEventVisitor*,
        calin::iact_data::event_visitor::TelescopeEventVisitor*>&
      antecedent_visitors = { });

  virtual bool visit_telescope_run(
    const calin::ix::iact_data::telescope_run_configuration::
      TelescopeRunConfiguration* run_config);
  virtual bool leave_telescope_run();

  virtual bool visit_telescope_event(uint64_t seq_index,
    calin::ix::iact_data::telescope_event::TelescopeEvent* event);
  virtual bool leave_telescope_event();

  virtual bool visit_waveform(unsigned ichan,
    calin::ix::iact_data::telescope_event::ChannelWaveform* high_gain,
    calin::ix::iact_data::telescope_event::ChannelWaveform* low_gain);

  virtual bool merge_results();
};

class ParallelEventVisitor
{
public:
  virtual ~ParallelEventVisitor();

  virtual ParallelEventVisitor* new_sub_visitor(
    std::map<calin::iact_data::event_visitor::ParallelEventVisitor*,
        calin::iact_data::event_visitor::ParallelEventVisitor*>
      antecedent_visitors = { } /* note pass by value ! */ );

  virtual bool visit_telescope_run(
    const calin::ix::iact_data::telescope_run_configuration::
      TelescopeRunConfiguration* run_config,
      EventLifetimeManager* event_lifetime_manager);
  virtual bool leave_telescope_run();

  virtual bool visit_telescope_event(uint64_t seq_index,
    calin::ix::iact_data::telescope_event::TelescopeEvent* event);

  virtual bool merge_results();
};

} } } // namespace calin::iact_data::event_visitor
