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
#include <provenance/chronicle.hpp>

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
    EventLifetimeManager* event_lifetime_manager,
    calin::ix::provenance::chronicle::ProcessingRecord* processing_record = nullptr);
  virtual bool leave_telescope_run(
    calin::ix::provenance::chronicle::ProcessingRecord* processing_record = nullptr);

  virtual bool visit_telescope_event(uint64_t seq_index,
    calin::ix::iact_data::telescope_event::TelescopeEvent* event);

  virtual bool merge_results();
};

class FilteredDelegatingParallelEventVisitor: public ParallelEventVisitor
{
public:
  FilteredDelegatingParallelEventVisitor();

  virtual ~FilteredDelegatingParallelEventVisitor();

  FilteredDelegatingParallelEventVisitor* new_sub_visitor(
    std::map<calin::iact_data::event_visitor::ParallelEventVisitor*,
        calin::iact_data::event_visitor::ParallelEventVisitor*>
      antecedent_visitors = { } /* note pass by value ! */ ) override;

  bool visit_telescope_run(
    const calin::ix::iact_data::telescope_run_configuration::
      TelescopeRunConfiguration* run_config,
      EventLifetimeManager* event_lifetime_manager,
      calin::ix::provenance::chronicle::ProcessingRecord* processing_record = nullptr) override;
  bool leave_telescope_run(
      calin::ix::provenance::chronicle::ProcessingRecord* processing_record = nullptr) override;

  bool visit_telescope_event(uint64_t seq_index,
    calin::ix::iact_data::telescope_event::TelescopeEvent* event) override;

  bool merge_results() override;

  void add_visitor(ParallelEventVisitor* visitor, const std::string& processing_record_comment = {}, bool adopt_visitor = false);

  void add_trigger_type_filtered_visitor(ParallelEventVisitor* visitor,
    calin::ix::iact_data::telescope_event::TriggerType trigger_type,
    const std::string& processing_record_comment = {},
    bool adopt_visitor = false);

  void add_physics_trigger_visitor(ParallelEventVisitor* visitor, const std::string& processing_record_comment = {}, bool adopt_visitor = false) {
    add_trigger_type_filtered_visitor(visitor, calin::ix::iact_data::telescope_event::TRIGGER_PHYSICS, processing_record_comment, adopt_visitor);
  }
  void add_software_trigger_visitor(ParallelEventVisitor* visitor, const std::string& processing_record_comment = {}, bool adopt_visitor = false) {
    add_trigger_type_filtered_visitor(visitor, calin::ix::iact_data::telescope_event::TRIGGER_SOFTWARE, processing_record_comment, adopt_visitor);
  }
  void add_pedestal_trigger_visitor(ParallelEventVisitor* visitor, const std::string& processing_record_comment = {}, bool adopt_visitor = false) {
    add_trigger_type_filtered_visitor(visitor, calin::ix::iact_data::telescope_event::TRIGGER_PEDESTAL, processing_record_comment, adopt_visitor);
  }
  void add_external_flasher_trigger_visitor(ParallelEventVisitor* visitor, const std::string& processing_record_comment = {}, bool adopt_visitor = false) {
    add_trigger_type_filtered_visitor(visitor, calin::ix::iact_data::telescope_event::TRIGGER_EXTERNAL_FLASHER, processing_record_comment, adopt_visitor);
  }
  void add_internal_flasher_trigger_visitor(ParallelEventVisitor* visitor, const std::string& processing_record_comment = {}, bool adopt_visitor = false) {
    add_trigger_type_filtered_visitor(visitor, calin::ix::iact_data::telescope_event::TRIGGER_INTERNAL_FLASHER, processing_record_comment, adopt_visitor);
  }

  bool visitor_saw_event(ParallelEventVisitor* visitor) const;

protected:
#ifndef SWIG
  struct DelegatedVisitor {
    DelegatedVisitor(ParallelEventVisitor* _visitor, bool _adopt_visitor,
        bool _unfiltered, calin::ix::iact_data::telescope_event::TriggerType _trigger_type,
        const std::string& _processing_record_comment):
      visitor(_visitor), adopt_visitor(_adopt_visitor), unfiltered(_unfiltered),
      trigger_type(_trigger_type), processing_record_comment(_processing_record_comment)
    { /* nothing to see here */ }

    ParallelEventVisitor* visitor;
    bool adopt_visitor;
    bool unfiltered;
    calin::ix::iact_data::telescope_event::TriggerType trigger_type;
    bool visitor_saw_event = false;
    std::string processing_record_comment;
    calin::ix::provenance::chronicle::ProcessingRecord* subprocessing_record = nullptr;
  };

  FilteredDelegatingParallelEventVisitor* parent_ = nullptr;
  std::map<calin::iact_data::event_visitor::ParallelEventVisitor*,
      calin::iact_data::event_visitor::ParallelEventVisitor*> parent_visitors_;

  std::vector<DelegatedVisitor> delegates_;
  std::map<calin::iact_data::event_visitor::ParallelEventVisitor*, int> visitor_delegates_;
#endif
};

} } } // namespace calin::iact_data::event_visitor
