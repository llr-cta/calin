/*

   calin/iact_data/event_visitor.cpp -- Stephen Fegan -- 2016-02-10

   A visitor of run and event data

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

#include <stdexcept>

#include <iact_data/event_visitor.hpp>
#include <util/log.hpp>

using namespace calin::iact_data::event_visitor;
using namespace calin::ix::iact_data::telescope_event;
using namespace calin::ix::iact_data::telescope_run_configuration;
using namespace calin::util::log;

EventLifetimeManager::~EventLifetimeManager()
{
  // nothing to see here
}

ParallelEventVisitor::~ParallelEventVisitor()
{
  // nothing to see here
}

ParallelEventVisitor* ParallelEventVisitor::new_sub_visitor(
  std::map<calin::iact_data::event_visitor::ParallelEventVisitor*,
    calin::iact_data::event_visitor::ParallelEventVisitor*> antecedent_visitors)
{
  return nullptr;
}

bool ParallelEventVisitor::visit_telescope_run(
  const calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration* run_config,
  EventLifetimeManager* event_lifetime_manager,
  calin::ix::provenance::chronicle::ProcessingRecord* processing_record)
{
  return true;
}

bool ParallelEventVisitor::leave_telescope_run(
  calin::ix::provenance::chronicle::ProcessingRecord* processing_record)
{
  return true;
}

bool ParallelEventVisitor::visit_telescope_event(uint64_t seq_index,
  calin::ix::iact_data::telescope_event::TelescopeEvent* event)
{
  return true;
}

bool ParallelEventVisitor::merge_results()
{
  return true;
}

FilteredDelegatingParallelEventVisitor::
FilteredDelegatingParallelEventVisitor():
  ParallelEventVisitor()
{
  // nothing to see here
}

FilteredDelegatingParallelEventVisitor::
~FilteredDelegatingParallelEventVisitor()
{
  for(auto ivisitor : delegates_) {
    if(ivisitor.adopt_visitor)delete ivisitor.visitor;
  }
}

FilteredDelegatingParallelEventVisitor*
FilteredDelegatingParallelEventVisitor::new_sub_visitor(
  std::map<calin::iact_data::event_visitor::ParallelEventVisitor*,
    calin::iact_data::event_visitor::ParallelEventVisitor*> antecedent_visitors)
{
  auto* sub_visitor = new FilteredDelegatingParallelEventVisitor();
  sub_visitor->parent_ = this;
  for(auto ivisitor : delegates_) {
    auto* dsv = ivisitor.visitor->new_sub_visitor(antecedent_visitors);
    if(dsv != nullptr) {
      if(ivisitor.unfiltered) {
        sub_visitor->add_visitor(dsv, ivisitor.processing_record_comment, /* adopt_visitor= */ true);
      } else {
        sub_visitor->add_trigger_type_filtered_visitor(dsv, ivisitor.trigger_type, ivisitor.processing_record_comment, /* adopt_visitor= */ true);
      }
      sub_visitor->parent_visitors_[dsv] = ivisitor.visitor;
    }
    antecedent_visitors[ivisitor.visitor] = dsv;
  }
  return sub_visitor;
}

bool FilteredDelegatingParallelEventVisitor::visit_telescope_run(
  const calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration* run_config,
  EventLifetimeManager* event_lifetime_manager,
  calin::ix::provenance::chronicle::ProcessingRecord* processing_record)
{
  bool good = true;
  for(auto ivisitor : delegates_) {
    ivisitor.sub_processing_record = nullptr;
    if(processing_record) {
      ivisitor.sub_processing_record = calin::provenance::chronicle::register_processing_start(
        __PRETTY_FUNCTION__, ivisitor.processing_record_comment);
      for(const auto& input : processing_record.primary_inputs()) {
        ivisitor.sub_processing_record->add_primary_inputs(input);
      }
    }
    good &= ivisitor.visitor->visit_telescope_run(run_config, event_lifetime_manager, ivisitor.sub_processing_record);
  }
  return good;
}

bool FilteredDelegatingParallelEventVisitor::leave_telescope_run(
  calin::ix::provenance::chronicle::ProcessingRecord* processing_record)
{
  bool good = true;
  for(auto ivisitor : delegates_) {
    if((processing_record == nullptr) ^ (ivisitor.sub_processing_record == nullptr)) {
      throw std::logic_error("FilteredDelegatingParallelEventVisitor::leave_telescope_run: inconsistent sub_processing_record");
    }
    good &= ivisitor.visitor->leave_telescope_run(ivisitor.sub_processing_record);
    if(ivisitor.sub_processing_record) {
      calin::provenance::chronicle::register_processing_finish(ivisitor.sub_processing_record);
      ivisitor.sub_processing_record = nullptr;
    }
  }
  return good;
}

bool FilteredDelegatingParallelEventVisitor::visit_telescope_event(uint64_t seq_index,
  calin::ix::iact_data::telescope_event::TelescopeEvent* event)
{
  bool good = true;
  for(auto& ivisitor : delegates_) {
    if(ivisitor.unfiltered or event->trigger_type() == ivisitor.trigger_type) {
      good &= ivisitor.visitor->visit_telescope_event(seq_index, event);
      ivisitor.visitor_saw_event = true;
    }
  }
  return good;
}

bool FilteredDelegatingParallelEventVisitor::merge_results()
{
  bool good = true;
  for(auto ivisitor : delegates_) {
    good &= ivisitor.visitor->merge_results();
    parent_->delegates_[parent_->visitor_delegates_.at(
      parent_visitors_.at(ivisitor.visitor))].visitor_saw_event |= ivisitor.visitor_saw_event;
  }
  return good;
}

void FilteredDelegatingParallelEventVisitor::add_visitor(
  ParallelEventVisitor* visitor, const std::string& processing_record_comment,
  bool adopt_visitor)
{
  delegates_.emplace_back(visitor, adopt_visitor, /*unfiltered=*/ true,
    calin::ix::iact_data::telescope_event::TRIGGER_UNKNOWN, processing_record_comment);
  visitor_delegates_[visitor] = delegates_.size() - 1;
}

void FilteredDelegatingParallelEventVisitor::
add_trigger_type_filtered_visitor(ParallelEventVisitor* visitor,
  calin::ix::iact_data::telescope_event::TriggerType trigger_type,
  const std::string& processing_record_comment,
  bool adopt_visitor)
{
  delegates_.emplace_back(visitor, adopt_visitor, /*unfiltered=*/ false, trigger_type, processing_record_comment);
  visitor_delegates_[visitor] = delegates_.size() - 1;
}

bool FilteredDelegatingParallelEventVisitor::visitor_saw_event(ParallelEventVisitor* visitor) const
{
  if(visitor_delegates_.count(visitor) == 0)return false;
  return delegates_[visitor_delegates_.at(visitor)].visitor_saw_event;
}
