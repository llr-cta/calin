/*

   calin/diagnostics/run_info.hpp -- Stephen Fegan -- 2018-10-26

   Run info diagnostics

   Copyright 2018, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, Ecole Polytechnique, CNRS/IN2P3

   This file is part of "calin"

   "calin" is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 or
   later, as published by the Free Software Foundation.

   "calin" is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

*/

#include <util/log.hpp>
#include <diagnostics/run_info.hpp>
#include <diagnostics/range.hpp>

using namespace calin::util::log;
using namespace calin::diagnostics::run_info;

RunInfoDiagnosticsVisitor::RunInfoDiagnosticsVisitor():
  calin::iact_data::event_visitor::ParallelEventVisitor(),
  arena_(new google::protobuf::Arena),
  results_(google::protobuf::Arena::CreateMessage<calin::ix::diagnostics::run_info::RunInfo>(arena_)),
  partials_(google::protobuf::Arena::CreateMessage<calin::ix::diagnostics::run_info::PartialRunInfo>(arena_)),
  run_config_(google::protobuf::Arena::CreateMessage<calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration>(arena_))
{
  // nothing to see here
}

RunInfoDiagnosticsVisitor::~RunInfoDiagnosticsVisitor()
{
  delete arena_;
}

RunInfoDiagnosticsVisitor* RunInfoDiagnosticsVisitor::new_sub_visitor(
  const std::map<calin::iact_data::event_visitor::ParallelEventVisitor*,
      calin::iact_data::event_visitor::ParallelEventVisitor*>&
    antecedent_visitors)
{
  RunInfoDiagnosticsVisitor* child = new RunInfoDiagnosticsVisitor();
  child->parent_ = this;
  return child;
}

bool RunInfoDiagnosticsVisitor::visit_telescope_run(
  const calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration* run_config,
  calin::iact_data::event_visitor::EventLifetimeManager* event_lifetime_manager)
{
  results_->Clear();
  partials_->Clear();
  run_config_->Clear();
  event_number_hist_.clear();
  elapsed_time_hist_.clear();
  run_config_->CopyFrom(*run_config);
  for(unsigned imod=0;imod<run_config->configured_module_id_size();imod++) {
    partials_->add_module();
  }
  return true;
}

bool RunInfoDiagnosticsVisitor::leave_telescope_run()
{
  return true;
}

bool RunInfoDiagnosticsVisitor::visit_telescope_event(uint64_t seq_index,
  calin::ix::iact_data::telescope_event::TelescopeEvent* event)
{
  partials_->increment_num_events_found();
  event_number_hist_.insert(event->local_event_number());
  elapsed_time_hist_.insert(event->elapsed_event_time().time_ns()*1e-9);

  partials_->add_event_number_sequence(event->local_event_number());

  // UCTS
  partials_->increment_num_events_missing_cdts_if(!event->has_cdts_data());
  calin::diagnostics::range::encode_value(
    partials_->mutable_cdts_presence(), event->has_cdts_data());

  // TIB
  partials_->increment_num_events_missing_tib_if(!event->has_tib_data());
  calin::diagnostics::range::encode_value(
    partials_->mutable_tib_presence(), event->has_tib_data());

  // SWAT
  partials_->increment_num_events_missing_swat_if(!event->has_swat_data());
  calin::diagnostics::range::encode_value(
    partials_->mutable_swat_presence(), event->has_swat_data());

  // ALL CHANNELS
  partials_->increment_num_events_missing_modules_if(!event->all_modules_present());
  calin::diagnostics::range::encode_value(
    partials_->mutable_all_channels_presence(), event->all_modules_present());

  // MODULES
  for(unsigned imodindex=0; imodindex<event->module_index_size(); imodindex++) {
    auto* mod = partials_->mutable_module(imodindex);
    if(event->module_index(imodindex)>=0) {
      mod->increment_num_events_present();
      calin::diagnostics::range::encode_value(mod->mutable_module_presence(), true);
    } else {
      calin::diagnostics::range::encode_value(mod->mutable_module_presence(), false);
    }
  }

  return true;
}

const calin::ix::diagnostics::run_info::RunInfo& RunInfoDiagnosticsVisitor::run_info()
{
  results_->Clear();
  for(unsigned imod=0;imod<run_config_->configured_module_id_size();imod++) {
    auto* mod = results_->add_module();
    mod->set_configured_module_rank(imod);
    mod->set_camera_module_id(run_config_->configured_module_id(imod));
  }
  integrate_histograms();
  integrate_partials();
  return *results_;
}

const calin::ix::diagnostics::run_info::PartialRunInfo& RunInfoDiagnosticsVisitor::partial_run_info()
{
  return *partials_;
}

bool RunInfoDiagnosticsVisitor::merge_results()
{
  integrate_histograms();
  parent_->partials_->IntegrateFrom(*partials_);
  partials_->Clear();
  return true;
}

void RunInfoDiagnosticsVisitor::integrate_histograms()
{
  auto* event_number_hist_data = event_number_hist_.dump_as_proto();
  partials_->mutable_event_number_histogram()->IntegrateFrom(*event_number_hist_data);
  delete event_number_hist_data;
  event_number_hist_.clear();

  auto* elapsed_time_hist_data = elapsed_time_hist_.dump_as_proto();
  partials_->mutable_elapsed_time_histogram()->IntegrateFrom(*elapsed_time_hist_data);
  delete elapsed_time_hist_data;
  elapsed_time_hist_.clear();
}

namespace {

  template<typename Iterator>
  void make_index_range_from_conditional_bool_rle(
    const Iterator& begin, const Iterator& end,
    const calin::ix::diagnostics::range::RunLengthEncodingBool& rle,
    calin::ix::diagnostics::range::IndexRange* range,
    bool value = false)
  {
    Iterator ievent = begin;
    std::vector<uint64_t> event_list;
    unsigned nrange = rle.count_size();
    for(unsigned irange=0;irange<nrange;irange++) {
      unsigned count = rle.count(irange);
      if(rle.value(irange) == value) {
        while(count--)event_list.push_back(*ievent++);
      } else {
        while(count--)ievent++;
      }
    }
    std::sort(event_list.begin(), event_list.end());
    calin::diagnostics::range::make_index_range(
      event_list.begin(), event_list.end(), range);
  }
}

void RunInfoDiagnosticsVisitor::integrate_partials()
{
  results_->increment_num_events_found(partials_->num_events_found());
  results_->increment_num_events_missing_cdts(partials_->num_events_missing_cdts());
  results_->increment_num_events_missing_tib(partials_->num_events_missing_tib());
  results_->increment_num_events_missing_swat(partials_->num_events_missing_swat());
  results_->increment_num_events_missing_modules(partials_->num_events_missing_modules());

  results_->mutable_event_number_histogram()->IntegrateFrom(partials_->event_number_histogram());
  results_->mutable_elapsed_time_histogram()->IntegrateFrom(partials_->elapsed_time_histogram());

  if(partials_->event_number_sequence_size() > 0)
  {
    std::vector<uint64_t> event_list;
    calin::ix::diagnostics::range::IndexRange range;
    event_list.reserve(partials_->event_number_sequence_size());

    for(auto ievent: partials_->event_number_sequence()) {
      event_list.push_back(ievent);
    }
    std::sort(event_list.begin(), event_list.end());
    unsigned last_event_number = event_list.front();
    for(unsigned iel=1;iel<event_list.size();iel++) {
      if(event_list[iel] == last_event_number) {
        results_->increment_num_duplicate_event_numbers();
        calin::diagnostics::range::encode_value(
          results_->mutable_duplicate_event_numbers(), last_event_number);
      }
      last_event_number = event_list[iel];
    }
    calin::diagnostics::range::make_index_range(
      event_list.begin(), event_list.end(), &range);
    results_->mutable_event_numbers_found()->IntegrateFrom(range);

    range.Clear();
    make_index_range_from_conditional_bool_rle(
      partials_->event_number_sequence().begin(), partials_->event_number_sequence().end(),
      partials_->cdts_presence(),
      &range, false);
    if(range.begin_index_size())
      results_->mutable_events_missing_cdts()->IntegrateFrom(range);

    range.Clear();
    make_index_range_from_conditional_bool_rle(
      partials_->event_number_sequence().begin(), partials_->event_number_sequence().end(),
      partials_->tib_presence(),
      &range, false);
    if(range.begin_index_size())
      results_->mutable_events_missing_tib()->IntegrateFrom(range);

    range.Clear();
    make_index_range_from_conditional_bool_rle(
      partials_->event_number_sequence().begin(), partials_->event_number_sequence().end(),
      partials_->swat_presence(),
      &range, false);
    if(range.begin_index_size())
      results_->mutable_events_missing_swat()->IntegrateFrom(range);

    range.Clear();
    make_index_range_from_conditional_bool_rle(
      partials_->event_number_sequence().begin(), partials_->event_number_sequence().end(),
      partials_->all_channels_presence(),
      &range, false);
    if(range.begin_index_size())
      results_->mutable_events_missing_modules()->IntegrateFrom(range);


    for(unsigned imod=0;imod<partials_->module_size();imod++) {
      range.Clear();
      make_index_range_from_conditional_bool_rle(
        partials_->event_number_sequence().begin(), partials_->event_number_sequence().end(),
        partials_->module(imod).module_presence(),
        &range, false);
      if(range.begin_index_size())
        results_->mutable_module(imod)->mutable_events_missing()->IntegrateFrom(range);
    }
  }
}
