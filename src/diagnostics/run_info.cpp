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

ModuleCounterProcessor::~ModuleCounterProcessor()
{
  // nothing to see here
}

DirectModuleCounterProcessor::~DirectModuleCounterProcessor()
{
  // nothing to see here
}

int64_t DirectModuleCounterProcessor::
processCounters(std::vector<int64_t>& values, int64_t event_number)
{
  return 0;
}

RelativeToEventNumberModuleCounterProcessor::
~RelativeToEventNumberModuleCounterProcessor()
{
  // nothing to see here
}

int64_t RelativeToEventNumberModuleCounterProcessor::
processCounters(std::vector<int64_t>& values, int64_t event_number)
{
  for(auto& iv: values) { iv -= event_number; }
  return event_number;
}

RelativeToMedianModuleCounterProcessor::
~RelativeToMedianModuleCounterProcessor()
{
  // nothing to see here
}

int64_t RelativeToMedianModuleCounterProcessor::
processCounters(std::vector<int64_t>& values, int64_t event_number)
{
  median_find_.clear();
  for(const auto& iv: values) { median_find_[iv]++; }
  unsigned count = 0;
  unsigned median_count = values.size()/2;
  auto imedian_find = median_find_.begin();
  while((count+=imedian_find->second) < median_count)++imedian_find;
  int64_t offset = imedian_find->first;
  for(auto& iv: values) { iv -= offset; }
  return offset;
}

RunInfoDiagnosticsVisitor::RunInfoDiagnosticsVisitor(
    const calin::ix::diagnostics::run_info::RunInfoConfig& config):
  calin::iact_data::event_visitor::ParallelEventVisitor(),
  config_(config),
  arena_(new google::protobuf::Arena),
  results_(google::protobuf::Arena::CreateMessage<calin::ix::diagnostics::run_info::RunInfo>(arena_)),
  partials_(google::protobuf::Arena::CreateMessage<calin::ix::diagnostics::run_info::PartialRunInfo>(arena_)),
  run_config_(google::protobuf::Arena::CreateMessage<calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration>(arena_))
{
  // nothing to see here
}

RunInfoDiagnosticsVisitor::~RunInfoDiagnosticsVisitor()
{
  for(auto* iproc: mod_counter_processor_)delete iproc;
  delete arena_;
}

RunInfoDiagnosticsVisitor* RunInfoDiagnosticsVisitor::new_sub_visitor(
  const std::map<calin::iact_data::event_visitor::ParallelEventVisitor*,
      calin::iact_data::event_visitor::ParallelEventVisitor*>&
    antecedent_visitors)
{
  RunInfoDiagnosticsVisitor* child = new RunInfoDiagnosticsVisitor(config_);
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

  calin::ix::diagnostics::run_info::RunInfoConfig config = config_;
  if(config.module_counter_test_id_size() == 0 and
      config.module_counter_test_mode_size() == 0) {
    if(run_config_->camera_layout().camera_type() ==
        calin::ix::iact_data::instrument_layout::CameraLayout::NECTARCAM) {
      *(config.mutable_module_counter_test_id()) = config.default_nectarcam_config().module_counter_test_id();
      *(config.mutable_module_counter_test_mode()) = config.default_nectarcam_config().module_counter_test_mode();
    } else if(run_config_->camera_layout().camera_type() ==
        calin::ix::iact_data::instrument_layout::CameraLayout::LSTCAM) {
      *(config.mutable_module_counter_test_id()) = config.default_lstcam_config().module_counter_test_id();
      *(config.mutable_module_counter_test_mode()) = config.default_lstcam_config().module_counter_test_mode();
    }
  }

  for(auto* iproc: mod_counter_processor_)delete iproc;
  mod_counter_id_.clear();
  mod_counter_mode_.clear();
  mod_counter_processor_.clear();
  for(int icounter=0;icounter<std::max(config.module_counter_test_id_size(),
      config.module_counter_test_mode_size());icounter++) {
    int counter_id = config.module_counter_test_id(icounter);
    if(counter_id<0 or counter_id>=run_config_->camera_layout().module_counter_name_size())continue;
    mod_counter_id_.push_back(counter_id);
    mod_counter_mode_.push_back(config.module_counter_test_mode(icounter));
    if(icounter<config.module_counter_test_mode_size()
        and config.module_counter_test_mode(icounter) == calin::ix::diagnostics::run_info::VALUE_RELATIVE_TO_MEDIAN) {
      mod_counter_processor_.push_back(new RelativeToMedianModuleCounterProcessor());
    } else if(icounter<config.module_counter_test_mode_size()
        and config.module_counter_test_mode(icounter) == calin::ix::diagnostics::run_info::VALUE_RELATIVE_TO_EVENT_NUMBER) {
      mod_counter_processor_.push_back(new RelativeToEventNumberModuleCounterProcessor());
    } else {
      mod_counter_processor_.push_back(new DirectModuleCounterProcessor());
    }
  }
  if(not mod_counter_id_.empty())
    mod_counter_values_.reserve(run_config->configured_module_id_size());

  for(unsigned imod=0;imod<run_config->configured_module_id_size();imod++) {
    auto* m = partials_->add_module();
    for(unsigned icounter=0;icounter<mod_counter_id_.size();icounter++) {
      m->add_counter_value();
    }
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
  for(int imodindex=0; imodindex<event->module_index_size(); imodindex++) {
    auto* mod = partials_->mutable_module(imodindex);
    if(event->module_index(imodindex)>=0) {
      mod->increment_num_events_present();
      calin::diagnostics::range::encode_value(mod->mutable_module_presence(), true);
    } else {
      calin::diagnostics::range::encode_value(mod->mutable_module_presence(), false);
    }
  }

  // MODULE COUNTERS
  for(unsigned icounter=0;icounter<mod_counter_id_.size(); icounter++) {
    int counter_id = mod_counter_id_[icounter];
    mod_counter_values_.clear();
    for(int imod=0; imod<event->module_counter_size(); imod++) {
      int64_t val = event->module_counter(imod).counter_value(counter_id);
      mod_counter_values_.push_back(val);
    }
    /* int64_t offset = */ mod_counter_processor_[icounter]->processCounters(
      mod_counter_values_, event->local_event_number());
    for(int imod=0; imod<event->module_counter_size(); imod++) {
      unsigned mod_id = event->module_counter(imod).module_id();
      auto* mod = partials_->mutable_module(mod_id);
      calin::diagnostics::range::encode_value(
        mod->mutable_counter_value(icounter), mod_counter_values_[imod]);
    }
  }

  return true;
}

calin::ix::diagnostics::run_info::RunInfoConfig RunInfoDiagnosticsVisitor::default_config()
{
  calin::ix::diagnostics::run_info::RunInfoConfig config;
  config.mutable_default_nectarcam_config()->add_module_counter_test_id(0);
  config.mutable_default_nectarcam_config()->add_module_counter_test_id(1);
  config.mutable_default_nectarcam_config()->add_module_counter_test_id(2);
  config.mutable_default_nectarcam_config()->add_module_counter_test_mode(
    calin::ix::diagnostics::run_info::VALUE_RELATIVE_TO_EVENT_NUMBER);
  config.mutable_default_nectarcam_config()->add_module_counter_test_mode(
    calin::ix::diagnostics::run_info::VALUE_RELATIVE_TO_MEDIAN);
  config.mutable_default_nectarcam_config()->add_module_counter_test_mode(
    calin::ix::diagnostics::run_info::VALUE_RELATIVE_TO_MEDIAN);

  config.mutable_default_lstcam_config()->add_module_counter_test_id(0);
  config.mutable_default_lstcam_config()->add_module_counter_test_id(1);
  config.mutable_default_lstcam_config()->add_module_counter_test_id(2);
  config.mutable_default_lstcam_config()->add_module_counter_test_mode(
    calin::ix::diagnostics::run_info::VALUE_RELATIVE_TO_MEDIAN);
  config.mutable_default_lstcam_config()->add_module_counter_test_mode(
    calin::ix::diagnostics::run_info::VALUE_RELATIVE_TO_MEDIAN);
  config.mutable_default_lstcam_config()->add_module_counter_test_mode(
    calin::ix::diagnostics::run_info::VALUE_RELATIVE_TO_MEDIAN);

  return config;
}

const calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration&
RunInfoDiagnosticsVisitor::run_config() const
{
  return *run_config_;
}

const calin::ix::diagnostics::run_info::RunInfo& RunInfoDiagnosticsVisitor::run_info()
{
  results_->Clear();
  for(int imod=0;imod<run_config_->configured_module_id_size();imod++) {
    auto* mod = results_->add_module();
    mod->set_configured_module_rank(imod);
    mod->set_camera_module_id(run_config_->configured_module_id(imod));
    for(unsigned icounter=0;icounter<mod_counter_id_.size(); icounter++) {
      auto* counter = mod->add_counter_value();
      counter->set_counter_id(mod_counter_id_[icounter]);
      counter->set_counter_name(run_config_->camera_layout().module_counter_name(
        mod_counter_id_[icounter]));
      counter->set_test_mode(
        calin::ix::diagnostics::run_info::CounterValueTestMode(mod_counter_mode_[icounter]));
    }
  }
  integrate_histograms();
  integrate_partials();
  return *results_;
}

const calin::ix::diagnostics::run_info::PartialRunInfo& RunInfoDiagnosticsVisitor::partial_run_info() const
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

  template<typename Iterable, typename PresenceRLE, typename ValueRLE, typename IndexValueRange>
  void make_index_range_from_value_rle(
    const Iterable& indexes, const PresenceRLE& has_rle, const ValueRLE& val_rle,
    IndexValueRange* ivr)
  {
    std::map<uint64_t, unsigned> index_and_value;
    unsigned val_irle = 0;
    unsigned val_nrle = val_rle.count(val_irle);
    unsigned has_irle = 0;
    unsigned has_nrle = has_rle.count(has_irle);
    for(auto index : indexes) {
      if(has_nrle==0) { ++has_irle; has_nrle = has_rle.count(has_irle); }
      if(val_nrle==0) { ++val_irle; val_nrle = val_rle.count(val_irle); }
      if(has_rle.value(has_irle)) {
        index_and_value[index] = val_irle;
        --val_nrle;
      }
      --has_nrle;
    }
    for(auto iv : index_and_value) {
      calin::diagnostics::range::encode_monotonic_index_and_value(
        ivr, iv.first, val_rle.value(iv.second));
    }
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

    for(int imod=0;imod<partials_->module_size();imod++) {
      auto* rimod = results_->mutable_module(imod);
      auto& pimod = partials_->module(imod);

      range.Clear();
      make_index_range_from_conditional_bool_rle(
        partials_->event_number_sequence().begin(), partials_->event_number_sequence().end(),
        pimod.module_presence(), &range, false);
      if(range.begin_index_size())
        results_->mutable_module(imod)->mutable_events_missing()->IntegrateFrom(range);
      rimod->set_num_events_present(pimod.num_events_present());

      for(unsigned icounter=0; icounter<mod_counter_id_.size(); icounter++) {
        make_index_range_from_value_rle(partials_->event_number_sequence(),
          pimod.module_presence(), pimod.counter_value(icounter),
          rimod->mutable_counter_value(icounter)->mutable_value_range());
      }
    }
  }
}
