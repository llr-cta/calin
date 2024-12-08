/*

   calin/diagnostics/run_info.hpp -- Stephen Fegan -- 2018-10-26

   Run info diagnostics

   Copyright 2018, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <util/log.hpp>
#include <util/algorithm.hpp>
#include <diagnostics/run_info.hpp>
#include <diagnostics/range.hpp>
#include <io/json.hpp>
#include <math/special.hpp>

using namespace calin::util::log;
using namespace calin::diagnostics::run_info;
using calin::math::special::SQR;

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

RunInfoDiagnosticsParallelEventVisitor::RunInfoDiagnosticsParallelEventVisitor(
    const calin::ix::diagnostics::run_info::RunInfoConfig& config):
  calin::iact_data::event_visitor::ParallelEventVisitor(),
  config_(config),
  arena_(new google::protobuf::Arena),
  results_(google::protobuf::Arena::CreateMessage<calin::ix::diagnostics::run_info::RunInfo>(arena_)),
  partials_(google::protobuf::Arena::CreateMessage<calin::ix::diagnostics::run_info::PartialRunInfo>(arena_)),
  run_config_(google::protobuf::Arena::CreateMessage<calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration>(arena_))
{
  muon_candidate_elapsed_time_hist_ = calin::math::histogram::SimpleHist(
    config_.event_time_histogram_resolution(),
    config_.event_time_histogram_min(), config_.event_time_histogram_max());
}

RunInfoDiagnosticsParallelEventVisitor::~RunInfoDiagnosticsParallelEventVisitor()
{
  for(auto* iproc: mod_counter_processor_)delete iproc;
  for(auto& mod_hists : mod_counter_hist_) {
    for(auto* hist : mod_hists)delete hist;
  }
  delete arena_;
}

RunInfoDiagnosticsParallelEventVisitor* RunInfoDiagnosticsParallelEventVisitor::new_sub_visitor(
  std::map<calin::iact_data::event_visitor::ParallelEventVisitor*,
      calin::iact_data::event_visitor::ParallelEventVisitor*>
    antecedent_visitors)
{
  RunInfoDiagnosticsParallelEventVisitor* child = new RunInfoDiagnosticsParallelEventVisitor(config_);
  child->parent_ = this;
  return child;
}

bool RunInfoDiagnosticsParallelEventVisitor::visit_telescope_run(
  const calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration* run_config,
  calin::iact_data::event_visitor::EventLifetimeManager* event_lifetime_manager,
  calin::ix::provenance::chronicle::ProcessingRecord* processing_record)
{
  if(processing_record) {
    processing_record->set_type("RunInfoDiagnosticsParallelEventVisitor");
    processing_record->set_description("Miscellaneous run information");
    auto* config_json = processing_record->add_config();
    config_json->set_type(config_.GetTypeName());
    config_json->set_json(calin::io::json::encode_protobuf_to_json_string(config_));
  }

  results_->Clear();
  partials_->Clear();
  run_config_->Clear();

  run_config_->CopyFrom(*run_config);

  results_->set_min_event_time(std::numeric_limits<int64_t>::max());
  results_->set_max_event_time(std::numeric_limits<int64_t>::min());

  partials_->set_min_event_time(std::numeric_limits<int64_t>::max());
  partials_->set_max_event_time(std::numeric_limits<int64_t>::min());

  for(int iclock=0; iclock<run_config_->camera_layout().camera_clock_name_size(); ++iclock) {
    results_->add_camera_clock_presence(0);
    results_->add_camera_clock_min_time(std::numeric_limits<int64_t>::max());
    results_->add_camera_clock_max_time(std::numeric_limits<int64_t>::min());
    partials_->add_camera_clock_presence(0);
    partials_->add_camera_clock_min_time(std::numeric_limits<int64_t>::max());
    partials_->add_camera_clock_max_time(std::numeric_limits<int64_t>::min());
  }

  trigger_type_code_hist_.clear();
  trigger_type_code_diff_hist_.clear();
  muon_candidate_elapsed_time_hist_.clear();

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
  for(auto& mod_hists : mod_counter_hist_) {
    for(auto* hist : mod_hists)delete hist;
  }
  mod_counter_hist_.clear();

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
    mod_counter_hist_.emplace_back(run_config->configured_module_id_size());
  }
  if(not mod_counter_id_.empty())
    mod_counter_values_.reserve(run_config->configured_module_id_size());

  for(int imod=0;imod<run_config->configured_module_id_size();imod++) {
    auto* m = partials_->add_module();
    for(unsigned icounter=0;icounter<mod_counter_id_.size();icounter++) {
      if(config_.enable_module_counter_test_value_range()) {
        m->add_counter_value();
      }
      m->add_counter_value_squared_sum_histogram();
    }
  }

  return true;
}

bool RunInfoDiagnosticsParallelEventVisitor::leave_telescope_run(
  calin::ix::provenance::chronicle::ProcessingRecord* processing_record)
{
  for(int imod=0;imod<run_config_->configured_module_id_size();imod++) {
    for(unsigned icounter=0;icounter<mod_counter_id_.size(); icounter++) {
      if(mod_counter_hist_[icounter][imod] != nullptr) {
        auto* pimod = partials_->mutable_module(imod);
        mod_counter_hist_[icounter][imod]->dump_as_proto(
          pimod->mutable_counter_value_squared_sum_histogram(icounter));
        delete mod_counter_hist_[icounter][imod];
      }
    }
  }
  mod_counter_hist_.clear();

  if(parent_ == nullptr) {
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
    integrate_partials();
  }
  return true;
}

bool RunInfoDiagnosticsParallelEventVisitor::visit_telescope_event(uint64_t seq_index,
  calin::ix::iact_data::telescope_event::TelescopeEvent* event)
{
  using namespace calin::ix::iact_data::telescope_event;

  partials_->increment_num_events_found();

  partials_->add_event_number_sequence(event->local_event_number());
  partials_->add_event_time_sequence(event->absolute_event_time().time_ns());
  partials_->add_event_type_sequence(event->trigger_type());

  // UCTS
  partials_->increment_num_events_missing_cdts_if(!event->has_cdts_data());
  calin::diagnostics::range::encode_value(
    partials_->mutable_cdts_presence(), event->has_cdts_data());

  // TIB
  partials_->increment_num_events_missing_tib_if(!event->has_tib_data());
  calin::diagnostics::range::encode_value(
    partials_->mutable_tib_presence(), event->has_tib_data());

  // UCTS & TIB
  partials_->increment_num_events_missing_tib_and_cdts_if(
    !event->has_tib_data() and !event->has_cdts_data());

  // SWAT
  partials_->increment_num_events_missing_swat_if(!event->has_swat_data());
  calin::diagnostics::range::encode_value(
    partials_->mutable_swat_presence(), event->has_swat_data());

  // ALL CHANNELS
  partials_->increment_num_events_missing_modules_if(!event->all_modules_present());
  calin::diagnostics::range::encode_value(
    partials_->mutable_all_channels_presence(), event->all_modules_present());

  // TIB TRIGGER BITS (from TIB or fallback to CDTS)
  if(event->has_tib_data()) {
    const auto& tib = event->tib_data();
    // partials_->increment_num_mono_trigger(tib.mono_trigger()); // NOT IN ICD
    // partials_->increment_num_stereo_trigger(tib.stereo_trigger()); // NOT IN ICD
    partials_->increment_num_external_calibration_trigger(tib.external_calibration_trigger());
    partials_->increment_num_internal_calibration_trigger(tib.internal_calibration_trigger());
    partials_->increment_num_ucts_aux_trigger(tib.ucts_aux_trigger());
    partials_->increment_num_pedestal_trigger(tib.pedestal_trigger());
    partials_->increment_num_slow_control_trigger(tib.slow_control_trigger());
    // partials_->increment_num_busy_trigger(tib.busy_trigger()); // NOT IN ICD
    trigger_type_code_hist_.insert(tib.trigger_type());
    if(event->has_cdts_data()) {
      const auto& cdts = event->cdts_data();
      partials_->increment_num_mono_trigger(cdts.mono_trigger());
      partials_->increment_num_stereo_trigger(cdts.stereo_trigger());
      partials_->increment_num_busy_trigger(cdts.busy_trigger());
      trigger_type_code_diff_hist_.insert(int(tib.trigger_type()) - int(cdts.trigger_type() & 0x7c));
      partials_->increment_num_tib_ucts_trigger_code_mismatch_if(tib.trigger_type() != (cdts.trigger_type() & 0x7c));
    } else {
      switch(event->trigger_type()) {
        case TRIGGER_PHYSICS: partials_->increment_num_mono_trigger(); break;
        case TRIGGER_FORCED_BY_ARRAY: partials_->increment_num_stereo_trigger(); break;
        default: break;
      }
    }
  } else if(event->has_cdts_data()) {
    const auto& cdts = event->cdts_data();
    partials_->increment_num_mono_trigger(cdts.mono_trigger());
    partials_->increment_num_stereo_trigger(cdts.stereo_trigger());
    partials_->increment_num_external_calibration_trigger(cdts.external_calibration_trigger());
    partials_->increment_num_internal_calibration_trigger(cdts.internal_calibration_trigger());
    partials_->increment_num_ucts_aux_trigger(cdts.ucts_aux_trigger());
    partials_->increment_num_pedestal_trigger(cdts.pedestal_trigger());
    partials_->increment_num_slow_control_trigger(cdts.slow_control_trigger());
    partials_->increment_num_busy_trigger(cdts.busy_trigger());
    trigger_type_code_hist_.insert(cdts.trigger_type());
  } else {
    switch(event->trigger_type()) {
      case TRIGGER_PHYSICS: partials_->increment_num_mono_trigger(); break;
      case TRIGGER_SOFTWARE: partials_->increment_num_slow_control_trigger(); break;
      case TRIGGER_PEDESTAL: partials_->increment_num_pedestal_trigger(); break;
      case TRIGGER_EXTERNAL_FLASHER: partials_->increment_num_external_calibration_trigger(); break;
      case TRIGGER_INTERNAL_FLASHER: partials_->increment_num_internal_calibration_trigger(); break;
      case TRIGGER_FORCED_BY_ARRAY: partials_->increment_num_stereo_trigger(); break;
      case TRIGGER_UCTS_AUX: partials_->increment_num_ucts_aux_trigger(); break;
      case TRIGGER_UNKNOWN:
      case TRIGGER_MULTIPLE:
      default: break;
    };
  }

  // MUON CANDIDATE
  partials_->increment_num_muon_candidate(event->is_muon_candidate());

  // EVENT TIME
  if(event->has_absolute_event_time() and event->absolute_event_time().time_ns()>0) {
    auto t = event->absolute_event_time().time_ns();
    partials_->set_max_event_time(std::max(partials_->max_event_time(), t));
    partials_->set_min_event_time(std::min(partials_->min_event_time(), t));
    if(event->is_muon_candidate()) {
      double elapsed_time_sec = (t-run_config_->run_start_time().time_ns())*1e-9;
      muon_candidate_elapsed_time_hist_.insert(elapsed_time_sec);
    }
  }

  // CAMERA CLOCK TIMES
  for(const auto& clock : event->camera_clock()) {
    if(clock.time_value()>0) {
      partials_->increment_camera_clock_presence(clock.clock_id());
      partials_->set_camera_clock_min_time(clock.clock_id(),
        std::min(partials_->camera_clock_min_time(clock.clock_id()), clock.time_value()));
      partials_->set_camera_clock_max_time(clock.clock_id(),
        std::max(partials_->camera_clock_max_time(clock.clock_id()), clock.time_value()));
    }
  }

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

      if(mod_counter_values_[imod]) {
        if(mod_counter_hist_[icounter][mod_id] == nullptr) {
          mod_counter_hist_[icounter][mod_id] =
            new calin::math::histogram::SimpleHist(config_.event_number_histogram_resolution());
        }
        mod_counter_hist_[icounter][mod_id]->insert(event->local_event_number(),
          SQR(mod_counter_values_[imod]));
      }

      if(config_.enable_module_counter_test_value_range()) {
        calin::diagnostics::range::encode_value(
          mod->mutable_counter_value(icounter), mod_counter_values_[imod]);
      }
    }
  }

  return true;
}

calin::ix::diagnostics::run_info::RunInfoConfig RunInfoDiagnosticsParallelEventVisitor::default_config()
{
  calin::ix::diagnostics::run_info::RunInfoConfig config;

  config.set_event_number_histogram_resolution(10000);
  config.set_event_time_histogram_resolution(1.0);
  config.set_event_time_histogram_max(7200.0);
  config.set_event_time_histogram_min(-60.0);  
  config.set_log10_delta_t_histogram_binsize(0.01);
  config.set_delta_t_timeslice(100.0);

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

calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration*
RunInfoDiagnosticsParallelEventVisitor::run_config(
  calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration* rc) const
{
  if(rc == nullptr)rc = run_config_->New();
  rc->CopyFrom(*run_config_);
  return rc;
}

calin::ix::diagnostics::run_info::RunInfo*
RunInfoDiagnosticsParallelEventVisitor::run_info(
  calin::ix::diagnostics::run_info::RunInfo* ri) const
{
  if(ri == nullptr)ri = results_->New();
  ri->CopyFrom(*results_);
  return ri;
}

const calin::ix::diagnostics::run_info::PartialRunInfo& RunInfoDiagnosticsParallelEventVisitor::partial_run_info() const
{
  return *partials_;
}

bool RunInfoDiagnosticsParallelEventVisitor::merge_results()
{
  if(parent_) {
    parent_->partials_->IntegrateFrom(*partials_);
    partials_->Clear();
    parent_->trigger_type_code_hist_.insert_hist(trigger_type_code_hist_);
    trigger_type_code_hist_.clear();
    parent_->trigger_type_code_diff_hist_.insert_hist(trigger_type_code_diff_hist_);
    trigger_type_code_diff_hist_.clear();
    parent_->muon_candidate_elapsed_time_hist_.insert_hist(muon_candidate_elapsed_time_hist_);
    muon_candidate_elapsed_time_hist_.clear();
  }
  return true;
}

namespace {

  template<typename Iterator>
  void make_index_range_from_conditional_bool_rle(
    const std::vector<size_t>& event_index,
    const Iterator& event_number_begin,
    const calin::ix::diagnostics::range::RunLengthEncodingBool& rle,
    std::vector<bool>& event_value,
    calin::ix::diagnostics::range::IndexRange* range,
    bool value = false)
  {
    event_value.resize(event_index.size());

    for(unsigned irange=0, iindex=0;irange<unsigned(rle.count_size());irange++) {
      uint64_t count = rle.count(irange);
      bool value =  rle.value(irange);
      while(count--) {
        event_value[iindex++] = value;
      }
    }
    for(unsigned iindex=0;iindex<event_index.size();iindex++) {
      if(event_value[event_index[iindex]] == value) {
        calin::diagnostics::range::encode_monotonic_index(range,
          *(event_number_begin + event_index[iindex]));
      }
    }
  }

  template<typename Iterator>
  void make_index_range_from_i64_value_rle(
    const std::vector<size_t>& event_index,
    const Iterator& event_number_begin,
    const calin::ix::diagnostics::range::RunLengthEncodingBool& has_rle,
    const calin::ix::diagnostics::range::RunLengthEncodingInt64& val_rle,
    std::vector<bool>& event_has,
    std::vector<int64_t>& event_val,
    calin::ix::diagnostics::range::IndexAndValueRangeInt64* range)
  {
    event_has.resize(event_index.size());
    event_val.resize(event_index.size());

    for(unsigned irange=0, iindex=0;irange<unsigned(has_rle.count_size());irange++) {
      uint64_t count = has_rle.count(irange);
      bool has = has_rle.value(irange);
      while(count--) {
        event_has[iindex++] = has;
      }
    }

    for(unsigned irange=0, iindex=0;irange<unsigned(val_rle.count_size());irange++) {
      uint64_t count = val_rle.count(irange);
      int64_t val = val_rle.value(irange);
      while(count) {
        event_val[iindex] = val;
        if(event_has[iindex])count--;
        iindex++;
      }
    }

    for(unsigned iindex=0;iindex<event_index.size();iindex++) {
      if(event_has[event_index[iindex]]) {
        calin::diagnostics::range::encode_monotonic_index_and_value(range,
          *(event_number_begin + event_index[iindex]), event_val[event_index[iindex]], true);
      }
    }
  }
}

void RunInfoDiagnosticsParallelEventVisitor::integrate_partials()
{
  results_->set_num_events_found(partials_->num_events_found());
  results_->set_num_events_missing_cdts(partials_->num_events_missing_cdts());
  results_->set_num_events_missing_tib(partials_->num_events_missing_tib());
  results_->set_num_events_missing_swat(partials_->num_events_missing_swat());
  results_->set_num_events_missing_modules(partials_->num_events_missing_modules());
  results_->set_num_events_missing_tib_and_cdts(partials_->num_events_missing_tib_and_cdts());

  results_->set_num_mono_trigger(partials_->num_mono_trigger());
  results_->set_num_stereo_trigger(partials_->num_stereo_trigger()) ;
  results_->set_num_external_calibration_trigger(partials_->num_external_calibration_trigger());
  results_->set_num_internal_calibration_trigger(partials_->num_internal_calibration_trigger());
  results_->set_num_ucts_aux_trigger(partials_->num_ucts_aux_trigger());
  results_->set_num_pedestal_trigger(partials_->num_pedestal_trigger());
  results_->set_num_slow_control_trigger(partials_->num_slow_control_trigger());
  results_->set_num_busy_trigger(partials_->num_busy_trigger());
  results_->set_num_muon_candidate(partials_->num_muon_candidate());
  results_->set_num_tib_ucts_trigger_code_mismatch(partials_->num_tib_ucts_trigger_code_mismatch());

  results_->set_min_event_time(partials_->min_event_time());
  results_->set_max_event_time(partials_->max_event_time());

  (*results_->mutable_camera_clock_presence()) = partials_->camera_clock_presence();
  (*results_->mutable_camera_clock_max_time()) = partials_->camera_clock_max_time();
  (*results_->mutable_camera_clock_min_time()) = partials_->camera_clock_min_time();

  trigger_type_code_hist_.dump_as_compactified_proto(0, 0,
    results_->mutable_trigger_code_histogram());
  trigger_type_code_diff_hist_.dump_as_compactified_proto(0, 0,
    results_->mutable_tib_ucts_trigger_code_diff_histogram());

  muon_candidate_elapsed_time_hist_.dump_as_proto(
    results_->mutable_elapsed_time_histogram_muon_candidate());

  if(partials_->event_number_sequence_size() > 0)
  {
    calin::math::histogram::Histogram1D event_number_hist {
      config_.event_number_histogram_resolution(), 0.0, 1.0e9,
      double(run_config_->camera_layout().first_event_number()) };

    double thistres = config_.event_time_histogram_resolution();
    double thistmin = config_.event_time_histogram_min();
    double thistmax = config_.event_time_histogram_max();

    // Check that no more than 5 percent of events are outside the window, else hardcode
    // histograms for 1 minute resolution over max of 48 hours. This does not work for
    // the muon histogram unfortunately
    unsigned num_event_outside_window = 0;
    for(int ievent=0; ievent<partials_->event_number_sequence_size(); ++ievent) {
        double elapsed_time_sec =
          (partials_->event_time_sequence(ievent)-run_config_->run_start_time().time_ns())*1e-9;
        if(elapsed_time_sec > thistmax)++num_event_outside_window;
        if(num_event_outside_window*20 > partials_->event_number_sequence_size()) {
          thistres = 60.0;
          thistmin = -360.0;
          thistmax = 172800.0;
          break;
        }
    }

    calin::math::histogram::Histogram1D elapsed_time_hist {
      thistres, thistmin, thistmax, 0.0 };
    calin::math::histogram::Histogram1D trigger_physics_elapsed_time_hist {
      thistres, thistmin, thistmax, 0.0 };
    calin::math::histogram::Histogram1D trigger_software_elapsed_time_hist {
      thistres, thistmin, thistmax, 0.0 };
    calin::math::histogram::Histogram1D trigger_pedestal_elapsed_time_hist {
      thistres, thistmin, thistmax, 0.0 };
    calin::math::histogram::Histogram1D trigger_external_flasher_elapsed_time_hist {
      thistres, thistmin, thistmax, 0.0 };
    calin::math::histogram::Histogram1D trigger_internal_flasher_elapsed_time_hist {
      thistres, thistmin, thistmax, 0.0 };
    calin::math::histogram::Histogram1D trigger_forced_array_elapsed_time_hist {
      thistres, thistmin, thistmax, 0.0 };
    calin::math::histogram::Histogram1D trigger_ucts_aux_elapsed_time_hist {
      thistres, thistmin, thistmax, 0.0 };

    for(int ievent=0; ievent<partials_->event_number_sequence_size(); ++ievent) {
      event_number_hist.insert(partials_->event_number_sequence(ievent));
      if(partials_->event_time_sequence(ievent) > 0) {
        double elapsed_time_sec =
          (partials_->event_time_sequence(ievent)-run_config_->run_start_time().time_ns())*1e-9;
        elapsed_time_hist.insert(elapsed_time_sec);
        switch(partials_->event_type_sequence(ievent)) {
        case calin::ix::iact_data::telescope_event::TRIGGER_PHYSICS:
          trigger_physics_elapsed_time_hist.insert(elapsed_time_sec); break;
        case calin::ix::iact_data::telescope_event::TRIGGER_SOFTWARE:
          trigger_software_elapsed_time_hist.insert(elapsed_time_sec); break;
        case calin::ix::iact_data::telescope_event::TRIGGER_PEDESTAL:
          trigger_pedestal_elapsed_time_hist.insert(elapsed_time_sec); break;
        case calin::ix::iact_data::telescope_event::TRIGGER_EXTERNAL_FLASHER:
          trigger_external_flasher_elapsed_time_hist.insert(elapsed_time_sec); break;
        case calin::ix::iact_data::telescope_event::TRIGGER_INTERNAL_FLASHER:
          trigger_internal_flasher_elapsed_time_hist.insert(elapsed_time_sec); break;
        case calin::ix::iact_data::telescope_event::TRIGGER_FORCED_BY_ARRAY:
          trigger_forced_array_elapsed_time_hist.insert(elapsed_time_sec); break;
        case calin::ix::iact_data::telescope_event::TRIGGER_UCTS_AUX:
          trigger_ucts_aux_elapsed_time_hist.insert(elapsed_time_sec); break;
        case calin::ix::iact_data::telescope_event::TRIGGER_UNKNOWN:
        case calin::ix::iact_data::telescope_event::TRIGGER_MULTIPLE:
        default:
          /* do nothing */
          break;
        };
      }
    }

    event_number_hist.dump_as_proto(results_->mutable_event_number_histogram());

    elapsed_time_hist.dump_as_proto(results_->mutable_elapsed_time_histogram());

    trigger_physics_elapsed_time_hist.dump_as_proto(
      results_->mutable_elapsed_time_histogram_trigger_physics());
    trigger_software_elapsed_time_hist.dump_as_proto(
      results_->mutable_elapsed_time_histogram_trigger_software());
    trigger_pedestal_elapsed_time_hist.dump_as_proto(
      results_->mutable_elapsed_time_histogram_trigger_pedestal());
    trigger_external_flasher_elapsed_time_hist.dump_as_proto(
      results_->mutable_elapsed_time_histogram_trigger_external_flasher());
    trigger_internal_flasher_elapsed_time_hist.dump_as_proto(
      results_->mutable_elapsed_time_histogram_trigger_internal_flasher());
    trigger_forced_array_elapsed_time_hist.dump_as_proto(
      results_->mutable_elapsed_time_histogram_trigger_forced_array());
    trigger_ucts_aux_elapsed_time_hist.dump_as_proto(
      results_->mutable_elapsed_time_histogram_trigger_ucts_aux());

    auto event_index = calin::util::algorithm::argsort(
      partials_->event_number_sequence().begin(), partials_->event_number_sequence().end());

    uint64_t first_event_number = partials_->event_number_sequence(event_index.front());
    uint64_t last_event_number = first_event_number;

    for(unsigned iel=1;iel<event_index.size();iel++) {
      uint64_t event_number = partials_->event_number_sequence(event_index[iel]);
      if(event_number == last_event_number) {
        results_->increment_num_duplicate_event_numbers();
        calin::diagnostics::range::encode_value(
          results_->mutable_duplicate_event_numbers(), event_number);
      } else if(event_number == last_event_number+1) {
        last_event_number = event_number;
      } else {
        results_->mutable_event_numbers_found()->add_begin_index(first_event_number);
        results_->mutable_event_numbers_found()->add_end_index(last_event_number + 1);
        first_event_number = last_event_number = event_number;
      }
    }

    results_->mutable_event_numbers_found()->add_begin_index(first_event_number);
    results_->mutable_event_numbers_found()->add_end_index(last_event_number + 1);

    std::vector<bool> event_value_bool(event_index.size());
    std::vector<int64_t> event_value_i64(event_index.size());

    make_index_range_from_conditional_bool_rle(
      event_index, partials_->event_number_sequence().begin(),
      partials_->cdts_presence(),
      event_value_bool, results_->mutable_events_missing_cdts(), false);

    make_index_range_from_conditional_bool_rle(
      event_index, partials_->event_number_sequence().begin(),
      partials_->tib_presence(),
      event_value_bool, results_->mutable_events_missing_tib(), false);

    make_index_range_from_conditional_bool_rle(
      event_index, partials_->event_number_sequence().begin(),
      partials_->swat_presence(),
      event_value_bool, results_->mutable_events_missing_swat(), false);

    make_index_range_from_conditional_bool_rle(
      event_index, partials_->event_number_sequence().begin(),
      partials_->all_channels_presence(),
      event_value_bool, results_->mutable_events_missing_modules(), false);

    for(int imod=0;imod<partials_->module_size();imod++) {
      auto* rimod = results_->mutable_module(imod);
      auto& pimod = partials_->module(imod);

      make_index_range_from_conditional_bool_rle(
        event_index, partials_->event_number_sequence().begin(),
        pimod.module_presence(),
        event_value_bool, results_->mutable_module(imod)->mutable_events_missing(), false);
      rimod->set_num_events_present(pimod.num_events_present());

      for(unsigned icounter=0; icounter<mod_counter_id_.size(); icounter++) {
        rimod->mutable_counter_value(icounter)->mutable_value_squared_sum_histogram()->IntegrateFrom(
          pimod.counter_value_squared_sum_histogram(icounter));
        if(config_.enable_module_counter_test_value_range()) {
          make_index_range_from_i64_value_rle(
            event_index, partials_->event_number_sequence().begin(),
            pimod.module_presence(), pimod.counter_value(icounter),
            event_value_bool, event_value_i64,
            rimod->mutable_counter_value(icounter)->mutable_value_range());
        }
      }
    }

    const double dthistmin = -9.0;
    const double dthistmax = 9.0;

    calin::math::histogram::Histogram1D log10_delta_t_hist {
      config_.log10_delta_t_histogram_binsize(), dthistmin, dthistmax, 0.0 };
    calin::math::histogram::Histogram1D log10_delta2_t_hist {
      config_.log10_delta_t_histogram_binsize(), dthistmin, dthistmax, 0.0 };

    calin::math::histogram::Histogram1D pt_log10_delta_t_hist {
      config_.log10_delta_t_histogram_binsize(), dthistmin, dthistmax, 0.0 };
    calin::math::histogram::Histogram1D pt_log10_delta2_t_hist {
      config_.log10_delta_t_histogram_binsize(), dthistmin, dthistmax, 0.0 };
    calin::math::histogram::Histogram1D pt2_log10_delta_t_hist {
      config_.log10_delta_t_histogram_binsize(), dthistmin, dthistmax, 0.0 };

    calin::math::histogram::Histogram1D rec_log10_delta_t_hist {
      config_.log10_delta_t_histogram_binsize(), dthistmin, dthistmax, 0.0 };

    uint64_t second_last_event_number = partials_->event_number_sequence(event_index[0]);
    int64_t  second_last_event_time = partials_->event_time_sequence(event_index[0]);

    last_event_number = partials_->event_number_sequence(event_index[0]);
    int64_t  last_event_time = partials_->event_time_sequence(event_index[0]);
    uint32_t last_event_type = partials_->event_type_sequence(event_index[0]);

    for(unsigned ievent=1; ievent<event_index.size(); ievent++) {
      uint64_t event_number = partials_->event_number_sequence(event_index[ievent]);
      int64_t  event_time = partials_->event_time_sequence(event_index[ievent]);
      uint32_t event_type = partials_->event_type_sequence(event_index[ievent]);

      double dt = double(event_time-last_event_time)*1e-9;
      double log10_dt = std::log10(std::max(dt, 0.0));

      if(dt>0) {
        rec_log10_delta_t_hist.insert(log10_dt);
      } else {
        results_->increment_num_delta_t_all_recorded_not_positive();
      }

      if(event_number == last_event_number+1 and event_time>0 and last_event_time>0) {
        double dt = double(event_time-last_event_time)*1e-9;
        if(dt>0) {
          log10_delta_t_hist.insert(log10_dt);
        }
        if(event_type == calin::ix::iact_data::telescope_event::TRIGGER_PHYSICS) {
          if(dt>0) {
            pt_log10_delta_t_hist.insert(log10_dt);
          }
        }

        if(event_number == second_last_event_number+2 and second_last_event_time>0) {
          double d2t = double(event_time-second_last_event_time)*1e-9;
          double log10_d2t = std::log10(std::max(d2t, 0.0));

          if(d2t>0) {
            log10_delta2_t_hist.insert(log10_d2t);
          }
          if(event_type == calin::ix::iact_data::telescope_event::TRIGGER_PHYSICS and
              last_event_type == calin::ix::iact_data::telescope_event::TRIGGER_PHYSICS) {
            if(dt>0) {
              pt2_log10_delta_t_hist.insert(log10_dt);
            }
            if(d2t>0) {
              pt_log10_delta2_t_hist.insert(log10_d2t);
            }
          }
        } else if(event_number == last_event_number+2 and event_time>0 and last_event_time>0) {
          double d2t = double(event_time-last_event_number)*1e-9;
          if(d2t>0) {
            log10_delta2_t_hist.insert(std::log10(d2t));
          }
        }
      }
      second_last_event_number = last_event_number;
      second_last_event_time = last_event_time;
      last_event_number = event_number;
      last_event_time = event_time;
      last_event_type = event_type;
    }

    log10_delta_t_hist.dump_as_proto(
      results_->mutable_log10_delta_t_histogram());
    log10_delta2_t_hist.dump_as_proto(
      results_->mutable_log10_delta2_t_histogram());
    pt_log10_delta_t_hist.dump_as_proto(
      results_->mutable_log10_delta_t_histogram_trigger_physics());
    pt_log10_delta2_t_hist.dump_as_proto(
      results_->mutable_log10_delta2_t_histogram_trigger_physics());
    pt2_log10_delta_t_hist.dump_as_proto(
      results_->mutable_log10_delta_t_histogram_2_trigger_physics());
    rec_log10_delta_t_hist.dump_as_proto(
      results_->mutable_log10_delta_t_histogram_all_recorded());
  }
}
