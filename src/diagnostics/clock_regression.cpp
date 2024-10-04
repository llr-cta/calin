/*

   calin/diagnostics/clock_regression.cpp -- Stephen Fegan -- 2020-05-25

   Clock regression visitor

   Copyright 2020, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <Eigen/Dense>
#include <diagnostics/clock_regression.hpp>
#include <io/json.hpp>
#include <util/log.hpp>

using namespace calin::util::log;

using namespace calin::diagnostics::clock_regression;

ClockRegressionParallelEventVisitor::
ClockRegressionParallelEventVisitor(
    const calin::ix::diagnostics::clock_regression::ClockRegressionConfig& config):
  calin::iact_data::event_visitor::ParallelEventVisitor(), config_(config)
{
  // nothing to see here
}

ClockRegressionParallelEventVisitor::~ClockRegressionParallelEventVisitor()
{
  // nothing to see here
}

ClockRegressionParallelEventVisitor*
ClockRegressionParallelEventVisitor::new_sub_visitor(
  std::map<calin::iact_data::event_visitor::ParallelEventVisitor*,
    calin::iact_data::event_visitor::ParallelEventVisitor*> antecedent_visitors)
{
  auto* sub_visior = new ClockRegressionParallelEventVisitor(config_);
  sub_visior->parent_ = this;
  return sub_visior;
}

namespace {
  int find_clock_id_by_name(const std::string& name, 
    const google::protobuf::RepeatedPtrField<std::string>& clock_names)
  {    
    for(int clock_id=0;clock_id<clock_names.size();++clock_id) {
      if(clock_names[clock_id] == name)return clock_id;
    }
    return -1;
  }
}

bool ClockRegressionParallelEventVisitor::visit_telescope_run(
  const calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration* run_config,
  calin::iact_data::event_visitor::EventLifetimeManager* event_lifetime_manager,
  calin::ix::provenance::chronicle::ProcessingRecord* processing_record)
{
  if(processing_record) {
    processing_record->set_type("ClockRegressionParallelEventVisitor");
    processing_record->set_description("Clock regression");
    auto* config_json = processing_record->add_config();
    config_json->set_type(config_.GetTypeName());
    config_json->set_json(calin::io::json::encode_protobuf_to_json_string(config_));
  }

  rebalance_ = config_.rebalance_nevent();

  camera_tests_.clear();
  module_tests_.clear();

  const google::protobuf::RepeatedPtrField<calin::ix::diagnostics::
    clock_regression::SingleClockRegressionConfig>* cam_tests = nullptr;
  const google::protobuf::RepeatedPtrField<calin::ix::diagnostics::
    clock_regression::SingleClockRegressionConfig>* mod_tests = nullptr;
  if(config_.camera_clocks_size() > 0) {
    cam_tests = &config_.camera_clocks();
  }
  if(config_.module_clocks_size() > 0) {
    mod_tests = &config_.module_clocks();
  }

  switch(run_config->camera_layout().camera_type()) {
  case calin::ix::iact_data::instrument_layout::CameraLayout::NECTARCAM:
  case calin::ix::iact_data::instrument_layout::CameraLayout::NECTARCAM_TESTBENCH_19CHANNEL:
  case calin::ix::iact_data::instrument_layout::CameraLayout::NECTARCAM_TESTBENCH_61CHANNEL:
    if(cam_tests == nullptr and config_.default_nectarcam_camera_clocks_size() > 0) {
      cam_tests = &config_.default_nectarcam_camera_clocks();
    }
    if(mod_tests == nullptr and config_.default_nectarcam_module_clocks_size() > 0) {
      mod_tests = &config_.default_nectarcam_module_clocks();
    }
    break;
  case calin::ix::iact_data::instrument_layout::CameraLayout::LSTCAM:
    if(cam_tests == nullptr and config_.default_lstcam_camera_clocks_size() > 0) {
      cam_tests = &config_.default_lstcam_camera_clocks();
    }
    if(mod_tests == nullptr and config_.default_lstcam_module_clocks_size() > 0) {
      mod_tests = &config_.default_lstcam_module_clocks();
    }
    break;
  default:
    break;
  }

  if(cam_tests) {
    for(const auto& cam_test : (*cam_tests)) {
      ClockTest ct;
      ct.config.CopyFrom(cam_test);
      if(not ct.config.clock_name().empty()) {
        ct.config.set_clock_id(find_clock_id_by_name(ct.config.clock_name(), 
          run_config->camera_layout().camera_clock_name()));
      } else if(ct.config.clock_id() < run_config->camera_layout().camera_clock_name_size()) {
        ct.config.set_clock_name(run_config->camera_layout().camera_clock_name(ct.config.clock_id()));
      }
      if(ct.config.clock_id()<0 
          or ct.config.clock_id()>=run_config->camera_layout().camera_clock_name_size())
        continue;
      camera_tests_.push_back(ct);
    }
  }

  if(mod_tests) {
    for(const auto& mod_test : (*mod_tests)) {
      ModuleClockTest mct;
      mct.config = mod_test;
      if(not mct.config.clock_name().empty()) {
        mct.config.set_clock_id(find_clock_id_by_name(mct.config.clock_name(), 
          run_config->camera_layout().module_clock_name()));
      } else if(mct.config.clock_id() < run_config->camera_layout().camera_clock_name_size()) {
        mct.config.set_clock_name(run_config->camera_layout().camera_clock_name(mct.config.clock_id()));
      }
      mct.modules.resize(run_config->configured_module_id_size());
      if(mct.config.clock_id()<0 
          or mct.config.clock_id()>=run_config->camera_layout().camera_clock_name_size())
        continue;
      // for(auto& ct : mct.modules) {
      //   ct.config = mct.config;
      // }
      module_tests_.push_back(mct);
    }
  }

  return true;
}

bool ClockRegressionParallelEventVisitor::leave_telescope_run(
  calin::ix::provenance::chronicle::ProcessingRecord* processing_record)
{
  return true;
}

namespace {
  const calin::ix::iact_data::telescope_event::Clock* find_clock(int id,
    const google::protobuf::RepeatedPtrField<calin::ix::iact_data::telescope_event::Clock>& clocks)
  {
    if(id < clocks.size() and clocks[id].clock_id() == id) {
      return &clocks[id];
    }
    for(auto& iclock : clocks) {
      if(iclock.clock_id() == id) {
        return &iclock;
      }
    }
    return nullptr;
  }

  void accumulate_clock(int64_t principal_time, uint64_t local_event_number,
    const calin::ix::iact_data::telescope_event::Clock& clock,
    const calin::ix::diagnostics::clock_regression::SingleClockRegressionConfig& config,
    std::map<int, calin::diagnostics::clock_regression::ClockRegressionParallelEventVisitor::RegressionAccumulator*>& bins,
    bool do_rebalance)
  {
    int ibin;
    switch(config.partition_mode()) {
    case calin::ix::diagnostics::clock_regression::PARTITION_BY_CLOCK_SEQUENCE_ID:
      ibin = clock.time_sequence_id();
      break;
    case calin::ix::diagnostics::clock_regression::PARTITION_BY_LOCAL_EVENT_NUMBER:
      ibin = local_event_number/config.partition_bin_size();
      break;
    case calin::ix::diagnostics::clock_regression::PARTITION_BY_MASTER_CLOCK:
      ibin = principal_time/config.partition_bin_size();
      break;
    case calin::ix::diagnostics::clock_regression::SINGLE_PARTITION:
    default:
      ibin = 0;
      break;
    }
    auto* accumulator = bins[ibin];
    if(accumulator == nullptr) {
      bins[ibin] = accumulator =
        new calin::diagnostics::clock_regression::ClockRegressionParallelEventVisitor::RegressionAccumulator();
    }
    accumulator->accumulate(principal_time, clock.time_value());
    if(do_rebalance) {
      accumulator->rebalance();
    }
  }
}

bool ClockRegressionParallelEventVisitor::visit_telescope_event(uint64_t seq_index,
  calin::ix::iact_data::telescope_event::TelescopeEvent* event)
{
  bool do_rebalance = false;
  if(rebalance_ > 0) {
    --rebalance_;
    if(rebalance_ == 0) {
      rebalance_ = config_.rebalance_nevent();
      do_rebalance = true;
    }
  }

  const auto* principal_clock = find_clock(config_.principal_clock_id(),event->camera_clock());
  if(principal_clock == nullptr) {
    return true;
  }

  if(principal_clock->time_value_may_be_suspect() and
      config_.include_possibly_suspect_principal_time_values() == false) {
    return true;
  }

  for(auto& ct : camera_tests_)
  {
    const auto* test_clock = find_clock(ct.config.clock_id(),event->camera_clock());
    if(test_clock) {
      if(test_clock->time_value_may_be_suspect() and
          ct.config.include_possibly_suspect_time_values() == false) {
        continue;
      }
      int64_t principal_time = principal_clock->time_value();
      if(ct.config.principal_clock_divisor() > 1) {
        principal_time /= ct.config.principal_clock_divisor();
      }
      accumulate_clock(principal_time, event->local_event_number(), *test_clock,
        ct.config, ct.bins, do_rebalance);
    }
  }

  for(auto& mt : module_tests_)
  {
    int64_t principal_time = principal_clock->time_value();
    if(mt.config.principal_clock_divisor() > 1) {
      principal_time /= mt.config.principal_clock_divisor();
    }

    for(auto& imod : event->module_clock()) {
      const auto* test_clock = find_clock(mt.config.clock_id(),imod.clock());
      if(test_clock) {
        if(test_clock->time_value_may_be_suspect() and
            mt.config.include_possibly_suspect_time_values() == false) {
          continue;
        }
        accumulate_clock(principal_time, event->local_event_number(), *test_clock,
          mt.config, mt.modules[imod.module_id()].bins, do_rebalance);
      }
    }
  }

  return true;
}

void ClockRegressionParallelEventVisitor::merge_into(ClockTest* to, const ClockTest& from)
{
  for(auto& ibin : from.bins) {
    if(to->bins[ibin.first] == nullptr) {
      to->bins[ibin.first] = new calin::diagnostics::clock_regression::ClockRegressionParallelEventVisitor::RegressionAccumulator();
    }
    ibin.second->integrate_into(*to->bins[ibin.first]);
    if(config_.rebalance_nevent() > 0) {
      to->bins[ibin.first]->rebalance();
    }
  }
}

bool ClockRegressionParallelEventVisitor::merge_results()
{
  if(parent_) {
    for(unsigned ict=0; ict<camera_tests_.size(); ict++) {
      merge_into(&parent_->camera_tests_[ict], camera_tests_[ict]);
    }
    for(unsigned imt=0; imt<module_tests_.size(); imt++) {
      for(unsigned imod=0; imod<module_tests_[imt].modules.size(); imod++) {
        merge_into(&parent_->module_tests_[imt].modules[imod],
          module_tests_[imt].modules[imod]);
      }
    }
  }
  return true;
}

calin::ix::diagnostics::clock_regression::ClockRegressionConfig
ClockRegressionParallelEventVisitor::default_config()
{
  calin::ix::diagnostics::clock_regression::ClockRegressionConfig config;
  config.set_principal_clock_id(0); // UCTS timestamp
  config.set_rebalance_nevent(1000);

  // NectarCAM

  auto* clock = config.add_default_nectarcam_camera_clocks();
  clock->set_clock_name("UCTS 10MHz counter"); // UCTS 10MHz
  clock->set_partition_mode(calin::ix::diagnostics::clock_regression::PARTITION_BY_CLOCK_SEQUENCE_ID);

  clock = config.add_default_nectarcam_camera_clocks();
  clock->set_clock_name("TIB 10MHz counter"); // TIB 10MHz
  clock->set_partition_mode(calin::ix::diagnostics::clock_regression::PARTITION_BY_CLOCK_SEQUENCE_ID);

  clock = config.add_default_nectarcam_camera_clocks();
  clock->set_clock_name("FEB local 2ns TDC counter sum"); // FEB local oscillator sum
  clock->set_partition_mode(calin::ix::diagnostics::clock_regression::PARTITION_BY_CLOCK_SEQUENCE_ID);

  clock = config.add_default_nectarcam_camera_clocks();
  clock->set_clock_name("FEB local 2ns TDC counter sum"); // FEB local oscillator sum
  clock->set_partition_mode(calin::ix::diagnostics::clock_regression::PARTITION_BY_MASTER_CLOCK);
  clock->set_partition_bin_size(250000000);

  clock = config.add_default_nectarcam_camera_clocks();
  clock->set_clock_name("UCTS timestamp"); // UCTS timestamp
  clock->set_partition_mode(calin::ix::diagnostics::clock_regression::PARTITION_BY_LOCAL_EVENT_NUMBER);
  clock->set_partition_bin_size(10000);

  clock = config.add_default_nectarcam_module_clocks();
  clock->set_clock_name("local 2ns TDC time"); // local ~2ns TDC time
  clock->set_partition_mode(calin::ix::diagnostics::clock_regression::PARTITION_BY_CLOCK_SEQUENCE_ID);

  // LSTCAM

  clock = config.add_default_lstcam_camera_clocks();
  clock->set_clock_name("UCTS 10MHz counter"); // UCTS 10MHz
  clock->set_partition_mode(calin::ix::diagnostics::clock_regression::PARTITION_BY_CLOCK_SEQUENCE_ID);

  clock = config.add_default_lstcam_camera_clocks();
  clock->set_clock_name("TIB 10MHz counter"); // TIB 10MHz
  clock->set_partition_mode(calin::ix::diagnostics::clock_regression::PARTITION_BY_CLOCK_SEQUENCE_ID);

  clock = config.add_default_lstcam_module_clocks();
  clock->set_clock_name("backplane 10MHz counter"); // backplane 10MHz counter
  clock->set_partition_mode(calin::ix::diagnostics::clock_regression::SINGLE_PARTITION);

  clock = config.add_default_lstcam_module_clocks();
  clock->set_clock_name("local 133MHz oscillator counter"); // local 133MHz TDC time
  clock->set_partition_mode(calin::ix::diagnostics::clock_regression::PARTITION_BY_CLOCK_SEQUENCE_ID);

  return config;
}

void ClockRegressionParallelEventVisitor::transfer_clock_results(
  calin::ix::diagnostics::clock_regression::SingleClockRegressionResults* res,
  const ClockTest& ct) const
{
  res->set_clock_name(ct.config.clock_name());
  res->set_clock_id(ct.config.clock_id());
  for(const auto& ibin : ct.bins) {
    const auto* reg = &ibin.second->accumulator();

    if(reg->num_entries()) {
      auto& reg_param = (*res->mutable_bins())[ibin.first];
      reg_param.set_x0(reg->x0());
      reg_param.set_y0(reg->y0());
      reg_param.set_num_entries(reg->num_entries());
      double a;
      double b;
      double d2;
      reg->fit_parameters_and_d2(a,b,d2);
      reg_param.set_a(a);
      reg_param.set_b(b);
      reg_param.set_d2(d2);
    }
  }
}

calin::ix::diagnostics::clock_regression::ClockRegressionResults*
ClockRegressionParallelEventVisitor::clock_regression(
  calin::ix::diagnostics::clock_regression::ClockRegressionResults* results) const
{
  if(results) {
    results->Clear();
  } else {
    results = new calin::ix::diagnostics::clock_regression::ClockRegressionResults();
  }

  results->set_principal_clock_name(config_.principal_clock_name());
  results->set_principal_clock_id(config_.principal_clock_id());

  for(const auto& ct : camera_tests_) {
    transfer_clock_results(results->add_camera_clock(), ct);
  }

  for(const auto& mt : module_tests_) {
    auto* mod_res = results->add_module_clock();
    for(const auto& ct : mt.modules) {
      transfer_clock_results(mod_res->add_modules(), ct);
    }
  }

  return results;
}
