/*

   calin/iact_data/parallel_event_dispatcher.hpp -- Stephen Fegan -- 2018-02-08

   A parallel-only dispatcher of run and event data

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

#pragma once

#include <string>
#include <atomic>
#include <vector>
#include <map>

#include <calin_global_definitions.hpp>
#include <calin_global_config.hpp>
#include <iact_data/telescope_data_source.hpp>
#include <iact_data/event_visitor.hpp>
#include <iact_data/cta_data_source.hpp>
#include <iact_data/cta_acada_event_decoder.hpp>
#include <iact_data/event_dispatcher.pb.h>

namespace calin { namespace iact_data { namespace event_dispatcher {

class ParallelEventDispatcher: protected
  calin::iact_data::event_visitor::EventLifetimeManager
{
public:
  CALIN_TYPEALIAS(config_type,
    calin::ix::iact_data::event_dispatcher::EventDispatcherConfig);

  ParallelEventDispatcher();

  ~ParallelEventDispatcher();

  void add_visitor(
    calin::iact_data::event_visitor::ParallelEventVisitor* visitor,
    const std::string& processing_record_comment,
    bool adopt_visitor = false);

  void add_visitor(
      calin::iact_data::event_visitor::ParallelEventVisitor* visitor,
      bool adopt_visitor = false) {
    add_visitor(visitor, std::string{}, adopt_visitor);
  }

  void process_run(calin::iact_data::telescope_data_source::
    TelescopeRandomAccessDataSourceWithRunConfig* src,
    unsigned log_frequency = 0, int nthread = 0);

  void process_run(std::vector<calin::iact_data::telescope_data_source::
    TelescopeDataSourceWithRunConfig*> src_list,
    unsigned log_frequency = 0);

  void process_run(std::vector<calin::iact_data::telescope_data_source::
    TelescopeRandomAccessDataSourceWithRunConfig*> src_list,
    unsigned log_frequency = 0);

  void process_run(calin::io::data_source::DataSource<
      calin::ix::iact_data::telescope_event::TelescopeEvent>* src,
    calin::ix::iact_data::
      telescope_run_configuration::TelescopeRunConfiguration* run_config,
    unsigned log_frequency = 0, int nthread = 0);

  void process_run(std::vector<calin::io::data_source::DataSource<
      calin::ix::iact_data::telescope_event::TelescopeEvent>*> src_list,
    calin::ix::iact_data::
      telescope_run_configuration::TelescopeRunConfiguration* merged_run_config,
    unsigned log_frequency = 0);

  void process_cta_zfits_run(const std::string& filename,
    const calin::ix::iact_data::event_dispatcher::EventDispatcherConfig& config = default_config());

#if 0
  void process_cta_zmq_run(const std::vector<std::string>& endpoints,
    const calin::ix::iact_data::event_dispatcher::EventDispatcherConfig& config = default_config());

  void process_cta_zmq_run(const std::string& endpoint,
    const calin::ix::iact_data::event_dispatcher::EventDispatcherConfig& config = default_config());
#endif

  static calin::ix::iact_data::event_dispatcher::EventDispatcherConfig default_config();

  static bool merge_run_config(calin::ix::iact_data::
      telescope_run_configuration::TelescopeRunConfiguration* to,
    const calin::ix::iact_data::
      telescope_run_configuration::TelescopeRunConfiguration& from);

private:
#ifndef SWIG
  struct DelegatedVisitor {
    DelegatedVisitor(calin::iact_data::event_visitor::ParallelEventVisitor* _visitor,
        const std::string& _processing_record_comment):
      visitor(_visitor), processing_record_comment(_processing_record_comment)
    { /* nothing to see here */ }

    calin::iact_data::event_visitor::ParallelEventVisitor* visitor;
    std::string processing_record_comment;
    calin::ix::provenance::chronicle::ProcessingRecord* processing_record = nullptr;
  };

  void dispatch_run_configuration(calin::ix::iact_data::
    telescope_run_configuration::TelescopeRunConfiguration* run_config, bool register_processor);
  void dispatch_event(uint64_t seq_index,
    calin::ix::iact_data::telescope_event::TelescopeEvent* event);
  void dispatch_leave_run();
  void dispatch_merge_results();

  void do_parallel_dispatcher_loops(
    calin::ix::iact_data::
      telescope_run_configuration::TelescopeRunConfiguration* run_config,
    calin::io::data_source::DataSourceFactory<
      calin::ix::iact_data::telescope_event::TelescopeEvent>* src_factory,
    unsigned nthread, unsigned log_frequency,
    const std::chrono::system_clock::time_point& start_time,
    std::atomic<uint_fast64_t>& ndispatched);

  void do_dispatcher_loop(
    calin::io::data_source::DataSource<
      calin::ix::iact_data::telescope_event::TelescopeEvent>* src,
    unsigned log_frequency, unsigned nevents_to_dispatch,
    const std::chrono::system_clock::time_point& start_time,
    std::atomic<uint_fast64_t>& ndispatched);

  void write_initial_log_message(calin::ix::iact_data::
    telescope_run_configuration::TelescopeRunConfiguration* run_config, int nthread);

  void write_final_log_message(
    calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration* run_config,
    const std::chrono::system_clock::time_point& start_time,
    std::atomic<uint_fast64_t>& ndispatched);

  std::vector<calin::iact_data::event_visitor::ParallelEventVisitor*> adopted_visitors_;
  std::vector<DelegatedVisitor> visitors_;

  struct managed_event {
    unsigned usage_count;
    uint64_t seq_index;
    google::protobuf::Arena* arena;
  };

  std::map<calin::ix::iact_data::telescope_event::TelescopeEvent*, managed_event>
    event_keep_;
#endif

  void add_event_to_keep(
    calin::ix::iact_data::telescope_event::TelescopeEvent* event,
    uint64_t seq_index, google::protobuf::Arena* arena);
  void keep_event(
    calin::ix::iact_data::telescope_event::TelescopeEvent* event) override;
  void release_event(
    calin::ix::iact_data::telescope_event::TelescopeEvent* event) override;
};

} } } // namespace calin::iact_data::event_dispatcher
