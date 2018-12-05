/*

   calin/iact_data/parallel_event_dispatcher.hpp -- Stephen Fegan -- 2018-02-08

   A parallel-only dispatcher of run and event data

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
#include <iact_data/cta_actl_event_decoder.hpp>
#include <iact_data/event_dispatcher.pb.h>
//#include <iact_data/cta_actl_event_decoder.hpp>

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
    bool adopt_visitor = false);

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

#ifdef CALIN_HAVE_CTA_CAMERASTOACTL
  void process_cta_zfits_run(const std::string& filename,
    const calin::ix::iact_data::event_dispatcher::EventDispatcherConfig& config = default_config());

  void process_cta_zmq_run(const std::vector<std::string>& endpoints,
    const calin::ix::iact_data::event_dispatcher::EventDispatcherConfig& config = default_config());

  void process_cta_zmq_run(const std::string& endpoint,
    const calin::ix::iact_data::event_dispatcher::EventDispatcherConfig& config = default_config());
#endif

#if 0
  void process_nectarcam_zfits_run(const std::string& filename,
    unsigned log_frequency = 0, int nthread = 0,
    const calin::ix::iact_data::nectarcam_data_source::NectarCamCameraEventDecoderConfig& decoder_config =
      calin::iact_data::nectarcam_data_source::NectarCamZFITSDataSource_L0::default_decoder_config(),
    const calin::ix::iact_data::zfits_data_source::ZFITSDataSourceConfig& zfits_config =
      calin::iact_data::nectarcam_data_source::NectarCamZFITSDataSource_L0::default_config());

  void process_nectarcam_zfits_run(const std::string& filename,
    unsigned log_frequency, int nthread,
    const calin::ix::iact_data::zfits_data_source::ZFITSDataSourceConfig& zfits_config,
    const calin::ix::iact_data::nectarcam_data_source::NectarCamCameraEventDecoderConfig& decoder_config =
      calin::iact_data::nectarcam_data_source::NectarCamZFITSDataSource_L0::default_decoder_config()) {
    return process_nectarcam_zfits_run(filename, log_frequency, nthread,
      decoder_config, zfits_config);
  }

  void process_nectarcam_zfits_run(const std::string& filename,
    const calin::ix::iact_data::nectarcam_data_source::NectarCamCameraEventDecoderConfig& decoder_config,
    unsigned log_frequency = 0, int nthread = 0,
    const calin::ix::iact_data::zfits_data_source::ZFITSDataSourceConfig& zfits_config =
      calin::iact_data::nectarcam_data_source::NectarCamZFITSDataSource_L0::default_config()) {
    return process_nectarcam_zfits_run(filename, log_frequency, nthread,
      decoder_config, zfits_config);
  }
#endif

  static calin::ix::iact_data::event_dispatcher::EventDispatcherConfig default_config();

  static bool merge_run_config(calin::ix::iact_data::
      telescope_run_configuration::TelescopeRunConfiguration* to,
    const calin::ix::iact_data::
      telescope_run_configuration::TelescopeRunConfiguration& from);

private:
  // These functions allow events to be passed on to the visitors - they
  // are not meant to be called directly as the visiors expect them to be
  // called in a specific order. They are liable to be made private.
  void dispatch_run_configuration(calin::ix::iact_data::
    telescope_run_configuration::TelescopeRunConfiguration* run_config,
    bool dispatch_only_to_adopted_visitors = false);
  void dispatch_event(uint64_t seq_index,
    calin::ix::iact_data::telescope_event::TelescopeEvent* event);
  void dispatch_leave_run(bool dispatch_only_to_adopted_visitors = false);
  void dispatch_merge_results(bool dispatch_only_to_adopted_visitors = true);

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
    unsigned log_frequency,
    const std::chrono::system_clock::time_point& start_time,
    std::atomic<uint_fast64_t>& ndispatched);

  void write_final_log_message(
    unsigned log_frequency, const std::chrono::system_clock::time_point& start_time,
    std::atomic<uint_fast64_t>& ndispatched);

  std::vector<calin::iact_data::event_visitor::ParallelEventVisitor*>
    adopted_visitors_;
  std::vector<calin::iact_data::event_visitor::ParallelEventVisitor*>
    visitors_;

  struct managed_event {
    unsigned usage_count;
    uint64_t seq_index;
    google::protobuf::Arena* arena;
  };

  void add_event_to_keep(
    calin::ix::iact_data::telescope_event::TelescopeEvent* event,
    uint64_t seq_index, google::protobuf::Arena* arena);
  void keep_event(
    calin::ix::iact_data::telescope_event::TelescopeEvent* event) override;
  void release_event(
    calin::ix::iact_data::telescope_event::TelescopeEvent* event) override;

  std::map<calin::ix::iact_data::telescope_event::TelescopeEvent*, managed_event>
    event_keep_;
};

} } } // namespace calin::iact_data::event_dispatcher
