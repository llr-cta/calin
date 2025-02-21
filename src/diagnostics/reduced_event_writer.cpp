/*

   calin/diagnostics/reduced_event_writer.cpp -- Stephen Fegan -- 2025-02-17

   Write charge data to a HDF5 file.

   Copyright 2025, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <algorithm>

#include <io/json.hpp>
#include <util/string.hpp>
#include <util/log.hpp>
#include <diagnostics/reduced_event_writer.hpp>

using namespace calin::util::log;
using namespace calin::ix::diagnostics::reduced_event_writer;
using namespace calin::iact_data::waveform_treatment_event_visitor;
using namespace calin::diagnostics::reduced_event_writer;

ReducedEventWriterParallelEventVisitor::
ReducedEventWriterParallelEventVisitor(
    calin::iact_data::waveform_treatment_event_visitor::OptimalWindowSumWaveformTreatmentParallelEventVisitor* gain1_visitor,
    calin::iact_data::waveform_treatment_event_visitor::OptimalWindowSumWaveformTreatmentParallelEventVisitor* gain2_visitor,
    calin::ix::diagnostics::reduced_event_writer::ReducedEventWriterConfig config,
    bool adopt_gain_visitors):
  ParallelEventVisitor(), config_(config), gain1_visitor_(gain1_visitor), gain2_visitor_(gain2_visitor),
  adopt_gain_visitors_(adopt_gain_visitors)
{
  //nothing to see here
}

ReducedEventWriterParallelEventVisitor::~ReducedEventWriterParallelEventVisitor()
{
  if(adopt_gain_visitors_) {
    delete gain1_visitor_;
    delete gain2_visitor_;
  }
}

ReducedEventWriterParallelEventVisitor* ReducedEventWriterParallelEventVisitor::new_sub_visitor(
  std::map<calin::iact_data::event_visitor::ParallelEventVisitor*,
    calin::iact_data::event_visitor::ParallelEventVisitor*> antecedent_visitors)
{
  auto* g1v = dynamic_cast<OptimalWindowSumWaveformTreatmentParallelEventVisitor*>(
    antecedent_visitors[gain1_visitor_]);
  auto* g2v = dynamic_cast<OptimalWindowSumWaveformTreatmentParallelEventVisitor*>(
    antecedent_visitors[gain2_visitor_]);
  auto* child = new ReducedEventWriterParallelEventVisitor(g1v, g2v, config_, false);
  child->parent_ = this;
  return child;
}

bool ReducedEventWriterParallelEventVisitor::visit_telescope_run(
  const calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration* run_config,
  calin::iact_data::event_visitor::EventLifetimeManager* event_lifetime_manager,
  calin::ix::provenance::chronicle::ProcessingRecord* processing_record)
{
  if(processing_record) {
    processing_record->set_type("ReducedEventWriterParallelEventVisitor");
    processing_record->set_description("Reduced event file writer");
    auto* config_json = processing_record->add_config();
    config_json->set_type(config_.GetTypeName());
    config_json->set_json(calin::io::json::encode_protobuf_to_json_string(config_));
    config_json = processing_record->add_config();
    std::vector<std::pair<std::string,std::string> > keyval;
    keyval.emplace_back("gain1WaveformSumInstance",
      calin::io::json::json_string_value(calin::util::string::instance_identifier(gain1_visitor_)));
    keyval.emplace_back("gain2WaveformSumInstance",
      calin::io::json::json_string_value(calin::util::string::instance_identifier(gain2_visitor_)));
    config_json->set_json(calin::io::json::json_for_dictionary(keyval));
  }

  run_config_.Clear();
  run_config_.CopyFrom(*run_config);

  if(parent_ == nullptr) {
    std::string filename;
    if(not config_.directory().empty()) {
      filename = config_.directory();
      if(filename.back() != '/') filename += '/';
    }
    if(not config_.file_prefix().empty()) {
      filename += config_.file_prefix();
      if(filename.back() != '_') filename += '_';
    }
    filename += calin::util::string::to_string(run_config->run_number()) + ".h5";

    auto* run_config_writer = run_config->NewHDFStreamWriter(filename, 
      config_.run_configuration_group(), config_.truncate());
    run_config_writer->write(*run_config);
    delete run_config_writer;

    zmq_push_pull_ = std::make_unique<calin::io::zmq_inproc::ZMQInprocPushPull>();

    writer_thread_ = std::make_unique<std::thread>([this,filename](){
      std::unique_ptr<calin::io::zmq_inproc::ZMQPuller> zmq_puller { 
        zmq_push_pull_->new_puller(calin::io::zmq_inproc::ZMQBindOrConnect::BIND) };
      std::unique_ptr<calin::ix::diagnostics::reduced_event::ReducedEvent_StreamWriter> writer { 
        calin::ix::diagnostics::reduced_event::ReducedEvent::NewHDFStreamWriter(filename, config_.event_group()) };
      while(zmq_puller->wait_for_data()) {
        std::pair<google::protobuf::Arena*, calin::ix::diagnostics::reduced_event::ReducedEvent*> arena_event;
        if(not zmq_puller->pull_assert_size(&arena_event, sizeof(arena_event),true))break;
        writer->write(*arena_event.second);
        delete(arena_event.first);
      }
    });

    zmq_pusher_.reset(zmq_push_pull_->new_pusher(calin::io::zmq_inproc::ZMQBindOrConnect::CONNECT));
  } else {
    zmq_pusher_.reset(parent_->zmq_push_pull_->new_pusher(calin::io::zmq_inproc::ZMQBindOrConnect::CONNECT));
  }

  return true;
}

bool ReducedEventWriterParallelEventVisitor::leave_telescope_run(
  calin::ix::provenance::chronicle::ProcessingRecord* processing_record)
{
  zmq_pusher_.reset();
  if(parent_ == nullptr) {
    zmq_push_pull_.reset();
    writer_thread_->join();
    writer_thread_.reset();
  }
  return true;
}

namespace {
  void copy_gain(calin::ix::diagnostics::reduced_event::OneGainChargeSums* dest,
    calin::iact_data::waveform_treatment_event_visitor::OptimalWindowSumWaveformTreatmentParallelEventVisitor* gain_visitor,
    ReducedEventWriterConfig& config) 
  {
    unsigned nchan = gain_visitor->nchan();
    dest->mutable_signal_type()->Reserve(nchan);
    dest->mutable_signal_type()->Add(
      gain_visitor->array_chan_signal_type(), gain_visitor->array_chan_signal_type()+nchan);
    if(config.write_max_sample()) {
      dest->mutable_max_sample()->Reserve(nchan);
      dest->mutable_max_sample()->Add(
        gain_visitor->array_chan_max(), gain_visitor->array_chan_max()+nchan);
    }
    if(config.write_max_index()) {
      dest->mutable_max_index()->Reserve(nchan);
      dest->mutable_max_index()->Add(
        gain_visitor->array_chan_max_index(), gain_visitor->array_chan_max_index()+nchan);
    }
    if(config.write_bkg_qsum()) {
      dest->mutable_bkg_qsum()->Reserve(nchan);
      dest->mutable_bkg_qsum()->Add(
        gain_visitor->array_chan_bkg_win_sum(), gain_visitor->array_chan_bkg_win_sum()+nchan);
    }
    if(config.write_sig_qsum()) {
      dest->mutable_sig_qsum()->Reserve(nchan);
      dest->mutable_sig_qsum()->Add(
        gain_visitor->array_chan_sig_win_sum(), gain_visitor->array_chan_sig_win_sum()+nchan);
    }
    if(config.write_opt_qsum()) {
      dest->mutable_opt_qsum()->Reserve(nchan);
      dest->mutable_opt_qsum()->Add(
        gain_visitor->array_chan_opt_win_sum(), gain_visitor->array_chan_opt_win_sum()+nchan);
    }
    if(config.write_opt_qtsum()) {
      dest->mutable_opt_qtsum()->Reserve(nchan);
      dest->mutable_opt_qtsum()->Add(
        gain_visitor->array_chan_opt_win_sum_qt(), gain_visitor->array_chan_opt_win_sum_qt()+nchan);
    }
    if(config.write_opt_index()) {
      dest->mutable_opt_index()->Reserve(nchan);
      dest->mutable_opt_index()->Add(
        gain_visitor->array_chan_opt_win_index(), gain_visitor->array_chan_opt_win_index()+nchan);
    }
    if(config.write_all_qsum()) {
      dest->mutable_all_qsum()->Reserve(nchan);
      dest->mutable_all_qsum()->Add(
        gain_visitor->array_chan_all_sum(), gain_visitor->array_chan_all_sum()+nchan);
    }
    if(config.write_opt_bkg_qsum_diff()) {
      dest->mutable_opt_bkg_qsum_diff()->Reserve(nchan);
      for(unsigned ichan=0;ichan<nchan;ichan++) {
        dest->add_opt_bkg_qsum_diff(gain_visitor->array_chan_opt_win_sum()[ichan] - gain_visitor->array_chan_bkg_win_sum()[ichan]);
      }
    }
  }
}

bool ReducedEventWriterParallelEventVisitor::visit_telescope_event(uint64_t seq_index,
  calin::ix::iact_data::telescope_event::TelescopeEvent* event)
{
  auto* arena = new google::protobuf::Arena();
  auto* reduced_event = google::protobuf::Arena::CreateMessage<calin::ix::diagnostics::reduced_event::ReducedEvent>(arena);

  reduced_event->set_local_event_number(event->local_event_number());
  reduced_event->set_trigger_type(
    static_cast<calin::ix::diagnostics::reduced_event::TriggerType>(event->trigger_type()));
  reduced_event->set_absolute_event_time_ns(event->absolute_event_time().time_ns());
  if(config_.write_gain1() and gain1_visitor_ and gain1_visitor_->is_same_event(seq_index)) {
    copy_gain(reduced_event->mutable_gain1(), gain1_visitor_, config_);
  }
  if(config_.write_gain2() and gain2_visitor_ and gain2_visitor_->is_same_event(seq_index)) {
    copy_gain(reduced_event->mutable_gain2(), gain2_visitor_, config_);
  }
  if(config_.write_l0_trigger_map() and event->has_trigger_map() and event->trigger_map().trigger_image_size() > 0) {
    auto* l0map = reduced_event->mutable_l0_trigger_map();
    l0map->mutable_trigger_hit()->Add(
      event->trigger_map().trigger_image().begin(), event->trigger_map().trigger_image().end());
  }

  std::pair<google::protobuf::Arena*, calin::ix::diagnostics::reduced_event::ReducedEvent*> arena_event{arena, reduced_event};
  zmq_pusher_->push(&arena_event, sizeof(arena_event));

  return true;
}

bool ReducedEventWriterParallelEventVisitor::merge_results()
{
  return true;
}

calin::ix::diagnostics::reduced_event_writer::ReducedEventWriterConfig ReducedEventWriterParallelEventVisitor::default_config()
{
  calin::ix::diagnostics::reduced_event_writer::ReducedEventWriterConfig config;
  config.set_write_gain1(true);
  config.set_write_gain2(true);
  config.set_write_l0_trigger_map(true);
  
  config.set_write_max_sample(true);
  config.set_write_max_index(true);
  config.set_write_bkg_qsum(true);
  config.set_write_sig_qsum(false);
  config.set_write_opt_qsum(false);
  config.set_write_opt_qtsum(false);
  config.set_write_opt_index(true);
  config.set_write_all_qsum(false);
  config.set_write_opt_bkg_qsum_diff(true);
  
  config.set_run_configuration_group("run_configuration");
  config.set_event_group("events");
  config.set_truncate(true);
  config.set_file_prefix("reduced_event"); 
  return config;
}
