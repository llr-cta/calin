/*

   calin/diagnostics/reduced_file_writer.cpp -- Stephen Fegan -- 2025-02-17

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
using namespace calin::diagnostics::reduced_file_writer;

ReducedFileWriterParallelEventVisitor::
ReducedFileWriterParallelEventVisitor(
    calin::iact_data::waveform_treatment_event_visitor::OptimalWindowSumWaveformTreatmentParallelEventVisitor* gain1_visitor,
    calin::iact_data::waveform_treatment_event_visitor::OptimalWindowSumWaveformTreatmentParallelEventVisitor* gain2_visitor,
    calin::ix::diagnostics::reduced_event_writer::ReducedEventWriterConfig config,
    bool adopt_gain_visitors):
  ParallelEventVisitor(), config_(config), gain1_visitor_(gain1_visitor), gain2_visitor_(gain2_visitor),
  adopt_gain_visitors_(adopt_gain_visitors)
{
  //nothing to see here
}

ReducedFileWriterParallelEventVisitor::~ReducedFileWriterParallelEventVisitor()
{
  if(adopt_gain_visitors_) {
    delete gain1_visitor_;
    delete gain2_visitor_;
  }
}

ReducedFileWriterParallelEventVisitor* ReducedFileWriterParallelEventVisitor::new_sub_visitor(
  std::map<calin::iact_data::event_visitor::ParallelEventVisitor*,
    calin::iact_data::event_visitor::ParallelEventVisitor*> antecedent_visitors)
{
  auto* g1v = dynamic_cast<OptimalWindowSumWaveformTreatmentParallelEventVisitor*>(
    antecedent_visitors[gain1_visitor_]);
  auto* g2v = dynamic_cast<OptimalWindowSumWaveformTreatmentParallelEventVisitor*>(
    antecedent_visitors[gain2_visitor_]);
  auto* child = new ReducedFileWriterParallelEventVisitor(g1v, g2v, config_, false);
  child->parent_ = this;
  return child;
}

bool ReducedFileWriterParallelEventVisitor::visit_telescope_run(
  const calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration* run_config,
  calin::iact_data::event_visitor::EventLifetimeManager* event_lifetime_manager,
  calin::ix::provenance::chronicle::ProcessingRecord* processing_record)
{
  if(processing_record) {
    processing_record->set_type("ReducedFileWriterParallelEventVisitor");
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

    auto* run_config_writer = run_config->NewHDFStreamWriter(filename, 
      config_.run_configuration_group(), config_.truncate());
    run_config_writer->write(*run_config);
    delete run_config_writer;

    event_writer_.reset(calin::ix::diagnostics::reduced_event::ReducedEvent::NewHDFStreamWriter(filename, config_.event_group()));
  }

  return true;
}

bool ReducedFileWriterParallelEventVisitor::leave_telescope_run(
  calin::ix::provenance::chronicle::ProcessingRecord* processing_record)
{
  event_writer_.reset();
  return true;
}

namespace {
  void copy_gain(calin::ix::diagnostics::reduced_event::OneGainChargeSums* dest,
    calin::iact_data::waveform_treatment_event_visitor::OptimalWindowSumWaveformTreatmentParallelEventVisitor* gain_visitor,
    ReducedEventWriterConfig& config) 
  {
    dest->mutable_signal_type()->Assign(
      gain_visitor->chan_signal_type().begin(), gain_visitor->chan_signal_type().end());
    if(config.write_max_sample()) {
      dest->mutable_max_sample()->Assign(
        gain_visitor->chan_max().begin(), gain_visitor->chan_max().end());
    }
    if(config.write_max_index()) {
      dest->mutable_max_index()->Assign(
        gain_visitor->chan_max_index().begin(), gain_visitor->chan_max_index().end());
    }
    if(config.write_bkg_win_qsum()) {
      dest->mutable_bkg_win_qsum()->Assign(
        gain_visitor->chan_bkg_win_sum().begin(), gain_visitor->chan_bkg_win_sum().end());
    }
    if(config.write_sig_win_qsum()) {
      dest->mutable_sig_win_qsum()->Assign(
        gain_visitor->chan_sig_win_sum().begin(), gain_visitor->chan_sig_win_sum().end());
    }
    if(config.write_opt_win_qsum()) {
      dest->mutable_opt_win_qsum()->Assign(
        gain_visitor->chan_opt_win_sum().begin(), gain_visitor->chan_opt_win_sum().end());
    }
    if(config.write_opt_win_qtsum()) {
      dest->mutable_opt_win_qtsum()->Assign(
        gain_visitor->chan_opt_win_sum_qt().begin(), gain_visitor->chan_opt_win_sum_qt().end());
    }
    if(config.write_opt_win_index()) {
      dest->mutable_opt_win_index()->Assign(
        gain_visitor->chan_opt_win_index().begin(), gain_visitor->chan_opt_win_index().end());
    }
    if(config.write_all_win_qsum()) {
      dest->mutable_all_win_qsum()->Assign(
        gain_visitor->chan_all_sum().begin(), gain_visitor->chan_all_sum().end());
    }
  }
}

bool ReducedFileWriterParallelEventVisitor::visit_telescope_event(uint64_t seq_index,
  calin::ix::iact_data::telescope_event::TelescopeEvent* event)
{
  calin::ix::diagnostics::reduced_event::ReducedEvent reduced_event;
  reduced_event.set_local_event_number(event->local_event_number());
  reduced_event.set_trigger_type(
    static_cast<calin::ix::diagnostics::reduced_event::TriggerType>(event->trigger_type()));
  reduced_event.set_absolute_event_time_ns(event->absolute_event_time().time_ns());
  if(config_.write_gain1() and gain1_visitor_ and gain1_visitor_->is_same_event(seq_index)) {
    copy_gain(reduced_event.mutable_gain1(), gain1_visitor_, config_);
  }
  if(config_.write_gain2() and gain2_visitor_ and gain2_visitor_->is_same_event(seq_index)) {
    copy_gain(reduced_event.mutable_gain2(), gain2_visitor_, config_);
  }

  if(parent_) {
    std::lock_guard<std::mutex> lock(event_writer_mutex_);
    parent_->event_writer_->write(reduced_event);
  } else {
    event_writer_->write(reduced_event);
  }

  return true;
}

bool ReducedFileWriterParallelEventVisitor::merge_results()
{
  return true;
}
