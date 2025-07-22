/*

   calin/diagnostics/reduced_event_writer.hpp -- Stephen Fegan -- 2025-02-20

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

#pragma once

#include <thread>
#include <memory>

#include <io/zmq_inproc_push_pull.hpp>
#include <iact_data/event_visitor.hpp>
#include <iact_data/waveform_treatment_event_visitor.hpp>

#include <diagnostics/reduced_event.pb.h>
#include <diagnostics/reduced_event_writer.pb.h>

namespace calin { namespace diagnostics { namespace reduced_event_writer {

class ReducedEventWriterParallelEventVisitor:
  public calin::iact_data::event_visitor::ParallelEventVisitor
{
public:
  ReducedEventWriterParallelEventVisitor(
    calin::iact_data::waveform_treatment_event_visitor::OptimalWindowSumWaveformTreatmentParallelEventVisitor* gain1_visitor,
    calin::iact_data::waveform_treatment_event_visitor::OptimalWindowSumWaveformTreatmentParallelEventVisitor* gain2_visitor = nullptr,
    calin::ix::diagnostics::reduced_event_writer::ReducedEventWriterConfig config = default_config(),
    bool adopt_gain_visitors = false);

  virtual ~ReducedEventWriterParallelEventVisitor();

  ReducedEventWriterParallelEventVisitor* new_sub_visitor(
    std::map<calin::iact_data::event_visitor::ParallelEventVisitor*,
        calin::iact_data::event_visitor::ParallelEventVisitor*>
      antecedent_visitors = { }) override;

  bool visit_telescope_run(
    const calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration* run_config,
    calin::iact_data::event_visitor::EventLifetimeManager* event_lifetime_manager,
    calin::ix::provenance::chronicle::ProcessingRecord* processing_record = nullptr) override;

  bool leave_telescope_run(
    calin::ix::provenance::chronicle::ProcessingRecord* processing_record = nullptr) override;

  bool visit_telescope_event(uint64_t seq_index,
    calin::ix::iact_data::telescope_event::TelescopeEvent* event) override;

  bool merge_results() override;

  const calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration& run_config() const {
    return run_config_;
  }

  calin::ix::diagnostics::reduced_event_writer::ReducedEventWriterConfig config() const { return config_; }
  static calin::ix::diagnostics::reduced_event_writer::ReducedEventWriterConfig default_config();

private:
#ifndef SWIG
  calin::ix::diagnostics::reduced_event_writer::ReducedEventWriterConfig config_;
  ReducedEventWriterParallelEventVisitor* parent_ = nullptr;

  calin::iact_data::waveform_treatment_event_visitor::OptimalWindowSumWaveformTreatmentParallelEventVisitor* gain1_visitor_ = nullptr;
  calin::iact_data::waveform_treatment_event_visitor::OptimalWindowSumWaveformTreatmentParallelEventVisitor* gain2_visitor_ = nullptr;
  bool adopt_gain_visitors_ = false;

  std::unique_ptr<std::thread> writer_thread_;
  std::unique_ptr<calin::io::zmq_inproc::ZMQPusher> zmq_pusher_;
  std::unique_ptr<calin::io::zmq_inproc::ZMQInprocPushPull> zmq_push_pull_;
  std::unique_ptr<calin::ix::diagnostics::reduced_event::ReducedEvent_StreamWriter> event_writer_;

  calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration run_config_;
#endif
};

} } } // namespace calin::diagnostics::reduced_event_writer
