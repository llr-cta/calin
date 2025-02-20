/*

   calin/diagnostics/reduced_file_writer.hpp -- Stephen Fegan -- 2025-02-20

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

#include <map>
#include <Eigen/Dense>

#include <iact_data/event_visitor.hpp>
#include <iact_data/waveform_treatment_event_visitor.hpp>

namespace calin { namespace diagnostics { namespace reduced_file_writer {

class SimpleChargeCaptureParallelEventVisitor:
  public calin::iact_data::event_visitor::ParallelEventVisitor
{
public:
  SimpleChargeCaptureParallelEventVisitor(
    calin::iact_data::waveform_treatment_event_visitor::OptimalWindowSumWaveformTreatmentParallelEventVisitor* waveform_sum_visitor,
    ChargeCaptureDatum capture_datum = CCD_OPT_WINDOW_SUM,
    unsigned max_event_number = 0, unsigned min_event_number = 0,
    bool adopt_waveform_sum_visitor = false);

  virtual ~SimpleChargeCaptureParallelEventVisitor();

  SimpleChargeCaptureParallelEventVisitor* new_sub_visitor(
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
  uint64_t size() const { return captured_data_.size(); }
  std::vector<uint64_t> keys() const;
  bool has(uint64_t event_number) const;
  int64_t time(uint64_t event_number) const;
  int32_t trigger_code(uint64_t event_number) const;
  Eigen::VectorXi signal_types(uint64_t event_number) const;
  Eigen::VectorXi values(uint64_t event_number) const;

private:
#ifndef SWIG
  struct CapturedEventData {
    uint64_t event_number;
    int64_t time;
    int32_t trigger_code;
    std::vector<int32_t> values;
  };

  SimpleChargeCaptureParallelEventVisitor* parent_ = nullptr;

  calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration run_config_;

  calin::iact_data::waveform_treatment_event_visitor::OptimalWindowSumWaveformTreatmentParallelEventVisitor* waveform_sum_visitor_ = nullptr;
  bool adopt_waveform_sum_visitor_ = false;
  std::map<uint64_t, CapturedEventData*> captured_data_;
  ChargeCaptureDatum capture_datum_;
  unsigned max_event_number_;
  unsigned min_event_number_;
#endif
};

} } } // namespace calin::diagnostics::simple_charge_capture
