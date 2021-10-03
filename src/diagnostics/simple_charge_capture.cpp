/*

   calin/diagnostics/simple_charge_stats.cpp -- Stephen Fegan -- 2020-03-22

   Channel info diagnostics

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

#include <algorithm>

#include <util/log.hpp>
#include <diagnostics/simple_charge_capture.hpp>
#include <io/json.hpp>
#include <util/string.hpp>

using namespace calin::util::log;
using namespace calin::diagnostics::simple_charge_capture;
using namespace calin::iact_data::waveform_treatment_event_visitor;

SimpleChargeCaptureParallelEventVisitor::
SimpleChargeCaptureParallelEventVisitor(
    calin::iact_data::waveform_treatment_event_visitor::OptimalWindowSumWaveformTreatmentParallelEventVisitor* waveform_sum_visitor,
    ChargeCaptureDatum capture_datum, unsigned max_event_number, unsigned min_event_number, bool adopt_waveform_sum_visitor):
  ParallelEventVisitor(), waveform_sum_visitor_(waveform_sum_visitor), adopt_waveform_sum_visitor_(adopt_waveform_sum_visitor),
  capture_datum_(capture_datum), max_event_number_(max_event_number),
  min_event_number_(min_event_number)
{
  //nothing to see here
}

SimpleChargeCaptureParallelEventVisitor::~SimpleChargeCaptureParallelEventVisitor()
{
  for(auto& icaptured : captured_data_) {
    delete icaptured.second;
  }
  if(adopt_waveform_sum_visitor_) {
    delete waveform_sum_visitor_;
  }
}

SimpleChargeCaptureParallelEventVisitor* SimpleChargeCaptureParallelEventVisitor::new_sub_visitor(
  std::map<calin::iact_data::event_visitor::ParallelEventVisitor*,
    calin::iact_data::event_visitor::ParallelEventVisitor*> antecedent_visitors)
{
  auto* v = dynamic_cast<OptimalWindowSumWaveformTreatmentParallelEventVisitor*>(
    antecedent_visitors[waveform_sum_visitor_]);
  auto* child = new SimpleChargeCaptureParallelEventVisitor(v,
    capture_datum_, max_event_number_, min_event_number_, false);
  child->parent_ = this;
  return child;
}

bool SimpleChargeCaptureParallelEventVisitor::visit_telescope_run(
  const calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration* run_config,
  calin::iact_data::event_visitor::EventLifetimeManager* event_lifetime_manager,
  calin::ix::provenance::chronicle::ProcessingRecord* processing_record)
{
  if(processing_record) {
    processing_record->set_type("SimpleChargeCaptureParallelEventVisitor");
    processing_record->set_description("Waveform charge capture");
    auto* config_json = processing_record->add_config();
    config_json = processing_record->add_config();
    std::vector<std::pair<std::string,std::string> > keyval;
    keyval.emplace_back("waveformSumVisitor",
      calin::io::json::json_string_value(calin::util::string::instance_identifier(waveform_sum_visitor_)));
    switch(capture_datum_) {
#define CASE_DATUM(x) case x : keyval.emplace_back("captureDatum", calin::io::json::json_string_value(#x)); break
    CASE_DATUM(CCD_MAX_SAMPLE_VALUE);
    CASE_DATUM(CCD_MAX_SAMPLE_INDEX);
    CASE_DATUM(CCD_BKG_WINDOW_SUM);
    CASE_DATUM(CCD_SIG_WINDOW_SUM);
    CASE_DATUM(CCD_OPT_WINDOW_SUM);
    CASE_DATUM(CCD_OPT_WINDOW_SUM_QT);
    CASE_DATUM(CCD_OPT_WINDOW_INDEX);
    CASE_DATUM(CCD_WAVEFORM_SUM);
    }
    keyval.emplace_back("maxEventNumber",
      calin::io::json::json_value(max_event_number_));
    keyval.emplace_back("minEventNumber",
      calin::io::json::json_value(min_event_number_));
    config_json->set_json(calin::io::json::json_for_dictionary(keyval));
  }
  run_config_.Clear();
  run_config_.CopyFrom(*run_config);
  for(auto& icaptured : captured_data_) {
    delete icaptured.second;
  }
  captured_data_.clear();
  return true;
}

bool SimpleChargeCaptureParallelEventVisitor::leave_telescope_run(
  calin::ix::provenance::chronicle::ProcessingRecord* processing_record)
{
  return true;
}

bool SimpleChargeCaptureParallelEventVisitor::visit_telescope_event(uint64_t seq_index,
  calin::ix::iact_data::telescope_event::TelescopeEvent* event)
{
  if(event->local_event_number()<min_event_number_ or (max_event_number_>0 and
      event->local_event_number()>=max_event_number_)) {
    return true;
  }
  unsigned nchan = waveform_sum_visitor_->nchan();
  auto* data = new CapturedEventData;
  data->event_number = event->local_event_number();
  data->time = event->absolute_event_time().time_ns();
  data->trigger_code = event->trigger_type();
  data->values.resize(nchan);
  const int32_t*__restrict__ array = nullptr;
  switch(capture_datum_) {
  case CCD_MAX_SAMPLE_VALUE:
    array = waveform_sum_visitor_->array_chan_max(); break;
  case CCD_MAX_SAMPLE_INDEX:
    array = waveform_sum_visitor_->array_chan_max_index(); break;
  case CCD_BKG_WINDOW_SUM:
    array = waveform_sum_visitor_->array_chan_bkg_win_sum(); break;
  case CCD_SIG_WINDOW_SUM:
    array = waveform_sum_visitor_->array_chan_sig_win_sum(); break;
  case CCD_OPT_WINDOW_SUM:
    array = waveform_sum_visitor_->array_chan_opt_win_sum(); break;
  case CCD_OPT_WINDOW_SUM_QT:
    array = waveform_sum_visitor_->array_chan_opt_win_sum_qt(); break;
  case CCD_OPT_WINDOW_INDEX:
    array = waveform_sum_visitor_->array_chan_opt_win_index(); break;
  case CCD_WAVEFORM_SUM:
    array = waveform_sum_visitor_->array_chan_all_sum(); break;
  }
  for(unsigned ichan=0;ichan<nchan;ichan++) {
    data->values[ichan] =
      (waveform_sum_visitor_->array_chan_signal_type()[ichan]&0x3) | (array[ichan] << 2);
  }
  captured_data_[data->event_number] = data;
  return true;
}

bool SimpleChargeCaptureParallelEventVisitor::merge_results()
{
  if(parent_) {
    for(auto& icaptured : captured_data_) {
      if(parent_->captured_data_[icaptured.first]) {
        delete icaptured.second;
      } else {
        parent_->captured_data_[icaptured.first] = icaptured.second;
      }
    }
    captured_data_.clear();
  }
  return true;
}

std::vector<uint64_t> SimpleChargeCaptureParallelEventVisitor::keys() const
{
  std::vector<uint64_t> kv;
  kv.reserve(captured_data_.size());
  for(auto& icaptured : captured_data_) {
    kv.push_back(icaptured.first);
  }
  return kv;
}

bool SimpleChargeCaptureParallelEventVisitor::has(uint64_t event_number) const
{
  return captured_data_.count(event_number);
}

int64_t SimpleChargeCaptureParallelEventVisitor::time(uint64_t event_number) const
{
  auto icd = captured_data_.find(event_number);
  if(icd == captured_data_.end()) {
    throw std::range_error("Event not found : " + std::to_string(event_number));
  }
  return icd->second->time;
}

int32_t SimpleChargeCaptureParallelEventVisitor::trigger_code(uint64_t event_number) const
{
  auto icd = captured_data_.find(event_number);
  if(icd == captured_data_.end()) {
    throw std::range_error("Event not found : " + std::to_string(event_number));
  }
  return icd->second->trigger_code;
}

Eigen::VectorXi SimpleChargeCaptureParallelEventVisitor::signal_types(uint64_t event_number) const
{
  auto icd = captured_data_.find(event_number);
  if(icd == captured_data_.end()) {
    throw std::range_error("Event not found : " + std::to_string(event_number));
  }
  const auto* cd = icd->second;
  unsigned nchan = cd->values.size();
  Eigen::VectorXi values(nchan);
  for(unsigned ichan=0;ichan<nchan;++ichan) {
    values[ichan] = cd->values[ichan] & 0x03;
  }
  return values;
}

Eigen::VectorXi SimpleChargeCaptureParallelEventVisitor::values(uint64_t event_number) const
{
  auto icd = captured_data_.find(event_number);
  if(icd == captured_data_.end()) {
    throw std::range_error("Event not found : " + std::to_string(event_number));
  }
  const auto* cd = icd->second;
  unsigned nchan = cd->values.size();
  Eigen::VectorXi values(nchan);
  for(unsigned ichan=0;ichan<nchan;++ichan) {
    values[ichan] = cd->values[ichan] >> 2;
  }
  return values;
}
