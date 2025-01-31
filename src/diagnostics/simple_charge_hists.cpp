/*

   calin/diagnostics/simple_charge_hists.cpp -- Stephen Fegan -- 2020-04-18

   Various histograms computed from channel data

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

#include <util/log.hpp>
#include <math/special.hpp>
#include <util/algorithm.hpp>
#include <util/string.hpp>
#include <io/json.hpp>
#include <diagnostics/simple_charge_hists.hpp>
#include <math/covariance_calc.hpp>

using namespace calin::util::log;
using namespace calin::diagnostics::simple_charge_hists;
using namespace calin::ix::diagnostics::simple_charge_hists;
using namespace calin::iact_data::waveform_treatment_event_visitor;

using calin::math::special::SQR;
using calin::math::covariance_calc::cov_i64_gen;
using calin::math::covariance_calc::cov_double_gen;

SimpleChargeHistsParallelEventVisitor::
SimpleChargeHistsParallelEventVisitor(
    OptimalWindowSumWaveformTreatmentParallelEventVisitor* high_gain_visitor,
    OptimalWindowSumWaveformTreatmentParallelEventVisitor* low_gain_visitor,
    const SimpleChargeHistsConfig& config,
    SimpleChargeHistsFilter* filter, bool adopt_filter):
  ParallelEventVisitor(), config_(config),
  high_gain_visitor_(high_gain_visitor), low_gain_visitor_(low_gain_visitor),
  filter_(filter), adopt_filter_(adopt_filter)
{
  // nothing to see here
}

SimpleChargeHistsParallelEventVisitor::~SimpleChargeHistsParallelEventVisitor()
{
  for(auto* h : chan_hists_)delete h;
  delete cam_hists_high_gain_;
  delete cam_hists_low_gain_;
  delete cam_hists_dual_gain_;
  if(adopt_filter_)delete filter_;
}

SimpleChargeHistsParallelEventVisitor* SimpleChargeHistsParallelEventVisitor::new_sub_visitor(
  std::map<calin::iact_data::event_visitor::ParallelEventVisitor*,
    calin::iact_data::event_visitor::ParallelEventVisitor*> antecedent_visitors)
{
  auto* hgv = dynamic_cast<OptimalWindowSumWaveformTreatmentParallelEventVisitor*>(
    antecedent_visitors[high_gain_visitor_]);
  auto* lgv = dynamic_cast<OptimalWindowSumWaveformTreatmentParallelEventVisitor*>(
    antecedent_visitors[low_gain_visitor_]); // good in case of nullptr also
  auto* child = new SimpleChargeHistsParallelEventVisitor(hgv, lgv, config_,
    filter_==nullptr?nullptr:filter_->clone(), filter_==nullptr?false:true);
  child->parent_ = this;
  return child;
}

bool SimpleChargeHistsParallelEventVisitor::visit_telescope_run(
  const calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration* run_config,
  calin::iact_data::event_visitor::EventLifetimeManager* event_lifetime_manager,
  calin::ix::provenance::chronicle::ProcessingRecord* processing_record)
{
  if(processing_record) {
    processing_record->set_type("SimpleChargeHistsParallelEventVisitor");
    processing_record->set_description("Per-channel waveform sum histograms");
    auto* config_json = processing_record->add_config();
    config_json->set_type(config_.GetTypeName());
    config_json->set_json(calin::io::json::encode_protobuf_to_json_string(config_));
    config_json = processing_record->add_config();
    std::vector<std::pair<std::string,std::string> > keyval;
    keyval.emplace_back("highGainWaveformSumInstance",
      calin::io::json::json_string_value(calin::util::string::instance_identifier(high_gain_visitor_)));
    keyval.emplace_back("lowGainWaveformSumInstance",
      calin::io::json::json_string_value(calin::util::string::instance_identifier(low_gain_visitor_)));
    keyval.emplace_back("filterInstance",
      calin::io::json::json_string_value(calin::util::string::instance_identifier(filter_)));
    config_json->set_json(calin::io::json::json_for_dictionary(keyval));
  }

  has_dual_gain_ = (run_config->camera_layout().adc_gains() !=
    calin::ix::iact_data::instrument_layout::CameraLayout::SINGLE_GAIN);

  for(auto* h :chan_hists_)delete h;
  chan_hists_.resize(run_config->configured_channel_id_size());
  for(auto*& h : chan_hists_) {
    h = new ChannelHists(has_dual_gain_, config_);
  }
  delete cam_hists_high_gain_;
  delete cam_hists_low_gain_;
  delete cam_hists_dual_gain_;
  if(config_.has_high_gain() and config_.high_gain().enable_hists()) {
    cam_hists_high_gain_ = new SingleGainCameraHists(config_.high_gain());
  } else {
    cam_hists_high_gain_ = nullptr;
  }
  if(has_dual_gain_ and config_.has_low_gain() and config_.low_gain().enable_hists()) {
    cam_hists_low_gain_ = new SingleGainCameraHists(config_.low_gain());
  } else {
    cam_hists_low_gain_ = nullptr;
  }
  if(has_dual_gain_ and config_.has_dual_gain() and config_.dual_gain().enable_hists()) {
    cam_hists_dual_gain_ = new DualGainCameraHists(config_.dual_gain());
  } else {
    cam_hists_dual_gain_ = nullptr;
  }
  return true;
}

bool SimpleChargeHistsParallelEventVisitor::leave_telescope_run(
  calin::ix::provenance::chronicle::ProcessingRecord* processing_record)
{
  return true;
}

void SimpleChargeHistsParallelEventVisitor::record_one_visitor_data(
  uint64_t seq_index, const calin::ix::iact_data::telescope_event::TelescopeEvent* event,
  const calin::iact_data::waveform_treatment_event_visitor::OptimalWindowSumWaveformTreatmentParallelEventVisitor* sum_visitor)
{
  if(not sum_visitor or not sum_visitor->is_same_event(seq_index))
    return;

  bool has_rel_qsum = false;
  for(unsigned ichan=0; ichan<sum_visitor->nchan(); ichan++) {
    SingleGainChannelHists* hists = nullptr;
    SingleGainCameraHists* cam_hists = nullptr;
    switch(sum_visitor->array_chan_signal_type()[ichan]) {
    case calin::ix::iact_data::telescope_event::SIGNAL_UNIQUE_GAIN:
    case calin::ix::iact_data::telescope_event::SIGNAL_HIGH_GAIN:
      if(filter_==nullptr or filter_->event_should_be_accepted(event, ichan, /*low_gain=*/ false)) {
        hists = chan_hists_[ichan]->high_gain;
        cam_hists = cam_hists_high_gain_;
      }
      break;
    case calin::ix::iact_data::telescope_event::SIGNAL_LOW_GAIN:
      if(filter_==nullptr or filter_->event_should_be_accepted(event, ichan, /*low_gain=*/ true)) {
        hists = chan_hists_[ichan]->low_gain;
        cam_hists = cam_hists_low_gain_;
      }
      break;
    case calin::ix::iact_data::telescope_event::SIGNAL_NONE:
    default:
      // do nothing
      break;
    }

    if(hists) {
      if(hists->full_wf_qsum)
        hists->full_wf_qsum->insert(sum_visitor->array_chan_all_sum()[ichan]);
      if(hists->full_wf_max)
        hists->full_wf_max->insert(sum_visitor->array_chan_max()[ichan]);
      if(hists->full_wf_max_index)
        hists->full_wf_max_index->insert(sum_visitor->array_chan_max_index()[ichan]);
      if(hists->opt_win_qsum)
        hists->opt_win_qsum->insert(sum_visitor->array_chan_opt_win_sum()[ichan]);
      if(hists->opt_win_qtsum)
        hists->opt_win_qtsum->insert(sum_visitor->array_chan_opt_win_sum_qt()[ichan]);
      if(hists->opt_win_index)
        hists->opt_win_index->insert(sum_visitor->array_chan_opt_win_index()[ichan]);
      if(hists->opt_ped_qsum_diff)
        hists->opt_ped_qsum_diff->insert(sum_visitor->array_chan_opt_win_sum()[ichan]
          - sum_visitor->array_chan_bkg_win_sum()[ichan]);
      if(hists->ped_win_qsum)
        hists->ped_win_qsum->insert(sum_visitor->array_chan_bkg_win_sum()[ichan]);
      if(hists->sig_win_qsum)
        hists->sig_win_qsum->insert(sum_visitor->array_chan_sig_win_sum()[ichan]);
      if(hists->sig_ped_qsum_diff)
        hists->sig_ped_qsum_diff->insert(sum_visitor->array_chan_sig_win_sum()[ichan]
          - sum_visitor->array_chan_bkg_win_sum()[ichan]);
      if(hists->rel_qsum) {
        has_rel_qsum = true;
      }
    }

    if(cam_hists) {
      ++cam_hists->event_nchan_present_sum;
      cam_hists->event_full_wf_qsum += sum_visitor->array_chan_all_sum()[ichan];
      cam_hists->event_opt_win_qsum += sum_visitor->array_chan_opt_win_sum()[ichan];
      cam_hists->event_opt_ped_win_qsum_diff += sum_visitor->array_chan_opt_win_sum()[ichan]
        - sum_visitor->array_chan_bkg_win_sum()[ichan];
      cam_hists->event_ped_win_qsum += sum_visitor->array_chan_bkg_win_sum()[ichan];
      cam_hists->event_sig_win_qsum += sum_visitor->array_chan_sig_win_sum()[ichan];
      cam_hists->event_sig_ped_win_qsum_diff += sum_visitor->array_chan_sig_win_sum()[ichan]
        - sum_visitor->array_chan_bkg_win_sum()[ichan];
    }
  }

  if(!has_rel_qsum) return;

  for(unsigned ichan=0; ichan<sum_visitor->nchan(); ichan++) {
    SingleGainChannelHists* hists = nullptr;
    SingleGainCameraHists* cam_hists = nullptr;
    switch(sum_visitor->array_chan_signal_type()[ichan]) {
    case calin::ix::iact_data::telescope_event::SIGNAL_UNIQUE_GAIN:
    case calin::ix::iact_data::telescope_event::SIGNAL_HIGH_GAIN:
      if(filter_==nullptr or filter_->event_should_be_accepted(event, ichan, /*low_gain=*/ false)) {
        hists = chan_hists_[ichan]->high_gain;
        cam_hists = cam_hists_high_gain_;
      }
      break;
    case calin::ix::iact_data::telescope_event::SIGNAL_LOW_GAIN:
      if(filter_==nullptr or filter_->event_should_be_accepted(event, ichan, /*low_gain=*/ true)) {
        hists = chan_hists_[ichan]->low_gain;
        cam_hists = cam_hists_low_gain_;
      }
    case calin::ix::iact_data::telescope_event::SIGNAL_NONE:
    default:
      // do nothing
      break;
    }

    if(hists and cam_hists and hists->rel_qsum) {
      double chan_val;
      double cam_val;
      switch(config_.rel_qsum_source()) {
      case OPT_PED_DIFF:
      default:
        chan_val = sum_visitor->array_chan_opt_win_sum()[ichan] - sum_visitor->array_chan_bkg_win_sum()[ichan];
        cam_val = cam_hists->event_opt_ped_win_qsum_diff;
        break;
      case SIG_PED_DIFF:
        chan_val = sum_visitor->array_chan_sig_win_sum()[ichan] - sum_visitor->array_chan_bkg_win_sum()[ichan];
        cam_val = cam_hists->event_sig_ped_win_qsum_diff;
        break;
      case OPT_WINDOW:
        chan_val = sum_visitor->array_chan_opt_win_sum()[ichan];
        cam_val = cam_hists->event_opt_win_qsum;
        break;
      case SIG_WINDOW:
        chan_val = sum_visitor->array_chan_sig_win_sum()[ichan];
        cam_val = cam_hists->event_sig_win_qsum;
        break;
      }
      hists->rel_qsum->insert(chan_val/cam_val*cam_hists->event_nchan_present_sum);
    }
  }
}

bool SimpleChargeHistsParallelEventVisitor::visit_telescope_event(uint64_t seq_index,
  calin::ix::iact_data::telescope_event::TelescopeEvent* event)
{
  if(cam_hists_high_gain_) {
    cam_hists_high_gain_->reset_sums();
  }
  if(cam_hists_low_gain_) {
    cam_hists_low_gain_->reset_sums();
  }
  if(high_gain_visitor_) {
    record_one_visitor_data(seq_index, event, high_gain_visitor_);
    if(cam_hists_high_gain_) {
      if(cam_hists_high_gain_->nchan_present) {
        cam_hists_high_gain_->nchan_present->insert(cam_hists_high_gain_->event_nchan_present_sum);
      }
      if(cam_hists_high_gain_->event_nchan_present_sum) {
        double scale = double(high_gain_visitor_->nchan())/double(cam_hists_high_gain_->event_nchan_present_sum);
        if(cam_hists_high_gain_->full_wf_qsum) {
          cam_hists_high_gain_->full_wf_qsum->insert(scale*cam_hists_high_gain_->event_full_wf_qsum);
        }
        if(cam_hists_high_gain_->opt_win_qsum) {
          cam_hists_high_gain_->opt_win_qsum->insert(scale*cam_hists_high_gain_->event_opt_win_qsum);
        }
        if(cam_hists_high_gain_->opt_ped_qsum_diff) {
          cam_hists_high_gain_->opt_ped_qsum_diff->insert(scale*cam_hists_high_gain_->event_opt_ped_win_qsum_diff);
        }
        if(cam_hists_high_gain_->ped_win_qsum) {
          cam_hists_high_gain_->ped_win_qsum->insert(scale*cam_hists_high_gain_->event_ped_win_qsum);
        }
        if(cam_hists_high_gain_->sig_win_qsum) {
          cam_hists_high_gain_->sig_win_qsum->insert(scale*cam_hists_high_gain_->event_sig_win_qsum);
        }
        if(cam_hists_high_gain_->sig_ped_qsum_diff) {
          cam_hists_high_gain_->sig_ped_qsum_diff->insert(scale*cam_hists_high_gain_->event_sig_ped_win_qsum_diff);
        }
      }
    }
  }
  if(low_gain_visitor_) {
    record_one_visitor_data(seq_index, event, low_gain_visitor_);
    if(cam_hists_low_gain_) {
      if(cam_hists_low_gain_->nchan_present) {
        cam_hists_low_gain_->nchan_present->insert(cam_hists_low_gain_->event_nchan_present_sum);
      }
      if(cam_hists_low_gain_->event_nchan_present_sum) {
        double scale = double(low_gain_visitor_->nchan())/double(cam_hists_low_gain_->event_nchan_present_sum);
        if(cam_hists_low_gain_->full_wf_qsum) {
          cam_hists_low_gain_->full_wf_qsum->insert(scale*cam_hists_low_gain_->event_full_wf_qsum);
        }
        if(cam_hists_low_gain_->opt_win_qsum) {
          cam_hists_low_gain_->opt_win_qsum->insert(scale*cam_hists_low_gain_->event_opt_win_qsum);
        }
        if(cam_hists_low_gain_->opt_ped_qsum_diff) {
          cam_hists_low_gain_->opt_ped_qsum_diff->insert(scale*cam_hists_low_gain_->event_opt_ped_win_qsum_diff);
        }
        if(cam_hists_low_gain_->ped_win_qsum) {
          cam_hists_low_gain_->ped_win_qsum->insert(scale*cam_hists_low_gain_->event_ped_win_qsum);
        }
        if(cam_hists_low_gain_->sig_win_qsum) {
          cam_hists_low_gain_->sig_win_qsum->insert(scale*cam_hists_low_gain_->event_sig_win_qsum);
        }
        if(cam_hists_low_gain_->sig_ped_qsum_diff) {
          cam_hists_low_gain_->sig_ped_qsum_diff->insert(scale*cam_hists_low_gain_->event_sig_ped_win_qsum_diff);
        }
      }
    }
  }
  if((high_gain_visitor_ or low_gain_visitor_) and cam_hists_dual_gain_ and cam_hists_dual_gain_->nchan_present) {
    unsigned nchan_present = 0;
    if(cam_hists_high_gain_) {
      nchan_present += cam_hists_high_gain_->event_nchan_present_sum;
    }
    if(cam_hists_low_gain_) {
      nchan_present += cam_hists_low_gain_->event_nchan_present_sum;
    }
    cam_hists_dual_gain_->nchan_present->insert(nchan_present);
  }
  return true;
}

void SimpleChargeHistsParallelEventVisitor::merge_one_gain_hists(SingleGainChannelHists* into,
  const SingleGainChannelHists* from)
{
  if(from->full_wf_qsum)into->full_wf_qsum->insert_hist(*from->full_wf_qsum);
  if(from->full_wf_max)into->full_wf_max->insert_hist(*from->full_wf_max);
  if(from->full_wf_max_index)into->full_wf_max_index->insert_hist(*from->full_wf_max_index);
  if(from->opt_win_qsum)into->opt_win_qsum->insert_hist(*from->opt_win_qsum);
  if(from->opt_win_qtsum)into->opt_win_qtsum->insert_hist(*from->opt_win_qtsum);
  if(from->opt_win_index)into->opt_win_index->insert_hist(*from->opt_win_index);
  if(from->ped_win_qsum)into->ped_win_qsum->insert_hist(*from->ped_win_qsum);
  if(from->sig_win_qsum)into->sig_win_qsum->insert_hist(*from->sig_win_qsum);
  if(from->opt_ped_qsum_diff)into->opt_ped_qsum_diff->insert_hist(*from->opt_ped_qsum_diff);
  if(from->sig_ped_qsum_diff)into->sig_ped_qsum_diff->insert_hist(*from->sig_ped_qsum_diff);
  if(from->rel_qsum)into->rel_qsum->insert_hist(*from->rel_qsum);
}

void SimpleChargeHistsParallelEventVisitor::merge_one_gain_cam_hists(SingleGainCameraHists* into,
  const SingleGainCameraHists* from)
{
  if(from->nchan_present)into->nchan_present->insert_hist(*from->nchan_present);
  if(from->full_wf_qsum)into->full_wf_qsum->insert_hist(*from->full_wf_qsum);
  if(from->opt_win_qsum)into->opt_win_qsum->insert_hist(*from->opt_win_qsum);
  if(from->ped_win_qsum)into->ped_win_qsum->insert_hist(*from->ped_win_qsum);
  if(from->sig_win_qsum)into->sig_win_qsum->insert_hist(*from->sig_win_qsum);
  if(from->opt_ped_qsum_diff)into->opt_ped_qsum_diff->insert_hist(*from->opt_ped_qsum_diff);
  if(from->sig_ped_qsum_diff)into->sig_ped_qsum_diff->insert_hist(*from->sig_ped_qsum_diff);
}

void SimpleChargeHistsParallelEventVisitor::merge_dual_gain_cam_hists(DualGainCameraHists* into,
  const DualGainCameraHists* from)
{
  if(from->nchan_present)into->nchan_present->insert_hist(*from->nchan_present);
}

bool SimpleChargeHistsParallelEventVisitor::merge_results()
{
  if(parent_) {
    for(unsigned ichan=0; ichan<chan_hists_.size(); ichan++) {
      if(chan_hists_[ichan]->high_gain) {
        merge_one_gain_hists(parent_->chan_hists_[ichan]->high_gain, chan_hists_[ichan]->high_gain);
      }
      if(chan_hists_[ichan]->low_gain) {
        merge_one_gain_hists(parent_->chan_hists_[ichan]->low_gain, chan_hists_[ichan]->low_gain);
      }
    }
    if(cam_hists_high_gain_) {
      merge_one_gain_cam_hists(parent_->cam_hists_high_gain_, cam_hists_high_gain_);
    }
    if(cam_hists_low_gain_) {
      merge_one_gain_cam_hists(parent_->cam_hists_low_gain_, cam_hists_low_gain_);
    }
    if(cam_hists_dual_gain_) {
      merge_dual_gain_cam_hists(parent_->cam_hists_dual_gain_, cam_hists_dual_gain_);
    }
  }
  return true;
}

void SimpleChargeHistsParallelEventVisitor::extract_one_gain_hists(
  calin::ix::diagnostics::simple_charge_hists::OneGainSimpleChargeHists* into,
  const SingleGainChannelHists* from) const
{
  if(from->full_wf_qsum)from->full_wf_qsum->serialize(into->mutable_full_wf_qsum());
  if(from->full_wf_max)from->full_wf_max->serialize(into->mutable_full_wf_max());
  if(from->full_wf_max_index)from->full_wf_max_index->serialize(into->mutable_full_wf_max_index());
  if(from->opt_win_qsum)from->opt_win_qsum->serialize(into->mutable_opt_win_qsum());
  if(from->opt_win_qtsum)from->opt_win_qtsum->serialize(into->mutable_opt_win_qtsum());
  if(from->opt_win_index)from->opt_win_index->serialize(into->mutable_opt_win_index());
  if(from->ped_win_qsum)from->ped_win_qsum->serialize(into->mutable_ped_win_qsum());
  if(from->sig_win_qsum)from->sig_win_qsum->serialize(into->mutable_sig_win_qsum());
  if(from->opt_ped_qsum_diff)from->opt_ped_qsum_diff->serialize(into->mutable_opt_ped_qsum_diff());
  if(from->sig_ped_qsum_diff)from->sig_ped_qsum_diff->serialize(into->mutable_sig_ped_qsum_diff());
  if(from->rel_qsum) {
    from->rel_qsum->serialize(into->mutable_rel_qsum());
    into->mutable_rel_qsum()->set_name(into->rel_qsum().name() + " ("
      + RelativeQSumSource_Name(config_.rel_qsum_source()) + ")");
  }
}

void SimpleChargeHistsParallelEventVisitor::extract_one_gain_cam_hists(
  calin::ix::diagnostics::simple_charge_hists::OneGainSimpleChargeCameraHists* into,
  const SingleGainCameraHists* from) const
{
  if(from->nchan_present)from->nchan_present->serialize(into->mutable_nchan_present());
  if(from->full_wf_qsum)from->full_wf_qsum->serialize(into->mutable_full_wf_qsum());
  if(from->opt_win_qsum)from->opt_win_qsum->serialize(into->mutable_opt_win_qsum());
  if(from->ped_win_qsum)from->ped_win_qsum->serialize(into->mutable_ped_win_qsum());
  if(from->sig_win_qsum)from->sig_win_qsum->serialize(into->mutable_sig_win_qsum());
  if(from->opt_ped_qsum_diff)from->opt_ped_qsum_diff->serialize(into->mutable_opt_ped_qsum_diff());
  if(from->sig_ped_qsum_diff)from->sig_ped_qsum_diff->serialize(into->mutable_sig_ped_qsum_diff());
}

void SimpleChargeHistsParallelEventVisitor::extract_dual_gain_cam_hists(
  calin::ix::diagnostics::simple_charge_hists::DualGainSimpleChargeCameraHists* into,
  const DualGainCameraHists* from) const
{
  if(from->nchan_present)from->nchan_present->serialize(into->mutable_nchan_present());
}

calin::ix::diagnostics::simple_charge_hists::SimpleChargeHists*
SimpleChargeHistsParallelEventVisitor::simple_charge_hists(
  calin::ix::diagnostics::simple_charge_hists::SimpleChargeHists* hists) const
{
  if(hists == nullptr) {
    hists = new calin::ix::diagnostics::simple_charge_hists::SimpleChargeHists();
  } else {
    hists->Clear();
  }

  for(unsigned ichan=0; ichan<chan_hists_.size(); ichan++) {
    if(chan_hists_[ichan]->high_gain) {
      extract_one_gain_hists(hists->add_high_gain_channel(), chan_hists_[ichan]->high_gain);
    }
    if(chan_hists_[ichan]->low_gain) {
      extract_one_gain_hists(hists->add_low_gain_channel(), chan_hists_[ichan]->low_gain);
    }
  }

  if(cam_hists_high_gain_) {
    extract_one_gain_cam_hists(hists->mutable_high_gain_camera(), cam_hists_high_gain_);
  }
  if(cam_hists_low_gain_) {
    extract_one_gain_cam_hists(hists->mutable_low_gain_camera(), cam_hists_low_gain_);
  }
  if(cam_hists_dual_gain_) {
    extract_dual_gain_cam_hists(hists->mutable_dual_gain_camera(), cam_hists_dual_gain_);
  }
  return hists;
}

namespace {

  void pre_fill_hist_config(
    calin::ix::diagnostics::simple_charge_hists::SingleGainSimpleChargeHistsConfig* config,
    const std::string& trig, const std::string& gain,
    double qsum_dx, int qsum_nmax, int qsum_rebin,
    double max_dx, int max_nmax, int max_rebin)
  {
    config->Clear();
    config->set_enable_hists(true);

    auto* hist = config->mutable_full_wf_qsum();
    hist->set_enable(false);
    hist->set_dxval(qsum_dx);
    hist->set_xval_align(0.0);
    hist->set_name(trig + " events " + gain + "-gain full waveform charge sum");
    hist->set_xval_units("DC");
    hist->set_weight_units("events");
    hist->set_compactify_output(true);
    hist->set_max_dense_bins_in_output(qsum_nmax);
    hist->set_max_output_rebinning(qsum_rebin);

    hist = config->mutable_full_wf_max();
    hist->set_enable(false);
    hist->set_dxval(max_dx);
    hist->set_xval_align(0.0);
    hist->set_name(trig + " events " + gain + "-gain full waveform max sample");
    hist->set_xval_units("DC");
    hist->set_weight_units("events");
    hist->set_compactify_output(true);
    hist->set_max_dense_bins_in_output(max_nmax);
    hist->set_max_output_rebinning(max_rebin);

    hist = config->mutable_full_wf_max_index();
    hist->set_enable(false);
    hist->set_dxval(1.0);
    hist->set_xval_align(0.5);
    hist->set_name(trig + " events " + gain + "-gain full waveform index of max sample");
    hist->set_xval_units("samples");
    hist->set_weight_units("events");
    hist->set_compactify_output(false);
    hist->set_max_dense_bins_in_output(0);
    hist->set_max_output_rebinning(0);

    hist = config->mutable_opt_win_qsum();
    hist->set_enable(false);
    hist->set_dxval(qsum_dx);
    hist->set_xval_align(0.0);
    hist->set_name(trig + " events " + gain + "-gain optimal window charge sum");
    hist->set_xval_units("DC");
    hist->set_weight_units("events");
    hist->set_compactify_output(true);
    hist->set_max_dense_bins_in_output(qsum_nmax);
    hist->set_max_output_rebinning(qsum_rebin);

    hist = config->mutable_opt_win_qtsum();
    hist->set_enable(false);
    hist->set_dxval(qsum_dx);
    hist->set_xval_align(0.0);
    hist->set_name(trig + " events " + gain + "-gain optimal window charge x time sum");
    hist->set_xval_units("DC");
    hist->set_weight_units("events");
    hist->set_compactify_output(true);
    hist->set_max_dense_bins_in_output(qsum_nmax);
    hist->set_max_output_rebinning(qsum_rebin);

    hist = config->mutable_opt_win_index();
    hist->set_enable(false);
    hist->set_dxval(1.0);
    hist->set_xval_align(0.5);
    hist->set_name(trig + " events " + gain + "-gain optimal window index");
    hist->set_xval_units("samples");
    hist->set_weight_units("events");
    hist->set_compactify_output(false);
    hist->set_max_dense_bins_in_output(0);
    hist->set_max_output_rebinning(0);

    hist = config->mutable_ped_win_qsum();
    hist->set_enable(false);
    hist->set_dxval(qsum_dx);
    hist->set_xval_align(0.0);
    hist->set_name(trig + " events " + gain + "-gain fixed pedestal window charge sum");
    hist->set_xval_units("DC");
    hist->set_weight_units("events");
    hist->set_compactify_output(true);
    hist->set_max_dense_bins_in_output(qsum_nmax);
    hist->set_max_output_rebinning(qsum_rebin);

    hist = config->mutable_sig_win_qsum();
    hist->set_enable(false);
    hist->set_dxval(qsum_dx);
    hist->set_xval_align(0.0);
    hist->set_name(trig + " events " + gain + "-gain fixed pedestal window charge sum");
    hist->set_xval_units("DC");
    hist->set_weight_units("events");
    hist->set_compactify_output(true);
    hist->set_max_dense_bins_in_output(qsum_nmax);
    hist->set_max_output_rebinning(qsum_rebin);

    hist = config->mutable_opt_ped_qsum_diff();
    hist->set_enable(false);
    hist->set_dxval(qsum_dx);
    hist->set_xval_align(0.0);
    hist->set_name(trig + " events " + gain + "-gain difference between optimal and pedestal fixed-window charge sum");
    hist->set_xval_units("DC");
    hist->set_weight_units("events");
    hist->set_compactify_output(true);
    hist->set_max_dense_bins_in_output(qsum_nmax);
    hist->set_max_output_rebinning(qsum_rebin);

    hist = config->mutable_sig_ped_qsum_diff();
    hist->set_enable(false);
    hist->set_dxval(qsum_dx);
    hist->set_xval_align(0.0);
    hist->set_name(trig + " events " + gain + "-gain difference between signal and pedestal fixed-window charge sum");
    hist->set_xval_units("DC");
    hist->set_weight_units("events");
    hist->set_compactify_output(true);
    hist->set_max_dense_bins_in_output(qsum_nmax);
    hist->set_max_output_rebinning(qsum_rebin);

    hist = config->mutable_rel_qsum();
    hist->set_enable(false);
    hist->set_dxval(0.01);
    hist->set_xval_align(0.0);
    hist->set_name(trig + " events " + gain + "-gain charge sum relative to camera mean");
    hist->set_xval_units("1");
    hist->set_weight_units("events");
    hist->set_limited(true);
    hist->set_xval_limit_lo(0.0);
    hist->set_xval_limit_hi(10.0);
    hist->set_compactify_output(true);
    hist->set_max_dense_bins_in_output(200);
    hist->set_max_output_rebinning(5);

    hist = config->mutable_nchan_present();
    hist->set_enable(false);
    hist->set_dxval(1);
    hist->set_xval_align(0.5);
    hist->set_name(trig + " events " + gain + "-gain number of channels present");
    hist->set_xval_units("channels");
    hist->set_weight_units("events");
    hist->set_compactify_output(true);
    hist->set_max_dense_bins_in_output(0);
    hist->set_max_output_rebinning(0);
  }

  void pre_fill_dual_gain_hist_config(
    calin::ix::diagnostics::simple_charge_hists::DualGainSimpleChargeHistsConfig* config,
    const std::string& trig)
  {
    config->Clear();
    config->set_enable_hists(true);

    auto* hist = config->mutable_nchan_present();
    hist->set_enable(false);
    hist->set_dxval(1);
    hist->set_xval_align(0.5);
    hist->set_name(trig + " events dual-gain number of channels present");
    hist->set_xval_units("channels");
    hist->set_weight_units("events");
    hist->set_compactify_output(true);
    hist->set_max_dense_bins_in_output(0);
    hist->set_max_output_rebinning(0);
  }

} // anonymous namespace

calin::ix::diagnostics::simple_charge_hists::SimpleChargeHistsConfig
SimpleChargeHistsParallelEventVisitor::all_enabled_config(
  double qsum_dx, int qsum_nmax, int qsum_rebin,
  double max_dx, int max_nmax, int max_rebin)
{
  calin::ix::diagnostics::simple_charge_hists::SimpleChargeHistsConfig config;
  auto* hg_config = config.mutable_high_gain();
  pre_fill_hist_config(hg_config, "pedestal", "high",
    qsum_dx, qsum_nmax, qsum_rebin, max_dx, max_nmax, max_rebin);
  hg_config->mutable_full_wf_qsum()->set_enable(true);
  hg_config->mutable_full_wf_max()->set_enable(true);
  hg_config->mutable_full_wf_max_index()->set_enable(true);
  hg_config->mutable_opt_win_qsum()->set_enable(true);
  hg_config->mutable_opt_win_qtsum ()->set_enable(true);
  hg_config->mutable_opt_win_index()->set_enable(true);
  hg_config->mutable_ped_win_qsum()->set_enable(true);
  hg_config->mutable_sig_win_qsum()->set_enable(true);
  hg_config->mutable_sig_ped_qsum_diff()->set_enable(true);
  hg_config->mutable_nchan_present()->set_enable(true);

  auto* lg_config = config.mutable_low_gain();
  pre_fill_hist_config(lg_config, "pedestal", "low",
    qsum_dx, qsum_nmax, qsum_rebin, max_dx, max_nmax, max_rebin);
  lg_config->mutable_full_wf_qsum()->set_enable(true);
  lg_config->mutable_full_wf_max()->set_enable(true);
  lg_config->mutable_full_wf_max_index()->set_enable(true);
  lg_config->mutable_opt_win_qsum()->set_enable(true);
  lg_config->mutable_opt_win_qtsum ()->set_enable(true);
  lg_config->mutable_opt_win_index()->set_enable(true);
  lg_config->mutable_ped_win_qsum()->set_enable(true);
  lg_config->mutable_sig_win_qsum()->set_enable(true);
  lg_config->mutable_sig_ped_qsum_diff()->set_enable(true);
  lg_config->mutable_nchan_present()->set_enable(true);

  auto* g2_config = config.mutable_dual_gain();
  pre_fill_dual_gain_hist_config(g2_config, "pedestal");
  g2_config->mutable_nchan_present()->set_enable(true);

  return config;
}

calin::ix::diagnostics::simple_charge_hists::SimpleChargeHistsConfig
SimpleChargeHistsParallelEventVisitor::ped_trig_default_config()
{
  calin::ix::diagnostics::simple_charge_hists::SimpleChargeHistsConfig config;

  auto* hg_config = config.mutable_high_gain();
  pre_fill_hist_config(hg_config, "pedestal", "high", 1.0, 1000, 5, 1.0, 100, 5);
  hg_config->mutable_full_wf_qsum()->set_enable(true);
  hg_config->mutable_opt_win_qsum()->set_enable(true);
  hg_config->mutable_opt_win_index()->set_enable(true);
  hg_config->mutable_nchan_present()->set_enable(true);
  auto* lg_config = config.mutable_low_gain();
  pre_fill_hist_config(lg_config, "pedestal", "low", 1.0, 200, 5, 1.0, 100, 5);
  lg_config->mutable_full_wf_qsum()->set_enable(true);
  lg_config->mutable_opt_win_qsum()->set_enable(true);
  lg_config->mutable_opt_win_index()->set_enable(true);
  lg_config->mutable_nchan_present()->set_enable(true);
  auto* g2_config = config.mutable_dual_gain();
  pre_fill_dual_gain_hist_config(g2_config, "pedestal");
  g2_config->mutable_nchan_present()->set_enable(true);

  return config;
}

calin::ix::diagnostics::simple_charge_hists::SimpleChargeHistsConfig
SimpleChargeHistsParallelEventVisitor::phy_trig_default_config()
{
  calin::ix::diagnostics::simple_charge_hists::SimpleChargeHistsConfig config;

  auto* hg_config = config.mutable_high_gain();
  pre_fill_hist_config(hg_config, "physics", "high", 10.0, 100, 10, 1.0, 100, 5);
  hg_config->mutable_full_wf_max()->set_enable(true);
  hg_config->mutable_full_wf_max_index()->set_enable(true);
  hg_config->mutable_opt_win_qsum()->set_enable(true);
  hg_config->mutable_opt_win_index()->set_enable(true);
  hg_config->mutable_nchan_present()->set_enable(true);
  auto* lg_config = config.mutable_low_gain();
  pre_fill_hist_config(lg_config, "physics", "low", 1.0, 100, 10, 1.0, 100, 5);
  lg_config->mutable_full_wf_max()->set_enable(true);
  lg_config->mutable_full_wf_max_index()->set_enable(true);
  lg_config->mutable_opt_win_qsum()->set_enable(true);
  lg_config->mutable_opt_win_index()->set_enable(true);
  lg_config->mutable_nchan_present()->set_enable(true);
  auto* g2_config = config.mutable_dual_gain();
  pre_fill_dual_gain_hist_config(g2_config, "physics");
  g2_config->mutable_nchan_present()->set_enable(true);

  return config;
}

calin::ix::diagnostics::simple_charge_hists::SimpleChargeHistsConfig
SimpleChargeHistsParallelEventVisitor::ext_trig_default_config()
{
  calin::ix::diagnostics::simple_charge_hists::SimpleChargeHistsConfig config;

  auto* hg_config = config.mutable_high_gain();
  pre_fill_hist_config(hg_config, "external-flasher", "high", 10.0, 100, 10, 1.0, 100, 5);
  hg_config->mutable_full_wf_max()->set_enable(true);
  hg_config->mutable_full_wf_max_index()->set_enable(true);
  hg_config->mutable_opt_win_qsum()->set_enable(true);
  hg_config->mutable_opt_win_index()->set_enable(true);
  hg_config->mutable_nchan_present()->set_enable(true);
  hg_config->mutable_opt_ped_qsum_diff()->set_enable(true);
  hg_config->mutable_rel_qsum()->set_enable(true);
  auto* lg_config = config.mutable_low_gain();
  pre_fill_hist_config(lg_config, "external-flasher", "low", 1.0, 100, 10, 1.0, 100, 5);
  lg_config->mutable_full_wf_max()->set_enable(true);
  lg_config->mutable_full_wf_max_index()->set_enable(true);
  lg_config->mutable_opt_win_qsum()->set_enable(true);
  lg_config->mutable_opt_win_index()->set_enable(true);
  lg_config->mutable_nchan_present()->set_enable(true);
  lg_config->mutable_opt_ped_qsum_diff()->set_enable(true);
  lg_config->mutable_rel_qsum()->set_enable(true);
  auto* g2_config = config.mutable_dual_gain();
  pre_fill_dual_gain_hist_config(g2_config, "external-flasher");
  g2_config->mutable_nchan_present()->set_enable(true);

  return config;
}

calin::ix::diagnostics::simple_charge_hists::SimpleChargeHistsConfig
SimpleChargeHistsParallelEventVisitor::int_trig_default_config()
{
  calin::ix::diagnostics::simple_charge_hists::SimpleChargeHistsConfig config;

  auto* hg_config = config.mutable_high_gain();
  pre_fill_hist_config(hg_config, "internal-flasher", "high", 1.0, 1000, 5, 1.0, 100, 5);
  hg_config->mutable_sig_win_qsum()->set_enable(true);
  hg_config->mutable_ped_win_qsum()->set_enable(true);
  hg_config->mutable_opt_win_qsum()->set_enable(true);
  hg_config->mutable_opt_win_index()->set_enable(true);
  hg_config->mutable_sig_ped_qsum_diff()->set_enable(true);
  hg_config->mutable_nchan_present()->set_enable(true);
  auto* lg_config = config.mutable_low_gain();
  pre_fill_hist_config(lg_config, "internal-flasher", "low", 1.0, 100, 10, 1.0, 100, 5);
  lg_config->mutable_nchan_present()->set_enable(true);
  auto* g2_config = config.mutable_dual_gain();
  pre_fill_dual_gain_hist_config(g2_config, "internal-flasher");
  g2_config->mutable_nchan_present()->set_enable(true);

  return config;
}

calin::ix::diagnostics::simple_charge_hists::SimpleChargeHistsConfig
SimpleChargeHistsParallelEventVisitor::l0_trig_bits_default_config()
{
  calin::ix::diagnostics::simple_charge_hists::SimpleChargeHistsConfig config;
  auto* hg_config = config.mutable_high_gain();

  hg_config->Clear();
  hg_config->set_enable_hists(true);

  auto* hist = hg_config->mutable_full_wf_max();
  hist->set_enable(true);
  hist->set_dxval(1.0);
  hist->set_xval_align(0.0);
  hist->set_name("All events (L0 trigger-bit selected) high-gain max sample");
  hist->set_xval_units("DC");
  hist->set_weight_units("events");
  hist->set_compactify_output(true);

  hist = hg_config->mutable_opt_win_qsum();
  hist->set_enable(true);
  hist->set_dxval(10.0);
  hist->set_xval_align(0.0);
  hist->set_name("All events (L0 trigger-bit selected) high-gain optimal window charge sum");
  hist->set_xval_units("DC");
  hist->set_weight_units("events");
  hist->set_compactify_output(true);

  return config;
}

SimpleChargeHistsFilter::~SimpleChargeHistsFilter()
{
  // nothing to see here
}

SimpleChargeHistsTriggerBitFilter::
SimpleChargeHistsTriggerBitFilter(bool trigger_bit_status_required_for_accept):
  SimpleChargeHistsFilter(),
  trigger_bit_status_required_for_accept_(trigger_bit_status_required_for_accept)
{
  // nothing to see here
}

SimpleChargeHistsTriggerBitFilter::~SimpleChargeHistsTriggerBitFilter()
{
  // nothing to see here
}

bool SimpleChargeHistsTriggerBitFilter::
event_should_be_accepted(const calin::ix::iact_data::telescope_event::TelescopeEvent* event,
  unsigned ichan, bool low_gain)
{
  if(event->has_trigger_map()) {
    return (event->trigger_map().trigger_image(ichan)>0) ==
      trigger_bit_status_required_for_accept_;
  } else {
    return false;
  }
}

SimpleChargeHistsTriggerBitFilter* SimpleChargeHistsTriggerBitFilter::clone()
{
  return new SimpleChargeHistsTriggerBitFilter(trigger_bit_status_required_for_accept_);
}
