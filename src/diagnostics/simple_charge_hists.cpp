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
    const SimpleChargeHistsConfig& config):
  ParallelEventVisitor(), config_(config),
  high_gain_visitor_(high_gain_visitor), low_gain_visitor_(low_gain_visitor)
{
  // nothing to see here
}

SimpleChargeHistsParallelEventVisitor::~SimpleChargeHistsParallelEventVisitor()
{
  for(auto* h : chan_hists_)delete h;
}

SimpleChargeHistsParallelEventVisitor* SimpleChargeHistsParallelEventVisitor::new_sub_visitor(
  std::map<calin::iact_data::event_visitor::ParallelEventVisitor*,
    calin::iact_data::event_visitor::ParallelEventVisitor*> antecedent_visitors)
{
  auto* hgv = dynamic_cast<OptimalWindowSumWaveformTreatmentParallelEventVisitor*>(
    antecedent_visitors[high_gain_visitor_]);
  auto* lgv = dynamic_cast<OptimalWindowSumWaveformTreatmentParallelEventVisitor*>(
    antecedent_visitors[low_gain_visitor_]); // good in case of nullptr also
  auto* child = new SimpleChargeHistsParallelEventVisitor(hgv, lgv, config_);
  child->parent_ = this;
  return child;
}

bool SimpleChargeHistsParallelEventVisitor::visit_telescope_run(
  const calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration* run_config,
  calin::iact_data::event_visitor::EventLifetimeManager* event_lifetime_manager)
{
  has_dual_gain_ = (run_config->camera_layout().adc_gains() !=
    calin::ix::iact_data::instrument_layout::CameraLayout::SINGLE_GAIN);

  for(auto* h :chan_hists_)delete h;
  chan_hists_.resize(run_config->configured_channel_id_size());
  for(auto*& h : chan_hists_) {
    h = new ChannelHists(has_dual_gain_, config_);
  }

  return true;
}

bool SimpleChargeHistsParallelEventVisitor::leave_telescope_run()
{
  return true;
}

void SimpleChargeHistsParallelEventVisitor::record_one_visitor_data(
  uint64_t seq_index, const calin::ix::iact_data::telescope_event::TelescopeEvent* event,
  const calin::iact_data::waveform_treatment_event_visitor::OptimalWindowSumWaveformTreatmentParallelEventVisitor* sum_visitor)
{
  if(not sum_visitor or not sum_visitor->is_same_event(seq_index))
    return;

  for(unsigned ichan=0; ichan<sum_visitor->nchan(); ichan++) {
    SingleGainChannelHists* hists = nullptr;
    switch(sum_visitor->array_chan_signal_type()[ichan]) {
    case calin::ix::iact_data::telescope_event::SIGNAL_UNIQUE_GAIN:
    case calin::ix::iact_data::telescope_event::SIGNAL_HIGH_GAIN:
      hists = chan_hists_[ichan]->high_gain;
      break;
    case calin::ix::iact_data::telescope_event::SIGNAL_LOW_GAIN:
      hists = chan_hists_[ichan]->low_gain;
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
      if(hists->ped_win_qsum)
        hists->ped_win_qsum->insert(sum_visitor->array_chan_bkg_win_sum()[ichan]);
      if(hists->sig_win_qsum)
        hists->sig_win_qsum->insert(sum_visitor->array_chan_sig_win_sum()[ichan]);
    }
  }
}

bool SimpleChargeHistsParallelEventVisitor::visit_telescope_event(uint64_t seq_index,
  calin::ix::iact_data::telescope_event::TelescopeEvent* event)
{
  if(high_gain_visitor_) {
    record_one_visitor_data(seq_index, event, high_gain_visitor_);
  }
  if(low_gain_visitor_) {
    record_one_visitor_data(seq_index, event, low_gain_visitor_);
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

  return config;
}

calin::ix::diagnostics::simple_charge_hists::SimpleChargeHistsConfig
SimpleChargeHistsParallelEventVisitor::ped_trig_default_config()
{
  calin::ix::diagnostics::simple_charge_hists::SimpleChargeHistsConfig config;

  auto* hg_config = config.mutable_high_gain();
  pre_fill_hist_config(hg_config, "pedestal", "high", 1.0, 1000, 5, 1.0, 100, 5);
  hg_config->mutable_full_wf_qsum()->set_enable(true);
  auto* lg_config = config.mutable_low_gain();
  pre_fill_hist_config(lg_config, "pedestal", "low", 1.0, 200, 5, 1.0, 100, 5);
  lg_config->mutable_full_wf_qsum()->set_enable(true);

  return config;
}

calin::ix::diagnostics::simple_charge_hists::SimpleChargeHistsConfig
SimpleChargeHistsParallelEventVisitor::phy_trig_default_config()
{
  calin::ix::diagnostics::simple_charge_hists::SimpleChargeHistsConfig config;

  auto* hg_config = config.mutable_high_gain();
  pre_fill_hist_config(hg_config, "external-flasher", "high", 10.0, 100, 10, 1.0, 100, 5);
  hg_config->mutable_full_wf_max()->set_enable(true);
  hg_config->mutable_opt_win_qsum()->set_enable(true);
  hg_config->mutable_opt_win_index()->set_enable(true);
  auto* lg_config = config.mutable_low_gain();
  pre_fill_hist_config(lg_config, "external-flasher", "low", 1.0, 100, 10, 1.0, 100, 5);
  lg_config->mutable_full_wf_max()->set_enable(true);
  lg_config->mutable_opt_win_qsum()->set_enable(true);
  lg_config->mutable_opt_win_index()->set_enable(true);

  return config;
}

calin::ix::diagnostics::simple_charge_hists::SimpleChargeHistsConfig
SimpleChargeHistsParallelEventVisitor::ext_trig_default_config()
{
  calin::ix::diagnostics::simple_charge_hists::SimpleChargeHistsConfig config;
  return config;
}

calin::ix::diagnostics::simple_charge_hists::SimpleChargeHistsConfig
SimpleChargeHistsParallelEventVisitor::int_trig_default_config()
{
  calin::ix::diagnostics::simple_charge_hists::SimpleChargeHistsConfig config;
  return config;
}
