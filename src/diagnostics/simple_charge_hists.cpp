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
  // for(auto* h : chan_hists_)delete h;
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
  // partials_.Clear();
  // results_.Clear();
  //
  // has_dual_gain_ = (run_config->camera_layout().adc_gains() !=
  //   calin::ix::iact_data::instrument_layout::CameraLayout::SINGLE_GAIN);
  // for(unsigned ichan=0;ichan<run_config->configured_channel_id_size();++ichan) {
  //   partials_.add_channel();
  // }
  //
  // for(auto* h :chan_hists_)delete h;
  // chan_hists_.resize(run_config->configured_channel_id_size());
  // for(auto*& h : chan_hists_) {
  //   h = new ChannelHists(has_dual_gain_, config_.ped_time_hist_resolution(),
  //     config_.dark_hist_hg_resolution(), config_.bright_hist_hg_resolution(),
  //     config_.dark_hist_lg_resolution(), config_.bright_hist_lg_resolution());
  // }

  // if(high_gain_visitor_) {
  //   high_gain_results_->Clear();
  //   high_gain_results_->set_integration_n(high_gain_visitor_->window_n());
  //   high_gain_results_->set_bkg_integration_0(high_gain_visitor_->bkg_window_0());
  //   int sig_win_0 = -2;
  //   for(int sw0 : high_gain_visitor_->sig_window_0()) {
  //     if(sig_win_0==-2)sig_win_0 = sw0;
  //     else if(sig_win_0 != sw0)sig_win_0 = -1;
  //     high_gain_results_->add_chan_sig_integration_0(sw0);
  //   };
  //   high_gain_results_->sig_integration_0(std::max(sig_win_0,-1));
  // }
  // if(low_gain_visitor_) {
  //   low_gain_results_->Clear();
  //   low_gain_results_->set_integration_n(low_gain_visitor_->window_n());
  //   low_gain_results_->set_bkg_integration_0(low_gain_visitor_->bkg_window_0());
  //   int sig_win_0 = -2;
  //   for(int sw0 : low_gain_visitor_->sig_window_0()) {
  //     if(sig_win_0==-2)sig_win_0 = sw0;
  //     else if(sig_win_0 != sw0)sig_win_0 = -1;
  //     low_gain_results_->add_chan_sig_integration_0(sw0);
  //   };
  //   low_gain_results_->sig_integration_0(std::max(sig_win_0,-1));
  // }
  return true;
}

bool SimpleChargeHistsParallelEventVisitor::leave_telescope_run()
{
  return true;
}

#if 0
void SimpleChargeHistsParallelEventVisitor::record_one_gain_channel_data(
  const calin::ix::iact_data::telescope_event::TelescopeEvent* event,
  const calin::iact_data::waveform_treatment_event_visitor::OptimalWindowSumWaveformTreatmentParallelEventVisitor* sum_visitor,
  unsigned ichan, double elapsed_event_time,
  calin::ix::diagnostics::simple_charge_stats::PartialOneGainChannelSimpleChargeStats* one_gain_stats,
  SingleGainChannelHists* one_gain_hists)
{
  one_gain_stats->increment_all_trig_num_events();
  one_gain_stats->increment_all_trig_ped_win_sum(sum_visitor->array_chan_bkg_win_sum()[ichan]);
  one_gain_stats->increment_all_trig_ped_win_sumsq(SQR(sum_visitor->array_chan_bkg_win_sum()[ichan]));
  if(event->trigger_type() == calin::ix::iact_data::telescope_event::TRIGGER_PEDESTAL) {
    double wf_all_sum = sum_visitor->array_chan_all_sum()[ichan];
    double sqr_wf_all_sum = SQR(wf_all_sum);
    one_gain_stats->increment_ped_trig_num_events();
    one_gain_stats->increment_ped_trig_full_wf_sum(wf_all_sum);
    one_gain_stats->increment_ped_trig_full_wf_sumsq(sqr_wf_all_sum);
    one_gain_hists->ped_wf_q_sum->insert(wf_all_sum);
    one_gain_hists->ped_wf_1_sum_vs_time->insert(elapsed_event_time);
    one_gain_hists->ped_wf_q_sum_vs_time->insert(elapsed_event_time, wf_all_sum);
    one_gain_hists->ped_wf_q2_sum_vs_time->insert(elapsed_event_time, sqr_wf_all_sum);
  } else if(event->trigger_type() == calin::ix::iact_data::telescope_event::TRIGGER_EXTERNAL_FLASHER) {
    one_gain_hists->ext_opt_win_q_sum->insert(sum_visitor->array_chan_opt_win_sum()[ichan]);
    one_gain_hists->ext_opt_win_index->insert(sum_visitor->array_chan_opt_win_index()[ichan]);
    one_gain_hists->ext_opt_win_max->insert(sum_visitor->array_chan_max()[ichan]);
  } else if(event->trigger_type() == calin::ix::iact_data::telescope_event::TRIGGER_PHYSICS) {
    one_gain_hists->phys_opt_win_q_sum->insert(sum_visitor->array_chan_opt_win_sum()[ichan]);
    one_gain_hists->phys_opt_win_index->insert(sum_visitor->array_chan_opt_win_index()[ichan]);
  }
}

void SimpleChargeHistsParallelEventVisitor::record_one_visitor_data(
  uint64_t seq_index, const calin::ix::iact_data::telescope_event::TelescopeEvent* event,
  const calin::iact_data::waveform_treatment_event_visitor::OptimalWindowSumWaveformTreatmentParallelEventVisitor* sum_visitor,
  calin::ix::diagnostics::simple_charge_stats::PartialSimpleChargeStats* partials)
{
  double elapsed_event_time = event->elapsed_event_time().time_ns() * 1e-9;

  if(sum_visitor and sum_visitor->is_same_event(seq_index)) {
    for(unsigned ichan=0; ichan<sum_visitor->nchan(); ichan++) {
      auto* pc = partials->mutable_channel(ichan);
      switch(sum_visitor->array_chan_signal_type()[ichan]) {
      case calin::ix::iact_data::telescope_event::SIGNAL_UNIQUE_GAIN:
      case calin::ix::iact_data::telescope_event::SIGNAL_HIGH_GAIN:
        record_one_gain_channel_data(event, sum_visitor, ichan, elapsed_event_time,
          pc->mutable_high_gain(), chan_hists_[ichan]->high_gain);
        break;
      case calin::ix::iact_data::telescope_event::SIGNAL_LOW_GAIN:
        record_one_gain_channel_data(event, sum_visitor, ichan, elapsed_event_time,
          pc->mutable_low_gain(), chan_hists_[ichan]->low_gain);
        break;
      case calin::ix::iact_data::telescope_event::SIGNAL_NONE:
      default:
        // do nothing
        break;
      }
    }
  }
}
#endif

bool SimpleChargeHistsParallelEventVisitor::visit_telescope_event(uint64_t seq_index,
  calin::ix::iact_data::telescope_event::TelescopeEvent* event)
{
  if(high_gain_visitor_) {
    // record_one_visitor_data(seq_index, event, high_gain_visitor_, &partials_);
  }
  if(low_gain_visitor_) {
    // record_one_visitor_data(seq_index, event, low_gain_visitor_, &partials_);
  }
  return true;
}

bool SimpleChargeHistsParallelEventVisitor::merge_results()
{
  if(parent_) {
    // parent_->partials_.IntegrateFrom(partials_);
  }
  return true;
}

#if 0
calin::ix::diagnostics::simple_charge_stats::SimpleChargeStats*
SimpleChargeStatsParallelEventVisitor::simple_charge_stats(
  calin::ix::diagnostics::simple_charge_stats::SimpleChargeStats* stats) const
{
  if(stats == nullptr)stats = results_.New();
  stats->CopyFrom(results_);
  return stats;
}
#endif

calin::ix::diagnostics::simple_charge_hists::SimpleChargeHists*
SimpleChargeHistsParallelEventVisitor::simple_charge_hists(
  calin::ix::diagnostics::simple_charge_hists::SimpleChargeHists* stats) const
{
  return nullptr;
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
