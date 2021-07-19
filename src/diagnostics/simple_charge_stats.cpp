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

#include <util/log.hpp>
#include <math/special.hpp>
#include <util/algorithm.hpp>
#include <diagnostics/simple_charge_stats.hpp>
#include <math/covariance_calc.hpp>

using namespace calin::util::log;
using namespace calin::diagnostics::simple_charge_stats;
using namespace calin::ix::diagnostics::simple_charge_stats;
using namespace calin::iact_data::waveform_treatment_event_visitor;

using calin::math::special::SQR;
using calin::math::covariance_calc::cov_i64_gen;
using calin::math::covariance_calc::cov_double_gen;

SimpleChargeStatsParallelEventVisitor::
SimpleChargeStatsParallelEventVisitor(
    OptimalWindowSumWaveformTreatmentParallelEventVisitor* high_gain_visitor,
    OptimalWindowSumWaveformTreatmentParallelEventVisitor* low_gain_visitor,
    const SimpleChargeStatsConfig& config):
  ParallelEventVisitor(), config_(config),
  high_gain_visitor_(high_gain_visitor), low_gain_visitor_(low_gain_visitor)
{
  // nothing to see here
}

SimpleChargeStatsParallelEventVisitor::~SimpleChargeStatsParallelEventVisitor()
{
  for(auto* h : chan_hists_)delete h;
  delete camera_hists_;
}

SimpleChargeStatsParallelEventVisitor* SimpleChargeStatsParallelEventVisitor::new_sub_visitor(
  std::map<calin::iact_data::event_visitor::ParallelEventVisitor*,
    calin::iact_data::event_visitor::ParallelEventVisitor*> antecedent_visitors)
{
  auto* hgv = dynamic_cast<OptimalWindowSumWaveformTreatmentParallelEventVisitor*>(
    antecedent_visitors[high_gain_visitor_]);
  auto* lgv = dynamic_cast<OptimalWindowSumWaveformTreatmentParallelEventVisitor*>(
    antecedent_visitors[low_gain_visitor_]); // good in case of nullptr also
  auto* child = new SimpleChargeStatsParallelEventVisitor(hgv, lgv, config_);
  child->parent_ = this;
  return child;
}

bool SimpleChargeStatsParallelEventVisitor::visit_telescope_run(
  const calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration* run_config,
  calin::iact_data::event_visitor::EventLifetimeManager* event_lifetime_manager)
{
  partials_.Clear();
  results_.Clear();

  has_dual_gain_ = (run_config->camera_layout().adc_gains() !=
    calin::ix::iact_data::instrument_layout::CameraLayout::SINGLE_GAIN);
  for(int ichan=0;ichan<run_config->configured_channel_id_size();++ichan) {
    partials_.add_channel();
  }

  for(auto* h :chan_hists_)delete h;
  chan_hists_.resize(run_config->configured_channel_id_size());
  for(auto*& h : chan_hists_) {
    // For now on do not use per-channel low gain histograms
    h = new ChannelHists(/* has_dual_gain_ */ false, config_.ped_time_hist_resolution(), 3600.0);
  }

  delete camera_hists_;
  camera_hists_ = new ChannelHists(has_dual_gain_, 60.0, 86400.0);

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

void SimpleChargeStatsParallelEventVisitor::integrate_one_gain_partials(
  calin::ix::diagnostics::simple_charge_stats::OneGainSimpleChargeStats* results_g,
  const calin::ix::diagnostics::simple_charge_stats::PartialOneGainChannelSimpleChargeStats& partials_gc)
{
  results_g->add_all_trigger_event_count(partials_gc.all_trig_num_events());
  if(partials_gc.all_trig_num_events() > 0) {
    results_g->add_all_trigger_ped_win_mean(
      double(partials_gc.all_trig_ped_win_sum())/double(partials_gc.all_trig_num_events()));
    results_g->add_all_trigger_ped_win_var(cov_i64_gen(
      partials_gc.all_trig_ped_win_sumsq(), partials_gc.all_trig_num_events(),
      partials_gc.all_trig_ped_win_sum(), partials_gc.all_trig_num_events(),
      partials_gc.all_trig_ped_win_sum(), partials_gc.all_trig_num_events()));
  } else {
    results_g->add_all_trigger_ped_win_mean(0.0);
    results_g->add_all_trigger_ped_win_var(0.0);
  }

  if(partials_gc.has_all_trig_pedwin_vs_time_1_sum()) {
    auto* mean_hist = new calin::ix::math::histogram::Histogram1DData;
    auto* var_hist = new calin::ix::math::histogram::Histogram1DData;
    mean_hist->CopyFrom(partials_gc.all_trig_pedwin_vs_time_q_sum());
    var_hist->CopyFrom(partials_gc.all_trig_pedwin_vs_time_q2_sum());
    for(int ibin=0;ibin<partials_gc.all_trig_pedwin_vs_time_1_sum().bins_size();ibin++) {
      double count = partials_gc.all_trig_pedwin_vs_time_1_sum().bins(ibin);
      if(count>0) {
        mean_hist->set_bins(ibin, mean_hist->bins(ibin)/count);
        var_hist->set_bins(ibin, var_hist->bins(ibin)/count - SQR(mean_hist->bins(ibin)));
      }
    }
    calin::math::histogram::sparsify(partials_gc.all_trig_pedwin_vs_time_1_sum(),
      results_g->add_all_trigger_ped_win_count_vs_time());
    calin::math::histogram::sparsify(*mean_hist,
      results_g->add_all_trigger_ped_win_mean_vs_time());
    calin::math::histogram::sparsify(*var_hist,
      results_g->add_all_trigger_ped_win_var_vs_time());
    delete var_hist;
    delete mean_hist;
  }

  if(partials_gc.has_ped_trig_vs_time_1_sum()) {
    auto* mean_hist = new calin::ix::math::histogram::Histogram1DData;
    auto* var_hist = new calin::ix::math::histogram::Histogram1DData;
    mean_hist->CopyFrom(partials_gc.ped_trig_vs_time_q_sum());
    var_hist->CopyFrom(partials_gc.ped_trig_vs_time_q2_sum());
    for(int ibin=0;ibin<partials_gc.ped_trig_vs_time_1_sum().bins_size();ibin++) {
      double count = partials_gc.ped_trig_vs_time_1_sum().bins(ibin);
      if(count>0) {
        mean_hist->set_bins(ibin, mean_hist->bins(ibin)/count);
        var_hist->set_bins(ibin, var_hist->bins(ibin)/count - SQR(mean_hist->bins(ibin)));
      }
    }
    calin::math::histogram::sparsify(partials_gc.ped_trig_vs_time_1_sum(),
      results_g->add_ped_trigger_full_wf_count_vs_time());
    calin::math::histogram::sparsify(*mean_hist,
      results_g->add_ped_trigger_full_wf_mean_vs_time());
    calin::math::histogram::sparsify(*var_hist,
      results_g->add_ped_trigger_full_wf_var_vs_time());
    delete var_hist;
    delete mean_hist;
  }

  results_g->add_ped_trigger_event_count(partials_gc.ped_trig_num_events());
  if(partials_gc.ped_trig_num_events() > 0) {
    results_g->add_ped_trigger_full_wf_mean(
      double(partials_gc.ped_trig_full_wf_sum())/double(partials_gc.ped_trig_num_events()));
    results_g->add_ped_trigger_full_wf_var(cov_i64_gen(
      partials_gc.ped_trig_full_wf_sumsq(), partials_gc.ped_trig_num_events(),
      partials_gc.ped_trig_full_wf_sum(), partials_gc.ped_trig_num_events(),
      partials_gc.ped_trig_full_wf_sum(), partials_gc.ped_trig_num_events()));

    results_g->add_ped_trigger_ped_win_mean(
      double(partials_gc.ped_trig_ped_win_sum())/double(partials_gc.ped_trig_num_events()));
    results_g->add_ped_trigger_ped_win_var(cov_i64_gen(
      partials_gc.ped_trig_ped_win_sumsq(), partials_gc.ped_trig_num_events(),
      partials_gc.ped_trig_ped_win_sum(), partials_gc.ped_trig_num_events(),
      partials_gc.ped_trig_ped_win_sum(), partials_gc.ped_trig_num_events()));

    results_g->add_ped_trigger_sig_win_mean(
      double(partials_gc.ped_trig_sig_win_sum())/double(partials_gc.ped_trig_num_events()));
    results_g->add_ped_trigger_sig_win_var(cov_i64_gen(
      partials_gc.ped_trig_sig_win_sumsq(), partials_gc.ped_trig_num_events(),
      partials_gc.ped_trig_sig_win_sum(), partials_gc.ped_trig_num_events(),
      partials_gc.ped_trig_sig_win_sum(), partials_gc.ped_trig_num_events()));

    results_g->add_ped_trigger_opt_win_mean(
      double(partials_gc.ped_trig_opt_win_sum())/double(partials_gc.ped_trig_num_events()));
    results_g->add_ped_trigger_opt_win_var(cov_i64_gen(
      partials_gc.ped_trig_opt_win_sumsq(), partials_gc.ped_trig_num_events(),
      partials_gc.ped_trig_opt_win_sum(), partials_gc.ped_trig_num_events(),
      partials_gc.ped_trig_opt_win_sum(), partials_gc.ped_trig_num_events()));
  } else {
    results_g->add_ped_trigger_full_wf_mean(0.0);
    results_g->add_ped_trigger_full_wf_var(0.0);
    results_g->add_ped_trigger_ped_win_mean(0.0);
    results_g->add_ped_trigger_ped_win_var(0.0);
    results_g->add_ped_trigger_sig_win_mean(0.0);
    results_g->add_ped_trigger_sig_win_var(0.0);
    results_g->add_ped_trigger_opt_win_mean(0.0);
    results_g->add_ped_trigger_opt_win_var(0.0);
  }

  results_g->add_ext_trigger_event_count(partials_gc.ext_trig_num_events());
  if(partials_gc.ext_trig_num_events() > 0) {
    results_g->add_ext_trigger_sig_win_mean(
      double(partials_gc.ext_trig_sig_win_sum())/double(partials_gc.ext_trig_num_events()));
    results_g->add_ext_trigger_sig_win_var(cov_i64_gen(
      partials_gc.ext_trig_sig_win_sumsq(), partials_gc.ext_trig_num_events(),
      partials_gc.ext_trig_sig_win_sum(), partials_gc.ext_trig_num_events(),
      partials_gc.ext_trig_sig_win_sum(), partials_gc.ext_trig_num_events()));

    results_g->add_ext_trigger_opt_win_mean(
      double(partials_gc.ext_trig_opt_win_sum())/double(partials_gc.ext_trig_num_events()));
    results_g->add_ext_trigger_opt_win_var(cov_i64_gen(
      partials_gc.ext_trig_opt_win_sumsq(), partials_gc.ext_trig_num_events(),
      partials_gc.ext_trig_opt_win_sum(), partials_gc.ext_trig_num_events(),
      partials_gc.ext_trig_opt_win_sum(), partials_gc.ext_trig_num_events()));
  } else {
    results_g->add_ext_trigger_sig_win_mean(0.0);
    results_g->add_ext_trigger_sig_win_var(0.0);
    results_g->add_ext_trigger_opt_win_mean(0.0);
    results_g->add_ext_trigger_opt_win_var(0.0);
  }
}

void SimpleChargeStatsParallelEventVisitor::integrate_one_gain_camera_partials(
  calin::ix::diagnostics::simple_charge_stats::OneGainSimpleChargeStats* results_g,
  const calin::ix::diagnostics::simple_charge_stats::PartialOneGainCameraSimpleChargeStats& partials_g)
{
  results_g->set_ext_trigger_all_channel_count(partials_g.ext_trig_all_num_events());
  if(partials_g.ext_trig_all_num_events() > 0) {
    results_g->set_ext_trigger_all_channel_sig_win_mean(
      double(partials_g.ext_trig_all_sig_win_sum())/double(partials_g.ext_trig_all_num_events()));
    results_g->set_ext_trigger_all_channel_sig_win_var(cov_double_gen(
      partials_g.ext_trig_all_sig_win_sumsq(), partials_g.ext_trig_all_num_events(),
      partials_g.ext_trig_all_sig_win_sum(), partials_g.ext_trig_all_num_events(),
      partials_g.ext_trig_all_sig_win_sum(), partials_g.ext_trig_all_num_events()));

    results_g->set_ext_trigger_all_channel_opt_win_mean(
      double(partials_g.ext_trig_all_opt_win_sum())/double(partials_g.ext_trig_all_num_events()));
    results_g->set_ext_trigger_all_channel_opt_win_var(cov_double_gen(
      partials_g.ext_trig_all_opt_win_sumsq(), partials_g.ext_trig_all_num_events(),
      partials_g.ext_trig_all_opt_win_sum(), partials_g.ext_trig_all_num_events(),
      partials_g.ext_trig_all_opt_win_sum(), partials_g.ext_trig_all_num_events()));
  } else {
    results_g->add_ext_trigger_sig_win_mean(0.0);
    results_g->add_ext_trigger_sig_win_var(0.0);
    results_g->add_ext_trigger_opt_win_mean(0.0);
    results_g->add_ext_trigger_opt_win_var(0.0);
  }

  if(partials_g.has_all_trig_pedwin_vs_time_1_sum()) {
    auto* mean_hist = new calin::ix::math::histogram::Histogram1DData;
    auto* var_hist = new calin::ix::math::histogram::Histogram1DData;
    mean_hist->CopyFrom(partials_g.all_trig_pedwin_vs_time_q_sum());
    var_hist->CopyFrom(partials_g.all_trig_pedwin_vs_time_q2_sum());
    for(int ibin=0;ibin<partials_g.all_trig_pedwin_vs_time_1_sum().bins_size();ibin++) {
      double count = partials_g.all_trig_pedwin_vs_time_1_sum().bins(ibin);
      if(count>0) {
        mean_hist->set_bins(ibin, mean_hist->bins(ibin)/count);
        var_hist->set_bins(ibin, var_hist->bins(ibin)/count - SQR(mean_hist->bins(ibin)));
      }
    }
    calin::math::histogram::sparsify(partials_g.all_trig_pedwin_vs_time_1_sum(),
      results_g->mutable_camera_all_trigger_ped_win_count_vs_time());
    calin::math::histogram::sparsify(*mean_hist,
      results_g->mutable_camera_all_trigger_ped_win_mean_vs_time());
    calin::math::histogram::sparsify(*var_hist,
      results_g->mutable_camera_all_trigger_ped_win_var_vs_time());
    delete var_hist;
    delete mean_hist;
  }

  if(partials_g.has_ped_trig_vs_time_1_sum()) {
    auto* mean_hist = new calin::ix::math::histogram::Histogram1DData;
    auto* var_hist = new calin::ix::math::histogram::Histogram1DData;
    mean_hist->CopyFrom(partials_g.ped_trig_vs_time_q_sum());
    var_hist->CopyFrom(partials_g.ped_trig_vs_time_q2_sum());
    for(int ibin=0;ibin<partials_g.ped_trig_vs_time_1_sum().bins_size();ibin++) {
      double count = partials_g.ped_trig_vs_time_1_sum().bins(ibin);
      if(count>0) {
        mean_hist->set_bins(ibin, mean_hist->bins(ibin)/count);
        var_hist->set_bins(ibin, var_hist->bins(ibin)/count - SQR(mean_hist->bins(ibin)));
      }
    }
    calin::math::histogram::sparsify(partials_g.ped_trig_vs_time_1_sum(),
      results_g->mutable_camera_ped_trigger_full_wf_count_vs_time());
    calin::math::histogram::sparsify(*mean_hist,
      results_g->mutable_camera_ped_trigger_full_wf_mean_vs_time());
    calin::math::histogram::sparsify(*var_hist,
      results_g->mutable_camera_ped_trigger_full_wf_var_vs_time());
    delete var_hist;
    delete mean_hist;
  }

}

void SimpleChargeStatsParallelEventVisitor::dump_single_gain_channel_hists_to_partials(
  const SingleGainChannelHists& hists,
  calin::ix::diagnostics::simple_charge_stats::PartialOneGainChannelSimpleChargeStats* partials)
{
  auto* hp = hists.all_pedwin_1_sum_vs_time->dump_as_proto();
  partials->mutable_all_trig_pedwin_vs_time_1_sum()->IntegrateFrom(*hp);

  hists.all_pedwin_q_sum_vs_time->dump_as_proto(hp);
  partials->mutable_all_trig_pedwin_vs_time_q_sum()->IntegrateFrom(*hp);

  hists.all_pedwin_q2_sum_vs_time->dump_as_proto(hp);
  partials->mutable_all_trig_pedwin_vs_time_q2_sum()->IntegrateFrom(*hp);

  hists.ped_wf_1_sum_vs_time->dump_as_proto(hp);
  partials->mutable_ped_trig_vs_time_1_sum()->IntegrateFrom(*hp);

  hists.ped_wf_q_sum_vs_time->dump_as_proto(hp);
  partials->mutable_ped_trig_vs_time_q_sum()->IntegrateFrom(*hp);

  hists.ped_wf_q2_sum_vs_time->dump_as_proto(hp);
  partials->mutable_ped_trig_vs_time_q2_sum()->IntegrateFrom(*hp);

  delete hp;
}

void SimpleChargeStatsParallelEventVisitor::dump_single_gain_camera_hists_to_partials(
  const SingleGainChannelHists& hists,
  calin::ix::diagnostics::simple_charge_stats::PartialOneGainCameraSimpleChargeStats* partials)
{
  auto* hp = hists.all_pedwin_1_sum_vs_time->dump_as_proto();
  partials->mutable_all_trig_pedwin_vs_time_1_sum()->IntegrateFrom(*hp);

  hists.all_pedwin_q_sum_vs_time->dump_as_proto(hp);
  partials->mutable_all_trig_pedwin_vs_time_q_sum()->IntegrateFrom(*hp);

  hists.all_pedwin_q2_sum_vs_time->dump_as_proto(hp);
  partials->mutable_all_trig_pedwin_vs_time_q2_sum()->IntegrateFrom(*hp);

  hists.ped_wf_1_sum_vs_time->dump_as_proto(hp);
  partials->mutable_ped_trig_vs_time_1_sum()->IntegrateFrom(*hp);

  hists.ped_wf_q_sum_vs_time->dump_as_proto(hp);
  partials->mutable_ped_trig_vs_time_q_sum()->IntegrateFrom(*hp);

  hists.ped_wf_q2_sum_vs_time->dump_as_proto(hp);
  partials->mutable_ped_trig_vs_time_q2_sum()->IntegrateFrom(*hp);

  delete hp;
}

bool SimpleChargeStatsParallelEventVisitor::leave_telescope_run()
{
  for(int ichan = 0; ichan<partials_.channel_size(); ichan++) {
    dump_single_gain_channel_hists_to_partials(*chan_hists_[ichan]->high_gain,
      partials_.mutable_channel(ichan)->mutable_high_gain());
    delete chan_hists_[ichan]->high_gain;
    chan_hists_[ichan]->high_gain = nullptr;
    if(has_dual_gain_ and chan_hists_[ichan]->low_gain) {
      dump_single_gain_channel_hists_to_partials(*chan_hists_[ichan]->low_gain,
        partials_.mutable_channel(ichan)->mutable_low_gain());
      delete chan_hists_[ichan]->low_gain;
      chan_hists_[ichan]->low_gain = nullptr;
    }
  }

  dump_single_gain_camera_hists_to_partials(*camera_hists_->high_gain,
    partials_.mutable_camera()->mutable_high_gain());
  if(has_dual_gain_ and camera_hists_->low_gain) {
    dump_single_gain_camera_hists_to_partials(*camera_hists_->low_gain,
      partials_.mutable_camera()->mutable_low_gain());
  }

  if(parent_)return true;

  for(int ichan = 0; ichan<partials_.channel_size(); ichan++) {
    auto& partials_chan = partials_.channel(ichan);
    integrate_one_gain_partials(results_.mutable_high_gain(), partials_chan.high_gain());
    if(has_dual_gain_) {
      integrate_one_gain_partials(results_.mutable_low_gain(), partials_chan.low_gain());
    }
    if(partials_.camera().num_event_trigger_hitmap_found() > 0) {
      results_.add_channel_triggered_count(partials_chan.all_trig_num_events_triggered());
    }
  }
  integrate_one_gain_camera_partials(results_.mutable_high_gain(), partials_.camera().high_gain());
  if(has_dual_gain_) {
    integrate_one_gain_camera_partials(results_.mutable_low_gain(), partials_.camera().low_gain());
  }
  partials_.Clear();
  return true;
}

void SimpleChargeStatsParallelEventVisitor::record_one_gain_channel_data(
  const calin::ix::iact_data::telescope_event::TelescopeEvent* event,
  const calin::iact_data::waveform_treatment_event_visitor::OptimalWindowSumWaveformTreatmentParallelEventVisitor* sum_visitor,
  unsigned ichan, double elapsed_event_time,
  calin::ix::diagnostics::simple_charge_stats::PartialOneGainChannelSimpleChargeStats* one_gain_stats,
  SingleGainChannelHists* one_gain_hists,
  unsigned& nsum, int64_t& opt_sum, int64_t& sig_sum, int64_t& bkg_sum, int64_t& wf_sum)
{
  one_gain_stats->increment_all_trig_num_events();

  int64_t ped_win_sum = sum_visitor->array_chan_bkg_win_sum()[ichan];
  int64_t sqr_ped_win_sum = SQR(ped_win_sum);
  one_gain_stats->increment_all_trig_ped_win_sum(ped_win_sum);
  one_gain_stats->increment_all_trig_ped_win_sumsq(sqr_ped_win_sum);
  if(one_gain_hists) {
    one_gain_hists->all_pedwin_1_sum_vs_time->insert(elapsed_event_time);
    one_gain_hists->all_pedwin_q_sum_vs_time->insert(elapsed_event_time, ped_win_sum);
    one_gain_hists->all_pedwin_q2_sum_vs_time->insert(elapsed_event_time, sqr_ped_win_sum);
  }
  int64_t wf_all_sum = sum_visitor->array_chan_all_sum()[ichan];
  int64_t sig_win_sum = sum_visitor->array_chan_sig_win_sum()[ichan];
  int64_t bkg_win_sum = sum_visitor->array_chan_bkg_win_sum()[ichan];
  int64_t opt_win_sum = sum_visitor->array_chan_opt_win_sum()[ichan];
  if(event->trigger_type() == calin::ix::iact_data::telescope_event::TRIGGER_PEDESTAL) {
    int64_t sqr_wf_all_sum = SQR(wf_all_sum);
    one_gain_stats->increment_ped_trig_num_events();
    one_gain_stats->increment_ped_trig_full_wf_sum(wf_all_sum);
    one_gain_stats->increment_ped_trig_full_wf_sumsq(sqr_wf_all_sum);
    one_gain_stats->increment_ped_trig_ped_win_sum(ped_win_sum);
    one_gain_stats->increment_ped_trig_ped_win_sumsq(sqr_ped_win_sum);
    int64_t sqr_sig_win_sum = SQR(sig_win_sum);
    one_gain_stats->increment_ped_trig_sig_win_sum(sig_win_sum);
    one_gain_stats->increment_ped_trig_sig_win_sumsq(sqr_sig_win_sum);
    int64_t sqr_opt_win_sum = SQR(opt_win_sum);
    one_gain_stats->increment_ped_trig_opt_win_sum(opt_win_sum);
    one_gain_stats->increment_ped_trig_opt_win_sumsq(sqr_opt_win_sum);
    if(one_gain_hists) {
      one_gain_hists->ped_wf_1_sum_vs_time->insert(elapsed_event_time);
      one_gain_hists->ped_wf_q_sum_vs_time->insert(elapsed_event_time, wf_all_sum);
      one_gain_hists->ped_wf_q2_sum_vs_time->insert(elapsed_event_time, sqr_wf_all_sum);
    }
  } else if(event->trigger_type() == calin::ix::iact_data::telescope_event::TRIGGER_EXTERNAL_FLASHER) {
    int64_t sqr_sig_win_sum = SQR(sig_win_sum);
    one_gain_stats->increment_ext_trig_num_events();
    one_gain_stats->increment_ext_trig_sig_win_sum(sig_win_sum);
    one_gain_stats->increment_ext_trig_sig_win_sumsq(sqr_sig_win_sum);
    int64_t sqr_opt_win_sum = SQR(opt_win_sum);
    one_gain_stats->increment_ext_trig_opt_win_sum(opt_win_sum);
    one_gain_stats->increment_ext_trig_opt_win_sumsq(sqr_opt_win_sum);
  } else if(event->trigger_type() == calin::ix::iact_data::telescope_event::TRIGGER_PHYSICS) {
  }
  ++nsum;
  opt_sum += opt_win_sum;
  sig_sum += sig_win_sum;
  bkg_sum += bkg_win_sum;
  wf_sum += wf_all_sum;
}

void SimpleChargeStatsParallelEventVisitor::record_one_visitor_data(
  uint64_t seq_index, const calin::ix::iact_data::telescope_event::TelescopeEvent* event,
  const calin::iact_data::waveform_treatment_event_visitor::OptimalWindowSumWaveformTreatmentParallelEventVisitor* sum_visitor,
  calin::ix::diagnostics::simple_charge_stats::PartialSimpleChargeStats* partials)
{
  double elapsed_event_time = event->elapsed_event_time().time_ns() * 1e-9;

  if(sum_visitor and sum_visitor->is_same_event(seq_index)) {
    unsigned nsum_hg = 0;
    int64_t opt_sum_hg = 0;
    int64_t sig_sum_hg = 0;
    int64_t bkg_sum_hg = 0;
    int64_t all_sum_hg = 0;

    unsigned nsum_lg = 0;
    int64_t opt_sum_lg = 0;
    int64_t sig_sum_lg = 0;
    int64_t bkg_sum_lg = 0;
    int64_t all_sum_lg = 0;
    for(unsigned ichan=0; ichan<sum_visitor->nchan(); ichan++) {
      auto* pc = partials->mutable_channel(ichan);
      switch(sum_visitor->array_chan_signal_type()[ichan]) {
      case calin::ix::iact_data::telescope_event::SIGNAL_UNIQUE_GAIN:
      case calin::ix::iact_data::telescope_event::SIGNAL_HIGH_GAIN:
        record_one_gain_channel_data(event, sum_visitor, ichan, elapsed_event_time,
          pc->mutable_high_gain(), chan_hists_[ichan]->high_gain,
          nsum_hg, opt_sum_hg, sig_sum_hg, bkg_sum_hg, all_sum_hg);
        break;
      case calin::ix::iact_data::telescope_event::SIGNAL_LOW_GAIN:
        record_one_gain_channel_data(event, sum_visitor, ichan, elapsed_event_time,
          pc->mutable_low_gain(), chan_hists_[ichan]->low_gain,
          nsum_lg, opt_sum_lg, sig_sum_lg, bkg_sum_lg, all_sum_lg);
        break;
      case calin::ix::iact_data::telescope_event::SIGNAL_NONE:
      default:
        // do nothing
        break;
      }
    }
    if(nsum_hg == sum_visitor->nchan()) {
      camera_hists_->high_gain->all_pedwin_1_sum_vs_time->insert(elapsed_event_time);
      camera_hists_->high_gain->all_pedwin_q_sum_vs_time->insert(elapsed_event_time, bkg_sum_hg);
      camera_hists_->high_gain->all_pedwin_q2_sum_vs_time->insert(elapsed_event_time, bkg_sum_hg*bkg_sum_hg);
      if(event->trigger_type() == calin::ix::iact_data::telescope_event::TRIGGER_EXTERNAL_FLASHER) {
        auto* pcam = partials->mutable_camera()->mutable_high_gain();
        pcam->increment_ext_trig_all_num_events();
        pcam->increment_ext_trig_all_sig_win_sum(sig_sum_hg);
        pcam->increment_ext_trig_all_sig_win_sumsq(SQR(sig_sum_hg));
        pcam->increment_ext_trig_all_opt_win_sum(opt_sum_hg);
        pcam->increment_ext_trig_all_opt_win_sumsq(SQR(opt_sum_hg));
      } else if(event->trigger_type() == calin::ix::iact_data::telescope_event::TRIGGER_PEDESTAL) {
        camera_hists_->high_gain->ped_wf_1_sum_vs_time->insert(elapsed_event_time);
        camera_hists_->high_gain->ped_wf_q_sum_vs_time->insert(elapsed_event_time, all_sum_hg);
        camera_hists_->high_gain->ped_wf_q2_sum_vs_time->insert(elapsed_event_time, all_sum_hg*all_sum_hg);
      }
    }
    if(nsum_lg == sum_visitor->nchan()) {
      camera_hists_->low_gain->all_pedwin_1_sum_vs_time->insert(elapsed_event_time);
      camera_hists_->low_gain->all_pedwin_q_sum_vs_time->insert(elapsed_event_time, bkg_sum_lg);
      camera_hists_->low_gain->all_pedwin_q2_sum_vs_time->insert(elapsed_event_time, bkg_sum_lg*bkg_sum_lg);
      if(event->trigger_type() == calin::ix::iact_data::telescope_event::TRIGGER_EXTERNAL_FLASHER) {
        auto* pcam = partials->mutable_camera()->mutable_low_gain();
        pcam->increment_ext_trig_all_num_events();
        pcam->increment_ext_trig_all_sig_win_sum(sig_sum_lg);
        pcam->increment_ext_trig_all_sig_win_sumsq(SQR(sig_sum_lg));
        pcam->increment_ext_trig_all_opt_win_sum(opt_sum_lg);
        pcam->increment_ext_trig_all_opt_win_sumsq(SQR(opt_sum_lg));
      } else if(event->trigger_type() == calin::ix::iact_data::telescope_event::TRIGGER_PEDESTAL) {
        camera_hists_->low_gain->ped_wf_1_sum_vs_time->insert(elapsed_event_time);
        camera_hists_->low_gain->ped_wf_q_sum_vs_time->insert(elapsed_event_time, all_sum_lg);
        camera_hists_->low_gain->ped_wf_q2_sum_vs_time->insert(elapsed_event_time, all_sum_lg*all_sum_lg);
      }
    }
  }
}

bool SimpleChargeStatsParallelEventVisitor::visit_telescope_event(uint64_t seq_index,
  calin::ix::iact_data::telescope_event::TelescopeEvent* event)
{

  if(high_gain_visitor_) {
    record_one_visitor_data(seq_index, event, high_gain_visitor_, &partials_);
  }
  if(low_gain_visitor_) {
    record_one_visitor_data(seq_index, event, low_gain_visitor_, &partials_);
  }
  if(event->has_trigger_map()) {
    if(event->trigger_map().hit_channel_id_size()>0) {
      partials_.mutable_camera()->increment_num_event_trigger_hitmap_found();
    }
    for(auto ichan : event->trigger_map().hit_channel_id()) {
      partials_.mutable_channel(ichan)->increment_all_trig_num_events_triggered();
    }
  }
  return true;
}

bool SimpleChargeStatsParallelEventVisitor::merge_results()
{
  if(parent_) {
    parent_->partials_.IntegrateFrom(partials_);
  }
  return true;
}

calin::ix::diagnostics::simple_charge_stats::SimpleChargeStats*
SimpleChargeStatsParallelEventVisitor::simple_charge_stats(
  calin::ix::diagnostics::simple_charge_stats::SimpleChargeStats* stats) const
{
  if(stats == nullptr)stats = results_.New();
  stats->CopyFrom(results_);
  return stats;
}

calin::ix::diagnostics::simple_charge_stats::SimpleChargeStatsConfig
SimpleChargeStatsParallelEventVisitor::default_config()
{
  calin::ix::diagnostics::simple_charge_stats::SimpleChargeStatsConfig config;
  config.set_ped_time_hist_resolution(5.0);
  return config;
}
