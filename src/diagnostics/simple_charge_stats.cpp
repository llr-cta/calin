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
  for(unsigned ichan=0;ichan<run_config->configured_channel_id_size();++ichan) {
    partials_.add_channel();
  }

  for(auto* h :chan_hists_)delete h;
  chan_hists_.resize(run_config->configured_channel_id_size());
  for(auto*& h : chan_hists_) {
    h = new ChannelHists(has_dual_gain_, config_.ped_time_hist_resolution());
  }

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

namespace {
  void transfer_histogram_with_rebin_if_necessary(
    const calin::ix::math::histogram::Histogram1DData& from_hist,
    calin::ix::math::histogram::Histogram1DData* to_hist,
    int maximum_size, int maximum_rebin)
  {
    auto* hs = calin::math::histogram::sparsify(from_hist);
    int rebin = 1;
    if(maximum_size>0 and hs->bins_size()>maximum_size) {
      rebin = (hs->bins_size()+maximum_size-1)/maximum_size;
      if(maximum_rebin>0) {
        rebin = std::min(rebin, maximum_rebin);
      }
    }
    if(rebin) {
      calin::math::histogram::rebin(*hs, rebin, to_hist);
    } else {
      to_hist->CopyFrom(*hs);
    }
    delete hs;
  }
} // anonymous namespace

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

  if(partials_gc.has_ped_trig_full_wf_hist()) {
    transfer_histogram_with_rebin_if_necessary(partials_gc.ped_trig_full_wf_hist(),
      results_g->add_ped_trigger_full_wf_hist(),
      config_.max_ped_hist_bins(), config_.max_ped_hist_rebin_factor());
  }

  if(partials_gc.has_ext_trig_opt_sum_hist()) {
    transfer_histogram_with_rebin_if_necessary(partials_gc.ext_trig_opt_sum_hist(),
      results_g->add_ext_trig_opt_sum_hist(),
      config_.max_ped_hist_bins(), config_.max_ped_hist_rebin_factor());
  }

  if(partials_gc.has_ext_trig_opt_max_hist()) {
    transfer_histogram_with_rebin_if_necessary(partials_gc.ext_trig_opt_max_hist(),
      results_g->add_ext_trig_opt_max_hist(),
      config_.max_ped_hist_bins(), config_.max_ped_hist_rebin_factor());
  }

  if(partials_gc.has_ext_trig_opt_index_hist()) {
    results_g->add_ext_trig_opt_index_hist()->CopyFrom(
      partials_gc.ext_trig_opt_index_hist());
  }

  if(partials_gc.has_ped_trig_vs_time_1_sum()) {
    auto* mean_hist = new calin::ix::math::histogram::Histogram1DData;
    auto* var_hist = new calin::ix::math::histogram::Histogram1DData;
    mean_hist->CopyFrom(partials_gc.ped_trig_vs_time_q_sum());
    var_hist->CopyFrom(partials_gc.ped_trig_vs_time_q2_sum());
    for(unsigned ibin=0;ibin<partials_gc.ped_trig_vs_time_1_sum().bins_size();ibin++) {
      double count = partials_gc.ped_trig_vs_time_1_sum().bins(ibin);
      if(count>0) {
        mean_hist->set_bins(ibin, mean_hist->bins(ibin)/count);
        var_hist->set_bins(ibin, var_hist->bins(ibin)/count - SQR(mean_hist->bins(ibin)));
      }
    }
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
  } else {
    results_g->add_ped_trigger_full_wf_mean(0.0);
    results_g->add_ped_trigger_full_wf_var(0.0);
  }
}

void SimpleChargeStatsParallelEventVisitor::dump_single_gain_channel_hists_to_partials(
  const SingleGainChannelHists& hists,
  calin::ix::diagnostics::simple_charge_stats::PartialOneGainChannelSimpleChargeStats* partials)
{
  auto* hp = hists.ped_wf_q_sum->dump_as_proto();
  partials->mutable_ped_trig_full_wf_hist()->IntegrateFrom(*hp);

  hists.ext_opt_q_sum->dump_as_proto(hp);
  partials->mutable_ext_trig_opt_sum_hist()->IntegrateFrom(*hp);

  hists.ext_opt_index->dump_as_proto(hp);
  partials->mutable_ext_trig_opt_index_hist()->IntegrateFrom(*hp);

  hists.ext_opt_max->dump_as_proto(hp);
  partials->mutable_ext_trig_opt_max_hist()->IntegrateFrom(*hp);

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
    if(has_dual_gain_) {
      dump_single_gain_channel_hists_to_partials(*chan_hists_[ichan]->low_gain,
        partials_.mutable_channel(ichan)->mutable_low_gain());
    }
  }

  if(parent_)return true;

  for(int ichan = 0; ichan<partials_.channel_size(); ichan++) {
    auto& partials_chan = partials_.channel(ichan);
    integrate_one_gain_partials(results_.mutable_high_gain(), partials_chan.high_gain());
    if(has_dual_gain_) {
      integrate_one_gain_partials(results_.mutable_low_gain(), partials_chan.low_gain());
    }
  }

  partials_.Clear();
  return true;
}

void SimpleChargeStatsParallelEventVisitor::record_one_gain_channel_data(
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
    one_gain_hists->ext_opt_q_sum->insert(sum_visitor->array_chan_opt_win_sum()[ichan]);
    one_gain_hists->ext_opt_index->insert(sum_visitor->array_chan_opt_win_index()[ichan]);
    one_gain_hists->ext_opt_max->insert(sum_visitor->array_chan_max()[ichan]);
  }
}

void SimpleChargeStatsParallelEventVisitor::record_one_visitor_data(
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

bool SimpleChargeStatsParallelEventVisitor::visit_telescope_event(uint64_t seq_index,
  calin::ix::iact_data::telescope_event::TelescopeEvent* event)
{
  if(high_gain_visitor_) {
    record_one_visitor_data(seq_index, event, high_gain_visitor_, &partials_);
  }
  if(low_gain_visitor_) {
    record_one_visitor_data(seq_index, event, low_gain_visitor_, &partials_);
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
  config.set_max_ped_hist_bins(1000);
  config.set_max_ped_hist_rebin_factor(0);
  config.set_ped_time_hist_resolution(30.0);
  return config;
}
