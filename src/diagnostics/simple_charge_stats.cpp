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

SimpleChargeStatsVisitor::
SimpleChargeStatsVisitor(
    OptimalWindowSumWaveformTreatmentEventVisitor* high_gain_visitor,
    OptimalWindowSumWaveformTreatmentEventVisitor* low_gain_visitor,
    const SimpleChargeStatsConfig& config):
  ParallelEventVisitor(), config_(config),
  high_gain_visitor_(high_gain_visitor), low_gain_visitor_(low_gain_visitor)
{
  // nothing to see here
}

SimpleChargeStatsVisitor::~SimpleChargeStatsVisitor()
{
  // nothing to see here
}

SimpleChargeStatsVisitor* SimpleChargeStatsVisitor::new_sub_visitor(
  std::map<calin::iact_data::event_visitor::ParallelEventVisitor*,
    calin::iact_data::event_visitor::ParallelEventVisitor*> antecedent_visitors)
{
  auto* hgv = dynamic_cast<OptimalWindowSumWaveformTreatmentEventVisitor*>(
    antecedent_visitors[high_gain_visitor_]);
  auto* lgv = dynamic_cast<OptimalWindowSumWaveformTreatmentEventVisitor*>(
    antecedent_visitors[low_gain_visitor_]); // good in case of nullptr also
  auto* child = new SimpleChargeStatsVisitor(hgv, lgv, config_);
  child->parent_ = this;
  return child;
}

bool SimpleChargeStatsVisitor::visit_telescope_run(
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
  void integrate_one_gain_partials(
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
}

bool SimpleChargeStatsVisitor::leave_telescope_run()
{
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

namespace {
  void record_one_gain_channel_data(const calin::ix::iact_data::telescope_event::TelescopeEvent* event,
    const calin::iact_data::waveform_treatment_event_visitor::OptimalWindowSumWaveformTreatmentEventVisitor* sum_visitor,
    unsigned ichan,
    calin::ix::diagnostics::simple_charge_stats::PartialOneGainChannelSimpleChargeStats* one_gain_stats)
  {
    one_gain_stats->increment_all_trig_num_events();
    one_gain_stats->increment_all_trig_ped_win_sum(sum_visitor->array_chan_bkg_win_sum()[ichan]);
    one_gain_stats->increment_all_trig_ped_win_sumsq(SQR(sum_visitor->array_chan_bkg_win_sum()[ichan]));
    if(event->trigger_type() == calin::ix::iact_data::telescope_event::TRIGGER_PEDESTAL) {
      one_gain_stats->increment_ped_trig_num_events();
      one_gain_stats->increment_ped_trig_full_wf_sum(sum_visitor->array_chan_all_sum()[ichan]);
      one_gain_stats->increment_ped_trig_full_wf_sumsq(SQR(sum_visitor->array_chan_all_sum  ()[ichan]));
    }
  }

  void record_one_visitor_data(uint64_t seq_index, const calin::ix::iact_data::telescope_event::TelescopeEvent* event,
    const calin::iact_data::waveform_treatment_event_visitor::OptimalWindowSumWaveformTreatmentEventVisitor* sum_visitor,
    calin::ix::diagnostics::simple_charge_stats::PartialSimpleChargeStats* partials)
  {
    if(sum_visitor and sum_visitor->is_same_event(seq_index)) {
      for(unsigned ichan=0; ichan<sum_visitor->nchan(); ichan++) {
        auto* pc = partials->mutable_channel(ichan);
        switch(sum_visitor->array_chan_signal_type()[ichan]) {
        case calin::ix::iact_data::telescope_event::SIGNAL_UNIQUE_GAIN:
        case calin::ix::iact_data::telescope_event::SIGNAL_HIGH_GAIN:
          record_one_gain_channel_data(event, sum_visitor, ichan, pc->mutable_high_gain());
          break;
        case calin::ix::iact_data::telescope_event::SIGNAL_LOW_GAIN:
          record_one_gain_channel_data(event, sum_visitor, ichan, pc->mutable_low_gain());
          break;
        case calin::ix::iact_data::telescope_event::SIGNAL_NONE:
        default:
          // do nothing
          break;
        }
      }
    }
  }
} // anonymous namespace

bool SimpleChargeStatsVisitor::visit_telescope_event(uint64_t seq_index,
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

bool SimpleChargeStatsVisitor::merge_results()
{
  if(parent_) {
    parent_->partials_.IntegrateFrom(partials_);
  }
  return true;
}

calin::ix::diagnostics::simple_charge_stats::SimpleChargeStats*
SimpleChargeStatsVisitor::simple_charge_stats() const
{
  auto* res = results_.New();
  res->CopyFrom(results_);
  return res;
}

calin::ix::diagnostics::simple_charge_stats::SimpleChargeStatsConfig
SimpleChargeStatsVisitor::default_config()
{
  calin::ix::diagnostics::simple_charge_stats::SimpleChargeStatsConfig config;
  return config;
}
