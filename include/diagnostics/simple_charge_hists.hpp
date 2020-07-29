/*

   calin/diagnostics/simple_charge_hists.hpp -- Stephen Fegan -- 2020-04-17

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

#pragma once

#include <math/histogram.hpp>
#include <iact_data/event_visitor.hpp>
#include <iact_data/waveform_treatment_event_visitor.hpp>
#include <diagnostics/simple_charge_hists.pb.h>

namespace calin { namespace diagnostics { namespace simple_charge_hists {

class SimpleChargeHistsParallelEventVisitor:
  public calin::iact_data::event_visitor::ParallelEventVisitor
{
public:
  SimpleChargeHistsParallelEventVisitor(
    calin::iact_data::waveform_treatment_event_visitor::OptimalWindowSumWaveformTreatmentParallelEventVisitor* high_gain_visitor,
    calin::iact_data::waveform_treatment_event_visitor::OptimalWindowSumWaveformTreatmentParallelEventVisitor* low_gain_visitor,
    const calin::ix::diagnostics::simple_charge_hists::SimpleChargeHistsConfig& config = ped_trig_default_config());

  SimpleChargeHistsParallelEventVisitor(
      calin::iact_data::waveform_treatment_event_visitor::OptimalWindowSumWaveformTreatmentParallelEventVisitor* mixed_or_unique_gain_visitor,
      const calin::ix::diagnostics::simple_charge_hists::SimpleChargeHistsConfig& config = ped_trig_default_config()):
    SimpleChargeHistsParallelEventVisitor(mixed_or_unique_gain_visitor, nullptr, config)
  { /* nothing to see here */ }

  virtual ~SimpleChargeHistsParallelEventVisitor();

  SimpleChargeHistsParallelEventVisitor* new_sub_visitor(
    std::map<calin::iact_data::event_visitor::ParallelEventVisitor*,
        calin::iact_data::event_visitor::ParallelEventVisitor*>
      antecedent_visitors = { }) override;

  bool visit_telescope_run(
    const calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration* run_config,
    calin::iact_data::event_visitor::EventLifetimeManager* event_lifetime_manager) override;

  bool leave_telescope_run() override;

  bool visit_telescope_event(uint64_t seq_index,
    calin::ix::iact_data::telescope_event::TelescopeEvent* event) override;

  bool merge_results() override;

#ifndef SWIG
  calin::ix::diagnostics::simple_charge_hists::SimpleChargeHists* simple_charge_hists(
    calin::ix::diagnostics::simple_charge_hists::SimpleChargeHists* stats = nullptr) const;
#else
  calin::ix::diagnostics::simple_charge_hists::SimpleChargeHists* simple_charge_hists() const;
  void simple_charge_hists(calin::ix::diagnostics::simple_charge_hists::SimpleChargeHists* stats) const;
#endif

  static calin::ix::diagnostics::simple_charge_hists::SimpleChargeHistsConfig all_enabled_config(
    double qsum_dx=10.0, int qsum_nmax=100, int qsum_rebin=5,
    double max_dx=10.0, int max_nmax=100, int max_rebin=5);
  static calin::ix::diagnostics::simple_charge_hists::SimpleChargeHistsConfig ped_trig_default_config();
  static calin::ix::diagnostics::simple_charge_hists::SimpleChargeHistsConfig phy_trig_default_config();
  static calin::ix::diagnostics::simple_charge_hists::SimpleChargeHistsConfig ext_trig_default_config();
  static calin::ix::diagnostics::simple_charge_hists::SimpleChargeHistsConfig int_trig_default_config();

  // const calin::ix::diagnostics::simple_charge_stats::PartialSimpleChargeStats& partials() const { return partials_; }

private:
  struct SingleGainChannelHists {
    SingleGainChannelHists(const
        calin::ix::diagnostics::simple_charge_hists::SingleGainSimpleChargeHistsConfig& config):
      full_wf_qsum(calin::math::histogram::new_histogram_if_enabled(config.full_wf_qsum())),
      full_wf_max(calin::math::histogram::new_histogram_if_enabled(config.full_wf_max())),
      full_wf_max_index(calin::math::histogram::new_histogram_if_enabled(config.full_wf_max_index())),
      opt_win_qsum(calin::math::histogram::new_histogram_if_enabled(config.opt_win_qsum())),
      opt_win_qtsum(calin::math::histogram::new_histogram_if_enabled(config.opt_win_qtsum())),
      opt_win_index(calin::math::histogram::new_histogram_if_enabled(config.opt_win_index())),
      ped_win_qsum(calin::math::histogram::new_histogram_if_enabled(config.ped_win_qsum())),
      sig_win_qsum(calin::math::histogram::new_histogram_if_enabled(config.sig_win_qsum()))
    {
      /* nothing to see here */
    }

    ~SingleGainChannelHists() {
      delete full_wf_qsum;
      delete full_wf_max;
      delete full_wf_max_index;
      delete opt_win_qsum;
      delete opt_win_qtsum;
      delete opt_win_index;
      delete ped_win_qsum;
      delete sig_win_qsum;
    }

    calin::math::histogram::Histogram1D* full_wf_qsum = nullptr;
    calin::math::histogram::Histogram1D* full_wf_max = nullptr;
    calin::math::histogram::Histogram1D* full_wf_max_index = nullptr;
    calin::math::histogram::Histogram1D* opt_win_qsum = nullptr;
    calin::math::histogram::Histogram1D* opt_win_qtsum = nullptr;
    calin::math::histogram::Histogram1D* opt_win_index = nullptr;
    calin::math::histogram::Histogram1D* ped_win_qsum = nullptr;
    calin::math::histogram::Histogram1D* sig_win_qsum = nullptr;
  };

  struct DualGainChannelHists {
    DualGainChannelHists(const calin::ix::diagnostics::simple_charge_hists::SimpleChargeHistsConfig& config)
    { /* nothing to see here */ }

    ~DualGainChannelHists() {
    }
  };

  struct SingleGainCameraHists {
    SingleGainCameraHists(const
        calin::ix::diagnostics::simple_charge_hists::SingleGainSimpleChargeHistsConfig& config):
      nchan_present(calin::math::histogram::new_histogram_if_enabled(config.nchan_present()))
    {
      /* nothing to see here */
    }

    ~SingleGainCameraHists() {
      delete nchan_present;
    }

    calin::math::histogram::Histogram1D* nchan_present = nullptr;
  };

  struct DualGainCameraHists {
    DualGainCameraHists(const
        calin::ix::diagnostics::simple_charge_hists::DualGainSimpleChargeHistsConfig& config):
      nchan_present(calin::math::histogram::new_histogram_if_enabled(config.nchan_present()))
    {
      /* nothing to see here */
    }

    ~DualGainCameraHists() {
      delete nchan_present;
    }

    calin::math::histogram::Histogram1D* nchan_present = nullptr;
  };

  struct ChannelHists {
    ChannelHists(bool has_dual_gain_,
        const calin::ix::diagnostics::simple_charge_hists::SimpleChargeHistsConfig& config):
      high_gain((config.has_high_gain() and config.high_gain().enable_hists()) ?
        new SingleGainChannelHists(config.high_gain()) : nullptr),
      low_gain((has_dual_gain_ and config.has_low_gain() and config.low_gain().enable_hists()) ?
        new SingleGainChannelHists(config.low_gain()) : nullptr),
      dual_gain(has_dual_gain_ ? new DualGainChannelHists(config) : nullptr)
    { /* nothing to see here */ }

    ~ChannelHists() {
      delete high_gain;
      delete low_gain;
      delete dual_gain;
    }

    SingleGainChannelHists* high_gain = nullptr;
    SingleGainChannelHists* low_gain = nullptr;
    DualGainChannelHists* dual_gain = nullptr;
  };

  void merge_one_gain_hists(SingleGainChannelHists* into, const SingleGainChannelHists* from);
  void merge_one_gain_cam_hists(SingleGainCameraHists* into, const SingleGainCameraHists* from);
  void merge_dual_gain_cam_hists(DualGainCameraHists* into, const DualGainCameraHists* from);

  void extract_one_gain_hists(calin::ix::diagnostics::simple_charge_hists::OneGainSimpleChargeHists* into,
    const SingleGainChannelHists* from) const;
  void extract_one_gain_cam_hists(calin::ix::diagnostics::simple_charge_hists::OneGainSimpleChargeCameraHists* into,
    const SingleGainCameraHists* from) const;
  void extract_dual_gain_cam_hists(calin::ix::diagnostics::simple_charge_hists::DualGainSimpleChargeCameraHists* into,
    const DualGainCameraHists* from) const;

  void record_one_visitor_data(uint64_t seq_index, const calin::ix::iact_data::telescope_event::TelescopeEvent* event,
    const calin::iact_data::waveform_treatment_event_visitor::OptimalWindowSumWaveformTreatmentParallelEventVisitor* sum_visitor,
    unsigned& high_gain_nchan_presence, unsigned& low_gain_nchan_presence);

  SimpleChargeHistsParallelEventVisitor* parent_ = nullptr;
  calin::ix::diagnostics::simple_charge_hists::SimpleChargeHistsConfig config_;
  calin::iact_data::waveform_treatment_event_visitor::OptimalWindowSumWaveformTreatmentParallelEventVisitor* high_gain_visitor_ = nullptr;
  calin::iact_data::waveform_treatment_event_visitor::OptimalWindowSumWaveformTreatmentParallelEventVisitor* low_gain_visitor_ = nullptr;

  bool has_dual_gain_ = false;

  std::vector<ChannelHists*> chan_hists_;
  SingleGainCameraHists* cam_hists_high_gain_ = nullptr;
  SingleGainCameraHists* cam_hists_low_gain_ = nullptr;
  DualGainCameraHists* cam_hists_dual_gain_ = nullptr;
};

} } } // namespace calin::diagnostics::simple_charge_stats
