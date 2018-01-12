/*

   calin/iact_data/waveform_treatment_event_visitor.hpp -- Stephen Fegan -- 2018-01-11

   Waveform treatment event data visitor - process waveforms

   Copyright 2018, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, Ecole Polytechnique, CNRS/IN2P3

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

#include <string>

#include <iact_data/event_visitor.hpp>
#include <math/simd.hpp>
#include <iact_data/waveform_treatment_event_visitor.pb.h>

namespace calin { namespace iact_data { namespace waveform_treatment_event_visitor {

class SingleGainDualWindowWaveformTreatmentEventVisitor:
  public calin::iact_data::event_visitor::TelescopeEventVisitor
{
public:
  SingleGainDualWindowWaveformTreatmentEventVisitor(
    calin::ix::iact_data::waveform_treatment_event_visitor::
      SingleGainDualWindowWaveformTreatmentEventVisitorConfig config = default_config(),
    bool treat_high_gain = true);
  virtual ~SingleGainDualWindowWaveformTreatmentEventVisitor();

  virtual bool demand_waveforms();
  virtual bool is_parallelizable();
  virtual SingleGainDualWindowWaveformTreatmentEventVisitor* new_sub_visitor(
    const std::map<calin::iact_data::event_visitor::TelescopeEventVisitor*,
        calin::iact_data::event_visitor::TelescopeEventVisitor*>&
      antecedent_visitors = { });

  virtual bool visit_telescope_run(
    const calin::ix::iact_data::telescope_run_configuration::
      TelescopeRunConfiguration* run_config);
  virtual bool leave_telescope_run();

  virtual bool visit_telescope_event(uint64_t seq_index,
    calin::ix::iact_data::telescope_event::TelescopeEvent* event);
  virtual bool leave_telescope_event();

  // This function allows testing of the code
  static void analyze_waveforms(const uint16_t* data, unsigned nchan, unsigned nsamp,
    float* ped, float ped_iir_old, float ped_iir_new,
    int* chan_max_index, int* chan_max, int* chan_bkg_win_sum, int* chan_sig_win_sum,
    int* chan_all_sum_q, int* chan_all_sum_qt, float* all_sig, float* all_mean_t);

  static calin::ix::iact_data::waveform_treatment_event_visitor::
    SingleGainDualWindowWaveformTreatmentEventVisitorConfig default_config()
  {
    calin::ix::iact_data::waveform_treatment_event_visitor::
      SingleGainDualWindowWaveformTreatmentEventVisitorConfig cfg;
    cfg.set_integration_n(16);
    return cfg;
  }

protected:
  calin::ix::iact_data::waveform_treatment_event_visitor::
    SingleGainDualWindowWaveformTreatmentEventVisitorConfig config_;
  bool treat_high_gain_ = true;
  unsigned nchan_ = 0;
  unsigned nsamp_ = 0;
  unsigned window_n_;
  int bkg_window_0_;
  std::vector<int> sig_window_0_;

  std::vector<float> chan_ped_est_;
  float ped_iir_old_;
  float ped_iir_new_;

  std::vector<int> chan_max_index_;
  std::vector<int> chan_max_;
  std::vector<int> chan_bkg_win_sum_;
  std::vector<int> chan_sig_win_sum_;
  std::vector<int> chan_all_sum_q_;
  std::vector<int> chan_all_sum_qt_;

  std::vector<float> all_sig_;
  std::vector<float> all_mean_t_;
};

} } } // namespace calin::iact_data::event_visitor
