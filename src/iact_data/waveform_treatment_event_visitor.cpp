/*

   calin/iact_data/waveform_treatment_event_visitor.cpp -- Stephen Fegan -- 2018-01-11

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

#include <stdexcept>

#include <iact_data/waveform_treatment_event_visitor.hpp>

using namespace calin::ix::iact_data::waveform_treatment_event_visitor;
using namespace calin::iact_data::waveform_treatment_event_visitor;
using namespace calin::ix::iact_data::telescope_run_configuration;
using namespace calin::ix::iact_data::telescope_event;

SingleGainDualWindowWaveformTreatmentEventVisitor::
SingleGainDualWindowWaveformTreatmentEventVisitor(
    calin::ix::iact_data::waveform_treatment_event_visitor::
      SingleGainDualWindowWaveformTreatmentEventVisitorConfig config,
    bool treat_high_gain):
  TelescopeEventVisitor(), config_(config), treat_high_gain_(treat_high_gain)
{
  // nothing to see here
}

SingleGainDualWindowWaveformTreatmentEventVisitor::
~SingleGainDualWindowWaveformTreatmentEventVisitor()
{
  // nothing to see here
}

bool SingleGainDualWindowWaveformTreatmentEventVisitor::demand_waveforms()
{
  return false;
}

bool SingleGainDualWindowWaveformTreatmentEventVisitor::is_parallelizable()
{
  return true;
}

SingleGainDualWindowWaveformTreatmentEventVisitor* SingleGainDualWindowWaveformTreatmentEventVisitor::new_sub_visitor(
  const std::map<TelescopeEventVisitor*,TelescopeEventVisitor*>&
    antecedent_visitors)
{
  return new SingleGainDualWindowWaveformTreatmentEventVisitor(config_);
}

bool SingleGainDualWindowWaveformTreatmentEventVisitor::
visit_telescope_run(const TelescopeRunConfiguration* run_config)
{
  if(config_.integration_n()==0)window_n_ = run_config->num_samples();
  else window_n_ = config_.integration_n();
  if(window_n_ > run_config->num_samples())
    throw std::out_of_range("SingleGainDualWindowWaveformTreatmentEventVisitor: "
      "requested background window larger than number of samples: "
      + std::to_string(window_n_) + ">"
      + std::to_string(run_config->num_samples()));

  if(config_.bkg_integration_0() < -int(run_config->num_samples()-window_n_)
      or config_.bkg_integration_0() > int(run_config->num_samples()-window_n_))
    throw std::out_of_range("SingleGainDualWindowWaveformTreatmentEventVisitor: "
      "requested background window start "
      + std::to_string(config_.bkg_integration_0())
      + " out of allowed range: ["
      + std::to_string(-int(run_config->num_samples()-window_n_)) + ", "
      + std::to_string(run_config->num_samples()-window_n_) + "]");

  if(config_.bkg_integration_0()<0)
    bkg_window_0_ = config_.bkg_integration_0() + run_config->num_samples()
      - window_n_;
  else
    bkg_window_0_ = config_.bkg_integration_0();


  nchan_ = run_config->configured_channel_id_size();
  nsamp_ = run_config->num_samples();

  sig_window_0_.resize(nchan_);
  if(config_.chan_sig_integration_0_size() == 0)
    std::fill(sig_window_0_.begin(), sig_window_0_.end(), config_.sig_integration_0());
  else if(config_.chan_sig_integration_0_size() == nchan_)
    std::copy(config_.chan_sig_integration_0().begin(), config_.chan_sig_integration_0().end(), sig_window_0_.begin());
  else
    throw std::out_of_range("SingleGainDualWindowWaveformTreatmentEventVisitor: "
      "size of chan_sig_integration_0 vector must be either zero or number of configured channels, "
      + std::to_string(config_.pedestal_size()) + " != 0 or "
      + std::to_string(nchan_));

  for(unsigned ichan=0; ichan<nchan_; ichan++)
  {
    if(sig_window_0_[ichan] < -int(run_config->num_samples()-window_n_)
        or sig_window_0_[ichan] > int(run_config->num_samples()-window_n_))
      throw std::out_of_range("SingleGainDualWindowWaveformTreatmentEventVisitor: "
        "requested start of signal window for channel " + std::to_string(ichan) + ", "
        + std::to_string(sig_window_0_[ichan])
        + ", out of allowed range: ["
        + std::to_string(-int(run_config->num_samples()-window_n_)) + ", "
        + std::to_string(run_config->num_samples()-window_n_) + "]");

    if(sig_window_0_[ichan]<0)
      sig_window_0_[ichan] += run_config->num_samples() - window_n_;
  }

  chan_ped_est_.resize(nchan_);
  if(config_.pedestal_size() == 0)
    std::fill(chan_ped_est_.begin(), chan_ped_est_.end(), -1.0f);
  else if(config_.pedestal_size() == nchan_)
    std::copy(config_.pedestal().begin(), config_.pedestal().end(), chan_ped_est_.begin());
  else
    throw std::out_of_range("SingleGainDualWindowWaveformTreatmentEventVisitor: "
      "size of pedestal vector must be either zero or number of configured channels, "
      + std::to_string(config_.pedestal_size()) + " != 0 or "
      + std::to_string(nchan_));

  ped_iir_old_ = std::min(std::max(config_.pedestal_filter_constant(), 0.0f), 1.0f);
  ped_iir_new_ = (1.0f - ped_iir_old_);

  chan_max_index_.resize(nchan_);
  chan_max_.resize(nchan_);
  chan_bkg_win_sum_.resize(nchan_);
  chan_sig_win_sum_.resize(nchan_);
  chan_sig_max_sum_.resize(nchan_);
  chan_sig_max_sum_index_.resize(nchan_);
  chan_all_sum_q_.resize(nchan_);
  chan_all_sum_qt_.resize(nchan_);

  chan_sig_.resize(nchan_);
  chan_mean_t_.resize(nchan_);

  return true;
}

bool SingleGainDualWindowWaveformTreatmentEventVisitor::leave_telescope_run()
{
  return true;
}

bool SingleGainDualWindowWaveformTreatmentEventVisitor::
visit_telescope_event(uint64_t seq_index, TelescopeEvent* event)
{
  const Waveforms* wf = nullptr;
  if(treat_high_gain_) {
    if(event->has_high_gain_image() and
        event->high_gain_image().has_camera_waveforms()) {
      wf = &event->high_gain_image().camera_waveforms();
    }
  } else if(event->has_low_gain_image() and
      event->low_gain_image().has_camera_waveforms()) {
    wf = &event->low_gain_image().camera_waveforms();
  }
  if(wf == nullptr)return true;
  const uint16_t* data = reinterpret_cast<const uint16_t*>(
    wf->raw_samples_array().data() + wf->raw_samples_array_start());
  analyze_waveforms(data, nchan_, nsamp_,
    window_n_, bkg_window_0_, sig_window_0_.data(),
    chan_ped_est_.data(), ped_iir_old_, ped_iir_new_,
    chan_max_index_.data(), chan_max_.data(),
    chan_bkg_win_sum_.data(), chan_sig_win_sum_.data(),
    chan_sig_max_sum_.data(), chan_sig_max_sum_index_.data(),
    chan_all_sum_q_.data(), chan_all_sum_qt_.data(), chan_sig_.data(), chan_mean_t_.data());
  return true;
}

bool SingleGainDualWindowWaveformTreatmentEventVisitor::leave_telescope_event()
{
  return true;
}

void SingleGainDualWindowWaveformTreatmentEventVisitor::
analyze_waveforms(const uint16_t* data, unsigned nchan, unsigned nsamp,
  unsigned window_n, int bkg_window_0, const int* sig_window_0,
  float* ped, float ped_iir_old, float ped_iir_new,
  int* chan_max_index, int* chan_max,
  int* chan_bkg_win_sum, int* chan_sig_win_sum,
  int* chan_sig_max_sum, int* chan_sig_max_sum_index,
  int* chan_all_sum_q, int* chan_all_sum_qt, float* chan_sig, float* chan_mean_t)
{
  int imax = 0;
  int max = 0;
  int bkg = 0;
  int sig = 0;
  int sig_max = 0;
  int isig_max = 0;
  int sum_q = 0;
  int sum_qt = 0;
  int win = 0;
  int isamp = 0;
  int samp[nsamp];

  for(unsigned ichan=0;ichan<nchan;ichan++)
  {
    samp[0] = data[ichan*nsamp];
    imax = 0;
    max = samp[0];
    win = max;
    sum_qt = 0;
    for(isamp = 1;isamp<window_n;isamp++) {
      const unsigned _samp = data[ichan*nsamp+isamp];
      samp[isamp] = _samp;
      win += _samp;
      if(_samp > max) {
        imax = isamp;
        max = _samp;
      }
      sum_qt += _samp*isamp;
    }
    sig_max = win;
    isig_max = 0;
    sum_q = bkg;
    while(isamp<nsamp) {
      unsigned iss = isamp-16;
      if(bkg_window_0 == iss)bkg = win;
      if(sig_window_0[ichan] == iss)sig = win;
      const unsigned _samp = data[ichan*nsamp+isamp];
      samp[isamp] = _samp;
      sum_q += _samp;
      sum_qt += _samp*isamp;
      win += _samp - samp[iss];
      if(win>sig_max) {
        sig_max = win;
        isig_max = iss;
      }
      if(_samp > max) {
        imax = isamp;
        max = _samp;
      }
      ++isamp;
    }
    if(bkg_window_0 == nsamp-window_n)bkg = win;
    if(sig_window_0[ichan] == nsamp-window_n)sig = win;

    chan_max_index[ichan] = imax;
    chan_max[ichan] = max;
    chan_bkg_win_sum[ichan] = bkg;
    chan_sig_win_sum[ichan] = sig;
    chan_sig_max_sum[ichan] = sig_max;
    chan_sig_max_sum_index[ichan] = isig_max;
    chan_all_sum_q[ichan] = sum_q;
    chan_all_sum_qt[ichan] = sum_qt;
    if(ped[ichan] >= 0) {
      ped[ichan] = ped_iir_old*ped[ichan] + ped_iir_new*float(bkg);
    } else {
      ped[ichan] = float(bkg);
    }
    chan_sig[ichan] = float(sig) - ped[ichan];
    chan_mean_t[ichan] =
      (double(sum_qt) - double(bkg*nsamp*(nsamp-1)/2)/double(window_n))/
        (double(sum_q) - double(bkg*nsamp)/double(window_n));
  }
}
