/*

   calin/iact_data/waveform_treatment_event_visitor.cpp -- Stephen Fegan -- 2018-01-11

   Waveform treatment event data visitor - process waveforms

   Copyright 2018, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <stdexcept>

#include <util/memory.hpp>
#include <iact_data/waveform_treatment_event_visitor.hpp>
#include <provenance/system_info.hpp>

using namespace calin::ix::iact_data::waveform_treatment_event_visitor;
using namespace calin::iact_data::waveform_treatment_event_visitor;
using namespace calin::ix::iact_data::telescope_run_configuration;
using namespace calin::ix::iact_data::telescope_event;
using namespace calin::iact_data::event_visitor;
using calin::util::memory::safe_aligned_recalloc;

SingleGainDualWindowWaveformTreatmentEventVisitor::
SingleGainDualWindowWaveformTreatmentEventVisitor(
    calin::ix::iact_data::waveform_treatment_event_visitor::
      SingleGainDualWindowWaveformTreatmentEventVisitorConfig config,
    bool treat_high_gain):
  ParallelEventVisitor(), config_(config), treat_high_gain_(treat_high_gain)
{
  // nothing to see here
}

SingleGainDualWindowWaveformTreatmentEventVisitor::
~SingleGainDualWindowWaveformTreatmentEventVisitor()
{
  free(sig_window_0_);
  free(chan_ped_est_);
  free(chan_max_index_);
  free(chan_max_);
  free(chan_bkg_win_sum_);
  free(chan_sig_win_sum_);
  free(chan_sig_max_sum_);
  free(chan_sig_max_sum_index_);
  free(chan_all_sum_q_);
  free(chan_all_sum_qt_);
  free(chan_sig_);
  free(chan_mean_t_);
}

SingleGainDualWindowWaveformTreatmentEventVisitor* SingleGainDualWindowWaveformTreatmentEventVisitor::new_sub_visitor(
  const std::map<ParallelEventVisitor*,ParallelEventVisitor*>&
    antecedent_visitors)
{
  return new SingleGainDualWindowWaveformTreatmentEventVisitor(config_);
}

bool SingleGainDualWindowWaveformTreatmentEventVisitor::
visit_telescope_run(const TelescopeRunConfiguration* run_config,
  EventLifetimeManager* event_lifetime_manager)
{
  reconfigure(run_config->configured_channel_id_size(), run_config->num_samples());

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

  if(config_.chan_sig_integration_0_size() == 0)
    std::fill(sig_window_0_, sig_window_0_+nchan_, config_.sig_integration_0());
  else if(config_.chan_sig_integration_0_size() == int(nchan_))
    std::copy(config_.chan_sig_integration_0().begin(), config_.chan_sig_integration_0().end(), sig_window_0_);
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

  if(config_.pedestal_size() == 0)
    std::fill(chan_ped_est_, chan_ped_est_+nchan_, -1.0f);
  else if(config_.pedestal_size() == int(nchan_))
    std::copy(config_.pedestal().begin(), config_.pedestal().end(), chan_ped_est_);
  else
    throw std::out_of_range("SingleGainDualWindowWaveformTreatmentEventVisitor: "
      "size of pedestal vector must be either zero or number of configured channels, "
      + std::to_string(config_.pedestal_size()) + " != 0 or "
      + std::to_string(nchan_));

  ped_iir_old_ = std::min(std::max(config_.pedestal_filter_constant(), 0.0f), 1.0f);
  ped_iir_new_ = (1.0f - ped_iir_old_);

  return true;
}

void SingleGainDualWindowWaveformTreatmentEventVisitor::
reconfigure(unsigned nchan, unsigned nsamp)
{
  if(nchan != nchan_){
    auto* host_info = calin::provenance::system_info::the_host_info();
    nchan_ = nchan;
    nsamp_ = nsamp;
    unsigned nalloc = ((nchan_+15)/16)*16;
    safe_aligned_recalloc(sig_window_0_, nalloc, host_info->log2_simd_vec_size());
    safe_aligned_recalloc(chan_ped_est_, nalloc, host_info->log2_simd_vec_size());
    safe_aligned_recalloc(chan_max_index_, nalloc, host_info->log2_simd_vec_size());
    safe_aligned_recalloc(chan_max_, nalloc, host_info->log2_simd_vec_size());
    safe_aligned_recalloc(chan_bkg_win_sum_, nalloc, host_info->log2_simd_vec_size());
    safe_aligned_recalloc(chan_sig_win_sum_, nalloc, host_info->log2_simd_vec_size());
    safe_aligned_recalloc(chan_sig_max_sum_, nalloc, host_info->log2_simd_vec_size());
    safe_aligned_recalloc(chan_sig_max_sum_index_, nalloc, host_info->log2_simd_vec_size());
    safe_aligned_recalloc(chan_all_sum_q_, nalloc, host_info->log2_simd_vec_size());
    safe_aligned_recalloc(chan_all_sum_qt_, nalloc, host_info->log2_simd_vec_size());
    safe_aligned_recalloc(chan_sig_, nalloc, host_info->log2_simd_vec_size());
    safe_aligned_recalloc(chan_mean_t_, nalloc, host_info->log2_simd_vec_size());
  }
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
    wf->raw_samples_array().data());
  scalar_analyze_waveforms(data);
  return true;
}

#if 0
AVX2_SingleGainDualWindowWaveformTreatmentEventVisitor::
~AVX2_SingleGainDualWindowWaveformTreatmentEventVisitor()
{
#if defined(__AVX2__) and defined(__FMA__)
  free(samples_);
  free(q_l_);
  free(q_u_);
  free(qt_l_);
  free(qt_u_);
#else
  // throw std::runtime_error("AVX2_SingleGainDualWindowWaveformTreatmentEventVisitor: AVX2 or FMA not available at compile time");
#endif
}

AVX2_SingleGainDualWindowWaveformTreatmentEventVisitor*
AVX2_SingleGainDualWindowWaveformTreatmentEventVisitor::new_sub_visitor(
  const std::map<ParallelEventVisitor*,ParallelEventVisitor*>&
    antecedent_visitors)
{
  return new AVX2_SingleGainDualWindowWaveformTreatmentEventVisitor(config_);
}

bool AVX2_SingleGainDualWindowWaveformTreatmentEventVisitor::
visit_telescope_run(const TelescopeRunConfiguration* run_config,
  EventLifetimeManager* event_lifetime_manager)
{
#if defined(__AVX2__) and defined(__FMA__)
  bool old_nsamp = nsamp_;
  bool good = SingleGainDualWindowWaveformTreatmentEventVisitor::visit_telescope_run(run_config, event_lifetime_manager);
  if(nsamp_!=old_nsamp) {
    auto* host_info = calin::provenance::system_info::the_host_info();
    const unsigned nv_samp = (nsamp_+15)/16;
    const unsigned nv_block = nv_samp*16;
    safe_aligned_recalloc(samples_, nv_block, host_info->log2_simd_vec_size());
    safe_aligned_recalloc(q_l_, nsamp_, host_info->log2_simd_vec_size());
    safe_aligned_recalloc(q_u_, nsamp_, host_info->log2_simd_vec_size());
    safe_aligned_recalloc(qt_l_, nsamp_, host_info->log2_simd_vec_size());
    safe_aligned_recalloc(qt_u_, nsamp_, host_info->log2_simd_vec_size());

    qt_l_[0] = qt_u_[0] = _mm256_setzero_si256();
  }
  return good;
#else
  throw std::runtime_error("AVX2_SingleGainDualWindowWaveformTreatmentEventVisitor: AVX2 or FMA not available at compile time");
#endif
}

bool AVX2_SingleGainDualWindowWaveformTreatmentEventVisitor::
visit_telescope_event(uint64_t seq_index,
  calin::ix::iact_data::telescope_event::TelescopeEvent* event)
{
#if defined(__AVX2__) and defined(__FMA__)
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
    wf->raw_samples_array().data());
  avx2_analyze_waveforms(data);
  return true;
#else
  throw std::runtime_error("AVX2_SingleGainDualWindowWaveformTreatmentEventVisitor: AVX2 or FMA not available at compile time");
#endif
}
#endif
