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

#include <util/memory.hpp>
#include <iact_data/waveform_treatment_event_visitor.hpp>
#include <provenance/system_info.hpp>
#include <math/simd.hpp>

using namespace calin::ix::iact_data::waveform_treatment_event_visitor;
using namespace calin::iact_data::waveform_treatment_event_visitor;
using namespace calin::ix::iact_data::telescope_run_configuration;
using namespace calin::ix::iact_data::telescope_event;
using calin::util::memory::safe_aligned_recalloc;

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
    safe_aligned_recalloc(sig_window_0_, nchan_, host_info->log2_simd_vec_size());
    safe_aligned_recalloc(chan_ped_est_, nchan_, host_info->log2_simd_vec_size());
    safe_aligned_recalloc(chan_max_index_, nchan_, host_info->log2_simd_vec_size());
    safe_aligned_recalloc(chan_max_, nchan_, host_info->log2_simd_vec_size());
    safe_aligned_recalloc(chan_bkg_win_sum_, nchan_, host_info->log2_simd_vec_size());
    safe_aligned_recalloc(chan_sig_win_sum_, nchan_, host_info->log2_simd_vec_size());
    safe_aligned_recalloc(chan_sig_max_sum_, nchan_, host_info->log2_simd_vec_size());
    safe_aligned_recalloc(chan_sig_max_sum_index_, nchan_, host_info->log2_simd_vec_size());
    safe_aligned_recalloc(chan_all_sum_q_, nchan_, host_info->log2_simd_vec_size());
    safe_aligned_recalloc(chan_all_sum_qt_, nchan_, host_info->log2_simd_vec_size());
    safe_aligned_recalloc(chan_sig_, nchan_, host_info->log2_simd_vec_size());
    safe_aligned_recalloc(chan_mean_t_, nchan_, host_info->log2_simd_vec_size());
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
    wf->raw_samples_array().data() + wf->raw_samples_array_start());
  scalar_analyze_waveforms(data);
  return true;
}

AVX2_SingleGainDualWindowWaveformTreatmentEventVisitor::
~AVX2_SingleGainDualWindowWaveformTreatmentEventVisitor()
{
#if defined(__AVX2__) and defined(__FMA__)
  delete[] samples_;
  delete[] q_l_;
  delete[] q_u_;
  delete[] qt_l_;
  delete[] qt_u_;
#else
  throw std::runtime_error("AVX2_SingleGainDualWindowWaveformTreatmentEventVisitor: AVX2 or FMA not available at compile time");
#endif
}

AVX2_SingleGainDualWindowWaveformTreatmentEventVisitor*
AVX2_SingleGainDualWindowWaveformTreatmentEventVisitor::new_sub_visitor(
  const std::map<TelescopeEventVisitor*,TelescopeEventVisitor*>&
    antecedent_visitors)
{
  return new AVX2_SingleGainDualWindowWaveformTreatmentEventVisitor(config_);
}

bool AVX2_SingleGainDualWindowWaveformTreatmentEventVisitor::
visit_telescope_run(const TelescopeRunConfiguration* run_config)
{
  bool old_nsamp = nsamp_;
  bool good = SingleGainDualWindowWaveformTreatmentEventVisitor::visit_telescope_run(run_config);
  if(nsamp_!=old_nsamp) {
    delete[] samples_;
    delete[] q_l_;
    delete[] q_u_;
    delete[] qt_l_;
    delete[] qt_u_;
    const unsigned nv_samp = (nsamp_+15)/16;
    const unsigned nv_block = nv_samp*16;
    samples_ = new __m256i[nv_block];
    q_l_ = new __m256i[nsamp_];
    q_u_ = new __m256i[nsamp_];
    qt_l_ = new __m256i[nsamp_];
    qt_u_ = new __m256i[nsamp_];
  }
  return good;
}

bool AVX2_SingleGainDualWindowWaveformTreatmentEventVisitor::
visit_telescope_event(uint64_t seq_index,
  calin::ix::iact_data::telescope_event::TelescopeEvent* event)
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
  avx2_analyze_waveforms(data);
  return true;
}
