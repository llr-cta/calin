/*

   calin/diagnostics/waveform.cpp -- Stephen Fegan -- 2016-04-13

   Waveform diagnostics visitor

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <math/special.hpp>
#include <util/log.hpp>
#include <util/memory.hpp>
#include <diagnostics/waveform.hpp>
#include <math/covariance_calc.hpp>
#include <math/fftw_util.hpp>

using namespace calin::util::log;
using namespace calin::diagnostics::waveform;
using namespace calin::math::fftw_util;
using calin::math::special::SQR;
using calin::math::covariance_calc::cov_i64_gen;
using calin::util::memory::safe_aligned_recalloc_and_fill;

WaveformSumParallelVisitor::WaveformSumParallelVisitor(): ParallelEventVisitor()
{
  // nothing to see here
}

WaveformSumParallelVisitor::~WaveformSumParallelVisitor()
{
  ::free(high_gain_count_i32_);
  ::free(low_gain_count_i32_);
  ::free(high_gain_wf_sum_i32_);
  ::free(low_gain_wf_sum_i32_);
  ::free(high_gain_count_i64_);
  ::free(low_gain_count_i64_);
  ::free(high_gain_wf_sum_i64_);
  ::free(low_gain_wf_sum_i64_);
}

WaveformSumParallelVisitor* WaveformSumParallelVisitor::new_sub_visitor(
  std::map<calin::iact_data::event_visitor::ParallelEventVisitor*,
    calin::iact_data::event_visitor::ParallelEventVisitor*> antecedent_visitors)
{
  auto* sub_visitor = new WaveformSumParallelVisitor();
  sub_visitor->parent_ = this;
  return sub_visitor;
}

bool WaveformSumParallelVisitor::visit_telescope_run(
  const calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration* run_config,
  calin::iact_data::event_visitor::EventLifetimeManager* event_lifetime_manager)
{
  has_dual_gain_ = (run_config->camera_layout().adc_gains() !=
    calin::ix::iact_data::instrument_layout::CameraLayout::SINGLE_GAIN);

  nchan_ = run_config->configured_channel_id_size();
  nsamp_ = run_config->num_samples();

  num_entries_i32_ = 0;

  safe_aligned_recalloc_and_fill(high_gain_count_i32_, nchan_);
  safe_aligned_recalloc_and_fill(high_gain_wf_sum_i32_, nchan_ * nsamp_);
  safe_aligned_recalloc_and_fill(high_gain_count_i64_, nchan_);
  safe_aligned_recalloc_and_fill(high_gain_wf_sum_i64_, nchan_ * nsamp_);

  if(has_dual_gain_) {
    safe_aligned_recalloc_and_fill(low_gain_count_i32_, nchan_);
    safe_aligned_recalloc_and_fill(low_gain_wf_sum_i32_, nchan_ * nsamp_);
    safe_aligned_recalloc_and_fill(low_gain_count_i64_, nchan_);
    safe_aligned_recalloc_and_fill(low_gain_wf_sum_i64_, nchan_ * nsamp_);
  }

  return true;
}

namespace {
  void integrate_i32_to_i64(unsigned nchan, unsigned nsamp,
    uint32_t*__restrict__ count_i32, int32_t*__restrict__ wf_sum_i32,
    uint64_t*__restrict__ count_i64, int64_t*__restrict__ wf_sum_i64)
  {
    for(unsigned ichan = 0; ichan<nchan; ichan++) {
      count_i64[ichan] += count_i32[ichan];
      count_i32[ichan] = 0;
    }

    nsamp *= nchan;
    for(unsigned isamp = 0; isamp<nsamp; isamp++) {
      wf_sum_i64[isamp] += wf_sum_i32[isamp];
      wf_sum_i32[isamp] = 0;
    }
  }

  void integrate_i64_to_parent(unsigned nchan, unsigned nsamp,
    uint64_t*__restrict__ count_i64c, int64_t*__restrict__ wf_sum_i64c,
    uint64_t*__restrict__ count_i64p, int64_t*__restrict__ wf_sum_i64p)
  {
    for(unsigned ichan = 0; ichan<nchan; ichan++) {
      count_i64p[ichan] += count_i64c[ichan];
    }

    nsamp *= nchan;
    for(unsigned isamp = 0; isamp<nsamp; isamp++) {
      wf_sum_i64p[isamp] += wf_sum_i64c[isamp];
    }
  }
}

bool WaveformSumParallelVisitor::leave_telescope_run()
{
  num_entries_i32_ = 0;
  integrate_i32_to_i64(nchan_, nsamp_, high_gain_count_i32_, high_gain_wf_sum_i32_,
    high_gain_count_i64_, high_gain_wf_sum_i64_);
  if(has_dual_gain_) {
    integrate_i32_to_i64(nchan_, nsamp_, low_gain_count_i32_, low_gain_wf_sum_i32_,
      low_gain_count_i64_, low_gain_wf_sum_i64_);
  }
  return true;
}

namespace {
  void perform_single_gain_waveform_sum(unsigned nchan, unsigned nsamp,
    const calin::ix::iact_data::telescope_event::Waveforms& wf,
    uint32_t*__restrict__ count_i32, int32_t*__restrict__ wf_sum_i32)
  {
    for(unsigned ichan = 0; ichan<nchan; ichan++, count_i32++) {
      if(wf.channel_signal_type(ichan) != calin::ix::iact_data::telescope_event::SIGNAL_NONE) {
        ++(*count_i32);
      }
    }
    nsamp *= nchan;
    const int16_t*__restrict__ wf_sum_i16 = reinterpret_cast<const int16_t*>(wf.raw_samples_array().data());
    for(unsigned isamp = 0; isamp<nsamp; isamp++) {
      // Here we use the fact that the waveform data for missing channels has
      // been zeroed out by the decoder
      *(wf_sum_i32++) += *(wf_sum_i16++);
    }
  }

  void perform_mixed_gain_waveform_sum(unsigned nchan, unsigned nsamp,
    const calin::ix::iact_data::telescope_event::Waveforms& wf,
    uint32_t*__restrict__ hg_count_i32, int32_t*__restrict__ hg_wf_sum_i32,
    uint32_t*__restrict__ lg_count_i32, int32_t*__restrict__ lg_wf_sum_i32)
  {
    const int16_t*__restrict__ wf_sum_i16 = reinterpret_cast<const int16_t*>(wf.raw_samples_array().data());
    for(unsigned ichan = 0; ichan<nchan; ichan++, hg_count_i32++, lg_count_i32++) {
      if(wf.channel_signal_type(ichan) != calin::ix::iact_data::telescope_event::SIGNAL_LOW_GAIN) {
        ++(*lg_count_i32);
        for(unsigned isamp = 0; isamp<nsamp; isamp++, hg_wf_sum_i32++) {
          *(lg_wf_sum_i32++) += *(wf_sum_i16++);
        }
      } else if(wf.channel_signal_type(ichan) != calin::ix::iact_data::telescope_event::SIGNAL_NONE) {
        ++(*hg_count_i32);
        for(unsigned isamp = 0; isamp<nsamp; isamp++, lg_wf_sum_i32++) {
          *(hg_wf_sum_i32++) += *(wf_sum_i16++);
        }
      } else {
        wf_sum_i16 += nsamp;
        hg_wf_sum_i32 += nsamp;
        lg_wf_sum_i32 += nsamp;
      }
    }
  }

}

bool WaveformSumParallelVisitor::visit_telescope_event(uint64_t seq_index,
  calin::ix::iact_data::telescope_event::TelescopeEvent* event)
{
  if(!has_dual_gain_) {
    if(event->has_image() and event->image().has_camera_waveforms()) {
      num_entries_i32_++;
      perform_single_gain_waveform_sum(nchan_, nsamp_,
        event->image().camera_waveforms(), high_gain_count_i32_, high_gain_wf_sum_i32_);
    }
  } else {
    if(event->has_image() and event->image().has_camera_waveforms()) {
      num_entries_i32_++;
      perform_mixed_gain_waveform_sum(nchan_, nsamp_,
        event->image().camera_waveforms(),
        high_gain_count_i32_, high_gain_wf_sum_i32_,
        low_gain_count_i32_, low_gain_wf_sum_i32_);
    } else {
      if(event->has_high_gain_image() and event->high_gain_image().has_camera_waveforms()) {
        num_entries_i32_++;
        perform_single_gain_waveform_sum(nchan_, nsamp_,
          event->high_gain_image().camera_waveforms(), high_gain_count_i32_, high_gain_wf_sum_i32_);
        if(event->has_low_gain_image() and event->low_gain_image().has_camera_waveforms()) {
          perform_single_gain_waveform_sum(nchan_, nsamp_,
            event->low_gain_image().camera_waveforms(), low_gain_count_i32_, low_gain_wf_sum_i32_);
        }
      } else if(event->has_low_gain_image() and event->low_gain_image().has_camera_waveforms()) {
        num_entries_i32_++;
        perform_single_gain_waveform_sum(nchan_, nsamp_,
          event->low_gain_image().camera_waveforms(), low_gain_count_i32_, low_gain_wf_sum_i32_);
      }
    }
  }

  if(num_entries_i32_ == partial_max_num_entries_) {
    num_entries_i32_ = 0;
    integrate_i32_to_i64(nchan_, nsamp_, high_gain_count_i32_, high_gain_wf_sum_i32_,
      high_gain_count_i64_, high_gain_wf_sum_i64_);
    if(has_dual_gain_) {
      integrate_i32_to_i64(nchan_, nsamp_, low_gain_count_i32_, low_gain_wf_sum_i32_,
        low_gain_count_i64_, low_gain_wf_sum_i64_);
    }
  }
  return true;
}

bool WaveformSumParallelVisitor::merge_results()
{
  if(parent_) {
    integrate_i64_to_parent(nchan_, nsamp_, high_gain_count_i64_, high_gain_wf_sum_i64_,
      parent_->high_gain_count_i64_, parent_->high_gain_wf_sum_i64_);
    if(has_dual_gain_) {
      integrate_i64_to_parent(nchan_, nsamp_, low_gain_count_i64_, low_gain_wf_sum_i64_,
        parent_->low_gain_count_i64_, parent_->low_gain_wf_sum_i64_);
    }
  }
  return true;
}

calin::ix::diagnostics::waveform::WaveformMean*
WaveformSumParallelVisitor::mean_waveforms() const
{
  auto* results = new calin::ix::diagnostics::waveform::WaveformMean;
  if(high_gain_count_i64_) {
    for(unsigned ichan=0;ichan<nchan_;ichan++) {
      if(high_gain_count_i64_[ichan]) {
        double count = high_gain_count_i64_[ichan];
        auto* wf = results->add_high_gain();
        wf->set_num_entries(high_gain_count_i64_[ichan]);
        for(unsigned isamp=0;isamp<nsamp_;isamp++) {
          wf->add_mean_waveform(
            double(high_gain_wf_sum_i64_[ichan*nsamp_ + isamp])/count);
        }
      }
    }
  }
  if(low_gain_count_i64_) {
    for(unsigned ichan=0;ichan<nchan_;ichan++) {
      if(low_gain_count_i64_[ichan]) {
        double count = low_gain_count_i64_[ichan];
        auto* wf = results->add_low_gain();
        wf->set_num_entries(low_gain_count_i64_[ichan]);
        for(unsigned isamp=0;isamp<nsamp_;isamp++) {
          wf->add_mean_waveform(
            double(low_gain_wf_sum_i64_[ichan*nsamp_ + isamp])/count);
        }
      }
    }
  }
  return results;
}

void WaveformSumParallelVisitor::high_gain_mean_wf(Eigen::MatrixXd& wf_matrix_out)
{
  wf_matrix_out.resize(nchan_, nsamp_);
  for(unsigned ichan=0;ichan<nchan_;ichan++) {
    double count = high_gain_count_i64_[ichan];
    for(unsigned isamp=0;isamp<nsamp_;isamp++) {
      wf_matrix_out(ichan, isamp) =
        double(high_gain_wf_sum_i64_[ichan*nsamp_ + isamp])/count;
    }
  }
}

void WaveformSumParallelVisitor::low_gain_mean_wf(Eigen::MatrixXd& wf_matrix_out)
{
  wf_matrix_out.resize(nchan_, nsamp_);
  for(unsigned ichan=0;ichan<nchan_;ichan++) {
    double count = low_gain_count_i64_[ichan];
    for(unsigned isamp=0;isamp<nsamp_;isamp++) {
      wf_matrix_out(ichan, isamp) =
        double(low_gain_wf_sum_i64_[ichan*nsamp_ + isamp])/count;
    }
  }
}

void WaveformSumParallelVisitor::high_gain_event_count(Eigen::VectorXi& count_out)
{
  count_out.resize(nchan_);
  for(unsigned ichan=0;ichan<nchan_;ichan++) {
    count_out(ichan) = high_gain_count_i64_[ichan];
  }
}

void WaveformSumParallelVisitor::low_gain_event_count(Eigen::VectorXi& count_out)
{
  count_out.resize(nchan_);
  for(unsigned ichan=0;ichan<nchan_;ichan++) {
    count_out(ichan) = low_gain_count_i64_[ichan];
  }
}



WaveformStatsParallelVisitor::WaveformStatsParallelVisitor(
    bool calculate_psd, bool calculate_covariance):
  ParallelEventVisitor(),
  calculate_psd_(calculate_psd), calculate_covariance_(calculate_covariance)
{
  // nothing to see here
}

WaveformStatsParallelVisitor::~WaveformStatsParallelVisitor()
{
  // nothing to see here
}

WaveformStatsParallelVisitor* WaveformStatsParallelVisitor::new_sub_visitor(
  std::map<calin::iact_data::event_visitor::ParallelEventVisitor*,
      calin::iact_data::event_visitor::ParallelEventVisitor*>
    antecedent_visitors)
{
  auto* sub_visitor = new WaveformStatsParallelVisitor(calculate_covariance_);
  sub_visitor->parent_ = this;
  return sub_visitor;
}

bool WaveformStatsParallelVisitor::visit_telescope_run(
  const calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration* run_config,
  calin::iact_data::event_visitor::EventLifetimeManager* event_lifetime_manager)
{
  run_config_ = run_config;
  results_.Clear();
  unsigned N = run_config->num_samples();
  for(int ichan = 0; ichan<run_config->configured_channel_id_size(); ichan++)
  {
    auto* hg_wf = results_.add_high_gain();
    hg_wf->mutable_sum()->Resize(N,0);
    hg_wf->mutable_sum_squared()->Resize(N,0);
    if(calculate_covariance_)
      hg_wf->mutable_sum_product()->Resize(N*(N-1)/2,0);
    auto* lg_wf = results_.add_low_gain();
    lg_wf->mutable_sum()->Resize(N,0);
    lg_wf->mutable_sum_squared()->Resize(N,0);
    if(calculate_covariance_)
      lg_wf->mutable_sum_product()->Resize(N*(N-1)/2,0);

    auto* phg_wf = partial_.add_high_gain();
    phg_wf->mutable_sum()->Resize(N,0);
    phg_wf->mutable_sum_squared()->Resize(N,0);
    if(calculate_covariance_)
      phg_wf->mutable_sum_product()->Resize(N*(N-1)/2,0);
    auto* plg_wf = partial_.add_low_gain();
    plg_wf->mutable_sum()->Resize(N,0);
    plg_wf->mutable_sum_squared()->Resize(N,0);
    if(calculate_covariance_)
      plg_wf->mutable_sum_product()->Resize(N*(N-1)/2,0);
  }

  if(calculate_psd_)
  {
    psd_results_.Clear();
    unsigned nfreq = 1 + N/2; // 1->1 2->2 3->2 4->3 5->3 etc
    free(waveform_t_);
    free(waveform_f_);
    waveform_t_ = fftwf_alloc_real(N);
    waveform_f_ = fftwf_alloc_real(N);
    assert(waveform_t_);
    assert(waveform_f_);
    fftw_plan_fwd_ =
      fftwf_plan_r2r_1d(N, waveform_t_, waveform_f_, FFTW_R2HC, FFTW_MEASURE|FFTW_DESTROY_INPUT);
    assert(fftw_plan_fwd_);
    fftw_plan_bwd_ =
      fftwf_plan_r2r_1d(N, waveform_f_, waveform_t_, FFTW_HC2R, FFTW_MEASURE|FFTW_DESTROY_INPUT);
    assert(fftw_plan_bwd_);
    for(int ichan = 0; ichan<run_config->configured_channel_id_size(); ichan++)
    {
      auto* hg_wf = psd_results_.add_high_gain();
      hg_wf->mutable_psd_sum()->Resize(nfreq,0);
      hg_wf->mutable_psd_sum_squared()->Resize(nfreq,0);
      hg_wf->mutable_corr_sum()->Resize(N,0);
      hg_wf->mutable_corr_sum_squared()->Resize(N,0);
      auto* lg_wf = psd_results_.add_low_gain();
      lg_wf->mutable_psd_sum()->Resize(nfreq,0);
      lg_wf->mutable_psd_sum_squared()->Resize(nfreq,0);
      lg_wf->mutable_corr_sum()->Resize(N,0);
      lg_wf->mutable_corr_sum_squared()->Resize(N,0);
    }
  }

  return true;
}

bool WaveformStatsParallelVisitor::leave_telescope_run()
{
  for(int ichan=0; ichan<results_.high_gain_size(); ichan++)
    merge_partial(partial_.mutable_high_gain(ichan),
      results_.mutable_high_gain(ichan));
  for(int ichan=0; ichan<results_.low_gain_size(); ichan++)
    merge_partial(partial_.mutable_low_gain(ichan),
      results_.mutable_low_gain(ichan));
  run_config_ = nullptr;
  if(calculate_psd_)
  {
    fftwf_destroy_plan(fftw_plan_bwd_);
    fftw_plan_bwd_ = nullptr;
    fftwf_destroy_plan(fftw_plan_fwd_);
    fftw_plan_fwd_ = nullptr;
    fftwf_free(waveform_t_);
    waveform_t_ = nullptr;
    fftwf_free(waveform_f_);
    waveform_f_ = nullptr;
  }
  return true;
}

bool WaveformStatsParallelVisitor::visit_telescope_event(uint64_t seq_index,
  calin::ix::iact_data::telescope_event::TelescopeEvent* event)
{
  const int nchan = run_config_->configured_channel_id_size();
  const int nsamp = run_config_->num_samples();

  if(event->has_high_gain_image() and event->high_gain_image().has_camera_waveforms())
  {
    const ix::iact_data::telescope_event::Waveforms* wf =
      &event->high_gain_image().camera_waveforms();
    const uint16_t*__restrict__ wf_data = reinterpret_cast<const uint16_t*__restrict__>(
        wf->raw_samples_array().data());
    for(int ichan = 0; ichan<nchan; ichan++) {
      calin::ix::diagnostics::waveform::WaveformRawPSD* psd = nullptr;
      if(calculate_psd_)psd = psd_results_.mutable_high_gain(ichan);
      process_one_waveform(wf_data, partial_.mutable_high_gain(ichan),
        results_.mutable_high_gain(ichan), psd);
      wf_data += nsamp;
    }
  }

  if(event->has_low_gain_image() and event->low_gain_image().has_camera_waveforms())
  {
    const ix::iact_data::telescope_event::Waveforms* wf =
      &event->low_gain_image().camera_waveforms();
    const uint16_t*__restrict__ wf_data = reinterpret_cast<const uint16_t*__restrict__>(
        wf->raw_samples_array().data());
    for(int ichan = 0; ichan<nchan; ichan++) {
      calin::ix::diagnostics::waveform::WaveformRawPSD* psd = nullptr;
      if(calculate_psd_)psd = psd_results_.mutable_low_gain(ichan);
      process_one_waveform(wf_data, partial_.mutable_low_gain(ichan),
        results_.mutable_low_gain(ichan), psd);
      wf_data += nsamp;
    }
  }

  // nothing to see here
  return true;
}

void WaveformStatsParallelVisitor::
process_one_waveform(const uint16_t*__restrict__ wf,
  ix::diagnostics::waveform::PartialWaveformRawStats* p_stat,
  ix::diagnostics::waveform::WaveformRawStats* r_stat,
  calin::ix::diagnostics::waveform::WaveformRawPSD* psd)
{
  const unsigned nsample = run_config_->num_samples();
  p_stat->set_num_entries(p_stat->num_entries()+1);
  auto*__restrict__ sum = p_stat->mutable_sum()->mutable_data();
  auto*__restrict__ sum_squared = p_stat->mutable_sum_squared()->mutable_data();
  if(psd) {
    for(unsigned isample=0; isample<nsample; isample++) {
      uint32_t sample = wf[isample];
      waveform_t_[isample] = float(sample);
      sum[isample] += sample;
      sum_squared[isample] += sample*sample;
    }
    fftwf_execute(fftw_plan_fwd_);
    psd->set_num_entries(psd->num_entries()+1);
    auto* psd_sum = psd->mutable_psd_sum()->mutable_data();
    auto* psd_sum_squared = psd->mutable_psd_sum_squared()->mutable_data();
    const float*__restrict__ ri = waveform_f_;
    const float*__restrict__ ci = waveform_f_ + nsample-1;
    double psdi = SQR(*ri++);
    (*psd_sum++) += psdi;
    (*psd_sum_squared++) += SQR(psdi);
    while(ri < ci)
    {
      psdi = SQR(*ri++) + SQR(*ci--);
      (*psd_sum++) += psdi;
      (*psd_sum_squared++) += SQR(psdi);
    }
    if(ri==ci)
    {
      double psdi = SQR(*ri);
      *psd_sum += psdi;
      *psd_sum_squared += SQR(psdi);
    }
    hcvec_scale_and_multiply_conj(waveform_f_, waveform_f_, waveform_f_, nsample);
    fftwf_execute(fftw_plan_bwd_);
    auto* corr_sum = psd->mutable_corr_sum()->mutable_data();
    for(unsigned isample=0; isample<nsample; isample++)
      corr_sum[isample] += waveform_t_[isample];
    auto* corr_sum_squared = psd->mutable_corr_sum_squared()->mutable_data();
    for(unsigned isample=0; isample<nsample; isample++)
      corr_sum_squared[isample] += SQR(waveform_t_[isample]);
  } else {
    for(unsigned isample=0; isample<nsample; isample++) {
      uint32_t sample = wf[isample];
      sum[isample] += sample;
      sum_squared[isample] += sample*sample;
    }
  }

  if(calculate_covariance_)
  {
    auto*__restrict__ sum_product = p_stat->mutable_sum_product()->mutable_data();
    const uint16_t*__restrict__ wf_j = wf;
    unsigned msample = nsample;
    for(unsigned isample=0; isample<nsample; isample++)
    {
      uint32_t sample_i = *wf_j;
      ++wf_j;
      --msample;
      for(unsigned jsample=0; jsample<msample; jsample++)
        sum_product[jsample] += sample_i * uint32_t(wf_j[jsample]);
      sum_product += msample;
    }
  }
  if(p_stat->num_entries() == partial_max_num_entries_)
    merge_partial(p_stat, r_stat);
}

bool WaveformStatsParallelVisitor::merge_results()
{
  if(parent_) {
    parent_->results_.IntegrateFrom(results_);
    if(calculate_psd_)parent_->psd_results_.IntegrateFrom(psd_results_);
  }
  return true;
}

namespace {

template<typename T1, typename T2>
void transfer_partial_array(T1* to, T2* from)
{
  int n = from->size();
  assert(to->size() == n);
  auto* to_data = to->mutable_data();
  auto* from_data = from->mutable_data();
  for(int i=0; i<n; i++)to_data[i] += from_data[i];
  for(int i=0; i<n; i++)from_data[i] = 0;
}

} // anonymous namespace

void WaveformStatsParallelVisitor::merge_partial(
  ix::diagnostics::waveform::PartialWaveformRawStats* p_stat,
  ix::diagnostics::waveform::WaveformRawStats* r_stat)
{
  r_stat->set_num_entries(r_stat->num_entries() + p_stat->num_entries());
  p_stat->set_num_entries(0);
  transfer_partial_array(r_stat->mutable_sum(),
    p_stat->mutable_sum());
  transfer_partial_array(r_stat->mutable_sum_squared(),
    p_stat->mutable_sum_squared());
  transfer_partial_array(r_stat->mutable_sum_product(),
    p_stat->mutable_sum_product());
}

Eigen::VectorXd WaveformStatsParallelVisitor::waveform_mean(
  const ix::diagnostics::waveform::WaveformRawStats* stat)
{
  const int N = stat->sum_size();
  Eigen::VectorXd m(N);
  const double one_over_n = 1.0/double(stat->num_entries());
  for(int i=0; i<N; i++)
    m(i) = double(stat->sum(i)) * one_over_n;
  return m;
}

Eigen::VectorXd WaveformStatsParallelVisitor::waveform_var(
  const ix::diagnostics::waveform::WaveformRawStats* stat)
{
  const int N = stat->sum_size();
  Eigen::VectorXd v(N);
  for(int i=0; i<N; i++)
    v(i) = cov_i64_gen(stat->sum_squared(i), stat->num_entries(),
      stat->sum(i), stat->num_entries(), stat->sum(i), stat->num_entries());
  return v;
}

Eigen::MatrixXd WaveformStatsParallelVisitor::waveform_cov(
  const ix::diagnostics::waveform::WaveformRawStats* stat)
{
  const int N = stat->sum_size();
  Eigen::MatrixXd c(N,N);
  for(int i=0; i<N; i++)
    c(i,i) = cov_i64_gen(stat->sum_squared(i), stat->num_entries(),
      stat->sum(i), stat->num_entries(), stat->sum(i), stat->num_entries());
  for(int i=0; i<N; i++)
    for(int j=i+1; j<N; j++)
    {
      double cij =
        cov_i64_gen(stat->sum_product(N*(N-1)/2-(N-i)*(N-i-1)/2 + j - i - 1),
          stat->num_entries(),
          stat->sum(i), stat->num_entries(),
          stat->sum(j), stat->num_entries());
      c(i,j) = cij;
      c(j,i) = cij;
    }
  return c;
}

Eigen::MatrixXd WaveformStatsParallelVisitor::waveform_cov_frac(
  const ix::diagnostics::waveform::WaveformRawStats* stat)
{
  Eigen::MatrixXd c = waveform_cov(stat);
  const int N = stat->sum_size();
  for(int i=0; i<N; i++) {
    double scale = 1.0/std::sqrt(c(i,i));
    for(int j=0; j<N; j++){
      c(i,j) *= scale;
      c(j,i) *= scale;
    }
  }
  return c;
}

// *****************************************************************************
// *****************************************************************************
// *****************************************************************************

// OBSOLETE version that uses TelescopeEventVisitor interface

// *****************************************************************************
// *****************************************************************************
// *****************************************************************************


WaveformStatsVisitor::WaveformStatsVisitor(bool calculate_covariance):
  TelescopeEventVisitor(), calculate_covariance_(calculate_covariance)
{
  // nothing to see here
}

WaveformStatsVisitor::~WaveformStatsVisitor()
{
  // nothing to see here
}

bool WaveformStatsVisitor::demand_waveforms()
{
  return true;
}

bool WaveformStatsVisitor::is_parallelizable()
{
  return true;
}

WaveformStatsVisitor* WaveformStatsVisitor::new_sub_visitor(
  const std::map<calin::iact_data::event_visitor::TelescopeEventVisitor*,
      calin::iact_data::event_visitor::TelescopeEventVisitor*>&
    antecedent_visitors)
{
  auto* sub_visitor = new WaveformStatsVisitor;
  sub_visitor->parent_ = this;
  return sub_visitor;
}

bool WaveformStatsVisitor::visit_telescope_run(
  const calin::ix::iact_data::telescope_run_configuration::
    TelescopeRunConfiguration* run_config)
{
  run_config_ = run_config;
  results_.Clear();
  unsigned N = run_config->num_samples();
  for(int ichan = 0; ichan<run_config->configured_channel_id_size(); ichan++)
  {
    auto* hg_wf = results_.add_high_gain();
    hg_wf->mutable_sum()->Resize(N,0);
    hg_wf->mutable_sum_squared()->Resize(N,0);
    if(calculate_covariance_)
      hg_wf->mutable_sum_product()->Resize(N*(N-1)/2,0);
    auto* lg_wf = results_.add_low_gain();
    lg_wf->mutable_sum()->Resize(N,0);
    lg_wf->mutable_sum_squared()->Resize(N,0);
    if(calculate_covariance_)
      lg_wf->mutable_sum_product()->Resize(N*(N-1)/2,0);

    auto* phg_wf = partial_.add_high_gain();
    phg_wf->mutable_sum()->Resize(N,0);
    phg_wf->mutable_sum_squared()->Resize(N,0);
    if(calculate_covariance_)
      phg_wf->mutable_sum_product()->Resize(N*(N-1)/2,0);
    auto* plg_wf = partial_.add_low_gain();
    plg_wf->mutable_sum()->Resize(N,0);
    plg_wf->mutable_sum_squared()->Resize(N,0);
    if(calculate_covariance_)
      plg_wf->mutable_sum_product()->Resize(N*(N-1)/2,0);
  }
  return true;
}

bool WaveformStatsVisitor::leave_telescope_run()
{
  for(int ichan=0; ichan<results_.high_gain_size(); ichan++)
    merge_partial(partial_.mutable_high_gain(ichan),
      results_.mutable_high_gain(ichan));
  for(int ichan=0; ichan<results_.low_gain_size(); ichan++)
    merge_partial(partial_.mutable_low_gain(ichan),
      results_.mutable_low_gain(ichan));
  run_config_ = nullptr;
  return true;
}

bool WaveformStatsVisitor::visit_telescope_event(uint64_t seq_index,
  calin::ix::iact_data::telescope_event::TelescopeEvent* event)
{
  // nothing to see here
  return true;
}

bool WaveformStatsVisitor::visit_waveform(unsigned ichan,
  calin::ix::iact_data::telescope_event::ChannelWaveform* high_gain,
  calin::ix::iact_data::telescope_event::ChannelWaveform* low_gain)
{
  const int index = ichan; //run_config_->configured_channel_index(ichan);
  if(high_gain)
    process_one_waveform(high_gain, partial_.mutable_high_gain(index),
      results_.mutable_high_gain(index));
  if(low_gain)
    process_one_waveform(low_gain, partial_.mutable_low_gain(index),
      results_.mutable_low_gain(index));
  return true;
}

void WaveformStatsVisitor::
process_one_waveform(
  const calin::ix::iact_data::telescope_event::ChannelWaveform* wf,
  ix::diagnostics::waveform::PartialWaveformRawStats* p_stat,
  ix::diagnostics::waveform::WaveformRawStats* r_stat)
{
  const unsigned nsample = run_config_->num_samples();
  assert(wf->samples_size() == int(nsample));
  const auto* sample = wf->samples().data();
  p_stat->set_num_entries(p_stat->num_entries()+1);
  auto* sum = p_stat->mutable_sum()->mutable_data();
  for(unsigned isample=0; isample<nsample; isample++)
    sum[isample] += sample[isample];
  auto* sum_squared = p_stat->mutable_sum_squared()->mutable_data();
  for(unsigned isample=0; isample<nsample; isample++)
    sum_squared[isample] += sample[isample] * sample[isample];
  if(calculate_covariance_)
  {
    auto* sum_product = p_stat->mutable_sum_product()->mutable_data();
    const auto* sample_j = sample;
    unsigned msample = nsample;
    for(unsigned isample=0; isample<nsample; isample++)
    {
      uint32_t sample_i = *sample_j;
      ++sample_j;
      --msample;
      for(unsigned jsample=0; jsample<msample; jsample++)
        sum_product[jsample] += sample_i * sample_j[jsample];
      sum_product += msample;
    }
  }
  if(p_stat->num_entries() == partial_max_num_entries_)
    merge_partial(p_stat, r_stat);
}

bool WaveformStatsVisitor::merge_results()
{
  if(parent_)parent_->results_.IntegrateFrom(results_);
  return true;
}

void WaveformStatsVisitor::merge_partial(
  ix::diagnostics::waveform::PartialWaveformRawStats* p_stat,
  ix::diagnostics::waveform::WaveformRawStats* r_stat)
{
  r_stat->set_num_entries(r_stat->num_entries() + p_stat->num_entries());
  p_stat->set_num_entries(0);
  transfer_partial_array(r_stat->mutable_sum(),
    p_stat->mutable_sum());
  transfer_partial_array(r_stat->mutable_sum_squared(),
    p_stat->mutable_sum_squared());
  transfer_partial_array(r_stat->mutable_sum_product(),
    p_stat->mutable_sum_product());
}

Eigen::VectorXd WaveformStatsVisitor::waveform_mean(
  const ix::diagnostics::waveform::WaveformRawStats* stat)
{
  const int N = stat->sum_size();
  Eigen::VectorXd m(N);
  const double one_over_n = 1.0/double(stat->num_entries());
  for(int i=0; i<N; i++)
    m(i) = double(stat->sum(i)) * one_over_n;
  return m;
}

Eigen::VectorXd WaveformStatsVisitor::waveform_var(
  const ix::diagnostics::waveform::WaveformRawStats* stat)
{
  const int N = stat->sum_size();
  Eigen::VectorXd v(N);
  for(int i=0; i<N; i++)
    v(i) = cov_i64_gen(stat->sum_squared(i), stat->num_entries(),
      stat->sum(i), stat->num_entries(), stat->sum(i), stat->num_entries());
  return v;
}

Eigen::MatrixXd WaveformStatsVisitor::waveform_cov(
  const ix::diagnostics::waveform::WaveformRawStats* stat)
{
  const int N = stat->sum_size();
  Eigen::MatrixXd c(N,N);
  for(int i=0; i<N; i++)
    c(i,i) = cov_i64_gen(stat->sum_squared(i), stat->num_entries(),
      stat->sum(i), stat->num_entries(), stat->sum(i), stat->num_entries());
  for(int i=0; i<N; i++)
    for(int j=i+1; j<N; j++)
    {
      double cij =
        cov_i64_gen(stat->sum_product(N*(N-1)/2-(N-i)*(N-i-1)/2 + j - i - 1),
          stat->num_entries(),
          stat->sum(i), stat->num_entries(),
          stat->sum(j), stat->num_entries());
      c(i,j) = cij;
      c(j,i) = cij;
    }
  return c;
}

Eigen::MatrixXd WaveformStatsVisitor::waveform_cov_frac(
  const ix::diagnostics::waveform::WaveformRawStats* stat)
{
  Eigen::MatrixXd c = waveform_cov(stat);
  const int N = stat->sum_size();
  for(int i=0; i<N; i++) {
    double scale = 1.0/std::sqrt(c(i,i));
    for(int j=0; j<N; j++){
      c(i,j) *= scale;
      c(j,i) *= scale;
    }
  }
  return c;
}
