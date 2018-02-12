/*

   calin/diagnostics/waveform.cpp -- Stephen Fegan -- 2016-04-13

   Waveform diagnostics visitor

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <math/special.hpp>
#include <util/log.hpp>
#include <diagnostics/waveform.hpp>
#include <math/covariance_calc.hpp>
#include <util/memory.hpp>

using namespace calin::util::log;
using namespace calin::diagnostics::waveform;
using calin::util::memory::safe_aligned_recalloc;

AVX2_Unroll8_WaveformStatsVisitor::
AVX2_Unroll8_WaveformStatsVisitor(bool high_gain, bool calculate_covariance):
  ParallelEventVisitor(),
  high_gain_(high_gain), calculate_covariance_(calculate_covariance)
{
#if defined(__AVX2__)
  // nothing to see here
#else // defined(__AVX2__)
  throw std::runtime_error("AVX2_Unroll8_WaveformStatsVisitor: AVX2 not supported at compile time");
#endif // defined(__AVX2__)
}

AVX2_Unroll8_WaveformStatsVisitor::~AVX2_Unroll8_WaveformStatsVisitor()
{
  for(unsigned i=0;i<8;i++)free(samples_[i]);
  free(partial_chan_nevent_);
  free(partial_chan_sum_);
  free(partial_chan_sum_squared_);
  if(calculate_covariance_)free(partial_chan_sum_cov_);
}

AVX2_Unroll8_WaveformStatsVisitor* AVX2_Unroll8_WaveformStatsVisitor::
new_sub_visitor(
  const std::map<calin::iact_data::event_visitor::ParallelEventVisitor*,
      calin::iact_data::event_visitor::ParallelEventVisitor*>&
    antecedent_visitors)
{
#if defined(__AVX2__)
  AVX2_Unroll8_WaveformStatsVisitor* child =
    new AVX2_Unroll8_WaveformStatsVisitor(high_gain_, calculate_covariance_);
  child->parent_ = this;
  return child;
#else // defined(__AVX2__)
  throw std::runtime_error("AVX2_Unroll8_WaveformStatsVisitor: AVX2 not supported at compile time");
#endif // defined(__AVX2__)
}

bool AVX2_Unroll8_WaveformStatsVisitor::visit_telescope_run(
  const calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration* run_config,
  calin::iact_data::event_visitor::EventLifetimeManager* event_lifetime_manager)
{
#if defined(__AVX2__)
  unsigned old_nsamp = nsamp_;
  unsigned old_nchan = nchan_;

  nsamp_ = run_config->num_samples();
  nchan_ = run_config->configured_channel_id_size();

  results_.Clear();
  for(unsigned ichan = 0; ichan<nchan_; ichan++)
  {
    if(high_gain_) {
      auto* hg_wf = results_.add_high_gain();
      hg_wf->mutable_sum()->Resize(nsamp_,0);
      hg_wf->mutable_sum_squared()->Resize(nsamp_,0);
      if(calculate_covariance_)
        hg_wf->mutable_sum_product()->Resize(nsamp_*(nsamp_-1)/2,0);
    } else {
      auto* lg_wf = results_.add_low_gain();
      lg_wf->mutable_sum()->Resize(nsamp_,0);
      lg_wf->mutable_sum_squared()->Resize(nsamp_,0);
      if(calculate_covariance_)
        lg_wf->mutable_sum_product()->Resize(nsamp_*(nsamp_-1)/2,0);
    }
  }

  unsigned nsamp_block = (nsamp_ + 15)/16; // number of uint16_t vectors in nsamp
  unsigned nchan_block = (nchan_ + 15)/16; // number of uint16_t vectors in nchan

  nkept_events_ = 0;

  if(nsamp_ != old_nsamp or samples_[0] == nullptr) {
    for(unsigned i=0;i<8;i++)safe_aligned_recalloc(samples_[i], nsamp_block);
  }

  partial_num_entries_ = 0;
  if(nsamp_ != old_nchan or partial_chan_nevent_ == nullptr) {
    safe_aligned_recalloc(partial_chan_nevent_, nchan_);
    safe_aligned_recalloc(partial_chan_sum_, nchan_block*nsamp_);
    safe_aligned_recalloc(partial_chan_sum_squared_, nchan_block*nsamp_);
    if(calculate_covariance_)
      safe_aligned_recalloc(partial_chan_sum_cov_, nchan_block*nsamp_*(nsamp_-1)/2);
  }

  std::fill(partial_chan_nevent_, partial_chan_nevent_+nchan_, 0);
  std::fill(partial_chan_sum_, partial_chan_sum_+nchan_block*nsamp_, _mm256_setzero_si256());
  std::fill(partial_chan_sum_squared_, partial_chan_sum_squared_+nchan_block*nsamp_, _mm256_setzero_si256());
  std::fill(partial_chan_sum_cov_, partial_chan_sum_cov_+nchan_block*nsamp_*(nsamp_-1)/2, _mm256_setzero_si256());

  event_lifetime_manager_ = event_lifetime_manager;
  return true;
#else // defined(__AVX2__)
  throw std::runtime_error("AVX2_Unroll8_WaveformStatsVisitor: AVX2 not supported at compile time");
#endif // defined(__AVX2__)
}

bool AVX2_Unroll8_WaveformStatsVisitor::leave_telescope_run()
{
#if defined(__AVX2__)
  process_8_events();
  merge_partials();
  event_lifetime_manager_ = nullptr;
  return true;
#else // defined(__AVX2__)
  throw std::runtime_error("AVX2_Unroll8_WaveformStatsVisitor: AVX2 not supported at compile time");
#endif // defined(__AVX2__)
}

bool AVX2_Unroll8_WaveformStatsVisitor::visit_telescope_event(uint64_t seq_index,
  calin::ix::iact_data::telescope_event::TelescopeEvent* event)
{
#if defined(__AVX2__)
  bool keep_this_event = false;
  if(high_gain_) {
    if(event->has_high_gain_image() and
        event->high_gain_image().has_camera_waveforms()) {
      keep_this_event = true;
    }
  } else if(event->has_low_gain_image() and
      event->low_gain_image().has_camera_waveforms()) {
    keep_this_event = true;
  }

  if(keep_this_event) {
    event_lifetime_manager_->keep_event(event);
    kept_events_[nkept_events_++] = event;
    if(nkept_events_ == 8)process_8_events();
    if(partial_num_entries_ == partial_max_num_entries_)merge_partials();
  }
  return true;
#else // defined(__AVX2__)
  throw std::runtime_error("AVX2_Unroll8_WaveformStatsVisitor: AVX2 not supported at compile time");
#endif // defined(__AVX2__)
}

bool AVX2_Unroll8_WaveformStatsVisitor::merge_results()
{
#if defined(__AVX2__)
  const unsigned nsamp_block = (nsamp_ + 15)/16; // number of uint16_t vectors in nsamp
  const unsigned nchan_block = (nchan_ + 15)/16; // number of uint16_t vectors in nchan

  std::fill(partial_chan_nevent_, partial_chan_nevent_+nchan_, 0);
  std::fill(partial_chan_sum_, partial_chan_sum_+nchan_block*nsamp_, _mm256_setzero_si256());
  std::fill(partial_chan_sum_squared_, partial_chan_sum_squared_+nchan_block*nsamp_, _mm256_setzero_si256());
  std::fill(partial_chan_sum_cov_, partial_chan_sum_cov_+nchan_block*nsamp_*(nsamp_-1)/2, _mm256_setzero_si256());
  return true;
#else // defined(__AVX2__)
  throw std::runtime_error("AVX2_Unroll8_WaveformStatsVisitor: AVX2 not supported at compile time");
#endif // defined(__AVX2__)
}

void AVX2_Unroll8_WaveformStatsVisitor::process_8_events()
{
#if defined(__AVX2__)
  const unsigned nsamp_block = (nsamp_ + 15)/16; // number of uint16_t vectors in nsamp
  const unsigned nchan_block = (nchan_ + 15)/16; // number of uint16_t vectors in nchan

  for(unsigned ievent=nkept_events_; ievent<8; ievent++) {
    __m256i*__restrict__ samples = samples_[ievent];
    __m256i*__restrict__ vp = samples;
    for(unsigned isamp_block=0; isamp_block<nsamp_block; isamp_block++) {
      for(unsigned ichan_block=0; ichan_block<nchan_block; ichan_block++) {
        *(vp++) = _mm256_setzero_si256();
      }
    }
  }

  __m256i*__restrict__ sum_base = partial_chan_sum_;
  __m256i*__restrict__ ssq_base = partial_chan_sum_squared_;
  __m256i*__restrict__ cov_base = partial_chan_sum_cov_;

  const __m256i*__restrict__ samples0 = samples_[0];
  const __m256i*__restrict__ samples1 = samples_[1];
  const __m256i*__restrict__ samples2 = samples_[2];
  const __m256i*__restrict__ samples3 = samples_[3];
  const __m256i*__restrict__ samples4 = samples_[4];
  const __m256i*__restrict__ samples5 = samples_[5];
  const __m256i*__restrict__ samples6 = samples_[6];
  const __m256i*__restrict__ samples7 = samples_[7];

  for(unsigned iblock=0;iblock<nchan_block;iblock++)
  {
    const unsigned nchan_block = std::min(16U, nchan_-iblock*16);
    for(unsigned ievent=0; ievent<nkept_events_; ievent++) {
      const ix::iact_data::telescope_event::Waveforms* wf = nullptr;
      if(ievent <= nkept_events_) {
        if(high_gain_) {
          wf = &kept_events_[ievent]->high_gain_image().camera_waveforms();
        } else {
          wf = &kept_events_[ievent]->low_gain_image().camera_waveforms();
        }
      }

      const uint16_t*__restrict__ data =
        reinterpret_cast<const uint16_t*__restrict__>(
          wf->raw_samples_array().data() + wf->raw_samples_array_start());
      const uint16_t*__restrict__ base = data + iblock*nsamp_*16;
      __m256i*__restrict__ samples = samples_[ievent];
      __m256i*__restrict__ vp = samples;
      for(unsigned isamp_block=0; isamp_block<nsamp_block; isamp_block++) {
        for(unsigned ichan_block=0; ichan_block<nchan_block; ichan_block++) {
          *(vp++) = _mm256_loadu_si256((__m256i*)(base + isamp_block*16 + nsamp_*ichan_block));
        }
        calin::math::simd::avx2_m256_swizzle_u16(samples + isamp_block*16);
      }
    }

    for(unsigned isamp=0;isamp<nsamp_;isamp++)
    {
      const __m256i sampi0 = samples0[isamp];
      const __m256i sampi1 = samples1[isamp];
      const __m256i sampi2 = samples2[isamp];
      const __m256i sampi3 = samples3[isamp];
      const __m256i sampi4 = samples4[isamp];
      const __m256i sampi5 = samples5[isamp];
      const __m256i sampi6 = samples6[isamp];
      const __m256i sampi7 = samples7[isamp];

      __m256i sum_l;
      __m256i sum_u;

      sum_l = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(sampi0,0));
      sum_u = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(sampi0,1));
      sum_l = _mm256_add_epi32(sum_l, _mm256_cvtepu16_epi32(_mm256_extracti128_si256(sampi1,0)));
      sum_u = _mm256_add_epi32(sum_l, _mm256_cvtepu16_epi32(_mm256_extracti128_si256(sampi1,1)));
      sum_l = _mm256_add_epi32(sum_l, _mm256_cvtepu16_epi32(_mm256_extracti128_si256(sampi2,0)));
      sum_u = _mm256_add_epi32(sum_l, _mm256_cvtepu16_epi32(_mm256_extracti128_si256(sampi2,1)));
      sum_l = _mm256_add_epi32(sum_l, _mm256_cvtepu16_epi32(_mm256_extracti128_si256(sampi3,0)));
      sum_u = _mm256_add_epi32(sum_l, _mm256_cvtepu16_epi32(_mm256_extracti128_si256(sampi3,1)));
      sum_l = _mm256_add_epi32(sum_l, _mm256_cvtepu16_epi32(_mm256_extracti128_si256(sampi4,0)));
      sum_u = _mm256_add_epi32(sum_l, _mm256_cvtepu16_epi32(_mm256_extracti128_si256(sampi4,1)));
      sum_l = _mm256_add_epi32(sum_l, _mm256_cvtepu16_epi32(_mm256_extracti128_si256(sampi5,0)));
      sum_u = _mm256_add_epi32(sum_l, _mm256_cvtepu16_epi32(_mm256_extracti128_si256(sampi5,1)));
      sum_l = _mm256_add_epi32(sum_l, _mm256_cvtepu16_epi32(_mm256_extracti128_si256(sampi6,0)));
      sum_u = _mm256_add_epi32(sum_l, _mm256_cvtepu16_epi32(_mm256_extracti128_si256(sampi6,1)));
      sum_l = _mm256_add_epi32(sum_l, _mm256_cvtepu16_epi32(_mm256_extracti128_si256(sampi7,0)));
      sum_u = _mm256_add_epi32(sum_l, _mm256_cvtepu16_epi32(_mm256_extracti128_si256(sampi7,1)));

      *sum_base = _mm256_add_epi32(*sum_base, sum_l);
      sum_base++;
      *sum_base = _mm256_add_epi32(*sum_base, sum_u);
      sum_base++;

      __m256i prod_lo;
      __m256i prod_hi;

      prod_lo = _mm256_mullo_epi16(sampi0, sampi0);
      prod_hi = _mm256_mulhi_epu16(sampi0, sampi0);
      sum_l = _mm256_unpacklo_epi16(prod_lo, prod_hi);
      sum_u = _mm256_unpackhi_epi16(prod_lo, prod_hi);
      prod_lo = _mm256_mullo_epi16(sampi1, sampi1);
      prod_hi = _mm256_mulhi_epu16(sampi1, sampi1);
      sum_l = _mm256_add_epi32(sum_l, _mm256_unpacklo_epi16(prod_lo, prod_hi));
      sum_u = _mm256_add_epi32(sum_u, _mm256_unpackhi_epi16(prod_lo, prod_hi));
      prod_lo = _mm256_mullo_epi16(sampi2, sampi2);
      prod_hi = _mm256_mulhi_epu16(sampi2, sampi2);
      sum_l = _mm256_add_epi32(sum_l, _mm256_unpacklo_epi16(prod_lo, prod_hi));
      sum_u = _mm256_add_epi32(sum_u, _mm256_unpackhi_epi16(prod_lo, prod_hi));
      prod_lo = _mm256_mullo_epi16(sampi3, sampi3);
      prod_hi = _mm256_mulhi_epu16(sampi3, sampi3);
      sum_l = _mm256_add_epi32(sum_l, _mm256_unpacklo_epi16(prod_lo, prod_hi));
      sum_u = _mm256_add_epi32(sum_u, _mm256_unpackhi_epi16(prod_lo, prod_hi));
      prod_lo = _mm256_mullo_epi16(sampi4, sampi4);
      prod_hi = _mm256_mulhi_epu16(sampi4, sampi4);
      sum_l = _mm256_add_epi32(sum_l, _mm256_unpacklo_epi16(prod_lo, prod_hi));
      sum_u = _mm256_add_epi32(sum_u, _mm256_unpackhi_epi16(prod_lo, prod_hi));
      prod_lo = _mm256_mullo_epi16(sampi5, sampi5);
      prod_hi = _mm256_mulhi_epu16(sampi5, sampi5);
      sum_l = _mm256_add_epi32(sum_l, _mm256_unpacklo_epi16(prod_lo, prod_hi));
      sum_u = _mm256_add_epi32(sum_u, _mm256_unpackhi_epi16(prod_lo, prod_hi));
      prod_lo = _mm256_mullo_epi16(sampi6, sampi6);
      prod_hi = _mm256_mulhi_epu16(sampi6, sampi6);
      sum_l = _mm256_add_epi32(sum_l, _mm256_unpacklo_epi16(prod_lo, prod_hi));
      sum_u = _mm256_add_epi32(sum_u, _mm256_unpackhi_epi16(prod_lo, prod_hi));
      prod_lo = _mm256_mullo_epi16(sampi7, sampi7);
      prod_hi = _mm256_mulhi_epu16(sampi7, sampi7);
      sum_l = _mm256_add_epi32(sum_l, _mm256_unpacklo_epi16(prod_lo, prod_hi));
      sum_u = _mm256_add_epi32(sum_u, _mm256_unpackhi_epi16(prod_lo, prod_hi));

      *ssq_base = _mm256_add_epi32(*ssq_base, sum_l);
      ssq_base++;
      *ssq_base = _mm256_add_epi32(*ssq_base, sum_u);
      ssq_base++;

      if(calculate_covariance_)
      {
        for(unsigned jsamp=isamp+1;jsamp<nsamp_;jsamp++)
        {
          __m256i sampj;

          sampj = samples0[jsamp];
          prod_lo = _mm256_mullo_epi16(sampi0, sampj);
          prod_hi = _mm256_mulhi_epu16(sampi0, sampj);
          sum_l = _mm256_unpacklo_epi16(prod_lo, prod_hi);
          sum_u = _mm256_unpackhi_epi16(prod_lo, prod_hi);

          sampj = samples1[jsamp];
          prod_lo = _mm256_mullo_epi16(sampi1, sampj);
          prod_hi = _mm256_mulhi_epu16(sampi1, sampj);
          sum_l = _mm256_add_epi32(sum_l, _mm256_unpacklo_epi16(prod_lo, prod_hi));
          sum_u = _mm256_add_epi32(sum_u, _mm256_unpackhi_epi16(prod_lo, prod_hi));

          sampj = samples2[jsamp];
          prod_lo = _mm256_mullo_epi16(sampi2, sampj);
          prod_hi = _mm256_mulhi_epu16(sampi2, sampj);
          sum_l = _mm256_add_epi32(sum_l, _mm256_unpacklo_epi16(prod_lo, prod_hi));
          sum_u = _mm256_add_epi32(sum_u, _mm256_unpackhi_epi16(prod_lo, prod_hi));

          sampj = samples3[jsamp];
          prod_lo = _mm256_mullo_epi16(sampi3, sampj);
          prod_hi = _mm256_mulhi_epu16(sampi3, sampj);
          sum_l = _mm256_add_epi32(sum_l, _mm256_unpacklo_epi16(prod_lo, prod_hi));
          sum_u = _mm256_add_epi32(sum_u, _mm256_unpackhi_epi16(prod_lo, prod_hi));

          sampj = samples4[jsamp];
          prod_lo = _mm256_mullo_epi16(sampi4, sampj);
          prod_hi = _mm256_mulhi_epu16(sampi4, sampj);
          sum_l = _mm256_add_epi32(sum_l, _mm256_unpacklo_epi16(prod_lo, prod_hi));
          sum_u = _mm256_add_epi32(sum_u, _mm256_unpackhi_epi16(prod_lo, prod_hi));

          sampj = samples5[jsamp];
          prod_lo = _mm256_mullo_epi16(sampi5, sampj);
          prod_hi = _mm256_mulhi_epu16(sampi5, sampj);
          sum_l = _mm256_add_epi32(sum_l, _mm256_unpacklo_epi16(prod_lo, prod_hi));
          sum_u = _mm256_add_epi32(sum_u, _mm256_unpackhi_epi16(prod_lo, prod_hi));

          sampj = samples6[jsamp];
          prod_lo = _mm256_mullo_epi16(sampi6, sampj);
          prod_hi = _mm256_mulhi_epu16(sampi6, sampj);
          sum_l = _mm256_add_epi32(sum_l, _mm256_unpacklo_epi16(prod_lo, prod_hi));
          sum_u = _mm256_add_epi32(sum_u, _mm256_unpackhi_epi16(prod_lo, prod_hi));

          sampj = samples7[jsamp];
          prod_lo = _mm256_mullo_epi16(sampi7, sampj);
          prod_hi = _mm256_mulhi_epu16(sampi7, sampj);
          sum_l = _mm256_add_epi32(sum_l, _mm256_unpacklo_epi16(prod_lo, prod_hi));
          sum_u = _mm256_add_epi32(sum_u, _mm256_unpackhi_epi16(prod_lo, prod_hi));

          *cov_base = _mm256_add_epi32(*cov_base, sum_l);
          cov_base++;
          *cov_base = _mm256_add_epi32(*cov_base, sum_u);
          cov_base++;
        }
      }
    }
  }

  for(unsigned ievent=0; ievent<nkept_events_; ievent++) {
    event_lifetime_manager_->release_event(kept_events_[ievent]);
    kept_events_[ievent] = nullptr;
  }
  nkept_events_ = 0;
#else // defined(__AVX2__)
  throw std::runtime_error("AVX2_Unroll8_WaveformStatsVisitor: AVX2 not supported at compile time");
#endif // defined(__AVX2__)
}

void AVX2_Unroll8_WaveformStatsVisitor::merge_partials()
{
#if defined(__AVX2__)

#else // defined(__AVX2__)
  throw std::runtime_error("AVX2_Unroll8_WaveformStatsVisitor: AVX2 not supported at compile time");
#endif // defined(__AVX2__)
}
