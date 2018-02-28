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

AVX2_Unroll8_WaveformStatsParallelVisitor::
AVX2_Unroll8_WaveformStatsParallelVisitor(bool high_gain, bool calculate_covariance):
  ParallelEventVisitor(),
  high_gain_(high_gain), calculate_covariance_(calculate_covariance)
{
#if defined(__AVX2__)
  // nothing to see here
#else // defined(__AVX2__)
  throw std::runtime_error("AVX2_Unroll8_WaveformStatsParallelVisitor: AVX2 not supported at compile time");
#endif // defined(__AVX2__)
}

AVX2_Unroll8_WaveformStatsParallelVisitor::~AVX2_Unroll8_WaveformStatsParallelVisitor()
{
#if defined(__AVX2__)
  for(unsigned i=0;i<8;i++)free(samples_[i]);
  free(partial_chan_nevent_);
  free(partial_chan_sum_);
  free(partial_chan_sum_squared_);
  if(calculate_covariance_)free(partial_chan_sum_cov_);
#else // defined(__AVX2__)
  // throw std::runtime_error("AVX2_Unroll8_WaveformStatsParallelVisitor: AVX2 not supported at compile time");
#endif // defined(__AVX2__)
}

AVX2_Unroll8_WaveformStatsParallelVisitor* AVX2_Unroll8_WaveformStatsParallelVisitor::
new_sub_visitor(
  const std::map<calin::iact_data::event_visitor::ParallelEventVisitor*,
      calin::iact_data::event_visitor::ParallelEventVisitor*>&
    antecedent_visitors)
{
#if defined(__AVX2__)
  AVX2_Unroll8_WaveformStatsParallelVisitor* child =
    new AVX2_Unroll8_WaveformStatsParallelVisitor(high_gain_, calculate_covariance_);
  child->parent_ = this;
  return child;
#else // defined(__AVX2__)
  throw std::runtime_error("AVX2_Unroll8_WaveformStatsParallelVisitor: AVX2 not supported at compile time");
#endif // defined(__AVX2__)
}

bool AVX2_Unroll8_WaveformStatsParallelVisitor::visit_telescope_run(
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
    calin::ix::diagnostics::waveform::WaveformRawStats* stats;
    if(high_gain_) {
      stats = results_.add_high_gain();
    } else {
      stats = results_.add_low_gain();
    }

    stats->mutable_sum()->Resize(nsamp_,0);
    stats->mutable_sum_squared()->Resize(nsamp_,0);
    if(calculate_covariance_)
      stats->mutable_sum_product()->Resize(nsamp_*(nsamp_-1)/2,0);

    stats->set_num_entries(0);
    std::fill(stats->mutable_sum()->begin(), stats->mutable_sum()->end(), 0);
    std::fill(stats->mutable_sum_squared()->begin(), stats->mutable_sum_squared()->end(), 0);
    if(calculate_covariance_)
      std::fill(stats->mutable_sum_product()->begin(), stats->mutable_sum_product()->end(), 0);
  }

  unsigned nsamp_block = (nsamp_ + 15)/16; // number of uint16_t vectors in nsamp
  unsigned nchan_block = (nchan_ + 15)/16; // number of uint16_t vectors in nchan

  nkept_events_ = 0;

  if(nsamp_ != old_nsamp or samples_[0] == nullptr) {
    for(unsigned i=0;i<8;i++)safe_aligned_recalloc(samples_[i], nsamp_block*16);
  }

  partial_num_entries_ = 0;
  if(nsamp_ != old_nchan or partial_chan_nevent_ == nullptr) {
    safe_aligned_recalloc(partial_chan_nevent_, nchan_);
    safe_aligned_recalloc(partial_chan_sum_, nchan_block*nsamp_*2);
    safe_aligned_recalloc(partial_chan_sum_squared_, nchan_block*nsamp_*2);
    if(calculate_covariance_)
      safe_aligned_recalloc(partial_chan_sum_cov_, nchan_block*nsamp_*(nsamp_-1));
  }

  std::fill(partial_chan_nevent_, partial_chan_nevent_+nchan_, 0);
  std::fill(partial_chan_sum_, partial_chan_sum_+nchan_block*nsamp_*2, _mm256_setzero_si256());
  std::fill(partial_chan_sum_squared_, partial_chan_sum_squared_+nchan_block*nsamp_*2, _mm256_setzero_si256());
  std::fill(partial_chan_sum_cov_, partial_chan_sum_cov_+nchan_block*nsamp_*(nsamp_-1), _mm256_setzero_si256());

  event_lifetime_manager_ = event_lifetime_manager;
  return true;
#else // defined(__AVX2__)
  throw std::runtime_error("AVX2_Unroll8_WaveformStatsParallelVisitor: AVX2 not supported at compile time");
#endif // defined(__AVX2__)
}

bool AVX2_Unroll8_WaveformStatsParallelVisitor::leave_telescope_run()
{
#if defined(__AVX2__)
  if(nkept_events_)process_8_events();
  merge_partials();
  event_lifetime_manager_ = nullptr;
  return true;
#else // defined(__AVX2__)
  throw std::runtime_error("AVX2_Unroll8_WaveformStatsParallelVisitor: AVX2 not supported at compile time");
#endif // defined(__AVX2__)
}

bool AVX2_Unroll8_WaveformStatsParallelVisitor::visit_telescope_event(uint64_t seq_index,
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
    if(partial_num_entries_ >= partial_max_num_entries_)merge_partials();
  }
  return true;
#else // defined(__AVX2__)
  throw std::runtime_error("AVX2_Unroll8_WaveformStatsParallelVisitor: AVX2 not supported at compile time");
#endif // defined(__AVX2__)
}

bool AVX2_Unroll8_WaveformStatsParallelVisitor::merge_results()
{
#if defined(__AVX2__)
  if(parent_)parent_->results_.IntegrateFrom(results_);
  return true;
#else // defined(__AVX2__)
  throw std::runtime_error("AVX2_Unroll8_WaveformStatsParallelVisitor: AVX2 not supported at compile time");
#endif // defined(__AVX2__)
}

#if defined(__AVX2__)
void AVX2_Unroll8_WaveformStatsParallelVisitor::process_8_events()
{
  const unsigned nsamp_block = (nsamp_ + 15)/16; // number of uint16_t vectors in nsamp
  const unsigned nchan_block = (nchan_ + 15)/16; // number of uint16_t vectors in nchan

  for(unsigned ievent=nkept_events_; ievent<8; ievent++) {
    __m256i*__restrict__ samples = samples_[ievent];
    __m256i*__restrict__ vp = samples;
    for(unsigned isamp_block=0; isamp_block<nsamp_block; isamp_block++) {
      for(unsigned ichan_block=0; ichan_block<16; ichan_block++) {
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
    const unsigned nvec = std::min(16U, nchan_-iblock*16);
    for(unsigned ievent=0; ievent<nkept_events_; ievent++) {
      const ix::iact_data::telescope_event::Waveforms* wf = nullptr;
      if(ievent <= nkept_events_) {
        if(high_gain_) {
          wf = &kept_events_[ievent]->high_gain_image().camera_waveforms();
        } else {
          wf = &kept_events_[ievent]->low_gain_image().camera_waveforms();
        }
      }

#if 0
      std::cerr << "process_8_events loop : " << ievent << ' '
        << kept_events_[ievent] << ' ' << wf << ' '
        << reinterpret_cast<const void*>(wf->raw_samples_array().data()) << ' '
        << wf->raw_samples_array().size() << ' '
        << wf->raw_samples_array_start() << ' '
        << iblock*nsamp_*16 << '\n';
#endif

      const uint16_t*__restrict__ data =
        reinterpret_cast<const uint16_t*__restrict__>(
          wf->raw_samples_array().data() + wf->raw_samples_array_start());
      const uint16_t*__restrict__ base = data + iblock*nsamp_*16;
      __m256i*__restrict__ samples = samples_[ievent];
      __m256i*__restrict__ vp = samples;
      for(unsigned isamp_block=0; isamp_block<nsamp_block; isamp_block++) {
        for(unsigned ichan_block=0; ichan_block<nvec; ichan_block++) {
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
      sum_u = _mm256_add_epi32(sum_u, _mm256_cvtepu16_epi32(_mm256_extracti128_si256(sampi1,1)));
      sum_l = _mm256_add_epi32(sum_l, _mm256_cvtepu16_epi32(_mm256_extracti128_si256(sampi2,0)));
      sum_u = _mm256_add_epi32(sum_u, _mm256_cvtepu16_epi32(_mm256_extracti128_si256(sampi2,1)));
      sum_l = _mm256_add_epi32(sum_l, _mm256_cvtepu16_epi32(_mm256_extracti128_si256(sampi3,0)));
      sum_u = _mm256_add_epi32(sum_u, _mm256_cvtepu16_epi32(_mm256_extracti128_si256(sampi3,1)));
      sum_l = _mm256_add_epi32(sum_l, _mm256_cvtepu16_epi32(_mm256_extracti128_si256(sampi4,0)));
      sum_u = _mm256_add_epi32(sum_u, _mm256_cvtepu16_epi32(_mm256_extracti128_si256(sampi4,1)));
      sum_l = _mm256_add_epi32(sum_l, _mm256_cvtepu16_epi32(_mm256_extracti128_si256(sampi5,0)));
      sum_u = _mm256_add_epi32(sum_u, _mm256_cvtepu16_epi32(_mm256_extracti128_si256(sampi5,1)));
      sum_l = _mm256_add_epi32(sum_l, _mm256_cvtepu16_epi32(_mm256_extracti128_si256(sampi6,0)));
      sum_u = _mm256_add_epi32(sum_u, _mm256_cvtepu16_epi32(_mm256_extracti128_si256(sampi6,1)));
      sum_l = _mm256_add_epi32(sum_l, _mm256_cvtepu16_epi32(_mm256_extracti128_si256(sampi7,0)));
      sum_u = _mm256_add_epi32(sum_u, _mm256_cvtepu16_epi32(_mm256_extracti128_si256(sampi7,1)));

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

  unsigned*__restrict__ chan_nevent = partial_chan_nevent_;
  const unsigned nchan = nchan_;
  const unsigned nkept = nkept_events_;
  for(unsigned ichan=0; ichan<nchan; ichan++) {
    *(chan_nevent++) += nkept;
  }

  partial_num_entries_ += nkept_events_;

  for(unsigned ievent=0; ievent<nkept_events_; ievent++) {
    event_lifetime_manager_->release_event(kept_events_[ievent]);
    kept_events_[ievent] = nullptr;
  }
  nkept_events_ = 0;
}

void AVX2_Unroll8_WaveformStatsParallelVisitor::merge_partials()
{
  const unsigned nchan = nchan_;
  const unsigned nsamp = nsamp_;
  const unsigned nchan_block = (nchan + 15)/16; // number of uint16_t vectors in nchan

  __m256i*__restrict__ sum_base = partial_chan_sum_;
  __m256i*__restrict__ ssq_base = partial_chan_sum_squared_;
  __m256i*__restrict__ cov_base = partial_chan_sum_cov_;

  for(unsigned ichan=0,iblock=0,iel=0; ichan<nchan; ichan++,iel++) {
    if(iel==16)iblock++, iel=0;

    calin::ix::diagnostics::waveform::WaveformRawStats* stats;
    if(high_gain_)stats = results_.mutable_high_gain(ichan);
    else stats = results_.mutable_low_gain(ichan);

    stats->set_num_entries(stats->num_entries() + partial_chan_nevent_[ichan]);
    for(unsigned isamp=0,ijsamp=0; isamp<nsamp; isamp++) {
      stats->set_sum(isamp, stats->sum(isamp) +
        reinterpret_cast<uint32_t*>(&sum_base[(iblock*nsamp+isamp)*2])[iel]);
      stats->set_sum_squared(isamp, stats->sum_squared(isamp) +
        reinterpret_cast<uint32_t*>(&ssq_base[(iblock*nsamp+isamp)*2])[iel]);
      if(calculate_covariance_)
      {
        for(unsigned jsamp=isamp+1;jsamp<nsamp_;jsamp++,ijsamp++) {
          stats->set_sum_product(ijsamp, stats->sum_product(ijsamp) +
            reinterpret_cast<uint32_t*>(&cov_base[iblock*nsamp*(nsamp-1) + ijsamp*2])[iel]);
        }
      }
    }
  }

  partial_num_entries_ = 0;

  std::fill(partial_chan_nevent_, partial_chan_nevent_+nchan_, 0);
  std::fill(partial_chan_sum_, partial_chan_sum_+nchan_block*nsamp*2, _mm256_setzero_si256());
  std::fill(partial_chan_sum_squared_, partial_chan_sum_squared_+nchan_block*nsamp*2, _mm256_setzero_si256());
  std::fill(partial_chan_sum_cov_, partial_chan_sum_cov_+nchan_block*nsamp*(nsamp-1), _mm256_setzero_si256());
}
#endif // defined(__AVX2__)
