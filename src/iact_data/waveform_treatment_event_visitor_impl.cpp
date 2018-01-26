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
#include <algorithm>

#include <util/memory.hpp>
#include <iact_data/waveform_treatment_event_visitor.hpp>
#include <math/simd.hpp>

using namespace calin::iact_data::waveform_treatment_event_visitor;
using calin::util::memory::safe_aligned_recalloc;

void calin::iact_data::waveform_treatment_event_visitor::
SingleGainDualWindowWaveformTreatmentEventVisitor::
scalar_analyze_waveforms(const uint16_t*__restrict__ data)
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
  unsigned isamp = 0;
  int samp[nsamp_];

  for(unsigned ichan=0;ichan<nchan_;ichan++)
  {
    samp[0] = data[ichan*nsamp_];
    imax = 0;
    max = samp[0];
    win = max;
    sum_qt = 0;
    for(isamp = 1;isamp<window_n_;isamp++) {
      const int _samp = data[ichan*nsamp_+isamp];
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
    sum_q = win;
    while(isamp<nsamp_) {
      int iss = isamp-16;
      if(bkg_window_0_ == iss)bkg = win;
      if(sig_window_0_[ichan] == iss)sig = win;
      const int _samp = data[ichan*nsamp_+isamp];
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
    if(bkg_window_0_ == nsamp_-window_n_)bkg = win;
    if(sig_window_0_[ichan] == nsamp_-window_n_)sig = win;

    chan_max_index_[ichan] = imax;
    chan_max_[ichan] = max;
    chan_bkg_win_sum_[ichan] = bkg;
    chan_sig_win_sum_[ichan] = sig;
    chan_sig_max_sum_[ichan] = sig_max;
    chan_sig_max_sum_index_[ichan] = isig_max;
    chan_all_sum_q_[ichan] = sum_q;
    chan_all_sum_qt_[ichan] = sum_qt;
    if(chan_ped_est_[ichan] >= 0) {
      chan_ped_est_[ichan] = ped_iir_old_*chan_ped_est_[ichan] + ped_iir_new_*float(bkg);
    } else {
      chan_ped_est_[ichan] = float(bkg);
    }
    chan_sig_[ichan] = float(sig) - chan_ped_est_[ichan];
    chan_mean_t_[ichan] =
      (double(window_n_*sum_qt) - double(chan_ped_est_[ichan]*nsamp_*(nsamp_-1)/2))/
        (double(window_n_*sum_q) - double(chan_ped_est_[ichan]*nsamp_));
  }
}

void calin::iact_data::waveform_treatment_event_visitor::
AVX2_SingleGainDualWindowWaveformTreatmentEventVisitor::
avx2_analyze_waveforms(const uint16_t*__restrict__ data)
{
#if defined(__AVX2__) and defined(__FMA__)
  const unsigned nv_samp = (nsamp_+15)/16;
  //const unsigned nv_block = nv_samp*16;

  __m256i*__restrict__ samples = samples_;

  const __m256 mean_t_c1 = _mm256_set1_ps(float(nsamp_*(nsamp_-1)/2)/float(window_n_));
  const __m256 mean_t_c2 = _mm256_set1_ps(float(nsamp_)/float(window_n_));

  const unsigned nblock = (nchan_+15)/16;
  unsigned iresvec = 0;

  __m256i visamp;
  __m256i mask;
  __m256i samp;

  for(unsigned iblock=0;iblock<nblock;iblock++)
  {
    const unsigned nvec = std::min(16U, nchan_-iblock*16);
    const uint16_t* base = data + iblock*nsamp_*16;
    __m256i* vp = samples;
    for(unsigned iv_samp=0; iv_samp<nv_samp; iv_samp++) {
      for(unsigned ivec=0; ivec<nvec; ivec++) {
        *(vp++) = _mm256_loadu_si256((__m256i*)(base + iv_samp*16 + nsamp_*ivec));
      }
      calin::math::simd::avx2_m256_swizzle_u16(samples + iv_samp*16);
    }

    __m256i max = samples[0];
    __m256i imax = _mm256_set1_epi16(0);

    __m256i win_l = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(samples[0],0));
    __m256i win_u = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(samples[0],1));

    __m256i sum_qt_l = _mm256_set1_epi32(0);
    __m256i sum_qt_u = _mm256_set1_epi32(0);

    int isamp;
    for(isamp = 1;isamp<window_n_;++isamp) {
      visamp = _mm256_set1_epi16(isamp);
      mask = _mm256_cmpgt_epi16(samples[isamp], max);
      max = _mm256_blendv_epi8(max, samples[isamp], mask);
      imax = _mm256_blendv_epi8(imax, visamp, mask);

      visamp = _mm256_set1_epi32(isamp);

      samp = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(samples[isamp],0));
      win_l = _mm256_add_epi32(win_l, samp);
      sum_qt_l = _mm256_add_epi32(sum_qt_l, _mm256_mullo_epi32(visamp, samp));

      samp = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(samples[isamp],1));
      win_u = _mm256_add_epi32(win_u, samp);
      sum_qt_u = _mm256_add_epi32(sum_qt_u, _mm256_mullo_epi32(visamp, samp));
    }

    __m256i sig_max_l = win_l;
    __m256i sig_max_u = win_u;
    __m256i isig_max_l = _mm256_set1_epi32(0);
    __m256i isig_max_u = _mm256_set1_epi32(0);

    __m256i sig_l = win_l;
    __m256i sig_u = win_u;
    __m256i bkg_l = win_l;
    __m256i bkg_u = win_u;
    __m256i sum_q_l = win_l;
    __m256i sum_q_u = win_u;

    int iss = 0;

    for(;isamp<nsamp_;++isamp)
    {
      visamp = _mm256_set1_epi16(isamp);
      mask = _mm256_cmpgt_epi16(samples[isamp], max);
      max = _mm256_blendv_epi8(max, samples[isamp], mask);
      imax = _mm256_blendv_epi8(imax, visamp, mask);

      visamp = _mm256_set1_epi32(isamp);

      samp = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(samples[isamp],0));
      win_l = _mm256_add_epi32(win_l, samp);
      sum_q_l = _mm256_add_epi32(sum_q_l, samp);
      sum_qt_l = _mm256_add_epi32(sum_qt_l, _mm256_mullo_epi32(visamp, samp));
      samp = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(samples[iss],0));
      win_l = _mm256_sub_epi32(win_l, samp);

      samp = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(samples[isamp],1));
      win_u = _mm256_add_epi32(win_u, samp);
      sum_q_u = _mm256_add_epi32(sum_q_u, samp);
      sum_qt_u = _mm256_add_epi32(sum_qt_u, _mm256_mullo_epi32(visamp, samp));
      samp = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(samples[iss],1));
      win_u = _mm256_sub_epi32(win_u, samp);

      ++iss;
      if(bkg_window_0_ == iss){
        bkg_l = win_l;
        bkg_u = win_u;
      }

      visamp = _mm256_set1_epi32(iss);

      mask = _mm256_cmpgt_epi32(win_l, sig_max_l);
      sig_max_l = _mm256_blendv_epi8(sig_max_l, win_l, mask);
      isig_max_l = _mm256_blendv_epi8(isig_max_l, visamp, mask);
      mask = _mm256_cmpeq_epi32(*(__m256i*)(sig_window_0_ + iblock*16), visamp);
      sig_l = _mm256_blendv_epi8(sig_l, win_l, mask);

      mask = _mm256_cmpgt_epi32(win_u, sig_max_u);
      sig_max_u = _mm256_blendv_epi8(sig_max_u, win_u, mask);
      isig_max_u = _mm256_blendv_epi8(isig_max_u, visamp, mask);
      mask = _mm256_cmpeq_epi32(*(__m256i*)(sig_window_0_ + iblock*16 + 8), visamp);
      sig_u = _mm256_blendv_epi8(sig_u, win_u, mask);
    }

    _mm256_store_si256((__m256i*)chan_max_index_+iresvec, _mm256_cvtepu16_epi32(_mm256_extracti128_si256(imax,0)));
    _mm256_store_si256((__m256i*)chan_max_index_+iresvec+1, _mm256_cvtepu16_epi32(_mm256_extracti128_si256(imax,1)));

    _mm256_store_si256((__m256i*)chan_max_+iresvec, _mm256_cvtepu16_epi32(_mm256_extracti128_si256(max,0)));
    _mm256_store_si256((__m256i*)chan_max_+iresvec+1, _mm256_cvtepu16_epi32(_mm256_extracti128_si256(max,1)));

    _mm256_store_si256((__m256i*)chan_sig_win_sum_+iresvec, sig_l);
    _mm256_store_si256((__m256i*)chan_sig_win_sum_+iresvec+1, sig_u);

    _mm256_store_si256((__m256i*)chan_sig_max_sum_+iresvec, sig_max_l);
    _mm256_store_si256((__m256i*)chan_sig_max_sum_+iresvec+1, sig_max_u);

    _mm256_store_si256((__m256i*)chan_sig_max_sum_index_+iresvec, isig_max_l);
    _mm256_store_si256((__m256i*)chan_sig_max_sum_index_+iresvec+1, isig_max_u);

    _mm256_store_si256((__m256i*)chan_bkg_win_sum_+iresvec, bkg_l);
    _mm256_store_si256((__m256i*)chan_bkg_win_sum_+iresvec+1, bkg_u);

    __m256 fbkg_l = _mm256_cvtepi32_ps(bkg_l);
    __m256 fbkg_u = _mm256_cvtepi32_ps(bkg_u);

    __m256 ped_l = _mm256_load_ps(chan_ped_est_+iresvec*8);
    __m256 ped_n = _mm256_mul_ps(ped_l, _mm256_set1_ps(ped_iir_old_));
    ped_n = _mm256_fmadd_ps(fbkg_l, _mm256_set1_ps(ped_iir_new_), ped_n);
    // This tricky blend uses the sign bit of each pedestal estimate to either
    // select the background sum directly (if pedestal<0) or the filtered
    // update (if pedestal>=0) - the initial pedestal estimates can therefore
    // be set to a negative value to initialize the filter
    ped_l = _mm256_blendv_ps(ped_n, fbkg_l, ped_l);
    _mm256_store_ps(chan_ped_est_+iresvec*8, ped_l);

    __m256 ped_u = _mm256_load_ps(chan_ped_est_+iresvec*8+8);
    ped_n = _mm256_mul_ps(ped_u, _mm256_set1_ps(ped_iir_old_));
    ped_n = _mm256_fmadd_ps(fbkg_u, _mm256_set1_ps(ped_iir_new_), ped_n);
    ped_u = _mm256_blendv_ps(ped_n, fbkg_u, ped_u);
    _mm256_store_ps(chan_ped_est_+iresvec*8+8, ped_u);

    _mm256_store_si256((__m256i*)chan_all_sum_q_+iresvec, sum_q_l);
    _mm256_store_si256((__m256i*)chan_all_sum_q_+iresvec+1, sum_q_u);

    _mm256_store_si256((__m256i*)chan_all_sum_qt_+iresvec, sum_qt_l);
    _mm256_store_si256((__m256i*)chan_all_sum_qt_+iresvec+1, sum_qt_u);

    _mm256_store_ps(chan_sig_+iresvec*8, _mm256_sub_ps(_mm256_cvtepi32_ps(sig_l), ped_l));
    _mm256_store_ps(chan_sig_+iresvec*8+8, _mm256_sub_ps(_mm256_cvtepi32_ps(sig_u), ped_u));

    __m256 nom = _mm256_fmsub_ps(mean_t_c1, ped_l, _mm256_cvtepi32_ps(sum_qt_l));
    __m256 denom = _mm256_fmsub_ps(mean_t_c2, ped_l, _mm256_cvtepi32_ps(sum_q_l));
    _mm256_store_ps(chan_mean_t_+iresvec*8, _mm256_div_ps(nom, denom));

    nom = _mm256_fmsub_ps(mean_t_c1, ped_u, _mm256_cvtepi32_ps(sum_qt_u));
    denom = _mm256_fmsub_ps(mean_t_c2, ped_u, _mm256_cvtepi32_ps(sum_q_u));
    _mm256_store_ps(chan_mean_t_+iresvec*8+8, _mm256_div_ps(nom, denom));

    iresvec += 2;
  }

#else  // defined(__AVX2__) and defined(__FMA__)
  throw std::runtime_error("AVX2_SingleGainDualWindowWaveformTreatmentEventVisitor: AVX2 or FMA not available at compile time");
#endif // defined(__AVX2__) and defined(__FMA__)
}

void calin::iact_data::waveform_treatment_event_visitor::
AVX2_SingleGainDualWindowWaveformTreatmentEventVisitor::
avx2_analyze_waveforms_v2(const uint16_t*__restrict__ data)
{
#if defined(__AVX2__) and defined(__FMA__)
  const unsigned nv_samp = (nsamp_+15)/16;

#if 0
  const unsigned nv_block = nv_samp*16;

  __m256i* samples = new __m256i[nv_block];
  __m256i* q_l = new __m256i[nsamp_];
  __m256i* q_u = new __m256i[nsamp_];
  __m256i* qt_l = new __m256i[nsamp_];
  __m256i* qt_u = new __m256i[nsamp_];
#endif
  __m256i*__restrict__ samples = samples_;
  __m256i*__restrict__ q_l = q_l_;
  __m256i*__restrict__ q_u = q_u_;
  __m256i*__restrict__ qt_l = qt_l_;
  __m256i*__restrict__ qt_u = qt_u_;

  __m256i*__restrict__ sig_window_0_l = (__m256i*)sig_window_0_;
  __m256i*__restrict__ sig_window_0_u = (__m256i*)sig_window_0_ + 1;
  __m256i*__restrict__ chan_max_index = (__m256i*)chan_max_index_;
  __m256i*__restrict__ chan_max = (__m256i*)chan_max_;
  __m256i*__restrict__ chan_sig_win_sum = (__m256i*)chan_sig_win_sum_;
  __m256i*__restrict__ chan_sig_max_sum = (__m256i*)chan_sig_max_sum_;
  __m256i*__restrict__ chan_sig_max_sum_index = (__m256i*)chan_sig_max_sum_index_;
  __m256i*__restrict__ chan_bkg_win_sum = (__m256i*)chan_bkg_win_sum_;
  float*__restrict__ chan_ped_est = chan_ped_est_;
  __m256i*__restrict__ chan_all_sum_q = (__m256i*)chan_all_sum_q_;
  __m256i*__restrict__ chan_all_sum_qt = (__m256i*)chan_all_sum_qt_;
  float*__restrict__ chan_sig = chan_sig_;
  float*__restrict__ chan_mean_t = chan_mean_t_;

  qt_l[0] = qt_u[0] = _mm256_setzero_si256();

  const __m256 mean_t_c1 = _mm256_set1_ps(float(nsamp_*(nsamp_-1)/2)/float(window_n_));
  const __m256 mean_t_c2 = _mm256_set1_ps(float(nsamp_)/float(window_n_));

  const unsigned nblock = (nchan_+15)/16;

  __m256i visamp;
  __m256i viss;
  __m256i mask;
  __m256i one;

  for(unsigned iblock=0;iblock<nblock;iblock++)
  {
    const unsigned nvec = std::min(16U, nchan_-iblock*16);
    const uint16_t*__restrict__ base = data + iblock*nsamp_*16;
    __m256i*__restrict__ vp = samples;
    for(unsigned iv_samp=0; iv_samp<nv_samp; iv_samp++) {
      for(unsigned ivec=0; ivec<nvec; ivec++) {
        *(vp++) = _mm256_loadu_si256((__m256i*)(base + iv_samp*16 + nsamp_*ivec));
      }
      calin::math::simd::avx2_m256_swizzle_u16(samples + iv_samp*16);
    }

    int isamp;
    __m256i max = samples[0];
    __m256i imax = _mm256_setzero_si256();
    visamp = one = _mm256_set1_epi16(1);
    for(isamp = 1;isamp<nsamp_;++isamp) {
      mask = _mm256_cmpgt_epi16(samples[isamp], max);
      max = _mm256_blendv_epi8(max, samples[isamp], mask);
      imax = _mm256_blendv_epi8(imax, visamp, mask);
      visamp = _mm256_add_epi16(visamp, one);
    }

    q_l[0] = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(samples[0],0));
    q_u[0] = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(samples[0],1));
    visamp = one = _mm256_set1_epi32(1);
#if 1
    for(isamp=1; isamp<nsamp_; isamp++) {
      const __m256i samp = samples[isamp];
      q_l[isamp] = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(samp,0));
      qt_l[isamp] = _mm256_mullo_epi32(visamp, q_l[isamp]);
      q_u[isamp] = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(samp,1));
      qt_u[isamp] = _mm256_mullo_epi32(visamp, q_u[isamp]);
      visamp = _mm256_add_epi32(visamp, one);
    }
#else
    for(isamp=1; isamp<nsamp_; isamp++) {
      q_l[isamp] = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(samples[isamp],0));
      qt_l[isamp] = _mm256_mullo_epi32(visamp, q_l[isamp]);
      q_u[isamp] = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(samples[isamp],1));
      qt_u[isamp] = _mm256_mullo_epi32(visamp, q_u[isamp]);
      visamp = _mm256_add_epi32(visamp, one);
    }
#endif

    __m256i win_l = q_l[0];
    __m256i win_u = q_u[0];

    __m256i sum_qt_l = qt_l[0];
    __m256i sum_qt_u = qt_u[0];

    for(isamp = 1;isamp<window_n_;++isamp) {
      win_l = _mm256_add_epi32(win_l, q_l[isamp]);
      sum_qt_l = _mm256_add_epi32(sum_qt_l, qt_l[isamp]);
      win_u = _mm256_add_epi32(win_u, q_u[isamp]);
      sum_qt_u = _mm256_add_epi32(sum_qt_u, qt_u[isamp]);
    }

    __m256i sig_max_l = win_l;
    __m256i sig_max_u = win_u;
    __m256i isig_max_l = _mm256_setzero_si256();
    __m256i isig_max_u = _mm256_setzero_si256();

    __m256i sig_l = win_l;
    __m256i sig_u = win_u;
    __m256i sum_q_l = win_l;
    __m256i sum_q_u = win_u;

    int iss = 0;
    viss = one;

    for(;iss<bkg_window_0_;++iss,++isamp)
    {
      win_l = _mm256_add_epi32(win_l, q_l[isamp]);
      sum_q_l = _mm256_add_epi32(sum_q_l, q_l[isamp]);
      sum_qt_l = _mm256_add_epi32(sum_qt_l, qt_l[isamp]);
      win_l = _mm256_sub_epi32(win_l, q_l[iss]);

      win_u = _mm256_add_epi32(win_u, q_u[isamp]);
      sum_q_u = _mm256_add_epi32(sum_q_u, q_u[isamp]);
      sum_qt_u = _mm256_add_epi32(sum_qt_u, qt_u[isamp]);
      win_u = _mm256_sub_epi32(win_u, q_u[iss]);

      mask = _mm256_cmpgt_epi32(win_l, sig_max_l);
      sig_max_l = _mm256_blendv_epi8(sig_max_l, win_l, mask);
      isig_max_l = _mm256_blendv_epi8(isig_max_l, viss, mask);
      mask = _mm256_cmpeq_epi32(*sig_window_0_l, viss);
      sig_l = _mm256_blendv_epi8(sig_l, win_l, mask);

      mask = _mm256_cmpgt_epi32(win_u, sig_max_u);
      sig_max_u = _mm256_blendv_epi8(sig_max_u, win_u, mask);
      isig_max_u = _mm256_blendv_epi8(isig_max_u, viss, mask);
      mask = _mm256_cmpeq_epi32(*sig_window_0_u, viss);
      sig_u = _mm256_blendv_epi8(sig_u, win_u, mask);

      // visamp = _mm256_add_epi32(visamp, one);
      viss = _mm256_add_epi32(viss, one);
    }

    __m256i bkg_l = win_l;
    __m256i bkg_u = win_u;

    for(;isamp<nsamp_;++iss,++isamp)
    {
      win_l = _mm256_add_epi32(win_l, q_l[isamp]);
      sum_q_l = _mm256_add_epi32(sum_q_l, q_l[isamp]);
      sum_qt_l = _mm256_add_epi32(sum_qt_l, qt_l[isamp]);
      win_l = _mm256_sub_epi32(win_l, q_l[iss]);

      win_u = _mm256_add_epi32(win_u, q_u[isamp]);
      sum_q_u = _mm256_add_epi32(sum_q_u, q_u[isamp]);
      sum_qt_u = _mm256_add_epi32(sum_qt_u, qt_u[isamp]);
      win_u = _mm256_sub_epi32(win_u, q_u[iss]);

      mask = _mm256_cmpgt_epi32(win_l, sig_max_l);
      sig_max_l = _mm256_blendv_epi8(sig_max_l, win_l, mask);
      isig_max_l = _mm256_blendv_epi8(isig_max_l, viss, mask);
      mask = _mm256_cmpeq_epi32(*(__m256i*)(sig_window_0_ + iblock*16), viss);
      sig_l = _mm256_blendv_epi8(sig_l, win_l, mask);

      mask = _mm256_cmpgt_epi32(win_u, sig_max_u);
      sig_max_u = _mm256_blendv_epi8(sig_max_u, win_u, mask);
      isig_max_u = _mm256_blendv_epi8(isig_max_u, viss, mask);
      mask = _mm256_cmpeq_epi32(*(__m256i*)(sig_window_0_ + iblock*16 + 8), viss);
      sig_u = _mm256_blendv_epi8(sig_u, win_u, mask);

      // visamp = _mm256_add_epi32(visamp, one);
      viss = _mm256_add_epi32(viss, one);
    }

    _mm256_store_si256(chan_max_index++, _mm256_cvtepu16_epi32(_mm256_extracti128_si256(imax,0)));
    _mm256_store_si256(chan_max_index++, _mm256_cvtepu16_epi32(_mm256_extracti128_si256(imax,1)));

    _mm256_store_si256(chan_max++, _mm256_cvtepu16_epi32(_mm256_extracti128_si256(max,0)));
    _mm256_store_si256(chan_max++, _mm256_cvtepu16_epi32(_mm256_extracti128_si256(max,1)));

    _mm256_store_si256(chan_sig_win_sum++, sig_l);
    _mm256_store_si256(chan_sig_win_sum++, sig_u);

    _mm256_store_si256(chan_sig_max_sum++, sig_max_l);
    _mm256_store_si256(chan_sig_max_sum++, sig_max_u);

    _mm256_store_si256(chan_sig_max_sum_index++, isig_max_l);
    _mm256_store_si256(chan_sig_max_sum_index++, isig_max_u);

    _mm256_store_si256(chan_bkg_win_sum++, bkg_l);
    _mm256_store_si256(chan_bkg_win_sum++, bkg_u);

    __m256 fbkg_l = _mm256_cvtepi32_ps(bkg_l);
    __m256 fbkg_u = _mm256_cvtepi32_ps(bkg_u);

    __m256 ped_l = _mm256_load_ps(chan_ped_est);
    __m256 ped_n = _mm256_mul_ps(ped_l, _mm256_set1_ps(ped_iir_old_));
    ped_n = _mm256_fmadd_ps(fbkg_l, _mm256_set1_ps(ped_iir_new_), ped_n);
    // This tricky blend uses the sign bit of each pedestal estimate to either
    // select the background sum directly (if pedestal<0) or the filtered
    // update (if pedestal>=0) - the initial pedestal estimates can therefore
    // be set to a negative value to initialize the filter
    ped_l = _mm256_blendv_ps(ped_n, fbkg_l, ped_l);
    _mm256_store_ps(chan_ped_est, ped_l);
    chan_ped_est+=8;

    __m256 ped_u = _mm256_load_ps(chan_ped_est);
    ped_n = _mm256_mul_ps(ped_u, _mm256_set1_ps(ped_iir_old_));
    ped_n = _mm256_fmadd_ps(fbkg_u, _mm256_set1_ps(ped_iir_new_), ped_n);
    ped_u = _mm256_blendv_ps(ped_n, fbkg_u, ped_u);
    _mm256_store_ps(chan_ped_est, ped_u);
    chan_ped_est+=8;

    _mm256_store_si256(chan_all_sum_q++, sum_q_l);
    _mm256_store_si256(chan_all_sum_q++, sum_q_u);

    _mm256_store_si256(chan_all_sum_qt++, sum_qt_l);
    _mm256_store_si256(chan_all_sum_qt++, sum_qt_u);

    _mm256_store_ps(chan_sig, _mm256_sub_ps(_mm256_cvtepi32_ps(sig_l), ped_l));
    chan_sig+=8;
    _mm256_store_ps(chan_sig, _mm256_sub_ps(_mm256_cvtepi32_ps(sig_u), ped_u));
    chan_sig+=8;

    __m256 nom = _mm256_fmsub_ps(mean_t_c1, ped_l, _mm256_cvtepi32_ps(sum_qt_l));
    __m256 denom = _mm256_fmsub_ps(mean_t_c2, ped_l, _mm256_cvtepi32_ps(sum_q_l));
    _mm256_store_ps(chan_mean_t, _mm256_div_ps(nom, denom));
    chan_mean_t+=8;

    nom = _mm256_fmsub_ps(mean_t_c1, ped_u, _mm256_cvtepi32_ps(sum_qt_u));
    denom = _mm256_fmsub_ps(mean_t_c2, ped_u, _mm256_cvtepi32_ps(sum_q_u));
    _mm256_store_ps(chan_mean_t, _mm256_div_ps(nom, denom));
    chan_mean_t+=8;

    sig_window_0_l += 2;
    sig_window_0_u += 2;
  }

#if 0
  delete[] samples;
  delete[] q_l;
  delete[] q_u;
  delete[] qt_l;
  delete[] qt_u;
#endif

#else  // defined(__AVX2__) and defined(__FMA__)
  throw std::runtime_error("AVX2_SingleGainDualWindowWaveformTreatmentEventVisitor: AVX2 or FMA not available at compile time");
#endif // defined(__AVX2__) and defined(__FMA__)
}








void calin::iact_data::waveform_treatment_event_visitor::
AVX2_SingleGainDualWindowWaveformTreatmentEventVisitor::
avx2_analyze_waveforms_v3(const uint16_t*__restrict__ data)
{
#if defined(__AVX2__) and defined(__FMA__)
  const unsigned nv_samp = (nsamp_+15)/16;

#if 0
  const unsigned nv_block = nv_samp*16;

  __m256i* samples = new __m256i[nv_block];
  __m256i* q_l = new __m256i[nsamp_];
  __m256i* q_u = new __m256i[nsamp_];
  __m256i* qt_l = new __m256i[nsamp_];
  __m256i* qt_u = new __m256i[nsamp_];
#endif
  __m256i*__restrict__ samples = samples_;
  __m256i*__restrict__ q_l = q_l_;
  __m256i*__restrict__ q_u = q_u_;

  __m256i*__restrict__ sig_window_0_l = (__m256i*)sig_window_0_;
  __m256i*__restrict__ sig_window_0_u = (__m256i*)sig_window_0_ + 1;
  __m256i*__restrict__ chan_max_index = (__m256i*)chan_max_index_;
  __m256i*__restrict__ chan_max = (__m256i*)chan_max_;
  __m256i*__restrict__ chan_sig_win_sum = (__m256i*)chan_sig_win_sum_;
  __m256i*__restrict__ chan_sig_max_sum = (__m256i*)chan_sig_max_sum_;
  __m256i*__restrict__ chan_sig_max_sum_index = (__m256i*)chan_sig_max_sum_index_;
  __m256i*__restrict__ chan_bkg_win_sum = (__m256i*)chan_bkg_win_sum_;
  float*__restrict__ chan_ped_est = chan_ped_est_;
  __m256i*__restrict__ chan_all_sum_q = (__m256i*)chan_all_sum_q_;
  __m256i*__restrict__ chan_all_sum_qt = (__m256i*)chan_all_sum_qt_;
  float*__restrict__ chan_sig = chan_sig_;
  float*__restrict__ chan_mean_t = chan_mean_t_;

  const __m256 mean_t_c1 = _mm256_set1_ps(float(nsamp_*(nsamp_-1)/2)/float(window_n_));
  const __m256 mean_t_c2 = _mm256_set1_ps(float(nsamp_)/float(window_n_));

  const unsigned nblock = (nchan_+15)/16;

  __m256i visamp;
  __m256i viss;
  __m256i mask;
  __m256i one;

  for(unsigned iblock=0;iblock<nblock;iblock++)
  {
    const unsigned nvec = std::min(16U, nchan_-iblock*16);
    const uint16_t*__restrict__ base = data + iblock*nsamp_*16;
    __m256i*__restrict__ vp = samples;
    for(unsigned iv_samp=0; iv_samp<nv_samp; iv_samp++) {
      for(unsigned ivec=0; ivec<nvec; ivec++) {
        *(vp++) = _mm256_loadu_si256((__m256i*)(base + iv_samp*16 + nsamp_*ivec));
      }
      calin::math::simd::avx2_m256_swizzle_u16(samples + iv_samp*16);
    }

    int isamp;
    __m256i max = samples[0];
    __m256i imax = _mm256_setzero_si256();
    visamp = one = _mm256_set1_epi16(1);
    for(isamp = 1;isamp<nsamp_;++isamp) {
      mask = _mm256_cmpgt_epi16(samples[isamp], max);
      max = _mm256_blendv_epi8(max, samples[isamp], mask);
      imax = _mm256_blendv_epi8(imax, visamp, mask);
      visamp = _mm256_add_epi16(visamp, one);
    }

    __m256i sum_qt_l = _mm256_setzero_si256();
    __m256i sum_qt_u = _mm256_setzero_si256();

    q_l[0] = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(samples[0],0));
    q_u[0] = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(samples[0],1));
    visamp = one = _mm256_set1_epi32(1);
    if(nsamp_ & 0x1) {
      __m256i sum_qt_l2 = _mm256_setzero_si256();
      __m256i sum_qt_u2 = _mm256_setzero_si256();

      for(isamp=1; isamp<nsamp_; isamp++) {
        __m256i samp16 = samples[isamp];
        __m256i samp32 = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(samp16,0));
        q_l[isamp] = samp32;
        sum_qt_l = _mm256_add_epi32(sum_qt_l, _mm256_mullo_epi32(visamp, samp32));
        samp32 = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(samp16,1));
        q_u[isamp] = samp32;
        sum_qt_u = _mm256_add_epi32(sum_qt_u, _mm256_mullo_epi32(visamp, samp32));
        visamp = _mm256_add_epi32(visamp, one);

        samp16 = samples[++isamp];
        samp32 = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(samp16,0));
        q_l[isamp] = samp32;
        sum_qt_l2 = _mm256_add_epi32(sum_qt_l2, _mm256_mullo_epi32(visamp, samp32));
        samp32 = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(samp16,1));
        q_u[isamp] = samp32;
        sum_qt_u2 = _mm256_add_epi32(sum_qt_u2, _mm256_mullo_epi32(visamp, samp32));
        visamp = _mm256_add_epi32(visamp, one);
      }
      sum_qt_l = _mm256_add_epi32(sum_qt_l, sum_qt_l2);
      sum_qt_u = _mm256_add_epi32(sum_qt_u, sum_qt_u2);
    } else {
      __m256i samp16 = samples[1];
      __m256i samp32 = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(samp16,0));
      q_l[1] = samp32;
      __m256i sum_qt_l2 = _mm256_mullo_epi32(visamp, samp32);
      samp32 = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(samp16,1));
      q_u[1] = samp32;
      __m256i sum_qt_u2 = _mm256_mullo_epi32(visamp, samp32);
      visamp = _mm256_add_epi32(visamp, one);
      for(isamp=2; isamp<nsamp_; isamp++) {
        samp16 = samples[isamp];
        samp32 = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(samp16,0));
        q_l[isamp] = samp32;
        sum_qt_l = _mm256_add_epi32(sum_qt_l, _mm256_mullo_epi32(visamp, samp32));
        samp32 = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(samp16,1));
        q_u[isamp] = samp32;
        sum_qt_u = _mm256_add_epi32(sum_qt_u, _mm256_mullo_epi32(visamp, samp32));
        visamp = _mm256_add_epi32(visamp, one);

        samp16 = samples[++isamp];
        samp32 = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(samp16,0));
        q_l[isamp] = samp32;
        sum_qt_l2 = _mm256_add_epi32(sum_qt_l2, _mm256_mullo_epi32(visamp, samp32));
        samp32 = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(samp16,1));
        q_u[isamp] = samp32;
        sum_qt_u2 = _mm256_add_epi32(sum_qt_u2, _mm256_mullo_epi32(visamp, samp32));
        visamp = _mm256_add_epi32(visamp, one);
      }
      sum_qt_l = _mm256_add_epi32(sum_qt_l, sum_qt_l2);
      sum_qt_u = _mm256_add_epi32(sum_qt_u, sum_qt_u2);
    }

    q_l[0] = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(samples[0],0));
    q_u[0] = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(samples[0],1));
    visamp = one = _mm256_set1_epi32(1);
#if 1
    for(isamp=1; isamp<nsamp_; isamp++) {
      const __m256i samp = samples[isamp];
      q_l[isamp] = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(samp,0));
      q_u[isamp] = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(samp,1));
      visamp = _mm256_add_epi32(visamp, one);
    }
#else
    for(isamp=1; isamp<nsamp_; isamp++) {
      q_l[isamp] = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(samples[isamp],0));
      qt_l[isamp] = _mm256_mullo_epi32(visamp, q_l[isamp]);
      q_u[isamp] = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(samples[isamp],1));
      qt_u[isamp] = _mm256_mullo_epi32(visamp, q_u[isamp]);
      visamp = _mm256_add_epi32(visamp, one);
    }
#endif

    __m256i win_l = q_l[0];
    __m256i win_u = q_u[0];

    for(isamp = 1;isamp<window_n_;++isamp) {
      win_l = _mm256_add_epi32(win_l, q_l[isamp]);
      win_u = _mm256_add_epi32(win_u, q_u[isamp]);
    }

    __m256i sig_max_l = win_l;
    __m256i sig_max_u = win_u;
    __m256i isig_max_l = _mm256_setzero_si256();
    __m256i isig_max_u = _mm256_setzero_si256();

    __m256i sig_l = win_l;
    __m256i sig_u = win_u;
    __m256i sum_q_l = win_l;
    __m256i sum_q_u = win_u;

    int iss = 0;
    viss = one;

    for(;iss<bkg_window_0_;++iss,++isamp)
    {
      win_l = _mm256_add_epi32(win_l, q_l[isamp]);
      sum_q_l = _mm256_add_epi32(sum_q_l, q_l[isamp]);
      win_l = _mm256_sub_epi32(win_l, q_l[iss]);

      win_u = _mm256_add_epi32(win_u, q_u[isamp]);
      sum_q_u = _mm256_add_epi32(sum_q_u, q_u[isamp]);
      win_u = _mm256_sub_epi32(win_u, q_u[iss]);

      mask = _mm256_cmpgt_epi32(win_l, sig_max_l);
      sig_max_l = _mm256_blendv_epi8(sig_max_l, win_l, mask);
      isig_max_l = _mm256_blendv_epi8(isig_max_l, viss, mask);
      mask = _mm256_cmpeq_epi32(*sig_window_0_l, viss);
      sig_l = _mm256_blendv_epi8(sig_l, win_l, mask);

      mask = _mm256_cmpgt_epi32(win_u, sig_max_u);
      sig_max_u = _mm256_blendv_epi8(sig_max_u, win_u, mask);
      isig_max_u = _mm256_blendv_epi8(isig_max_u, viss, mask);
      mask = _mm256_cmpeq_epi32(*sig_window_0_u, viss);
      sig_u = _mm256_blendv_epi8(sig_u, win_u, mask);

      // visamp = _mm256_add_epi32(visamp, one);
      viss = _mm256_add_epi32(viss, one);
    }

    __m256i bkg_l = win_l;
    __m256i bkg_u = win_u;

    for(;isamp<nsamp_;++iss,++isamp)
    {
      win_l = _mm256_add_epi32(win_l, q_l[isamp]);
      sum_q_l = _mm256_add_epi32(sum_q_l, q_l[isamp]);
      win_l = _mm256_sub_epi32(win_l, q_l[iss]);

      win_u = _mm256_add_epi32(win_u, q_u[isamp]);
      sum_q_u = _mm256_add_epi32(sum_q_u, q_u[isamp]);
      win_u = _mm256_sub_epi32(win_u, q_u[iss]);

      mask = _mm256_cmpgt_epi32(win_l, sig_max_l);
      sig_max_l = _mm256_blendv_epi8(sig_max_l, win_l, mask);
      isig_max_l = _mm256_blendv_epi8(isig_max_l, viss, mask);
      mask = _mm256_cmpeq_epi32(*(__m256i*)(sig_window_0_ + iblock*16), viss);
      sig_l = _mm256_blendv_epi8(sig_l, win_l, mask);

      mask = _mm256_cmpgt_epi32(win_u, sig_max_u);
      sig_max_u = _mm256_blendv_epi8(sig_max_u, win_u, mask);
      isig_max_u = _mm256_blendv_epi8(isig_max_u, viss, mask);
      mask = _mm256_cmpeq_epi32(*(__m256i*)(sig_window_0_ + iblock*16 + 8), viss);
      sig_u = _mm256_blendv_epi8(sig_u, win_u, mask);

      // visamp = _mm256_add_epi32(visamp, one);
      viss = _mm256_add_epi32(viss, one);
    }

    _mm256_store_si256(chan_max_index++, _mm256_cvtepu16_epi32(_mm256_extracti128_si256(imax,0)));
    _mm256_store_si256(chan_max_index++, _mm256_cvtepu16_epi32(_mm256_extracti128_si256(imax,1)));

    _mm256_store_si256(chan_max++, _mm256_cvtepu16_epi32(_mm256_extracti128_si256(max,0)));
    _mm256_store_si256(chan_max++, _mm256_cvtepu16_epi32(_mm256_extracti128_si256(max,1)));

    _mm256_store_si256(chan_sig_win_sum++, sig_l);
    _mm256_store_si256(chan_sig_win_sum++, sig_u);

    _mm256_store_si256(chan_sig_max_sum++, sig_max_l);
    _mm256_store_si256(chan_sig_max_sum++, sig_max_u);

    _mm256_store_si256(chan_sig_max_sum_index++, isig_max_l);
    _mm256_store_si256(chan_sig_max_sum_index++, isig_max_u);

    _mm256_store_si256(chan_bkg_win_sum++, bkg_l);
    _mm256_store_si256(chan_bkg_win_sum++, bkg_u);

    __m256 fbkg_l = _mm256_cvtepi32_ps(bkg_l);
    __m256 fbkg_u = _mm256_cvtepi32_ps(bkg_u);

    __m256 ped_l = _mm256_load_ps(chan_ped_est);
    __m256 ped_n = _mm256_mul_ps(ped_l, _mm256_set1_ps(ped_iir_old_));
    ped_n = _mm256_fmadd_ps(fbkg_l, _mm256_set1_ps(ped_iir_new_), ped_n);
    // This tricky blend uses the sign bit of each pedestal estimate to either
    // select the background sum directly (if pedestal<0) or the filtered
    // update (if pedestal>=0) - the initial pedestal estimates can therefore
    // be set to a negative value to initialize the filter
    ped_l = _mm256_blendv_ps(ped_n, fbkg_l, ped_l);
    _mm256_store_ps(chan_ped_est, ped_l);
    chan_ped_est+=8;

    __m256 ped_u = _mm256_load_ps(chan_ped_est);
    ped_n = _mm256_mul_ps(ped_u, _mm256_set1_ps(ped_iir_old_));
    ped_n = _mm256_fmadd_ps(fbkg_u, _mm256_set1_ps(ped_iir_new_), ped_n);
    ped_u = _mm256_blendv_ps(ped_n, fbkg_u, ped_u);
    _mm256_store_ps(chan_ped_est, ped_u);
    chan_ped_est+=8;

    _mm256_store_si256(chan_all_sum_q++, sum_q_l);
    _mm256_store_si256(chan_all_sum_q++, sum_q_u);

    _mm256_store_si256(chan_all_sum_qt++, sum_qt_l);
    _mm256_store_si256(chan_all_sum_qt++, sum_qt_u);

    _mm256_store_ps(chan_sig, _mm256_sub_ps(_mm256_cvtepi32_ps(sig_l), ped_l));
    chan_sig+=8;
    _mm256_store_ps(chan_sig, _mm256_sub_ps(_mm256_cvtepi32_ps(sig_u), ped_u));
    chan_sig+=8;

    __m256 nom = _mm256_fmsub_ps(mean_t_c1, ped_l, _mm256_cvtepi32_ps(sum_qt_l));
    __m256 denom = _mm256_fmsub_ps(mean_t_c2, ped_l, _mm256_cvtepi32_ps(sum_q_l));
    _mm256_store_ps(chan_mean_t, _mm256_div_ps(nom, denom));
    chan_mean_t+=8;

    nom = _mm256_fmsub_ps(mean_t_c1, ped_u, _mm256_cvtepi32_ps(sum_qt_u));
    denom = _mm256_fmsub_ps(mean_t_c2, ped_u, _mm256_cvtepi32_ps(sum_q_u));
    _mm256_store_ps(chan_mean_t, _mm256_div_ps(nom, denom));
    chan_mean_t+=8;

    sig_window_0_l += 2;
    sig_window_0_u += 2;
  }

#if 0
  delete[] samples;
  delete[] q_l;
  delete[] q_u;
  delete[] qt_l;
  delete[] qt_u;
#endif

#else  // defined(__AVX2__) and defined(__FMA__)
  throw std::runtime_error("AVX2_SingleGainDualWindowWaveformTreatmentEventVisitor: AVX2 or FMA not available at compile time");
#endif // defined(__AVX2__) and defined(__FMA__)
}
