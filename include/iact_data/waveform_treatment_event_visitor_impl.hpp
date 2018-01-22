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

#include <math/simd.hpp>

inline void calin::iact_data::waveform_treatment_event_visitor::
SingleGainDualWindowWaveformTreatmentEventVisitor::
analyze_waveforms(
  const uint16_t*__restrict__ data, unsigned nchan, int nsamp,
  int window_n, int bkg_window_0, const int* sig_window_0,
  float*__restrict__ ped, float ped_iir_old, float ped_iir_new,
  int*__restrict__ chan_max_index, int*__restrict__ chan_max,
  int*__restrict__ chan_bkg_win_sum, int* chan_sig_win_sum,
  int*__restrict__ chan_sig_max_sum, int*__restrict__ chan_sig_max_sum_index,
  int*__restrict__ chan_all_sum_q, int*__restrict__ chan_all_sum_qt,
  float*__restrict__ chan_sig, float*__restrict__ chan_mean_t)
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
  int samp[nsamp];

  for(unsigned ichan=0;ichan<nchan;ichan++)
  {
    samp[0] = data[ichan*nsamp];
    imax = 0;
    max = samp[0];
    win = max;
    sum_qt = 0;
    for(isamp = 1;isamp<window_n;isamp++) {
      const int _samp = data[ichan*nsamp+isamp];
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
    while(isamp<nsamp) {
      int iss = isamp-16;
      if(bkg_window_0 == iss)bkg = win;
      if(sig_window_0[ichan] == iss)sig = win;
      const int _samp = data[ichan*nsamp+isamp];
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
      (double(window_n*sum_qt) - double(ped[ichan]*nsamp*(nsamp-1)/2))/
        (double(window_n*sum_q) - double(ped[ichan]*nsamp));
  }
}

inline void calin::iact_data::waveform_treatment_event_visitor::
AVX2_SingleGainDualWindowWaveformTreatmentEventVisitor::
avx2_analyze_waveforms(
  const uint16_t*__restrict__ data, unsigned nchan, int nsamp,
  int window_n, int bkg_window_0, const int* sig_window_0,
  float*__restrict__ ped, float ped_iir_old, float ped_iir_new,
  int*__restrict__ chan_max_index, int*__restrict__ chan_max,
  int*__restrict__ chan_bkg_win_sum, int* chan_sig_win_sum,
  int*__restrict__ chan_sig_max_sum, int*__restrict__ chan_sig_max_sum_index,
  int*__restrict__ chan_all_sum_q, int*__restrict__ chan_all_sum_qt,
  float*__restrict__ chan_sig, float*__restrict__ chan_mean_t)
{
#if defined(__AVX2__) and defined(__FMA__)
  const unsigned nv_samp = (nsamp+15)/16;
  const unsigned nv_block = nv_samp*16;

  __m256i* samples = new __m256i[nv_block];

  __m256 mean_t_c1 = _mm256_set1_ps(float(nsamp*(nsamp-1)/2)/float(window_n));
  __m256 mean_t_c2 = _mm256_set1_ps(float(nsamp)/float(window_n));

  const unsigned nblock = (nchan+15)/16;
  for(unsigned iblock=0;iblock<nblock;iblock++)
  {
    const uint16_t* base = data + iblock*nsamp*16;
    __m256i* vp = samples;
    for(unsigned iv_samp=0; iv_samp<nv_samp; iv_samp++) {
      for(unsigned ivec=0; ivec<16; ivec++) {
        *(vp++) = _mm256_loadu_si256((__m256i*)(base + iv_samp*16 + nsamp*ivec));
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
    for(isamp = 1;isamp<window_n;++isamp) {
      __m256i visamp = _mm256_set1_epi16(isamp);
      __m256i mask = _mm256_cmpgt_epi16(samples[isamp], max);
      max = _mm256_blendv_epi8(max, samples[isamp], mask);
      imax = _mm256_blendv_epi8(imax, visamp, mask);

      __m256i samp_l = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(samples[isamp],0));
      __m256i samp_u = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(samples[isamp],1));

      win_l = _mm256_add_epi32(win_l, samp_l);
      win_u = _mm256_add_epi32(win_u, samp_u);

      visamp = _mm256_set1_epi32(isamp);
      sum_qt_l = _mm256_add_epi32(sum_qt_l, _mm256_mullo_epi32(visamp, samp_l));
      sum_qt_u = _mm256_add_epi32(sum_qt_u, _mm256_mullo_epi32(visamp, samp_u));
    }

    __m256i sig_max_l = win_l;
    __m256i sig_max_u = win_u;
    __m256i isig_max_l = _mm256_set1_epi32(0);
    __m256i isig_max_u = _mm256_set1_epi32(0);

    __m256i sig_l;
    __m256i sig_u;
    __m256i bkg_l;
    __m256i bkg_u;
    __m256i sum_q_l = win_l;
    __m256i sum_q_u = win_u;

    for(isamp=window_n;isamp<nsamp;++isamp)
    {
      int iss = isamp-16;

      if(bkg_window_0 == iss){
        bkg_l = win_l;
        bkg_u = win_u;
      }

      __m256i visamp = _mm256_set1_epi32(iss);
      __m256i mask = _mm256_cmpeq_epi32(*(__m256i*)(sig_window_0 + iblock*16), visamp);
      sig_l = _mm256_blendv_epi8(sig_l, win_l, mask);
      mask = _mm256_cmpeq_epi32(*(__m256i*)(sig_window_0 + iblock*16 + 8), visamp);
      sig_u = _mm256_blendv_epi8(sig_u, win_u, mask);

      visamp = _mm256_set1_epi16(isamp);
      mask = _mm256_cmpgt_epi16(samples[isamp], max);
      max = _mm256_blendv_epi8(max, samples[isamp], mask);
      imax = _mm256_blendv_epi8(imax, visamp, mask);

      __m256i samp_l = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(samples[isamp],0));
      __m256i samp_u = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(samples[isamp],1));

      win_l = _mm256_add_epi32(win_l, samp_l);
      win_u = _mm256_add_epi32(win_u, samp_u);

      sum_q_l = _mm256_add_epi32(sum_q_l, samp_l);
      sum_q_u = _mm256_add_epi32(sum_q_u, samp_u);

      visamp = _mm256_set1_epi32(isamp);
      sum_qt_l = _mm256_add_epi32(sum_qt_l, _mm256_mullo_epi32(visamp, samp_l));
      sum_qt_u = _mm256_add_epi32(sum_qt_u, _mm256_mullo_epi32(visamp, samp_u));

      samp_l = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(samples[iss],0));
      samp_u = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(samples[iss],1));

      win_l = _mm256_sub_epi32(win_l, samp_l);
      win_u = _mm256_sub_epi32(win_u, samp_u);

      visamp = _mm256_set1_epi32(iss+1);
      mask = _mm256_cmpgt_epi16(win_l, sig_max_l);
      sig_max_l = _mm256_blendv_epi8(sig_max_l, win_l, mask);
      isig_max_l = _mm256_blendv_epi8(isig_max_l, visamp, mask);
      mask = _mm256_cmpgt_epi16(win_u, sig_max_u);
      sig_max_u = _mm256_blendv_epi8(sig_max_u, win_u, mask);
      isig_max_u = _mm256_blendv_epi8(isig_max_u, visamp, mask);
    }

    _mm256_store_si256((__m256i*)chan_max_index, _mm256_cvtepu16_epi32(_mm256_extracti128_si256(imax,0)));
    chan_max_index+=8;
    _mm256_store_si256((__m256i*)chan_max_index, _mm256_cvtepu16_epi32(_mm256_extracti128_si256(imax,1)));
    chan_max_index+=8;

    _mm256_store_si256((__m256i*)chan_max, _mm256_cvtepu16_epi32(_mm256_extracti128_si256(max,0)));
    chan_max+=8;
    _mm256_store_si256((__m256i*)chan_max, _mm256_cvtepu16_epi32(_mm256_extracti128_si256(max,1)));
    chan_max+=8;

    _mm256_store_si256((__m256i*)chan_sig_win_sum, sig_l);
    chan_sig_win_sum+=8;
    _mm256_store_si256((__m256i*)chan_sig_win_sum, sig_u);
    chan_sig_win_sum+=8;

    _mm256_store_si256((__m256i*)chan_sig_max_sum, sig_max_l);
    chan_sig_max_sum+=8;
    _mm256_store_si256((__m256i*)chan_sig_max_sum, sig_max_u);
    chan_sig_max_sum+=8;

    _mm256_store_si256((__m256i*)chan_sig_max_sum_index, isig_max_l);
    chan_sig_max_sum_index+=8;
    _mm256_store_si256((__m256i*)chan_sig_max_sum_index, isig_max_u);
    chan_sig_max_sum_index+=8;

    _mm256_store_si256((__m256i*)chan_bkg_win_sum, bkg_l);
    chan_bkg_win_sum+=8;
    _mm256_store_si256((__m256i*)chan_bkg_win_sum, bkg_u);
    chan_bkg_win_sum+=8;

    __m256 fbkg_l = _mm256_cvtepi32_ps(bkg_l);
    __m256 fbkg_u = _mm256_cvtepi32_ps(bkg_u);

    __m256 ped_l = _mm256_load_ps(ped);
    __m256 ped_n = _mm256_mul_ps(ped_l, _mm256_set1_ps(ped_iir_old));
    ped_n = _mm256_fmadd_ps(fbkg_l, _mm256_set1_ps(ped_iir_new), ped_n);
    // This tricky blend uses the sign bit to
    ped_l = _mm256_blendv_ps(ped_n, fbkg_l, ped_l);
    _mm256_store_ps(ped, ped_l);
    ped+=8;

    __m256 ped_u = _mm256_load_ps(ped);
    ped_n = _mm256_mul_ps(ped_u, _mm256_set1_ps(ped_iir_old));
    ped_n = _mm256_fmadd_ps(fbkg_u, _mm256_set1_ps(ped_iir_new), ped_n);
    ped_u = _mm256_blendv_ps(ped_n, fbkg_u, ped_u);
    _mm256_store_ps(ped, ped_u);
    ped+=8;

    _mm256_store_si256((__m256i*)chan_all_sum_q, sum_q_l);
    chan_all_sum_q+=8;
    _mm256_store_si256((__m256i*)chan_all_sum_q, sum_q_u);
    chan_all_sum_q+=8;

    _mm256_store_si256((__m256i*)chan_all_sum_qt, sum_qt_l);
    chan_all_sum_qt+=8;
    _mm256_store_si256((__m256i*)chan_all_sum_qt, sum_qt_u);
    chan_all_sum_qt+=8;

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
  }
#else  // defined(__AVX2__) and defined(__FMA__)
  throw std::runtime_error("AVX2_SingleGainDualWindowWaveformTreatmentEventVisitor: AVX2 or FMA not available at compile time");
#endif // defined(__AVX2__)
}
