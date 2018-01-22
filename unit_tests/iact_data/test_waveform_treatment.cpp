/*

   calin/unit_tests/math/test_waveform_treatment.cpp -- Stephen Fegan -- 2018-01-19

   Unit tests for waveform treatment

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

#include <iostream>
#include <iomanip>
#include <gtest/gtest.h>
#include <vector>
#include <tuple>
#include <cmath>

#include <util/memory.hpp>
#include <math/rng.hpp>
#include <math/simd.hpp>
#include <iact_data/waveform_treatment_event_visitor.hpp>
#include <iact_data/waveform_treatment_event_visitor_impl.hpp>
#include <provenance/system_info.hpp>

using calin::math::rng::NR3_AVX2_RNGCore;
using calin::util::memory::safe_aligned_calloc;
using namespace calin::math::simd;
using namespace calin::iact_data::waveform_treatment_event_visitor;

static constexpr unsigned NSIM_TRACEANAL = 10000;

uint16_t* samples_data = nullptr;
int* sig_window_0_ = nullptr;
float* chan_ped_est_ = nullptr;
int* chan_max_ = nullptr;
int* chan_max_index_ = nullptr;
int* chan_bkg_win_sum_ = nullptr;
int* chan_sig_win_sum_ = nullptr;
int* chan_sig_max_sum_ = nullptr;
int* chan_sig_max_sum_index_ = nullptr;
int* chan_all_sum_q_ = nullptr;
int* chan_all_sum_qt_ = nullptr;
float* chan_sig_ = nullptr;
float* chan_mean_t_ = nullptr;

namespace {

void free_all_arrays()
{
  free(samples_data);
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

void alloc_all_arrays(unsigned nchan_, unsigned nsamp_)
{
  auto* host_info = calin::provenance::system_info::the_host_info();
  safe_aligned_calloc(samples_data, nchan_*nsamp_, host_info->log2_simd_vec_size());
  safe_aligned_calloc(sig_window_0_, nchan_, host_info->log2_simd_vec_size());
  safe_aligned_calloc(chan_ped_est_, nchan_, host_info->log2_simd_vec_size());
  safe_aligned_calloc(chan_max_index_, nchan_, host_info->log2_simd_vec_size());
  safe_aligned_calloc(chan_max_, nchan_, host_info->log2_simd_vec_size());
  safe_aligned_calloc(chan_bkg_win_sum_, nchan_, host_info->log2_simd_vec_size());
  safe_aligned_calloc(chan_sig_win_sum_, nchan_, host_info->log2_simd_vec_size());
  safe_aligned_calloc(chan_sig_max_sum_, nchan_, host_info->log2_simd_vec_size());
  safe_aligned_calloc(chan_sig_max_sum_index_, nchan_, host_info->log2_simd_vec_size());
  safe_aligned_calloc(chan_all_sum_q_, nchan_, host_info->log2_simd_vec_size());
  safe_aligned_calloc(chan_all_sum_qt_, nchan_, host_info->log2_simd_vec_size());
  safe_aligned_calloc(chan_sig_, nchan_, host_info->log2_simd_vec_size());
  safe_aligned_calloc(chan_mean_t_, nchan_, host_info->log2_simd_vec_size());
}

}

constexpr int nchan = 1024;
constexpr int nsamp = 60;

#if defined(__AVX2__) and defined(__FMA__)

TEST(TestWaveformTreatment, TraceAnalysis)
{
  alloc_all_arrays(nchan, nsamp);
  std::fill(sig_window_0_,sig_window_0_+nsamp,24);
  NR3_AVX2_RNGCore core(12345);
  for(unsigned iloop=0;iloop<NSIM_TRACEANAL;iloop++)
  {
    const __m256i mask_12bit = _mm256_set1_epi16((1<<12)-1);
    for(unsigned i=0;i<nchan*nsamp/16;i++) {
      __m256i x = _mm256_and_si256(core.uniform_uivec256(), mask_12bit);
      _mm256_store_si256((__m256i*)(samples_data+i*16), x);
    }

    SingleGainDualWindowWaveformTreatmentEventVisitor::
    analyze_waveforms(samples_data, nchan, nsamp, 16, 0, sig_window_0_,
      chan_ped_est_, 0.0, 1.0,
      chan_max_index_, chan_max_,
      chan_bkg_win_sum_, chan_sig_win_sum_,
      chan_sig_max_sum_, chan_sig_max_sum_index_,
      chan_all_sum_q_, chan_all_sum_qt_,
      chan_sig_, chan_mean_t_);

#if 0
    if(iloop==0) {
      for(unsigned ichan=0; ichan<32; ichan++) {
        std::cout << ichan << ' ' << all_imax[ichan] << ' ' << all_max[ichan] << ' '
          << all_bkg[ichan] << ' ' << all_sig[ichan] << ' '
          << all_sum_q[ichan] << ' ' << all_sum_qt[ichan] << ' '
          << all_mean_t[ichan] << '\n';
#endif
  }
  free_all_arrays();
}

TEST(TestWaveformTreatment, AVX_TraceAnalysis)
{
  alloc_all_arrays(nchan, nsamp);
  std::fill(sig_window_0_,sig_window_0_+nsamp,24);
  NR3_AVX2_RNGCore core(12345);
  for(unsigned iloop=0;iloop<NSIM_TRACEANAL;iloop++)
  {
    const __m256i mask_12bit = _mm256_set1_epi16((1<<12)-1);
    for(unsigned i=0;i<nchan*nsamp/16;i++) {
      __m256i x = _mm256_and_si256(core.uniform_uivec256(), mask_12bit);
      _mm256_store_si256((__m256i*)(samples_data+i*16), x);
    }

    AVX2_SingleGainDualWindowWaveformTreatmentEventVisitor::
    avx2_analyze_waveforms(samples_data, nchan, nsamp, 16, 0, sig_window_0_,
      chan_ped_est_, 0.0, 1.0,
      chan_max_index_, chan_max_,
      chan_bkg_win_sum_, chan_sig_win_sum_,
      chan_sig_max_sum_, chan_sig_max_sum_index_,
      chan_all_sum_q_, chan_all_sum_qt_,
      chan_sig_, chan_mean_t_);

#if 0
    if(iloop==0) {
      for(unsigned ichan=0; ichan<32; ichan++) {
        std::cout << ichan << ' ' << all_imax[ichan] << ' ' << all_max[ichan] << ' '
          << all_bkg[ichan] << ' ' << all_sig[ichan] << ' '
          << all_sum_q[ichan] << ' ' << all_sum_qt[ichan] << ' '
          << all_mean_t[ichan] << '\n';
#endif
  }
  free_all_arrays();
}

#endif // defined(__AVX2__) and defined(__FMA__)

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
