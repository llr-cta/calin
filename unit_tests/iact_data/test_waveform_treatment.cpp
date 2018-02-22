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
#include <provenance/system_info.hpp>

using calin::math::rng::NR3_AVX2_RNGCore;
using calin::util::memory::safe_aligned_calloc;
using namespace calin::math::simd;
using namespace calin::iact_data::waveform_treatment_event_visitor;
using namespace calin::ix::iact_data::waveform_treatment_event_visitor;

static constexpr unsigned NSIM_TRACEANAL = 10000;

constexpr int nchan = 1024;
constexpr int nsamp = 60;

#if defined(__AVX2__) and defined(__FMA__)

TEST(TestWaveformTreatment, TraceAnalysis)
{
  SingleGainDualWindowWaveformTreatmentEventVisitorConfig config =
    SingleGainDualWindowWaveformTreatmentEventVisitor::default_config();
  config.set_sig_integration_0(24);
  SingleGainDualWindowWaveformTreatmentEventVisitor wfev(config);

  calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration run_config;
  run_config.set_num_samples(nsamp);
  for(unsigned ichan=0;ichan<nchan;ichan++)run_config.add_configured_channel_id(ichan);
  wfev.visit_telescope_run(&run_config, nullptr);

  uint16_t* samples_data;
  safe_aligned_calloc(samples_data, nchan*nsamp);

  NR3_AVX2_RNGCore core(12345);
  for(unsigned iloop=0;iloop<NSIM_TRACEANAL;iloop++)
  {
    const __m256i mask_12bit = _mm256_set1_epi16((1<<12)-1);
    for(unsigned i=0;i<nchan*nsamp/16;i++) {
      __m256i x = _mm256_and_si256(core.uniform_uivec256(), mask_12bit);
      _mm256_store_si256((__m256i*)(samples_data+i*16), x);
    }

    wfev.scalar_analyze_waveforms(samples_data);

#if 0
    if(iloop==0) {
      for(unsigned ichan=0; ichan<32; ichan++) {
        std::cout << ichan << ' '
          << wfev.chan_max_index()[ichan] << ' ' << wfev.chan_max()[ichan] << ' '
          << wfev.chan_bkg_win_sum()[ichan] << ' ' << wfev.chan_sig_win_sum()[ichan] << ' '
          << wfev.chan_all_sum_q()[ichan] << ' ' << wfev.chan_all_sum_qt()[ichan] << ' '
          << wfev.chan_ped()[ichan] << ' ' << wfev.chan_sig()[ichan] << ' '
          << wfev.chan_mean_t()[ichan] << ' '
          << wfev.chan_sig_max_sum()[ichan] << ' ' << wfev.chan_sig_max_sum_index()[ichan] << ' '
          << '\n';
      }
    }
#endif
  }
  free(samples_data);
}

TEST(TestWaveformTreatment, AVX_TraceAnalysis)
{
  SingleGainDualWindowWaveformTreatmentEventVisitorConfig config =
    SingleGainDualWindowWaveformTreatmentEventVisitor::default_config();
  config.set_sig_integration_0(24);
  AVX2_SingleGainDualWindowWaveformTreatmentEventVisitor wfev(config);

  calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration run_config;
  run_config.set_num_samples(nsamp);
  for(unsigned ichan=0;ichan<nchan;ichan++)run_config.add_configured_channel_id(ichan);
  wfev.visit_telescope_run(&run_config, nullptr);

  uint16_t* samples_data;
  safe_aligned_calloc(samples_data, nchan*nsamp);

  NR3_AVX2_RNGCore core(12345);
  for(unsigned iloop=0;iloop<NSIM_TRACEANAL;iloop++)
  {
    const __m256i mask_12bit = _mm256_set1_epi16((1<<12)-1);
    for(unsigned i=0;i<nchan*nsamp/16;i++) {
      __m256i x = _mm256_and_si256(core.uniform_uivec256(), mask_12bit);
      _mm256_store_si256((__m256i*)(samples_data+i*16), x);
    }

    wfev.avx2_analyze_waveforms(samples_data);

#if 0
    if(iloop==0) {
      for(unsigned ichan=0; ichan<32; ichan++) {
        std::cout << ichan << ' '
          << wfev.chan_max_index()[ichan] << ' ' << wfev.chan_max()[ichan] << ' '
          << wfev.chan_bkg_win_sum()[ichan] << ' ' << wfev.chan_sig_win_sum()[ichan] << ' '
          << wfev.chan_all_sum_q()[ichan] << ' ' << wfev.chan_all_sum_qt()[ichan] << ' '
          << wfev.chan_ped()[ichan] << ' ' << wfev.chan_sig()[ichan] << ' '
          << wfev.chan_mean_t()[ichan] << ' '
          << wfev.chan_sig_max_sum()[ichan] << ' ' << wfev.chan_sig_max_sum_index()[ichan] << ' '
          << '\n';
      }
    }
#endif
  }
  free(samples_data);
}

TEST(TestWaveformTreatment, AVX_TraceAnalysis_V2)
{
  SingleGainDualWindowWaveformTreatmentEventVisitorConfig config =
    SingleGainDualWindowWaveformTreatmentEventVisitor::default_config();
  config.set_sig_integration_0(24);
  AVX2_SingleGainDualWindowWaveformTreatmentEventVisitor wfev(config);

  calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration run_config;
  run_config.set_num_samples(nsamp);
  for(unsigned ichan=0;ichan<nchan;ichan++)run_config.add_configured_channel_id(ichan);
  wfev.visit_telescope_run(&run_config, nullptr);

  uint16_t* samples_data;
  safe_aligned_calloc(samples_data, nchan*nsamp);

  NR3_AVX2_RNGCore core(12345);
  for(unsigned iloop=0;iloop<NSIM_TRACEANAL;iloop++)
  {
    const __m256i mask_12bit = _mm256_set1_epi16((1<<12)-1);
    for(unsigned i=0;i<nchan*nsamp/16;i++) {
      __m256i x = _mm256_and_si256(core.uniform_uivec256(), mask_12bit);
      _mm256_store_si256((__m256i*)(samples_data+i*16), x);
    }

    wfev.avx2_analyze_waveforms_v2(samples_data);

#if 0
    if(iloop==0) {
      for(unsigned ichan=0; ichan<32; ichan++) {
        std::cout << ichan << ' '
          << wfev.chan_max_index()[ichan] << ' ' << wfev.chan_max()[ichan] << ' '
          << wfev.chan_bkg_win_sum()[ichan] << ' ' << wfev.chan_sig_win_sum()[ichan] << ' '
          << wfev.chan_all_sum_q()[ichan] << ' ' << wfev.chan_all_sum_qt()[ichan] << ' '
          << wfev.chan_ped()[ichan] << ' ' << wfev.chan_sig()[ichan] << ' '
          << wfev.chan_mean_t()[ichan] << ' '
          << wfev.chan_sig_max_sum()[ichan] << ' ' << wfev.chan_sig_max_sum_index()[ichan] << ' '
          << '\n';
      }
    }
#endif
  }
  free(samples_data);
}


TEST(TestWaveformTreatment, AVX_TraceAnalysis_V3)
{
  SingleGainDualWindowWaveformTreatmentEventVisitorConfig config =
    SingleGainDualWindowWaveformTreatmentEventVisitor::default_config();
  config.set_sig_integration_0(24);
  AVX2_SingleGainDualWindowWaveformTreatmentEventVisitor wfev(config);

  calin::ix::iact_data::telescope_run_configuration::TelescopeRunConfiguration run_config;
  run_config.set_num_samples(nsamp);
  for(unsigned ichan=0;ichan<nchan;ichan++)run_config.add_configured_channel_id(ichan);
  wfev.visit_telescope_run(&run_config, nullptr);

  uint16_t* samples_data;
  safe_aligned_calloc(samples_data, nchan*nsamp);

  NR3_AVX2_RNGCore core(12345);
  for(unsigned iloop=0;iloop<NSIM_TRACEANAL;iloop++)
  {
    const __m256i mask_12bit = _mm256_set1_epi16((1<<12)-1);
    for(unsigned i=0;i<nchan*nsamp/16;i++) {
      __m256i x = _mm256_and_si256(core.uniform_uivec256(), mask_12bit);
      _mm256_store_si256((__m256i*)(samples_data+i*16), x);
    }

    wfev.avx2_analyze_waveforms_v3(samples_data);

#if 0
    if(iloop==0) {
      for(unsigned ichan=0; ichan<32; ichan++) {
        std::cout << ichan << ' '
          << wfev.chan_max_index()[ichan] << ' ' << wfev.chan_max()[ichan] << ' '
          << wfev.chan_bkg_win_sum()[ichan] << ' ' << wfev.chan_sig_win_sum()[ichan] << ' '
          << wfev.chan_all_sum_q()[ichan] << ' ' << wfev.chan_all_sum_qt()[ichan] << ' '
          << wfev.chan_ped()[ichan] << ' ' << wfev.chan_sig()[ichan] << ' '
          << wfev.chan_mean_t()[ichan] << ' '
          << wfev.chan_sig_max_sum()[ichan] << ' ' << wfev.chan_sig_max_sum_index()[ichan] << ' '
          << '\n';
      }
    }
#endif
  }
  free(samples_data);
}

#endif // defined(__AVX2__) and defined(__FMA__)

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
