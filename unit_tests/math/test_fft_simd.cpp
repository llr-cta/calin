/*

   calin/unit_tests/math/test_simd.cpp -- Stephen Fegan -- 2018-01-05

   Unit tests for SIMD FFT class

   Copyright 2018, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, CNRS, Ecole Polytechnique, Institut Polytechnique de Paris

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

#include <math/rng.hpp>
#include <math/simd.hpp>
#include <math/special.hpp>
#include <util/memory.hpp>
#include <math/fft_simd.hpp>

using calin::math::rng::NR3_AVX2_RNGCore;
using calin::math::rng::NR3RNGCore;
using calin::math::special::SQR;
using namespace calin::math::simd;
using namespace calin::math::fft_simd;

#ifdef CALIN_HAS_NR3_AVX2_RNGCORE

static constexpr unsigned NSIM_COMPARE = 256;
static constexpr unsigned NSIM_SPEED = 1024*1024;

class Compare : public ::testing::TestWithParam<unsigned>
{
  // nothing to see here
};

TEST_P(Compare, AVX_R2C_Equals_FFT_Float)
{
  unsigned n = GetParam();
  NR3_AVX2_RNGCore core(12345);
  auto* codelet = new_m256_codelet_r2c_dft(n);
  auto* fftw = new_m256_fftw_r2c_dft(n);
  __m256* xt = codelet->alloc_real_array();
  __m256* xf1 = codelet->alloc_complex_array();
  __m256* xf2 = fftw->alloc_complex_array();
  for(unsigned isim=0; isim<NSIM_COMPARE; isim++) {
    for(unsigned i=0;i<n;i++)xt[i] = core.uniform_zc_m256();
    codelet->r2c(xt, xf1);
    fftw->r2c(xt, xf2);
    for(unsigned i=0;i<codelet->complex_array_size();i++)
      for(unsigned j=0;j<8;j++)
        EXPECT_NEAR(xf1[i][j], xf2[i][j], 1.5e-6)
          << "With i=" << i << " and j=" << j;
  }
  free(xf1);
  free(xf2);
  free(xt);
  delete codelet;
  delete fftw;
}

TEST_P(Compare, AVX_R2HC_Equals_FFT_Float)
{
  unsigned n = GetParam();
  NR3_AVX2_RNGCore core(12345);
  auto* codelet = new_m256_codelet_r2hc_dft(n);
  auto* fftw = new_m256_fftw_r2hc_dft(n);
  __m256* xt = codelet->alloc_real_array();
  __m256* xf1 = codelet->alloc_half_complex_array();
  __m256* xf2 = fftw->alloc_half_complex_array();
  for(unsigned isim=0; isim<NSIM_COMPARE; isim++) {
    for(unsigned i=0;i<n;i++)xt[i] = core.uniform_zc_m256();
    codelet->r2hc(xt, xf1);
    fftw->r2hc(xt, xf2);
    for(unsigned i=0;i<codelet->half_complex_array_size();i++)
      for(unsigned j=0;j<8;j++)
        EXPECT_NEAR(xf1[i][j], xf2[i][j], 1.5e-6)
          << "With i=" << i << " and j=" << j;
  }
  free(xf1);
  free(xf2);
  free(xt);
  delete codelet;
  delete fftw;
}

TEST_P(Compare, AVX_C2R_Equals_FFT_Float)
{
  unsigned n = GetParam();
  NR3_AVX2_RNGCore core(12345);
  auto* codelet = new_m256_codelet_r2c_dft(n);
  auto* fftw = new_m256_fftw_r2c_dft(n);
  __m256* xf1 = codelet->alloc_complex_array();
  __m256* xf2 = fftw->alloc_complex_array();
  __m256* xt1 = codelet->alloc_real_array();
  __m256* xt2 = fftw->alloc_real_array();
  for(unsigned isim=0; isim<NSIM_COMPARE; isim++) {
    for(unsigned i=0;i<codelet->complex_array_size();i++)
      xf1[i] = xf2[i] = core.uniform_zc_m256();
    codelet->c2r(xt1, xf1);
    fftw->c2r(xt2, xf2);
    for(unsigned i=0;i<n;i++)
      for(unsigned j=0;j<8;j++)
        EXPECT_NEAR(xt1[i][j], xt2[i][j], 3e-6)
          << "With i=" << i << " and j=" << j;
  }
  free(xf1);
  free(xf2);
  free(xt1);
  free(xt2);
  delete codelet;
  delete fftw;
}

TEST_P(Compare, AVX_HC2R_Equals_FFT_Float)
{
  unsigned n = GetParam();
  NR3_AVX2_RNGCore core(12345);
  auto* codelet = new_m256_codelet_r2hc_dft(n);
  auto* fftw = new_m256_fftw_r2hc_dft(n);
  __m256* xf1 = codelet->alloc_half_complex_array();
  __m256* xf2 = fftw->alloc_half_complex_array();
  __m256* xt1 = codelet->alloc_real_array();
  __m256* xt2 = fftw->alloc_real_array();
  for(unsigned isim=0; isim<NSIM_COMPARE; isim++) {
    for(unsigned i=0;i<n;i++)
      xf1[i] = xf2[i] = core.uniform_zc_m256();
    codelet->hc2r(xt1, xf1);
    fftw->hc2r(xt2, xf2);
    for(unsigned i=0;i<n;i++)
      for(unsigned j=0;j<8;j++)
        EXPECT_NEAR(xt1[i][j], xt2[i][j], 3e-6)
          << "With i=" << i << " and j=" << j;
  }
  free(xf1);
  free(xf2);
  free(xt1);
  free(xt2);
  delete codelet;
  delete fftw;
}

class Speed_1M : public ::testing::TestWithParam<unsigned>
{
  // nothing to see here
};

TEST_P(Speed_1M, FFTW_R2HC)
{
  unsigned n = GetParam();
  NR3_AVX2_RNGCore core(12345);
  auto* fftw = new_m256_fftw_r2hc_dft(n);
  __m256* xt = fftw->alloc_real_array();
  __m256* xf = fftw->alloc_half_complex_array();
  for(unsigned i=0;i<n;i++)xt[i] = core.uniform_zc_m256();
  for(unsigned isim=0; isim<NSIM_SPEED; isim++) {
    fftw->r2hc(xt, xf);
  }
  free(xf);
  free(xt);
  delete fftw;
}

TEST_P(Speed_1M, AVXCodelet_R2HC)
{
  unsigned n = GetParam();
  NR3_AVX2_RNGCore core(12345);
  auto* codelet = new_m256_codelet_r2hc_dft(n);
  __m256* xt = codelet->alloc_real_array();
  __m256* xf = codelet->alloc_half_complex_array();
  for(unsigned i=0;i<n;i++)xt[i] = core.uniform_zc_m256();
  for(unsigned isim=0; isim<NSIM_SPEED; isim++) {
    codelet->r2hc(xt, xf);
  }
  free(xf);
  free(xt);
  delete codelet;
}

INSTANTIATE_TEST_CASE_P(TestFFTSIMD,
                        Compare,
                        ::testing::ValuesIn(list_available_m256_codelets()),
                        ::testing::PrintToStringParamName());

INSTANTIATE_TEST_CASE_P(TestFFTSIMD,
                        Speed_1M,
                        ::testing::Values(16,60,64),
                        ::testing::PrintToStringParamName());

#endif // defined CALIN_HAS_NR3_AVX2_RNGCORE

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
