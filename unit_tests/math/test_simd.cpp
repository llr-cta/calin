/*

   calin/unit_tests/math/test_simd.cpp -- Stephen Fegan -- 2018-01-05

   Unit tests for SIMD

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

#include <fftw3.h>

#include <math/rng.hpp>
#include <math/simd.hpp>
#include <math/special.hpp>
#include <util/memory.hpp>

using calin::math::rng::NR3_AVX2_RNGCore;
using calin::math::rng::NR3RNGCore;
using calin::math::special::SQR;
using namespace calin::math::simd;

static constexpr uint64_t NSIM_RANDSINCOS = 100000000ULL;
#if 0
#ifdef CALIN_HAS_NR3_AVX2_RNGCORE
static constexpr unsigned NSIM_TRACEPSD = 4096;
#endif
#endif

TEST(TestSIMD, SpeedTest100M_Random64SinCos32)
{
  float sum_s = 0;
  float sum_c = 0;
  uint64_t seed = calin::math::rng::RNG::uint64_from_random_device();
  NR3RNGCore core(seed);
  for(unsigned i=0;i<NSIM_RANDSINCOS/2;i++) {
    uint64_t x64 = core.uniform_uint64();
    float x = 2.328306437e-10 * 2.0 * M_PI * unsigned(x64);
    sum_s += ::sinf(x);
    sum_c += ::cosf(x);
    x = 2.328306437e-10 * 2.0 * M_PI * unsigned(x64>>32);
    sum_s += ::sinf(x);
    sum_c += ::cosf(x);
  }

  EXPECT_GE(sum_s, -float(NSIM_RANDSINCOS));
  EXPECT_LE(sum_s, float(NSIM_RANDSINCOS));
  EXPECT_GE(sum_c, -float(NSIM_RANDSINCOS));
  EXPECT_LE(sum_c, float(NSIM_RANDSINCOS));
}

#ifdef CALIN_HAS_NR3_AVX2_RNGCORE

TEST(TestSIMD, SpeedTest100M_Random256SinCos32)
{
  float sum_s = 0;
  float sum_c = 0;
  uint64_t seed = calin::math::rng::RNG::uint64_from_random_device();
  NR3_AVX2_RNGCore core(seed);
  for(unsigned i=0;i<NSIM_RANDSINCOS/8;i++) {
    float x[8] __attribute__((aligned(32)));
    __m256 vx = core.uniform_zc_m256(2*M_PI);
    _mm256_store_ps(x, vx);
    sum_s += ::sinf(x[0]);
    sum_c += ::cosf(x[0]);
    sum_s += ::sinf(x[1]);
    sum_c += ::cosf(x[1]);
    sum_s += ::sinf(x[2]);
    sum_c += ::cosf(x[2]);
    sum_s += ::sinf(x[3]);
    sum_c += ::cosf(x[3]);
    sum_s += ::sinf(x[4]);
    sum_c += ::cosf(x[4]);
    sum_s += ::sinf(x[5]);
    sum_c += ::cosf(x[5]);
    sum_s += ::sinf(x[6]);
    sum_c += ::cosf(x[6]);
    sum_s += ::sinf(x[7]);
    sum_c += ::cosf(x[7]);
  }

  EXPECT_GE(sum_s, -float(NSIM_RANDSINCOS));
  EXPECT_LE(sum_s, float(NSIM_RANDSINCOS));
  EXPECT_GE(sum_c, -float(NSIM_RANDSINCOS));
  EXPECT_LE(sum_c, float(NSIM_RANDSINCOS));
}

TEST(TestSIMD, SpeedTest100M_Random256SinCos256)
{
  __m256 sum_s = _mm256_setzero_ps();
  __m256 sum_c = _mm256_setzero_ps();
  uint64_t seed = calin::math::rng::RNG::uint64_from_random_device();
  NR3_AVX2_RNGCore core(seed);
  for(unsigned i=0;i<NSIM_RANDSINCOS/8;i++) {
    __m256 x = core.uniform_zc_m256(8.0);
    __m256 s;
    __m256 c;
    avx2_sincosf_domain_pi_poly3(x, s, c);
    sum_s = _mm256_add_ps(sum_s, s);
    sum_c = _mm256_add_ps(sum_c, c);
  }
  EXPECT_GE(reinterpret_cast<float*>(&sum_s)[0], -float(NSIM_RANDSINCOS));
  EXPECT_LE(reinterpret_cast<float*>(&sum_s)[0], float(NSIM_RANDSINCOS));
  EXPECT_GE(reinterpret_cast<float*>(&sum_c)[0], -float(NSIM_RANDSINCOS));
  EXPECT_LE(reinterpret_cast<float*>(&sum_c)[0], float(NSIM_RANDSINCOS));
}

TEST(TestSIMD, TransposeM256)
{
  __m256 a = _mm256_set_ps(81,71,61,51,41,31,21,11);
  __m256 b = _mm256_set_ps(82,72,62,52,42,32,22,12);
  __m256 c = _mm256_set_ps(83,73,63,53,43,33,23,13);
  __m256 d = _mm256_set_ps(84,74,64,54,44,34,24,14);
  __m256 e = _mm256_set_ps(85,75,65,55,45,35,25,15);
  __m256 f = _mm256_set_ps(86,76,66,56,46,36,26,16);
  __m256 g = _mm256_set_ps(87,77,67,57,47,37,27,17);
  __m256 h = _mm256_set_ps(88,78,68,58,48,38,28,18);
  calin::math::simd::transpose(a,b,c,d,e,f,g,h);
  for(unsigned i=0;i<8;i++)EXPECT_EQ(a[i], float(11+i));
  for(unsigned i=0;i<8;i++)EXPECT_EQ(b[i], float(21+i));
  for(unsigned i=0;i<8;i++)EXPECT_EQ(c[i], float(31+i));
  for(unsigned i=0;i<8;i++)EXPECT_EQ(d[i], float(41+i));
  for(unsigned i=0;i<8;i++)EXPECT_EQ(e[i], float(51+i));
  for(unsigned i=0;i<8;i++)EXPECT_EQ(f[i], float(61+i));
  for(unsigned i=0;i<8;i++)EXPECT_EQ(g[i], float(71+i));
  for(unsigned i=0;i<8;i++)EXPECT_EQ(h[i], float(81+i));
}

TEST(TestSIMD, TransposeM256D)
{
  __m256d a = _mm256_set_pd(41,31,21,11);
  __m256d b = _mm256_set_pd(42,32,22,12);
  __m256d c = _mm256_set_pd(43,33,23,13);
  __m256d d = _mm256_set_pd(44,34,24,14);
  calin::math::simd::transpose(a,b,c,d);
  for(unsigned i=0;i<4;i++)EXPECT_EQ(a[i], double(11+i));
  for(unsigned i=0;i<4;i++)EXPECT_EQ(b[i], double(21+i));
  for(unsigned i=0;i<4;i++)EXPECT_EQ(c[i], double(31+i));
  for(unsigned i=0;i<4;i++)EXPECT_EQ(d[i], double(41+i));
}


#if 0
TEST(TestTracePSD, Scalar)
{
  NR3_AVX2_RNGCore core(12345);
  uint16_t* hg = new uint16_t[nchan*nsamp];
  float* xt = fftwf_alloc_real(nsamp);
  float* xf = fftwf_alloc_real(2*(nsamp/2 + 1));
  double* psd_sum = fftw_alloc_real(nchan*(nsamp/2 + 1));
  double* psd_sumsq = fftw_alloc_real(nchan*(nsamp/2 + 1));
  std::fill(psd_sum, psd_sum+nsamp/2+1, 0.0);
  std::fill(psd_sumsq, psd_sum+nsamp/2+1, 0.0);
  fftwf_plan plan = fftwf_plan_dft_r2c_1d(nsamp, xt, (fftwf_complex*)xf, FFTW_MEASURE);
  for(unsigned iloop=0;iloop<NSIM_TRACEPSD;iloop++)
  {
    const __m256i mask_12bit = _mm256_set1_epi16((1<<12)-1);
    for(unsigned i=0;i<nchan*nsamp/16;i++) {
      __m256i x = _mm256_and_si256(core.uniform_m256i(), mask_12bit);
      _mm256_storeu_si256((__m256i*)(hg+i*16), x);
    }

    const uint16_t*__restrict__ samp = hg;
    double*__restrict__ ipsd_sum = psd_sum;
    double*__restrict__ ipsd_sumsq = psd_sumsq;
    for(unsigned ichan=0;ichan<nchan;ichan++) {
      float*__restrict__ ixt = xt;
      for(unsigned isamp=0;isamp<nsamp;isamp++) {
        *ixt++ = float(*samp++);
      }

      fftwf_execute(plan);

      const float*__restrict__ ri = xf;
      const float*__restrict__ ci = xf + nsamp - 1;
      double psdi = SQR(*ri++);
      *ipsd_sum++ += psdi;
      *ipsd_sumsq++ += SQR(psdi);
      while(ri < ci)
      {
        psdi = SQR(*ri++) + SQR(*ci--);
        *ipsd_sum++ += psdi;
        *ipsd_sumsq++ += SQR(psdi);
      }
      if(ri==ci)
      {
        double psdi = SQR(*ri);
        *ipsd_sum++ += psdi;
        *ipsd_sumsq++ += SQR(psdi);
      }
    }
  }

  for(unsigned ifreq=0;ifreq<nsamp/2+1;ifreq++)
    std::cout << psd_sum[ifreq] << ' ';
  std::cout << '\n';

  fftw_free(psd_sum);
  fftw_free(psd_sumsq);
  fftwf_destroy_plan(plan);
  fftwf_free(xf);
  fftwf_free(xt);
  delete[] hg;
}

TEST(TestTracePSD, AVX2)
{
  NR3_AVX2_RNGCore core(12345);
  uint16_t* hg;
  calin::util::memory::safe_aligned_calloc(hg, nchan*nsamp);
  std::fill(cov, cov+nchan*nsamp*(nsamp+1)/2, 0);
  const unsigned nv_trace = (nsamp+15)/16;
  const unsigned nv_block = nv_trace*16;
  __m256i* samples;
  calin::util::memory::safe_aligned_calloc(samples, nv_block);

  __m256* xt;
  __m256* xf;
  calin::util::memory::safe_aligned_calloc(xt, 2*nsamp);
  calin::util::memory::safe_aligned_calloc(xf, 2*(nsamp/2 + 1));

  const unsigned nblock = nchan/16;

  __m256d* psd_sum;
  __m256d* psd_sumsq;
  calin::util::memory::safe_aligned_calloc(psd_sum, 4*nblock*(nsamp/2 + 1));
  calin::util::memory::safe_aligned_calloc(psd_sumsq, 4*nblock*(nsamp/2 + 1));
  std::fill(psd_sum, psd_sum+4*nblock*(nsamp/2 + 1), _mm256_setzero_pd());
  std::fill(psd_sumsq, psd_sumsq+4*nblock*(nsamp/2 + 1), _mm256_setzero_pd());

  int n = nsamp;
  fftwf_plan plan = fftwf_plan_many_dft_r2c(1, &n, 16,
                            (float*)xt, nullptr, 16, 1,
                            (fftwf_complex*)xf, nullptr, 16, 1,
                            FFTW_MEASURE);

  for(unsigned iloop=0;iloop<NSIM_TRACEPSD;iloop++)
  {
    const __m256i mask_12bit = _mm256_set1_epi16((1<<12)-1);
    for(unsigned i=0;i<nchan*nsamp/16;i++) {
      __m256i x = _mm256_and_si256(core.uniform_m256i(), mask_12bit);
      _mm256_storeu_si256((__m256i*)(hg+i*16), x);
    }

    __m256d*__restrict__ ipsd_sum = psd_sum;
    __m256d*__restrict__ ipsd_sumsq = psd_sumsq;

    for(unsigned iblock=0;iblock<nblock;iblock++)
    {
      uint16_t* base = hg + iblock*nsamp*16;
      __m256i*__restrict__ vp = samples;
      for(unsigned iv_trace=0; iv_trace<nv_trace; iv_trace++) {
        for(unsigned ivec=0; ivec<16; ivec++) {
          *(vp++) = _mm256_loadu_si256((__m256i*)(base + iv_trace*16 + nsamp*ivec));
        }
        calin::math::simd::avx2_m256_swizzle_u16(samples + iv_trace*16);
      }

      vp = samples;
      __m256*__restrict__ ixt = xt;
      for(unsigned isamp=0;isamp<nsamp;isamp++) {
        __m256i samp = *vp++;
        *ixt++ = _mm256_cvtepi32_ps(_mm256_cvtepu16_epi32(_mm256_extracti128_si256(samp,0)));
        *ixt++ = _mm256_cvtepi32_ps(_mm256_cvtepu16_epi32(_mm256_extracti128_si256(samp,1)));
      }

      fftwf_execute(plan);

      const __m256*__restrict__ ri = xf;
      const __m256*__restrict__ ci = ri + 2*(nsamp - 1);

      __m256 ri_f;
      __m256 ci_f;

      __m256d val_d;

      ri_f = *ri++;
      val_d = _mm256_cvtps_pd(_mm256_extractf128_ps(ri_f, 0));
      val_d = _mm256_mul_pd(val_d, val_d);
      *ipsd_sum = _mm256_add_pd(val_d, *ipsd_sum);
      ipsd_sum++;
      *ipsd_sumsq = _mm256_fmadd_pd(val_d, val_d, *ipsd_sumsq);
      ipsd_sumsq++;

      val_d = _mm256_cvtps_pd(_mm256_extractf128_ps(ri_f, 1));
      val_d = _mm256_mul_pd(val_d, val_d);
      *ipsd_sum = _mm256_add_pd(val_d, *ipsd_sum);
      ipsd_sum++;
      *ipsd_sumsq = _mm256_fmadd_pd(val_d, val_d, *ipsd_sumsq);
      ipsd_sumsq++;

      ri_f = *ri++;
      val_d = _mm256_cvtps_pd(_mm256_extractf128_ps(ri_f, 0));
      val_d = _mm256_mul_pd(val_d, val_d);
      *ipsd_sum = _mm256_add_pd(val_d, *ipsd_sum);
      ipsd_sum++;
      *ipsd_sumsq = _mm256_fmadd_pd(val_d, val_d, *ipsd_sumsq);
      ipsd_sumsq++;

      val_d = _mm256_cvtps_pd(_mm256_extractf128_ps(ri_f, 1));
      val_d = _mm256_mul_pd(val_d, val_d);
      *ipsd_sum = _mm256_add_pd(val_d, *ipsd_sum);
      ipsd_sum++;
      *ipsd_sumsq = _mm256_fmadd_pd(val_d, val_d, *ipsd_sumsq);
      ipsd_sumsq++;

      while(ri < ci)
      {
        __m256d cval_d;

        ri_f = *ri++;
        ci_f = *ci++;

        val_d = _mm256_cvtps_pd(_mm256_extractf128_ps(ri_f, 0));
        val_d = _mm256_mul_pd(val_d, val_d);
        cval_d = _mm256_cvtps_pd(_mm256_extractf128_ps(ci_f, 0));
        val_d = _mm256_fmadd_pd(cval_d, cval_d, val_d);
        *ipsd_sum = _mm256_add_pd(val_d, *ipsd_sum);
        ipsd_sum++;
        *ipsd_sumsq = _mm256_fmadd_pd(val_d, val_d, *ipsd_sumsq);
        ipsd_sumsq++;

        val_d = _mm256_cvtps_pd(_mm256_extractf128_ps(ri_f, 1));
        val_d = _mm256_mul_pd(val_d, val_d);
        cval_d = _mm256_cvtps_pd(_mm256_extractf128_ps(ci_f, 1));
        val_d = _mm256_fmadd_pd(cval_d, cval_d, val_d);
        *ipsd_sum = _mm256_add_pd(val_d, *ipsd_sum);
        ipsd_sum++;
        *ipsd_sumsq = _mm256_fmadd_pd(val_d, val_d, *ipsd_sumsq);
        ipsd_sumsq++;

        ri_f = *ri++;
        ci_f = *ci;
        ci -= 3;

        val_d = _mm256_cvtps_pd(_mm256_extractf128_ps(ri_f, 0));
        val_d = _mm256_mul_pd(val_d, val_d);
        cval_d = _mm256_cvtps_pd(_mm256_extractf128_ps(ci_f, 0));
        val_d = _mm256_fmadd_pd(cval_d, cval_d, val_d);
        *ipsd_sum = _mm256_add_pd(val_d, *ipsd_sum);
        ipsd_sum++;
        *ipsd_sumsq = _mm256_fmadd_pd(val_d, val_d, *ipsd_sumsq);
        ipsd_sumsq++;

        val_d = _mm256_cvtps_pd(_mm256_extractf128_ps(ri_f, 1));
        val_d = _mm256_mul_pd(val_d, val_d);
        cval_d = _mm256_cvtps_pd(_mm256_extractf128_ps(ci_f, 1));
        val_d = _mm256_fmadd_pd(cval_d, cval_d, val_d);
        *ipsd_sum = _mm256_add_pd(val_d, *ipsd_sum);
        ipsd_sum++;
        *ipsd_sumsq = _mm256_fmadd_pd(val_d, val_d, *ipsd_sumsq);
        ipsd_sumsq++;
      }
      if(ri==ci)
      {
        ri_f = *ri++;
        val_d = _mm256_cvtps_pd(_mm256_extractf128_ps(ri_f, 0));
        val_d = _mm256_mul_pd(val_d, val_d);
        *ipsd_sum = _mm256_add_pd(val_d, *ipsd_sum);
        ipsd_sum++;
        *ipsd_sumsq = _mm256_fmadd_pd(val_d, val_d, *ipsd_sumsq);
        ipsd_sumsq++;

        val_d = _mm256_cvtps_pd(_mm256_extractf128_ps(ri_f, 1));
        val_d = _mm256_mul_pd(val_d, val_d);
        *ipsd_sum = _mm256_add_pd(val_d, *ipsd_sum);
        ipsd_sum++;
        *ipsd_sumsq = _mm256_fmadd_pd(val_d, val_d, *ipsd_sumsq);
        ipsd_sumsq++;

        ri_f = *ri++;
        val_d = _mm256_cvtps_pd(_mm256_extractf128_ps(ri_f, 0));
        val_d = _mm256_mul_pd(val_d, val_d);
        *ipsd_sum = _mm256_add_pd(val_d, *ipsd_sum);
        ipsd_sum++;
        *ipsd_sumsq = _mm256_fmadd_pd(val_d, val_d, *ipsd_sumsq);
        ipsd_sumsq++;

        val_d = _mm256_cvtps_pd(_mm256_extractf128_ps(ri_f, 1));
        val_d = _mm256_mul_pd(val_d, val_d);
        *ipsd_sum = _mm256_add_pd(val_d, *ipsd_sum);
        ipsd_sum++;
        *ipsd_sumsq = _mm256_fmadd_pd(val_d, val_d, *ipsd_sumsq);
        ipsd_sumsq++;
      }
    }
  }

  for(unsigned ifreq=0;ifreq<nsamp/2+1;ifreq++)
    std::cout << psd_sum[ifreq*4][0] << ' ';
  std::cout << '\n';

  free(psd_sum);
  free(psd_sumsq);
  fftwf_destroy_plan(plan);
  free(xf);
  free(xt);
  free(hg);
  free(samples);
}

#endif // defined CALIN_HAS_NR3_AVX2_RNGCORE
#endif

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
