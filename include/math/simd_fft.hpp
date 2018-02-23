/*

   calin/math/simd_fft.hpp -- Stephen Fegan -- 2018-02-21

   SIMD FFT functions using codelets from FFTW/genfft

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

#include <vector>
#include <fftw3.h>

#if defined(__AVX__)
#include <immintrin.h>
#endif // defined(__AVX__)

#include <util/memory.hpp>

namespace calin { namespace math { namespace simd_fft {

template<typename T> void shuffle_rc(T& r, T&c) { }

#if defined(__AVX__)
void shuffle_c_v2s(__m256& r, __m256&c)
{
  __m256 t;
  t = _mm256_set_m128(_mm256_extractf128_ps(r,0), _mm256_extractf128_ps(c,0));
  c = _mm256_set_m128(_mm256_extractf128_ps(r,1), _mm256_extractf128_ps(c,1));
  c = _mm256_permute_ps(c, 0b10110001);
  r = _mm256_blend_ps(t, c, 0b10101010);
  c = _mm256_blend_ps(t, c, 0x01010101);
  r = _mm256_permute_ps(r, 0b11011000);
  c = _mm256_permute_ps(c, 0b10001101);
}

void shuffle_c_s2v(__m256& r, __m256&c)
{
  __m256 t;
  r = _mm256_permute_ps(r, 0b11011000);
  c = _mm256_permute_ps(c, 0b01110010);
  t = _mm256_blend_ps(r, c, 0b10101010);
  c = _mm256_blend_ps(r, c, 0x01010101);
  c = _mm256_permute_ps(c, 0b10110001);
  r = _mm256_set_m128(_mm256_extractf128_ps(t,0), _mm256_extractf128_ps(c,0));
  c = _mm256_set_m128(_mm256_extractf128_ps(t,1), _mm256_extractf128_ps(c,1));
}

void shuffle_c_v2s(__m256d& r, __m256d&c)
{
  __m256d t;
  t = _mm256_set_m128d(_mm256_extractf128_pd(r,0), _mm256_extractf128_pd(c,0));
  c = _mm256_set_m128d(_mm256_extractf128_pd(r,1), _mm256_extractf128_pd(c,1));
  r = _mm256_unpacklo_pd(t, c);
  c = _mm256_unpackhi_pd(t, c);
}

void shuffle_c_s2v(__m256d& r, __m256d&c)
{
  __m256d t;
  t = _mm256_unpacklo_pd(r, c);
  c = _mm256_unpackhi_pd(r, c);
  r = _mm256_set_m128d(_mm256_extractf128_pd(t,0), _mm256_extractf128_pd(c,0));
  c = _mm256_set_m128d(_mm256_extractf128_pd(t,1), _mm256_extractf128_pd(c,1));
}
#endif

template<typename T> class FixedSizeRealToComplexDFT
{
public:
  FixedSizeRealToComplexDFT(unsigned N, unsigned real_stride = 1, unsigned complex_stride = 1):
    N_(N), rs_(real_stride), cs_(complex_stride) { assert(rs_ != 0); assert(cs_ != 0); };
  virtual ~FixedSizeRealToComplexDFT() { }
  virtual void r2c(T* r_in, T* c_out) = 0;
  virtual void c2r(T* r_out, T* c_in) = 0;

  // Utility functions
  int real_stride() const { return rs_; }
  int complex_stride() const { return cs_; }
  unsigned num_real() const { return N_; }
  unsigned num_complex() const { return N_/2+1; }
  unsigned real_array_size() const { return num_real()*rs_; }
  unsigned complex_array_size() const { return 2*num_complex()*cs_; }
  T* alloc_real_array() const {
    return calin::util::memory::aligned_calloc<T>(real_array_size()); }
  T* alloc_complex_array() const {
    return calin::util::memory::aligned_calloc<T>(complex_array_size()); }
protected:
  unsigned N_;
  unsigned rs_;
  unsigned cs_;
};

#if 0
template<typename T> class FixedSizeRealToHalfComplexDFT
{
public:
  FixedSizeRealToHalfComplexDFT(unsigned N, unsigned real_stride = 1, unsigned half_complex_stride = 1):
    N_(N), rs_(real_stride), hcs_(half_complex_stride) { assert(rs_ != 0); assert(hcs_ != 0); };
  virtual ~FixedSizeRealToHalfComplexDFT() { }
  virtual void r2hc(T* r_in, T* hc_out) = 0;
  virtual void hc2r(T* r_out, T* hc_in) = 0;

  // Utility functions
  int real_stride() const { return rs_; }
  int half_complex_stride() const { return hcs_; }
  unsigned num_real() const { return N_; }
  unsigned num_half_complex() const { return N_; }
  unsigned real_array_size() const { return num_real()*rs_; }
  unsigned half_complex_array_size() const { return num_complex()*hcs_; }
  T* alloc_real_array() const {
    return calin::util::memory::aligned_calloc<T>(real_array_size()); }
  T* alloc_half_complex_array() const {
    return calin::util::memory::aligned_calloc<T>(complex_array_size()); }
protected:
  unsigned N_;
  unsigned rs_;
  unsigned hcs_;
};
#endif

template<typename T> class FFTW_FixedSizeRealToComplexDFT:
  public FixedSizeRealToComplexDFT<T>
{
public:
  FFTW_FixedSizeRealToComplexDFT(unsigned N, unsigned real_stride = 1, unsigned complex_stride = 1, unsigned flags = FFTW_MEASURE):
      FixedSizeRealToComplexDFT<T>(N, real_stride, complex_stride) {
    T* xt = FixedSizeRealToComplexDFT<T>::alloc_real_array();
    T* xf = FixedSizeRealToComplexDFT<T>::alloc_complex_array();
    int howmany = sizeof(T)/sizeof(double);
    int n = N;
    plan_r2c_ = fftw_plan_many_dft_r2c(1, &n, howmany,
      (double*)xt, nullptr, FixedSizeRealToComplexDFT<T>::rs_*howmany, 1,
      (fftw_complex*)xf, nullptr, FixedSizeRealToComplexDFT<T>::cs_*howmany, 1,
      flags);
    plan_c2r_ = fftw_plan_many_dft_c2r(1, &n, howmany,
      (fftw_complex*)xf, nullptr, FixedSizeRealToComplexDFT<T>::cs_*howmany, 1,
      (double*)xt, nullptr, FixedSizeRealToComplexDFT<T>::rs_*howmany, 1,
      flags);
    free(xf);
    free(xt);
  }
  virtual ~FFTW_FixedSizeRealToComplexDFT() {
    fftw_destroy_plan(plan_r2c_);
    fftw_destroy_plan(plan_c2r_);
  }
  void r2c(T* r_in, T* c_out) override {
    fftw_execute_dft_r2c(plan_r2c_, (double*)r_in, (fftw_complex*)c_out);
    unsigned nc = FixedSizeRealToComplexDFT<T>::num_complex();
    for(unsigned ic=0; ic<nc; ic++) {
      shuffle_c_v2s(c_out[2*FixedSizeRealToComplexDFT<T>::cs_*ic],
        c_out[2*FixedSizeRealToComplexDFT<T>::cs_*ic+1]);
    }
  }
  void c2r(T* r_out, T* c_in) override {
    unsigned nc = FixedSizeRealToComplexDFT<T>::num_complex();
    for(unsigned ic=0; ic<nc; ic++) {
      shuffle_c_s2v(c_in[2*FixedSizeRealToComplexDFT<T>::cs_*ic],
        c_in[2*FixedSizeRealToComplexDFT<T>::cs_*ic+1]);
    }
    fftw_execute_dft_c2r(plan_c2r_, (fftw_complex*)c_in, (double*)r_out);
  }
protected:
  fftw_plan plan_r2c_;
  fftw_plan plan_c2r_;
};

template<typename T> class FFTWF_FixedSizeRealToComplexDFT:
  public FixedSizeRealToComplexDFT<T>
{
public:
  FFTWF_FixedSizeRealToComplexDFT(unsigned N, unsigned real_stride = 1, unsigned complex_stride = 1, unsigned flags = FFTW_MEASURE):
      FixedSizeRealToComplexDFT<T>(N, real_stride, complex_stride) {
    T* xt = FixedSizeRealToComplexDFT<T>::alloc_real_array();
    T* xf = FixedSizeRealToComplexDFT<T>::alloc_complex_array();
    int howmany = sizeof(T)/sizeof(float);
    int n = N;
    plan_r2c_ = fftwf_plan_many_dft_r2c(1, &n, howmany,
      (float*)xt, nullptr, FixedSizeRealToComplexDFT<T>::rs_*howmany, 1,
      (fftwf_complex*)xf, nullptr, FixedSizeRealToComplexDFT<T>::cs_*howmany, 1,
      flags);
    plan_c2r_ = fftwf_plan_many_dft_c2r(1, &n, howmany,
      (fftwf_complex*)xf, nullptr, FixedSizeRealToComplexDFT<T>::cs_*howmany, 1,
      (float*)xt, nullptr, FixedSizeRealToComplexDFT<T>::rs_*howmany, 1,
      flags);
    free(xf);
    free(xt);
  }
  virtual ~FFTWF_FixedSizeRealToComplexDFT() {
    fftwf_destroy_plan(plan_r2c_);
    fftwf_destroy_plan(plan_c2r_);
  }
  void r2c(T* r_in, T* c_out) override {
    fftwf_execute_dft_r2c(plan_r2c_, (float*)r_in, (fftwf_complex*)c_out);
    unsigned nc = FixedSizeRealToComplexDFT<T>::num_complex();
    for(unsigned ic=0; ic<nc; ic++) {
      shuffle_c_v2s(c_out[2*FixedSizeRealToComplexDFT<T>::cs_*ic],
        c_out[2*FixedSizeRealToComplexDFT<T>::cs_*ic+1]);
    }
  }
  void c2r(T* r_out, T* c_in) override {
    unsigned nc = FixedSizeRealToComplexDFT<T>::num_complex();
    for(unsigned ic=0; ic<nc; ic++) {
      shuffle_c_s2v(c_in[2*FixedSizeRealToComplexDFT<T>::cs_*ic],
        c_in[2*FixedSizeRealToComplexDFT<T>::cs_*ic+1]);
    }
    fftwf_execute_dft_c2r(plan_c2r_, (fftwf_complex*)c_in, (float*)r_out);
  }
protected:
  fftwf_plan plan_r2c_;
  fftwf_plan plan_c2r_;
};

#if defined(__AVX__)
FixedSizeRealToComplexDFT<__m256>* new_m256_r2c_dft(unsigned n,
  unsigned real_stride = 1, unsigned complex_stride = 1);
FixedSizeRealToComplexDFT<__m256d>* new_m256d_r2c_dft(unsigned n,
  unsigned real_stride = 1, unsigned complex_stride = 1);

FixedSizeRealToComplexDFT<__m256>* new_m256_codelet_r2c_dft(unsigned n,
  unsigned real_stride = 1, unsigned complex_stride = 1);
FixedSizeRealToComplexDFT<__m256d>* new_m256d_codelet_r2c_dft(unsigned n,
  unsigned real_stride = 1, unsigned complex_stride = 1);
#endif // defined(__AVX__)

std::vector<unsigned> list_available_m256_codelets();
std::vector<unsigned> list_available_m256d_codelets();

std::vector<float> test_m256_r2c_dft(const std::vector<float>& data);
std::vector<float> test_m256_c2r_dft(const std::vector<float>& fft, unsigned n);
std::vector<float> test_fftw_m256_r2c_dft(const std::vector<float>& data);
std::vector<float> test_fftw_m256_c2r_dft(const std::vector<float>& fft, unsigned n);

std::vector<double> test_m256d_r2c_dft(const std::vector<double>& data);
std::vector<double> test_m256d_c2r_dft(const std::vector<double>& fft, unsigned n);
std::vector<double> test_fftw_m256d_r2c_dft(const std::vector<double>& data);
std::vector<double> test_fftw_m256d_c2r_dft(const std::vector<double>& fft, unsigned n);


} } } // namespace calin::math::simd_fft
