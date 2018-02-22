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
  }
  void c2r(T* r_out, T* c_in) override {
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
  }
  void c2r(T* r_out, T* c_in) override {
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

} } } // namespace calin::math::simd_fft
