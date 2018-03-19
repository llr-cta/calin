/*

   calin/math/fft_simd.cpp -- Stephen Fegan -- 2018-02-21

   SIMD FFT functions

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

#include <string>
#include <stdexcept>

#include <immintrin.h>

#include <math/fft_simd.hpp>

namespace calin { namespace math { namespace fft_simd {

#if defined(__AVX__)
namespace m256 {

using float_type = __m256;
using E = float_type;
using R = float_type;
using INT = int;
using stride = int;

inline int WS(const stride s, const stride i) { return s*i; }

inline E ADD(const E a, const E b) { return _mm256_add_ps(a,b); }
inline E SUB(const E a, const E b) { return _mm256_sub_ps(a,b); }
inline E MUL(const E a, const E b) { return _mm256_mul_ps(a,b); }

//inline E NEG(const E a) { return _mm256_sub_ps(_mm256_setzero_ps(),a); }
inline E NEG(const E a) { return _mm256_xor_ps(a, _mm256_set1_ps(-0.0)); }

#if defined(__FMA__)
inline E FMA(const E a, const E b, const E c) { return _mm256_fmadd_ps(a,b,c); }
inline E FMS(const E a, const E b, const E c) { return _mm256_fmsub_ps(a,b,c); }
// Note: inconsistency between FFTW and Intel intrinsics definitions of FNMA/S
inline E FNMA(const E a, const E b, const E c) { return _mm256_fnmsub_ps(a,b,c); }
inline E FNMS(const E a, const E b, const E c) { return _mm256_fnmadd_ps(a,b,c); }
#else
inline E FMA(const E a, const E b, const E c) { return _mm256_add_ps(_mm256_mul_ps(a,b),c); }
inline E FMS(const E a, const E b, const E c) { return _mm256_sub_ps(_mm256_mul_ps(a,b),c); }
inline E FNMA(const E a, const E b, const E c) { return _mm256_sub_ps(NEG(_mm256_mul_ps(a,b)),c); }
inline E FNMS(const E a, const E b, const E c) { return _mm256_add_ps(NEG(_mm256_mul_ps(a,b)),c); }
#endif

inline void MAKE_VOLATILE_STRIDE(int a, int b) { }

#define DK(name, val) \
  static const E name = { (val),(val),(val),(val),(val),(val),(val),(val) }

inline E ZERO() { return _mm256_setzero_ps(); }

#include "genfft_codelets/define_all_dft_codelets.cpp_incl"

#undef DK

} // namespace calin::math::fft_simd::m256

namespace m256d {

using float_type = __m256d;
using E = float_type;
using R = float_type;
using INT = int;
using stride = int;

inline int WS(const stride s, const stride i) { return s*i; }

inline E ADD(const E a, const E b) { return _mm256_add_pd(a,b); }
inline E SUB(const E a, const E b) { return _mm256_sub_pd(a,b); }
inline E MUL(const E a, const E b) { return _mm256_mul_pd(a,b); }

//inline E NEG(const E a) { return _mm256_sub_ps(_mm256_setzero_ps(),a); }
inline E NEG(const E a) { return _mm256_xor_pd(a, _mm256_set1_pd(-0.0)); }

#if defined(__FMA__)
inline E FMA(const E a, const E b, const E c) { return _mm256_fmadd_pd(a,b,c); }
inline E FMS(const E a, const E b, const E c) { return _mm256_fmsub_pd(a,b,c); }
// Note: inconsistency between FFTW and Intel intrinsics definitions of FNMA/S
inline E FNMA(const E a, const E b, const E c) { return _mm256_fnmsub_pd(a,b,c); }
inline E FNMS(const E a, const E b, const E c) { return _mm256_fnmadd_pd(a,b,c); }
#else
inline E FMA(const E a, const E b, const E c) { return _mm256_add_pd(_mm256_mul_pd(a,b),c); }
inline E FMS(const E a, const E b, const E c) { return _mm256_sub_pd(_mm256_mul_pd(a,b),c); }
inline E FNMA(const E a, const E b, const E c) { return _mm256_sub_pd(NEG(_mm256_mul_pd(a,b)),c); }
inline E FNMS(const E a, const E b, const E c) { return _mm256_add_pd(NEG(_mm256_mul_pd(a,b)),c); }
#endif

inline void MAKE_VOLATILE_STRIDE(int a, int b) { }

#define DK(name, val) \
  static const E name = { (val),(val),(val),(val) }

inline E ZERO() { return _mm256_setzero_pd(); }

#include "genfft_codelets/define_all_dft_codelets.cpp_incl"

#undef DK

} // namespace calin::math::fft_simd::m256d

FixedSizeRealToComplexDFT<__m256>* new_m256_r2c_dft(unsigned n,
  unsigned real_stride, unsigned complex_stride)
{
  FixedSizeRealToComplexDFT<__m256>* ptr
    = m256::new_codelet_r2c_dft(n, real_stride, complex_stride);
  if(ptr == nullptr)
    ptr = new FFTWF_FixedSizeRealToComplexDFT<__m256>(n, real_stride, complex_stride);
  return ptr;
}

FixedSizeRealToComplexDFT<__m256d>* new_m256d_r2c_dft(unsigned n,
  unsigned real_stride, unsigned complex_stride)
{
  FixedSizeRealToComplexDFT<__m256d>* ptr
    = m256d::new_codelet_r2c_dft(n, real_stride, complex_stride);
  if(ptr == nullptr)
    ptr = new FFTW_FixedSizeRealToComplexDFT<__m256d>(n, real_stride, complex_stride);
  return ptr;
}

FixedSizeRealToComplexDFT<__m256>* new_m256_codelet_r2c_dft(unsigned n,
  unsigned real_stride, unsigned complex_stride)
{
  return m256::new_codelet_r2c_dft(n, real_stride, complex_stride);
}

FixedSizeRealToComplexDFT<__m256d>* new_m256d_codelet_r2c_dft(unsigned n,
  unsigned real_stride, unsigned complex_stride)
{
  return m256d::new_codelet_r2c_dft(n, real_stride, complex_stride);
}

FixedSizeRealToComplexDFT<__m256>* new_m256_fftw_r2c_dft(unsigned n,
  unsigned real_stride, unsigned complex_stride)
{
  return new FFTWF_FixedSizeRealToComplexDFT<__m256>(n, real_stride, complex_stride);
}

FixedSizeRealToComplexDFT<__m256d>* new_m256d_fftw_r2c_dft(unsigned n,
  unsigned real_stride, unsigned complex_stride)
{
  return new FFTW_FixedSizeRealToComplexDFT<__m256d>(n, real_stride, complex_stride);
}

FixedSizeRealToHalfComplexDFT<__m256>* new_m256_r2hc_dft(unsigned n,
  unsigned real_stride, unsigned half_complex_stride)
{
  FixedSizeRealToHalfComplexDFT<__m256>* ptr
    = m256::new_codelet_r2hc_dft(n, real_stride, half_complex_stride);
  if(ptr == nullptr)
    ptr = new FFTWF_FixedSizeRealToHalfComplexDFT<__m256>(n, real_stride, half_complex_stride);
  return ptr;
}

FixedSizeRealToHalfComplexDFT<__m256d>* new_m256d_r2hc_dft(unsigned n,
  unsigned real_stride, unsigned half_complex_stride)
{
  FixedSizeRealToHalfComplexDFT<__m256d>* ptr
    = m256d::new_codelet_r2hc_dft(n, real_stride, half_complex_stride);
  if(ptr == nullptr)
    ptr = new FFTW_FixedSizeRealToHalfComplexDFT<__m256d>(n, real_stride, half_complex_stride);
  return ptr;
}

FixedSizeRealToHalfComplexDFT<__m256>* new_m256_codelet_r2hc_dft(unsigned n,
  unsigned real_stride, unsigned half_complex_stride)
{
  return m256::new_codelet_r2hc_dft(n, real_stride, half_complex_stride);
}

FixedSizeRealToHalfComplexDFT<__m256d>* new_m256d_codelet_r2hc_dft(unsigned n,
  unsigned real_stride, unsigned half_complex_stride)
{
  return m256d::new_codelet_r2hc_dft(n, real_stride, half_complex_stride);
}

FixedSizeRealToHalfComplexDFT<__m256>* new_m256_fftw_r2hc_dft(unsigned n,
  unsigned real_stride, unsigned half_complex_stride)
{
  return new FFTWF_FixedSizeRealToHalfComplexDFT<__m256>(n, real_stride, half_complex_stride);
}

FixedSizeRealToHalfComplexDFT<__m256d>* new_m256d_fftw_r2hc_dft(unsigned n,
  unsigned real_stride, unsigned half_complex_stride)
{
  return new FFTW_FixedSizeRealToHalfComplexDFT<__m256d>(n, real_stride, half_complex_stride);
}

std::vector<unsigned> list_available_m256_codelets()
{
  return m256::list_available_codelets();
}

std::vector<unsigned> list_available_m256d_codelets()
{
  return m256d::list_available_codelets();
}

#endif // defined(__AVX__)

// *****************************************************************************
// ******************* FLOAT AVX REAL <-> COMPLEX TEST CODE ********************
// *****************************************************************************

std::vector<float> test_m256_r2c_dft(const std::vector<float>& data)
{
#if defined(__AVX__)
  auto* dft = new_m256_r2c_dft(data.size());
  std::vector<float> fft(dft->complex_array_size());
  auto* xt = dft->alloc_real_array();
  auto* xf = dft->alloc_complex_array();
  std::transform(data.begin(), data.end(), xt, [](float x){ return _mm256_set_ps(0,0,0,0,0,0,0,x); });
  dft->r2c(xt, xf);
  std::transform(xf, xf+fft.size(), fft.begin(), [](__m256 x){ return /* _mm256_cvtss_f32(x) */ (float)x[0]; });
  free(xf);
  free(xt);
  delete dft;
  return fft;
#else // defined(__AVX__)
  throw std::runtime_error("AVX not present at compile type");
#endif // defined(__AVX__)
}

std::vector<float> test_m256_c2r_dft(const std::vector<float>& fft, unsigned n)
{
#if defined(__AVX__)
  if(fft.size() == 0)throw std::runtime_error("Input DFT empty");
  if(fft.size() != 2*(n/2+1))throw std::runtime_error("Input DFT must have size " + std::to_string(2*(n/2+1)));
  auto* dft = new_m256_r2c_dft(n);
  std::vector<float> data(dft->real_array_size());
  auto* xt = dft->alloc_real_array();
  auto* xf = dft->alloc_complex_array();
  std::transform(fft.begin(), fft.end(), xf, [](float x){ return _mm256_set_ps(0,0,0,0,0,0,0,x); });
  dft->c2r(xt, xf);
  std::transform(xt, xt+data.size(), data.begin(), [](__m256 x){ return /* _mm256_cvtss_f32(x) */ (float)x[0]; });
  free(xf);
  free(xt);
  delete dft;
  return data;
#else // defined(__AVX__)
  throw std::runtime_error("AVX not present at compile type");
#endif // defined(__AVX__)
}

std::vector<float> test_fftw_m256_r2c_dft(const std::vector<float>& data)
{
#if defined(__AVX__)
  auto* dft = new FFTWF_FixedSizeRealToComplexDFT<__m256>(data.size());
  std::vector<float> fft(dft->complex_array_size());
  auto* xt = dft->alloc_real_array();
  auto* xf = dft->alloc_complex_array();
  std::transform(data.begin(), data.end(), xt, [](float x){ return _mm256_set_ps(0,0,0,0,0,0,0,x); });
  dft->r2c(xt, xf);
  std::transform(xf, xf+fft.size(), fft.begin(), [](__m256 x){ return /* _mm256_cvtss_f32(x) */ x[0]; });
  free(xf);
  free(xt);
  delete dft;
  return fft;
#else // defined(__AVX__)
  throw std::runtime_error("AVX not present at compile type");
#endif // defined(__AVX__)
}

std::vector<float> test_fftw_m256_c2r_dft(const std::vector<float>& fft, unsigned n)
{
#if defined(__AVX__)
  if(fft.size() == 0)throw std::runtime_error("Input DFT empty");
  if(fft.size() != 2*(n/2+1))throw std::runtime_error("Input DFT must have size " + std::to_string(2*(n/2+1)));
  auto* dft = new FFTWF_FixedSizeRealToComplexDFT<__m256>(n);
  std::vector<float> data(dft->real_array_size());
  auto* xt = dft->alloc_real_array();
  auto* xf = dft->alloc_complex_array();
  std::transform(fft.begin(), fft.end(), xf, [](float x){ return _mm256_set_ps(0,0,0,0,0,0,0,x); });
  dft->c2r(xt, xf);
  std::transform(xt, xt+data.size(), data.begin(), [](__m256 x){ return /* _mm256_cvtss_f32(x) */ x[0]; });
  free(xf);
  free(xt);
  delete dft;
  return data;
#else // defined(__AVX__)
  throw std::runtime_error("AVX not present at compile type");
#endif // defined(__AVX__)
}

// *****************************************************************************
// ******************* DOUBLE AVX REAL <-> COMPLEX TEST CODE *******************
// *****************************************************************************

std::vector<double> test_m256d_r2c_dft(const std::vector<double>& data)
{
#if defined(__AVX__)
  auto* dft = new_m256d_r2c_dft(data.size());
  std::vector<double> fft(dft->complex_array_size());
  auto* xt = dft->alloc_real_array();
  auto* xf = dft->alloc_complex_array();
  std::transform(data.begin(), data.end(), xt, [](double x){ return _mm256_set_pd(0,0,0,x); });
  dft->r2c(xt, xf);
  std::transform(xf, xf+fft.size(), fft.begin(), [](__m256d x){ return /* _mm256_cvtsd_f64(x) */ x[0]; });
  free(xf);
  free(xt);
  delete dft;
  return fft;
#else // defined(__AVX__)
  throw std::runtime_error("AVX not present at compile type");
#endif // defined(__AVX__)
}

std::vector<double> test_m256d_c2r_dft(const std::vector<double>& fft, unsigned n)
{
#if defined(__AVX__)
  if(fft.size() == 0)throw std::runtime_error("Input DFT empty");
  if(fft.size() != 2*(n/2+1))throw std::runtime_error("Input DFT must have size " + std::to_string(2*(n/2+1)));
  auto* dft = new_m256d_r2c_dft(n);
  std::vector<double> data(dft->real_array_size());
  auto* xt = dft->alloc_real_array();
  auto* xf = dft->alloc_complex_array();
  std::transform(fft.begin(), fft.end(), xf, [](double x){ return _mm256_set_pd(0,0,0,x); });
  dft->c2r(xt, xf);
  std::transform(xt, xt+data.size(), data.begin(), [](__m256d x){ return /* _mm256_cvtsd_f64(x) */ x[0]; });
  free(xf);
  free(xt);
  delete dft;
  return data;
#else // defined(__AVX__)
  throw std::runtime_error("AVX not present at compile type");
#endif // defined(__AVX__)
}

std::vector<double> test_fftw_m256d_r2c_dft(const std::vector<double>& data)
{
#if defined(__AVX__)
  auto* dft = new FFTW_FixedSizeRealToComplexDFT<__m256d>(data.size());
  std::vector<double> fft(dft->complex_array_size());
  auto* xt = dft->alloc_real_array();
  auto* xf = dft->alloc_complex_array();
  std::transform(data.begin(), data.end(), xt, [](double x){ return _mm256_set_pd(0,0,0,x); });
  dft->r2c(xt, xf);
  std::transform(xf, xf+fft.size(), fft.begin(), [](__m256d x){ return /* _mm256_cvtsd_f64(x) */ x[0]; });
  free(xf);
  free(xt);
  delete dft;
  return fft;
#else // defined(__AVX__)
  throw std::runtime_error("AVX not present at compile type");
#endif // defined(__AVX__)
}

std::vector<double> test_fftw_m256d_c2r_dft(const std::vector<double>& fft, unsigned n)
{
#if defined(__AVX__)
  if(fft.size() == 0)throw std::runtime_error("Input DFT empty");
  if(fft.size() != 2*(n/2+1))throw std::runtime_error("Input DFT must have size " + std::to_string(2*(n/2+1)));
  auto* dft = new FFTW_FixedSizeRealToComplexDFT<__m256d>(n);
  std::vector<double> data(dft->real_array_size());
  auto* xt = dft->alloc_real_array();
  auto* xf = dft->alloc_complex_array();
  std::transform(fft.begin(), fft.end(), xf, [](double x){ return _mm256_set_pd(0,0,0,x); });
  dft->c2r(xt, xf);
  std::transform(xt, xt+data.size(), data.begin(), [](__m256d x){ return /* _mm256_cvtsd_f64(x) */ x[0]; });
  free(xf);
  free(xt);
  delete dft;
  return data;
#else // defined(__AVX__)
  throw std::runtime_error("AVX not present at compile type");
#endif // defined(__AVX__)
}

// *****************************************************************************
// ***************** FLOAT AVX REAL <-> HALF COMPLEX TEST CODE *****************
// *****************************************************************************

std::vector<float> test_m256_r2hc_dft(const std::vector<float>& data)
{
#if defined(__AVX__)
  auto* dft = new_m256_r2hc_dft(data.size());
  std::vector<float> fft(dft->half_complex_array_size());
  auto* xt = dft->alloc_real_array();
  auto* xf = dft->alloc_half_complex_array();
  std::transform(data.begin(), data.end(), xt, [](float x){ return _mm256_set_ps(0,0,0,0,0,0,0,x); });
  dft->r2hc(xt, xf);
  std::transform(xf, xf+fft.size(), fft.begin(), [](__m256 x){ return /* _mm256_cvtss_f32(x) */ x[0]; });
  free(xf);
  free(xt);
  delete dft;
  return fft;
#else // defined(__AVX__)
  throw std::runtime_error("AVX not present at compile type");
#endif // defined(__AVX__)
}

std::vector<float> test_m256_hc2r_dft(const std::vector<float>& fft)
{
#if defined(__AVX__)
  auto* dft = new_m256_r2hc_dft(fft.size());
  std::vector<float> data(dft->real_array_size());
  auto* xt = dft->alloc_real_array();
  auto* xf = dft->alloc_half_complex_array();
  std::transform(fft.begin(), fft.end(), xf, [](float x){ return _mm256_set_ps(0,0,0,0,0,0,0,x); });
  dft->hc2r(xt, xf);
  std::transform(xt, xt+data.size(), data.begin(), [](__m256 x){ return /* _mm256_cvtss_f32(x) */ x[0]; });
  free(xf);
  free(xt);
  delete dft;
  return data;
#else // defined(__AVX__)
  throw std::runtime_error("AVX not present at compile type");
#endif // defined(__AVX__)
}

std::vector<float> test_fftw_m256_r2hc_dft(const std::vector<float>& data)
{
#if defined(__AVX__)
  auto* dft = new FFTWF_FixedSizeRealToHalfComplexDFT<__m256>(data.size());
  std::vector<float> fft(dft->half_complex_array_size());
  auto* xt = dft->alloc_real_array();
  auto* xf = dft->alloc_half_complex_array();
  std::transform(data.begin(), data.end(), xt, [](float x){ return _mm256_set_ps(0,0,0,0,0,0,0,x); });
  dft->r2hc(xt, xf);
  std::transform(xf, xf+fft.size(), fft.begin(), [](__m256 x){ return /* _mm256_cvtss_f32(x) */ x[0]; });
  free(xf);
  free(xt);
  delete dft;
  return fft;
#else // defined(__AVX__)
  throw std::runtime_error("AVX not present at compile type");
#endif // defined(__AVX__)
}

std::vector<float> test_fftw_m256_hc2r_dft(const std::vector<float>& fft)
{
#if defined(__AVX__)
  auto* dft = new FFTWF_FixedSizeRealToHalfComplexDFT<__m256>(fft.size());
  std::vector<float> data(dft->real_array_size());
  auto* xt = dft->alloc_real_array();
  auto* xf = dft->alloc_half_complex_array();
  std::transform(fft.begin(), fft.end(), xf, [](float x){ return _mm256_set_ps(0,0,0,0,0,0,0,x); });
  dft->hc2r(xt, xf);
  std::transform(xt, xt+data.size(), data.begin(), [](__m256 x){ return /* _mm256_cvtss_f32(x) */ x[0]; });
  free(xf);
  free(xt);
  delete dft;
  return data;
#else // defined(__AVX__)
  throw std::runtime_error("AVX not present at compile type");
#endif // defined(__AVX__)
}

// *****************************************************************************
// ***************** DOUBLE AVX REAL <-> HALF COMPLEX TEST CODE ****************
// *****************************************************************************

std::vector<double> test_m256d_r2hc_dft(const std::vector<double>& data)
{
#if defined(__AVX__)
  auto* dft = new_m256d_r2hc_dft(data.size());
  std::vector<double> fft(dft->half_complex_array_size());
  auto* xt = dft->alloc_real_array();
  auto* xf = dft->alloc_half_complex_array();
  std::transform(data.begin(), data.end(), xt, [](double x){ return _mm256_set_pd(0,0,0,x); });
  dft->r2hc(xt, xf);
  std::transform(xf, xf+fft.size(), fft.begin(), [](__m256d x){ return /* _mm256_cvtsd_f64(x) */ x[0]; });
  free(xf);
  free(xt);
  delete dft;
  return fft;
#else // defined(__AVX__)
  throw std::runtime_error("AVX not present at compile type");
#endif // defined(__AVX__)
}

std::vector<double> test_m256d_hc2r_dft(const std::vector<double>& fft)
{
#if defined(__AVX__)
  auto* dft = new_m256d_r2hc_dft(fft.size());
  std::vector<double> data(dft->real_array_size());
  auto* xt = dft->alloc_real_array();
  auto* xf = dft->alloc_half_complex_array();
  std::transform(fft.begin(), fft.end(), xf, [](double x){ return _mm256_set_pd(0,0,0,x); });
  dft->hc2r(xt, xf);
  std::transform(xt, xt+data.size(), data.begin(), [](__m256d x){ return /* _mm256_cvtsd_f64(x) */ x[0]; });
  free(xf);
  free(xt);
  delete dft;
  return data;
#else // defined(__AVX__)
  throw std::runtime_error("AVX not present at compile type");
#endif // defined(__AVX__)
}

std::vector<double> test_fftw_m256d_r2hc_dft(const std::vector<double>& data)
{
#if defined(__AVX__)
  auto* dft = new FFTW_FixedSizeRealToHalfComplexDFT<__m256d>(data.size());
  std::vector<double> fft(dft->half_complex_array_size());
  auto* xt = dft->alloc_real_array();
  auto* xf = dft->alloc_half_complex_array();
  std::transform(data.begin(), data.end(), xt, [](double x){ return _mm256_set_pd(0,0,0,x); });
  dft->r2hc(xt, xf);
  std::transform(xf, xf+fft.size(), fft.begin(), [](__m256d x){ return /* _mm256_cvtsd_f64(x) */ x[0]; });
  free(xf);
  free(xt);
  delete dft;
  return fft;
#else // defined(__AVX__)
  throw std::runtime_error("AVX not present at compile type");
#endif // defined(__AVX__)
}

std::vector<double> test_fftw_m256d_hc2r_dft(const std::vector<double>& fft)
{
#if defined(__AVX__)
  auto* dft = new FFTW_FixedSizeRealToHalfComplexDFT<__m256d>(fft.size());
  std::vector<double> data(dft->real_array_size());
  auto* xt = dft->alloc_real_array();
  auto* xf = dft->alloc_half_complex_array();
  std::transform(fft.begin(), fft.end(), xf, [](double x){ return _mm256_set_pd(0,0,0,x); });
  dft->hc2r(xt, xf);
  std::transform(xt, xt+data.size(), data.begin(), [](__m256d x){ return /* _mm256_cvtsd_f64(x) */ x[0]; });
  free(xf);
  free(xt);
  delete dft;
  return data;
#else // defined(__AVX__)
  throw std::runtime_error("AVX not present at compile type");
#endif // defined(__AVX__)
}

} } } // namespace calin::math::fft_simd
