/*

   calin/math/simd_fft.cpp -- Stephen Fegan -- 2018-02-21

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

#include <math/simd_fft.hpp>

namespace calin { namespace math { namespace simd_fft {

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

} // namespace calin::math::simd_fft::m256

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

} // namespace calin::math::simd_fft::m256d

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

std::vector<unsigned> list_available_m256_codelets()
{
  return m256::list_available_codelets();
}

std::vector<unsigned> list_available_m256d_codelets()
{
  return m256d::list_available_codelets();
}

#endif // defined(__AVX__)


} } } // namespace calin::math::simd_fft
