/*

   calin/util/vcl.cpp -- Stephen Fegan -- 2018-08-14

   Vector class library helpers

   Copyright 2018, Stephen Fegan <sfegan@llr.in2p3.fr>
   Laboratoire Leprince-Ringuet, CNRS/IN2P3, Ecole Polytechnique, Institut Polytechnique de Paris

   This file is part of "calin"

   "calin" is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 or
   later, as published by the Free Software Foundation.

   "calin" is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

*/

#include <util/vcl.hpp>

using namespace calin::util::vcl;

constexpr char VCL128Architecture::architecture_name[];
constexpr unsigned VCL128Architecture::vec_bits;
constexpr unsigned VCL128Architecture::vec_bytes;
constexpr unsigned VCL128Architecture::num_int8;
constexpr unsigned VCL128Architecture::num_uint8;
constexpr unsigned VCL128Architecture::num_int16;
constexpr unsigned VCL128Architecture::num_uint16;
constexpr unsigned VCL128Architecture::num_int32;
constexpr unsigned VCL128Architecture::num_uint32;
constexpr unsigned VCL128Architecture::num_int64;
constexpr unsigned VCL128Architecture::num_uint64;
constexpr unsigned VCL128Architecture::num_float;
constexpr unsigned VCL128Architecture::num_double;

constexpr char VCL256Architecture::architecture_name[];
constexpr unsigned VCL256Architecture::vec_bits;
constexpr unsigned VCL256Architecture::vec_bytes;
constexpr unsigned VCL256Architecture::num_int8;
constexpr unsigned VCL256Architecture::num_uint8;
constexpr unsigned VCL256Architecture::num_int16;
constexpr unsigned VCL256Architecture::num_uint16;
constexpr unsigned VCL256Architecture::num_int32;
constexpr unsigned VCL256Architecture::num_uint32;
constexpr unsigned VCL256Architecture::num_int64;
constexpr unsigned VCL256Architecture::num_uint64;
constexpr unsigned VCL256Architecture::num_float;
constexpr unsigned VCL256Architecture::num_double;

constexpr char VCL512Architecture::architecture_name[];
constexpr unsigned VCL512Architecture::vec_bits;
constexpr unsigned VCL512Architecture::vec_bytes;
constexpr unsigned VCL512Architecture::num_int8;
constexpr unsigned VCL512Architecture::num_uint8;
constexpr unsigned VCL512Architecture::num_int16;
constexpr unsigned VCL512Architecture::num_uint16;
constexpr unsigned VCL512Architecture::num_int32;
constexpr unsigned VCL512Architecture::num_uint32;
constexpr unsigned VCL512Architecture::num_int64;
constexpr unsigned VCL512Architecture::num_uint64;
constexpr unsigned VCL512Architecture::num_float;
constexpr unsigned VCL512Architecture::num_double;

namespace {
  template<typename T> inline void do_one_128_swizzle_16(T& a, T& b) {
    __m128i tmp = _mm_unpackhi_epi16(a, b);
    a = _mm_unpacklo_epi16(a, b);
    b = tmp;
  }

  template<typename T> inline void do_one_128_swizzle_32(T& a, T& b) {
    __m128i tmp = _mm_unpackhi_epi32(a, b);
    a = _mm_unpacklo_epi32(a, b);
    b = tmp;
  }

  template<typename T> inline void do_one_128_swizzle_64(T& a, T& b) {
    __m128i tmp = _mm_unpackhi_epi64(a, b);
    a = _mm_unpacklo_epi64(a, b);
    b = tmp;
  }

  inline void do_one_128_swizzle_ps(Vec4f& a, Vec4f& b) {
    __m128 tmp = _mm_unpackhi_ps(a, b);
    a = _mm_unpacklo_ps(a, b);
    b = tmp;
  }

  inline void do_one_128_swizzle_pd(Vec4f& a, Vec4f& b) {
    __m128d tmp = _mm_unpackhi_pd(vcl::reinterpret_d(a), vcl::reinterpret_d(b));
    a = _mm_castpd_ps(_mm_unpacklo_pd(vcl::reinterpret_d(a), vcl::reinterpret_d(b)));
    b = _mm_castpd_ps(tmp);
  }

  inline void do_one_128_swizzle_pd(Vec2d& a, Vec2d& b) {
    __m128d tmp = _mm_unpackhi_pd(a, b);
    a = _mm_unpacklo_pd(a, b);
    b = tmp;
  }

}

void calin::util::vcl::transpose(Vec8s* x)
{
  do_one_128_swizzle_16(x[0], x[1]);
  do_one_128_swizzle_16(x[2], x[3]);
  do_one_128_swizzle_16(x[4], x[5]);
  do_one_128_swizzle_16(x[6], x[7]);

  do_one_128_swizzle_32(x[0], x[2]);
  do_one_128_swizzle_32(x[1], x[3]);
  do_one_128_swizzle_32(x[4], x[6]);
  do_one_128_swizzle_32(x[5], x[7]);

  do_one_128_swizzle_64(x[0], x[4]);
  do_one_128_swizzle_64(x[1], x[5]);
  do_one_128_swizzle_64(x[2], x[6]);
  do_one_128_swizzle_64(x[3], x[7]);

  std::swap(x[1], x[4]);
  std::swap(x[3], x[6]);
}

void calin::util::vcl::transpose(Vec8us* x)
{
  do_one_128_swizzle_16(x[0], x[1]);
  do_one_128_swizzle_16(x[2], x[3]);
  do_one_128_swizzle_16(x[4], x[5]);
  do_one_128_swizzle_16(x[6], x[7]);

  do_one_128_swizzle_32(x[0], x[2]);
  do_one_128_swizzle_32(x[1], x[3]);
  do_one_128_swizzle_32(x[4], x[6]);
  do_one_128_swizzle_32(x[5], x[7]);

  do_one_128_swizzle_64(x[0], x[4]);
  do_one_128_swizzle_64(x[1], x[5]);
  do_one_128_swizzle_64(x[2], x[6]);
  do_one_128_swizzle_64(x[3], x[7]);

  std::swap(x[1], x[4]);
  std::swap(x[3], x[6]);
}

void calin::util::vcl::transpose(Vec4i* x)
{
  do_one_128_swizzle_32(x[0], x[1]);
  do_one_128_swizzle_32(x[2], x[3]);

  do_one_128_swizzle_64(x[0], x[2]);
  do_one_128_swizzle_64(x[1], x[3]);

  std::swap(x[1], x[2]);
}

void calin::util::vcl::transpose(Vec4ui* x)
{
  do_one_128_swizzle_32(x[0], x[1]);
  do_one_128_swizzle_32(x[2], x[3]);

  do_one_128_swizzle_64(x[0], x[2]);
  do_one_128_swizzle_64(x[1], x[3]);

  std::swap(x[1], x[2]);
}

void calin::util::vcl::transpose(Vec2q* x)
{
  do_one_128_swizzle_64(x[0], x[1]);
}

void calin::util::vcl::transpose(Vec2uq* x)
{
  do_one_128_swizzle_64(x[0], x[1]);
}

void calin::util::vcl::transpose(Vec4f* x)
{
  do_one_128_swizzle_ps(x[0], x[1]);
  do_one_128_swizzle_ps(x[2], x[3]);

  do_one_128_swizzle_pd(x[0], x[2]);
  do_one_128_swizzle_pd(x[1], x[3]);

  std::swap(x[1], x[2]);
}

void calin::util::vcl::transpose(Vec2d* x)
{
  do_one_128_swizzle_pd(x[0], x[1]);
}

#if MAX_VECTOR_SIZE >= 256

namespace {
#if INSTRSET >= 8
  template<typename T> inline void do_one_256_swizzle_16(T& a, T& b) {
    __m256i tmp = _mm256_unpackhi_epi16(a, b);
    a = _mm256_unpacklo_epi16(a, b);
    b = tmp;
  }

  template<typename T> inline void do_one_256_swizzle_32(T& a, T& b) {
    __m256i tmp = _mm256_unpackhi_epi32(a, b);
    a = _mm256_unpacklo_epi32(a, b);
    b = tmp;
  }

  template<typename T> inline void do_one_256_swizzle_64(T& a, T& b) {
    __m256i tmp = _mm256_unpackhi_epi64(a, b);
    a = _mm256_unpacklo_epi64(a, b);
    b = tmp;
  }

  template<typename T> inline void do_one_256_swizzle_128(T& a, T& b) {
    __m256i tmp =  _mm256_permute2x128_si256(a, b, 0x31);
    a = _mm256_permute2x128_si256(a, b, 0x20);
    b = tmp;
  }

  template<typename T> inline void do_one_256_swizzle_ps(T& a, T& b) {
    T tmp = _mm256_unpackhi_ps(a, b);
    a = _mm256_unpacklo_ps(a, b);
    b = tmp;
  }

  inline void do_one_256_swizzle_pd(Vec8f& a, Vec8f& b) {
    __m256d tmp = _mm256_unpackhi_pd(vcl::reinterpret_d(a), vcl::reinterpret_d(b));
    a = _mm256_castpd_ps(_mm256_unpacklo_pd(vcl::reinterpret_d(a), vcl::reinterpret_d(b)));
    b = _mm256_castpd_ps(tmp);
  }

  inline void do_one_256_swizzle_pd(Vec4d& a, Vec4d& b) {
    __m256d tmp = _mm256_unpackhi_pd(a, b);
    a = _mm256_unpacklo_pd(a, b);
    b = tmp;
  }

  inline void do_one_256_swizzle_flt128(Vec8f& a, Vec8f& b) {
    __m256 tmp =  _mm256_permute2f128_ps(a, b, 0x31);
    a = _mm256_permute2f128_ps(a, b, 0x20);
    b = tmp;
  }

  inline void do_one_256_swizzle_flt128(Vec4d& a, Vec4d& b) {
    __m256d tmp =  _mm256_permute2f128_pd(a, b, 0x31);
    a = _mm256_permute2f128_pd(a, b, 0x20);
    b = tmp;
  }

#else // INSTRSET >= 8
  template<typename T> inline void do_one_256_swizzle_16(T& a, T& b) {
    T tmp(_mm_unpackhi_epi16(a.get_low(), b.get_low()),
          _mm_unpackhi_epi16(a.get_high(), b.get_high()));
    a = T(_mm_unpacklo_epi16(a.get_low(), b.get_low()),
          _mm_unpacklo_epi16(a.get_high(), b.get_high()));
    b = tmp;
  }

  template<typename T> inline void do_one_256_swizzle_32(T& a, T& b) {
    T tmp(_mm_unpackhi_epi32(a.get_low(), b.get_low()),
          _mm_unpackhi_epi32(a.get_high(), b.get_high()));
    a = T(_mm_unpacklo_epi32(a.get_low(), b.get_low()),
          _mm_unpacklo_epi32(a.get_high(), b.get_high()));
    b = tmp;
  }

  template<typename T> inline void do_one_256_swizzle_64(T& a, T& b) {
    T tmp(_mm_unpackhi_epi64(a.get_low(), b.get_low()),
          _mm_unpackhi_epi64(a.get_high(), b.get_high()));
    a = T(_mm_unpacklo_epi64(a.get_low(), b.get_low()),
          _mm_unpacklo_epi64(a.get_high(), b.get_high()));
    b = tmp;
  }

  template<typename T> inline void do_one_256_swizzle_128(T& a, T& b) {
    T tmp(a.get_high(), b.get_high());
    a = T(a.get_low(), b.get_low());
    b = tmp;
  }

  inline void do_one_256_swizzle_ps(Vec8f& a, Vec8f& b) {
    Vec8f tmp(_mm_unpackhi_ps(a.get_low(), b.get_low()),
          _mm_unpackhi_ps(a.get_high(), b.get_high()));
    a = Vec8f(_mm_unpacklo_ps(a.get_low(), b.get_low()),
          _mm_unpacklo_ps(a.get_high(), b.get_high()));
    b = tmp;
  }

  inline void do_one_256_swizzle_pd(Vec8f& a, Vec8f& b) {
    Vec8f tmp(_mm_castpd_ps(_mm_unpackhi_pd(vcl::reinterpret_d(a.get_low()), vcl::reinterpret_d(b.get_low()))),
              _mm_castpd_ps(_mm_unpackhi_pd(vcl::reinterpret_d(a.get_high()), vcl::reinterpret_d(b.get_high()))));
    a = Vec8f(_mm_castpd_ps(_mm_unpacklo_pd(vcl::reinterpret_d(a.get_low()), vcl::reinterpret_d(b.get_low()))),
              _mm_castpd_ps(_mm_unpacklo_pd(vcl::reinterpret_d(a.get_high()), vcl::reinterpret_d(b.get_high()))));
    b = tmp;
  }

  inline void do_one_256_swizzle_pd(Vec4d& a, Vec4d& b) {
    Vec4d tmp(_mm_unpackhi_pd(a.get_low(), b.get_low()),
              _mm_unpackhi_pd(a.get_high(), b.get_high()));
    a = Vec4d(_mm_unpacklo_pd(a.get_low(), b.get_low()),
              _mm_unpacklo_pd(a.get_high(), b.get_high()));
    b = tmp;
  }



  template<typename T> inline void do_one_256_swizzle_flt128(T& a, T& b) {
    T tmp(a.get_high(), b.get_high());
    a = T(a.get_low(), b.get_low());
    b = tmp;
  }

#endif // INSTRSET >= 8
}

void calin::util::vcl::transpose(Vec16s* x)
{
  do_one_256_swizzle_16(x[0],   x[1]);
  do_one_256_swizzle_16(x[2],   x[3]);
  do_one_256_swizzle_16(x[4],   x[5]);
  do_one_256_swizzle_16(x[6],   x[7]);
  do_one_256_swizzle_16(x[8],   x[9]);
  do_one_256_swizzle_16(x[10],  x[11]);
  do_one_256_swizzle_16(x[12],  x[13]);
  do_one_256_swizzle_16(x[14],  x[15]);

  do_one_256_swizzle_32(x[0],   x[2]);
  do_one_256_swizzle_32(x[1],   x[3]);
  do_one_256_swizzle_32(x[4],   x[6]);
  do_one_256_swizzle_32(x[5],   x[7]);
  do_one_256_swizzle_32(x[8],   x[10]);
  do_one_256_swizzle_32(x[9],   x[11]);
  do_one_256_swizzle_32(x[12],  x[14]);
  do_one_256_swizzle_32(x[13],  x[15]);

  do_one_256_swizzle_64(x[0],   x[4]);
  do_one_256_swizzle_64(x[1],   x[5]);
  do_one_256_swizzle_64(x[2],   x[6]);
  do_one_256_swizzle_64(x[3],   x[7]);
  do_one_256_swizzle_64(x[8],   x[12]);
  do_one_256_swizzle_64(x[9],   x[13]);
  do_one_256_swizzle_64(x[10],  x[14]);
  do_one_256_swizzle_64(x[11],  x[15]);

  do_one_256_swizzle_128(x[0],  x[8]);
  do_one_256_swizzle_128(x[1],  x[9]);
  do_one_256_swizzle_128(x[2],  x[10]);
  do_one_256_swizzle_128(x[3],  x[11]);
  do_one_256_swizzle_128(x[4],  x[12]);
  do_one_256_swizzle_128(x[5],  x[13]);
  do_one_256_swizzle_128(x[6],  x[14]);
  do_one_256_swizzle_128(x[7],  x[15]);

  std::swap(x[1],  x[4]);
  std::swap(x[3],  x[6]);
  std::swap(x[9],  x[12]);
  std::swap(x[11], x[14]);
}

void calin::util::vcl::transpose(Vec16us* x)
{
  do_one_256_swizzle_16(x[0],   x[1]);
  do_one_256_swizzle_16(x[2],   x[3]);
  do_one_256_swizzle_16(x[4],   x[5]);
  do_one_256_swizzle_16(x[6],   x[7]);
  do_one_256_swizzle_16(x[8],   x[9]);
  do_one_256_swizzle_16(x[10],  x[11]);
  do_one_256_swizzle_16(x[12],  x[13]);
  do_one_256_swizzle_16(x[14],  x[15]);

  do_one_256_swizzle_32(x[0],   x[2]);
  do_one_256_swizzle_32(x[1],   x[3]);
  do_one_256_swizzle_32(x[4],   x[6]);
  do_one_256_swizzle_32(x[5],   x[7]);
  do_one_256_swizzle_32(x[8],   x[10]);
  do_one_256_swizzle_32(x[9],   x[11]);
  do_one_256_swizzle_32(x[12],  x[14]);
  do_one_256_swizzle_32(x[13],  x[15]);

  do_one_256_swizzle_64(x[0],   x[4]);
  do_one_256_swizzle_64(x[1],   x[5]);
  do_one_256_swizzle_64(x[2],   x[6]);
  do_one_256_swizzle_64(x[3],   x[7]);
  do_one_256_swizzle_64(x[8],   x[12]);
  do_one_256_swizzle_64(x[9],   x[13]);
  do_one_256_swizzle_64(x[10],  x[14]);
  do_one_256_swizzle_64(x[11],  x[15]);

  do_one_256_swizzle_128(x[0],  x[8]);
  do_one_256_swizzle_128(x[1],  x[9]);
  do_one_256_swizzle_128(x[2],  x[10]);
  do_one_256_swizzle_128(x[3],  x[11]);
  do_one_256_swizzle_128(x[4],  x[12]);
  do_one_256_swizzle_128(x[5],  x[13]);
  do_one_256_swizzle_128(x[6],  x[14]);
  do_one_256_swizzle_128(x[7],  x[15]);

  std::swap(x[1],  x[4]);
  std::swap(x[3],  x[6]);
  std::swap(x[9],  x[12]);
  std::swap(x[11], x[14]);
}

void calin::util::vcl::transpose(Vec8i* x)
{
  do_one_256_swizzle_32(x[0],   x[1]);
  do_one_256_swizzle_32(x[2],   x[3]);
  do_one_256_swizzle_32(x[4],   x[5]);
  do_one_256_swizzle_32(x[6],   x[7]);

  do_one_256_swizzle_64(x[0],   x[2]);
  do_one_256_swizzle_64(x[1],   x[3]);
  do_one_256_swizzle_64(x[4],   x[6]);
  do_one_256_swizzle_64(x[5],   x[7]);

  do_one_256_swizzle_128(x[0],  x[4]);
  do_one_256_swizzle_128(x[1],  x[5]);
  do_one_256_swizzle_128(x[2],  x[6]);
  do_one_256_swizzle_128(x[3],  x[7]);

  std::swap(x[1],  x[2]);
  std::swap(x[5],  x[6]);
}

void calin::util::vcl::transpose(Vec8ui* x)
{
  do_one_256_swizzle_32(x[0],   x[1]);
  do_one_256_swizzle_32(x[2],   x[3]);
  do_one_256_swizzle_32(x[4],   x[5]);
  do_one_256_swizzle_32(x[6],   x[7]);

  do_one_256_swizzle_64(x[0],   x[2]);
  do_one_256_swizzle_64(x[1],   x[3]);
  do_one_256_swizzle_64(x[4],   x[6]);
  do_one_256_swizzle_64(x[5],   x[7]);

  do_one_256_swizzle_128(x[0],  x[4]);
  do_one_256_swizzle_128(x[1],  x[5]);
  do_one_256_swizzle_128(x[2],  x[6]);
  do_one_256_swizzle_128(x[3],  x[7]);

  std::swap(x[1],  x[2]);
  std::swap(x[5],  x[6]);
}

void calin::util::vcl::transpose(Vec4q* x)
{
  do_one_256_swizzle_64(x[0],   x[1]);
  do_one_256_swizzle_64(x[2],   x[3]);

  do_one_256_swizzle_128(x[0],  x[2]);
  do_one_256_swizzle_128(x[1],  x[3]);
}

void calin::util::vcl::transpose(Vec4uq* x)
{
  do_one_256_swizzle_64(x[0],   x[1]);
  do_one_256_swizzle_64(x[2],   x[3]);

  do_one_256_swizzle_128(x[0],  x[2]);
  do_one_256_swizzle_128(x[1],  x[3]);
}

void calin::util::vcl::transpose(Vec8f* x)
{
  do_one_256_swizzle_ps(x[0],   x[1]);
  do_one_256_swizzle_ps(x[2],   x[3]);
  do_one_256_swizzle_ps(x[4],   x[5]);
  do_one_256_swizzle_ps(x[6],   x[7]);

  do_one_256_swizzle_pd(x[0],   x[2]);
  do_one_256_swizzle_pd(x[1],   x[3]);
  do_one_256_swizzle_pd(x[4],   x[6]);
  do_one_256_swizzle_pd(x[5],   x[7]);

  do_one_256_swizzle_flt128(x[0],  x[4]);
  do_one_256_swizzle_flt128(x[1],  x[5]);
  do_one_256_swizzle_flt128(x[2],  x[6]);
  do_one_256_swizzle_flt128(x[3],  x[7]);

  std::swap(x[1],  x[2]);
  std::swap(x[5],  x[6]);
}

void calin::util::vcl::transpose(Vec4d* x)
{
  do_one_256_swizzle_pd(x[0],   x[1]);
  do_one_256_swizzle_pd(x[2],   x[3]);

  do_one_256_swizzle_flt128(x[0],  x[2]);
  do_one_256_swizzle_flt128(x[1],  x[3]);
}

#endif // MAX_VECTOR_SIZE >= 256

#if MAX_VECTOR_SIZE >= 512

namespace {
#if INSTRSET >= 9
  template<typename T> inline void do_one_512_swizzle_16(T& a, T& b) {
    __m512i tmp = _mm512_unpackhi_epi16(a, b);
    a = _mm512_unpacklo_epi16(a, b);
    b = tmp;
  }

  template<typename T> inline void do_one_512_swizzle_32(T& a, T& b) {
    __m512i tmp = _mm512_unpackhi_epi32(a, b);
    a = _mm512_unpacklo_epi32(a, b);
    b = tmp;
  }

  template<typename T> inline void do_one_512_swizzle_64(T& a, T& b) {
    __m512i tmp = _mm512_unpackhi_epi64(a, b);
    a = _mm512_unpacklo_epi64(a, b);
    b = tmp;
  }

  template<typename T> inline void do_one_512_swizzle_128(T& a, T& b) {
    __m512i tmp =  _mm512_shuffle_i32x4(a, b, 0xdd);
    a = _mm512_shuffle_i32x4(a, b, 0x88);
    b = tmp;
  }

  // template<typename T> inline void do_one_512_swizzle_ps(T& a, T& b) {
  //   T tmp = _mm512_unpackhi_ps(a, b);
  //   a = _mm512_unpacklo_ps(a, b);
  //   b = tmp;
  // }
  //
  // inline void do_one_512_swizzle_pd(Vec16f& a, Vec16f& b) {
  //   __m512d tmp = _mm512_unpackhi_pd(vcl::reinterpret_d(a), vcl::reinterpret_d(b));
  //   a = _mm512_castpd_ps(_mm512_unpacklo_pd(vcl::reinterpret_d(a), vcl::reinterpret_d(b)));
  //   b = _mm512_castpd_ps(tmp);
  // }
  //
  // inline void do_one_512_swizzle_pd(Vec8d& a, Vec8d& b) {
  //   __m512d tmp = _mm512_unpackhi_pd(a, b);
  //   a = _mm512_unpacklo_pd(a, b);
  //   b = tmp;
  // }
  //
  // inline void do_one_512_swizzle_flt128(Vec16f& a, Vec16f& b) {
  //   __m512 tmp =  _mm512_shuffle_f32x4(a, b, 0xdd);
  //   a = _mm512_shuffle_f32x4(a, b, 0x88);
  //   b = tmp;
  // }
  //
  // inline void do_one_512_swizzle_flt128(Vec8d& a, Vec8d& b) {
  //   __m512 tmp =  _mm512_shuffle_d64x2(a, b, 0xdd);
  //   a = _mm512_shuffle_d64x2(a, b, 0x88);
  //   b = tmp;
  // }

#else // INSTRSET >= 9
  template<typename T> inline void do_one_512_swizzle_16(T& a, T& b) {
    auto a_lo = a.get_low();
    auto b_lo = b.get_low();
    auto a_hi = a.get_high();
    auto b_hi = b.get_high();
    do_one_256_swizzle_16(a_lo,b_lo);
    do_one_256_swizzle_16(a_hi,b_hi);
    a = T(a_lo, a_hi);
    b = T(b_lo, b_hi);
  }

  template<typename T> inline void do_one_512_swizzle_32(T& a, T& b) {
    auto a_lo = a.get_low();
    auto b_lo = b.get_low();
    auto a_hi = a.get_high();
    auto b_hi = b.get_high();
    do_one_256_swizzle_32(a_lo,b_lo);
    do_one_256_swizzle_32(a_hi,b_hi);
    a = T(a_lo, a_hi);
    b = T(b_lo, b_hi);
  }

  template<typename T> inline void do_one_512_swizzle_64(T& a, T& b) {
    auto a_lo = a.get_low();
    auto b_lo = b.get_low();
    auto a_hi = a.get_high();
    auto b_hi = b.get_high();
    do_one_256_swizzle_64(a_lo,b_lo);
    do_one_256_swizzle_64(a_hi,b_hi);
    a = T(a_lo, a_hi);
    b = T(b_lo, b_hi);
  }

  template<typename T> inline void do_one_512_swizzle_128(T& a, T& b) {
    auto a_lo = a.get_low();
    auto b_lo = b.get_low();
    auto a_hi = a.get_high();
    auto b_hi = b.get_high();
    do_one_256_swizzle_128(a_lo,a_hi);
    do_one_256_swizzle_128(b_lo,b_hi);
    a = T(a_lo, b_lo);
    b = T(a_hi, b_hi);
  }

  // inline void do_one_512_swizzle_ps(Vec8f& a, Vec8f& b) {
  //   Vec8f tmp(_mm_unpackhi_ps(a.get_low(), b.get_low()),
  //         _mm_unpackhi_ps(a.get_high(), b.get_high()));
  //   a = Vec8f(_mm_unpacklo_ps(a.get_low(), b.get_low()),
  //         _mm_unpacklo_ps(a.get_high(), b.get_high()));
  //   b = tmp;
  // }
  //
  // inline void do_one_512_swizzle_pd(Vec8f& a, Vec8f& b) {
  //   Vec8f tmp(_mm_castpd_ps(_mm_unpackhi_pd(vcl::reinterpret_d(a.get_low()), vcl::reinterpret_d(b.get_low()))),
  //             _mm_castpd_ps(_mm_unpackhi_pd(vcl::reinterpret_d(a.get_high()), vcl::reinterpret_d(b.get_high()))));
  //   a = Vec8f(_mm_castpd_ps(_mm_unpacklo_pd(vcl::reinterpret_d(a.get_low()), vcl::reinterpret_d(b.get_low()))),
  //             _mm_castpd_ps(_mm_unpacklo_pd(vcl::reinterpret_d(a.get_high()), vcl::reinterpret_d(b.get_high()))));
  //   b = tmp;
  // }
  //
  // inline void do_one_512_swizzle_pd(Vec4d& a, Vec4d& b) {
  //   Vec4d tmp(_mm_unpackhi_pd(a.get_low(), b.get_low()),
  //             _mm_unpackhi_pd(a.get_high(), b.get_high()));
  //   a = Vec4d(_mm_unpacklo_pd(a.get_low(), b.get_low()),
  //             _mm_unpacklo_pd(a.get_high(), b.get_high()));
  //   b = tmp;
  // }
  //
  // template<typename T> inline void do_one_512_swizzle_flt128(T& a, T& b) {
  //   T tmp(a.get_high(), b.get_high());
  //   a = T(a.get_low(), b.get_low());
  //   b = tmp;
  // }
#endif // INSTRSET >= 9
}

void calin::util::vcl::transpose(Vec32s* x)
{
  do_one_512_swizzle_16(x[ 0], x[ 1]);
  do_one_512_swizzle_16(x[ 2], x[ 3]);
  do_one_512_swizzle_16(x[ 4], x[ 5]);
  do_one_512_swizzle_16(x[ 6], x[ 7]);
  do_one_512_swizzle_16(x[ 8], x[ 9]);
  do_one_512_swizzle_16(x[10], x[11]);
  do_one_512_swizzle_16(x[12], x[13]);
  do_one_512_swizzle_16(x[14], x[15]);
  do_one_512_swizzle_16(x[16], x[17]);
  do_one_512_swizzle_16(x[18], x[19]);
  do_one_512_swizzle_16(x[20], x[21]);
  do_one_512_swizzle_16(x[22], x[23]);
  do_one_512_swizzle_16(x[24], x[25]);
  do_one_512_swizzle_16(x[26], x[27]);
  do_one_512_swizzle_16(x[28], x[29]);
  do_one_512_swizzle_16(x[30], x[31]);

  do_one_512_swizzle_32(x[ 0], x[ 2]);
  do_one_512_swizzle_32(x[ 1], x[ 3]);
  do_one_512_swizzle_32(x[ 4], x[ 6]);
  do_one_512_swizzle_32(x[ 5], x[ 7]);
  do_one_512_swizzle_32(x[ 8], x[10]);
  do_one_512_swizzle_32(x[ 9], x[11]);
  do_one_512_swizzle_32(x[12], x[14]);
  do_one_512_swizzle_32(x[13], x[15]);
  do_one_512_swizzle_32(x[16], x[18]);
  do_one_512_swizzle_32(x[17], x[19]);
  do_one_512_swizzle_32(x[20], x[22]);
  do_one_512_swizzle_32(x[21], x[23]);
  do_one_512_swizzle_32(x[24], x[26]);
  do_one_512_swizzle_32(x[25], x[27]);
  do_one_512_swizzle_32(x[28], x[30]);
  do_one_512_swizzle_32(x[29], x[31]);

  do_one_512_swizzle_64(x[ 0], x[ 4]);
  do_one_512_swizzle_64(x[ 1], x[ 5]);
  do_one_512_swizzle_64(x[ 2], x[ 6]);
  do_one_512_swizzle_64(x[ 3], x[ 7]);
  do_one_512_swizzle_64(x[ 8], x[12]);
  do_one_512_swizzle_64(x[ 9], x[13]);
  do_one_512_swizzle_64(x[10], x[14]);
  do_one_512_swizzle_64(x[11], x[15]);
  do_one_512_swizzle_64(x[16], x[20]);
  do_one_512_swizzle_64(x[17], x[21]);
  do_one_512_swizzle_64(x[18], x[22]);
  do_one_512_swizzle_64(x[19], x[23]);
  do_one_512_swizzle_64(x[24], x[28]);
  do_one_512_swizzle_64(x[25], x[29]);
  do_one_512_swizzle_64(x[26], x[30]);
  do_one_512_swizzle_64(x[27], x[31]);

  do_one_512_swizzle_128(x[ 0], x[ 8]);
  do_one_512_swizzle_128(x[ 1], x[ 9]);
  do_one_512_swizzle_128(x[ 2], x[10]);
  do_one_512_swizzle_128(x[ 3], x[11]);
  do_one_512_swizzle_128(x[ 4], x[12]);
  do_one_512_swizzle_128(x[ 5], x[13]);
  do_one_512_swizzle_128(x[ 6], x[14]);
  do_one_512_swizzle_128(x[ 7], x[15]);
  do_one_512_swizzle_128(x[16], x[24]);
  do_one_512_swizzle_128(x[17], x[25]);
  do_one_512_swizzle_128(x[18], x[26]);
  do_one_512_swizzle_128(x[19], x[27]);
  do_one_512_swizzle_128(x[20], x[28]);
  do_one_512_swizzle_128(x[21], x[29]);
  do_one_512_swizzle_128(x[22], x[30]);
  do_one_512_swizzle_128(x[23], x[31]);

  do_one_512_swizzle_128(x[ 0], x[16]);
  do_one_512_swizzle_128(x[ 1], x[17]);
  do_one_512_swizzle_128(x[ 2], x[18]);
  do_one_512_swizzle_128(x[ 3], x[19]);
  do_one_512_swizzle_128(x[ 4], x[20]);
  do_one_512_swizzle_128(x[ 5], x[21]);
  do_one_512_swizzle_128(x[ 6], x[22]);
  do_one_512_swizzle_128(x[ 7], x[23]);
  do_one_512_swizzle_128(x[ 8], x[24]);
  do_one_512_swizzle_128(x[ 9], x[25]);
  do_one_512_swizzle_128(x[10], x[26]);
  do_one_512_swizzle_128(x[11], x[27]);
  do_one_512_swizzle_128(x[12], x[28]);
  do_one_512_swizzle_128(x[13], x[29]);
  do_one_512_swizzle_128(x[14], x[30]);
  do_one_512_swizzle_128(x[15], x[31]);

  std::swap(x[ 1], x[ 4]);
  std::swap(x[ 3], x[ 6]);
  std::swap(x[ 9], x[12]);
  std::swap(x[11], x[14]);
  std::swap(x[17], x[20]);
  std::swap(x[19], x[22]);
  std::swap(x[25], x[28]);
  std::swap(x[27], x[30]);
}

void calin::util::vcl::transpose(Vec32us* x)
{
  do_one_512_swizzle_16(x[ 0], x[ 1]);
  do_one_512_swizzle_16(x[ 2], x[ 3]);
  do_one_512_swizzle_16(x[ 4], x[ 5]);
  do_one_512_swizzle_16(x[ 6], x[ 7]);
  do_one_512_swizzle_16(x[ 8], x[ 9]);
  do_one_512_swizzle_16(x[10], x[11]);
  do_one_512_swizzle_16(x[12], x[13]);
  do_one_512_swizzle_16(x[14], x[15]);
  do_one_512_swizzle_16(x[16], x[17]);
  do_one_512_swizzle_16(x[18], x[19]);
  do_one_512_swizzle_16(x[20], x[21]);
  do_one_512_swizzle_16(x[22], x[23]);
  do_one_512_swizzle_16(x[24], x[25]);
  do_one_512_swizzle_16(x[26], x[27]);
  do_one_512_swizzle_16(x[28], x[29]);
  do_one_512_swizzle_16(x[30], x[31]);

  do_one_512_swizzle_32(x[ 0], x[ 2]);
  do_one_512_swizzle_32(x[ 1], x[ 3]);
  do_one_512_swizzle_32(x[ 4], x[ 6]);
  do_one_512_swizzle_32(x[ 5], x[ 7]);
  do_one_512_swizzle_32(x[ 8], x[10]);
  do_one_512_swizzle_32(x[ 9], x[11]);
  do_one_512_swizzle_32(x[12], x[14]);
  do_one_512_swizzle_32(x[13], x[15]);
  do_one_512_swizzle_32(x[16], x[18]);
  do_one_512_swizzle_32(x[17], x[19]);
  do_one_512_swizzle_32(x[20], x[22]);
  do_one_512_swizzle_32(x[21], x[23]);
  do_one_512_swizzle_32(x[24], x[26]);
  do_one_512_swizzle_32(x[25], x[27]);
  do_one_512_swizzle_32(x[28], x[30]);
  do_one_512_swizzle_32(x[29], x[31]);

  do_one_512_swizzle_64(x[ 0], x[ 4]);
  do_one_512_swizzle_64(x[ 1], x[ 5]);
  do_one_512_swizzle_64(x[ 2], x[ 6]);
  do_one_512_swizzle_64(x[ 3], x[ 7]);
  do_one_512_swizzle_64(x[ 8], x[12]);
  do_one_512_swizzle_64(x[ 9], x[13]);
  do_one_512_swizzle_64(x[10], x[14]);
  do_one_512_swizzle_64(x[11], x[15]);
  do_one_512_swizzle_64(x[16], x[20]);
  do_one_512_swizzle_64(x[17], x[21]);
  do_one_512_swizzle_64(x[18], x[22]);
  do_one_512_swizzle_64(x[19], x[23]);
  do_one_512_swizzle_64(x[24], x[28]);
  do_one_512_swizzle_64(x[25], x[29]);
  do_one_512_swizzle_64(x[26], x[30]);
  do_one_512_swizzle_64(x[27], x[31]);

  do_one_512_swizzle_128(x[ 0], x[ 8]);
  do_one_512_swizzle_128(x[ 1], x[ 9]);
  do_one_512_swizzle_128(x[ 2], x[10]);
  do_one_512_swizzle_128(x[ 3], x[11]);
  do_one_512_swizzle_128(x[ 4], x[12]);
  do_one_512_swizzle_128(x[ 5], x[13]);
  do_one_512_swizzle_128(x[ 6], x[14]);
  do_one_512_swizzle_128(x[ 7], x[15]);
  do_one_512_swizzle_128(x[16], x[24]);
  do_one_512_swizzle_128(x[17], x[25]);
  do_one_512_swizzle_128(x[18], x[26]);
  do_one_512_swizzle_128(x[19], x[27]);
  do_one_512_swizzle_128(x[20], x[28]);
  do_one_512_swizzle_128(x[21], x[29]);
  do_one_512_swizzle_128(x[22], x[30]);
  do_one_512_swizzle_128(x[23], x[31]);

  do_one_512_swizzle_128(x[ 0], x[16]);
  do_one_512_swizzle_128(x[ 1], x[17]);
  do_one_512_swizzle_128(x[ 2], x[18]);
  do_one_512_swizzle_128(x[ 3], x[19]);
  do_one_512_swizzle_128(x[ 4], x[20]);
  do_one_512_swizzle_128(x[ 5], x[21]);
  do_one_512_swizzle_128(x[ 6], x[22]);
  do_one_512_swizzle_128(x[ 7], x[23]);
  do_one_512_swizzle_128(x[ 8], x[24]);
  do_one_512_swizzle_128(x[ 9], x[25]);
  do_one_512_swizzle_128(x[10], x[26]);
  do_one_512_swizzle_128(x[11], x[27]);
  do_one_512_swizzle_128(x[12], x[28]);
  do_one_512_swizzle_128(x[13], x[29]);
  do_one_512_swizzle_128(x[14], x[30]);
  do_one_512_swizzle_128(x[15], x[31]);

  std::swap(x[ 1], x[ 4]);
  std::swap(x[ 3], x[ 6]);
  std::swap(x[ 9], x[12]);
  std::swap(x[11], x[14]);
  std::swap(x[17], x[20]);
  std::swap(x[19], x[22]);
  std::swap(x[25], x[28]);
  std::swap(x[27], x[30]);
}

// void transpose(Vec16i* x);
// void transpose(Vec16ui* x);
// void transpose(Vec8q* x);
// void transpose(Vec8uq* x);
// void transpose(Vec16f* x);
// void transpose(Vec8d* x);

#endif // MAX_VECTOR_SIZE >= 512
