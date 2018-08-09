/*

   calin/util/vcl.hpp -- Stephen Fegan -- 2018-08-08

   Calin interface to Agner Fog's Vector Class Library (VCL)

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

#include <ostream>

#define MAX_VECTOR_SIZE 512
#define VCL_NAMESPACE vcl

#include <VCL/vectorclass.h>
#include <VCL/vectormath_trig.h>
#include <VCL/vectormath_exp.h>
#include <VCL/vectormath_hyp.h>

#include <Eigen/Core>

namespace calin { namespace util { namespace vcl {

  using namespace ::vcl;

  struct VCL128Architecture
  {
    constexpr static unsigned vec_size   = 128;
    constexpr static unsigned num_int8   = vec_size/sizeof(int8_t);
    constexpr static unsigned num_uint8  = vec_size/sizeof(uint8_t);
    constexpr static unsigned num_int16  = vec_size/sizeof(int16_t);
    constexpr static unsigned num_uint16 = vec_size/sizeof(uint16_t);
    constexpr static unsigned num_int32  = vec_size/sizeof(int32_t);
    constexpr static unsigned num_uint32 = vec_size/sizeof(uint32_t);
    constexpr static unsigned num_int64  = vec_size/sizeof(int64_t);
    constexpr static unsigned num_uint64 = vec_size/sizeof(uint64_t);
    constexpr static unsigned num_float  = vec_size/sizeof(float);
    constexpr static unsigned num_double = vec_size/sizeof(double);

    typedef Vec128b bool_vec_type;
    typedef Vec16c  int8_vec_type;
    typedef Vec16uc uint8_vec_type;
    typedef Vec8s   int16_vec_type;
    typedef Vec8us  uint16_vec_type;
    typedef Vec4i   int32_vec_type;
    typedef Vec4ui  uint32_vec_type;
    typedef Vec2q   int64_vec_type;
    typedef Vec2uq  uint64_vec_type;
    typedef Vec4f   float_vec_type;
    typedef Vec2d   double_vec_type;

    // template<int32_t x> static int32_type constant_int32() {
    //   return vcl::constant4i<x,x,x,x>(); }
    // template<uint32_t x> static uint32_type constant_uint32() {
    //   return vcl::constant4ui<x,x,x,x>(); }
    // template<int64_t x> static int64_type constant_int64() {
    //   return vcl::constant2i<x,x>(); }
    // template<uint64_t x> static uint64_type constant_uint64() {
    //   return vcl::constant2ui<x,x>(); }
    // template<float x> static float_type constant_float() {
    //   return vcl::constant4f<x,x,x,x>(); }
    // template<double x> static double_type constant_double() {
    //   return vcl::constant2d<x,x>(); }
  };

  struct VCL256Architecture
  {
    constexpr static unsigned vec_size   = 256;
    constexpr static unsigned num_int8   = vec_size/sizeof(int8_t);
    constexpr static unsigned num_uint8  = vec_size/sizeof(uint8_t);
    constexpr static unsigned num_int16  = vec_size/sizeof(int16_t);
    constexpr static unsigned num_uint16 = vec_size/sizeof(uint16_t);
    constexpr static unsigned num_int32  = vec_size/sizeof(int32_t);
    constexpr static unsigned num_uint32 = vec_size/sizeof(uint32_t);
    constexpr static unsigned num_int64  = vec_size/sizeof(int64_t);
    constexpr static unsigned num_uint64 = vec_size/sizeof(uint64_t);
    constexpr static unsigned num_float  = vec_size/sizeof(float);
    constexpr static unsigned num_double = vec_size/sizeof(double);

    typedef Vec256b bool_vec_type;
    typedef Vec32c  int8_vec_type;
    typedef Vec32uc uint8_vec_type;
    typedef Vec16s  int16_vec_type;
    typedef Vec16us uint16_vec_type;
    typedef Vec8i   int32_vec_type;
    typedef Vec8ui  uint32_vec_type;
    typedef Vec4q   int64_vec_type;
    typedef Vec4uq  uint64_vec_type;
    typedef Vec8f   float_vec_type;
    typedef Vec4d   double_vec_type;

    // template<int32_t x> static int32_type constant_int32() {
    //   return vcl::constant8i<x,x,x,x,x,x,x,x>(); }
    // template<uint32_t x> static uint32_type constant_uint32() {
    //   return vcl::constant8ui<x,x,x,x,x,x,x,x>(); }
    // template<int64_t x> static int64_type constant_int64() {
    //   return vcl::constant4i<x,x,x,x>(); }
    // template<uint64_t x> static uint64_type constant_uint64() {
    //   return vcl::constant4ui<x,x,x,x>(); }
    // template<float x> static float_type constant_float() {
    //   return vcl::constant8f<x,x,x,x,x,x,x,x>(); }
    // template<double x> static double_type constant_double() {
    //   return vcl::constant4d<x,x,x,x>(); }
  };

  struct VCL512Architecture
  {
    constexpr static unsigned vec_size   = 512;
    constexpr static unsigned num_int8   = vec_size/sizeof(int8_t);
    constexpr static unsigned num_uint8  = vec_size/sizeof(uint8_t);
    constexpr static unsigned num_int16  = vec_size/sizeof(int16_t);
    constexpr static unsigned num_uint16 = vec_size/sizeof(uint16_t);
    constexpr static unsigned num_int32  = vec_size/sizeof(int32_t);
    constexpr static unsigned num_uint32 = vec_size/sizeof(uint32_t);
    constexpr static unsigned num_int64  = vec_size/sizeof(int64_t);
    constexpr static unsigned num_uint64 = vec_size/sizeof(uint64_t);
    constexpr static unsigned num_float  = vec_size/sizeof(float);
    constexpr static unsigned num_double = vec_size/sizeof(double);

    // typedef Vec512b bool_vec_type;
    // typedef Vec64c  int8_vec_type;
    // typedef Vec64uc uint8_vec_type;
    // typedef Vec32s  int16_vec_type;
    // typedef Vec32us uint16_vec_type;
    typedef void    bool_vec_type;
    typedef void    int8_vec_type;
    typedef void    uint8_vec_type;
    typedef void    int16_vec_type;
    typedef void    uint16_vec_type;
    typedef Vec16i  int32_vec_type;
    typedef Vec16ui uint32_vec_type;
    typedef Vec8q   int64_vec_type;
    typedef Vec8uq  uint64_vec_type;
    typedef Vec16f  float_vec_type;
    typedef Vec8d   double_vec_type;
  };

  template<typename T> struct vcl_type
  {
    typedef void scalar_type;
    typedef void vcl_architecture;
  };

#define DEFINE_VCL_TYPE(Vec,Scalar,Arch) \
  template<> struct vcl_type<Vec> \
  { \
    typedef Scalar scalar_type; \
    typedef Arch vcl_architecture; \
  }

  DEFINE_VCL_TYPE(Vec16c,  int8_t,   VCL128Architecture);
  DEFINE_VCL_TYPE(Vec16uc, uint8_t,  VCL128Architecture);
  DEFINE_VCL_TYPE(Vec8s,   int16_t,  VCL128Architecture);
  DEFINE_VCL_TYPE(Vec8us,  uint16_t, VCL128Architecture);
  DEFINE_VCL_TYPE(Vec4i,   int32_t,  VCL128Architecture);
  DEFINE_VCL_TYPE(Vec4ui,  uint32_t, VCL128Architecture);
  DEFINE_VCL_TYPE(Vec2q,   int64_t,  VCL128Architecture);
  DEFINE_VCL_TYPE(Vec2uq,  uint64_t, VCL128Architecture);
  DEFINE_VCL_TYPE(Vec4f,   float,    VCL128Architecture);
  DEFINE_VCL_TYPE(Vec2d,   double,   VCL128Architecture);

  DEFINE_VCL_TYPE(Vec32c,  int8_t,   VCL256Architecture);
  DEFINE_VCL_TYPE(Vec32uc, uint8_t,  VCL256Architecture);
  DEFINE_VCL_TYPE(Vec16s,  int16_t,  VCL256Architecture);
  DEFINE_VCL_TYPE(Vec16us, uint16_t, VCL256Architecture);
  DEFINE_VCL_TYPE(Vec8i,   int32_t,  VCL256Architecture);
  DEFINE_VCL_TYPE(Vec8ui,  uint32_t, VCL256Architecture);
  DEFINE_VCL_TYPE(Vec4q,   int64_t,  VCL256Architecture);
  DEFINE_VCL_TYPE(Vec4uq,  uint64_t, VCL256Architecture);
  DEFINE_VCL_TYPE(Vec8f,   float,    VCL256Architecture);
  DEFINE_VCL_TYPE(Vec4d,   double,   VCL256Architecture);

  // DEFINE_VCL_TYPE(Vec64c,  int8_t,   VCL512Architecture);
  // DEFINE_VCL_TYPE(Vec64uc, uint8_t,  VCL512Architecture);
  // DEFINE_VCL_TYPE(Vec32s,  int16_t,  VCL512Architecture);
  // DEFINE_VCL_TYPE(Vec32us, uint16_t, VCL512Architecture);
  DEFINE_VCL_TYPE(Vec16i,  int32_t,  VCL512Architecture);
  DEFINE_VCL_TYPE(Vec16ui, uint32_t, VCL512Architecture);
  DEFINE_VCL_TYPE(Vec8q,   int64_t,  VCL512Architecture);
  DEFINE_VCL_TYPE(Vec8uq,  uint64_t, VCL512Architecture);
  DEFINE_VCL_TYPE(Vec16f,  float,    VCL512Architecture);
  DEFINE_VCL_TYPE(Vec8d,   double,   VCL512Architecture);

#undef DEFINE_VCL_TYPE

  template<typename Vec> void print_vec(std::ostream& s, const Vec& v)
  {
    s << v[0];
    for(unsigned i=1;i<v.size();i++)
      s << '_' << v[i];
  }

  inline Vec2uq multiply_low_32bit(const Vec2uq& a, const Vec2uq& b) {
    return _mm_mul_epu32(a, b);
  }

#if MAX_VECTOR_SIZE >= 256
#if INSTRSET >= 8
  inline Vec4uq multiply_low_32bit(const Vec4uq& a, const Vec4uq& b) {
    return _mm256_mul_epu32(a, b);
  }
#else // INSTRSET < 8
  inline Vec4uq multiply_low_32bit(const Vec4uq& a, const Vec4uq& b) {
    return Vec4uq(multiply_low_32bit(a.get_low(), b.get_low()),
                  multiply_low_32bit(a.get_high(), b.get_high()));
  }
#endif // INSTRSET < 8
#endif // MAX_VECTOR_SIZE >= 256

#if MAX_VECTOR_SIZE >= 512
#if INSTRSET >= 9
  inline Vec8uq multiply_low_32bit(const Vec8uq& a, const Vec8uq& b) {
    return _mm512_mul_epu32(a, b);
  }
#else // INSTRSET < 9
  inline Vec8uq multiply_low_32bit(const Vec8uq& a, const Vec8uq& b) {
    return Vec8uq(multiply_low_32bit(a.get_low(), b.get_low()),
                  multiply_low_32bit(a.get_high(), b.get_high()));
  }
#endif // INSTRSET < 9
#endif // MAX_VECTOR_SIZE >= 512

} } } // namespace calin::util::vcl

#define ADD_OSTREAM_OPERATOR(Vec) \
  std::ostream& operator<<(std::ostream& s, const Vec& v) { \
    calin::util::vcl::print_vec(s, v); \
    return s; \
  }

ADD_OSTREAM_OPERATOR(calin::util::vcl::Vec16c);
ADD_OSTREAM_OPERATOR(calin::util::vcl::Vec16uc);
ADD_OSTREAM_OPERATOR(calin::util::vcl::Vec8s);
ADD_OSTREAM_OPERATOR(calin::util::vcl::Vec8us);
ADD_OSTREAM_OPERATOR(calin::util::vcl::Vec4i);
ADD_OSTREAM_OPERATOR(calin::util::vcl::Vec4ui);
ADD_OSTREAM_OPERATOR(calin::util::vcl::Vec2q);
ADD_OSTREAM_OPERATOR(calin::util::vcl::Vec2uq);
ADD_OSTREAM_OPERATOR(calin::util::vcl::Vec4f);
ADD_OSTREAM_OPERATOR(calin::util::vcl::Vec2d);

ADD_OSTREAM_OPERATOR(calin::util::vcl::Vec32c);
ADD_OSTREAM_OPERATOR(calin::util::vcl::Vec32uc);
ADD_OSTREAM_OPERATOR(calin::util::vcl::Vec16s);
ADD_OSTREAM_OPERATOR(calin::util::vcl::Vec16us);
ADD_OSTREAM_OPERATOR(calin::util::vcl::Vec8i);
ADD_OSTREAM_OPERATOR(calin::util::vcl::Vec8ui);
ADD_OSTREAM_OPERATOR(calin::util::vcl::Vec4q);
ADD_OSTREAM_OPERATOR(calin::util::vcl::Vec4uq);
ADD_OSTREAM_OPERATOR(calin::util::vcl::Vec8f);
ADD_OSTREAM_OPERATOR(calin::util::vcl::Vec4d);

// ADD_OSTREAM_OPERATOR(calin::util::vcl::Vec64c);
// ADD_OSTREAM_OPERATOR(calin::util::vcl::Vec64uc);
// ADD_OSTREAM_OPERATOR(calin::util::vcl::Vec32s);
// ADD_OSTREAM_OPERATOR(calin::util::vcl::Vec32us);
ADD_OSTREAM_OPERATOR(calin::util::vcl::Vec16i);
ADD_OSTREAM_OPERATOR(calin::util::vcl::Vec16ui);
ADD_OSTREAM_OPERATOR(calin::util::vcl::Vec8q);
ADD_OSTREAM_OPERATOR(calin::util::vcl::Vec8uq);
ADD_OSTREAM_OPERATOR(calin::util::vcl::Vec16f);
ADD_OSTREAM_OPERATOR(calin::util::vcl::Vec8d);

#undef ADD_OSTREAM_OPERATOR

namespace calin { namespace util { namespace vcl {



} } } // namespace calin::util::vcl
