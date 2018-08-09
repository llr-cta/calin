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

#define VCL_NAMESPACE vcl
#include <VCL/vectorclass.h>
#include <VCL/vectormath_trig.h>
#include <VCL/vectormath_exp.h>
#include <VCL/vectormath_hyp.h>

#include <Eigen/Core>

namespace calin { namespace util { namespace vcl {

  using namespace ::vcl;

  template<unsigned VEC_SIZE> struct vcl_architecture
  {
    constexpr static unsigned vec_size   = VEC_SIZE;
    constexpr static unsigned num_int8   = VEC_SIZE/sizeof(int8_t);
    constexpr static unsigned num_uint8  = VEC_SIZE/sizeof(uint8_t);
    constexpr static unsigned num_int16  = VEC_SIZE/sizeof(int16_t);
    constexpr static unsigned num_uint16 = VEC_SIZE/sizeof(uint16_t);
    constexpr static unsigned num_int32  = VEC_SIZE/sizeof(int32_t);
    constexpr static unsigned num_uint32 = VEC_SIZE/sizeof(uint32_t);
    constexpr static unsigned num_int64  = VEC_SIZE/sizeof(int64_t);
    constexpr static unsigned num_uint64 = VEC_SIZE/sizeof(uint64_t);
    constexpr static unsigned num_float  = VEC_SIZE/sizeof(float);
    constexpr static unsigned num_double = VEC_SIZE/sizeof(double);

    typedef void bool_type;
    typedef void int8_type;
    typedef void uint8_type;
    typedef void int16_type;
    typedef void uint16_type;
    typedef void int32_type;
    typedef void uint32_type;
    typedef void int64_type;
    typedef void uint64_type;
    typedef void float_type;
    typedef void double_type;

    // template<int32_t x> static int32_type constant_int32() { }
    // template<uint32_t x> static uint32_type constant_uint32() { }
    // template<int64_t x> static int64_type constant_int64() { }
    // template<uint64_t x> static uint64_type constant_uint64() { }
    // template<float x> static float_type constant_float() { }
    // template<double x> static double_type constant_double() { }
  };

  template<> struct vcl_architecture<128>
  {
    constexpr static unsigned vec_size   = 128;
    constexpr static unsigned num_int8   = 128/sizeof(int8_t);
    constexpr static unsigned num_uint8  = 128/sizeof(uint8_t);
    constexpr static unsigned num_int16  = 128/sizeof(int16_t);
    constexpr static unsigned num_uint16 = 128/sizeof(uint16_t);
    constexpr static unsigned num_int32  = 128/sizeof(int32_t);
    constexpr static unsigned num_uint32 = 128/sizeof(uint32_t);
    constexpr static unsigned num_int64  = 128/sizeof(int64_t);
    constexpr static unsigned num_uint64 = 128/sizeof(uint64_t);
    constexpr static unsigned num_float  = 128/sizeof(float);
    constexpr static unsigned num_double = 128/sizeof(double);

    typedef Vec128b bool_type;
    typedef Vec16c  int8_type;
    typedef Vec16uc uint8_type;
    typedef Vec8s   int16_type;
    typedef Vec8us  uint16_type;
    typedef Vec4i   int32_type;
    typedef Vec4ui  uint32_type;
    typedef Vec2q   int64_type;
    typedef Vec2uq  uint64_type;
    typedef Vec4f   float_type;
    typedef Vec2d   double_type;

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

  template<> struct vcl_architecture<256>
  {
    constexpr static unsigned vec_size   = 256;
    constexpr static unsigned num_int8   = 256/sizeof(int8_t);
    constexpr static unsigned num_uint8  = 256/sizeof(uint8_t);
    constexpr static unsigned num_int16  = 256/sizeof(int16_t);
    constexpr static unsigned num_uint16 = 256/sizeof(uint16_t);
    constexpr static unsigned num_int32  = 256/sizeof(int32_t);
    constexpr static unsigned num_uint32 = 256/sizeof(uint32_t);
    constexpr static unsigned num_int64  = 256/sizeof(int64_t);
    constexpr static unsigned num_uint64 = 256/sizeof(uint64_t);
    constexpr static unsigned num_float  = 256/sizeof(float);
    constexpr static unsigned num_double = 256/sizeof(double);

    typedef Vec256b bool_type;
    typedef Vec32c  int8_type;
    typedef Vec32uc uint8_type;
    typedef Vec16s  int16_type;
    typedef Vec16us uint16_type;
    typedef Vec8i   int32_type;
    typedef Vec8ui  uint32_type;
    typedef Vec4q   int64_type;
    typedef Vec4uq  uint64_type;
    typedef Vec8f   float_type;
    typedef Vec4d   double_type;

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

  DEFINE_VCL_TYPE(Vec16c,  int8_t,   vcl_architecture<128>);
  DEFINE_VCL_TYPE(Vec16uc, uint8_t,  vcl_architecture<128>);
  DEFINE_VCL_TYPE(Vec8s,   int16_t,  vcl_architecture<128>);
  DEFINE_VCL_TYPE(Vec8us,  uint16_t, vcl_architecture<128>);
  DEFINE_VCL_TYPE(Vec4i,   int32_t,  vcl_architecture<128>);
  DEFINE_VCL_TYPE(Vec4ui,  uint32_t, vcl_architecture<128>);
  DEFINE_VCL_TYPE(Vec2q,   int64_t,  vcl_architecture<128>);
  DEFINE_VCL_TYPE(Vec2uq,  uint64_t, vcl_architecture<128>);
  DEFINE_VCL_TYPE(Vec4f,   float,    vcl_architecture<128>);
  DEFINE_VCL_TYPE(Vec2d,   double,   vcl_architecture<128>);

  DEFINE_VCL_TYPE(Vec32c,  int8_t,   vcl_architecture<256>);
  DEFINE_VCL_TYPE(Vec32uc, uint8_t,  vcl_architecture<256>);
  DEFINE_VCL_TYPE(Vec16s,  int16_t,  vcl_architecture<256>);
  DEFINE_VCL_TYPE(Vec16us, uint16_t, vcl_architecture<256>);
  DEFINE_VCL_TYPE(Vec8i,   int32_t,  vcl_architecture<256>);
  DEFINE_VCL_TYPE(Vec8ui,  uint32_t, vcl_architecture<256>);
  DEFINE_VCL_TYPE(Vec4q,   int64_t,  vcl_architecture<256>);
  DEFINE_VCL_TYPE(Vec4uq,  uint64_t, vcl_architecture<256>);
  DEFINE_VCL_TYPE(Vec8f,   float,    vcl_architecture<256>);
  DEFINE_VCL_TYPE(Vec4d,   double,   vcl_architecture<256>);

#undef DEFINE_VCL_TYPE

  template<typename Vec> void print_vec(std::ostream& s, const Vec& v)
  {
    s << v[0];
    for(unsigned i=1;i<v.size();i++)
      s << '_' << v[i];
  }

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

#undef ADD_OSTREAM_OPERATOR

namespace calin { namespace util { namespace vcl {



} } } // namespace calin::util::vcl
