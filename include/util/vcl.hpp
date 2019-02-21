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
#include <string>

#define MAX_VECTOR_SIZE 512
#define VCL_NAMESPACE vcl

#include <VCL/vectorclass.h>
#include <VCL/vectormath_trig.h>
#include <VCL/vectormath_exp.h>
#include <VCL/vectormath_hyp.h>

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Dense>

namespace Eigen {

template<> struct NumTraits<::vcl::Vec4f>: NumTraits<float>
{
  typedef ::vcl::Vec4f Real;
  typedef ::vcl::Vec4f NonInteger;
  typedef ::vcl::Vec4f Nested;
  typedef float Literal;
  enum {
    IsComplex = 0,
    IsInteger = 0,
    IsSigned = 1,
    RequireInitialization = 0,
    ReadCost = 1,
    AddCost = 3,
    MulCost = 3
  };
};

template<> struct NumTraits<::vcl::Vec8f>: NumTraits<float>
{
  typedef ::vcl::Vec8f Real;
  typedef ::vcl::Vec8f NonInteger;
  typedef ::vcl::Vec8f Nested;
  typedef float Literal;
  enum {
    IsComplex = 0,
    IsInteger = 0,
    IsSigned = 1,
    RequireInitialization = 0,
    ReadCost = 1,
    AddCost = 3,
    MulCost = 3
  };
};

template<> struct NumTraits<::vcl::Vec16f>: NumTraits<float>
{
  typedef ::vcl::Vec16f Real;
  typedef ::vcl::Vec16f NonInteger;
  typedef ::vcl::Vec16f Nested;
  typedef float Literal;
  enum {
    IsComplex = 0,
    IsInteger = 0,
    IsSigned = 1,
    RequireInitialization = 0,
    ReadCost = 1,
    AddCost = 3,
    MulCost = 3
  };
};

template<> struct NumTraits<::vcl::Vec2d>: NumTraits<double>
{
  typedef ::vcl::Vec2d Real;
  typedef ::vcl::Vec2d NonInteger;
  typedef ::vcl::Vec2d Nested;
  typedef double Literal;
  enum {
    IsComplex = 0,
    IsInteger = 0,
    IsSigned = 1,
    RequireInitialization = 0,
    ReadCost = 1,
    AddCost = 3,
    MulCost = 3
  };
};

template<> struct NumTraits<::vcl::Vec4d>: NumTraits<double>
{
  typedef ::vcl::Vec4d Real;
  typedef ::vcl::Vec4d NonInteger;
  typedef ::vcl::Vec4d Nested;
  typedef double Literal;
  enum {
    IsComplex = 0,
    IsInteger = 0,
    IsSigned = 1,
    RequireInitialization = 0,
    ReadCost = 1,
    AddCost = 3,
    MulCost = 3
  };
};

template<> struct NumTraits<::vcl::Vec8d>: NumTraits<double>
{
  typedef ::vcl::Vec8d Real;
  typedef ::vcl::Vec8d NonInteger;
  typedef ::vcl::Vec8d Nested;
  typedef double Literal;
  enum {
    IsComplex = 0,
    IsInteger = 0,
    IsSigned = 1,
    RequireInitialization = 0,
    ReadCost = 1,
    AddCost = 3,
    MulCost = 3
  };
};

typedef Eigen::Matrix< ::vcl::Vec4f , 3 , 1> Vector3_4f;
typedef Eigen::Matrix< ::vcl::Vec4f , 3 , 3> Matrix3_4f;
typedef Eigen::Matrix< ::vcl::Vec8f , 3 , 1> Vector3_8f;
typedef Eigen::Matrix< ::vcl::Vec8f , 3 , 3> Matrix3_8f;
typedef Eigen::Matrix< ::vcl::Vec16f , 3 , 1> Vector3_16f;
typedef Eigen::Matrix< ::vcl::Vec16f , 3 , 3> Matrix3_16f;

typedef Eigen::Matrix< ::vcl::Vec2d , 3 , 1> Vector3_2d;
typedef Eigen::Matrix< ::vcl::Vec2d , 3 , 3> Matrix3_2d;
typedef Eigen::Matrix< ::vcl::Vec4d , 3 , 1> Vector3_4d;
typedef Eigen::Matrix< ::vcl::Vec4d , 3 , 3> Matrix3_4d;
typedef Eigen::Matrix< ::vcl::Vec8d , 3 , 1> Vector3_8d;
typedef Eigen::Matrix< ::vcl::Vec8d , 3 , 3> Matrix3_8d;

} // namespace Eigen

namespace vcl {

  static inline Vec4d to_double_low(Vec8i const & a) {
    return to_double(a.get_low());
  }

  static inline Vec4d to_double_high(Vec8i const & a) {
    return to_double(a.get_high());
  }

  static inline Vec8d to_double_low(Vec16i const & a) {
    return to_double(a.get_low());
  }

  static inline Vec8d to_double_high(Vec16i const & a) {
    return to_double(a.get_high());
  }

  static inline Vec8i truncate_to_int(Vec4d const & a, Vec4d const & b) {
    Vec4i t1 = truncate_to_int(a);
    Vec4i t2 = truncate_to_int(b);
    return Vec8i(t1,t2);
  }

  static inline Vec16i truncate_to_int(Vec8d const & a, Vec8d const & b) {
    Vec8i t1 = truncate_to_int(a);
    Vec8i t2 = truncate_to_int(b);
    return Vec16i(t1,t2);
  }

}

namespace calin { namespace util { namespace vcl {

using namespace ::vcl;

template<typename VCLArchitecture> struct VCLFloatReal
{
  constexpr static unsigned num_real              = VCLArchitecture::num_float;
  typedef VCLArchitecture                         architecture;

  typedef int32_t                                 int_t;
  typedef uint32_t                                uint_t;
  typedef float                                   real_t;
  typedef Eigen::Vector3f                         vec3_t;
  typedef Eigen::Matrix3f                         mat3_t;

  typedef typename VCLArchitecture::int32_at      int_at;
  typedef typename VCLArchitecture::uint32_at     uint_at;
  typedef typename VCLArchitecture::float_at      real_at;

  typedef typename VCLArchitecture::int32_vt      int_vt;
  typedef typename VCLArchitecture::uint32_vt     uint_vt;
  typedef typename VCLArchitecture::int32_bvt     bool_int_vt;
  typedef typename VCLArchitecture::uint32_bvt    bool_uint_vt;

  typedef typename VCLArchitecture::float_vt      real_vt;
  typedef typename VCLArchitecture::float_bvt     bool_vt;

  typedef typename VCLArchitecture::Vector3f_vt   vec3_vt;
  typedef typename VCLArchitecture::Matrix3f_vt   mat3_vt;

  static inline int_vt truncate_to_int_limited(real_vt x) {
    return vcl::truncate_to_int(x);
  }

  static inline int_vt round_to_int_limited(real_vt x) {
    return vcl::round_to_int(x);
  }

  static inline int_vt truncate_to_int(real_vt x) {
    return vcl::truncate_to_int(x);
  }

  static inline int_vt round_to_int(real_vt x) {
    return vcl::round_to_int(x);
  }
};

template<typename VCLArchitecture> struct VCLDoubleReal
{
  constexpr static unsigned num_real              = VCLArchitecture::num_double;
  typedef VCLArchitecture                         architecture;

  typedef int64_t                                 int_t;
  typedef uint64_t                                uint_t;
  typedef double                                  real_t;
  typedef Eigen::Vector3d                         vec3_t;
  typedef Eigen::Matrix3d                         mat3_t;

  typedef typename VCLArchitecture::int64_at      int_at;
  typedef typename VCLArchitecture::uint64_at     uint_at;
  typedef typename VCLArchitecture::double_at     real_at;

  typedef typename VCLArchitecture::int64_vt      int_vt;
  typedef typename VCLArchitecture::uint64_vt     uint_vt;
  typedef typename VCLArchitecture::int64_bvt     bool_int_vt;
  typedef typename VCLArchitecture::uint64_bvt    bool_uint_vt;

  typedef typename VCLArchitecture::double_vt     real_vt;
  typedef typename VCLArchitecture::double_bvt    bool_vt;

  typedef typename VCLArchitecture::Vector3d_vt   vec3_vt;
  typedef typename VCLArchitecture::Matrix3d_vt   mat3_vt;

  static inline int_vt truncate_to_int_limited(real_vt x) {
    return vcl::truncate_to_int64_limited(x);
  }

  static inline int_vt round_to_int_limited(real_vt x) {
    return vcl::round_to_int64_limited(x);
  }

  static inline int_vt truncate_to_int(real_vt x) {
    return vcl::truncate_to_int64(x);
  }

  static inline int_vt round_to_int(real_vt x) {
    return vcl::round_to_int64(x);
  }
};

struct VCL128Architecture
{
  constexpr static char architecture_name[] = "VCL128Architecture";

  constexpr static unsigned vec_bits   = 128;
  constexpr static unsigned vec_bytes  = vec_bits/8;
  constexpr static unsigned num_int8   = vec_bytes/sizeof(int8_t);
  constexpr static unsigned num_uint8  = vec_bytes/sizeof(uint8_t);
  constexpr static unsigned num_int16  = vec_bytes/sizeof(int16_t);
  constexpr static unsigned num_uint16 = vec_bytes/sizeof(uint16_t);
  constexpr static unsigned num_int32  = vec_bytes/sizeof(int32_t);
  constexpr static unsigned num_uint32 = vec_bytes/sizeof(uint32_t);
  constexpr static unsigned num_int64  = vec_bytes/sizeof(int64_t);
  constexpr static unsigned num_uint64 = vec_bytes/sizeof(uint64_t);
  constexpr static unsigned num_float  = vec_bytes/sizeof(float);
  constexpr static unsigned num_double = vec_bytes/sizeof(double);

  typedef Vec128b bool_vt;
  typedef Vec16c  int8_vt;
  typedef Vec16uc uint8_vt;
  typedef Vec8s   int16_vt;
  typedef Vec8us  uint16_vt;
  typedef Vec4i   int32_vt;
  typedef Vec4ui  uint32_vt;
  typedef Vec2q   int64_vt;
  typedef Vec2uq  uint64_vt;
  typedef Vec4f   float_vt;
  typedef Vec2d   double_vt;

  typedef Vec16cb int8_bvt;
  typedef Vec16cb uint8_bvt;
  typedef Vec8sb  int16_bvt;
  typedef Vec8sb  uint16_bvt;
  typedef Vec4ib  int32_bvt;
  typedef Vec4ib  uint32_bvt;
  typedef Vec2qb  int64_bvt;
  typedef Vec2qb  uint64_bvt;
  typedef Vec4fb  float_bvt;
  typedef Vec2db  double_bvt;

  typedef int32_t  int32_at[num_int32] __attribute((aligned(vec_bytes)));
  typedef uint32_t uint32_at[num_int32] __attribute((aligned(vec_bytes)));
  typedef int64_t  int64_at[num_int64] __attribute((aligned(vec_bytes)));
  typedef uint64_t uint64_at[num_int64] __attribute((aligned(vec_bytes)));
  typedef float    float_at[num_float] __attribute((aligned(vec_bytes)));
  typedef double   double_at[num_double] __attribute((aligned(vec_bytes)));

  typedef Eigen::Vector3_4f Vector3f_vt;
  typedef Eigen::Matrix3_4f Matrix3f_vt;

  typedef Eigen::Vector3_2d  Vector3d_vt;
  typedef Eigen::Matrix3_2d  Matrix3d_vt;

  typedef VCLFloatReal<VCL128Architecture> float_real;
  typedef VCLDoubleReal<VCL128Architecture> double_real;
};

struct VCL256Architecture
{
  constexpr static char architecture_name[] = "VCL256Architecture";

  constexpr static unsigned vec_bits   = 256;
  constexpr static unsigned vec_bytes  = vec_bits/8;
  constexpr static unsigned num_int8   = vec_bytes/sizeof(int8_t);
  constexpr static unsigned num_uint8  = vec_bytes/sizeof(uint8_t);
  constexpr static unsigned num_int16  = vec_bytes/sizeof(int16_t);
  constexpr static unsigned num_uint16 = vec_bytes/sizeof(uint16_t);
  constexpr static unsigned num_int32  = vec_bytes/sizeof(int32_t);
  constexpr static unsigned num_uint32 = vec_bytes/sizeof(uint32_t);
  constexpr static unsigned num_int64  = vec_bytes/sizeof(int64_t);
  constexpr static unsigned num_uint64 = vec_bytes/sizeof(uint64_t);
  constexpr static unsigned num_float  = vec_bytes/sizeof(float);
  constexpr static unsigned num_double = vec_bytes/sizeof(double);

  typedef Vec256b bool_vt;
  typedef Vec32c  int8_vt;
  typedef Vec32uc uint8_vt;
  typedef Vec16s  int16_vt;
  typedef Vec16us uint16_vt;
  typedef Vec8i   int32_vt;
  typedef Vec8ui  uint32_vt;
  typedef Vec4q   int64_vt;
  typedef Vec4uq  uint64_vt;
  typedef Vec8f   float_vt;
  typedef Vec4d   double_vt;

  typedef Vec32cb int8_bvt;
  typedef Vec32cb uint8_bvt;
  typedef Vec16sb int16_bvt;
  typedef Vec16sb uint16_bvt;
  typedef Vec8ib  int32_bvt;
  typedef Vec8ib  uint32_bvt;
  typedef Vec4qb  int64_bvt;
  typedef Vec4qb  uint64_bvt;
  typedef Vec8fb  float_bvt;
  typedef Vec4db  double_bvt;

  typedef int32_t  int32_at[num_int32] __attribute((aligned(vec_bytes)));
  typedef uint32_t uint32_at[num_int32] __attribute((aligned(vec_bytes)));
  typedef int64_t  int64_at[num_int64] __attribute((aligned(vec_bytes)));
  typedef uint64_t uint64_at[num_int64] __attribute((aligned(vec_bytes)));
  typedef float    float_at[num_float] __attribute((aligned(vec_bytes)));
  typedef double   double_at[num_double] __attribute((aligned(vec_bytes)));

  typedef Eigen::Vector3_8f Vector3f_vt;
  typedef Eigen::Matrix3_8f Matrix3f_vt;

  typedef Eigen::Vector3_4d  Vector3d_vt;
  typedef Eigen::Matrix3_4d  Matrix3d_vt;

  typedef VCLFloatReal<VCL256Architecture> float_real;
  typedef VCLDoubleReal<VCL256Architecture> double_real;
};

struct VCL512Architecture
{
  constexpr static char architecture_name[] = "VCL512Architecture";

  constexpr static unsigned vec_bits   = 512;
  constexpr static unsigned vec_bytes  = vec_bits/8;
  constexpr static unsigned num_int8   = vec_bytes/sizeof(int8_t);
  constexpr static unsigned num_uint8  = vec_bytes/sizeof(uint8_t);
  constexpr static unsigned num_int16  = vec_bytes/sizeof(int16_t);
  constexpr static unsigned num_uint16 = vec_bytes/sizeof(uint16_t);
  constexpr static unsigned num_int32  = vec_bytes/sizeof(int32_t);
  constexpr static unsigned num_uint32 = vec_bytes/sizeof(uint32_t);
  constexpr static unsigned num_int64  = vec_bytes/sizeof(int64_t);
  constexpr static unsigned num_uint64 = vec_bytes/sizeof(uint64_t);
  constexpr static unsigned num_float  = vec_bytes/sizeof(float);
  constexpr static unsigned num_double = vec_bytes/sizeof(double);

  typedef Vec16i  int32_vt;
  typedef Vec16ui uint32_vt;
  typedef Vec8q   int64_vt;
  typedef Vec8uq  uint64_vt;
  typedef Vec16f  float_vt;
  typedef Vec8d   double_vt;

  typedef Vec16ib int32_bvt;
  typedef Vec16ib uint32_bvt;
  typedef Vec8qb  int64_bvt;
  typedef Vec8qb  uint64_bvt;
  typedef Vec16fb float_bvt;
  typedef Vec8db  double_bvt;

  typedef int32_t  int32_at[num_int32] __attribute((aligned(vec_bytes)));
  typedef uint32_t uint32_at[num_int32] __attribute((aligned(vec_bytes)));
  typedef int64_t  int64_at[num_int64] __attribute((aligned(vec_bytes)));
  typedef uint64_t uint64_at[num_int64] __attribute((aligned(vec_bytes)));
  typedef float    float_at[num_float] __attribute((aligned(vec_bytes)));
  typedef double   double_at[num_double] __attribute((aligned(vec_bytes)));

  typedef Eigen::Vector3_16f Vector3f_vt;
  typedef Eigen::Matrix3_16f Matrix3f_vt;

  typedef Eigen::Vector3_8d  Vector3d_vt;
  typedef Eigen::Matrix3_8d  Matrix3d_vt;

  typedef VCLFloatReal<VCL512Architecture> float_real;
  typedef VCLDoubleReal<VCL512Architecture> double_real;
};

typedef VCLFloatReal<VCL128Architecture> VCL128FloatReal;
typedef VCLFloatReal<VCL256Architecture> VCL256FloatReal;
typedef VCLFloatReal<VCL512Architecture> VCL512FloatReal;

typedef VCLDoubleReal<VCL128Architecture> VCL128DoubleReal;
typedef VCLDoubleReal<VCL256Architecture> VCL256DoubleReal;
typedef VCLDoubleReal<VCL512Architecture> VCL512DoubleReal;

template<typename VCLArchitecture> std::string templated_class_name(
  const std::string& class_name)
{
  return class_name + "<" + VCLArchitecture::architecture_name + ">";
}

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
  for(int i=1;i<v.size();i++)
    s << '_' << v[i];
}

template<typename VCLReal> inline void
insert_into_vec3_with_mask(typename VCLReal::vec3_vt& vv,
    const typename VCLReal::vec3_t& vs, const typename VCLReal::bool_vt& mask) {
  vv.x() = vcl::select(mask, vs.x(), vv.x());
  vv.y() = vcl::select(mask, vs.y(), vv.y());
  vv.z() = vcl::select(mask, vs.z(), vv.z());
}

template<typename VCLReal> inline void
insert_into_with_mask(typename VCLReal::real_vt& v,
    const typename VCLReal::real_t& s, const typename VCLReal::bool_vt& mask) {
  v = vcl::select(mask, s, v);
}

typedef Eigen::Matrix< ::vcl::Vec4f , 3 , 1> Vector3_4f;
typedef Eigen::Matrix< ::vcl::Vec4f , 3 , 3> Matrix3_4f;
typedef Eigen::Matrix< ::vcl::Vec8f , 3 , 1> Vector3_8f;
typedef Eigen::Matrix< ::vcl::Vec8f , 3 , 3> Matrix3_8f;
typedef Eigen::Matrix< ::vcl::Vec16f , 3 , 1> Vector3_16f;
typedef Eigen::Matrix< ::vcl::Vec16f , 3 , 3> Matrix3_16f;

typedef Eigen::Matrix< ::vcl::Vec2d , 3 , 1> Vector3_2d;
typedef Eigen::Matrix< ::vcl::Vec2d , 3 , 3> Matrix3_2d;
typedef Eigen::Matrix< ::vcl::Vec4d , 3 , 1> Vector3_4d;
typedef Eigen::Matrix< ::vcl::Vec4d , 3 , 3> Matrix3_4d;
typedef Eigen::Matrix< ::vcl::Vec8d , 3 , 1> Vector3_8d;
typedef Eigen::Matrix< ::vcl::Vec8d , 3 , 3> Matrix3_8d;


inline Vec2uq mul_low32_packed64(const Vec2uq& a, const Vec2uq& b) {
  return _mm_mul_epu32(a, b);
}

#if (defined (__AVX512DQ__) && defined (__AVX512VL__)) || INSTRSET < 5
  inline Vec2uq mul_64(const Vec2uq& a, const Vec2uq& b) {
    return a*b;
  }
#else
  inline Vec2uq mul_64(const Vec2uq& a, const Vec2uq& b) {
    __m128i prod = _mm_mul_epu32(_mm_shuffle_epi32(a, 0xB1), b);
    __m128i tmp = _mm_mul_epu32(a, _mm_shuffle_epi32(b, 0xB1));
    prod = _mm_add_epi64(prod, tmp);
    prod = _mm_slli_epi64(prod, 32);
    tmp = _mm_mul_epu32(a, b);
    prod = _mm_add_epi64(prod, tmp);
    return prod;
  }
#endif

#if MAX_VECTOR_SIZE >= 256
#if INSTRSET >= 8
  inline Vec4uq mul_low32_packed64(const Vec4uq& a, const Vec4uq& b) {
    return _mm256_mul_epu32(a, b);
  }

#if defined (__AVX512DQ__) && defined (__AVX512VL__)
  inline Vec4uq mul_64(const Vec4uq& a, const Vec4uq& b) {
    return a*b;
  }
#else
  inline Vec4uq mul_64(const Vec4uq& a, const Vec4uq& b) {
    __m256i prod = _mm256_mul_epu32(_mm256_shuffle_epi32(a, 0xB1), b);
    __m256i tmp = _mm256_mul_epu32(a, _mm256_shuffle_epi32(b, 0xB1));
    prod = _mm256_add_epi64(prod, tmp);
    prod = _mm256_slli_epi64(prod, 32);
    tmp = _mm256_mul_epu32(a, b);
    prod = _mm256_add_epi64(prod, tmp);
    return prod;
  }
#endif // defined (__AVX512DQ__) && defined (__AVX512VL__)
#else // INSTRSET < 8
  inline Vec4uq mul_low32_packed64(const Vec4uq& a, const Vec4uq& b) {
    return Vec4uq(mul_low32_packed64(a.get_low(), b.get_low()),
                  mul_low32_packed64(a.get_high(), b.get_high()));
  }

  inline Vec4uq mul_64(const Vec4uq& a, const Vec4uq& b) {
    return a*b;
  }
#endif // INSTRSET < 8
#endif // MAX_VECTOR_SIZE >= 256

#if MAX_VECTOR_SIZE >= 512
#if INSTRSET >= 9
  inline Vec8uq mul_low32_packed64(const Vec8uq& a, const Vec8uq& b) {
    return _mm512_mul_epu32(a, b);
  }
#else // INSTRSET < 9
  inline Vec8uq mul_low32_packed64(const Vec8uq& a, const Vec8uq& b) {
    return Vec8uq(mul_low32_packed64(a.get_low(), b.get_low()),
                  mul_low32_packed64(a.get_high(), b.get_high()));
  }
#endif // INSTRSET < 9

  inline Vec8uq mul_64(const Vec8uq& a, const Vec8uq& b) {
    return a*b;
  }

#endif // MAX_VECTOR_SIZE >= 512

} } } // namespace calin::util::vcl

#define ADD_OSTREAM_OPERATOR(Vec) \
  inline std::ostream& operator<<(std::ostream& s, const Vec& v) { \
    calin::util::vcl::print_vec(s, v); \
    return s; \
  }

namespace vcl {

ADD_OSTREAM_OPERATOR(Vec16c);
ADD_OSTREAM_OPERATOR(Vec16uc);
ADD_OSTREAM_OPERATOR(Vec8s);
ADD_OSTREAM_OPERATOR(Vec8us);
ADD_OSTREAM_OPERATOR(Vec4i);
ADD_OSTREAM_OPERATOR(Vec4ui);
ADD_OSTREAM_OPERATOR(Vec2q);
ADD_OSTREAM_OPERATOR(Vec2uq);
ADD_OSTREAM_OPERATOR(Vec4f);
ADD_OSTREAM_OPERATOR(Vec2d);

ADD_OSTREAM_OPERATOR(Vec16cb);
ADD_OSTREAM_OPERATOR(Vec8sb);
ADD_OSTREAM_OPERATOR(Vec4ib);
ADD_OSTREAM_OPERATOR(Vec2qb);
ADD_OSTREAM_OPERATOR(Vec4fb);
ADD_OSTREAM_OPERATOR(Vec2db);

ADD_OSTREAM_OPERATOR(Vec32c);
ADD_OSTREAM_OPERATOR(Vec32uc);
ADD_OSTREAM_OPERATOR(Vec16s);
ADD_OSTREAM_OPERATOR(Vec16us);
ADD_OSTREAM_OPERATOR(Vec8i);
ADD_OSTREAM_OPERATOR(Vec8ui);
ADD_OSTREAM_OPERATOR(Vec4q);
ADD_OSTREAM_OPERATOR(Vec4uq);
ADD_OSTREAM_OPERATOR(Vec8f);
ADD_OSTREAM_OPERATOR(Vec4d);

ADD_OSTREAM_OPERATOR(Vec32cb);
ADD_OSTREAM_OPERATOR(Vec16sb);
ADD_OSTREAM_OPERATOR(Vec8ib);
ADD_OSTREAM_OPERATOR(Vec4qb);
ADD_OSTREAM_OPERATOR(Vec8fb);
ADD_OSTREAM_OPERATOR(Vec4db);

// ADD_OSTREAM_OPERATOR(Vec64c);
// ADD_OSTREAM_OPERATOR(Vec64uc);
// ADD_OSTREAM_OPERATOR(Vec32s);
// ADD_OSTREAM_OPERATOR(Vec32us);
ADD_OSTREAM_OPERATOR(Vec16i);
ADD_OSTREAM_OPERATOR(Vec16ui);
ADD_OSTREAM_OPERATOR(Vec8q);
ADD_OSTREAM_OPERATOR(Vec8uq);
ADD_OSTREAM_OPERATOR(Vec16f);
ADD_OSTREAM_OPERATOR(Vec8d);

ADD_OSTREAM_OPERATOR(Vec16ib);
ADD_OSTREAM_OPERATOR(Vec8qb);
ADD_OSTREAM_OPERATOR(Vec16fb);
ADD_OSTREAM_OPERATOR(Vec8db);
}

#undef ADD_OSTREAM_OPERATOR

namespace calin { namespace util { namespace vcl {

//template

} } } // namespace calin::util::vcl
