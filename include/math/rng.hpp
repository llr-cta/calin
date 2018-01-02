/*

   calin/math/rng.hpp -- Stephen Fegan -- 2015-11-19

   Random number generators. Based on RandomNumbers_TNG written
   by the author while at UCLA in 2007.

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
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

/*! \file RandomNumbrs_TNG.hpp

  The next generation random number generator. Features NR3 generators,
  RanLux and old NR2 generator.

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    $Revision: 1.5 $
  \date       10/31/2007

  $Id: RandomNumbers_TNG.hpp 5429 2013-06-27 13:40:55Z sfegan $

*/

#pragma once

#include <cstdint>
#include <cmath>
#include <random>
#include <utility>

#if defined(__AVX2__)
#include <immintrin.h>
#endif

#include <calin_global_definitions.hpp>
#include <math/rng.pb.h>

namespace calin { namespace math { namespace rng {

// Decision : should this be a template or a base class. As a base
// class the code is simpler and more flexable. As a template it's
// faster. In most cases I guess the amount of time spent generating
// and processing of the deviates by the caller should outweigh the
// time spent doing the indirect call so we'll go with a base class
// and potentially reevaluate if profiling shows this assumption is
// wrong.

class RNGCore
{
public:
  virtual ~RNGCore();
  virtual uint64_t uniform_uint64() = 0;
  virtual void save_to_proto(ix::math::rng::RNGData* proto) const = 0;
  ix::math::rng::RNGData* as_proto() const {
    auto* proto = new ix::math::rng::RNGData;
    save_to_proto(proto); return proto; }
  static RNGCore* create_from_proto(const ix::math::rng::RNGData& proto,
                                    bool restore_state = false);
};

class RNG
{
public:
  enum class CoreType { NR3, RANLUX48, MT19937 /*, DEV_RANDOM*/ };
  RNG(CoreType core_type = CoreType::NR3, const std::string& created_by = "");
  RNG(const std::string& created_by, CoreType core_type = CoreType::NR3):
    RNG(core_type, created_by) { /* nothing to see here */ }
  RNG(uint64_t seed, CoreType core_type = CoreType::NR3, const std::string& created_by = "");
  RNG(uint64_t seed, const std::string& created_by, CoreType core_type = CoreType::NR3):
    RNG(seed, core_type, created_by) { /* nothing to see here */ }
  RNG(RNGCore* core, bool adopt_core = false);
  RNG(const ix::math::rng::RNGData& proto, bool restore_state = false,
    const std::string& created_by = "");

  ~RNG();

  void save_to_proto(ix::math::rng::RNGData* proto) const;

  ix::math::rng::RNGData* as_proto() const {
    auto* proto = new ix::math::rng::RNGData;
    save_to_proto(proto); return proto; }

  uint64_t uniform_uint64() { return core_->uniform_uint64(); }

  uint32_t uniform_uint32() {
    if(dev32_hascached_)
    {
      dev32_hascached_ = false;
      return dev32_cachedval_&0xFFFFFFFF;
    }
    else
    {
      dev32_hascached_ = true;
      dev32_cachedval_ = uniform_uint64();
      uint32_t ret = dev32_cachedval_ & 0xFFFFFFFF;
      dev32_cachedval_ >>= 32;
      return ret;
    }
  }

  double uniform_double() {
    return 5.42101086242752217E-20 * double(uniform_uint64());
  }

  float uniform_float() {
    return 2.328306437e-10 * float(uniform_uint32());
  }

  void uniform_by_type(uint64_t& x) { x = uniform_uint64(); }
  void uniform_by_type(uint32_t& x) { x = uniform_uint32(); }
  void uniform_by_type(double& x) { x = uniform_double(); }
  void uniform_by_type(float& x) { x = uniform_float(); }

  double uniform() { return uniform_double(); }
  double exponential() { return -std::log(uniform()); }
  double exponential(double mean) { return -mean*std::log(uniform()); }
  double normal();
  double normal(double mean, double sigma) { return mean+normal()*sigma; }
  double gamma_by_alpha_and_beta(double alpha, double beta);
  double gamma_by_mean_and_sigma(const double mean, const double sigma) {
    double b = mean/(sigma*sigma);
    return gamma_by_alpha_and_beta(mean*b,b); }
  int poisson(double lambda);
  int polya(double mean, double non_poisson_sigma) {
    return poisson(gamma_by_mean_and_sigma(mean, non_poisson_sigma)); }
  int binomial(double pp, int n);

  static uint32_t uint32_from_random_device();
  static uint64_t uint64_from_random_device();

  static uint64_t nonzero_uint64_from_random_device() { uint64_t x;
    do x = uint64_from_random_device(); while(x==0ULL); return x; }
  static uint32_t nonzero_uint32_from_random_device() { uint32_t x;
    do x = uint32_from_random_device(); while(x==0U); return x; }

  double inverse_cdf(const std::vector<std::pair<double,double>> &inv_cdf);
  static void generate_inverse_cdf(std::vector<std::pair<double,double>> &cdf,
    unsigned nbins=0);

  static constexpr uint64_t std_test_seed = 12939; // essential supply

private:
  RNGCore* core_;
  bool adopt_core_;

  // Cached values
  bool bm_hascached_ = false;
  double bm_cachedval_ = 0.0;

  bool dev32_hascached_ = false;
  uint64_t dev32_cachedval_ = 0;

  // Speedup caches
  double   poi_lambdaexp_ = 1.0/M_E;
  double   poi_lambdasrt_ = 1.0;
  double   poi_lambdalog_ = 0.0;
  double   poi_lambdaold_ = 1.0;

  int      bin_nold_      = -1;
  double   bin_pold_      = -1;
  double   bin_pc_        = 0.0;
  double   bin_plog_      = 0.0;
  double   bin_pclog_     = 0.0;
  unsigned bin_en_        = 0;
  double   bin_oldg_      = 0.0;
};

// =============================================================================
//
// Specific random number generators
//
// =============================================================================

class NR3RNGCore: public RNGCore
{
 public:
  typedef calin::ix::math::rng::NR3RNGCoreData ix_core_data_type;

  NR3RNGCore(uint64_t seed = 0):
      RNGCore(),
      seed_(seed>0 ? seed : RNG::nonzero_uint64_from_random_device()),
      u_(UINT64_C(0)), v_(UINT64_C(4101842887655102017)), w_(UINT64_C(1))
  {
    u_ = seed_^v_; uniform_uint64();
    v_ = u_; uniform_uint64();
    w_ = v_; uniform_uint64();
    calls_ = 0;
  }
  NR3RNGCore(const ix::math::rng::NR3RNGCoreData& proto,
             bool restore_state = false);
  ~NR3RNGCore();

  uint64_t uniform_uint64() override {
    calls_++;
    u_ = u_*UINT64_C(2862933555777941757) + UINT64_C(7046029254386353087);
    v_ ^= v_ >> 17;
    v_ ^= v_ << 31;
    v_ ^= v_ >> 8;
    w_ = 4294957665U*(w_ & 0xFFFFFFFF) + (w_ >> 32);
    uint64_t x = u_ ^ (u_ << 21);
    x ^= x >> 35;
    x ^= x << 4;
    return (x+v_)^w_;
  }

  void save_to_proto(ix::math::rng::NR3RNGCoreData* data) const;
  void save_to_proto(ix::math::rng::RNGData* proto) const override;

  static ix_core_data_type* mutable_core_data(ix::math::rng::RNGData* proto) {
    return proto->mutable_nr3_core(); }
  static const ix_core_data_type& core_data(const ix::math::rng::RNGData& proto)
  { return proto.nr3_core(); }

 private:
  uint64_t seed_;
  uint64_t calls_ = 0;
  uint64_t u_;
  uint64_t v_;
  uint64_t w_;
};

class Ranlux48RNGCore: public RNGCore
{
 public:
  typedef calin::ix::math::rng::Ranlux48RNGCoreData ix_core_data_type;

  Ranlux48RNGCore(uint64_t seed = 0):
      RNGCore(),
      gen_seed_(seed>0 ? seed : RNG::nonzero_uint64_from_random_device()),
      gen_(gen_seed_) { }
  Ranlux48RNGCore(const ix::math::rng::Ranlux48RNGCoreData& proto,
                  bool restore_state = false);
  ~Ranlux48RNGCore();

  uint64_t uniform_uint64() override {
    uint64_t res = gen_(); // res has 48 bits
    gen_calls_++;
    res <<= 16;
    switch(dev_blocks_)
    {
      case 0:
        dev_ = gen_();
        gen_calls_++;
        res |= dev_&0x000000000000FFFFULL;
        dev_ >>= 16;
        dev_blocks_ = 2;
        break;
      case 1:
        res |= dev_&0x000000000000FFFFULL;
        dev_blocks_ = 0;
        break;
      case 2:
        res |= dev_&0x000000000000FFFFULL;
        dev_ >>= 16;
        dev_blocks_ = 1;
        break;
      default:
        assert(0);
    }
    return res;
  }

  void save_to_proto(ix::math::rng::RNGData* proto) const override;

  static ix_core_data_type* mutable_core_data(ix::math::rng::RNGData* proto) {
    return proto->mutable_ranlux48_core(); }
  static const ix_core_data_type& core_data(const ix::math::rng::RNGData& proto)
  { return proto.ranlux48_core(); }

 private:
  uint64_t gen_seed_;
  std::ranlux48 gen_;
  uint64_t gen_calls_ = 0;
  uint_fast64_t dev_ = 0;
  uint_fast32_t dev_blocks_ = 0;
};

class MT19937RNGCore: public RNGCore
{
 public:
  typedef calin::ix::math::rng::STLRNGCoreData ix_core_data_type;

  MT19937RNGCore(uint64_t seed = 0):
      RNGCore(),
      gen_seed_(seed>0 ? seed : RNG::nonzero_uint64_from_random_device()),
      gen_(gen_seed_) { }
  MT19937RNGCore(const ix::math::rng::STLRNGCoreData& proto,
                 bool restore_state = false);
  ~MT19937RNGCore();

  uint64_t uniform_uint64() override { gen_calls_++; return gen_(); }

  void save_to_proto(ix::math::rng::RNGData* proto) const override;

  static ix_core_data_type* mutable_core_data(ix::math::rng::RNGData* proto) {
    return proto->mutable_mt19937_core(); }
  static const ix_core_data_type& core_data(const ix::math::rng::RNGData& proto)
  { return proto.mt19937_core(); }

 private:
  uint64_t gen_seed_;
  std::mt19937_64 gen_;
  uint64_t gen_calls_ = 0;
};

template<unsigned NSTREAM>
class NR3_EmulateSIMD_RNGCore: public RNGCore
{
public:
  typedef calin::ix::math::rng::NR3_SIMD_RNGCoreData ix_core_data_type;

  NR3_EmulateSIMD_RNGCore(uint64_t seed):
      RNGCore(),
      seed_(seed>0 ? seed : RNG::nonzero_uint64_from_random_device())
  {
    std::mt19937_64 gen(seed_);
    for(unsigned i=0;i<NSTREAM;i++)core_[i] = new NR3RNGCore(gen());
  }

  NR3_EmulateSIMD_RNGCore(uint64_t seeds[NSTREAM]):
      RNGCore(), seed_(0)
  {
    for(unsigned i=0;i<NSTREAM;i++)core_[i] = new NR3RNGCore(seeds[i]);
  }

  NR3_EmulateSIMD_RNGCore(const ix::math::rng::NR3_SIMD_RNGCoreData& proto,
      bool restore_state = false):
    NR3_EmulateSIMD_RNGCore(proto.seed())
  {
    if(restore_state and proto.state_saved() and proto.vec_stream_seed_size()!=0)
    {
      calls_ = proto.calls();
      if(proto.vec_stream_seed_size() != NSTREAM)
        throw std::runtime_error("NR3_EmulateSIMD_RNGCore: need " +
          std::to_string(NSTREAM) + " seeds to restore state.");
      for(unsigned i=0; i<NSTREAM; i++)
      {
        calin::ix::math::rng::NR3RNGCoreData core_data;
        core_data.set_seed(proto.vec_stream_seed(i));
        core_data.set_calls(calls_);
        core_data.set_state_saved(true);
        core_data.set_u(proto.vec_u(i));
        core_data.set_v(proto.vec_v(i));
        core_data.set_w(proto.vec_w(i));
        delete core_[NSTREAM-i-1];
        core_[NSTREAM-i-1] = new NR3RNGCore(core_data, restore_state);
      }
      ndev_ = NSTREAM;
      for(unsigned i=0; i<proto.dev_size(); i++) {
        vec_dev_[--ndev_] = proto.dev(i);
      }
    }
  }

  ~NR3_EmulateSIMD_RNGCore()
  {
    for(unsigned i=0;i<NSTREAM;i++)delete core_[i];
  }

  uint64_t uniform_uint64() override
  {
    if(ndev_ >= NSTREAM)
    {
      for(unsigned i=0;i<NSTREAM;i++)vec_dev_[i] = core_[i]->uniform_uint64();
      ndev_ = 0;
    }
    return vec_dev_[ndev_++];
  }

  void save_to_proto(ix::math::rng::RNGData* proto) const override
  {
    auto* data = proto->mutable_nr3_simd_emu4_core();
    data->set_seed(seed_);
    data->set_calls(calls_);
    data->set_state_saved(true);
    for(unsigned i=NSTREAM; i>0;)
    {
      calin::ix::math::rng::NR3RNGCoreData core_data;
      core_[--i]->save_to_proto(&core_data);
      data->add_vec_stream_seed(core_data.seed());
      data->add_vec_u(core_data.u());
      data->add_vec_v(core_data.v());
      data->add_vec_w(core_data.w());
    }
    for(unsigned i=NSTREAM; i>ndev_;)data->add_dev(vec_dev_[--i]);
  }

  static ix_core_data_type* mutable_core_data(ix::math::rng::RNGData* proto) {
    return proto->mutable_nr3_simd_emu4_core(); }
  static const ix_core_data_type& core_data(const ix::math::rng::RNGData& proto)
  { return proto.nr3_simd_emu4_core(); }

private:
  NR3RNGCore* core_[NSTREAM];
  uint64_t seed_ = 0;
  uint64_t calls_ = 0;
  uint64_t vec_dev_[NSTREAM];
  unsigned ndev_ = NSTREAM;
};

class NR3_AVX2_RNGCore: public RNGCore
{
private:
#if defined(__AVX2__)
  void test_cpu() const;

  void init(uint64_t seed0, uint64_t seed1, uint64_t seed2, uint64_t seed3)
  {
    test_cpu();

    stream_seed0_ = seed0;
    stream_seed1_ = seed1;
    stream_seed2_ = seed2;
    stream_seed3_ = seed3;

    vec_u_ = _mm256_setzero_si256();
    vec_v_ = _mm256_set1_epi64x(UINT64_C(4101842887655102017));
    vec_w_ = _mm256_set1_epi64x(UINT64_C(1));

    __m256i vec_seed = _mm256_set_epi64x(seed0,seed1,seed2,seed3);

    vec_u_ = _mm256_xor_si256(vec_v_, vec_seed);
    uniform_uivec256();
    vec_v_ = vec_u_;
    uniform_uivec256();
    vec_w_ = vec_v_;
    uniform_uivec256();
    calls_ = 0;
  }
#endif

 public:
  typedef calin::ix::math::rng::NR3_SIMD_RNGCoreData ix_core_data_type;

  NR3_AVX2_RNGCore(uint64_t seed = 0):
    RNGCore(), seed_(seed>0 ? seed : RNG::nonzero_uint64_from_random_device())
  {
#if defined(__AVX2__)
    std::mt19937_64 gen(seed_);
    init(gen(), gen(), gen(), gen());
#else
    throw std::runtime_error("NR3_AVX2_RNGCore: AVX2 not present at compile time.");
#endif
  }

  NR3_AVX2_RNGCore(uint64_t seed0, uint64_t seed1, uint64_t seed2, uint64_t seed3):
    RNGCore(), seed_(0)
  {
#if defined(__AVX2__)
    init(seed0, seed1, seed2, seed3);
#else
    throw std::runtime_error("NR3_AVX2_RNGCore: AVX2 not present at compile time.");
#endif
  }

  NR3_AVX2_RNGCore(const ix::math::rng::NR3_SIMD_RNGCoreData& proto,
                   bool restore_state = false);

  virtual ~NR3_AVX2_RNGCore();

#if defined(__AVX2__)
  __m256i uniform_uivec256()
  {
    constexpr uint64_t CU1 = UINT64_C(2862933555777941757);
    constexpr uint64_t CU2 = UINT64_C(7046029254386353087);

    // AVX2 doesn't have 64bit->64bit multiply, so instead we do
    // three 32bit->64bit multiplies, two 64bit adds and one 64bit
    // shift.  u*C = ( u_hi*C_lo + u_lo*C_hi ) << 32 + u_lo*C_lo
    const __m256i vec_cu1_lo = _mm256_set1_epi64x(CU1);
    const __m256i vec_cu1_hi = _mm256_set1_epi64x(CU1>>32);
    const __m256i vec_u_lo = vec_u_;
    const __m256i vec_u_hi = _mm256_shuffle_epi32(vec_u_, 0xB1);
    vec_u_ = _mm256_mul_epu32(vec_u_hi, vec_cu1_lo);
    vec_u_ = _mm256_add_epi64(vec_u_, _mm256_mul_epu32(vec_u_lo, vec_cu1_hi));
    vec_u_ = _mm256_slli_epi64(vec_u_, 32);
    vec_u_ = _mm256_add_epi64(vec_u_, _mm256_mul_epu32(vec_u_lo, vec_cu1_lo));
    const __m256i vec_cu2 = _mm256_set1_epi64x(CU2);
    vec_u_ =  _mm256_add_epi64(vec_u_, vec_cu2);

    vec_v_ = _mm256_xor_si256(vec_v_, _mm256_srli_epi64(vec_v_, 17));
    vec_v_ = _mm256_xor_si256(vec_v_, _mm256_slli_epi64(vec_v_, 31));
    vec_v_ = _mm256_xor_si256(vec_v_, _mm256_srli_epi64(vec_v_, 8));

    constexpr uint64_t CW = UINT64_C(4294957665);
    const __m256i vec_cw =  _mm256_set1_epi64x(CW);
    vec_w_ = _mm256_add_epi64(_mm256_mul_epu32(vec_w_, vec_cw),
                              _mm256_srli_epi64(vec_w_, 32));

    __m256i vec_x =  _mm256_xor_si256(vec_u_, _mm256_slli_epi64(vec_u_, 21));
    vec_x = _mm256_xor_si256(vec_x, _mm256_srli_epi64(vec_x, 35));
    vec_x = _mm256_xor_si256(vec_x, _mm256_slli_epi64(vec_x, 4));

    ++calls_;

    return _mm256_xor_si256(_mm256_add_epi64(vec_x, vec_v_), vec_w_);
  }

  __m256 uniform_psvec256()
  {
    __m256i vec_ui = uniform_uivec256();
    __m256 vec_ps = _mm256_cvtepi32_ps(vec_ui);
    vec_ps = _mm256_mul_ps(vec_ps, _mm256_set1_ps(2.328306437e-10));
    return _mm256_add_ps(vec_ps, _mm256_set1_ps(0.5));
  }
#endif

  uint64_t uniform_uint64() override
  {
#if defined(__AVX2__)
    if(ndev_ == 0)
    {
      vec_dev_ = uniform_uivec256();
      ndev_ = 4;
    }
    return reinterpret_cast<uint64_t*>(&vec_dev_)[--ndev_];
#else
    throw std::runtime_error("NR3_AVX2_RNGCore: AVX2 not present at compile time.");
#endif
  }

  void save_to_proto(ix::math::rng::RNGData* proto) const override;

  static ix_core_data_type* mutable_core_data(ix::math::rng::RNGData* proto) {
    return proto->mutable_nr3_avx2_core(); }
  static const ix_core_data_type& core_data(const ix::math::rng::RNGData& proto)
  { return proto.nr3_avx2_core(); }

private:
  uint64_t seed_;
#if defined(__AVX2__)
  uint64_t stream_seed0_;
  uint64_t stream_seed1_;
  uint64_t stream_seed2_;
  uint64_t stream_seed3_;
  uint64_t calls_ = 0;
  __m256i vec_u_;
  __m256i vec_v_;
  __m256i vec_w_;
  __m256i vec_dev_;
  unsigned ndev_ = 0;
#endif
};

} } } // namespace calin::math::rng
