/* 

   calin/math/rng.cpp -- Stephen Fegan -- 2015-11-19

   Random number generators. Based on RandomNumbers_TNG written
   by the author while at UCLA in 2007.

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, Ecole polytechnique, CNRS/IN2P3, Universite Paris-Saclay

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

#include <proto/math_rng.pb.h>

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
  virtual double uniform_double() = 0;
  virtual uint64_t uniform_uint64() = 0;
  virtual uint32_t uniform_uint32() = 0;
  virtual void save_to_proto(ix::math::rng::RNGData* proto) const = 0;
  static RNGCore* create_from_proto(const ix::math::rng::RNGData& proto,
                                    bool restore_state = false);
  void uniform_by_type(uint64_t& x) { x = uniform_uint64(); }
  void uniform_by_type(uint32_t& x) { x = uniform_uint32(); }
  void uniform_by_type(double& x) { x = uniform_double(); }
};

class RNG
{
 public:
  enum class CoreType { NR3, RANLUX48, MT19937 /*, DEV_RANDOM*/ };
  RNG(CoreType core_type = CoreType::NR3);
  RNG(uint64_t seed, CoreType core_type = CoreType::NR3);
  RNG(RNGCore* core, bool adopt_core = false);
  RNG(const ix::math::rng::RNGData& proto, bool restore_state = false);

  ~RNG();

  void save_to_proto(ix::math::rng::RNGData* proto);

  double uniform_double() { return core_->uniform_double(); }
  double uniform_uint64() { return core_->uniform_uint64(); }
  double uniform_uint32() { return core_->uniform_uint32(); }
  
  double uniform() { return core_->uniform_double(); }
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
  
#if 0
  inline double inverse_cdf(const std::vector< Pair > &inv_cdf);
  static void generate_inverse_cdf(std::vector< Pair > &cdf, unsigned nbins=0);
#endif
  
 private:
  RNGCore* core_;
  bool adopt_core_;
  
  bool bm_hascached_ = false;
  double bm_cachedval_ = 0.0;

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
      RNGCore(), seed_(seed>0 ? seed : RNG::nonzero_uint64_from_random_device()),
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

  uint32_t uniform_uint32() override {
    return uint32_t(uniform_uint64()&0xFFFFFFFF);
  }

  double uniform_double() override {
    return 5.42101086242752217E-20 * double(uniform_uint64());
  }
  
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

  uint32_t uniform_uint32() override {
    uint32_t res = 0;
    switch(dev_blocks_)
    {
      case 0:
        dev_ = gen_();
        gen_calls_++;
        res = dev_&0x00000000FFFFFFFFULL;
        dev_ >>= 32;
        dev_blocks_ = 1;
        break;
      case 1:
        res = dev_&0x000000000000FFFFULL;
        res <<= 16;
        dev_ = gen_();
        gen_calls_++;
        res |= dev_&0x000000000000FFFFULL;
        dev_ >>= 16;
        dev_blocks_ = 2;
        break;
      case 2:
        res = dev_&0x00000000FFFFFFFF;
        dev_blocks_ = 0;
        break;
      default:
        assert(0);
    }
    return res;
  }

  double uniform_double() override {
    return 5.421001086242752217E-20 * double(uniform_uint64());
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
  uint32_t uniform_uint32() override { gen_calls_++;
    return uint32_t(gen_()&0xFFFFFFFF); }
  double uniform_double() override {
    return 5.421001086242752217E-20 * double(uniform_uint64());
  }
  
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
  
} } } // namespace calin::math::rng
