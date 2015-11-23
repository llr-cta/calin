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
  virtual void save_to_proto() const = 0;
  static RNGCore* create_from_proto();
};

class RNG
{
 public:
  RNG(uint64_t seed = 0);
  RNG(RNGCore* core, bool adopt_core = false);
  ~RNG() { if(adopt_core_)delete core_; }
  
  double uniform_double() { return core_->uniform_double(); }
  double uniform_uint64() { return core_->uniform_uint64(); }
  double uniform_uint32() { return core_->uniform_uint32(); }
  
  double uniform() { return core_->uniform_double(); }
  double exponential() { return -std::log(uniform()); }
  double exponential(double mean) { return -mean*std::log(uniform()); }
  double normal();
  double normal(double mean, double sigma) { return mean+normal()*sigma; }
  double gamma(double alpha, double beta=1.0);
  double gamma_by_mean_and_sigma(const double mean, const double sigma) {
    double b = mean/sigma/sigma;
    return gamma(mean*b,b); }
  int poisson(double lambda);
  int polya(double mean, double non_poisson_sigma) {
    return poisson(gamma_by_mean_and_sigma(mean, non_poisson_sigma)); }
  int binomial(double pp, int n);

  static uint32_t uint32_from_random_device();
  static uint64_t uint64_from_random_device();
  
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
  NR3RNGCore(uint64_t seed):
      m_u(UINT64_C(0)), m_v(UINT64_C(4101842887655102017)), m_w(UINT64_C(1))
  {
    m_u = seed^m_v; uniform_uint64();
    m_v = m_u; uniform_uint64();
    m_w = m_v; uniform_uint64();
  }

  uint64_t uniform_uint64() override {
    m_u = m_u*UINT64_C(2862933555777941757) + UINT64_C(7046029254386353087);
    m_v ^= m_v >> 17; 
    m_v ^= m_v << 31;
    m_v ^= m_v >> 8;
    m_w = 4294957665U*(m_w & 0xFFFFFFFF) + (m_w >> 32);
    uint64_t x = m_u ^ (m_u << 21);
    x ^= x >> 35;
    x ^= x << 4;
    return (x+m_v)^m_w;
  }

  uint32_t uniform_uint32() override {
    return uint32_t(uniform_uint64()&0xFFFFFFFF);
  }

  double uniform_double() override {
    return 5.42101086242752217E-20 * double(uniform_uint64());
  }
  
  void save_to_proto() const override;
  static RNGCore* create_from_proto();
  
 private:
  uint64_t m_u;
  uint64_t m_v;
  uint64_t m_w;
};

} } } // namespace calin::math::rng
