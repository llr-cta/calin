/*

   calin/math/rng.hpp -- Stephen Fegan -- 2015-11-19

   Random number generators. Based on RandomNumbers_TNG written
   by the author while at UCLA in 2007.

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
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
#include <cstddef>

#include <calin_global_definitions.hpp>
#include <math/rng.pb.h>
#include <provenance/chronicle.pb.h>

namespace calin { namespace math { namespace rng {

// Decision : should this be a template or a base class. As a base
// class the code is simpler and more flexable. As a template it's
// faster. In most cases I guess the amount of time spent generating
// and processing of the deviates by the caller should outweigh the
// time spent doing the indirect call so we'll go with a base class
// and potentially reevaluate if profiling shows this assumption is
// wrong.

constexpr float C_U32_TO_FLT = 2.328306437E-10;
constexpr double C_U64_TO_DBL = 5.42101086242752217E-20;

class RNGCore
{
public:
  virtual ~RNGCore();
  virtual uint64_t uniform_uint64() = 0;
  virtual void bulk_uniform_uint64(void* buffer, std::size_t nbytes);
  virtual void bulk_uniform_uint64_with_mask(void* buffer, std::size_t nbytes,
    uint64_t mask = 0xFFFFFFFFFFFFFFFFU);
  virtual void save_to_proto(ix::math::rng::RNGCoreData* proto) const = 0;
  ix::math::rng::RNGCoreData* as_proto() const {
    auto* proto = new ix::math::rng::RNGCoreData;
    save_to_proto(proto); return proto; }
  static RNGCore* create_from_proto(const ix::math::rng::RNGCoreData& proto,
    bool restore_state = false,
    const std::string& created_by = "", const std::string& comment = "");

  std::vector<uint64_t> vec_uniform_uint64(std::size_t nelements);
  std::vector<uint64_t> vec_uniform_uint64_with_mask(std::size_t nelements,
    uint64_t mask = 0xFFFFFFFFFFFFFFFFU);

protected:
  void write_provenance(const std::string& created_by, const std::string& comment = "");
  calin::ix::provenance::chronicle::RNGRecord* chronicle_record_ = nullptr;
};

class RNG
{
public:
  enum class CoreType { NR3, RANLUX48, MT19937 /*, DEV_RANDOM*/ };
  RNG(CoreType core_type = CoreType::NR3, const std::string& created_by = "",
    const std::string& comment = "");
  RNG(const std::string& created_by, const std::string& comment = "",
      CoreType core_type = CoreType::NR3):
    RNG(core_type, created_by, comment) { /* nothing to see here */ }
  RNG(uint64_t seed, CoreType core_type = CoreType::NR3, const std::string& created_by = "",
    const std::string& comment = "");
  RNG(uint64_t seed, const std::string& created_by, const std::string& comment = "",
      CoreType core_type = CoreType::NR3):
    RNG(seed, core_type, created_by, comment) { /* nothing to see here */ }
  RNG(RNGCore* core, bool adopt_core = false,
    const std::string& created_by = "", const std::string& comment = "");
  RNG(const ix::math::rng::RNGData& proto, bool restore_state = false,
    const std::string& created_by = "", const std::string& comment = "");

  ~RNG();

  void save_to_proto(ix::math::rng::RNGData* proto) const;

  ix::math::rng::RNGData* as_proto() const {
    auto* proto = new ix::math::rng::RNGData;
    save_to_proto(proto); return proto; }

  uint64_t uniform_uint64() { return core_->uniform_uint64(); }

  uint32_t uniform_uint32() { return uint32_t(uniform_uint64()); }

  double uniform_double() {
    return C_U64_TO_DBL * double(uniform_uint64());
  }

  float uniform_float() {
    return C_U32_TO_FLT * float(uniform_uint32());
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
  calin::ix::provenance::chronicle::RNGRecord* chronicle_record_ = nullptr;

  // Cached values
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

constexpr uint64_t C_NR3_U_INIT = UINT64_C(0);
constexpr uint64_t C_NR3_V_INIT = UINT64_C(4101842887655102017);
constexpr uint64_t C_NR3_W_INIT = UINT64_C(1);

constexpr uint64_t C_NR3_U_MUL = UINT64_C(2862933555777941757);
constexpr uint64_t C_NR3_U_ADD = UINT64_C(7046029254386353087);

constexpr int C_NR3_V_SHIFT1 = 17;
constexpr int C_NR3_V_SHIFT2 = 31;
constexpr int C_NR3_V_SHIFT3 = 8;

constexpr uint64_t C_NR3_W_MUL = UINT64_C(4294957665);
constexpr int C_NR3_W_SHIFT1 = 32;

constexpr int C_NR3_X_SHIFT1 = 21;
constexpr int C_NR3_X_SHIFT2 = 35;
constexpr int C_NR3_X_SHIFT3 = 4;

class NR3RNGCore: public RNGCore
{
 public:
  typedef calin::ix::math::rng::NR3RNGCoreData ix_core_data_type;

  NR3RNGCore(uint64_t seed = 0,
    const std::string& created_by = "", const std::string& comment = "");
  NR3RNGCore(const std::string& created_by, const std::string& comment = ""):
    NR3RNGCore(0, created_by, comment) { /* nothing to see here */ }
  NR3RNGCore(const ix::math::rng::NR3RNGCoreData& proto, bool restore_state = false,
    const std::string& created_by = "", const std::string& comment = "");
  ~NR3RNGCore();

  uint64_t uniform_uint64() override {
    calls_++;
    u_ = u_*C_NR3_U_MUL + C_NR3_U_ADD;
    v_ ^= v_ >> C_NR3_V_SHIFT1;
    v_ ^= v_ << C_NR3_V_SHIFT2;
    v_ ^= v_ >> C_NR3_V_SHIFT3;
    // Slight inefficiency in next statement - 32bit multiplies would be OK
    w_ = C_NR3_W_MUL*(w_ & 0xFFFFFFFF) + (w_ >> C_NR3_W_SHIFT1);
    uint64_t x = u_ ^ (u_ << C_NR3_X_SHIFT1);
    x ^= x >> C_NR3_X_SHIFT2;
    x ^= x << C_NR3_X_SHIFT3;
    return (x+v_)^w_;
  }

  void save_to_proto(ix::math::rng::NR3RNGCoreData* data) const;
  void save_to_proto(ix::math::rng::RNGCoreData* proto) const override;

  static ix_core_data_type* mutable_core_data(ix::math::rng::RNGCoreData* proto) {
    return proto->mutable_nr3_core(); }
  static const ix_core_data_type& core_data(const ix::math::rng::RNGCoreData& proto)
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

  Ranlux48RNGCore(uint64_t seed = 0,
    const std::string& created_by = "", const std::string& comment = "");
  Ranlux48RNGCore(const std::string& created_by, const std::string& comment = ""):
    Ranlux48RNGCore(0, created_by, comment) { /* nothing to see here */ }
  Ranlux48RNGCore(const ix::math::rng::Ranlux48RNGCoreData& proto,
    bool restore_state = false,
    const std::string& created_by = "", const std::string& comment = "");
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

  void save_to_proto(ix::math::rng::Ranlux48RNGCoreData* proto) const;
  void save_to_proto(ix::math::rng::RNGCoreData* proto) const override;

  static ix_core_data_type* mutable_core_data(ix::math::rng::RNGCoreData* proto) {
    return proto->mutable_ranlux48_core(); }
  static const ix_core_data_type& core_data(const ix::math::rng::RNGCoreData& proto)
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

  MT19937RNGCore(uint64_t seed = 0,
    const std::string& created_by = "", const std::string& comment = "");
  MT19937RNGCore(const std::string& created_by, const std::string& comment = ""):
    MT19937RNGCore(0, created_by, comment) { /* nothing to see here */ }
  MT19937RNGCore(const ix::math::rng::STLRNGCoreData& proto,
    bool restore_state = false,
    const std::string& created_by = "", const std::string& comment = "");
  ~MT19937RNGCore();

  uint64_t uniform_uint64() override { gen_calls_++; return gen_(); }

  void save_to_proto(ix::math::rng::STLRNGCoreData* proto) const;
  void save_to_proto(ix::math::rng::RNGCoreData* proto) const override;

  static ix_core_data_type* mutable_core_data(ix::math::rng::RNGCoreData* proto) {
    return proto->mutable_mt19937_core(); }
  static const ix_core_data_type& core_data(const ix::math::rng::RNGCoreData& proto)
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

  NR3_EmulateSIMD_RNGCore(uint64_t seed,
      const std::string& created_by = "", const std::string& comment = ""):
    RNGCore(),
    seed_(seed>0 ? seed : RNG::nonzero_uint64_from_random_device())
  {
    std::mt19937_64 gen(seed_);
    std::string sub_comment = comment;
    if(not created_by.empty()) {
      sub_comment = "Deligated from: " + created_by;
      if(not comment.empty()) {
        sub_comment += ". Comment: " + comment;
      }
    }
    for(unsigned i=0;i<NSTREAM;i++)
      core_[i] = new NR3RNGCore(gen(), __PRETTY_FUNCTION__, sub_comment);
    write_provenance(created_by, comment);
  }

  NR3_EmulateSIMD_RNGCore(uint64_t seeds[NSTREAM],
      const std::string& created_by = "", const std::string& comment = ""):
    RNGCore(), seed_(0)
  {
    std::string sub_comment = comment;
    if(not created_by.empty()) {
      sub_comment = "Deligated from: " + created_by;
      if(not comment.empty()) {
        sub_comment += ". Comment: " + comment;
      }
    }
    for(unsigned i=0;i<NSTREAM;i++)
      core_[i] = new NR3RNGCore(seeds[i], __PRETTY_FUNCTION__, sub_comment);
    write_provenance(created_by, comment);
  }

  NR3_EmulateSIMD_RNGCore(const ix::math::rng::NR3_SIMD_RNGCoreData& proto,
      bool restore_state = false,
      const std::string& created_by = "", const std::string& comment = ""):
    NR3_EmulateSIMD_RNGCore(proto.seed())
  {
    if(proto.vec_stream_seed_size()!=0)
    {
      if(proto.vec_stream_seed_size() != NSTREAM)
        throw std::runtime_error("NR3_EmulateSIMD_RNGCore: need " +
          std::to_string(NSTREAM) + " seeds to restore state.");
      if(restore_state and proto.state_saved()) {
        calls_ = proto.calls();
        if(proto.vec_u_size() != NSTREAM)
          throw std::runtime_error("NR3_EmulateSIMD_RNGCore: need " +
            std::to_string(NSTREAM) + " u elements to restore state.");
        if(proto.vec_v_size() != NSTREAM)
          throw std::runtime_error("NR3_EmulateSIMD_RNGCore: need " +
            std::to_string(NSTREAM) + " v elements to restore state.");
        if(proto.vec_w_size() != NSTREAM)
          throw std::runtime_error("NR3_EmulateSIMD_RNGCore: need " +
            std::to_string(NSTREAM) + " w elements to restore state.");
      }
      for(unsigned i=0; i<NSTREAM; i++)
      {
        calin::ix::math::rng::NR3RNGCoreData core_data;
        core_data.set_seed(proto.vec_stream_seed(i));
        if(restore_state and proto.state_saved()) {
          core_data.set_calls(calls_);
          core_data.set_state_saved(true);
          core_data.set_u(proto.vec_u(i));
          core_data.set_v(proto.vec_v(i));
          core_data.set_w(proto.vec_w(i));
        }
        delete core_[NSTREAM-i-1];
        core_[NSTREAM-i-1] = new NR3RNGCore(core_data, restore_state);
      }
      ndev_ = NSTREAM;
      if(restore_state and proto.state_saved()) {
        if(proto.dev_size() > int(NSTREAM))
          throw std::runtime_error("NR3_EmulateSIMD_RNGCore: need a most " +
            std::to_string(NSTREAM) + " saved deviates to restore state.");
        for(int i=0; i<proto.dev_size(); i++) {
          vec_dev_[--ndev_] = proto.dev(i);
        }
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

  void save_to_proto(ix::math::rng::NR3_SIMD_RNGCoreData* data) const
  {
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

  void save_to_proto(ix::math::rng::RNGCoreData* proto) const override
  {
    save_to_proto(proto->mutable_nr3_simd_emu4_core());
  }

  static ix_core_data_type* mutable_core_data(ix::math::rng::RNGCoreData* proto) {
    return proto->mutable_nr3_simd_emu4_core(); }
  static const ix_core_data_type& core_data(const ix::math::rng::RNGCoreData& proto)
  { return proto.nr3_simd_emu4_core(); }

private:
  NR3RNGCore* core_[NSTREAM];
  uint64_t seed_ = 0;
  uint64_t calls_ = 0;
  uint64_t vec_dev_[NSTREAM];
  unsigned ndev_ = NSTREAM;
};

} } } // namespace calin::math::rng
