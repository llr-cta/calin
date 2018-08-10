/*

   calin/math/rng_vcl.hpp -- Stephen Fegan -- 2018-08-09

   Random number generators for VCL vector types.

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

#include <math/rng.hpp>
#include <util/vcl.hpp>
#include <provenance/chronicle.hpp>

namespace calin { namespace math { namespace rng {

template<typename VCLArchitecture> class VCLRNGCore: public VCLArchitecture
{
public:
  virtual ~VCLRNGCore() { }

  virtual typename VCLArchitecture::uint64_vec_type uniform_uint64() = 0;

  virtual void save_to_proto(calin::ix::math::rng::VCLRNGCoreData* proto) const = 0;

  calin::ix::math::rng::VCLRNGCoreData* as_proto() const
  {
    auto* proto = new calin::ix::math::rng::VCLRNGCoreData;
    save_to_proto(proto);
    return proto;
  }

  // static RNGCore* create_from_proto(const ix::math::rng::RNGCoreData& proto,
  //   bool restore_state = false,
  //   const std::string& created_by = "", const std::string& comment = "")
  // {
  //
  // }

protected:
  void write_provenance(const std::string& created_by, const std::string& comment = "")
  {
    const ix::math::rng::VCLRNGCoreData* proto = this->as_proto();
    calin::provenance::chronicle::register_vcl_rng_core(*proto, created_by, comment);
    delete proto;
  }
};

// class RNG
// {
// public:
//   enum class CoreType { NR3, RANLUX48, MT19937 /*, DEV_RANDOM*/ };
//   RNG(CoreType core_type = CoreType::NR3, const std::string& created_by = "");
//   RNG(const std::string& created_by, CoreType core_type = CoreType::NR3):
//     RNG(core_type, created_by) { /* nothing to see here */ }
//   RNG(uint64_t seed, CoreType core_type = CoreType::NR3, const std::string& created_by = "");
//   RNG(uint64_t seed, const std::string& created_by, CoreType core_type = CoreType::NR3):
//     RNG(seed, core_type, created_by) { /* nothing to see here */ }
//   RNG(RNGCore* core, bool adopt_core = false);
//   RNG(const ix::math::rng::RNGData& proto, bool restore_state = false,
//     const std::string& created_by = "");
//
//   ~RNG();
//
//   void save_to_proto(ix::math::rng::RNGData* proto) const;
//
//   ix::math::rng::RNGData* as_proto() const {
//     auto* proto = new ix::math::rng::RNGData;
//     save_to_proto(proto); return proto; }
//
//   uint64_t uniform_uint64() { return core_->uniform_uint64(); }
//
//   uint32_t uniform_uint32() { return uint32_t(uniform_uint64()); }
//
//   double uniform_double() {
//     return C_U64_TO_DBL * double(uniform_uint64());
//   }
//
//   float uniform_float() {
//     return C_U32_TO_FLT * float(uniform_uint32());
//   }
//
//   void uniform_by_type(uint64_t& x) { x = uniform_uint64(); }
//   void uniform_by_type(uint32_t& x) { x = uniform_uint32(); }
//   void uniform_by_type(double& x) { x = uniform_double(); }
//   void uniform_by_type(float& x) { x = uniform_float(); }
//
//   double uniform() { return uniform_double(); }
//   double exponential() { return -std::log(uniform()); }
//   double exponential(double mean) { return -mean*std::log(uniform()); }
//   double normal();
//   double normal(double mean, double sigma) { return mean+normal()*sigma; }
//   double gamma_by_alpha_and_beta(double alpha, double beta);
//   double gamma_by_mean_and_sigma(const double mean, const double sigma) {
//     double b = mean/(sigma*sigma);
//     return gamma_by_alpha_and_beta(mean*b,b); }
//   int poisson(double lambda);
//   int polya(double mean, double non_poisson_sigma) {
//     return poisson(gamma_by_mean_and_sigma(mean, non_poisson_sigma)); }
//   int binomial(double pp, int n);
//
//   static uint32_t uint32_from_random_device();
//   static uint64_t uint64_from_random_device();
//
//   static uint64_t nonzero_uint64_from_random_device() { uint64_t x;
//     do x = uint64_from_random_device(); while(x==0ULL); return x; }
//   static uint32_t nonzero_uint32_from_random_device() { uint32_t x;
//     do x = uint32_from_random_device(); while(x==0U); return x; }
//
//   double inverse_cdf(const std::vector<std::pair<double,double>> &inv_cdf);
//   static void generate_inverse_cdf(std::vector<std::pair<double,double>> &cdf,
//     unsigned nbins=0);
//
//   static constexpr uint64_t std_test_seed = 12939; // essential supply
//
// private:
//   RNGCore* core_;
//   bool adopt_core_;
//
//   // Cached values
//   bool bm_hascached_ = false;
//   double bm_cachedval_ = 0.0;
//
//   // Speedup caches
//   double   poi_lambdaexp_ = 1.0/M_E;
//   double   poi_lambdasrt_ = 1.0;
//   double   poi_lambdalog_ = 0.0;
//   double   poi_lambdaold_ = 1.0;
//
//   int      bin_nold_      = -1;
//   double   bin_pold_      = -1;
//   double   bin_pc_        = 0.0;
//   double   bin_plog_      = 0.0;
//   double   bin_pclog_     = 0.0;
//   unsigned bin_en_        = 0;
//   double   bin_oldg_      = 0.0;
// };

// =============================================================================
//
// Specific random number generators
//
// =============================================================================

template<typename VCLArchitecture> class NR3_VCLRNGCore:
  public VCLRNGCore<VCLArchitecture>
{
public:
  typedef calin::ix::math::rng::NR3_SIMD_RNGCoreData ix_core_data_type;

  NR3_VCLRNGCore(uint64_t seed = 0,
      const std::string& created_by = "", const std::string& comment = ""):
    VCLRNGCore<VCLArchitecture>(),
    seed_(seed>0 ? seed : RNG::nonzero_uint64_from_random_device()),
    sequence_seeds_(), u_(C_NR3_U_INIT), v_(C_NR3_V_INIT), w_(C_NR3_W_INIT)
  {
    std::mt19937_64 gen(seed_);
    uint64_t seeds[VCLArchitecture::num_uint64];
    for(unsigned i=0; i<VCLArchitecture::num_uint64; i++)seeds[i] = gen();
    init(seeds);
    this->write_provenance(created_by, comment);
  }

  NR3_VCLRNGCore(uint64_t seeds[VCLArchitecture::num_uint64],
      const std::string& created_by = "", const std::string& comment = ""):
    VCLRNGCore<VCLArchitecture>(),
    seed_(0), sequence_seeds_(),
    u_(C_NR3_U_INIT), v_(C_NR3_V_INIT), w_(C_NR3_W_INIT)
  {
    init(seeds);
    this->write_provenance(created_by, comment);
  }

  NR3_VCLRNGCore(const ix::math::rng::NR3_SIMD_RNGCoreData& proto, bool restore_state = false,
      const std::string& created_by = "", const std::string& comment = ""):
    VCLRNGCore<VCLArchitecture>(),
    seed_(proto.seed())
  {
    // We allow either zero of num_uint64 seeds.. zero means we do "seed only" reinit
    if(proto.vec_stream_seed_size() == 0)
    {
      if(seed_ == 0)seed_ = RNG::nonzero_uint64_from_random_device();
      std::mt19937_64 gen(seed_);
      uint64_t seeds[VCLArchitecture::num_uint64];
      for(unsigned i=0; i<VCLArchitecture::num_uint64; i++)seeds[i] = gen();
      init(seeds);
    }
    else if(proto.vec_stream_seed_size() == VCLArchitecture::num_uint64)
    {
      init(proto.vec_stream_seed().data());
      if(restore_state and proto.state_saved())
      {
        if(proto.vec_u_size() != VCLArchitecture::num_uint64
            or proto.vec_v_size() != VCLArchitecture::num_uint64
            or proto.vec_w_size() != VCLArchitecture::num_uint64)
          throw std::runtime_error(std::string("NR3_AVX2_RNGCore<")
            + VCLArchitecture::architecture_name + ">: saved state vectors must have "
            + std::to_string(VCLArchitecture::num_uint64) + " elements");
        calls_ = proto.calls();
        u_.load(proto.vec_u().data());
        v_.load(proto.vec_v().data());
        w_.load(proto.vec_w().data());
      }
    }
    else
    {
      throw std::runtime_error(std::string("NR3_AVX2_RNGCore<")
        + VCLArchitecture::architecture_name + ">: saved seed vectors must have "
        + std::to_string(VCLArchitecture::num_uint64) + " elements");
    }
  }

  virtual ~NR3_VCLRNGCore() { }

  typename VCLArchitecture::uint64_vec_type uniform_uint64() override {
    calls_++;
    u_ = u_*C_NR3_U_MUL + C_NR3_U_ADD;
    v_ ^= v_ >> C_NR3_V_SHIFT1;
    v_ ^= v_ << C_NR3_V_SHIFT2;
    v_ ^= v_ >> C_NR3_V_SHIFT3;
    // This is inefficient as the full 64bit multiply is not needed
    // w_ = C_NR3_W_MUL*(w_ & 0xFFFFFFFF) + (w_ >> C_NR3_W_SHIFT1);
    w_ = calin::util::vcl::multiply_low_32bit(w_, C_NR3_W_MUL) + (w_ >> C_NR3_W_SHIFT1);
    typename VCLArchitecture::uint64_vec_type x = u_ ^ (u_ << C_NR3_X_SHIFT1);
    x ^= x >> C_NR3_X_SHIFT2;
    x ^= x << C_NR3_X_SHIFT3;
    return (x+v_)^w_;
  }

  void save_to_proto(calin::ix::math::rng::NR3_SIMD_RNGCoreData* data) const
  {
    data->set_seed(seed_);
    data->set_calls(calls_);
    data->set_state_saved(true);
    for(unsigned i=0; i<VCLArchitecture::num_uint64; i++) {
      data->add_vec_stream_seed(sequence_seeds_[i]);
      data->add_vec_u(u_[i]);
      data->add_vec_v(v_[i]);
      data->add_vec_w(w_[i]);
    }
  }

  void save_to_proto(calin::ix::math::rng::VCLRNGCoreData* proto) const override
  {
    proto->set_vcl_architecture(VCLArchitecture::vec_bits);
    auto* data = proto->mutable_nr3_vcl_core();
    save_to_proto(data);
  }

  // static ix_core_data_type* mutable_core_data(ix::math::rng::RNGCoreData* proto) {
  //   return proto->mutable_nr3_core(); }
  // static const ix_core_data_type& core_data(const ix::math::rng::RNGCoreData& proto)
  // { return proto.nr3_core(); }

#ifndef SWIG
  static void* operator new(size_t nbytes) {
    void* p = nullptr;
    if(::posix_memalign(&p, 32, nbytes)==0) {
      return p;
    }
    throw std::bad_alloc();
  }
  static void* operator new(size_t nbytes, void* p) {
    return p;
  }
  static void operator delete(void *p) {
    free(p);
  }
#endif

private:
  void init(const uint64_t* seeds)
  {
    sequence_seeds_.load(seeds);
    u_ = sequence_seeds_^v_; uniform_uint64();
    v_ = u_; uniform_uint64();
    w_ = v_; uniform_uint64();
  }

  uint64_t seed_;
  uint64_t calls_ = 0;
  typename VCLArchitecture::uint64_vec_type sequence_seeds_;
  typename VCLArchitecture::uint64_vec_type u_;
  typename VCLArchitecture::uint64_vec_type v_;
  typename VCLArchitecture::uint64_vec_type w_;
};

} } } // namespace calin::math::rng
