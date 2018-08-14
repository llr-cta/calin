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
#include <util/log.hpp>
#include <provenance/chronicle.hpp>

namespace calin { namespace math { namespace rng {

template<typename VCLArchitecture> class NR3_VCLRNGCore;

template<typename VCLArchitecture> class VCLRNGCore: public VCLArchitecture
{
public:
  using typename VCLArchitecture::uint64_vt;

  virtual ~VCLRNGCore() { }

  virtual uint64_vt uniform_uint64() = 0;

  virtual void save_to_proto(calin::ix::math::rng::VCLRNGCoreData* proto) const = 0;

  calin::ix::math::rng::VCLRNGCoreData* as_proto() const
  {
    auto* proto = new calin::ix::math::rng::VCLRNGCoreData;
    save_to_proto(proto);
    return proto;
  }

  static VCLRNGCore<VCLArchitecture>* create_from_proto(const ix::math::rng::VCLRNGCoreData& proto,
    bool restore_state = false,
    const std::string& created_by = "", const std::string& comment = "")
  {
    if(proto.vcl_architecture() != VCLArchitecture::vec_bits) {
      calin::util::log::LOG(calin::util::log::ERROR)
        << calin::util::vcl::templated_class_name<VCLArchitecture>("VCLRNGCore")
        << "Incorrect architecture : " << proto.vcl_architecture()
        << " != " << VCLArchitecture::vec_bits;
      return new NR3_VCLRNGCore<VCLArchitecture>(0, created_by, comment);
    }

    switch(proto.core_case()) {
      case ix::math::rng::VCLRNGCoreData::kNr3VclCore:
        return new NR3_VCLRNGCore<VCLArchitecture>(proto.nr3_vcl_core(), restore_state, created_by, comment);
      default:
        calin::util::log::LOG(calin::util::log::ERROR)
          << calin::util::vcl::templated_class_name<VCLArchitecture>("VCLRNGCore")
          << "Unrecognised RNG type case : " << proto.core_case();
      case 0:
        return new NR3_VCLRNGCore<VCLArchitecture>(0, created_by, comment);
    }
  }

protected:
  void write_provenance(const std::string& created_by, const std::string& comment = "")
  {
    const ix::math::rng::VCLRNGCoreData* proto = this->as_proto();
    calin::provenance::chronicle::register_vcl_rng_core(*proto, created_by, comment);
    delete proto;
  }
};

template<typename VCLArchitecture> class VCLRNG: public VCLArchitecture
{
public:
  using typename VCLArchitecture::uint32_vt;
  using typename VCLArchitecture::int32_vt;
  using typename VCLArchitecture::uint64_vt;
  using typename VCLArchitecture::int64_vt;
  using typename VCLArchitecture::float_vt;
  using typename VCLArchitecture::double_vt;
  using typename VCLArchitecture::float_bvt;
  using typename VCLArchitecture::double_bvt;

  enum class CoreType { NR3 };
  VCLRNG(uint64_t seed, CoreType core_type = CoreType::NR3,
    const std::string& created_by = "", const std::string& comment = "")
  {
    switch(core_type) {
    case CoreType::NR3:
      core_ = new NR3_VCLRNGCore<VCLArchitecture>(seed, created_by, comment);
      break;
    }
  }

  VCLRNG(CoreType core_type = CoreType::NR3, const std::string& created_by = "",
      const std::string& comment = ""):
    VCLRNG(0, core_type, created_by, comment) { /* nothing to see here */ }
  VCLRNG(const std::string& created_by, const std::string& comment = "",
      uint64_t seed = 0, CoreType core_type = CoreType::NR3):
    VCLRNG(seed, core_type, created_by, comment) { /* nothing to see here */ }
  VCLRNG(uint64_t seed, const std::string& created_by, const std::string& comment = "",
      CoreType core_type = CoreType::NR3):
    VCLRNG(seed, core_type, created_by, comment) { /* nothing to see here */ }

  VCLRNG(VCLRNGCore<VCLArchitecture>* core, bool adopt_core = false):
      core_(core), adopt_core_(adopt_core) { /* nothing to see here */ }

  VCLRNG(const calin::ix::math::rng::VCLRNGCoreData& proto, bool restore_state = false,
      const std::string& created_by = "", const std::string& comment = ""):
    core_(VCLRNGCore<VCLArchitecture>::create_from_proto(proto,
      restore_state, created_by, comment)) { /* nothing to see here */ }

  ~VCLRNG() {
    if(adopt_core_)delete core_;
  }

  //   void save_to_proto(ix::math::rng::RNGData* proto) const;
  //
  //   ix::math::rng::RNGData* as_proto() const {
  //     auto* proto = new ix::math::rng::RNGData;
  //     save_to_proto(proto); return proto; }

  uint64_vt uniform_uint64() { return core_->uniform_uint64(); }
  int64_vt uniform_int64() { return int64_vt(core_->uniform_uint64()); }

  uint32_vt uniform_uint32() { return uint32_vt(core_->uniform_uint64()); }
  int32_vt uniform_int32() { return int32_vt(core_->uniform_uint64()); }

  // Generate float/double deviates in the range [-0.5*scale+offset, 0.5*scale+offset)
  // The default values therefore create deviates in the "standard" range of [0,1)
  double_vt uniform_double_gen(const double scale = 1.0, const double offset = 0.5) {
    return mul_add(to_double(uniform_int64()), C_U64_TO_DBL*scale, offset);
  }

  float_vt uniform_float_gen(const float scale = 1.0f, const float offset = 0.5f) {
    return mul_add(to_float(uniform_int32()), C_U32_TO_FLT*scale, offset);
  }

  // Generate zero-centered float/double deviates in the range [-0.5*scale, 0.5*scale]
  // The default values therefore create deviates in the range of [-0.5,0.5]
  double_vt uniform_double_zc(const double scale = 1.0) {
    return to_double(uniform_int64()) * (C_U64_TO_DBL*scale);
  }

  float_vt uniform_float_zc(const float scale = 1.0f) {
    return to_float(uniform_int32()) * (C_U32_TO_FLT*scale);
  }

  // Generate standard float/double deviates in the range [0, scale]
  // The default values therefore create deviates in the range of [0,1]
  double_vt uniform_double(const double scale = 1.0) {
    const double multiplier = C_U64_TO_DBL*scale;
    double_vt x = to_double(uniform_int64());
    // return select(x<0, nmul_add(x,multiplier,0.5*scale), x*multiplier);
    return mul_add(x, multiplier, double_vt(scale) & (x<0));
}

  float_vt uniform_float(const float scale = 1.0f) {
    const float multiplier = C_U32_TO_FLT*scale;
    float_vt x = to_float(uniform_int32());
    // return select(x<0, nmul_add(x,multiplier,0.5*scale), x*multiplier);
    return mul_add(x, multiplier, float_vt(scale) & (x<0));
  }

  void uniform_by_type(uint64_vt& x) { x = uniform_uint64(); }
  void uniform_by_type(uint32_vt& x) { x = uniform_uint32(); }
  void uniform_by_type(double_vt& x) { x = uniform_double(); }
  // void uniform_by_type(double_vt& x, double scale) { x = uniform_double(scale); }
  // void uniform_by_type(double_vt& x, double scale, double offset) {
  //   x = uniform_double(scale,offset); }
  void uniform_by_type(float_vt& x) { x = uniform_float(); }
  // void uniform_by_type(float_vt& x, float scale) { x = uniform_float(scale); }
  // void uniform_by_type(float_vt& x, float scale, float offset) {
  //   x = uniform_float(scale, offset); }

  double_vt exponential_double() {
    return -calin::util::vcl::log(max(uniform_double(),DBL_MIN));
  }

  double_vt exponential_double(const double_vt& scale) {
    return scale * exponential_double();
  }

  float_vt exponential_float() {
    return -calin::util::vcl::log(max(uniform_float(),FLT_MIN));
  }

  float_vt exponential_float(const float_vt& scale) {
    return scale * exponential_float();
  }

  // This is a stripped-down version of "sincos_f" in VCL "vectormath_trig.h"
  // The code for range reduction and quadrant selection is removed since the
  // argument is in a controlled range.
  void sincos_float(float_vt& s, float_vt& c)
  {
    const float P0sinf = -1.6666654611E-1f;
    const float P1sinf = 8.3321608736E-3f;
    const float P2sinf = -1.9515295891E-4f;

    const float P0cosf = 4.166664568298827E-2f;
    const float P1cosf = -1.388731625493765E-3f;
    const float P2cosf = 2.443315711809948E-5f;

    // Sin and cos calculated from uniform 32-bit integers
    // Note : Lowest two bits used below to (1) flip cos sign bit and
    //        (2) exchange sin and cos polynomials. These bits are masked so
    //        that they are not used in the argument.
    uint32_vt uix = uniform_uint32();

    // Random deviate between -pi/4 and pi/4
    float_vt x = C_U32_TO_FLT*M_PI_2*to_float(int32_vt(uix & (~0x3)));

    // Polynomial for sin and cos valid between -pi/4 and pi/4
    float_vt x2 = x * x;
    s = polynomial_2(x2, P0sinf, P1sinf, P2sinf);
    c = polynomial_2(x2, P0cosf, P1cosf, P2cosf);
    s = mul_add(s, x*x2, x);
    c = mul_add(c, x2*x2, nmul_add(0.5f, x2, 1.0f));

    // Second lowest bit of original uint32 used to swap sin and cos.
    float_bvt swap_mask = float_bvt((uix & 0x2) == 0);

    // Lowest bit of original uint32 used to flip sign of cos - XOR the sign bit
    uix <<= 31;
    c ^= reinterpret_f(uix);

    float_vt c2 = select(swap_mask, c, s);
    s = select(swap_mask, s, c);
    c = c2;
  }

  void sincos_double(double_vt& s, double_vt& c)
  {
    // See comments from sincos_float

    const double P0sin = -1.66666666666666307295E-1;
    const double P1sin = 8.33333333332211858878E-3;
    const double P2sin = -1.98412698295895385996E-4;
    const double P3sin = 2.75573136213857245213E-6;
    const double P4sin = -2.50507477628578072866E-8;
    const double P5sin = 1.58962301576546568060E-10;

    const double P0cos = 4.16666666666665929218E-2;
    const double P1cos = -1.38888888888730564116E-3;
    const double P2cos = 2.48015872888517045348E-5;
    const double P3cos = -2.75573141792967388112E-7;
    const double P4cos = 2.08757008419747316778E-9;
    const double P5cos = -1.13585365213876817300E-11;

    uint64_vt uix = uniform_uint64();

    double_vt x = C_U64_TO_DBL*M_PI_2*to_double(int64_vt(uix & (~UINT64_C(0x3))));

    double_vt x2 = x * x;
    s = polynomial_5(x2, P0sin, P1sin, P2sin, P3sin, P4sin, P5sin);
    c = polynomial_5(x2, P0cos, P1cos, P2cos, P3cos, P4cos, P5cos);
    s = mul_add(x * x2, s, x);
    c = mul_add(x2 * x2, c, nmul_add(x2, 0.5, 1.0));

    double_bvt swap_mask = double_bvt((uix & UINT64_C(0x2)) == UINT64_C(0));

    uix <<= 63;
    c ^= reinterpret_d(uix);

    double_vt c2 = select(swap_mask, c, s);
    s = select(swap_mask, s, c);
    c = c2;
  }

  void normal_two_float_bm(float_vt& x1, float_vt& x2)
  {
    // double_vt r1 = sqrt(exponential_double(2.0));
    // double_vt r2 = sqrt(exponential_double(2.0));
    // float_vt r = compress(r1,r2);
    float_vt r = sqrt(exponential_float(2.0));
    float_vt s, c;
    sincos_float(s, c);
    x1 = r*c;
    x2 = r*s;
  }

  void normal_two_double_bm(double_vt& x1, double_vt& x2)
  {
    double_vt r = sqrt(exponential_double(2.0));
    double_vt s, c;
    sincos_double(s, c);
    x1 = r*c;
    x2 = r*s;
  }

  static uint64_vt uint64_from_seed(uint64_t seed = 0)
  {
    if(seed == 0)seed = RNG::uint64_from_random_device();
    std::mt19937_64 gen(seed);
    uint64_t u64[VCLArchitecture::num_uint64];
    for(unsigned i=0; i<VCLArchitecture::num_uint64; i++)u64[i] = gen();
    uint64_vt x;
    x.load(u64);
    return x;
  }

  static uint64_vt uint64_from_random_device()
  {
    std::random_device gen;
    static_assert(sizeof(std::random_device::result_type)==sizeof(uint32_t),
                  "std::random_device::result_type is not 32 bits");
    uint32_t u32[VCLArchitecture::num_int32];
    for(unsigned i=0;i<VCLArchitecture::num_int32;i++) {
      u32[i] = gen();
    }
    uint32_vt x;
    x.load(u32);
    return uint64_vt(x);
  }

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
private:
  VCLRNGCore<VCLArchitecture>* core_ = nullptr;
  bool adopt_core_ = true;
};

// =============================================================================
//
// Specific random number generators
//
// =============================================================================

template<typename VCLArchitecture> class NR3_VCLRNGCore:
  public VCLRNGCore<VCLArchitecture>
{
public:
  using typename VCLArchitecture::uint64_vt;

  typedef calin::ix::math::rng::NR3_SIMD_RNGCoreData ix_core_data_type;

  NR3_VCLRNGCore(uint64_t seed = 0,
      const std::string& created_by = "", const std::string& comment = ""):
    VCLRNGCore<VCLArchitecture>(),
    seed_(seed>0 ? seed : RNG::nonzero_uint64_from_random_device()),
    sequence_seeds_(VCLRNG<VCLArchitecture>::uint64_from_seed(seed_))
  {
    init();
    this->write_provenance(created_by, comment);
  }

  NR3_VCLRNGCore(uint64_t seeds[VCLArchitecture::num_uint64],
      const std::string& created_by = "", const std::string& comment = ""):
    VCLRNGCore<VCLArchitecture>(), seed_(0)
  {
    sequence_seeds_.load(seeds);
    init();
    this->write_provenance(created_by, comment);
  }

  NR3_VCLRNGCore(const ix::math::rng::NR3_SIMD_RNGCoreData& proto,
      bool restore_state = false,
      const std::string& created_by = "", const std::string& comment = ""):
    VCLRNGCore<VCLArchitecture>(), seed_(proto.seed())
  {
    // We allow either zero of num_uint64 seeds.. zero means we do "seed only" reinit
    if(proto.vec_stream_seed_size() == 0)
    {
      if(seed_ == 0)seed_ = RNG::nonzero_uint64_from_random_device();
      sequence_seeds_ = VCLRNG<VCLArchitecture>::uint64_from_seed(seed_);
      init();
    }
    else if(proto.vec_stream_seed_size() == VCLArchitecture::num_uint64)
    {
      sequence_seeds_.load(proto.vec_stream_seed().data());
      if(restore_state and proto.state_saved())
      {
        if(proto.vec_u_size() != VCLArchitecture::num_uint64
            or proto.vec_v_size() != VCLArchitecture::num_uint64
            or proto.vec_w_size() != VCLArchitecture::num_uint64)
          throw std::runtime_error(
            calin::util::vcl::templated_class_name<VCLArchitecture>("NR3_VCLRNGCore")
            + ": saved state vectors must have "
            + std::to_string(VCLArchitecture::num_uint64) + " elements");
        calls_ = proto.calls();
        u_.load(proto.vec_u().data());
        v_.load(proto.vec_v().data());
        w_.load(proto.vec_w().data());
      } else {
        init();
      }
    }
    else
    {
      throw std::runtime_error(
        calin::util::vcl::templated_class_name<VCLArchitecture>("NR3_VCLRNGCore")
        + ": saved seed vectors must have "
        + std::to_string(VCLArchitecture::num_uint64) + " elements");
    }
  }

  virtual ~NR3_VCLRNGCore() { }

  uint64_vt uniform_uint64() override {
    calls_++;
    u_ = u_*C_NR3_U_MUL + C_NR3_U_ADD;
    v_ ^= v_ >> C_NR3_V_SHIFT1;
    v_ ^= v_ << C_NR3_V_SHIFT2;
    v_ ^= v_ >> C_NR3_V_SHIFT3;
    // This is inefficient as the full 64bit multiply is not needed
    // w_ = C_NR3_W_MUL*(w_ & 0xFFFFFFFF) + (w_ >> C_NR3_W_SHIFT1);
    w_ = calin::util::vcl::multiply_low_32bit(w_, C_NR3_W_MUL) + (w_ >> C_NR3_W_SHIFT1);
    uint64_vt x = u_ ^ (u_ << C_NR3_X_SHIFT1);
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
  void init()
  {
    u_ = sequence_seeds_^v_; uniform_uint64();
    v_ = u_; uniform_uint64();
    w_ = v_; uniform_uint64();
  }

  uint64_t seed_;
  uint64_t calls_ = 0;
  uint64_vt sequence_seeds_;
  uint64_vt u_ = C_NR3_U_INIT;
  uint64_vt v_ = C_NR3_V_INIT;
  uint64_vt w_ = C_NR3_W_INIT;
};

} } } // namespace calin::math::rng
