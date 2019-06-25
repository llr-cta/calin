/*

   calin/math/rng_vcl.hpp -- Stephen Fegan -- 2018-08-09

   Random number generators for VCL vector types.

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
      return new NR3_VCLRNGCore<VCLArchitecture>(created_by, comment);
    }

    switch(proto.core_case()) {
      case ix::math::rng::VCLRNGCoreData::kNr3VclCore:
        return new NR3_VCLRNGCore<VCLArchitecture>(proto.nr3_vcl_core(), restore_state, created_by, comment);
      default:
        calin::util::log::LOG(calin::util::log::ERROR)
          << calin::util::vcl::templated_class_name<VCLArchitecture>("VCLRNGCore")
          << "Unrecognised RNG type case : " << proto.core_case();
      case 0:
        return new NR3_VCLRNGCore<VCLArchitecture>(created_by, comment);
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
  using typename VCLArchitecture::uint64_bvt;
  using typename VCLArchitecture::float_bvt;
  using typename VCLArchitecture::double_bvt;
  using typename VCLArchitecture::Vector3f_vt;
  using typename VCLArchitecture::Vector3d_vt;

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

  void uniform_uint(uint64_vt& x) { x = uniform_uint64(); }
  void uniform_uint(uint32_vt& x) { x = uniform_uint32(); }

  void uniform_int(int64_vt& x) { x = uniform_uint64(); }
  void uniform_int(int32_vt& x) { x = uniform_uint32(); }

  // Generate float/double deviates in the range [-0.5*scale+offset, 0.5*scale+offset)
  // The default values therefore create deviates in the "standard" range of [0,1)
  double_vt uniform_double_gen(const double_vt& scale = 1.0, const double_vt& offset = 0.5) {
    return mul_add(to_double(uniform_int64()), C_U64_TO_DBL*scale, offset);
  }

  float_vt uniform_float_gen(const float_vt& scale = 1.0f, const float_vt& offset = 0.5f) {
    return mul_add(to_float(uniform_int32()), C_U32_TO_FLT*scale, offset);
  }

  void uniform_real_gen(float_vt& x, const float_vt& scale = 1.0f, const float_vt& offset = 0.5f) {
    x = uniform_float_gen(scale, offset);
  }

  void uniform_real_gen(double_vt& x, const double_vt& scale = 1.0f, const double_vt& offset = 0.5f) {
    x = uniform_double_gen(scale, offset);
  }

  // Generate zero-centered float/double deviates in the range [-0.5*scale, 0.5*scale]
  // The default values therefore create deviates in the range of [-0.5,0.5]
  double_vt uniform_double_zc(const double_vt& scale = 1.0) {
    return to_double(uniform_int64()) * (C_U64_TO_DBL*scale);
  }

  float_vt uniform_float_zc(const float_vt& scale = 1.0f) {
    return to_float(uniform_int32()) * (C_U32_TO_FLT*scale);
  }

  void uniform_real_zc(float_vt& x, const float_vt& scale = 1.0f) {
    x = uniform_float_zc(scale);
  }

  void uniform_real_zc(double_vt& x, const double_vt& scale = 1.0f) {
    x = uniform_double_zc(scale);
  }

  // Generate standard float/double deviates in the range [0, scale]
  // The default values therefore create deviates in the range of [0,1]
  double_vt uniform_double(const double scale = 1.0) {
    const double multiplier = C_U64_TO_DBL*scale;
    double_vt x = to_double(uniform_int64());
    // return select(x<0, nmul_add(x,multiplier,0.5*scale), x*multiplier);
    return mul_add(x, multiplier, select(x<0, scale, 0.0));
  }

  double_vt uniform_double_53bit(const double scale = 1.0) {
    constexpr uint64_t MASK_HI = (1ULL<<52)-1;
    constexpr uint64_t EXP_HI = 1023ULL<<52;
    uint64_vt xui = uniform_uint64();
    double_vt xhi = vcl::reinterpret_d((xui&MASK_HI) | EXP_HI);
    return scale*(xhi-1.0);
  }

  double_vt uniform_double_alt(const double scale = 1.0) {
    constexpr uint64_t MASK_HI = (1ULL<<52)-1;
    constexpr uint64_t MASK_LO = ((1ULL<<12)-1)<<40;
    constexpr uint64_t EXP_HI = 1023ULL<<52;
    constexpr uint64_t EXP_LO = (1023ULL-52ULL)<<52;
    uint64_vt xui = uniform_uint64();
    double_vt xlo = vcl::reinterpret_d(((xui>>12)&MASK_LO) | EXP_LO);
    double_vt xhi = vcl::reinterpret_d((xui&MASK_HI) | EXP_HI);
    return scale*((xlo-2.2204460492503130808e-16) + (xhi-1.0));
  }

  float_vt uniform_float(const float scale = 1.0f) {
    const float multiplier = C_U32_TO_FLT*scale;
    float_vt x = to_float(uniform_int32());
    // return select(x<0, nmul_add(x,multiplier,0.5*scale), x*multiplier);
    return mul_add(x, multiplier, select(x<0, scale, 0));
  }

  void uniform_real(float_vt& x, const float_vt& scale = 1.0f) {
    x = uniform_float(scale);
  }

  void uniform_real(double_vt& x, const double_vt& scale = 1.0f) {
    x = uniform_double(scale);
  }


  // void uniform_by_type(uint64_vt& x) { x = uniform_uint64(); }
  // void uniform_by_type(uint32_vt& x) { x = uniform_uint32(); }
  // void uniform_by_type(double_vt& x) { x = uniform_double(); }
  // void uniform_by_type(double_vt& x, double scale) { x = uniform_double(scale); }
  // void uniform_by_type(double_vt& x, double scale, double offset) {
  //   x = uniform_double(scale,offset); }
  // void uniform_by_type(float_vt& x) { x = uniform_float(); }
  // void uniform_by_type(float_vt& x, float scale) { x = uniform_float(scale); }
  // void uniform_by_type(float_vt& x, float scale, float offset) {
  //   x = uniform_float(scale, offset); }

  double_vt exponential_double() {
    double_vt x = uniform_double();
    double_bvt xzero = x==0.0;
    while(horizontal_or(xzero)) {
      x = select(xzero, uniform_double(), x);
      xzero = x==0.0;
    }
    return -calin::util::vcl::log(x);
  }

  double_vt exponential_double(const double_vt& scale) {
    return scale * exponential_double();
  }

  float_vt exponential_float() {
    float_vt x = uniform_float();
    float_bvt xzero = x==0.0f;
    while(horizontal_or(xzero)) {
      x = select(xzero, uniform_float(), x);
      xzero = x==0.0f;
    }
    return -calin::util::vcl::log(x);
  }

  float_vt exponential_float(const float_vt& scale) {
    return scale * exponential_float();
  }

  void exponential_real(double_vt& x) {
    x = exponential_double();
  }

  void exponential_real(double_vt& x, const double_vt& scale) {
    x = exponential_double(scale);
  }

  void exponential_real(float_vt& x) {
    x = exponential_float();
  }

  void exponential_real(float_vt& x, const float_vt& scale) {
    x = exponential_float(scale);
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

  void sincos_real(float_vt& s, float_vt& c)
  {
    sincos_float(s, c);
  }

  void sincos_real(double_vt& s, double_vt& c)
  {
    sincos_double(s, c);
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

  void normal_two_real_bm(float_vt& x1, float_vt& x2)
  {
    normal_two_float_bm(x1, x2);
  }

  void normal_two_real_bm(double_vt& x1, double_vt& x2)
  {
    normal_two_double_bm(x1, x2);
  }

  void uniform_on_unit_sphere_float(float_vt& x, float_vt& y, float_vt& z)
  {
    z = uniform_float_zc(2.0f); // z=cos(theta) uniform from -1 to +1
    const float_vt rho = sqrt(nmul_add(z, z, 1.0f));
    sincos_float(x, y);
    x *= rho;
    y *= rho;
  }

  void uniform_on_unit_sphere_double(double_vt& x, double_vt& y, double_vt& z)
  {
    z = uniform_double_zc(2.0); // z=cos(theta) uniform from -1 to +1
    const double_vt rho = sqrt(nmul_add(z, z, 1.0));
    sincos_double(x, y);
    x *= rho;
    y *= rho;
  }

  void uniform_on_unit_sphere_real(float_vt& x, float_vt& y, float_vt& z)
  {
    uniform_on_unit_sphere_float(x, y, z);
  }

  void uniform_on_unit_sphere_real(double_vt& x, double_vt& y, double_vt& z)
  {
    uniform_on_unit_sphere_float(x, y, z);
  }

  Vector3f_vt uniform_on_unit_sphere_vec_float()
  {
    Vector3f_vt v;
    uniform_on_unit_sphere_float(v.x(), v.y(), v.z());
    return v;
  }

  Vector3d_vt uniform_on_unit_sphere_vec_double()
  {
    Vector3d_vt v;
    uniform_on_unit_sphere_double(v.x(), v.y(), v.z());
    return v;
  }

  void uniform_on_unit_sphere_vec_real(Vector3f_vt& v)
  {
    uniform_on_unit_sphere_float(v.x(), v.y(), v.z());
  }

  void uniform_on_unit_sphere_vec_real(Vector3d_vt& v)
  {
    uniform_on_unit_sphere_double(v.x(), v.y(), v.z());
  }

  float_vt from_inverse_cdf_float(const float* inverse_cdf, unsigned npoints)
  {
    float len = npoints-1;
    float_vt x = uniform_float_gen(len, 0.5*len);
    int32_vt ix = truncate_to_int(x);
    ix = max(ix, 0);
    ix = min(ix, npoints-2);
    float_vt dx = x - to_float(ix);
    float_vt y0 = vcl::lookup<0x40000000>(ix, inverse_cdf);
    float_vt y1 = vcl::lookup<0x40000000>(ix+1, inverse_cdf);
    // std::cout << len << ' ' << x[0] << ' ' << ix[0] << ' ' << dx[0]
    //   << ' ' << y0[0] << ' ' << y1[0] <<  '\n';
    return y0 + (y1-y0)*dx;
  }

  float_vt from_inverse_cdf_float(const std::vector<float> inverse_cdf)
  {
    return from_inverse_cdf_float(inverse_cdf.data(), inverse_cdf.size());
  }

  double_vt from_inverse_cdf_double(const double* inverse_cdf, unsigned npoints)
  {
    double len = npoints-1;
    double_vt x = uniform_double_gen(len, 0.5*len);
    int64_vt ix = truncate_to_int64(x);
    ix = max(ix, 0);
    ix = min(ix, npoints-2);
    double_vt dx = x - to_double(ix);
    double_vt y0 = vcl::lookup<0x40000000>(ix, inverse_cdf);
    double_vt y1 = vcl::lookup<0x40000000>(ix+1, inverse_cdf);
    return y0 + (y1-y0)*dx;
  }

  double_vt from_inverse_cdf_double(const std::vector<double> inverse_cdf)
  {
    return from_inverse_cdf_double(inverse_cdf.data(), inverse_cdf.size());
  }

  void from_inverse_cdf_real(float_vt& x, const float* inverse_cdf, unsigned npoints)
  {
    x = from_inverse_cdf_float(inverse_cdf, npoints);
  }

  void from_inverse_cdf_real(float_vt& x, const std::vector<float>& inverse_cdf)
  {
    x = from_inverse_cdf_float(inverse_cdf);
  }

  void from_inverse_cdf_real(double_vt& x, const double* inverse_cdf, unsigned npoints)
  {
    x = from_inverse_cdf_double(inverse_cdf, npoints);
  }

  void from_inverse_cdf_real(double_vt& x, const std::vector<double>& inverse_cdf)
  {
    x = from_inverse_cdf_double(inverse_cdf);
  }

  static uint64_vt uint64_from_seed(uint64_t seed = 0)
  {
    if(seed == 0)seed = RNG::uint64_from_random_device();
    std::mt19937_64 gen(seed);
    uint64_t u64[VCLArchitecture::num_uint64];
    for(unsigned i=0; i<VCLArchitecture::num_uint64; i++) {
      u64[i] = gen();
    }
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

  static uint64_vt nonzero_uint64_from_random_device()
  {
    uint64_vt x = uint64_from_random_device();
    uint64_bvt xzero = x == 0;
    while(horizontal_or(xzero)) {
      x = select(xzero, uint64_from_random_device(), x);
      xzero = x == 0;
    }
    return x;
  }

  static std::vector<double> generate_inverse_cdf_double(
    const std::vector<std::pair<double,double>> cdf, unsigned nbins)
  {
    if(cdf.size() < 2)
      throw std::runtime_error("generate_inverse_cdf_double : CDF must have at least 2 points.");
    if(cdf.front().second != 0.0)
      throw std::runtime_error("generate_inverse_cdf_double : First CDF bin must have probability zero.");
    if(cdf.back().second != 1.0)
      throw std::runtime_error("generate_inverse_cdf_double : Last CDF bin must have probability one.");
    double p = 0.0;
    for(auto icdf : cdf) {
      if(icdf.second < p)
        throw std::runtime_error("generate_inverse_cdf_double : CDF must be uniformly increasing.");
      p = icdf.second;
    }

    if(nbins < 2)nbins = cdf.size();
    nbins = ((nbins + VCLArchitecture::num_double - 1)/VCLArchitecture::num_double)
      * VCLArchitecture::num_double;

    std::vector<double> inverse_cdf(nbins);

    auto iter = std::upper_bound(cdf.begin(), cdf.end(), 0.0,
      [](double p, const std::pair<double,double>& x) { return p<x.second; });
    inverse_cdf[0] = (iter-1)->first;

    double n_inv = 1.0/double(nbins-1);
    for(unsigned ibin=1; ibin<nbins-1; ibin++) {
      double p = double(ibin)*n_inv;
      iter = std::upper_bound(cdf.begin(), cdf.end(), p,
        [](double p, const std::pair<double,double>& x) { return p<x.second; });
      double dx = (p-(iter-1)->second)/(iter->second-(iter-1)->second);
      inverse_cdf[ibin] = (iter-1)->first + dx*(iter->first - (iter-1)->first);
    }

    iter = std::lower_bound(cdf.begin(), cdf.end(), 1.0,
      [](const std::pair<double,double>& x, double p) { return x.second<p; });
    inverse_cdf[nbins-1] = iter->first;

    return inverse_cdf;
  }

  static std::vector<float> generate_inverse_cdf_float(
    const std::vector<std::pair<double,double>> cdf, unsigned nbins)
  {
    if(cdf.size() < 2)
      throw std::runtime_error("generate_inverse_cdf_double : CDF must have at least 2 points.");
    if(cdf.front().second != 0.0)
      throw std::runtime_error("generate_inverse_cdf_double : First CDF bin must have probability zero.");
    if(cdf.back().second != 1.0)
      throw std::runtime_error("generate_inverse_cdf_double : Last CDF bin must have probability one.");
    double p = 0.0;
    for(auto icdf : cdf) {
      if(icdf.second < p)
        throw std::runtime_error("generate_inverse_cdf_double : CDF must be uniformly increasing.");
      p = icdf.second;
    }

    if(nbins < 2)nbins = cdf.size();
    nbins = ((nbins + VCLArchitecture::num_double - 1)/VCLArchitecture::num_double)
      * VCLArchitecture::num_double;

    std::vector<float> inverse_cdf(nbins);

    auto iter = std::upper_bound(cdf.begin(), cdf.end(), 0.0,
      [](double p, const std::pair<double,double>& x) { return p<x.second; });
    inverse_cdf[0] = (iter-1)->first;

    double n_inv = 1.0/double(nbins-1);
    for(unsigned ibin=1; ibin<nbins-1; ibin++) {
      double p = double(ibin)*n_inv;
      iter = std::upper_bound(cdf.begin(), cdf.end(), p,
        [](double p, const std::pair<double,double>& x) { return p<x.second; });
      double dx = (p-(iter-1)->second)/(iter->second-(iter-1)->second);
      inverse_cdf[ibin] = (iter-1)->first + dx*(iter->first - (iter-1)->first);
    }

    iter = std::lower_bound(cdf.begin(), cdf.end(), 1.0,
      [](const std::pair<double,double>& x, double p) { return x.second<p; });
    inverse_cdf[nbins-1] = iter->first;

    return inverse_cdf;
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

template<typename VCLRealArch> class VCLRealRNG: public VCLRealArch
{
public:
  using typename VCLRealArch::architecture;
  using typename VCLRealArch::real_t;
  using typename VCLRealArch::real_vt;
  using typename VCLRealArch::int_vt;
  using typename VCLRealArch::uint_vt;
  using typename VCLRealArch::vec3_vt;

  VCLRealRNG(VCLRNG<architecture>* rng, bool adopt_rng = false):
    rng_(rng), adopt_rng_(adopt_rng)
  {
    // nothing to see here
  }

  VCLRealRNG(const std::string& created_by = "", const std::string& comment = "",
      uint64_t seed = 0, typename VCLRNG<architecture>::CoreType core_type = VCLRNG<architecture>::CoreType::NR3):
    rng_(new VCLRNG<architecture>(created_by, comment, seed, core_type)),
    adopt_rng_(true)
  {
    // nothing to see here
  }

  VCLRealRNG(uint64_t seed, const std::string& created_by = "", const std::string& comment = "",
      typename VCLRNG<architecture>::CoreType core_type = VCLRNG<architecture>::CoreType::NR3):
    rng_(new VCLRNG<architecture>(seed, created_by, comment, core_type)),
    adopt_rng_(true)
  {
    /* nothing to see here */
  }

  ~VCLRealRNG()
  {
    if(adopt_rng_)delete rng_;
  }

  int_vt uniform_int() { int_vt x; rng_->uniform_int(x); return x; }
  uint_vt uniform_uint() { uint_vt x; rng_->uniform_uint(x); return x; }

  real_vt uniform_gen(const real_vt& scale = 1.0, const real_vt& offset = 0.5) {
    real_vt x; rng_->uniform_real_gen(x, scale, offset); return x; }
  real_vt uniform_zc(const real_vt& scale = 1.0) {
    real_vt x; rng_->uniform_real_zc(x, scale); return x; }
  real_vt uniform(const real_vt& scale = 1.0) {
    real_vt x; rng_->uniform_real(x, scale); return x; }

  real_vt exponential() { real_vt x; rng_exponential_real(x); return x; }

  void normal_two_bm(real_vt& x1, real_vt& x2) {
    rng_->normal_two_real_bm(x1,x2); }
  real_vt normal() {
    real_vt x1, x2; rng_->normal_two_real_bm(x1,x2); return x1; }

  void uniform_on_unit_sphere(real_vt& x, real_vt& y, real_vt& z) {
    rng_->uniform_on_unit_sphere_real(x,y,z); }
  vec3_vt uniform_on_unit_sphere() {
    return rng_->uniform_on_unit_sphere_vec_real(); }

  real_vt from_inverse_cdf(const real_t* inverse_cdf, unsigned npoints) {
    return from_inverse_cdf_real(inverse_cdf, npoints); }
  real_vt from_inverse_cdf(const std::vector<real_t>& inverse_cdf) {
    return from_inverse_cdf_real(inverse_cdf); }

private:
  VCLRNG<architecture>* rng_;
  bool adopt_rng_;
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

  NR3_VCLRNGCore(
      const std::string& created_by = "", const std::string& comment = ""):
    VCLRNGCore<VCLArchitecture>(),
    seed_(RNG::nonzero_uint64_from_random_device()),
    sequence_seeds_(VCLRNG<VCLArchitecture>::uint64_from_seed(seed_))
  {
    init();
    this->write_provenance(created_by, comment);
  }

  NR3_VCLRNGCore(uint64_t seed,
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

  uint64_t calls() const { return calls_; }
  uint64_t seed() const { return seed_; }
  uint64_vt sequence_seeds() const { return sequence_seeds_; }

  uint64_vt uniform_uint64() override {
    calls_++;
    u_ = calin::util::vcl::mul_64(u_, C_NR3_U_MUL) + C_NR3_U_ADD;
    v_ ^= v_ >> C_NR3_V_SHIFT1;
    v_ ^= v_ << C_NR3_V_SHIFT2;
    v_ ^= v_ >> C_NR3_V_SHIFT3;
    // This is inefficient as the full 64bit multiply is not needed
    // w_ = C_NR3_W_MUL*(w_ & 0xFFFFFFFF) + (w_ >> C_NR3_W_SHIFT1);
    w_ = calin::util::vcl::mul_low32_packed64(w_, C_NR3_W_MUL) + (w_ >> C_NR3_W_SHIFT1);
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

  static ix_core_data_type* mutable_core_data(ix::math::rng::VCLRNGCoreData* proto) {
    return proto->mutable_nr3_vcl_core(); }
  static const ix_core_data_type& core_data(const ix::math::rng::VCLRNGCoreData& proto) {
    return proto.nr3_vcl_core(); }

#ifndef SWIG
  static void* operator new(size_t nbytes) {
    void* p = nullptr;
    if(::posix_memalign(&p, CALIN_NEW_ALIGN, nbytes)==0) {
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
