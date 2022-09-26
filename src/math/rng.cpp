/*

   calin/math/rng.cpp -- Stephen Fegan -- 2015-11-19

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

#include <sstream>
#include <stdexcept>
#include <cassert>
#include <random>

#include <math/rng.hpp>
#include <provenance/chronicle.hpp>
#include <provenance/system_info.hpp>
#include <util/log.hpp>
#include <util/memory.hpp>

using namespace calin::math::rng;
using namespace calin::util::log;

RNGCore::~RNGCore()
{
  if(chronicle_record_)
    calin::provenance::chronicle::register_rng_close(chronicle_record_);
}

void RNGCore::bulk_uniform_uint64(void* buffer, std::size_t nbytes)
{
  uint64_t* u64_ptr = reinterpret_cast<uint64_t*>(buffer);
  std::size_t n64 = nbytes >> 3;
  for(std::size_t i64=0; i64<n64; i64++) {
    u64_ptr[i64] = this->uniform_uint64();
  }
  u64_ptr += n64;
  nbytes -= n64*8;
  if(nbytes)
  {
    uint64_t u64 = this->uniform_uint64();
    std::copy(reinterpret_cast<uint8_t*>(u64_ptr),
      reinterpret_cast<uint8_t*>(u64_ptr)+nbytes, reinterpret_cast<uint8_t*>(&u64));
  }
}

void RNGCore::bulk_uniform_uint64_with_mask(void* buffer, std::size_t nbytes, uint64_t mask)
{
  uint64_t* u64_ptr = reinterpret_cast<uint64_t*>(buffer);
  std::size_t n64 = nbytes >> 3;
  for(std::size_t i64=0; i64<n64; i64++) {
    u64_ptr[i64] = this->uniform_uint64() & mask;
  }
  u64_ptr += n64;
  nbytes -= n64*8;
  if(nbytes)
  {
    uint64_t u64 = this->uniform_uint64() & mask;
    std::copy(reinterpret_cast<uint8_t*>(u64_ptr),
      reinterpret_cast<uint8_t*>(u64_ptr)+nbytes, reinterpret_cast<uint8_t*>(&u64));
  }
}

RNGCore* RNGCore::create_from_proto(const ix::math::rng::RNGCoreData& proto,
  bool restore_state, const std::string& created_by, const std::string& comment)
{
  switch(proto.core_case()) {
    case ix::math::rng::RNGCoreData::kRanlux48Core:
      return new Ranlux48RNGCore(proto.ranlux48_core(), restore_state, created_by, comment);
    case ix::math::rng::RNGCoreData::kMt19937Core:
      return new MT19937RNGCore(proto.mt19937_core(), restore_state, created_by, comment);
    case ix::math::rng::RNGCoreData::kNr3Core:
      return new NR3RNGCore(proto.nr3_core(), restore_state, created_by, comment);
    case ix::math::rng::RNGCoreData::kNr3SimdEmu4Core:
      return new NR3_EmulateSIMD_RNGCore<4>(proto.nr3_simd_emu4_core(), restore_state);
    default:
      LOG(ERROR) << "Unrecognised RNG type case : " << proto.core_case();
    case 0:
      return new NR3RNGCore();
  }
}

std::vector<uint64_t> RNGCore::vec_uniform_uint64(std::size_t nelements)
{
  std::vector<uint64_t> vec(nelements);
  bulk_uniform_uint64(vec.data(), nelements*sizeof(uint64_t));
  return vec;
}

std::vector<uint64_t> RNGCore::vec_uniform_uint64_with_mask(
  std::size_t nelements, uint64_t mask)
{
  std::vector<uint64_t> vec(nelements);
  bulk_uniform_uint64_with_mask(vec.data(), nelements*sizeof(uint64_t), mask);
  return vec;
}

void RNGCore::write_provenance(
  const std::string& created_by, const std::string& comment)
{
  const ix::math::rng::RNGCoreData* proto = this->as_proto();
  chronicle_record_ =
    calin::provenance::chronicle::register_rng_core_open(*proto, created_by, comment);
  delete proto;
}

RNG::RNG(CoreType core_type, const std::string& created_by,
  const std::string& comment): core_(nullptr), adopt_core_(true)
{
  switch(core_type) {
    case CoreType::NR3:
      core_ = new NR3RNGCore(created_by, comment);
      break;
    case CoreType::RANLUX48:
      core_ = new Ranlux48RNGCore(created_by, comment);
      break;
    case CoreType::MT19937:
      core_ = new MT19937RNGCore(created_by, comment);
      break;
    default:
      throw std::invalid_argument("Unsupported RNG core type: " + std::to_string(int(core_type)));
  }
  const ix::math::rng::RNGData* proto = this->as_proto();
  chronicle_record_ =
    calin::provenance::chronicle::register_calin_rng_open(*proto, created_by, comment);
  delete proto;
}

RNG::RNG(uint64_t seed, CoreType core_type, const std::string& created_by,
    const std::string& comment):
  core_(nullptr), adopt_core_(true)
{
  switch(core_type) {
    case CoreType::NR3:
      core_ = new NR3RNGCore(seed, created_by, comment);
      break;
    case CoreType::RANLUX48:
      core_ = new Ranlux48RNGCore(seed, created_by, comment);
      break;
    case CoreType::MT19937:
      core_ = new MT19937RNGCore(seed, created_by, comment);
      break;
    default:
      throw std::invalid_argument("Unsupported RNG core type: " + std::to_string(int(core_type)));
  }
  const ix::math::rng::RNGData* proto = this->as_proto();
  chronicle_record_ =
    calin::provenance::chronicle::register_calin_rng_open(*proto, created_by, comment);
  delete proto;
}

RNG::RNG(RNGCore* core, bool adopt_core,
    const std::string& created_by, const std::string& comment):
  core_(core), adopt_core_(adopt_core)
{
  const ix::math::rng::RNGData* proto = this->as_proto();
  chronicle_record_ =
    calin::provenance::chronicle::register_calin_rng_open(*proto, created_by, comment);
  delete proto;
}

RNG::RNG(const ix::math::rng::RNGData& proto, bool restore_state,
    const std::string& created_by, const std::string& comment):
  core_(RNGCore::create_from_proto(proto.core(), restore_state, created_by, comment)), adopt_core_(true)
{
  if(restore_state)
  {
    bm_hascached_ = proto.bm_hascached();
    bm_cachedval_ = proto.bm_cachedval();
  }
  const ix::math::rng::RNGData* proto2 = this->as_proto();
  chronicle_record_ =
    calin::provenance::chronicle::register_calin_rng_open(*proto2, created_by, comment);
  delete proto2;
}

RNG::~RNG()
{
  if(adopt_core_)delete core_;
  if(chronicle_record_)
    calin::provenance::chronicle::register_rng_close(chronicle_record_);
}

void RNG::save_to_proto(ix::math::rng::RNGData* proto) const
{
  core_->save_to_proto(proto->mutable_core());
  proto->set_bm_hascached(bm_hascached_);
  proto->set_bm_cachedval(bm_hascached_ ? bm_cachedval_ : 0.0);
}

uint64_t RNG::uint64_from_random_device()
{
  std::random_device gen;
  static_assert(sizeof(std::random_device::result_type)==sizeof(uint32_t),
                "std::random_device::result_type is not 32 bits");
  return uint64_t(gen())<<32 | uint64_t(gen());
}

uint32_t RNG::uint32_from_random_device()
{
  std::random_device gen;
  static_assert(sizeof(std::random_device::result_type)==sizeof(uint32_t),
                "std::random_device::result_type is not 32 bits");
  return gen();
}

double RNG::normal()
{
  if(bm_hascached_)
  {
    bm_hascached_ = false;
    return bm_cachedval_;
  }
  else
  {
    double v1;
    double v2;
    double rsq;
    do
  	{
  	  v1 = 2.0*uniform() - 1.0;
  	  v2 = 2.0*uniform() - 1.0;
  	  rsq = v1*v1 + v2*v2;
  	}while(rsq >= 1.0 || rsq == 0.0);
    const double fac = sqrt(-2.0*log(rsq)/rsq);
    bm_cachedval_ = v1*fac;
    bm_hascached_ = true;
    return v2*fac;
  }
}

void RNG::normal_two_bm(double &x, double& y)
{
  double v1;
  double v2;
  double rsq;
  do
	{
	  v1 = 2.0*uniform() - 1.0;
	  v2 = 2.0*uniform() - 1.0;
	  rsq = v1*v1 + v2*v2;
	}while(rsq >= 1.0 || rsq == 0.0);
  const double fac = sqrt(-2.0*log(rsq)/rsq);
  x = v1*fac;
  y = v2*fac;
}

namespace {
  inline double uint64_to_double(uint64_t u) {
    constexpr uint64_t MASK_HI = (1ULL<<52)-1;
    constexpr uint64_t EXP_HI = 1023ULL<<52;
    u &= MASK_HI;
    u |= EXP_HI;
    return *reinterpret_cast<double*>(&u) - 1.0;
  }
}

double RNG::ziggurat_normal()
{
  using namespace gaussian_ziggurat;
  constexpr uint64_t MASK_SIGN = 1ULL<<63;
  while(true) {
    uint64_t u0 = core_->uniform_uint64();
    uint64_t i = u0&0xFF;
    double x = xi[i+1]*uint64_to_double(u0>>8);
    if(x < xi[i]) {
      return (u0&MASK_SIGN) ? -x : x;
    } else if(i != 0xFF) {
      double fx = std::exp(-0.5*x*x);
      double y = uint64_to_double(core_->uniform_uint64());
      if(y*(fi[i]-fi[i+1]) < fx-fi[i+1]) {
        return (u0&MASK_SIGN) ? -x : x;
      }
    } else {
      double y;
      do
        {
          x = std::log(uint64_to_double(core_->uniform_uint64())) * r_inv;
          y = std::log(uint64_to_double(core_->uniform_uint64()));
        } while (-2 * y < x * x);
      return (u0&MASK_SIGN) ? x - r : r - x;
    }
  }
}

/**
 *  Returns a deviate distributed as a gamma distribution, i.e., a
 *  waiting time to the i'th event in a Poisson process of unit mean.
 *  pdf=beta^alpha * x^(alpha-1) * exp(-beta x) / Gamma(alpha)
 */
double RNG::gamma_by_alpha_and_beta(double alpha, double beta)
{
  double oalpha = alpha;
  if(alpha <= 0.0)throw(std::invalid_argument("RNG::gamma: negative alpha"));
  if(alpha <= 1.0)alpha += 1.0;
  const double a1 = alpha - 1.0/3.0;
  const double a2 = 1./std::sqrt(9.0*a1);
  double u;
  double v;
  double x,x2,x4;
  do
  {
    do
  	{
  	  x = normal();
  	  v = 1.0 + a2*x;
  	}while(v<=0.0);
    v = v*v*v;
    u = uniform();
    x2 = x*x;
    x4 = x2*x2;
  }while((u>1.0-0.331*x4)&&(std::log(u)>0.5*x2+a1*(1.0-v+std::log(v))));
  if(alpha == oalpha)return a1*v/beta;
  else
  {
    do u=uniform(); while(u==0);
    return std::pow(u,1.0/oalpha)*a1*v/beta;
  }
}

namespace {

double lfactorial(unsigned ix)
{
  static const double lfa[] = {
    0.00000000000000000, 0.00000000000000000, 0.69314718055994529,
    1.79175946922805496, 3.17805383034794575, 4.78749174278204581,
    6.57925121201010121, 8.52516136106541467, 10.60460290274525086,
    12.80182748008146909, 15.10441257307551588, 17.50230784587388655,
    19.98721449566188468, 22.55216385312342098, 25.19122118273868338,
    27.89927138384089389, 30.67186010608067548, 33.50507345013689076,
    36.39544520803305261, 39.33988418719949465 };
  if(ix<sizeof(lfa)/sizeof(*lfa))return lfa[ix];
  return std::lgamma(ix+1);
}

} // anonymous namespace

/**
 *  Returns random deviate drawn from a Poisson distribution of mean
 *  lambda, using CORE as a source of uniform random deviates.
 */
int RNG::poisson(double lambda)
{
  if(lambda<5.0)
  {
    if(lambda != poi_lambdaold_)
  	{
  	  poi_lambdaexp_ = std::exp(-lambda);
  	  poi_lambdaold_ = lambda;
  	}
    int k = 0;
#if 1
    double t = uniform();
    while(t>poi_lambdaexp_) { ++k; t*=uniform(); }
#else
    uint64_t ilexp = 18446744073709551616.0*poi_lambdaexp_;
    uint64_t t = uniform_uint64();
    while(t>ilexp) { ++k; t=(t>>32)*(uniform_uint64()>>32); }
#endif
    return k;
  }
  else
  {
    int k;
    if(lambda != poi_lambdaold_)
  	{
  	  poi_lambdasrt_ = std::sqrt(lambda);
  	  poi_lambdalog_ = std::log(lambda);
  	  poi_lambdaold_ = lambda;
  	}
    while(1)
  	{
  	  double u = 0.64*uniform();
  	  double v = -0.68 + 1.28*uniform();
  	  double u2;
  	  if(lambda>13.5)
	    {
	      double v2 = v*v;
	      if(v >= 0) { if(v2 > 6.5*u*(0.64-u)*(u+0.2))continue; }
	      else { if(v2 > 9.6*u*(0.66-u)*(u+0.07))continue; }

	      k = int(std::floor(poi_lambdasrt_*(v/u)+lambda+0.5));
	      if(k<0)continue;
	      u2 = u*u;

	      if(v >= 0){ if(v2 < 15.2*u2*(0.61-u)*(0.8-u))break; }
	      else { if(v2 < 6.76*u2*(0.62-u)*(1.4-u))break; }
	    }
  	  else
	    {
	      k = int(std::floor(poi_lambdasrt_*(v/u)+lambda+0.5));
	      if(k<0)continue;
	      u2 = u*u;
	    }

  	  double lfact = lfactorial(k);
  	  double p =
        poi_lambdasrt_*std::exp(-lambda + k*poi_lambdalog_ - lfact);
  	  if(u2 < p)break;
  	}
    return k;
  }
  assert(0);
  return 0;
}

/**
 *  Returns an integer value that is a random deviate drawn from a
 *  binomial distribution of n trials each of probability pp, using
 *  CORE as a source of uniform random deviates.
 */

int RNG::binomial(double pp, int n)
{
  double bnl;

  const double p = (pp<=0.5 ? pp : 1.0-pp);
  const double am = n*p;

  if(n < 25)               /* Use direct method */
  {
    int ibnl=0;
    for(int i=0;i<n;i++)if(uniform()<p)++ibnl;
    bnl = double(ibnl);
  }
  else if(am < 1.0)
  {
    double g = std::exp(-am);
    double t = 1.0;
    int j;
    for(j=0;j<=n;j++) { t *= uniform(); if(t<g)break; }
    bnl = (j<=n?j:n);
  }
  else                       /* Use rejection method */
  {
    if(n != bin_nold_)
  	{
  	  bin_en_    = n;
  	  bin_oldg_  = lfactorial(n);
  	  bin_nold_  = n;
  	}

    if (p != bin_pold_)
  	{
  	  bin_pc_    = 1.0-p;
  	  bin_plog_  = std::log(p);
  	  bin_pclog_ = std::log(bin_pc_);
  	  bin_pold_  = p;
  	}

    double sq=std::sqrt(2.0*am*bin_pc_);
    double t;
    double em;
    do
  	{
  	  double y;
  	  do
	    {
	      double angle = M_PI*uniform();
	      y  = std::tan(angle);
	      em = sq*y+am;
	    }while(em<0.0 || em>=(bin_en_+1));
  	  em = std::floor(em);
  	  t = 1.2*sq*(1.0+y*y)*std::exp(bin_oldg_
                                    - lfactorial(unsigned(em))
                                    - lfactorial(unsigned(bin_en_-em))
                                    + em*bin_plog_
                                    + (bin_en_-em)*bin_pclog_);
  	}while(uniform()>t);
    bnl = em;
  }
  if (p != pp)bnl=n-bnl;
  return int(bnl);
}

namespace {

/**
 *  Performs a linear interpolation at the coordinate x using two
 *  points along a vector of pairs.
 */
inline double interpolate(double x,
  const std::vector<std::pair<double,double>>::const_iterator& itr)
{
  return (itr-1)->second +
    (x-(itr-1)->first)*((itr)->second-(itr-1)->second)/
    ((itr)->first-(itr-1)->first);
}

} // anonympus namespace

/**
 *  Returns a continuous random deviate drawn from the inverse CDF
 *  provided in the argument.  Note that this function assumes that
 *  the inverse CDF has equidistant steps in cumulative probability.
 */
double RNG::inverse_cdf(const std::vector<std::pair<double,double>> &inv_cdf)
{
  double x = uniform();
  unsigned i = (unsigned)ceil(x*(double)(inv_cdf.size()-1));
  return interpolate(x,inv_cdf.begin()+i);
}

/**
 *  Overwrites the CDF provided in the argument with its inverse.  The
 *  second argument specifies the number of equidistant bins in
 *  cumulative probability that will be generated. The first entry
 *  should be zero, and the last should be one.
 */
void RNG::generate_inverse_cdf(std::vector<std::pair<double,double>> &cdf,
  unsigned nbins)
{
  if(nbins == 0)
    nbins = cdf.size();

  for(auto itr = cdf.begin(); itr != cdf.end(); ++itr)
    *itr = std::make_pair(itr->second,itr->first);

  std::remove_reference<decltype(cdf)>::type inv_cdf;
  auto itr = cdf.cbegin()+1;

  for(unsigned i = 0; i < nbins; i++)
  {
    double x;

    if(i==0)
      x = 0.0;
    else if(i == nbins-1)
      x = 1.0;
    else
      x = (double)i/(double)(nbins-1);

    while(x > itr->first && itr+1 != cdf.end())
      itr++;

    inv_cdf.emplace_back(x,interpolate(x,itr));
  }

  cdf = inv_cdf;
}


// *****************************************************************************
// *****************************************************************************
//
// CORES
//
// *****************************************************************************
// *****************************************************************************

NR3RNGCore::NR3RNGCore(uint64_t seed,
    const std::string& created_by, const std::string& comment):
  RNGCore(),
  seed_(seed>0 ? seed : RNG::nonzero_uint64_from_random_device()),
  u_(C_NR3_U_INIT), v_(C_NR3_V_INIT), w_(C_NR3_W_INIT)
{
  u_ = seed_^v_; uniform_uint64();
  v_ = u_; uniform_uint64();
  w_ = v_; uniform_uint64();
  calls_ = 0;
  write_provenance(created_by, comment);
}

NR3RNGCore::
NR3RNGCore(const ix::math::rng::NR3RNGCoreData& proto, bool restore_state,
    const std::string& created_by, const std::string& comment):
  RNGCore(),
  seed_(proto.seed()>0 ? proto.seed() : RNG::nonzero_uint64_from_random_device()),
  u_(C_NR3_U_INIT), v_(C_NR3_V_INIT), w_(C_NR3_W_INIT)
{
  if(proto.seed()>0 and restore_state and proto.state_saved())
  {
    u_ = proto.u();
    v_ = proto.v();
    w_ = proto.w();
    calls_ = proto.calls();
  } else {
    u_ = seed_^v_; uniform_uint64();
    v_ = u_; uniform_uint64();
    w_ = v_; uniform_uint64();
    calls_ = 0;
  }
  write_provenance(created_by, comment);
}

NR3RNGCore::~NR3RNGCore()
{
  // nothing to see here
}

void NR3RNGCore::save_to_proto(ix::math::rng::RNGCoreData* proto) const
{
  auto* data = proto->mutable_nr3_core();
  save_to_proto(data);
}

void NR3RNGCore::save_to_proto(ix::math::rng::NR3RNGCoreData* data) const
{
  data->set_seed(seed_);
  data->set_state_saved(true);
  data->set_calls(calls_);
  data->set_u(u_);
  data->set_v(v_);
  data->set_w(w_);
}

Ranlux48RNGCore::Ranlux48RNGCore(uint64_t seed,
    const std::string& created_by, const std::string& comment):
  RNGCore(),
  gen_seed_(seed>0 ? seed : RNG::nonzero_uint64_from_random_device()),
  gen_(gen_seed_)
{
  write_provenance(created_by, comment);
}

Ranlux48RNGCore::
Ranlux48RNGCore(const ix::math::rng::Ranlux48RNGCoreData& proto, bool restore_state,
    const std::string& created_by, const std::string& comment):
  RNGCore(),
  gen_seed_(proto.seed()>0 ? proto.seed() : RNG::nonzero_uint64_from_random_device()),
  gen_(gen_seed_)
{
  if(proto.seed()>0 and restore_state and proto.state_saved())
  {
    std::istringstream state(proto.state());
    state >> gen_;
    gen_calls_ = proto.calls();
    dev_ = proto.dev();
    dev_blocks_ = proto.dev_blocks();
  }
  write_provenance(created_by, comment);
}

Ranlux48RNGCore::~Ranlux48RNGCore()
{
  // nothing to see here
}

void Ranlux48RNGCore::save_to_proto(ix::math::rng::RNGCoreData* proto) const
{
  auto* data = proto->mutable_ranlux48_core();
  save_to_proto(data);
}

void Ranlux48RNGCore::save_to_proto(ix::math::rng::Ranlux48RNGCoreData* data) const
{
  data->set_seed(gen_seed_);
  data->set_calls(gen_calls_);
  data->set_state_saved(true);
  std::ostringstream state;
  state << gen_;
  data->set_state(state.str());
  data->set_dev(dev_blocks_>0?dev_:0ULL);
  data->set_dev_blocks(dev_blocks_);
}

MT19937RNGCore::MT19937RNGCore(uint64_t seed,
    const std::string& created_by, const std::string& comment):
  RNGCore(),
  gen_seed_(seed>0 ? seed : RNG::nonzero_uint64_from_random_device()),
  gen_(gen_seed_)
{
  write_provenance(created_by, comment);
}

MT19937RNGCore::
MT19937RNGCore(const ix::math::rng::STLRNGCoreData& proto, bool restore_state,
    const std::string& created_by, const std::string& comment):
  RNGCore(),
  gen_seed_(proto.seed()>0 ? proto.seed() : RNG::nonzero_uint64_from_random_device()),
  gen_(gen_seed_)
{
  if(proto.seed()>0 and restore_state and proto.state_saved())
  {
    std::istringstream state(proto.state());
    state >> gen_;
    gen_calls_ = proto.calls();
  }
  write_provenance(created_by, comment);
}

MT19937RNGCore::~MT19937RNGCore()
{
  // nothing to see here
}

void MT19937RNGCore::save_to_proto(ix::math::rng::RNGCoreData* proto) const
{
  auto* data = proto->mutable_mt19937_core();
  save_to_proto(data);
}

void MT19937RNGCore::save_to_proto(ix::math::rng::STLRNGCoreData* data) const
{
  data->set_seed(gen_seed_);
  data->set_calls(gen_calls_);
  data->set_state_saved(true);
  std::ostringstream state;
  state << gen_;
  data->set_state(state.str());
}
