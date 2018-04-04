/*

   calin/math/rng.cpp -- Stephen Fegan -- 2015-11-19

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
  // nothing to see here
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
                                    bool restore_state)
{
  switch(proto.core_case()) {
    case ix::math::rng::RNGCoreData::kRanlux48Core:
      return new Ranlux48RNGCore(proto.ranlux48_core(), restore_state);
    case ix::math::rng::RNGCoreData::kMt19937Core:
      return new MT19937RNGCore(proto.mt19937_core(), restore_state);
    case ix::math::rng::RNGCoreData::kNr3Core:
      return new NR3RNGCore(proto.nr3_core(), restore_state);
    case ix::math::rng::RNGCoreData::kNr3SimdEmu4Core:
      return new NR3_EmulateSIMD_RNGCore<4>(proto.nr3_simd_emu4_core(), restore_state);
    case ix::math::rng::RNGCoreData::kNr3Avx2Core:
      return new NR3_AVX2_RNGCore(proto.nr3_avx2_core(), restore_state);
    default:
      LOG(ERROR) << "Unrecognised RNG type case : " << proto.core_case();
    case 0:
      return new NR3RNGCore();
  }
}

void RNGCore::write_provenance(
  const std::string& created_by, const std::string& comment)
{
  const ix::math::rng::RNGCoreData* proto = this->as_proto();
  calin::provenance::chronicle::register_rng_core(*proto, created_by, comment);
  delete proto;
}

RNG::RNG(CoreType core_type, const std::string& created_by): core_(nullptr), adopt_core_(true)
{
  switch(core_type) {
    case CoreType::NR3:
      core_ = new NR3RNGCore;
      break;
    case CoreType::RANLUX48:
      core_ = new Ranlux48RNGCore;
      break;
    case CoreType::MT19937:
      core_ = new MT19937RNGCore;
      break;
    default:
      throw std::invalid_argument("Unsupported RNG core type: " + std::to_string(int(core_type)));
  }
  const ix::math::rng::RNGData* proto = this->as_proto();
  calin::provenance::chronicle::register_calin_rng(*proto, created_by);
  delete proto;
}

RNG::RNG(uint64_t seed, CoreType core_type, const std::string& created_by):
  core_(nullptr), adopt_core_(true)
{
  switch(core_type) {
    case CoreType::NR3:
      core_ = new NR3RNGCore(seed);
      break;
    case CoreType::RANLUX48:
      core_ = new Ranlux48RNGCore(seed);
      break;
    case CoreType::MT19937:
      core_ = new MT19937RNGCore(seed);
      break;
    default:
      throw std::invalid_argument("Unsupported RNG core type: " + std::to_string(int(core_type)));
  }
  const ix::math::rng::RNGData* proto = this->as_proto();
  calin::provenance::chronicle::register_calin_rng(*proto, created_by);
  delete proto;
}

RNG::RNG(RNGCore* core, bool adopt_core):
    core_(core), adopt_core_(adopt_core)
{
  // nothing to see here
}

RNG::RNG(const ix::math::rng::RNGData& proto, bool restore_state, const std::string& created_by):
    core_(RNGCore::create_from_proto(proto.core(), restore_state)), adopt_core_(true)
{
  if(restore_state)
  {
    bm_hascached_ = proto.bm_hascached();
    bm_cachedval_ = proto.bm_cachedval();
  }
  const ix::math::rng::RNGData* proto2 = this->as_proto();
  calin::provenance::chronicle::register_calin_rng(*proto2, created_by);
  delete proto2;
}

RNG::~RNG()
{
  if(adopt_core_)delete core_;
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


NR3RNGCore::
NR3RNGCore(const ix::math::rng::NR3RNGCoreData& proto, bool restore_state):
    NR3RNGCore(proto.seed())
{
  if(restore_state and proto.state_saved())
  {
    calls_ = proto.calls();
    u_ = proto.u();
    v_ = proto.v();
    w_ = proto.w();
  }
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

Ranlux48RNGCore::
Ranlux48RNGCore(const ix::math::rng::Ranlux48RNGCoreData& proto,
                bool restore_state):
    RNGCore(), gen_seed_(proto.seed()), gen_(gen_seed_)
{
  if(restore_state and proto.state_saved())
  {
    std::istringstream state(proto.state());
    state >> gen_;
    gen_calls_ = proto.calls();
    dev_ = proto.dev();
    dev_blocks_ = proto.dev_blocks();
  }
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

MT19937RNGCore::
MT19937RNGCore(const ix::math::rng::STLRNGCoreData& proto,
               bool restore_state):
  RNGCore(), gen_seed_(proto.seed()), gen_(gen_seed_)
{
  if(restore_state and proto.state_saved())
  {
    std::istringstream state(proto.state());
    state >> gen_;
    gen_calls_ = proto.calls();
  }
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

#if defined(CALIN_HAS_NR3_AVX2_RNGCORE)
void NR3_AVX2_RNGCore::test_cpu() const
{
  const auto* sysinfo = calin::provenance::system_info::the_host_info();
  if(not sysinfo->cpu_has_avx2())
    throw std::runtime_error("NR3_AVX2_RNGCore: CPU does not support AVX2 instructions.");
}
#endif

NR3_AVX2_RNGCore::~NR3_AVX2_RNGCore()
{
  // nothing to see here
}

NR3_AVX2_RNGCore::NR3_AVX2_RNGCore(const ix::math::rng::NR3_SIMD_RNGCoreData& proto,
    bool restore_state):
  NR3_AVX2_RNGCore(proto.seed())
{
#if defined(CALIN_HAS_NR3_AVX2_RNGCORE)
  test_cpu();
  if(proto.vec_stream_seed_size() != 0) {
    // We allow either zero of four seeds.. zero means we do "seed only" reinit
    if(proto.vec_stream_seed_size() != 4)
      throw std::runtime_error("NR3_AVX2_RNGCore: saved seed vectors must have 4 elements");
    init(proto.vec_stream_seed(3), proto.vec_stream_seed(2),
         proto.vec_stream_seed(1), proto.vec_stream_seed(0));
    if(restore_state and proto.state_saved())
    {
      if(proto.vec_u_size() != 4 or proto.vec_v_size() != 4 or proto.vec_w_size() != 4)
        throw std::runtime_error("NR3_AVX2_RNGCore: saved state vectors must have 4 elements");
      if(proto.dev_size() > 4)
        throw std::runtime_error("NR3_AVX2_RNGCore: saved deviate vector must have at most 4 elements");
      calls_ = proto.calls();
      vec_u_ = _mm256_set_epi64x(proto.vec_u(3), proto.vec_u(2), proto.vec_u(1), proto.vec_u(0));
      vec_v_ = _mm256_set_epi64x(proto.vec_v(3), proto.vec_v(2), proto.vec_v(1), proto.vec_v(0));
      vec_w_ = _mm256_set_epi64x(proto.vec_w(3), proto.vec_w(2), proto.vec_w(1), proto.vec_w(0));
      ndev_ = proto.dev_size();
      switch(ndev_) {
      case 0:
      default:
        vec_dev_ = _mm256_setzero_si256();
        break;
      case 1:
        vec_dev_ = _mm256_set_epi64x(0U, 0U, 0U, proto.dev(0));
        break;
      case 2:
        vec_dev_ = _mm256_set_epi64x(0U, 0U, proto.dev(1), proto.dev(0));
        break;
      case 3:
        vec_dev_ = _mm256_set_epi64x(0U, proto.dev(2), proto.dev(1), proto.dev(0));
        break;
      case 4:
        vec_dev_ = _mm256_set_epi64x(proto.dev(3), proto.dev(2), proto.dev(1), proto.dev(0));
        break;
      }
    }
  }
#else
  throw std::runtime_error("NR3_AVX2_RNGCore: AVX2 not present at compile time.");
#endif
}

void NR3_AVX2_RNGCore::save_to_proto(ix::math::rng::RNGCoreData* proto) const
{
  auto* data = proto->mutable_nr3_avx2_core();
  save_to_proto(data);
}

void NR3_AVX2_RNGCore::save_to_proto(ix::math::rng::NR3_SIMD_RNGCoreData* data) const
{
  data->set_seed(seed_);
#if defined(CALIN_HAS_NR3_AVX2_RNGCORE)
  data->set_calls(calls_);
  data->set_state_saved(true);
  data->add_vec_stream_seed(stream_seed3_);
  data->add_vec_stream_seed(stream_seed2_);
  data->add_vec_stream_seed(stream_seed1_);
  data->add_vec_stream_seed(stream_seed0_);

  data->add_vec_u(_mm256_extract_epi64(vec_u_,0));
  data->add_vec_u(_mm256_extract_epi64(vec_u_,1));
  data->add_vec_u(_mm256_extract_epi64(vec_u_,2));
  data->add_vec_u(_mm256_extract_epi64(vec_u_,3));

  data->add_vec_v(_mm256_extract_epi64(vec_v_,0));
  data->add_vec_v(_mm256_extract_epi64(vec_v_,1));
  data->add_vec_v(_mm256_extract_epi64(vec_v_,2));
  data->add_vec_v(_mm256_extract_epi64(vec_v_,3));

  data->add_vec_w(_mm256_extract_epi64(vec_w_,0));
  data->add_vec_w(_mm256_extract_epi64(vec_w_,1));
  data->add_vec_w(_mm256_extract_epi64(vec_w_,2));
  data->add_vec_w(_mm256_extract_epi64(vec_w_,3));

  if(ndev_>0U)data->add_dev(_mm256_extract_epi64(vec_dev_,0));
  if(ndev_>1U)data->add_dev(_mm256_extract_epi64(vec_dev_,1));
  if(ndev_>2U)data->add_dev(_mm256_extract_epi64(vec_dev_,2));
  if(ndev_>3U)data->add_dev(_mm256_extract_epi64(vec_dev_,3));
#else
  throw std::runtime_error("NR3_AVX2_RNGCore: AVX2 not present at compile time.");
#endif
}

void NR3_AVX2_RNGCore::bulk_uniform_uint64(void* buffer, std::size_t nbytes)
{
#if defined(CALIN_HAS_NR3_AVX2_RNGCORE)
  __m256i* u256_ptr = reinterpret_cast<__m256i*>(buffer);
  std::size_t n256 = nbytes >> 5;
  for(std::size_t i256=0; i256<n256; i256++) {
    __m256i u256 = uniform_m256i();
    _mm256_storeu_si256(u256_ptr+i256, u256);
  }
  u256_ptr += n256;
  nbytes -= n256*32;
  if(nbytes)
  {
    __m256i u256 = uniform_m256i();
    std::copy(reinterpret_cast<uint8_t*>(u256_ptr),
      reinterpret_cast<uint8_t*>(u256_ptr)+nbytes, reinterpret_cast<uint8_t*>(&u256));
  }
#else
  throw std::runtime_error("NR3_AVX2_RNGCore: AVX2 not present at compile time.");
#endif
}

void NR3_AVX2_RNGCore::
bulk_uniform_uint64_with_mask(void* buffer, std::size_t nbytes, uint64_t mask)
{
#if defined(CALIN_HAS_NR3_AVX2_RNGCORE)
  const __m256i mask256 = _mm256_set1_epi64x(mask);
  __m256i* u256_ptr = reinterpret_cast<__m256i*>(buffer);
  std::size_t n256 = nbytes >> 5;
  for(std::size_t i256=0; i256<n256; i256++) {
    __m256i u256 = _mm256_and_si256(uniform_m256i(), mask256);
    _mm256_storeu_si256(u256_ptr+i256, u256);
  }
  u256_ptr += n256;
  nbytes -= n256*32;
  if(nbytes)
  {
    __m256i u256 = _mm256_and_si256(uniform_m256i(), mask256);
    std::copy(reinterpret_cast<uint8_t*>(u256_ptr),
      reinterpret_cast<uint8_t*>(u256_ptr)+nbytes, reinterpret_cast<uint8_t*>(&u256));
  }
#else
  throw std::runtime_error("NR3_AVX2_RNGCore: AVX2 not present at compile time.");
#endif
}

#if 0

// *****************************************************************************
// *****************************************************************************
// *****************************************************************************
// *****************************************************************************
// *****************************************************************************
// *****************************************************************************

// LEGACY CODE THAT MIGHT BE USEFUL SOME TIME

// *****************************************************************************
// *****************************************************************************
// *****************************************************************************
// *****************************************************************************
// *****************************************************************************
// *****************************************************************************

namespace RNGCore
{

  /**
   *  Generator from NR3 which is recommended for everyday use. Faster than
   *  NR3Ran with period of 1.8E19.
   */

  class NR3Ranq1
  {
  public:
    typedef EmptyOptions Options;
    static Options defaultOptions() { return Options(); }

    NR3Ranq1(uint64_t seed, const Options& opt = Options()):
      m_v(UINT64_C(4101842887655102017))
    {
      m_v ^= seed;
      m_v = UInt64();
    }

    void Burn() { UInt64(); }

    uint64_t UInt64()
    {
      m_v ^= m_v>>21;
      m_v ^= m_v<<35;
      m_v ^= m_v>>4;
      return m_v*UINT64_C(2685821657736338717);
    }

    uint32_t UInt32() { return uint32_t(UInt64()); }
    double Double() { return 5.42101086242752217E-20 * double(UInt64()); }

    static std::string getCoreName() { return "NR3Ranq1"; }

    void saveCoreState(std::ostream& stream) const { stream << m_v << '\n'; }
    bool loadCoreState(std::istream& stream) { return stream >> m_v; }

  private:
    uint64_t m_v;
  };

  /**
   *  Generator from NR3 which is described as a backup for times when
   *  Ranq1 has too short a period and Ran is too slow. Period of 8.5E37.
   */

  class NR3Ranq2
  {
  public:
    typedef EmptyOptions Options;
    static Options defaultOptions() { return Options(); }

    NR3Ranq2(uint64_t seed, const Options& opt = Options()):
      m_v(UINT64_C(4101842887655102017)), m_w(UINT64_C(1))
    {
      m_v ^= seed;
      m_v = UInt64();
      m_w = UInt64();
    }

    void Burn() { UInt64(); }

    uint64_t UInt64()
    {
      m_v ^= m_v>>17;
      m_v ^= m_v<<31;
      m_v ^= m_v>>8;
      m_w = 4294957665U*(m_w & 0xFFFFFFFF) + (m_w >> 32);
      return m_v^m_w;
    }

    uint32_t UInt32() { return uint32_t(UInt64()); }
    double Double() { return 5.42101086242752217E-20 * double(UInt64()); }

    static std::string getCoreName() { return "NR3Ranq2"; }

    void saveCoreState(std::ostream& stream) const
    {
      stream << m_v << '\n' << m_w << '\n';
    }

    bool loadCoreState(std::istream& stream)
    {
      return stream >> m_v >> m_w;
    }

  private:
    uint64_t m_v;
    uint64_t m_w;
  };

  /**
   *  Generator from NR3 implementing Knuth's subtractive generator (aka
   *  Fibonacci generator), a fast double precision generator. NR say it
   *  is a "good but not great generator with speed as it principal
   *  recommendation".
   */

  class NR3Ranfib
  {
  public:
    typedef EmptyOptions Options;
    static Options defaultOptions() { return Options(); }

    NR3Ranfib(uint64_t seed, const Options& opt = Options()):
      m_dtab(), m_dd(), m_inext(0), m_inextp(31)
    {
      NR3Ran rng(seed);
      for(unsigned k=0;k<55;k++)m_dtab[k]=rng.Double();
    }

    void Burn() { Double(); }

    double Double()
    {
      if(++m_inext==55)m_inext=0;
      if(++m_inextp==55)m_inextp=0;
      m_dd = m_dtab[m_inext]-m_dtab[m_inextp];
      if(m_dd<0) m_dd+=1.0;
      return (m_dtab[m_inext]=m_dd);
    }

    static std::string getCoreName() { return "NR3Ranfib"; }

    void saveCoreState(std::ostream& stream) const
    {
      std::copy(m_dtab, m_dtab+55, std::ostream_iterator<double>(stream,"\n"));
      stream << m_dd << '\n' << m_inext << '\n' << m_inextp << '\n';
    }

    bool loadCoreState(std::istream& stream)
    {
      double dtab[55];
      double dd;
      int inext;
      int inextp;
      for(int iel=0;iel<55;iel++)if(!(stream >> dtab[iel]))return false;
      if(!(stream >> dd >> inext >> inextp))return false;
      for(int iel=0;iel<55;iel++)m_dtab[iel] = dtab[iel];
      m_dd = dd;
      m_inext = inext;
      m_inextp = inextp;
      return true;
    }

  private:
    double m_dtab[55];
    double m_dd;
    int m_inext;
    int m_inextp;
  };

  /**
   *  Generator from NR2, no longer recommended by NR3. Only supports
   *  double type. Do not try to use UInt32 and UInt64.
   */

  class NR2Ran2
  {
  public:
    typedef EmptyOptions Options;
    static Options defaultOptions() { return Options(); }

    NR2Ran2(uint64_t seed, const Options& opt = Options());

    void Burn() { ran2(); }
    double Double() { return ran2(); }

    static std::string getCoreName() { return "NR2Ran2"; }

    void saveCoreState(std::ostream& stream) const;
    bool loadCoreState(std::istream& stream);

  private:
    double ran2();
    void ran2_init(int32_t);

    static const unsigned RANDOMNUMBERS_NTAB = 32;

    static const int32_t RAN2_A1=40014;
    static const int32_t RAN2_A2=40692;
    static const int32_t RAN2_Q1=53668;
    static const int32_t RAN2_Q2=52774;
    static const int32_t RAN2_R1=12211;
    static const int32_t RAN2_R2=3791;
    static const int32_t RAN2_M1=2147483563;
    static const int32_t RAN2_M2=2147483399;
    static const double RNMX;

    int32_t ran2_idum1;
    int32_t ran2_idum2;
    int32_t ran2_iy;
    int32_t ran2_iv[RANDOMNUMBERS_NTAB];
  };
}

// ============================================================================
// ============================================================================
//
// FUNCTIONS FOR RANDOMNUMBERSBASE CLASS
//
// ============================================================================
// ============================================================================


// ============================================================================
// ============================================================================
//
// FUNCTIONS FOR RANDOMNUMBERS CLASS
//
// ============================================================================
// ============================================================================

#endif
