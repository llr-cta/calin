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

#include <stdexcept>
#include <cassert>

#include <math/rng.hpp>

using namespace calin::math::rng;

RNGCore::~RNGCore()
{
  // nothing to see here
}

RNG::RNG(uint64_t seed)
{
  
}

RNG::RNG(RNGCore* core, bool adopt_core):
    core_(core), adopt_core_(adopt_core)
{
  // nothing to see here
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
double RNG::gamma(double alpha, double beta)
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
      double t = uniform();
      while(t>poi_lambdaexp_) { ++k; t*=uniform(); }
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
      for(int i=0;i<=n;i++)if(uniform()<p)++ibnl;
      bnl = int(ibnl);
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

#if 0
namespace RNGCore
{
 
  /**
   *  Read random bytes from system random source (usually "/dev/random"). 
   *  This can supply high quality (non-pseudo) random numbers but is 
   *  usually VERY slow so it is advised that you only use this generator 
   *  to produce good quality seed numbers.
   */

  class DevRandom
  {
  public:
    typedef std::string Options;
    static Options defaultOptions() { return "/dev/random"; }
    
    DevRandom(const std::string& filename = defaultOptions()):
      m_file_name(filename), m_fp(fopen(filename.c_str(),"r")){ assert(m_fp); }
    DevRandom(uint64_t unused_seed, 
	      const std::string& filename = defaultOptions()):
      m_file_name(filename), m_fp(fopen(filename.c_str(),"r")){ assert(m_fp); }
    ~DevRandom() { fclose(m_fp); }

    void Burn() { }
    uint64_t UInt64() { uint64_t x; assert(fread(&x,8,1,m_fp)==1); return x; }
    uint32_t UInt32() { uint32_t x; assert(fread(&x,4,1,m_fp)==1); return x; }
    double Double() { return 5.42101086242752217E-20 * double(UInt64()); }

    static std::string getCoreName() { return "DevRandom"; }

    void saveCoreState(std::ostream& stream) const
    {
      stream << m_file_name << '\n';
    }

    bool loadCoreState(std::istream& stream)
    {
      std::string file_name;
      if(!(stream >> file_name))return false;
      if(file_name != m_file_name)
	{
	  fclose(m_fp);
	  m_file_name = file_name;
	  m_fp = fopen(file_name.c_str(),"r");
	  assert(m_fp);
	}
      return true;
    }

    static uint64_t oneUInt64() { DevRandom rng; return rng.UInt64(); }
    static uint64_t oneUInt32() { DevRandom rng; return rng.UInt32(); }
    static double oneDouble() { DevRandom rng; return rng.Double(); }
  private:
    std::string m_file_name;
    FILE*       m_fp;
  };  

  /**
   *  Highest quality random number generator from Numerical Recipes 3
   *  combines two unrelated generators to make generator with period
   *  of ~3.1E57.
   */

  class EmptyOptions {};

  class NR3Ran
  {
  public:
    typedef EmptyOptions Options;
    static Options defaultOptions() { return Options(); }

    NR3Ran(uint64_t seed, const Options& opt = Options()):
      m_u(UINT64_C(0)), m_v(UINT64_C(4101842887655102017)), m_w(UINT64_C(1))
    {
      m_u = seed^m_v; UInt64();
      m_v = m_u; UInt64();
      m_w = m_v; UInt64();
    }

    void Burn() { UInt64(); }

    uint64_t UInt64()
    {
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

    uint32_t UInt32() { return uint32_t(UInt64()); }
    double Double() { return 5.42101086242752217E-20 * double(UInt64()); }

    static std::string getCoreName() { return "NR3Ran"; }

    void saveCoreState(std::ostream& stream) const
    {
      stream << m_u << '\n' << m_v << '\n' << m_w << '\n';
    }

    bool loadCoreState(std::istream& stream) 
    {
      return stream >> m_u >> m_v >> m_w;
    }

  private:
    uint64_t m_u;
    uint64_t m_v;
    uint64_t m_w;
  };

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

  /**
   *  RanLux 3.2 generator. A particle physics favourite. 
   *  Copyright (C) 2005 Martin Luescher.
   *  http://luscher.web.cern.ch/luscher/ranlux/index.html
   *
   *  Integer option, which must be 1 (default) or 2, specifies
   *  "luxury" value, describing degree of residual correlation
   *  between subsequent deviates. The default should be sufficient
   *  for all applications.
   */

  class RanLuxV32
  {
  public:
    typedef int Options;
    static Options defaultOptions() { return 1; }

    RanLuxV32(uint64_t seed, const Options& opt = Options());

    void Burn() { double x; ranlxd(&x,1); }
    double Double() { double x; ranlxd(&x,1); return x; }

    static std::string getCoreName() { return "RanLuxV32"; }

    void saveCoreState(std::ostream& stream) const;
    bool loadCoreState(std::istream& stream);

  private:
    static const int BASE = 0x1000000;
    static const int MASK = 0xffffff;

    struct vec_t
    {
      int c1;
      int c2;
      int c3;
      int c4;
    };

    struct dble_vec_t
    {
      vec_t c1;
      vec_t c2;
    };
    
    int pr;
    int prm;
    int ir;
    int jr;
    int is;
    int is_old;
    int next[96];
    double one_bit;
    vec_t carry;

    union
    {
      dble_vec_t vec[12];
      int num[96];
    } x;

    void ranlxd(double r[],int n);
    void rlxd_init(int level,int seed);
    int rlxd_size(void) const;
    void rlxd_get(int state[]) const;
    void rlxd_reset(int state[]);

    void error(int no) const;
    void update();
    void define_constants();
  };

}

// ============================================================================
// ============================================================================
//
// FUNCTIONS FOR RANDOMNUMBERSBASE CLASS
//
// ============================================================================
// ============================================================================

/**
 *  Performs a linear interpolation at the coordinate x using two
 *  points along a vector of pairs.
 */ 
inline double RandomNumbersBase::
interpolate(double x, const std::vector<Pair>::const_iterator& itr)
{
  return (itr-1)->second + 
    (x-(itr-1)->first)*((itr)->second-(itr-1)->second)/
    ((itr)->first-(itr-1)->first);
}

// ============================================================================
// ============================================================================
//
// FUNCTIONS FOR RANDOMNUMBERS CLASS
//
// ============================================================================
// ============================================================================

// ----------------------------------------------------------------------------
// Constructors, destructors, and functions to load and save state
// ----------------------------------------------------------------------------

/** 
 *  Construct class from prescribed seed. Input and output of RNG state
 *  is disabled, so state will not be stored on destruction of the object.
 */
template<typename CORE> RandomNumbers_TNG<CORE>::
RandomNumbers_TNG(uint64_t seed, const CoreOptions& core_options): 
  RandomNumbersBase(), 
  CORE(seed>0?seed:RNGCore::DevRandom::oneUInt64(), core_options),
  m_bm_hascached(false), m_bm_cachedval(),
  m_poi_lambdaold(-1), m_poi_lambdaexp(), m_poi_lambdasrt(), m_poi_lambdalog(),
  m_bin_nold(-1), m_bin_pold(-1), m_bin_pc(), m_bin_plog(), 
  m_bin_pclog(), m_bin_en(), m_bin_oldg(),
  m_file_name()
{
  // nothing to see here
}

/** 
 *  Construct class from a stored state file. If state file does not exist 
 *  or is incompatible with this class then the generator is constructed 
 *  with a random seed and a state file constructed. The state file is locked
 *  for the duration of the existance of the generator.
 */
template<typename CORE> RandomNumbers_TNG<CORE>::
RandomNumbers_TNG(const std::string& state_filename, 
		  const CoreOptions& core_options):
  RandomNumbersBase(), 
  CORE(0/*RNGCore::DevRandom::oneUInt64()*/, core_options),
  m_bm_hascached(false), m_bm_cachedval(),
  m_poi_lambdaold(-1), m_poi_lambdaexp(), m_poi_lambdasrt(), m_poi_lambdalog(),
  m_bin_nold(-1), m_bin_pold(-1), m_bin_pc(), m_bin_plog(), 
  m_bin_pclog(), m_bin_en(), m_bin_oldg(),
  m_file_name(state_filename)
{
  lockState(m_file_name);                     // lock the lockfile
  if(!readState())                            // initialize seeds
    {
      uint64_t seed = RNGCore::DevRandom::oneUInt64();
      std::cerr 
	<< "RandomNumbers_TNG::readState(): could not read RNG state from file"
	<< std::endl
	<< "RandomNumbers_TNG::readState(): " << m_file_name << std::endl
	<< "RandomNumbers_TNG::readState(): initializing state from seed "
	<< seed << std::endl;
      *(CORE*)this = CORE(seed,core_options);
    }
  CORE::Burn();                               // generate a random number to
  writeState();                               // prevent repetitive crashes
}

#if 0
/** 
 *  Construct class from a stored state file. If state file does not exist 
 *  or is incompatible with this class then the generator is constructed 
 *  with a random seed and a state file constructed. The state file is locked
 *  for the duration of the existance of the generator.
 */
template<typename CORE> RandomNumbers_TNG<CORE>::
RandomNumbers_TNG(const char* state_filename, const CoreOptions& core_options):
  RandomNumbersBase(), 
  CORE(0/*RNGCore::DevRandom::oneUInt64()*/, core_options),
  m_bm_hascached(false), m_bm_cachedval(),
  m_poi_lambdaold(-1), m_poi_lambdaexp(), m_poi_lambdasrt(), m_poi_lambdalog(),
  m_bin_nold(-1), m_bin_pold(-1), m_bin_pc(), m_bin_plog(), 
  m_bin_pclog(), m_bin_en(), m_bin_oldg(),
  m_file_name(state_filename)
{
  assert(state_filename != 0);
  lockState(m_file_name);                     // lock the lockfile
  readState();                                // initialize seeds
  CORE::Burn();                               // generate a random number to
  writeState();                               // prevent repetitive crashes
}
#endif

/**
 *  Write the RNG state and unlock the state file.
 */

template<typename CORE> RandomNumbers_TNG<CORE>::~RandomNumbers_TNG()
{
  if(!m_file_name.empty())
    {
      writeState();                           // save seeds
      unlockState();                          // unlock the lockfile
    }
}

/**
 *  Return unique identifier string describing this class
 */
template<typename CORE> std::string RandomNumbers_TNG<CORE>::idTag()
{
  return 
    std::string("RandomNumbers_TNG<")+CORE::getCoreName()+std::string(">");
}

/**
 *  Load the generator state from a file
 */
template<typename CORE> bool RandomNumbers_TNG<CORE>::readState()
{
  std::ifstream fstream(m_file_name.c_str());
  if(fstream)
    {
      // Make sure state file is compatible with this class
      std::string idtag;
      std::getline(fstream,idtag);
      if((idtag == idTag())&&(CORE::loadCoreState(fstream)))
	{
	  fstream >> m_bm_hascached
		  >> m_bm_cachedval;
	  return true;
	}
    }

  return false;
}

/**
 *  Save the generator state to a file
 */
template<typename CORE> void RandomNumbers_TNG<CORE>::writeState()
{
  std::ofstream fstream(m_file_name.c_str());
  if(!fstream)
    {
      std::ostringstream stream;
      stream << "RandomNumbers_TNG::writeState(): could not open " 
	     << m_file_name << " for writing" << '\n'
	     << "RandomNumbers_TNG::writeState(): "
	     << strerror(errno) << '\n';
      throw stream.str();
    }

  fstream << std::scientific
	  << std::setprecision(std::numeric_limits<long double>::digits10)
	  << idTag() << '\n';
  CORE::saveCoreState(fstream);
  fstream << m_bm_hascached << '\n'
	  << m_bm_cachedval << '\n';
}

/**
 *  Return default state file name
 */
template<typename CORE> std::string RandomNumbers_TNG<CORE>::defaultFilename()
{
  std::ostringstream stream;
  const char* tmpdir = getenv("TMPDIR");
  if(tmpdir)
    {
      stream << tmpdir;
      if((*tmpdir!='\0')&&(tmpdir[strlen(tmpdir)-1]!='/'))stream << '/';
    }
  else
    {
      stream << "/tmp/";
    } 
  stream << "rng_" << CORE::getCoreName() << "_state_uid" << getuid()
	 << ".dat";
  return stream.str();
}

// ----------------------------------------------------------------------------
// Non-uniform distributions
// ----------------------------------------------------------------------------




/**
 *  Returns a continuous random deviate drawn from the inverse CDF
 *  provided in the argument.  Note that this function assumes that
 *  the inverse CDF has equidistant steps in cumulative probability.
 */
template<typename CORE> inline double RandomNumbers_TNG<CORE>::
InverseCDF(const std::vector<Pair> &inv_cdf)
{
  double x = CORE::Double();
  unsigned i = (unsigned)ceil(x*(double)(inv_cdf.size()-1));
  return interpolate(x,inv_cdf.begin()+i);
}

/**
 *  Overwrites the CDF provided in the argument with its inverse.  The
 *  second argument specifies the number of equidistant bins in
 *  cumulative probability that will be generated. The first entry
 *  should be zero, and the last should be one.
 */
template<typename CORE> void RandomNumbers_TNG<CORE>::
GenerateInverseCDF(std::vector<Pair> &cdf, unsigned nbins)
{
  if(nbins == 0)
    nbins = cdf.size();

  for(std::vector<Pair>::iterator itr = cdf.begin();
      itr != cdf.end(); ++itr)
    *itr = std::make_pair(itr->second,itr->first);

  std::vector<Pair> inv_cdf;
  std::vector<Pair>::const_iterator itr = cdf.begin()+1;

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

      inv_cdf.push_back(Pair(x,interpolate(x,itr)));
    }

  cdf = inv_cdf;
}

template<typename T> SimpleRNGAdapter_TNG<T>::~SimpleRNGAdapter_TNG()
{
  // nothing to see here
}

template<typename T> double SimpleRNGAdapter_TNG<T>::uniform()
{
  return m_rng->Double();
};

#endif
