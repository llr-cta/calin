/*

   calin/math/lomb_scargle.cpp -- Stephen Fegan -- 2017-12-05

   Calculate Lomb-Scargle periodogram

   Copyright 2017, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, CNRS, Ecole Polytechnique, Institut Polytechnique de Paris

   This file is part of "calin"

   "calin" is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 or
   later, as published by the Free Software Foundation.

   "calin" is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

*/

#include <stdexcept>
#include <cmath>

#include <provenance/system_info.hpp>
#include <math/lomb_scargle.hpp>
#include <util/memory.hpp>

namespace {

inline void validate(const Eigen::VectorXd& xi, const Eigen::VectorXd& ti)
{
  if(xi.size() != ti.size())
    throw std::runtime_error("Lomb-Scargle: xi and ti must be same size, "
      + std::to_string(xi.size()) + " != " + std::to_string(ti.size()));
}

} // anonymous namespace

std::pair<double, double> calin::math::lomb_scargle::
amplitudes(const Eigen::VectorXd& xi, const Eigen::VectorXd& ti, double freq)
{
  validate(xi,ti);

  unsigned ni = xi.size();
  double xmean = 0;
  for(unsigned i=0; i<ni; i++){
    xmean += xi[i];
  }
  xmean /= double(ni);

  const double omega = 2*M_PI*freq;
  double CC = 0;
  double SS = 0;
  double CS = 0;
  double XC = 0;
  double XS = 0;
  for(unsigned i=0; i<ni; i++) {
    const double theta = omega*ti[i];
    const double x = xi[i] - xmean;
    const double c = std::cos(theta);
    const double s = std::sin(theta);
    CC += c*c;
    SS += s*s;
    CS += c*s;
    XC += x*c;
    XS += x*s;
  }
  const double A = (XC*SS-XS*CS)/(SS*CC-CS*CS);
  const double B = (XS*CC-XC*CS)/(SS*CC-CS*CS);
  return {A,B};
}

double calin::math::lomb_scargle::
power(const Eigen::VectorXd& xi, const Eigen::VectorXd& ti, double freq)
{
  validate(xi,ti);

  unsigned ni = xi.size();
  double xmean = 0;
  for(unsigned i=0; i<ni; i++){
    xmean += xi[i];
  }
  xmean /= double(ni);

  const double omega = 2*M_PI*freq;
  double CC = 0;
  double SS = 0;
  double CS = 0;
  double XC = 0;
  double XS = 0;
  for(unsigned i=0; i<ni; i++) {
    const double theta = omega*ti[i];
    const double x = xi[i] - xmean;
    const double c = std::cos(theta);
    const double s = std::sin(theta);
    CC += c*c;
    SS += s*s;
    CS += c*s;
    XC += x*c;
    XS += x*s;
  }
  const double A = (XC*SS-XS*CS)/(SS*CC-CS*CS);
  const double B = (XS*CC-XC*CS)/(SS*CC-CS*CS);
  return XC*A+XS*B;
}

Eigen::VectorXd calin::math::lomb_scargle::
periodogram(const Eigen::VectorXd& xi, const Eigen::VectorXd& ti,
  double freq_lo, double freq_hi, double delta_freq)
{
#if defined(__AVX2__) && defined(__FMA__)
  const auto* sysinfo = calin::provenance::system_info::the_host_info();
  if(sysinfo->cpu_has_avx2() and sysinfo->cpu_has_fma3())
    return periodogram_avx2(xi, ti, freq_lo, freq_hi, delta_freq);
#endif
  return periodogram_fast(xi, ti, freq_lo, freq_hi, delta_freq);
}

Eigen::VectorXd calin::math::lomb_scargle::
periodogram_slow(const Eigen::VectorXd& xi, const Eigen::VectorXd& ti,
  double freq_lo, double freq_hi, double delta_freq)
{
  validate(xi,ti);
  unsigned nfreq = unsigned((freq_hi-freq_lo)/delta_freq);
  Eigen::VectorXd periodogram(nfreq);
  for(unsigned ifreq=0; ifreq<nfreq; ifreq++) {
    double freq = freq_lo + delta_freq*ifreq;
    periodogram[ifreq] = power(xi, ti, freq);
  }
  return periodogram;
}

Eigen::VectorXd calin::math::lomb_scargle::
periodogram_fast(const Eigen::VectorXd& xi, const Eigen::VectorXd& ti,
  double freq_lo, double freq_hi, double delta_freq, unsigned renormalize_nfreq)
{
  validate(xi,ti);
  const unsigned nfreq = unsigned((freq_hi-freq_lo)/delta_freq);
  const unsigned ni = xi.size();

  Eigen::VectorXd cdf(ni);
  Eigen::VectorXd sdf(ni);
  const double delta_omega = 2*M_PI*delta_freq;
  double xmean = 0;
  for(unsigned i=0; i<ni; i++) {
    xmean += xi[i];
    const double delta_theta = delta_omega*ti[i];
    cdf[i] = std::cos(delta_theta);
    sdf[i] = std::sin(delta_theta);
  }

  Eigen::VectorXd xi_centered(ni);
  xmean /= double(ni);
  for(unsigned i=0; i<ni; i++) {
    xi_centered[i] = xi[i]-xmean;
  }

  Eigen::VectorXd cf(ni);
  Eigen::VectorXd sf(ni);
  Eigen::VectorXd periodogram(nfreq);
  for(unsigned ifreq=0; ifreq<nfreq; ifreq++) {

    // Renormalize the cosine and sine vectors as requested
    if(ifreq==0 or (renormalize_nfreq!=0 and ifreq%renormalize_nfreq==0)) {
      const double omega = 2*M_PI*(freq_lo + delta_freq*ifreq);
      for(unsigned i=0; i<ni; i++) {
        const double theta = omega*ti[i];
        cf[i] = std::cos(theta);
        sf[i] = std::sin(theta);
      }
    }

    double CC = 0;
    double SS = 0;
    double CS = 0;
    double XC = 0;
    double XS = 0;
    for(unsigned i=0; i<ni; i++) {
      const double c = cf[i];
      const double s = sf[i];

      CC += c*c;
      SS += s*s;
      CS += c*s;
      const double x = xi_centered[i];
      XC += x*c;
      XS += x*s;

      cf[i] = c*cdf[i] - s*sdf[i];
      sf[i] = c*sdf[i] + s*cdf[i];
    }
    const double A = (XC*SS-XS*CS)/(SS*CC-CS*CS);
    const double B = (XS*CC-XC*CS)/(SS*CC-CS*CS);
    periodogram[ifreq] = XC*A+XS*B;
  }
  return periodogram;
}

Eigen::VectorXd calin::math::lomb_scargle::
frequencies(const Eigen::VectorXd& periodogram, double freq_lo, double delta_freq)
{
  unsigned nfreq = periodogram.size();
  Eigen::VectorXd freq(nfreq);
  for(unsigned ifreq=0; ifreq<nfreq; ifreq++)freq[ifreq] = freq_lo + delta_freq*ifreq;
  return freq;
}

#if defined(__AVX2__) && defined(__FMA__)

#include <immintrin.h>

Eigen::VectorXd calin::math::lomb_scargle::
periodogram_avx2(const Eigen::VectorXd& xi, const Eigen::VectorXd& ti,
  double freq_lo, double freq_hi, double delta_freq, unsigned renormalize_nfreq)
{
  // ===========================================================================
  //
  // CAUTION : No User Serviceable Parts Inside
  //
  // On my machine this function is twice as fast as the "_fast" version above.
  // I think this improvement comes only from the fact that it processes two
  // frequencies at once in the loop. The extra improvement of a factor of 4
  // from the vector operations doesn't appear. I think this is because the
  // compiler automatically vectorizes the "_fast" version also. This is a
  // useful lesson ! It probably wasn't worth bothering to do an "intrisics"
  // version. Even the additional factor of 2 could probably be achieved in
  // "_fast" above by restructuring the loop to process two frequencies at once
  // like here.
  //
  // ===========================================================================

  validate(xi,ti);

  const unsigned nfreq_block = (unsigned((freq_hi-freq_lo)/delta_freq)+1)/2;
  const unsigned nfreq = nfreq_block * 2;

  const unsigned ni_block = (xi.size()+3)/4;
  const unsigned ni = xi.size();

  double xmean = 0;
  for(unsigned i=0; i<ni; i++){
    xmean += xi[i];
  }
  xmean /= double(ni);

  Eigen::VectorXd periodogram(nfreq);

  __m256d* vx = calin::util::memory::aligned_calloc<__m256d>(ni_block);
  __m256d* vcf = calin::util::memory::aligned_calloc<__m256d>(ni_block);
  __m256d* vsf = calin::util::memory::aligned_calloc<__m256d>(ni_block);
  __m256d* vcdf = calin::util::memory::aligned_calloc<__m256d>(ni_block);
  __m256d* vsdf = calin::util::memory::aligned_calloc<__m256d>(ni_block);

  double delta_omega = 2*M_PI*delta_freq;

  for(unsigned iblock=0; iblock<ni_block; iblock++) {
    double x[4] = {0,0,0,0};
    double cdf[4] = {0,0,0,0};
    double sdf[4] = {0,0,0,0};

    for(unsigned i=0; i<4; i++) {
      unsigned ichan = iblock*4 + i;
      if(ichan >= ni)break;
      x[i] = xi[ichan]-xmean;
      cdf[i] = std::cos(delta_omega*ti[ichan]);
      sdf[i] = std::sin(delta_omega*ti[ichan]);
    }

    vx[iblock] = _mm256_set_pd(x[3], x[2], x[1], x[0]);
    vcdf[iblock] = _mm256_set_pd(cdf[3], cdf[2], cdf[1], cdf[0]);
    vsdf[iblock] = _mm256_set_pd(sdf[3], sdf[2], sdf[1], sdf[0]);
  }

  renormalize_nfreq /= 2;

  for(unsigned ifreq=0; ifreq<nfreq_block; ifreq++) {

    // Renormalize the cosine and sine vectors as requested
    if(ifreq==0 or (renormalize_nfreq!=0 and ifreq%renormalize_nfreq==0)) {
      const double omega = 2*M_PI*(freq_lo + delta_freq*(ifreq*2));
      for(unsigned iblock=0; iblock<ni_block; iblock++) {
        double cf[4] = {0,0,0,0};
        double sf[4] = {0,0,0,0};
        for(unsigned i=0; i<4; i++) {
          unsigned ichan = iblock*4 + i;
          if(ichan >= ni)break;
          cf[i] = std::cos(omega*ti[ichan]);
          sf[i] = std::sin(omega*ti[ichan]);
        }
        vcf[iblock] = _mm256_set_pd(cf[3], cf[2], cf[1], cf[0]);
        vsf[iblock] = _mm256_set_pd(sf[3], sf[2], sf[1], sf[0]);
      }
    }

    __m256d CC1 = _mm256_setzero_pd(); // YMM0 (notional register assignment)
    __m256d SS1 = _mm256_setzero_pd(); // YMM1 (in fact GCC decides which)
    __m256d CS1 = _mm256_setzero_pd(); // YMM2 (variable is assigned to which)
    __m256d XC1 = _mm256_setzero_pd(); // YMM3 (register.)
    __m256d XS1 = _mm256_setzero_pd(); // YMM4

    __m256d CC2 = _mm256_setzero_pd(); // YMM5
    __m256d SS2 = _mm256_setzero_pd(); // YMM6
    __m256d CS2 = _mm256_setzero_pd(); // YMM7
    __m256d XC2 = _mm256_setzero_pd(); // YMM8
    __m256d XS2 = _mm256_setzero_pd(); // YMM9

    for(unsigned ix=0; ix<ni_block; ix++) {
      const __m256d x = vx[ix];        // YMM10
      const __m256d cdf = vcdf[ix];    // YMM11
      const __m256d sdf = vsdf[ix];    // YMM12

      __m256d c = vcf[ix];             // YMM13
      __m256d s = vsf[ix];             // YMM14

      CC1 = _mm256_fmadd_pd(c,c,CC1);
      CS1 = _mm256_fmadd_pd(c,s,CS1);
      SS1 = _mm256_fmadd_pd(s,s,SS1);
      XC1 = _mm256_fmadd_pd(x,c,XC1);
      XS1 = _mm256_fmadd_pd(x,s,XS1);

      __m256d c_ = _mm256_mul_pd(s,sdf); // YMM15
      c_ = _mm256_fmsub_pd(c,cdf,c_);
      s = _mm256_mul_pd(s,cdf);
      s = _mm256_fmadd_pd(c,sdf,s);

      CC2 = _mm256_fmadd_pd(c_,c_,CC2);
      CS2 = _mm256_fmadd_pd(c_,s,CS2);
      SS2 = _mm256_fmadd_pd(s,s,SS2);
      XC2 = _mm256_fmadd_pd(x,c_,XC2);
      XS2 = _mm256_fmadd_pd(x,s,XS2);

      c = _mm256_mul_pd(s,sdf);
      c = _mm256_fmsub_pd(c_,cdf,c);
      s = _mm256_mul_pd(s,cdf);
      s = _mm256_fmadd_pd(c_,sdf,s);

      vcf[ix] = c;
      vsf[ix] = s;
    }

    CC1 = _mm256_hadd_pd(CC1, CC2);
    CS1 = _mm256_hadd_pd(CS1, CS2);
    SS1 = _mm256_hadd_pd(SS1, SS2);
    XC1 = _mm256_hadd_pd(XC1, XC2);
    XS1 = _mm256_hadd_pd(XS1, XS2);

    double x[4];

    _mm256_store_pd(x, CC1);
    const double cc1 = x[0]+x[2];
    const double cc2 = x[1]+x[3];

    _mm256_store_pd(x, CS1);
    const double cs1 = x[0]+x[2];
    const double cs2 = x[1]+x[3];

    _mm256_store_pd(x, SS1);
    const double ss1 = x[0]+x[2];
    const double ss2 = x[1]+x[3];

    _mm256_store_pd(x, XC1);
    const double xc1 = x[0]+x[2];
    const double xc2 = x[1]+x[3];

    _mm256_store_pd(x, XS1);
    const double xs1 = x[0]+x[2];
    const double xs2 = x[1]+x[3];

    const double A1 = (xc1*ss1-xs1*cs1)/(ss1*cc1-cs1*cs1);
    const double B1 = (xs1*cc1-xc1*cs1)/(ss1*cc1-cs1*cs1);
    periodogram[ifreq*2] = xc1*A1+xs1*B1;

    const double A2 = (xc2*ss2-xs2*cs2)/(ss2*cc2-cs2*cs2);
    const double B2 = (xs2*cc2-xc2*cs2)/(ss2*cc2-cs2*cs2);
    periodogram[ifreq*2+1] = xc2*A2+xs2*B2;
  }

  free(vx);
  free(vcf);
  free(vsf);
  free(vcdf);
  free(vsdf);

  return periodogram;
}

Eigen::VectorXd calin::math::lomb_scargle::
periodogram_avx2_float(const Eigen::VectorXd& xi, const Eigen::VectorXd& ti,
  double freq_lo, double freq_hi, double delta_freq, unsigned renormalize_nfreq)
{
  // ===========================================================================
  //
  // CAUTION : No User Serviceable Parts Inside
  //
  // On my machine this function is twice as fast as the double "_avx" version
  // above, which can be understood since we process 8 floats at once using
  // the vector intrinsics rather than 4 doubles.
  //
  // ===========================================================================

  validate(xi,ti);

  const unsigned nfreq_block = (unsigned((freq_hi-freq_lo)/delta_freq)+1)/2;
  const unsigned nfreq = nfreq_block * 2;

  const unsigned ni_block = (xi.size()+7)/8;
  const unsigned ni = xi.size();

  double xmean = 0;
  for(unsigned i=0; i<ni; i++){
    xmean += xi[i];
  }
  xmean /= double(ni);

  Eigen::VectorXd periodogram(nfreq);

  __m256* vx = calin::util::memory::aligned_calloc<__m256>(ni_block);
  __m256* vcf = calin::util::memory::aligned_calloc<__m256>(ni_block);
  __m256* vsf = calin::util::memory::aligned_calloc<__m256>(ni_block);
  __m256* vcdf = calin::util::memory::aligned_calloc<__m256>(ni_block);
  __m256* vsdf = calin::util::memory::aligned_calloc<__m256>(ni_block);

  double delta_omega = 2*M_PI*delta_freq;

  for(unsigned iblock=0; iblock<ni_block; iblock++) {
    float x[8] = {0,0,0,0,0,0,0,0};
    float cdf[8] = {0,0,0,0,0,0,0,0};
    float sdf[8] = {0,0,0,0,0,0,0,0};

    for(unsigned i=0; i<8; i++) {
      unsigned ichan = iblock*8 + i;
      if(ichan >= ni)break;
      x[i] = xi[ichan]-xmean;
      cdf[i] = std::cos(float(delta_omega*ti[ichan]));
      sdf[i] = std::sin(float(delta_omega*ti[ichan]));
    }

    vx[iblock] = _mm256_set_ps(x[7], x[6], x[5], x[4], x[3], x[2], x[1], x[0]);
    vcdf[iblock] = _mm256_set_ps(cdf[7], cdf[6], cdf[5], cdf[4], cdf[3], cdf[2], cdf[1], cdf[0]);
    vsdf[iblock] = _mm256_set_ps(sdf[7], sdf[6], sdf[5], sdf[4], sdf[3], sdf[2], sdf[1], sdf[0]);
  }

  renormalize_nfreq /= 2;

  for(unsigned ifreq=0; ifreq<nfreq_block; ifreq++) {

    // Renormalize the cosine and sine vectors as requested
    if(ifreq==0 or (renormalize_nfreq!=0 and ifreq%renormalize_nfreq==0)) {
      const double omega = 2*M_PI*(freq_lo + delta_freq*(ifreq*2));
      for(unsigned iblock=0; iblock<ni_block; iblock++) {
        float cf[8] = {0,0,0,0,0,0,0,0};
        float sf[8] = {0,0,0,0,0,0,0,0};
        for(unsigned i=0; i<8; i++) {
          unsigned ichan = iblock*8 + i;
          if(ichan >= ni)break;
          cf[i] = std::cos(float(omega*ti[ichan]));
          sf[i] = std::sin(float(omega*ti[ichan]));
        }
        vcf[iblock] = _mm256_set_ps(cf[7], cf[6], cf[5], cf[4], cf[3], cf[2], cf[1], cf[0]);
        vsf[iblock] = _mm256_set_ps(sf[7], sf[6], sf[5], sf[4], sf[3], sf[2], sf[1], sf[0]);
      }
    }

    __m256 CC1 = _mm256_setzero_ps(); // YMM0 (notional register assignment)
    __m256 SS1 = _mm256_setzero_ps(); // YMM1 (in fact GCC decides which)
    __m256 CS1 = _mm256_setzero_ps(); // YMM2 (variable is assigned to which)
    __m256 XC1 = _mm256_setzero_ps(); // YMM3 (register.)
    __m256 XS1 = _mm256_setzero_ps(); // YMM4

    __m256 CC2 = _mm256_setzero_ps(); // YMM5
    __m256 SS2 = _mm256_setzero_ps(); // YMM6
    __m256 CS2 = _mm256_setzero_ps(); // YMM7
    __m256 XC2 = _mm256_setzero_ps(); // YMM8
    __m256 XS2 = _mm256_setzero_ps(); // YMM9

    for(unsigned ix=0; ix<ni_block; ix++) {
      const __m256 x = vx[ix];        // YMM10
      const __m256 cdf = vcdf[ix];    // YMM11
      const __m256 sdf = vsdf[ix];    // YMM12

      __m256 c = vcf[ix];             // YMM13
      __m256 s = vsf[ix];             // YMM14

      CC1 = _mm256_fmadd_ps(c,c,CC1);
      CS1 = _mm256_fmadd_ps(c,s,CS1);
      SS1 = _mm256_fmadd_ps(s,s,SS1);
      XC1 = _mm256_fmadd_ps(x,c,XC1);
      XS1 = _mm256_fmadd_ps(x,s,XS1);

      __m256 c_ = _mm256_mul_ps(s,sdf); // YMM15
      c_ = _mm256_fmsub_ps(c,cdf,c_);
      s = _mm256_mul_ps(s,cdf);
      s = _mm256_fmadd_ps(c,sdf,s);

      CC2 = _mm256_fmadd_ps(c_,c_,CC2);
      CS2 = _mm256_fmadd_ps(c_,s,CS2);
      SS2 = _mm256_fmadd_ps(s,s,SS2);
      XC2 = _mm256_fmadd_ps(x,c_,XC2);
      XS2 = _mm256_fmadd_ps(x,s,XS2);

      c = _mm256_mul_ps(s,sdf);
      c = _mm256_fmsub_ps(c_,cdf,c);
      s = _mm256_mul_ps(s,cdf);
      s = _mm256_fmadd_ps(c_,sdf,s);

      vcf[ix] = c;
      vsf[ix] = s;
    }

    CC1 = _mm256_hadd_ps(CC1, CC2);
    CS1 = _mm256_hadd_ps(CS1, CS2);
    SS1 = _mm256_hadd_ps(SS1, SS2);
    XC1 = _mm256_hadd_ps(XC1, XC2);
    XS1 = _mm256_hadd_ps(XS1, XS2);

    float x[8];

    _mm256_store_ps(x, CC1);
    const float cc1 = x[0]+x[1]+x[4]+x[5];
    const float cc2 = x[2]+x[3]+x[6]+x[7];

    _mm256_store_ps(x, CS1);
    const float cs1 = x[0]+x[1]+x[4]+x[5];
    const float cs2 = x[2]+x[3]+x[6]+x[7];

    _mm256_store_ps(x, SS1);
    const float ss1 = x[0]+x[1]+x[4]+x[5];
    const float ss2 = x[2]+x[3]+x[6]+x[7];

    _mm256_store_ps(x, XC1);
    const float xc1 = x[0]+x[1]+x[4]+x[5];
    const float xc2 = x[2]+x[3]+x[6]+x[7];

    _mm256_store_ps(x, XS1);
    const float xs1 = x[0]+x[1]+x[4]+x[5];
    const float xs2 = x[2]+x[3]+x[6]+x[7];

    const double A1 = (xc1*ss1-xs1*cs1)/(ss1*cc1-cs1*cs1);
    const double B1 = (xs1*cc1-xc1*cs1)/(ss1*cc1-cs1*cs1);
    periodogram[ifreq*2] = xc1*A1+xs1*B1;

    const double A2 = (xc2*ss2-xs2*cs2)/(ss2*cc2-cs2*cs2);
    const double B2 = (xs2*cc2-xc2*cs2)/(ss2*cc2-cs2*cs2);
    periodogram[ifreq*2+1] = xc2*A2+xs2*B2;
  }

  free(vx);
  free(vcf);
  free(vsf);
  free(vcdf);
  free(vsdf);

  return periodogram;
}

#else // defined(__AVX2__) && defined(__FMA__)

Eigen::VectorXd calin::math::lomb_scargle::
periodogram_avx2(const Eigen::VectorXd& xi, const Eigen::VectorXd& ti,
  double freq_lo, double freq_hi, double delta_freq, unsigned renormalize_nfreq)
{
  throw std::runtime_error("periodogram_avx2: AVX2 support was not available at compile time");
}

Eigen::VectorXd calin::math::lomb_scargle::
periodogram_avx2_float(const Eigen::VectorXd& xi, const Eigen::VectorXd& ti,
  double freq_lo, double freq_hi, double delta_freq, unsigned renormalize_nfreq)
{
  throw std::runtime_error("periodogram_avx2_float: AVX2 support was not available at compile time");
}

#endif
