/*

   calin/math/lomb_scargle_vcl.hpp -- Stephen Fegan -- 2020-01-31

   Calculate Lomb-Scargle periodogram using VCL vectorization

   Copyright 2020, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <Eigen/Dense>
#include <util/vcl.hpp>

namespace calin { namespace math { namespace lomb_scargle {

template<typename VCLReal, unsigned UNROLL=2>
typename VCLReal::vecX_t periodogram_vcl(
  const typename VCLReal::vecX_t& xi, const typename VCLReal::vecX_t& ti,
  const typename VCLReal::real_t freq_lo, const typename VCLReal::real_t freq_hi,
  const typename VCLReal::real_t delta_freq, unsigned renormalize_nfreq = 0)
{
  const unsigned nfreq_block = (unsigned((freq_hi-freq_lo)/delta_freq)+UNROLL-1)/UNROLL;
  const unsigned nfreq = nfreq_block * UNROLL;

  const unsigned ni_block = (xi.size()+VCLReal::num_real-1)/VCLReal::num_real;
  const unsigned ni = xi.size();

  unsigned ix;
  unsigned iblock;

  renormalize_nfreq = (renormalize_nfreq+UNROLL-1)/UNROLL;

  // Calculate mean
  double xmean = xi.sum()/double(ni);

  // Allocate output vector
  typename VCLReal::vecX_t periodogram(nfreq);

  // Allocate temporary storage
  typename VCLReal::real_t* vx =
    calin::util::memory::aligned_calloc<typename VCLReal::real_t>(5*ni_block*VCLReal::num_real);
  typename VCLReal::real_t* vcf = vx + ni_block*VCLReal::num_real;
  typename VCLReal::real_t* vsf = vx + 2*ni_block*VCLReal::num_real;
  typename VCLReal::real_t* vcdf = vx + 3*ni_block*VCLReal::num_real;
  typename VCLReal::real_t* vsdf = vx + 4*ni_block*VCLReal::num_real;

  typename VCLReal::real_t delta_omega = 2*M_PI*delta_freq;
  typename VCLReal::real_t omega0 = 2*M_PI*freq_lo;

  // Initialize temporary storage
  for(iblock=0, ix=0; iblock<ni_block; iblock++, ix+=VCLReal::num_real) {
    typename VCLReal::real_vt x;
    typename VCLReal::real_vt t;
    typename VCLReal::real_vt mask(1.0);
    if(ix+VCLReal::num_real <= ni) {
      x.load(&xi[ix]);
      t.load(&ti[ix]);
    } else {
      x = 0;
      t = 0;
      x.load_partial(ni-ix, &xi[ix]);
      t.load_partial(ni-ix, &ti[ix]);
      mask.cutoff(ni-ix);
    }
    x -= xmean;
    x.store_a(vx + ix);
    typename VCLReal::real_vt cdf;
    typename VCLReal::real_vt sdf = vcl::sincos(&cdf, t * delta_omega);
    cdf *= mask;
    sdf *= mask;
    cdf.store_a(vcdf + ix);
    sdf.store_a(vsdf + ix);
    typename VCLReal::real_vt cf;
    typename VCLReal::real_vt sf = vcl::sincos(&cf, t * omega0);
    cf *= mask;
    sf *= mask;
    cf.store_a(vcf + ix);
    sf.store_a(vsf + ix);
  }

  // Main loop over frequency blocks
  for(unsigned ifreq=0; ifreq<nfreq_block; ifreq++) {

    // Renormalize the cosine and sine vectors as requested
    if(renormalize_nfreq!=0 and ifreq%renormalize_nfreq==0) {
      const double omega = 2*M_PI*(freq_lo + delta_freq*(ifreq*UNROLL));
      for(iblock=0, ix=0; iblock<ni_block; iblock++, ix+=VCLReal::num_real) {
        typename VCLReal::real_vt t;
        typename VCLReal::real_vt mask(1.0);
        if(ix+VCLReal::num_real <= ni) {
          t.load(&ti[ix]);
        } else {
          t = 0;
          t.load_partial(ni-ix, &ti[ix]);
          mask.cutoff(ni-ix);
        }
        typename VCLReal::real_vt cf;
        typename VCLReal::real_vt sf = vcl::sincos(&cf, t * omega);
        cf *= mask;
        sf *= mask;
        cf.store_a(vcf + ix);
        sf.store_a(vsf + ix);
      }
    }

    // Initialize the accumulators for these frequencies
    typename VCLReal::real_vt CC[UNROLL];
    typename VCLReal::real_vt SS[UNROLL];
    typename VCLReal::real_vt CS[UNROLL];
    typename VCLReal::real_vt XC[UNROLL];
    typename VCLReal::real_vt XS[UNROLL];

    for(unsigned iroll = 0; iroll<UNROLL; ++iroll) {
      CC[iroll] = 0;
      SS[iroll] = 0;
      CS[iroll] = 0;
      XC[iroll] = 0;
      XS[iroll] = 0;
    }

    // Accumulate the data, rotating frequencies
    for(unsigned iblock=0; iblock<ni_block; iblock++) {
      typename VCLReal::real_vt x;
      typename VCLReal::real_vt cdf;
      typename VCLReal::real_vt sdf;
      x.load_a(vx + iblock*VCLReal::num_real);
      cdf.load_a(vcdf + iblock*VCLReal::num_real);
      sdf.load_a(vsdf + iblock*VCLReal::num_real);

      typename VCLReal::real_vt c;
      typename VCLReal::real_vt s;
      c.load_a(vcf + iblock*VCLReal::num_real);
      s.load_a(vsf + iblock*VCLReal::num_real);

      for(unsigned iroll=0; iroll<UNROLL; ++iroll) {
        CC[iroll] = vcl::mul_add(c,c,CC[iroll]);
        CS[iroll] = vcl::mul_add(c,s,CS[iroll]);
        SS[iroll] = vcl::mul_add(s,s,SS[iroll]);
        XC[iroll] = vcl::mul_add(x,c,XC[iroll]);
        XS[iroll] = vcl::mul_add(x,s,XS[iroll]);

        typename VCLReal::real_vt c_ = vcl::mul_sub(c,cdf,s*sdf);
        s = vcl::mul_add(c,sdf,s*cdf);
        c = c_;
      }

      c.store_a(vcf + iblock*VCLReal::num_real);
      s.store_a(vsf + iblock*VCLReal::num_real);
    }

    // Calculate periodogram from accumulators

    for(unsigned iroll=0; iroll<UNROLL; ++iroll) {
      const typename VCLReal::real_t cc = vcl::horizontal_add(CC[iroll]);
      const typename VCLReal::real_t cs = vcl::horizontal_add(CS[iroll]);
      const typename VCLReal::real_t ss = vcl::horizontal_add(SS[iroll]);
      const typename VCLReal::real_t xc = vcl::horizontal_add(XC[iroll]);
      const typename VCLReal::real_t xs = vcl::horizontal_add(XS[iroll]);

      const typename VCLReal::real_t A = (xc*ss-xs*cs)/(ss*cc-cs*cs);
      const typename VCLReal::real_t B = (xs*cc-xc*cs)/(ss*cc-cs*cs);

      periodogram[ifreq*UNROLL + iroll] = xc*A+xs*B;
    }
  }

  free(vx);

  return periodogram;
}

template<typename VCLReal, unsigned UNROLL=2>
typename VCLReal::matX_t multi_periodogram_vcl(
  const typename VCLReal::matX_t& xi, const typename VCLReal::vecX_t& ti,
  const typename VCLReal::real_t freq_lo, const typename VCLReal::real_t freq_hi,
  const typename VCLReal::real_t delta_freq, unsigned renormalize_nfreq = 0)
{
  const unsigned nfreq_block = (unsigned((freq_hi-freq_lo)/delta_freq)+UNROLL-1)/UNROLL;
  const unsigned nfreq = nfreq_block * UNROLL;
  renormalize_nfreq = (renormalize_nfreq+UNROLL-1)/UNROLL;

  const unsigned ni_block = (xi.size()+VCLReal::num_real-1)/VCLReal::num_real;
  const unsigned ni = xi.rows();

  const unsigned nset = xi.cols();

  // Allocate output vector
  typename VCLReal::matX_t periodogram(nfreq,nset);

  // Allocate temporary storage
  typename VCLReal::real_t*__restrict__ vcf0 =
    calin::util::memory::aligned_calloc<typename VCLReal::real_t>(6*ni_block*VCLReal::num_real);
  typename VCLReal::real_t*__restrict__ vsf0 = vcf0 + 1*ni_block*VCLReal::num_real;
  typename VCLReal::real_t*__restrict__ vcdf = vcf0 + 2*ni_block*VCLReal::num_real;
  typename VCLReal::real_t*__restrict__ vsdf = vcf0 + 3*ni_block*VCLReal::num_real;
  typename VCLReal::real_t*__restrict__ vxcf = vcf0 + 4*ni_block*VCLReal::num_real;
  typename VCLReal::real_t*__restrict__ vxsf = vcf0 + 5*ni_block*VCLReal::num_real;

  typename VCLReal::real_t delta_omega = 2*M_PI*delta_freq;
  typename VCLReal::real_t omega0 = 2*M_PI*freq_lo;

  // Initialize frequency storage arrays
  for(unsigned iblock=0; iblock<ni_block; iblock++) {
    const unsigned ix = iblock * VCLReal::num_real;
    typename VCLReal::real_vt t;
    typename VCLReal::real_vt mask(1.0);
    if(ix+VCLReal::num_real <= ni) {
      t.load(&ti[ix]);
    } else {
      t = 0;
      t.load_partial(ni-ix, &ti[ix]);
      mask.cutoff(ni-ix);
    }
    typename VCLReal::real_vt cdf;
    typename VCLReal::real_vt sdf = vcl::sincos(&cdf, t * delta_omega);
    cdf *= mask;
    sdf *= mask;
    cdf.store_a(vcdf + ix);
    sdf.store_a(vsdf + ix);
    typename VCLReal::real_vt cf;
    typename VCLReal::real_vt sf = vcl::sincos(&cf, t * omega0);
    cf *= mask;
    sf *= mask;
    cf.store_a(vcf0 + ix);
    sf.store_a(vsf0 + ix);
  }

  // Precompute the frequency accumulator values
  typename VCLReal::real_t*__restrict__ vcc =
    calin::util::memory::aligned_calloc<typename VCLReal::real_t>(3*nfreq);
  typename VCLReal::real_t*__restrict__ vcs = vcc + 1*nfreq;
  typename VCLReal::real_t*__restrict__ vss = vcc + 2*nfreq;

  std::copy(vcf0, vcf0+ni_block*VCLReal::num_real, vxcf);
  std::copy(vsf0, vsf0+ni_block*VCLReal::num_real, vxsf);

  for(unsigned ifreq=0; ifreq<nfreq; ifreq++) {
    typename VCLReal::real_vt CC = 0;
    typename VCLReal::real_vt SS = 0;
    typename VCLReal::real_vt CS = 0;

    // Renormalize the cosine and sine vectors as requested
    if(ifreq!=0 and renormalize_nfreq!=0 and ifreq%(renormalize_nfreq*UNROLL)==0) {
      const double omega = 2*M_PI*(freq_lo + delta_freq*ifreq);
      for(unsigned iblock=0; iblock<ni_block; ++iblock) {
        unsigned ix = iblock*VCLReal::num_real;
        typename VCLReal::real_vt t;
        typename VCLReal::real_vt mask(1.0);
        if(ix+VCLReal::num_real <= ni) {
          t.load(&ti[ix]);
        } else {
          t = 0;
          t.load_partial(ni-ix, &ti[ix]);
          mask.cutoff(ni-ix);
        }
        typename VCLReal::real_vt cf;
        typename VCLReal::real_vt sf = vcl::sincos(&cf, t * omega);
        cf *= mask;
        sf *= mask;
        cf.store_a(vxcf + ix);
        sf.store_a(vxsf + ix);
      }
    }

    for(unsigned iblock=0; iblock<ni_block; iblock++) {
      typename VCLReal::real_vt cdf;
      typename VCLReal::real_vt sdf;
      cdf.load_a(vcdf + iblock*VCLReal::num_real);
      sdf.load_a(vsdf + iblock*VCLReal::num_real);

      typename VCLReal::real_vt c;
      typename VCLReal::real_vt s;
      c.load_a(vxcf + iblock*VCLReal::num_real);
      s.load_a(vxsf + iblock*VCLReal::num_real);

      CC = vcl::mul_add(c,c,CC);
      CS = vcl::mul_add(c,s,CS);
      SS = vcl::mul_add(s,s,SS);

      typename VCLReal::real_vt c_ = vcl::mul_sub(c,cdf,s*sdf);
      s = vcl::mul_add(c,sdf,s*cdf);
      c = c_;

      c.store_a(vxcf + iblock*VCLReal::num_real);
      s.store_a(vxsf + iblock*VCLReal::num_real);
    }

    vcc[ifreq] = vcl::horizontal_add(CC);
    vcs[ifreq] = vcl::horizontal_add(CS);
    vss[ifreq] = vcl::horizontal_add(SS);
  }

  // Compute periodograms

  for(unsigned iset=0; iset<nset; ++iset) {
    typename VCLReal::real_t xmean = xi.col(iset).sum()/double(ni);

    for(unsigned iblock=0; iblock<ni_block; ++iblock) {
      unsigned ix = iblock * VCLReal::num_real;
      typename VCLReal::real_vt x;
      if(ix+VCLReal::num_real <= ni) {
        x.load(&xi(ix,iset));
      } else {
        x = 0;
        x.load_partial(ni-ix, &xi(ix,iset));
      }
      x -= xmean;
      typename VCLReal::real_vt c;
      c.load_a(vcf0 + ix);
      c *= x;
      c.store_a(vxcf + ix);
      typename VCLReal::real_vt s;
      s.load_a(vsf0 + ix);
      s *= x;
      s.store_a(vxsf + ix);
    }

    typename VCLReal::real_vt XC[UNROLL];
    typename VCLReal::real_vt XS[UNROLL];

    for(unsigned ifreq_block=0; ifreq_block<nfreq_block; ifreq_block++) {
      // Renormalize the cosine and sine vectors as requested
      if(ifreq_block!=0 and renormalize_nfreq!=0 and ifreq_block%renormalize_nfreq==0) {
        const double omega = 2*M_PI*(freq_lo + delta_freq*ifreq_block*UNROLL);
        for(unsigned iblock=0; iblock<ni_block; ++iblock) {
          unsigned ix = iblock * VCLReal::num_real;
          typename VCLReal::real_vt t;
          typename VCLReal::real_vt x;

          if(ix+VCLReal::num_real <= ni) {
            t.load(&ti[ix]);
            x.load(&xi(ix,iset));
          } else {
            t = 0;
            t.load_partial(ni-ix, &ti[ix]);
            x = 0;
            x.load_partial(ni-ix, &xi(ix,iset));
          }
          typename VCLReal::real_vt cf;
          typename VCLReal::real_vt sf = vcl::sincos(&cf, t * omega);
          cf *= x;
          sf *= x;
          cf.store_a(vxcf + ix);
          sf.store_a(vxsf + ix);
        }
      }

      for(unsigned iroll=0;iroll<UNROLL;++iroll) {
        XC[iroll] = 0;
        XS[iroll] = 0;
      }

      for(unsigned iblock=0; iblock<ni_block; iblock++) {
        unsigned ix = iblock * VCLReal::num_real;

        typename VCLReal::real_vt cdf;
        typename VCLReal::real_vt sdf;
        cdf.load_a(vcdf + ix);
        sdf.load_a(vsdf + ix);

        typename VCLReal::real_vt xc;
        typename VCLReal::real_vt xs;
        xc.load_a(vxcf + ix);
        xs.load_a(vxsf + ix);

        for(unsigned iroll=0;iroll<UNROLL;++iroll) {
          XC[iroll] += xc;
          XS[iroll] += xs;

          typename VCLReal::real_vt xc_ = vcl::mul_sub(xc,cdf,xs*sdf);
          xs = vcl::mul_add(xc,sdf,xs*cdf);
          xc = xc_;
        }

        xc.store_a(vxcf + ix);
        xs.store_a(vxsf + ix);
      }

      for(unsigned iroll=0; iroll<UNROLL; ++iroll) {
        unsigned ifreq = ifreq_block*UNROLL + iroll;
        const typename VCLReal::real_t xc = vcl::horizontal_add(XC[iroll]);
        const typename VCLReal::real_t xs = vcl::horizontal_add(XS[iroll]);

        const typename VCLReal::real_t cc = vcc[ifreq];
        const typename VCLReal::real_t cs = vcs[ifreq];
        const typename VCLReal::real_t ss = vss[ifreq];

        const typename VCLReal::real_t A = (xc*ss-xs*cs)/(ss*cc-cs*cs);
        const typename VCLReal::real_t B = (xs*cc-xc*cs)/(ss*cc-cs*cs);

        periodogram(ifreq,iset) = xc*A+xs*B;
      }
    }
  }

  free(vcf0);
  free(vcc);

  return periodogram;
}

} } } // namespace calin::math::lomb_scargle
