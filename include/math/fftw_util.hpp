/*

   calin/math/fftw_util.hpp -- Stephen Fegan -- 2016-03-29

   Utility functions for fftw HC types

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include<cmath>
#include<algorithm>

#include<util/log.hpp>
#include<util/vcl.hpp>
#include<math/fftw_util.pb.h>
#include<math/special.hpp>
#include<fftw3.h>

namespace calin { namespace math { namespace fftw_util {

inline unsigned hcvec_num_real(unsigned nsample) {
  // 1 -> 1
  // 2 -> 2
  // 3 -> 2
  // 4 -> 3
  // 5 -> 3
  return std::min(nsample/2+1, nsample); // Min to handle zero case
}

inline unsigned hcvec_num_imag(unsigned nsample) {
  // 1 -> 0
  // 2 -> 0
  // 3 -> 1
  // 4 -> 1
  // 5 -> 2
  return (std::max(nsample,1U)-1)/2; // Max to handle zero case
}

#ifndef SWIG

template<typename T>
void hcvec_scale_and_multiply(T* ovec, const T* ivec1,
  const T* ivec2, unsigned nsample, T scale = 1.0)
{
  T *ro = ovec;
  T *co = ovec + nsample-1;
  const T *ri1 = ivec1;
  const T *ci1 = ivec1 + nsample-1;
  const T *ri2 = ivec2;
  const T *ci2 = ivec2 + nsample-1;

  (*ro++) = (*ri1++) * (*ri2++) * scale;
  if(ro==ri1 or ro==ri2)
  {
    while(ro < co)
    {
      T vri1 = *ri1++;
      T vci1 = *ci1--;
      T vri2 = *ri2++;
      T vci2 = *ci2--;
      (*ro++) = (vri1*vri2 - vci1*vci2)*scale;
      (*co--) = (vri1*vci2 + vci1*vri2)*scale;
    }
   }
  else
  {
    while(ro < co)
    {
      (*ro++) = ((*ri1)*(*ri2) - (*ci1)*(*ci2)) * scale;
      (*co--) = ((*ri1++)*(*ci2--) + (*ci1--)*(*ri2++)) * scale;
    }
  }
  if(ro==co)(*ro) = (*ri1) * (*ri2) * scale;
}

template<typename T>
void hcvec_scale_and_multiply_conj(T* ovec, const T* ivec1,
  const T* ivec2_conj, unsigned nsample, T scale = 1.0)
{
  T *ro = ovec;
  T *co = ovec + nsample-1;
  const T *ri1 = ivec1;
  const T *ci1 = ivec1 + nsample-1;
  const T *ri2 = ivec2_conj;
  const T *ci2 = ivec2_conj + nsample-1;

  (*ro++) = (*ri1++) * (*ri2++) * scale;
  if(ro==ri1 or ro==ri2)
  {
    while(ro < co)
    {
      T vri1 = *ri1++;
      T vci1 = *ci1--;
      T vri2 = *ri2++;
      T vci2 = *ci2--;
      (*ro++) = (vri1*vri2 + vci1*vci2)*scale;
      (*co--) = (vci1*vri2 - vri1*vci2)*scale;
    }
   }
  else
  {
    while(ro < co)
    {
      (*ro++) = ((*ri1)*(*ri2) + (*ci1)*(*ci2)) * scale;
      (*co--) = ((*ci1--)*(*ri2++) - (*ri1++)*(*ci2--)) * scale;
    }
  }
  if(ro==co)(*ro) = (*ri1) * (*ri2) * scale;
}

template<typename T>
void hcvec_scale(T* ovec, unsigned nsample,
  T scale = 1.0)
{
  for(unsigned i=0; i<nsample; ++i) {
    ovec[i] *= scale;
  }
}

template<typename T>
void hcvec_add_scaled(T* ovec, const T* ivec, unsigned nsample,
  T scale = 1.0)
{
  for(unsigned i=0; i<nsample; ++i) {
    ovec[i] += ivec[i]*scale;
  }
}

template<typename T>
void hcvec_set_real(T* ovec, T real_value, unsigned nsample)
{
  const unsigned nreal = hcvec_num_real(nsample);
  for(unsigned i=0; i<nreal; ++i) {
    ovec[i] = real_value;
  }
  for(unsigned i=nreal; i<nsample; ++i) {
    ovec[i] = T(0);
  }
}

template<typename T>
void hcvec_add_real(T* ovec, T real_addand, unsigned nsample)
{
  const unsigned nreal = hcvec_num_real(nsample);
  for(unsigned i=0; i<nreal; ++i) {
    ovec[i] += real_addand;
  }
}

template<typename T>
void hcvec_copy_with_scale_and_add_real(T* ovec, const T* ivec,
  T scale, T real_addand, unsigned nsample)
{
  const unsigned nreal = hcvec_num_real(nsample);
  for(unsigned i=0; i<nreal; ++i) {
    ovec[i] = ivec[i] * scale + real_addand;
  }
  for(unsigned i=nreal; i<nsample; ++i) {
    ovec[i] = ivec[i] * scale;
  }
}

template<typename T>
void hcvec_scale_and_add_real(T* ovec, T scale, T real_addand, unsigned nsample)
{
  const unsigned nreal = hcvec_num_real(nsample);
  for(unsigned i=0; i<nreal; ++i) {
    ovec[i] = ovec[i] * scale + real_addand;
  }
  for(unsigned i=nreal; i<nsample; ++i) {
    ovec[i] = ovec[i] * scale;
  }
}

template<typename T>
T hcvec_sum_real(const T* ivec, unsigned nsample)
{
  const unsigned nreal = hcvec_num_real(nsample);
  T sum_real = T(0);
  for(unsigned i=0; i<nreal; ++i) {
    sum_real += ivec[i];
  }
  return sum_real;
}

template<typename T>
T hcvec_avg_real(const T* ivec, unsigned nsample)
{
  const unsigned nreal = hcvec_num_real(nsample);
  T sum_real = T(0);
  for(unsigned i=0; i<nreal; ++i) {
    sum_real += ivec[i];
  }
  return sum_real/double(nreal);
}

template<typename T>
void hcvec_multiply_and_add_real(T* ovec, const T* ivec1,
  const T* ivec2, T real_addand, unsigned nsample)
{
  T *ro = ovec;
  T *co = ovec + nsample-1;
  const T *ri1 = ivec1;
  const T *ci1 = ivec1 + nsample-1;
  const T *ri2 = ivec2;
  const T *ci2 = ivec2 + nsample-1;

  (*ro++) = (*ri1++) * (*ri2++) + real_addand;
  if(ro==ri1 or ro==ri2)
  {
    while(ro < co)
    {
      T vri1 = *ri1++;
      T vci1 = *ci1--;
      T vri2 = *ri2++;
      T vci2 = *ci2--;
      (*ro++) = (vri1*vri2 - vci1*vci2) + real_addand;
      (*co--) = (vri1*vci2 + vci1*vri2);
    }
   }
  else
  {
    while(ro < co)
    {
      (*ro++) = ((*ri1)*(*ri2) - (*ci1)*(*ci2)) + real_addand;
      (*co--) = ((*ri1++)*(*ci2--) + (*ci1--)*(*ri2++));
    }
  }
  if(ro==co)(*ro) = (*ri1) * (*ri2) + real_addand;
}

template<typename VCLReal>
void hcvec_multiply_and_add_real_vcl(typename VCLReal::real_t* ovec,
  const typename VCLReal::real_t* ivec1, const typename VCLReal::real_t* ivec2,
  typename VCLReal::real_t real_addand, unsigned nsample)
{
  typename VCLReal::real_t* ro = ovec;
  typename VCLReal::real_t* co = ovec + nsample;
  const typename VCLReal::real_t* ri1 = ivec1;
  const typename VCLReal::real_t* ci1 = ivec1 + nsample;
  const typename VCLReal::real_t* ri2 = ivec2;
  const typename VCLReal::real_t* ci2 = ivec2 + nsample;

  (*ro++) = (*ri1++) * (*ri2++) + real_addand;

  while(co - ro >= 2*VCLReal::num_real)
  {
    co -= VCLReal::num_real;
    ci1 -= VCLReal::num_real;
    ci2 -= VCLReal::num_real;
    typename VCLReal::real_vt vri1; vri1.load(ri1);
    typename VCLReal::real_vt vci1; vci1.load(ci1);
    typename VCLReal::real_vt vri2; vri2.load(ri2);
    typename VCLReal::real_vt vci2; vci2.load(ci2);

    typename VCLReal::real_vt vro = (vri1*vri2 - calin::util::vcl::reverse(vci1*vci2)) + real_addand;
    typename VCLReal::real_vt vco = (calin::util::vcl::reverse(vri1)*vci2 + vci1*calin::util::vcl::reverse(vri2));

    vro.store(ro);
    vco.store(co);

    ro += VCLReal::num_real;
    ri1 += VCLReal::num_real;
    ri2 += VCLReal::num_real;
  }

  --co;
  --ci1;
  --ci2;
  while(ro < co)
  {
    double vri1 = *ri1++;
    double vci1 = *ci1--;
    double vri2 = *ri2++;
    double vci2 = *ci2--;
    (*ro++) = (vri1*vri2 - vci1*vci2) + real_addand;
    (*co--) = (vri1*vci2 + vci1*vri2);
  }
  if(ro==co)(*ro) = (*ri1) * (*ri2) + real_addand;
}

// NOTE : This function overrides the template for systems with AVX !!!
#if INSTRSET >= 7
void hcvec_multiply_and_add_real(double* ovec, const double* ivec1,
  const double* ivec2, double real_addand, unsigned nsample);
#endif

template<typename T>
void hcvec_polynomial(T* ovec, const T* ivec, const std::vector<T>& p, unsigned nsample)
{
  if(p.empty()) {
    hcvec_set_real(ovec, T(0), nsample);
    return;
  }

  T *ro = ovec;
  T *co = ovec + nsample;
  const T *ri = ivec;
  const T *ci = ivec + nsample;

  auto pi = p.end();
  --pi;
  T vro = *pi;
  T vco;
  T vri = *ri;
  T vci;
  while(pi != p.begin()) {
    --pi;
    vro = vro * vri + (*pi);
  }
  *ro = vro;

  ++ro, ++ri;
  --co, --ci;

  while(ro < co)
  {
    pi = p.end();
    --pi;
    vro = *pi;
    vco = T(0);

    vri = *ri;
    vci = *ci;

    while(pi != p.begin()) {
      --pi;
      T vvro = vro * vri - vco * vci + (*pi);
      vco    = vro * vci + vco * vri;
      vro   = vvro;
    }
    *ro = vro;
    *co = vco;

    ++ro, ++ri;
    --co, --ci;
  }

  if(ro==co) {
    pi = p.end();
    --pi;
    vro = *pi;
    vri = *ri;
    while(pi != p.begin()) {
      --pi;
      vro = vro * vri + (*pi);
    }
    *ro = vro;
  }
}

// More complex version with vectorization and unrolling of vector loop. In this
// version the order of the inner loops are inverted so that each frequency is
// passed through the full polynomial before moving to the next (to improve
// cache performance)
template<typename VCLReal>
void hcvec_polynomial_vcl(typename VCLReal::real_t* ovec,
  const typename VCLReal::real_t* ivec,
  const std::vector<typename VCLReal::real_t>& p, unsigned nsample)
{
  // No user servicable parts inside

  if(p.empty()) {
    hcvec_set_real(ovec, typename VCLReal::real_t(0), nsample);
    return;
  }

  typename VCLReal::real_t* ro = ovec;
  typename VCLReal::real_t* co = ovec + nsample;
  const typename VCLReal::real_t* ri = ivec;
  const typename VCLReal::real_t* ci = ivec + nsample;

  // Evaluate the zero frequency (real-only) component
  auto pi = p.end();
  --pi;
  *ro = *pi;
  while(pi != p.begin()) {
    --pi;
    *ro = (*ro) * (*ri) + (*pi);
  }

  ++ro, ++ri;

  // Evaluate two AVX vectors of real and complex compnents (i.e. 2 * num_real
  // frequencies) using vector types
  while(co - ro >= 4*VCLReal::num_real)
  {
    pi = p.end();
    typename VCLReal::real_vt vpi = *(--pi);

    typename VCLReal::real_vt vro_a = vpi;
    typename VCLReal::real_vt vco_a = typename VCLReal::real_t(0);
    typename VCLReal::real_vt vro_b = vpi;
    typename VCLReal::real_vt vco_b = typename VCLReal::real_t(0);

    typename VCLReal::real_vt vri_a; vri_a.load(ri);
    ri += VCLReal::num_real;
    typename VCLReal::real_vt vri_b; vri_b.load(ri);
    ri += VCLReal::num_real;

    ci -= VCLReal::num_real;
    typename VCLReal::real_vt vci_a; vci_a.load(ci);
    ci -= VCLReal::num_real;
    typename VCLReal::real_vt vci_b; vci_b.load(ci);

    while(pi != p.begin()) {
      vpi = *(--pi);

      typename VCLReal::real_vt vro_t;

      vro_t  = vro_a*vri_a - calin::util::vcl::reverse(vco_a*vci_a) + vpi;
      vco_a  = calin::util::vcl::reverse(vro_a)*vci_a + vco_a*calin::util::vcl::reverse(vri_a);
      vro_a  = vro_t;

      vro_t  = vro_b*vri_b - calin::util::vcl::reverse(vco_b*vci_b) + vpi;
      vco_b  = calin::util::vcl::reverse(vro_b)*vci_b + vco_b*calin::util::vcl::reverse(vri_b);
      vro_b  = vro_t;
    }

    vro_a.store(ro);
    ro += VCLReal::num_real;
    vro_b.store(ro);
    ro += VCLReal::num_real;

    co -= VCLReal::num_real;
    vco_a.store(co);
    co -= VCLReal::num_real;
    vco_b.store(co);
  }

  // Evaluate any remaining real & complex frequencies that don't fit into a vector
  --co, --ci;

  while(ro < co)
  {
    pi = p.end();
    --pi;
    typename VCLReal::real_t vro = *pi;
    typename VCLReal::real_t vco = typename VCLReal::real_t(0);

    typename VCLReal::real_t vri = *ri;
    typename VCLReal::real_t vci = *ci;

    while(pi != p.begin()) {
      --pi;
      typename VCLReal::real_t vrt;
      vrt = vro * vri - vco * vci + (*pi);
      vco = vro * vci + vco * vri;
      vro = vrt;
    }
    *ro = vro;
    *co = vco;

    ++ro, ++ri;
    --co, --ci;
  }

  // For even numbers of frequencies, finish with final real-only component
  if(ro==co) {
    pi = p.end();
    --pi;
    typename VCLReal::real_t vro = *pi;
    typename VCLReal::real_t vri = *ri;
    while(pi != p.begin()) {
      --pi;
      vro = vro * vri + (*pi);
    }
    *ro = vro;
  }
}

template<typename T>
void hcvec_polynomial_old(T* ovec, const T* ivec, const std::vector<T>& p, unsigned nsample)
{
  if(p.empty()) {
    hcvec_set_real(ovec, T(0), nsample);
    return;
  }

  auto pi = p.end();

  --pi;
  hcvec_set_real(ovec, *pi, nsample);

  while(pi != p.begin()) {
    --pi;
    hcvec_multiply_and_add_real(ovec, ovec, ivec, *pi, nsample);
  }
}

// NOTE : This function overrides the template for systems with AVX !!!
#if INSTRSET >= 7
void hcvec_polynomial(double* ovec, const double* ivec,
  const std::vector<double>& p, unsigned nsample);
#endif

// *****************************************************************************
// *****************************************************************************
//
// Multi stage polynomial
//
// *****************************************************************************
// *****************************************************************************

template<typename T>
void hcvec_multi_stage_polynomial(T* ovec, const T* ivec,
  const std::vector<const std::vector<T>*>& stage_p, unsigned nsample)
{
  T *ro = ovec;
  T *co = ovec + nsample;
  const T *ri = ivec;
  const T *ci = ivec + nsample;

  T vri = *ri;
  T vci;

  for(auto p: stage_p)
  {
    auto pi = p->end();
    --pi;
    T vro = *pi;
    while(pi != p->begin()) {
      --pi;
      vro = vro * vri + (*pi);
    }
    vri = vro;
  }
  *ro = vri;

  ++ro, ++ri;
  --co, --ci;

  while(ro < co)
  {
    vri = *ri;
    vci = *ci;

    for(auto p: stage_p)
    {
      auto pi = p->end();
      --pi;

      T vro = *pi;
      T vco = T(0);

      while(pi != p->begin()) {
        --pi;
        T vvro = vro * vri - vco * vci + (*pi);
        vco    = vro * vci + vco * vri;
        vro   = vvro;
      }
      vri = vro;
      vci = vco;
    }
    *ro = vri;
    *co = vci;

    ++ro, ++ri;
    --co, --ci;
  }

  if(ro==co) {
    vri = *ri;
    for(auto p: stage_p)
    {
      auto pi = p->end();
      --pi;
      T vro = *pi;
      while(pi != p->begin()) {
        --pi;
        vro = vro * vri + (*pi);
      }
      vri = vro;
    }
    *ro = vri;
  }
}

// Extremely complex version with vectorization and unrolling of vector loop. In
// this version the order of the inner loops are inverted so that each frequency
// is passed through all the polynomials before moving to the next (to improve
// cache performance)
template<typename VCLReal>
void hcvec_multi_stage_polynomial_vcl(typename VCLReal::real_t* ovec,
  const typename VCLReal::real_t* ivec,
  const std::vector<const std::vector<typename VCLReal::real_t>*>& stage_p, unsigned nsample)
{
  // No user servicable parts inside

  typename VCLReal::real_t* ro = ovec;
  typename VCLReal::real_t* co = ovec + nsample;
  const typename VCLReal::real_t* ri = ivec;
  const typename VCLReal::real_t* ci = ivec + nsample;

  // Evaluate the zero frequency (real-only) component

  typename VCLReal::real_t sri = *ri;
  typename VCLReal::real_t sci;
  for(auto p: stage_p)
  {
    auto pi = p->end();
    --pi;
    typename VCLReal::real_t sro = *pi;
    while(pi != p->begin()) {
      --pi;
      sro = sro * sri + (*pi);
    }
    sri = sro;
  }
  *ro = sri;
  ++ro, ++ri;

  // Evaluate two AVX vectors of real and complex compnents (i.e. 2 * num_real
  // frequencies) using vector types
  while(co - ro >= 6*VCLReal::num_real)
  {
    typename VCLReal::real_vt vri_a; vri_a.load(ri);
    ri += VCLReal::num_real;
    typename VCLReal::real_vt vri_b; vri_b.load(ri);
    ri += VCLReal::num_real;
    typename VCLReal::real_vt vri_c; vri_c.load(ri);
    ri += VCLReal::num_real;

    ci -= VCLReal::num_real;
    typename VCLReal::real_vt vci_a; vci_a.load(ci);
    ci -= VCLReal::num_real;
    typename VCLReal::real_vt vci_b; vci_b.load(ci);
    ci -= VCLReal::num_real;
    typename VCLReal::real_vt vci_c; vci_c.load(ci);

    for(auto p: stage_p)
    {
      auto pi = p->end();

      typename VCLReal::real_vt vpi = *(--pi);

      typename VCLReal::real_vt vro_a = vpi;
      typename VCLReal::real_vt vco_a = typename VCLReal::real_t(0);
      typename VCLReal::real_vt vro_b = vpi;
      typename VCLReal::real_vt vco_b = typename VCLReal::real_t(0);
      typename VCLReal::real_vt vro_c = vpi;
      typename VCLReal::real_vt vco_c = typename VCLReal::real_t(0);

      while(pi != p->begin()) {
        vpi = *(--pi);

        typename VCLReal::real_vt vro_t;

        vro_t  = vro_a*vri_a - calin::util::vcl::reverse(vco_a*vci_a) + vpi;
        vco_a  = calin::util::vcl::reverse(vro_a)*vci_a + vco_a*calin::util::vcl::reverse(vri_a);
        vro_a  = vro_t;

        vro_t  = vro_b*vri_b - calin::util::vcl::reverse(vco_b*vci_b) + vpi;
        vco_b  = calin::util::vcl::reverse(vro_b)*vci_b + vco_b*calin::util::vcl::reverse(vri_b);
        vro_b  = vro_t;

        vro_t  = vro_c*vri_c - calin::util::vcl::reverse(vco_c*vci_c) + vpi;
        vco_c  = calin::util::vcl::reverse(vro_c)*vci_c + vco_c*calin::util::vcl::reverse(vri_c);
        vro_c  = vro_t;
      }

      vri_a = vro_a;
      vci_a = vco_a;
      vri_b = vro_b;
      vci_b = vco_b;
      vri_c = vro_c;
      vci_c = vco_c;
    }

    vri_a.store(ro);
    ro += VCLReal::num_real;
    vri_b.store(ro);
    ro += VCLReal::num_real;
    vri_c.store(ro);
    ro += VCLReal::num_real;

    co -= VCLReal::num_real;
    vci_a.store(co);
    co -= VCLReal::num_real;
    vci_b.store(co);
    co -= VCLReal::num_real;
    vci_c.store(co);
  }

  // Evaluate any remaining real & complex frequencies that don't fit into a vector
  --co, --ci;

  while(ro < co)
  {
    sri = *ri;
    sci = *ci;

    for(auto p: stage_p)
    {
      auto pi = p->end();
      --pi;

      typename VCLReal::real_t sro = *pi;
      typename VCLReal::real_t sco = 0;

      while(pi != p->begin()) {
        --pi;
        typename VCLReal::real_t sro_t = sro * sri - sco * sci + (*pi);
        sco     = sro * sci + sco * sri;
        sro     = sro_t;
      }
      sri = sro;
      sci = sco;
    }
    *ro = sri;
    *co = sci;

    ++ro, ++ri;
    --co, --ci;
  }

  // For even numbers of frequencies, finish with final real-only component
  if(ro==co) {
    sri = *ri;
    for(auto p: stage_p)
    {
      auto pi = p->end();
      --pi;
      typename VCLReal::real_t sro = *pi;
      while(pi != p->begin()) {
        --pi;
        sro = sro * sri + (*pi);
      }
      sri = sro;
    }
    *ro = sri;
  }
}

// NOTE : This function overrides the template for systems with AVX !!!
#if INSTRSET >= 7
void hcvec_multi_stage_polynomial(double* ovec, double* ivec,
  const std::vector<const std::vector<double>*>& stage_p, unsigned nsample);
#endif

// *****************************************************************************
// *****************************************************************************
//
// Analytic DFT of Gaussian
//
// *****************************************************************************
// *****************************************************************************

template<typename T>
void hcvec_gaussian_dft(T* ovec, T mean, T sigma, unsigned nsample)
{
  T *ro = ovec;
  T *co = ovec + nsample-1;

  T nsample_inv = 1.0/T(nsample);
  T scale = 2.0*calin::math::special::SQR(M_PI*sigma*nsample_inv);
  T phase = 2*M_PI*mean*nsample_inv;

  (*ro++) = 1.0;

  while(ro < co)
  {
    T x = T(ro-ovec);
    T amp = std::exp(-calin::math::special::SQR(x) * scale);
    x *= phase;
    T c = std::cos(x);
    T s = std::sin(x);
    (*ro++) = amp*c;
    (*co--) = -amp*s;
  }

  if(ro==co) {
    T x = T(ro-ovec);
    T amp = std::exp(-calin::math::special::SQR(x) * scale);
    x *= phase;
    T c = std::cos(x);
    (*ro++) = amp*c;
  }
}

template<typename VCLReal>
void hcvec_gaussian_dft_vcl(typename VCLReal::real_t* ovec,
  typename VCLReal::real_t mean, typename VCLReal::real_t sigma, unsigned nsample)
{
  // No user servicable parts inside

  typename VCLReal::real_t* ro = ovec;
  typename VCLReal::real_t* co = ovec + nsample;

  typename VCLReal::real_t nsample_inv = 1.0/typename VCLReal::real_t(nsample);
  typename VCLReal::real_t scale = 2.0*calin::math::special::SQR(M_PI*sigma*nsample_inv);
  typename VCLReal::real_t phase = -2*M_PI*mean*nsample_inv;

  // Evaluate the zero frequency (real-only) component
  (*ro++) = 1.0;

  // Evaluate AVX vectors of real and complex compnents (i.e. num_real
  // frequencies) using vector types
  typename VCLReal::real_vt x = VCLReal::iota() + 1.0;
  while(co - ro >= 2*VCLReal::num_real)
  {
    typename VCLReal::real_vt amp = vcl::exp(-x*x*scale);
    typename VCLReal::real_vt c;
    typename VCLReal::real_vt s = vcl::sincos(&c, x*phase);

    c *= amp;
    s *= amp;

    s = calin::util::vcl::reverse(s);

    x += VCLReal::num_real;

    c.store(ro);
    ro += VCLReal::num_real;

    co -= VCLReal::num_real;
    s.store(co);
  }

  // Evaluate any remaining real & complex frequencies that don't fit into a vector
  --co;

  while(ro < co)
  {
    typename VCLReal::real_t x = typename VCLReal::real_t(ro-ovec);
    typename VCLReal::real_t amp = std::exp(-calin::math::special::SQR(x) * scale);
    x *= phase;
    typename VCLReal::real_t c = std::cos(x);
    typename VCLReal::real_t s = std::sin(x);
    (*ro++) = amp*c;
    (*co--) = amp*s;
  }

  if(ro==co) {
    typename VCLReal::real_t x = double(ro-ovec);
    typename VCLReal::real_t amp = std::exp(-calin::math::special::SQR(x) * scale);
    x *= phase;
    typename VCLReal::real_t c = std::cos(x);
    (*ro++) = amp*c;
  }
}

// NOTE : This function overrides the template for systems with AVX !!!
#if INSTRSET >= 7
void hcvec_gaussian_dft(double* ovec, double mean, double sigma, unsigned nsample);
#endif

// *****************************************************************************
// *****************************************************************************
//
// Analytic DFT of Delta(x-x0)
//
// *****************************************************************************
// *****************************************************************************

template<typename T>
void hcvec_delta_dft(T* ovec, T x0, unsigned nsample)
{
  T *ro = ovec;
  T *co = ovec + nsample-1;

  T nsample_inv = 1.0/T(nsample);
  T phase = 2*M_PI*x0*nsample_inv;

  (*ro++) = 1.0;

  while(ro < co)
  {
    T x = T(ro-ovec) * phase;
    T c = std::cos(x);
    T s = std::sin(x);
    (*ro++) = c;
    (*co--) = -s;
  }

  if(ro==co) {
    T x = T(ro-ovec) * phase;
    T c = std::cos(x);
    (*ro++) = c;
  }
}

template<typename VCLReal>
void hcvec_delta_dft_vcl(typename VCLReal::real_t* ovec, typename VCLReal::real_t x0, unsigned nsample)
{
  // No user servicable parts inside

  typename VCLReal::real_t* ro = ovec;
  typename VCLReal::real_t* co = ovec + nsample;

  typename VCLReal::real_t nsample_inv = 1.0/typename VCLReal::real_t(nsample);
  typename VCLReal::real_t phase = 2*M_PI*x0*nsample_inv;

  // Evaluate the zero frequency (real-only) component
  (*ro++) = 1.0;

  // Evaluate three AVX vectors of real and complex compnents (i.e. 3 * num_real
  // frequencies) using vector types
  typename VCLReal::real_vt x = VCLReal::iota() + 1.0;
  while(co - ro >= 2*VCLReal::num_real)
  {
    typename VCLReal::real_vt c;
    typename VCLReal::real_vt s = vcl::sincos(&c, x*phase);

    s = -calin::util::vcl::reverse(s);
    x += VCLReal::num_real;

    c.store(ro);
    ro += VCLReal::num_real;

    co -= VCLReal::num_real;
    s.store(co);
  }

  // Evaluate any remaining real & complex frequencies that don't fit into a vector
  --co;

  while(ro < co)
  {
    typename VCLReal::real_t x = typename VCLReal::real_t(ro-ovec) * phase;
    typename VCLReal::real_t c = std::cos(x);
    typename VCLReal::real_t s = std::sin(x);
    (*ro++) = c;
    (*co--) = -s;
  }

  if(ro==co) {
    typename VCLReal::real_t x = double(ro-ovec) * phase;
    typename VCLReal::real_t c = std::cos(x);
    (*ro++) = c;
  }
}

// NOTE : This function overrides the template for systems with AVX !!!
#if INSTRSET >= 7
void hcvec_delta_dft(double* ovec, double x0, unsigned nsample);
#endif

// *****************************************************************************
// *****************************************************************************
//
// Various other functions and SWIG definitions
//
// *****************************************************************************
// *****************************************************************************

using uptr_fftw_plan = std::unique_ptr<fftw_plan_s,void(*)(fftw_plan_s*)>;
using uptr_fftw_data = std::unique_ptr<double,void(*)(void*)>;

#endif // defined SWIG

// Expose some of these functions in more SWIG friendly way
double hcvec_sum_real(const Eigen::VectorXd& ivec);
double hcvec_avg_real(const Eigen::VectorXd& ivec);
Eigen::VectorXd hcvec_scale_and_add_real(const Eigen::VectorXd& ivec, double scale, double real_addand);
Eigen::VectorXd hcvec_gaussian_dft(double mean, double sigma, unsigned nsample, bool vcl = true);
Eigen::VectorXd hcvec_delta_dft(double x0, unsigned nsample, bool vcl = true);

Eigen::VectorXd fftw_r2hc(const Eigen::VectorXd& x,
  calin::ix::math::fftw_util::FFTWPlanningRigor fftw_rigor = calin::ix::math::fftw_util::ESTIMATE);
Eigen::VectorXd fftw_hc2r(const Eigen::VectorXd& f,
  calin::ix::math::fftw_util::FFTWPlanningRigor fftw_rigor = calin::ix::math::fftw_util::ESTIMATE);

int proto_planning_enum_to_fftw_flag(calin::ix::math::fftw_util::FFTWPlanningRigor x);

bool load_wisdom_from_file(std::string filename = "~/.calin_fft_wisdom");
bool load_wisdom_from_proto(const calin::ix::math::fftw_util::FFTWWisdom& proto);

bool save_wisdom_to_file(std::string filename = "~/.calin_fft_wisdom");
bool save_wisdom_to_proto(calin::ix::math::fftw_util::FFTWWisdom& proto);

} } } // namespace calin::math::fftw_util
