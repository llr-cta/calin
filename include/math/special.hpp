/*

   calin/math/special.hpp -- Stephen Fegan -- 2015-08-06

   Various special functions, which could be implemented as interfaces
   to GSL or directly

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

#pragma once

#include <gsl/gsl_sf_dawson.h>
#include <gsl/gsl_sf_erf.h>

#include <util/vcl.hpp>

namespace calin { namespace math { namespace special {

template<typename T> constexpr inline T SQR(T x) { return x*x; }
template<typename T> constexpr inline T CUBE(T x) { return x*x*x; }
template<typename T> constexpr inline T QUAD(T x) { return SQR(SQR(x)); }

inline double dawson(double x)
{
  return gsl_sf_dawson(x);
}

inline double lerfc(double x)
{
  return gsl_sf_log_erfc(x);
}

inline unsigned round_up_power_of_two(unsigned v)
{
  // see http://graphics.stanford.edu/~seander/bithacks.html#RoundUpPowerOf2
  v--;
  v |= v >> 1;
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;
  v++;
  return v;
}

#ifndef SWIG
template<typename T> void gaussian_scalar(T* vec, unsigned n, T mean, T sigma)
{
  const T scale = 0.5/SQR(sigma);
  const T norm = std::sqrt(M_1_PI*scale);
  for(unsigned i=0; i<n; i++) {
    T x2 = SQR(T(i) - mean)*scale;
    vec[i] = norm*std::exp(-x2);
  }
}

template<typename VCLReal>
void gaussian_vcl(typename VCLReal::real_t* vec,
  unsigned n, typename VCLReal::real_t mean, typename VCLReal::real_t sigma)
{
  typename VCLReal::real_t scale = 0.5/SQR(sigma);
  typename VCLReal::real_t norm = std::sqrt(M_1_PI*scale);

  int i = 0;
  const int nv = int(n) - VCLReal::num_real;
  typename VCLReal::real_vt x = VCLReal::iota();

  while(i <= nv) {
    typename VCLReal::real_vt x2 = SQR(x - mean)*scale;
    typename VCLReal::real_vt g = norm*vcl::exp(-x2);
    g.store(vec + i);
    x += VCLReal::num_real;
    i += VCLReal::num_real;
  }

  while(i<int(n)) {
    typename VCLReal::real_t x2 = SQR(typename VCLReal::real_t(i) - mean)*scale;
    vec[i] = norm*std::exp(-x2);
    ++i;
  }
}

void gaussian(double* vec, unsigned n, double mean, double sigma);

template<typename T> void two_gaussian_scalar(T* vec, unsigned n, T mean, T sigma, T split)
{
  const T mean_a = mean - 0.5*split;
  const T mean_b = mean + 0.5*split;
  const T scale = 0.5/SQR(sigma);
  const T norm = 0.5*std::sqrt(M_1_PI*scale);
  for(unsigned i=0; i<n; i++) {
    T x2_a = SQR(T(i) - mean_a)*scale;
    T x2_b = SQR(T(i) - mean_b)*scale;
    vec[i] = norm*(std::exp(-x2_a) + std::exp(-x2_b));
  }
}

template<typename VCLReal>
void two_gaussian_vcl(typename VCLReal::real_t* vec,
  unsigned n, typename VCLReal::real_t mean, typename VCLReal::real_t sigma,
  typename VCLReal::real_t split)
{
  typename VCLReal::real_t mean_a = mean - 0.5*split;
  typename VCLReal::real_t mean_b = mean + 0.5*split;
  typename VCLReal::real_t scale = 0.5/SQR(sigma);
  typename VCLReal::real_t norm = 0.5*std::sqrt(M_1_PI*scale);

  int i = 0;
  const int nv = int(n) - VCLReal::num_real;
  typename VCLReal::real_vt x = VCLReal::iota();

  while(i <= nv) {
    typename VCLReal::real_vt x2_a = SQR(x - mean_a)*scale;
    typename VCLReal::real_vt x2_b = SQR(x - mean_b)*scale;
    typename VCLReal::real_vt g = norm*(vcl::exp(-x2_a) + vcl::exp(-x2_b));
    g.store(vec + i);
    x += VCLReal::num_real;
    i += VCLReal::num_real;
  }

  while(i<int(n)) {
    typename VCLReal::real_t x2_a = SQR(typename VCLReal::real_t(i) - mean_a)*scale;
    typename VCLReal::real_t x2_b = SQR(typename VCLReal::real_t(i) - mean_b)*scale;
    vec[i] = norm*(std::exp(-x2_a) + std::exp(-x2_b));
    ++i;
  }
}

void two_gaussian(double* vec, unsigned n, double mean, double sigma, double split);
#endif

Eigen::VectorXd gaussian(unsigned n, double mean, double sigma, bool vcl = true);
Eigen::VectorXd two_gaussian(unsigned n, double mean, double sigma, double split, bool vcl = true);

} } } // namespace calin::math::special
