/*

   calin/math/b_spline.hpp -- Stephen Fegan -- 2018-10-22

   Class for evaluation and fitting of B splines

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

namespace calin { namespace math { namespace b_spline {

template<typename REAL> inline void basis0(REAL x, REAL& y00)
{
  y00 = 1.0;
}

template<typename REAL> inline void basis1(REAL x, REAL& y0, REAL& y1)
{
  y0 = x;
  y1 = 1.0-x;
}

template<typename REAL> inline void basis2(REAL x, REAL& y0, REAL& y1, REAL& y2)
{
  constexpr REAL S2 = 1.0/2.0;

  y1 = 1-x;
  y0 = x;

  y2 = S2 * y1*y1;
  y1 = S2 * ((1+x)*y1 + (2-x)*y0);
  y0 = S2 * x*y0;
}

template<typename REAL> inline void basis3(REAL x, REAL& y0, REAL& y1, REAL& y2, REAL& y3)
{
  constexpr REAL S2 = 1.0/2.0;
  constexpr REAL S3 = 1.0/3.0;

  y3 = 1-x;
  y0 = x;

  y2 = S2 * y3*y3;
  y1 = S2 * ((1+x)*y3 + (2-x)*y0);
  y0 = S2 * x*y0;

  y3 *= S3 * y2;
  y2 = S3 * ((2+x)*y2 + (2-x)*y1);
  y1 = S3 * ((1+x)*y1 + (3-x)*y0);
  y0 = S3 * x*y0;
}

template<typename REAL> inline REAL eval0(REAL x, REAL c0)
{
  return c0;
}

template<typename REAL> inline REAL eval1(REAL x, REAL c0, REAL c1)
{
  c0 = (1-x)*c1 + x*c0;
  return c0;
}

template<typename REAL> inline REAL eval2(REAL x, REAL c0, REAL c1, REAL c2)
{
  constexpr REAL S = 1.0/2.0;

  c0 = (2-x)*c1 + x*c0;
  c1 = (1-x)*c2 + (1+x)*c1;

  c0 = (1-x)*c1 + x*c0;
  return S*c0;
}

template<typename REAL> inline REAL eval3(REAL x, REAL c0, REAL c1, REAL c2, REAL c3)
{
  constexpr REAL S = 1.0/(2.0*3.0);

  c0 = (3-x)*c1 + x*c0;
  c1 = (2-x)*c2 + (1+x)*c1;
  c2 = (1-x)*c3 + (2+x)*c2;

  c0 = (2-x)*c1 + x*c0;
  c1 = (1-x)*c2 + (1+x)*c1;

  c0 = (1-x)*c1 + x*c0;
  return S*c0;
}

// See "The calculation of indefinite integrals of B-Spline" P.W. Gaffney (1974)
// but be careful of the different normalization we use here
template<typename REAL> inline REAL integrate0(REAL x, REAL c0)
{
  return x*c0;
}

template<typename REAL> inline REAL integrate1(REAL x, REAL c0, REAL c1)
{
  double I = 0;
  I += eval1(x, c0*x, c1*(1+x))/2.0;
  I += eval0(x, c1*x);
  return I;
}

template<typename REAL> inline REAL integrate2(REAL x, REAL c0, REAL c1, REAL c2)
{
  double I = 0;
  I += eval2(x, c0*x, c1*(1+x), c2*(2+x))/3.0;
  I += eval1(x, c1*x, c2*(1+x))/2.0;
  I += eval0(x, c2*x);
  return I;
}

template<typename REAL> inline REAL integrate3(REAL x, REAL c0, REAL c1, REAL c2, REAL c3)
{
  double I = 0;
  I += eval3(x, c0*x, c1*(1+x), c2*(2+x), c3*(3+x))/4.0;
  I += eval2(x, c1*x, c2*(1+x), c3*(2+x))/3.0;
  I += eval1(x, c2*x, c3*(1+x))/2.0;
  I += eval0(x, c3*x);
  return I;
}


} } } // namespace calin:math::b_spline
