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

namespace calin { namespace math { namespace special {

template<typename T> inline T SQR(T x) { return x*x; }
template<typename T> inline T CUBE(T x) { return x*x*x; }
template<typename T> inline T QUAD(T x) { return SQR(SQR(x)); }

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

} } } // namespace calin::math::special
