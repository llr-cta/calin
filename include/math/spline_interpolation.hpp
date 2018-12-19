/*

   calin/math/spline_interpolation.hpp -- Stephen Fegan -- 2018-12-17

   Spline interpolation functions.

   Some of this code from my Fermi tools code "RegularSpline.hpp",
   These portions are copyright Stephen Fegan, 2010.

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

#include <vector>

#include <math/special.hpp>

namespace calin { namespace math { namespace spline_interpolation {

struct InterpolationIntervals
{
  double xmin;
  double xmax;
  double regular_xmax;    // x at the start of the regularly sampled portion
  double regular;         // the interval of the regularly sampled portion
  double regular_dx_inv;  // one over the interval of the regularly sampled portion
  std::vector<double> x;  // x at the start / end of the intervals
};

struct CubicSplineIntervals: public InterpolationIntervals
{
  std::vector<double> y;     // y at the start and end of the intervals
  std::vector<double> dy_dx; // dy/dx at the start / end of the intervals
};

enum BoundaryConitions { BC_NATURAL, BC_NOT_A_KNOT, BC_CLAMPED_SLOPE };

std::vector<double> generate_cubic_spline_interpolation(
  const std::vector<double>& x, const std::vector<double>& y,
  BoundaryConitions bc_lhs = BC_NATURAL, double bc_lhs_val = 0.0,
  BoundaryConitions bc_rhs = BC_NATURAL, double bc_rhs_val = 0.0);

template<typename R> inline
R core_value(R x, R x0, R dx, R dx_inv, R y0, R y1, R D0, R D1)
{
  D0 *= dx;
  D1 *= dx;
  R t = (x-x0)*dx_inv;
  R dy = y1-y0;
  R c = 3*dy - (2*D0 + D1);
  R d = -2*dy + (D0 + D1);
  return y0 + t * (D0 + t * (c + t * d));
}

template<typename R> inline
R core_1st_derivative(R x, R x0, R dx, R dx_inv, R y0, R y1, R D0, R D1)
{
  D0 *= dx;
  D1 *= dx;
  R t = (x-x0)*dx_inv;
  R dy = y1-y0;
  R c = 3*dy - (2*D0 + D1);
  R d = -2*dy + (D0 + D1);
  return dx_inv * (D0 + t * (2*c + t * 3*d));
}

template<typename R> inline
R core_1st_derivative_and_value(R& value,
  R x, R x0, R dx, R dx_inv, R y0, R y1, R D0, R D1)
{
  D0 *= dx;
  D1 *= dx;
  double t = (x-x0)*dx_inv;
  double dy = y1-y0;
  double c = 3*dy - (2*D0 + D1);
  double d = -2*dy + (D0 + D1);
  value = y0 + t * (D0 + t * (c + t * d));
  return dx_inv * (D0 + t * (2*c + t * 3*d));
}

template<typename R> inline
R core_2nd_derivative(R x, R x0, R dx, R dx_inv, R y0, R y1, R D0, R D1)
{
  using calin::math::special::SQR;
  D0 *= dx;
  D1 *= dx;
  R t = (x-x0)*dx_inv;
  R dy = y1-y0;
  R c = 3*dy - (2*D0 + D1);
  R d = -2*dy + (D0 + D1);
  return SQR(dx_inv) * (2*c + t * 6*d);
}

template<typename R> inline
R core_3rd_derivative(R x, R x0, R dx, R dx_inv, R y0, R y1, R D0, R D1)
{
  using calin::math::special::CUBE;
  D0 *= dx;
  D1 *= dx;
  R t = (x-x0)*dx_inv;
  R dy = y1-y0;
  R d = -2*dy + (D0 + D1);
  return CUBE(dx_inv) * 6*d;
}

template<typename R> inline
R core_integral(R x, R x0, R dx, R dx_inv, R y0, R y1, R D0, R D1, R I0 = 0)
{
  D0 *= dx;
  D1 *= dx;
  R t = (x-x0)*dx_inv;
  R dy = y1-y0;
  R c = 3*dy - (2*D0 + D1);
  R d = -2*dy + (D0 + D1);
  return I0 + t * (y0 + t * ((1.0/2.0)*D0 + t * ((1.0/3.0)*c + t * (1.0/4.0) * d)));
}

} } } // namespace calin::math::spline_interpolation
