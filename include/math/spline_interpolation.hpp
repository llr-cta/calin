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

namespace calin { namespace math { namespace spline_interpolation {

struct InterpolatedIntervals
{
  double xmin;
  double xmax;
  double regular_xmax;    // x at the start of the regularly sampled portion
  double regular_dx_inv;  // one over the interval of the regularly sampled portion
  std::vector<double> x;  // x at the start / end of the intervals
  std::vector<double> y; // y at the start / end of the intervals
  std::vector<double> dx_inv;  // 1/dx for each intervals
};

struct CubicSplineIntervals: public InterpolatedIntervals
{
  std::vector<double> ypp_6; // (d2y/dx^2)/6 at the start / end of the intervals
};

enum BoundaryConitions { BC_NATURAL, BC_NOT_A_KNOT, BC_CLAMPED_SLOPE };

std::vector<double> generate_cubic_spline_interpolation(
  const std::vector<double>& x, const std::vector<double>& y,
  BoundaryConitions bc_lhs = BC_NATURAL, double bc_lhs_val = 0.0,
  BoundaryConitions bc_rhs = BC_NATURAL, double bc_rhs_val = 0.0);

} } } // namespace calin::math::spline_interpolation
