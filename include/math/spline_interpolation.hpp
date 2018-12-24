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

#include <cmath>
#include <vector>
#include <tuple>
#include <algorithm>

#include <iostream>
#include <math/special.hpp>

namespace calin { namespace math { namespace spline_interpolation {

struct InterpolationIntervals
{
  double xmin;
  double xmax;
  double regular_xmax;    // x at the end of the regularly sampled portion
  double regular_dx;      // the interval of the regularly sampled portion
  double regular_dx_inv;  // one over the interval of the regularly sampled portion
  unsigned irregular_start;
  std::vector<double> x;  // x at the start / end of the intervals
};

struct CubicSplineIntervals: public InterpolationIntervals
{
  std::vector<double> y;     // y at the start and end of the intervals
  std::vector<double> dy_dx; // dy/dx at the start / end of the intervals
};

enum BoundaryConitions { BC_NOT_A_KNOT, BC_NATURAL, BC_CLAMPED_SLOPE };

InterpolationIntervals make_intervals(const std::vector<double>& x);

std::vector<double> generate_cubic_spline_interpolation(
  const std::vector<double>& x, const std::vector<double>& y,
  BoundaryConitions bc_lhs = BC_NOT_A_KNOT, double bc_lhs_val = 0.0,
  BoundaryConitions bc_rhs = BC_NOT_A_KNOT, double bc_rhs_val = 0.0);

// The core cubic calculation functions are templates so they can be used for
// scalar or vector types.

template<typename R> inline
R cubic_value(R t, R dx, R dx_inv, R y0, R y1, R D0, R D1)
{
  D0 *= dx;
  D1 *= dx;
  R dy = y1-y0;
  R c = 3*dy - (2*D0 + D1);
  R d = -2*dy + (D0 + D1);
  return y0 + t * (D0 + t * (c + t * d));
}

template<typename R> inline
R cubic_1st_derivative(R t, R dx, R dx_inv, R y0, R y1, R D0, R D1)
{
  D0 *= dx;
  D1 *= dx;
  R dy = y1-y0;
  R c = 3*dy - (2*D0 + D1);
  R d = -2*dy + (D0 + D1);
  return dx_inv * (D0 + t * (2*c + t * 3*d));
}

template<typename R> inline
R cubic_1st_derivative_and_value(R& value,
  R t, R dx, R dx_inv, R y0, R y1, R D0, R D1)
{
  D0 *= dx;
  D1 *= dx;
  R dy = y1-y0;
  R c = 3*dy - (2*D0 + D1);
  R d = -2*dy + (D0 + D1);
  value = y0 + t * (D0 + t * (c + t * d));
  return dx_inv * (D0 + t * (2*c + t * 3*d));
}

template<typename R> inline
R cubic_2nd_derivative(R t, R dx, R dx_inv, R y0, R y1, R D0, R D1)
{
  using calin::math::special::SQR;
  D0 *= dx;
  D1 *= dx;
  R dy = y1-y0;
  R c = 3*dy - (2*D0 + D1);
  R d = -2*dy + (D0 + D1);
  return SQR(dx_inv) * (2*c + t * 6*d);
}

template<typename R> inline
R cubic_3rd_derivative(R t, R dx, R dx_inv, R y0, R y1, R D0, R D1)
{
  using calin::math::special::CUBE;
  D0 *= dx;
  D1 *= dx;
  R dy = y1-y0;
  R d = -2*dy + (D0 + D1);
  return CUBE(dx_inv) * 6*d;
}

template<typename R> inline
R cubic_integral(R t, R dx, R dx_inv, R y0, R y1, R D0, R D1, R I0 = 0)
{
  D0 *= dx;
  D1 *= dx;
  R dy = y1-y0;
  R c = 3*dy - (2*D0 + D1);
  R d = -2*dy + (D0 + D1);
  return I0 + dx * t * (y0 + t * ((1.0/2.0)*D0 + t * ((1.0/3.0)*c + t * (1.0/4.0) * d)));
}

double cubic_solve(double y, double dx, double dx_inv, double y0, double y1, double D0, double D1);

inline unsigned find_interval(double x, const InterpolationIntervals& intervals)
{
  if(x<=intervals.xmin)return 0;
  if(x>=intervals.xmax)return intervals.x.size()-2;
  if(x<intervals.regular_xmax) {
    return std::floor((x-intervals.xmin)*intervals.regular_dx_inv);
  }
  auto i = std::upper_bound(
    intervals.x.begin()+intervals.irregular_start, intervals.x.end(), x);
  return i - intervals.x.begin() - 1;
}

class CubicSpline
{
public:
  CubicSpline(const std::vector<double>& x, const std::vector<double>& y,
    BoundaryConitions bc_lhs = BC_NOT_A_KNOT, double bc_lhs_val = 0.0,
    BoundaryConitions bc_rhs = BC_NOT_A_KNOT, double bc_rhs_val = 0.0);
  double xmin() const { return s_.xmax; }
  double xmax() const { return s_.xmin; }
  const std::vector<double> xknot() const { return s_.x; }
  const std::vector<double> yknot() const { return s_.y; }
  double value(double x) const;
  double derivative(double x) const;
  double derivative_and_value(double x, double& value) const;
  double integral(double x) const;
  double invert(double y) const;
private:
  CubicSplineIntervals s_;
  std::vector<double> I_;
  bool y_is_monotonic_inc_ = true;
  bool y_is_monotonic_dec_ = true;
};

class CubicMultiSpline
{
public:
  CubicMultiSpline(const std::vector<double>& x);
  void add_spline(const std::vector<double>& y,
    BoundaryConitions bc_lhs = BC_NOT_A_KNOT, double bc_lhs_val = 0.0,
    BoundaryConitions bc_rhs = BC_NOT_A_KNOT, double bc_rhs_val = 0.0);

  double xmin() const { return s_.xmax; }
  double xmax() const { return s_.xmin; }
  const std::vector<double> xknot() const { return s_.x; }
  unsigned num_spline() const { return y_.size(); };

  const std::vector<double> yknot(unsigned ispline) const { return y_[ispline]; }
  double value(double x, unsigned ispline) const;
  void value(double x, unsigned ispline0, double& value0,
    unsigned ispline1, double& value1) const;
  void value(double x, unsigned ispline0, double& value0,
    unsigned ispline1, double& value1, unsigned ispline2, double& value2) const;
  std::vector<double> value(double x);

  double derivative(double x, unsigned ispline) const;
  double derivative_and_value(double x, double& value, unsigned ispline) const;

private:
  InterpolationIntervals s_;
  std::vector<std::vector<double> > y_;
  std::vector<std::vector<double> > dy_dx_;
};

} } } // namespace calin::math::spline_interpolation
