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

#include <stdexcept>
#include <iostream>

#include <math/spline_interpolation.hpp>
#include <math/special.hpp>

using namespace calin::math::spline_interpolation;

using calin::math::special::SQR;
using calin::math::special::CUBE;
using calin::math::special::QUAD;

InterpolationIntervals
calin::math::spline_interpolation::make_intervals(const std::vector<double>& x)
{
  if(x.size() < 2)
    throw std::runtime_error("make_intervals: need at least two x values");
  InterpolationIntervals intervals;
  intervals.xmin = x[0];

  intervals.regular_xmax = x[1];
  intervals.regular_dx = x[1] - x[0];
  intervals.regular_dx_inv = 1/intervals.regular_dx;
  intervals.irregular_start = 1;

  if(intervals.regular_dx <= 0)
    std::runtime_error("make_intervals: x values must be monotonously increasing");
  for(unsigned i=2; i<x.size(); i++) {
    double dx = x[i] - x[i-1];
    if(dx <= 0)
      std::runtime_error("make_intervals: x values must be monotonously increasing");
    if(intervals.irregular_start == (i-1) and
        abs(dx - intervals.regular_dx)/intervals.regular_dx < 1e-6) {
      intervals.irregular_start = i;
      intervals.regular_xmax = x[i];
    }
  }
  intervals.xmax = x.back();
  intervals.x = x;
  if(intervals.irregular_start == 1) {
    intervals.regular_xmax = x[0];
    intervals.regular_dx = 0;
    intervals.regular_dx_inv = 0;
    intervals.irregular_start = 0;
  }
  return intervals;
}

std::vector<double> calin::math::spline_interpolation::
generate_cubic_spline_interpolation(
  const std::vector<double>& x, const std::vector<double>& y,
  BoundaryConitions bc_lhs, double bc_lhs_val,
  BoundaryConitions bc_rhs, double bc_rhs_val)
{
  unsigned np = y.size();

  if(x.size() != np)
    throw std::runtime_error("generate_cubic_spline_interpolation: x and y must "
      "have equal number of points");

  if(np<3)
    throw std::runtime_error("generate_cubic_spline_interpolation: must "
      "have at least three points to interpolate");

  std::vector<double> dx_inv(np-1, 0); // 1/dx
  for(unsigned i=0; i<np-1; i++)dx_inv[i] = 1.0/(x[i+1]-x[i]);

  std::vector<double> dy(np-1, 0);
  for(unsigned i=0; i<np-1; i++)dy[i] = y[i+1]-y[i];

  // First solve the tridiagonal matrix equation for the slopes at
  // each of the points.. do this in a manner so that we can later
  // calculate the value of the spline and the derivative with respect
  // to the input points and slopes

  std::vector<double> a(np, 0); // Sub diagnoal elements   [1,np-1]
  std::vector<double> b(np, 0); // Diagonal elements       [0,np-1]
  std::vector<double> c(np, 0); // Super diagonal elements [0,np-2]
  std::vector<double> d(np, 0); // RHS of equation         [0,np-1]
  std::vector<double> m(np, 0); // Slopes we are seeking   [0,np-1]

  // Enter primary spline constraints into matrix
  for(unsigned i=1; i<np-1; i++) {
    b[i] = 2*(dx_inv[i]+dx_inv[i-1]);
    a[i] = dx_inv[i-1];
    c[i] = dx_inv[i];
    d[i] = 3*(dy[i]*SQR(dx_inv[i]) + dy[i-1]*SQR(dx_inv[i-1]));
  }

  // We cannot use the default NOT_A_KNOT algorithm with only 2 or 3
  // points (i.e. 1 or 2 spline segments), since there is no "knot" in
  // the first case, giving no constraints at all, and only one "knot"
  // in the second case giving a single constraint.

  // In the 2-point case, if either of the BCS is set to NOT_A_KNOT we
  // zero out the cubic term. If both are set to zero then we also set
  // the quadratic term to zero (by changing bcsN to BCS_NATURAL).
  // This has the pleasent effect of making the spline into a line
  // segment if either is NOT_A_KNOT and the other is either NATURAL
  // or NOT_A_KNOT. In the case of NOT_A_KNOT combined with CLAMPED
  // (preassigned slope) we get a quadratic.

  // In the 3-point case, if either of the BCS is set to NOT_A_KNOT we
  // set the quadratic term to be constant across the "knot" as
  // expected. The other constraint then depends on what the other
  // BCS is. For NATURAL and CLAMPED the other constraint is obvious.
  // If the other is also NOT_A_KNOT we set the cubic term to zero to
  // get a quadratic that fits the three points.

  switch(bc_lhs)
  {
  case BC_NATURAL:
    b[0] = 2;
    c[0] = 1;
    d[0] = 3*dy[0]*dx_inv[0];
    break;
  case BC_CLAMPED_SLOPE:
    b[0] = 1;
    c[0] = 0;
    d[0] = bc_lhs_val;
    break;
  case BC_NOT_A_KNOT:
    b[0] = dx_inv[0]*(dx_inv[1]+dx_inv[0]);
    c[0] = SQR(dx_inv[1]+dx_inv[0]);
    d[0] = dy[1]*CUBE(dx_inv[1]) + dy[0]*(3*dx_inv[1]+2*dx_inv[0])*SQR(dx_inv[0]);
    break;
  }

  switch(bc_rhs)
  {
  case BC_NATURAL:
    b[np-1] = 2;
    a[np-1] = 1;
    d[np-1] = 3*dy[np-2]*dx_inv[np-2];
    break;
  case BC_CLAMPED_SLOPE:
    b[np-1] = 1;
    a[np-1] = 0;
    d[0] = bc_rhs_val;
    break;
  case BC_NOT_A_KNOT:
    if((np==3)and(bc_lhs==BC_NOT_A_KNOT)) {
      b[2] = 1;
      a[2] = 1;
      d[2] = 2*dy[1]*dx_inv[1];
    } else {
      b[np-1] = dx_inv[np-2]*(dx_inv[np-2]+dx_inv[np-3]);
      a[np-1] = SQR(dx_inv[np-2]+dx_inv[np-3]);
      d[np-1] = dy[np-2]*(2*dx_inv[np-2]+3*dx_inv[np-3])*SQR(dx_inv[np-2]) + dy[np-3]*CUBE(dx_inv[np-3]);
    }
    break;
  }

  c[0] /= b[0];
  d[0] /= b[0];

  for(unsigned ip = 1; ip<np; ip++) {
    double id = 1.0/(b[ip] - c[ip-1]*a[ip]);
    c[ip] *= id;
    d[ip] = (d[ip] - d[ip-1] * a[ip])*id;
  }

  m[np-1] = d[np-1];
  for(int ip=np-2; ip>=0; ip--) {
    m[ip] = d[ip]-c[ip]*m[ip+1];
  }

  return m;

  // std::vector<double> ypp(np);
  // for(unsigned ip=0; ip<np-1; ip++) {
  //   ypp[ip] = 6*dy[ip]*SQR(dx_inv[ip]) - (4*m[ip]+2*m[ip+1])*dx_inv[ip];
  // }
  // ypp[np-1] = -6*dy[np-2]*SQR(dx_inv[np-2]) + (2*m[np-2]+4*m[np-1])*dx_inv[np-2];
  //
  // return ypp;
}

CubicSpline::
CubicSpline(const std::vector<double>& x, const std::vector<double>& y,
    BoundaryConitions bc_lhs, double bc_lhs_val,
    BoundaryConitions bc_rhs, double bc_rhs_val):
  s_(), I_()
{
  *static_cast<InterpolationIntervals*>(&s_) = make_intervals(x);
  s_.y = y;
  s_.dy_dx = generate_cubic_spline_interpolation(x, y, bc_lhs, bc_lhs_val, bc_rhs, bc_rhs_val);
  double I=0;
  I_.push_back(I);
  for(unsigned i=0; i<x.size()-1; i++) {
    double dx = s_.x[i+1]-s_.x[i];
    double dx_inv = 1.0/dx;
    I += cubic_integral(1.0, dx, dx_inv, s_.y[i], s_.y[i+1], s_.dy_dx[i], s_.dy_dx[i+1]);
    I_.push_back(I);
  }
}

double CubicSpline::value(double x) const
{
  unsigned i = find_interval(x, s_);
  double dx = s_.x[i+1]-s_.x[i];
  double dx_inv = 1.0/dx;
  double t = (x-s_.x[i])*dx_inv;
  return cubic_value(t, dx, dx_inv, s_.y[i], s_.y[i+1], s_.dy_dx[i], s_.dy_dx[i+1]);
}

double CubicSpline::derivative(double x) const
{
  unsigned i = find_interval(x, s_);
  double dx = s_.x[i+1]-s_.x[i];
  double dx_inv = 1.0/dx;
  double t = (x-s_.x[i])*dx_inv;
  return cubic_1st_derivative(t, dx, dx_inv, s_.y[i], s_.y[i+1], s_.dy_dx[i], s_.dy_dx[i+1]);
}

double CubicSpline::integral(double x) const
{
  unsigned i = find_interval(x, s_);
  double dx = s_.x[i+1]-s_.x[i];
  double dx_inv = 1.0/dx;
  double t = (x-s_.x[i])*dx_inv;
  return cubic_integral(t, dx, dx_inv, s_.y[i], s_.y[i+1], s_.dy_dx[i], s_.dy_dx[i+1], I_[i]);
}

double CubicSpline::invert(double y) const
{
  unsigned i;
  if(s_.y.back() > s_.y.front()) {
    auto ifind = std::upper_bound(s_.y.begin(), s_.y.end(), y);
    if(ifind == s_.y.end())return NAN;
    if(ifind == s_.y.begin())return NAN;
    i = ifind - s_.y.begin() - 1;
  } else {
    auto ifind = std::upper_bound(s_.y.begin(), s_.y.end(), y, std::greater<double>());
    if(ifind == s_.y.end())return NAN;
    if(ifind == s_.y.begin())return NAN;
    i = ifind - s_.y.begin() - 1;
  }
  double dx = s_.x[i+1]-s_.x[i];
  double dx_inv = 1.0/dx;
  double t = cubic_solve(y, dx, dx_inv, s_.y[i], s_.y[i+1], s_.dy_dx[i], s_.dy_dx[i+1]);
  return t*dx + s_.x[i];
}
