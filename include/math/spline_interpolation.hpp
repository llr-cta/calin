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
#include <util/vcl.hpp>

namespace calin { namespace math { namespace spline_interpolation {

struct InterpolationIntervals
{
  double xmin;
  double xmax;
  double regular_xmax;    // x at the end of the regularly sampled portion
  double regular_dx;      // the interval of the regularly sampled portion
  double regular_dx_inv;  // one over the interval of the regularly sampled portion
  unsigned irregular_begin;
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
R cubic_2nd_derivative_and_value(R& first_derivative, R& value,
  R t, R dx, R dx_inv, R y0, R y1, R D0, R D1)
{
  using calin::math::special::SQR;
  D0 *= dx;
  D1 *= dx;
  R dy = y1-y0;
  R c = 3*dy - (2*D0 + D1);
  R d = -2*dy + (D0 + D1);
  value = y0 + t * (D0 + t * (c + t * d));
  first_derivative = dx_inv * (D0 + t * (2*c + t * 3*d));
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
    intervals.x.begin()+intervals.irregular_begin, intervals.x.end(), x);
  return i - intervals.x.begin() - 1;
}

template<typename VCLReal> inline typename VCLReal::int_vt
vcl_find_interval(typename VCLReal::real_vt x, const InterpolationIntervals& intervals)
{
  typename VCLReal::int_vt ireg =
    VCLReal::round_to_int(vcl::floor((x-intervals.xmin)*intervals.regular_dx_inv));
  ireg = vcl::max(ireg, 0);
  unsigned iregmax = intervals.x.size()-2;
  if(intervals.irregular_begin <= iregmax)
  {
    if(vcl::horizontal_or(ireg > intervals.irregular_begin)) {
      // Yucky scalar code, sniff :'(
      typename VCLReal::real_at x_array;
      typename VCLReal::int_at ireg_array;
      x.store_a(x_array);
      ireg.store_a(ireg_array);
      for(unsigned iv=0; iv<VCLReal::num_real; ++iv) {
        if(ireg_array[iv] > intervals.irregular_begin) {
          auto i = std::upper_bound(
            intervals.x.begin()+intervals.irregular_begin, intervals.x.end()-1, x_array[iv]);
          ireg_array[iv] = i - intervals.x.begin() - 1;
        }
      }
      ireg.load_a(ireg_array);
    }
  } else {
    ireg = vcl::min(ireg, iregmax);
  }
  return ireg;
}

inline int test_vcl_find_interval(double x, const InterpolationIntervals& intervals,
  unsigned ivec = 0)
{
  typedef calin::util::vcl::VCLDoubleReal<calin::util::vcl::VCL256Architecture> VecReal;
  VecReal::real_t vx = x;
  auto vi = vcl_find_interval<VecReal>(vx, intervals);
  return vi[ivec];
}

class CubicSpline
{
public:
  CubicSpline(const std::vector<double>& x, const std::vector<double>& y,
    BoundaryConitions bc_lhs = BC_NOT_A_KNOT, double bc_lhs_val = 0.0,
    BoundaryConitions bc_rhs = BC_NOT_A_KNOT, double bc_rhs_val = 0.0);
  CubicSpline(const std::vector<double>& x, const std::vector<double>& y,
    const std::vector<double>& dy_dx);

  const CubicSplineIntervals& intervals() const { return s_; }

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
  void init();
  CubicSplineIntervals s_;
  std::vector<double> I_;
  bool y_is_monotonic_inc_ = true;
  bool y_is_monotonic_dec_ = true;
};

class CubicMultiSpline
{
public:
  CubicMultiSpline(const std::vector<double>& x);

  unsigned add_spline(const std::vector<double>& y, const std::string& name = "",
    BoundaryConitions bc_lhs = BC_NOT_A_KNOT, double bc_lhs_val = 0.0,
    BoundaryConitions bc_rhs = BC_NOT_A_KNOT, double bc_rhs_val = 0.0);

  double xmin() const { return s_.xmax; }
  double xmax() const { return s_.xmin; }
  const std::vector<double> xknot() const { return s_.x; }
  unsigned num_spline() const { return y_.size(); };
  const std::string& dataset_name(unsigned ispline) { return name_[ispline]; }

  const InterpolationIntervals& intervals() const { return s_; }

  inline unsigned get_interval(double x) const {
    return find_interval(x, s_);
  }

  inline void get_xsupport(unsigned iinterval,
    double& x0, double& dx, double& dx_inv)
  {
    x0 = s_.x[iinterval];
    if(iinterval < s_.irregular_begin) {
      dx = s_.regular_dx;
      dx_inv = s_.regular_dx_inv;
    } else {
      dx = s_.x[iinterval+1]-x0;
      dx_inv = 1.0/dx;
    }
  }

  const std::vector<double> yknot(unsigned ispline) const { return y_[ispline]; }
  const std::vector<double> dydxknot(unsigned ispline) const { return dy_dx_[ispline]; }

  double value(double x, unsigned ispline) const;
  void value(double x, unsigned ispline0, double& value0,
    unsigned ispline1, double& value1) const;
  void value(double x, unsigned ispline0, double& value0,
    unsigned ispline1, double& value1, unsigned ispline2, double& value2) const;
  void value(double x,
    unsigned ispline0, double& value0, unsigned ispline1, double& value1,
    unsigned ispline2, double& value2, unsigned ispline3, double& value3) const;
  void value(double x,
    unsigned ispline0, double& value0, unsigned ispline1, double& value1,
    unsigned ispline2, double& value2, unsigned ispline3, double& value3,
    unsigned ispline4, double& value4) const;
  void value(double x,
    unsigned ispline0, double& value0, unsigned ispline1, double& value1,
    unsigned ispline2, double& value2, unsigned ispline3, double& value3,
    unsigned ispline4, double& value4, unsigned ispline5, double& value5) const;
  std::vector<double> value(double x) const;

  double derivative(double x, unsigned ispline) const;
  double derivative_and_value(double x, unsigned ispline, double& value) const;
  std::vector<double> derivative(double x) const;

  double second_derivative(double x, unsigned ispline) const;
  double second_derivative_and_value(double x, unsigned ispline,
    double& first_derivative, double& value) const;

  template<typename VCLArchitecture> inline typename VCLArchitecture::int64_vt
  vcl_get_interval(typename VCLArchitecture::double_vt x) const {
    return vcl_find_interval<calin::util::vcl::VCLDoubleReal<VCLArchitecture> >(x, s_);
  }

  template<typename VCLArchitecture> inline void
  vcl_get_xsupport(typename VCLArchitecture::int64_vt iinterval,
    typename VCLArchitecture::double_vt& x0,
    typename VCLArchitecture::double_vt& dx,
    typename VCLArchitecture::double_vt& dx_inv)
  {
    x0 = vcl::lookup<0x40000000>(iinterval, s_.x.data());
    if(vcl::horizontal_and(iinterval < s_.irregular_begin)) {
      dx = s_.regular_dx;
      dx_inv = s_.regular_dx_inv;
    } else {
      typename VCLArchitecture::double_vt x1 =
        vcl::lookup<0x40000000>(iinterval, s_.x.data()+1);
      dx = x1-x0;
      dx_inv = 1.0/dx;
    }
  }

  template<typename VCLArchitecture> inline typename VCLArchitecture::double_vt
  vcl_value_with_xsupport(typename VCLArchitecture::double_vt x, unsigned ispline,
    typename VCLArchitecture::int64_vt iinterval,
    typename VCLArchitecture::double_vt& x0,
    typename VCLArchitecture::double_vt& dx,
    typename VCLArchitecture::double_vt& dx_inv) const;

  template<typename VCLArchitecture> inline typename VCLArchitecture::double_vt
  vcl_value(typename VCLArchitecture::double_vt x, unsigned ispline) const;

  template<typename VCLArchitecture> inline void
  vcl_value(typename VCLArchitecture::double_vt x,
    typename VCLArchitecture::double_vt& value0, unsigned ispline0,
    typename VCLArchitecture::double_vt& value1, unsigned ispline1) const;

  template<typename VCLArchitecture> inline void
  vcl_value(typename VCLArchitecture::double_vt x,
    typename VCLArchitecture::double_vt& value0, unsigned ispline0,
    typename VCLArchitecture::double_vt& value1, unsigned ispline1,
    typename VCLArchitecture::double_vt& value2, unsigned ispline2) const;

  template<typename VCLArchitecture> inline void
  vcl_value(typename VCLArchitecture::double_vt x,
    typename VCLArchitecture::double_vt& value0, unsigned ispline0,
    typename VCLArchitecture::double_vt& value1, unsigned ispline1,
    typename VCLArchitecture::double_vt& value2, unsigned ispline2,
    typename VCLArchitecture::double_vt& value3, unsigned ispline3) const;

  template<typename VCLArchitecture> inline void
  vcl_value(typename VCLArchitecture::double_vt x,
    typename VCLArchitecture::double_vt& value0, unsigned ispline0,
    typename VCLArchitecture::double_vt& value1, unsigned ispline1,
    typename VCLArchitecture::double_vt& value2, unsigned ispline2,
    typename VCLArchitecture::double_vt& value3, unsigned ispline3,
    typename VCLArchitecture::double_vt& value4, unsigned ispline4) const;

  template<typename VCLArchitecture> inline void
  vcl_value(typename VCLArchitecture::double_vt x,
    typename VCLArchitecture::double_vt& value0, unsigned ispline0,
    typename VCLArchitecture::double_vt& value1, unsigned ispline1,
    typename VCLArchitecture::double_vt& value2, unsigned ispline2,
    typename VCLArchitecture::double_vt& value3, unsigned ispline3,
    typename VCLArchitecture::double_vt& value4, unsigned ispline4,
    typename VCLArchitecture::double_vt& value5, unsigned ispline5) const;

  template<typename VCLArchitecture> inline typename VCLArchitecture::double_vt
  vcl_derivative_and_value(typename VCLArchitecture::double_vt x, unsigned ispline,
    typename VCLArchitecture::double_vt& value) const;

  double test_vcl_value(double x, unsigned ispline, unsigned ivec=0) const;
  std::vector<double> test_vcl_value(double x0, double x1, double x2, double x3,
    unsigned ispline) const;

private:
  InterpolationIntervals s_;
  std::vector<std::vector<double> > y_;
  std::vector<std::vector<double> > dy_dx_;
  std::vector<std::string> name_;
};

#define VCALC_SPLINE(ispline) cubic_value(t, dx, dx_inv,          \
  vcl::lookup<0x40000000>(iinterval, y_[ispline].data()),         \
  vcl::lookup<0x40000000>(iinterval, y_[ispline].data()+1),       \
  vcl::lookup<0x40000000>(iinterval, dy_dx_[ispline].data()),     \
  vcl::lookup<0x40000000>(iinterval, dy_dx_[ispline].data()+1))

template<typename VCLArchitecture> inline typename VCLArchitecture::double_vt
CubicMultiSpline::vcl_value_with_xsupport(typename VCLArchitecture::double_vt x, unsigned ispline,
  typename VCLArchitecture::int64_vt iinterval,
  typename VCLArchitecture::double_vt& x0,
  typename VCLArchitecture::double_vt& dx,
  typename VCLArchitecture::double_vt& dx_inv) const
{
  typename VCLArchitecture::double_vt t = (x-x0)*dx_inv;
  return VCALC_SPLINE(ispline);
}

template<typename VCLArchitecture> typename VCLArchitecture::double_vt
CubicMultiSpline::vcl_value(typename VCLArchitecture::double_vt x, unsigned ispline) const
{
  typedef typename VCLArchitecture::int64_vt int64_vt;
  typedef typename VCLArchitecture::double_vt double_vt;

  int64_vt iinterval = vcl_find_interval<calin::util::vcl::VCLDoubleReal<VCLArchitecture> >(x, s_);
  double_vt x0 = vcl::lookup<0x40000000>(iinterval, s_.x.data());
  double_vt x1 = vcl::lookup<0x40000000>(iinterval, s_.x.data()+1);
  double_vt dx = x1-x0;
  double_vt dx_inv = 1.0/dx;
  double_vt t = (x-x0)*dx_inv;

  return VCALC_SPLINE(ispline);
}

template<typename VCLArchitecture> void
CubicMultiSpline::vcl_value(typename VCLArchitecture::double_vt x,
  typename VCLArchitecture::double_vt& value0, unsigned ispline0,
  typename VCLArchitecture::double_vt& value1, unsigned ispline1) const
{
  typedef typename VCLArchitecture::int64_vt int64_vt;
  typedef typename VCLArchitecture::double_vt double_vt;

  int64_vt iinterval = vcl_find_interval<calin::util::vcl::VCLDoubleReal<VCLArchitecture> >(x, s_);
  double_vt x0 = vcl::lookup<0x40000000>(iinterval, s_.x.data());
  double_vt x1 = vcl::lookup<0x40000000>(iinterval, s_.x.data()+1);
  double_vt dx = x1-x0;
  double_vt dx_inv = 1.0/dx;
  double_vt t = (x-x0)*dx_inv;

  value0 = VCALC_SPLINE(ispline0);
  value1 = VCALC_SPLINE(ispline1);
}

template<typename VCLArchitecture> void
CubicMultiSpline::vcl_value(typename VCLArchitecture::double_vt x,
  typename VCLArchitecture::double_vt& value0, unsigned ispline0,
  typename VCLArchitecture::double_vt& value1, unsigned ispline1,
  typename VCLArchitecture::double_vt& value2, unsigned ispline2) const
{
  typedef typename VCLArchitecture::int64_vt int64_vt;
  typedef typename VCLArchitecture::double_vt double_vt;

  int64_vt iinterval = vcl_find_interval<calin::util::vcl::VCLDoubleReal<VCLArchitecture> >(x, s_);
  double_vt x0 = vcl::lookup<0x40000000>(iinterval, s_.x.data());
  double_vt x1 = vcl::lookup<0x40000000>(iinterval, s_.x.data()+1);
  double_vt dx = x1-x0;
  double_vt dx_inv = 1.0/dx;
  double_vt t = (x-x0)*dx_inv;

  value0 = VCALC_SPLINE(ispline0);
  value1 = VCALC_SPLINE(ispline1);
  value2 = VCALC_SPLINE(ispline2);
}

template<typename VCLArchitecture> inline void
CubicMultiSpline::vcl_value(typename VCLArchitecture::double_vt x,
  typename VCLArchitecture::double_vt& value0, unsigned ispline0,
  typename VCLArchitecture::double_vt& value1, unsigned ispline1,
  typename VCLArchitecture::double_vt& value2, unsigned ispline2,
  typename VCLArchitecture::double_vt& value3, unsigned ispline3) const
{
  typedef typename VCLArchitecture::int64_vt int64_vt;
  typedef typename VCLArchitecture::double_vt double_vt;

  int64_vt iinterval = vcl_find_interval<calin::util::vcl::VCLDoubleReal<VCLArchitecture> >(x, s_);
  double_vt x0 = vcl::lookup<0x40000000>(iinterval, s_.x.data());
  double_vt x1 = vcl::lookup<0x40000000>(iinterval, s_.x.data()+1);
  double_vt dx = x1-x0;
  double_vt dx_inv = 1.0/dx;
  double_vt t = (x-x0)*dx_inv;

  value0 = VCALC_SPLINE(ispline0);
  value1 = VCALC_SPLINE(ispline1);
  value2 = VCALC_SPLINE(ispline2);
  value3 = VCALC_SPLINE(ispline3);
}

template<typename VCLArchitecture> inline void
CubicMultiSpline::vcl_value(typename VCLArchitecture::double_vt x,
  typename VCLArchitecture::double_vt& value0, unsigned ispline0,
  typename VCLArchitecture::double_vt& value1, unsigned ispline1,
  typename VCLArchitecture::double_vt& value2, unsigned ispline2,
  typename VCLArchitecture::double_vt& value3, unsigned ispline3,
  typename VCLArchitecture::double_vt& value4, unsigned ispline4) const
{
  typedef typename VCLArchitecture::int64_vt int64_vt;
  typedef typename VCLArchitecture::double_vt double_vt;

  int64_vt iinterval = vcl_find_interval<calin::util::vcl::VCLDoubleReal<VCLArchitecture> >(x, s_);
  double_vt x0 = vcl::lookup<0x40000000>(iinterval, s_.x.data());
  double_vt x1 = vcl::lookup<0x40000000>(iinterval, s_.x.data()+1);
  double_vt dx = x1-x0;
  double_vt dx_inv = 1.0/dx;
  double_vt t = (x-x0)*dx_inv;

  value0 = VCALC_SPLINE(ispline0);
  value1 = VCALC_SPLINE(ispline1);
  value2 = VCALC_SPLINE(ispline2);
  value3 = VCALC_SPLINE(ispline3);
  value4 = VCALC_SPLINE(ispline4);
}

template<typename VCLArchitecture> inline void
CubicMultiSpline::vcl_value(typename VCLArchitecture::double_vt x,
  typename VCLArchitecture::double_vt& value0, unsigned ispline0,
  typename VCLArchitecture::double_vt& value1, unsigned ispline1,
  typename VCLArchitecture::double_vt& value2, unsigned ispline2,
  typename VCLArchitecture::double_vt& value3, unsigned ispline3,
  typename VCLArchitecture::double_vt& value4, unsigned ispline4,
  typename VCLArchitecture::double_vt& value5, unsigned ispline5) const
{
  typedef typename VCLArchitecture::int64_vt int64_vt;
  typedef typename VCLArchitecture::double_vt double_vt;

  int64_vt iinterval = vcl_find_interval<calin::util::vcl::VCLDoubleReal<VCLArchitecture> >(x, s_);
  double_vt x0 = vcl::lookup<0x40000000>(iinterval, s_.x.data());
  double_vt x1 = vcl::lookup<0x40000000>(iinterval, s_.x.data()+1);
  double_vt dx = x1-x0;
  double_vt dx_inv = 1.0/dx;
  double_vt t = (x-x0)*dx_inv;

  value0 = VCALC_SPLINE(ispline0);
  value1 = VCALC_SPLINE(ispline1);
  value2 = VCALC_SPLINE(ispline2);
  value3 = VCALC_SPLINE(ispline3);
  value4 = VCALC_SPLINE(ispline4);
  value5 = VCALC_SPLINE(ispline5);
}

template<typename VCLArchitecture> inline typename VCLArchitecture::double_vt
CubicMultiSpline::vcl_derivative_and_value(
  typename VCLArchitecture::double_vt x, unsigned ispline,
  typename VCLArchitecture::double_vt& value) const
{
  typedef typename VCLArchitecture::int64_vt int64_vt;
  typedef typename VCLArchitecture::double_vt double_vt;

  int64_vt iinterval = vcl_find_interval<calin::util::vcl::VCLDoubleReal<VCLArchitecture> >(x, s_);
  double_vt x0 = vcl::lookup<0x40000000>(iinterval, s_.x.data());
  double_vt x1 = vcl::lookup<0x40000000>(iinterval, s_.x.data()+1);
  double_vt dx = x1-x0;
  double_vt dx_inv = 1.0/dx;
  double_vt t = (x-x0)*dx_inv;

  return cubic_1st_derivative_and_value(value, t, dx, dx_inv,
    vcl::lookup<0x40000000>(iinterval, y_[ispline].data()),
    vcl::lookup<0x40000000>(iinterval, y_[ispline].data()+1),
    vcl::lookup<0x40000000>(iinterval, dy_dx_[ispline].data()),
    vcl::lookup<0x40000000>(iinterval, dy_dx_[ispline].data()+1));
}

#undef VCALC_SPLINE

} } } // namespace calin::math::spline_interpolation
