/*

   calin/math/spline_interpolation.hpp -- Stephen Fegan -- 2018-12-17

   Spline interpolation functions.

   Some of this code from my Fermi tools code "RegularSpline.hpp",
   These portions are copyright Stephen Fegan, 2010.

   Copyright 2018, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <stdexcept>
#include <iostream>

#include <math/spline_interpolation.hpp>
#include <math/special.hpp>
#include <util/log.hpp>

using namespace calin::util::log;
using namespace calin::math::spline_interpolation;

using calin::math::special::SQR;
using calin::math::special::CUBE;
using calin::math::special::QUAD;

TwoDimensionalCubicSpline::
TwoDimensionalCubicSpline(const Eigen::VectorXd& x, const Eigen::VectorXd& y,
  const Eigen::MatrixXd& z):
  sx_(make_intervals(calin::eigen_to_stdvec(x))),
  sy_(make_intervals(calin::eigen_to_stdvec(y))),
  u_(z), p_(z*0), q_(z*0), r_(z*0)
{
  // Implement algorithm from "Two Dimensional Spline Intepolation Algorithms"
  // H. Spath, 1995, Chapter

  if(z.cols() != x.size()) {
    throw std::range_error("Z matrix must have same number of columns as x vector");
  }
  if(z.rows() != y.size()) {
    throw std::range_error("Z matrix must have same number of rows as y vector");
  }

  // Solve equation 4.9 of Spath, 1995
  for(unsigned irow=0;irow<y.size();++irow) {
    p_.row(irow) = generate_cubic_spline_interpolation_eigen(x, u_.row(irow));
  }

  // Solve equation 4.10 of Spath, 1995
  for(unsigned icol=0;icol<x.size();++icol) {
    q_.col(icol) = generate_cubic_spline_interpolation_eigen(y, u_.col(icol));
  }

  // Solve equation 4.13 of Spath, 1995
  r_.row(0) = generate_cubic_spline_interpolation_eigen(x, q_.row(0));
  r_.row(y.size()-1) = generate_cubic_spline_interpolation_eigen(x, q_.row(y.size()-1));

  // Solve equation 4.14 of Spath, 1995
  for(unsigned icol=0;icol<x.size();++icol) {
    r_.col(icol) = generate_cubic_spline_interpolation_eigen(y, p_.col(icol),
      BC_CLAMPED_SLOPE, r_(0,icol), BC_CLAMPED_SLOPE, r_(y.size()-1,icol));
  }
}

double TwoDimensionalCubicSpline::value_old(double x, double y) const
{
  double dx;
  double dx_inv;
  unsigned ix = find_interval(x, sx_, dx, dx_inv);

  double dx_inv2 = dx_inv*dx_inv;
  double dx_inv3 = dx_inv2*dx_inv;
  Eigen::Matrix4d Vx_inv;
  Vx_inv << 1,                  0,          0,       0,
            0,                  1,          0,       0,
            -3*dx_inv2, -2*dx_inv,  3*dx_inv2, -dx_inv,
            2*dx_inv3,    dx_inv2, -2*dx_inv3, dx_inv2;

  // LOG(INFO) << "Vx_inv\n" << Vx_inv << '\n';

  double dy;
  double dy_inv;
  unsigned iy = find_interval(y, sy_, dy, dy_inv);

  double dy_inv2 = dy_inv*dy_inv;
  double dy_inv3 = dy_inv2*dy_inv;
  Eigen::Matrix4d Vy_inv;
  Vy_inv << 1,                  0,          0,       0,
            0,                  1,          0,       0,
            -3*dy_inv2, -2*dy_inv,  3*dy_inv2, -dy_inv,
            2*dy_inv3,    dy_inv2, -2*dy_inv3, dy_inv2;

  // LOG(INFO) << "Vy_inv\n" << Vy_inv << '\n';

  Eigen::Matrix4d C;
  C << u_(iy,ix),   q_(iy,ix),   u_(iy+1,ix),   q_(iy+1,ix),
       p_(iy,ix),   r_(iy,ix),   p_(iy+1,ix),   r_(iy+1,ix),
       u_(iy,ix+1), q_(iy,ix+1), u_(iy+1,ix+1), q_(iy+1,ix+1),
       p_(iy,ix+1), r_(iy,ix+1), p_(iy+1,ix+1), r_(iy+1,ix+1);

 // LOG(INFO) << "C\n" << C << '\n';

  x -= sx_.x[ix];
  y -= sy_.x[iy];

  double x2 = x*x;
  double x3 = x2*x;
  Eigen::Vector4d gx;
  gx << 1, x, x2, x3;

  // LOG(INFO) << "gx\n" << gx << '\n';

  double y2 = y*y;
  double y3 = y2*y;
  Eigen::Vector4d gy;
  gy << 1, y, y2, y3;

  // LOG(INFO) << "gy\n" << gx << '\n';

  // LOG(INFO) << "Vy_inv * gy\n" << Vy_inv.transpose() * gy << '\n';


  return gx.transpose() * Vx_inv * C * Vy_inv.transpose() * gy;
}

double TwoDimensionalCubicSpline::value(double x, double y) const
{
  double dx;
  double dx_inv;
  unsigned ix = find_interval(x, sx_, dx, dx_inv);

  x = (x - sx_.x[ix])*dx_inv;
  double xn = x;

  Eigen::Vector4d gx;
  gx(0) = 1;
  gx(1) = xn;
  xn *= x;
  gx(1) -= 2*xn;
  gx(2) = 3*xn;
  gx(3) = -xn;
  xn *= x;
  gx(1) += xn;
  gx(2) -= 2*xn;
  gx(3) += xn;
  gx(0) -= gx(2);
  gx(1) *= dx;
  gx(3) *= dx;

  double dy;
  double dy_inv;
  unsigned iy = find_interval(y, sy_, dy, dy_inv);

  y = (y - sy_.x[iy])*dy_inv;
  double yn = y;

  Eigen::Vector4d gy;
  gy(0) = 1;
  gy(1) = yn;
  yn *= y;
  gy(1) -= 2*yn;
  gy(2) = 3*yn;
  gy(3) = -yn;
  yn *= y;
  gy(1) += yn;
  gy(2) -= 2*yn;
  gy(3) += yn;
  gy(0) -= gy(2);
  gy(1) *= dy;
  gy(3) *= dy;

  Eigen::Matrix4d C;
  C << u_(iy,ix),   q_(iy,ix),   u_(iy+1,ix),   q_(iy+1,ix),
       p_(iy,ix),   r_(iy,ix),   p_(iy+1,ix),   r_(iy+1,ix),
       u_(iy,ix+1), q_(iy,ix+1), u_(iy+1,ix+1), q_(iy+1,ix+1),
       p_(iy,ix+1), r_(iy,ix+1), p_(iy+1,ix+1), r_(iy+1,ix+1);

  return gx.transpose() * C * gy;
}

double TwoDimensionalCubicSpline::test_value(unsigned n, double x, double y) const
{
  volatile double z;
  while(n--) {
    z = value(x,y);
  }
  return z;
}

double TwoDimensionalCubicSpline::test_value_old(unsigned n, double x, double y) const
{
  volatile double z;
  while(n--) {
    z = value_old(x,y);
  }
  return z;
}

double TwoDimensionalCubicSpline::test_vcl_value_128(double x, double y) const
{
  return vcl_value<calin::util::vcl::VCL128Architecture>(x, y)[0];
}

double TwoDimensionalCubicSpline::test_vcl_value_128(unsigned n, double x, double y) const
{
  typename calin::util::vcl::VCL128Architecture::double_vt z;
  while(n--) {
    z = vcl_value<calin::util::vcl::VCL128Architecture>(x, y);
  }
  return z[0];
}

double TwoDimensionalCubicSpline::test_vcl_value_256(double x, double y) const
{
  return vcl_value<calin::util::vcl::VCL256Architecture>(x, y)[0];
}

double TwoDimensionalCubicSpline::test_vcl_value_256(unsigned n, double x, double y) const
{
  typename calin::util::vcl::VCL256Architecture::double_vt z;
  while(n--) {
    z = vcl_value<calin::util::vcl::VCL256Architecture>(x, y);
  }
  return z[0];
}

double TwoDimensionalCubicSpline::test_vcl_value_512(double x, double y) const
{
  return vcl_value<calin::util::vcl::VCL512Architecture>(x, y)[0];
}

double TwoDimensionalCubicSpline::test_vcl_value_512(unsigned n, double x, double y) const
{
  typename calin::util::vcl::VCL512Architecture::double_vt z;
  while(n--) {
    z = vcl_value<calin::util::vcl::VCL512Architecture>(x, y);
  }
  return z[0];
}
