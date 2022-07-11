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

using namespace calin::math::spline_interpolation;

using calin::math::special::SQR;
using calin::math::special::CUBE;
using calin::math::special::QUAD;

TwoDimensionalCubicSpline::
TwoDimensionalCubicSpline(const Eigen::VectorXd& x, const Eigen::VectorXd& y,
  const Eigen::MatrixXd& z):
  sx_(make_intervals(calin::eigen_to_stdvec(x))),
  sy_(make_intervals(calin::eigen_to_stdvec(y))),
  z_(z), dz_dx_(z*0), dz_dy_(z*0), d2z_dx_dy_(z*0)
{
  if(z.cols() != x.size()) {
    throw std::range_error("Z matrix must have same number of columns as x vector");
  }
  if(z.rows() != y.size()) {
    throw std::range_error("Z matrix must have same number of rows as y vector");
  }

  for(unsigned irow=0;irow<y.size();++irow) {
    Eigen::VectorXd zrow = z_.row(irow);
    CubicSpline xspline(sx_.x, calin::eigen_to_stdvec(zrow));
    dz_dx_.row(irow) = xspline.dydxknot_as_eigen();
  }

  for(unsigned icol=0;icol<x.size();++icol) {
    Eigen::VectorXd zcol = z_.col(icol);
    CubicSpline yspline(sy_.x, calin::eigen_to_stdvec(zcol));
    dz_dy_.col(icol) = yspline.dydxknot_as_eigen();
  }

  for(unsigned icol=0;icol<x.size();++icol) {
    Eigen::VectorXd dz_dx_col = dz_dx_.col(icol);
    CubicSpline dz_dx_spline(sy_.x, calin::eigen_to_stdvec(dz_dx_col));
    d2z_dx_dy_.col(icol) = dz_dx_spline.dydxknot_as_eigen();
  }
}

double TwoDimensionalCubicSpline::value(double x, double y) const
{
  double dx;
  double dx_inv;
  unsigned ix = find_interval(x, sx_, dx, dx_inv);

  double dy;
  double dy_inv;
  unsigned iy = find_interval(y, sy_, dy, dy_inv);

  double tx = (x-sx_.x[ix])*dx_inv;
  double ty = (y-sy_.x[iy])*dy_inv;

  double z0 = cubic_value(ty, dy, dy_inv, z_(iy,ix), z_(iy+1,ix), dz_dy_(iy,ix), dz_dy_(iy+1,ix));
  double z1 = cubic_value(ty, dy, dy_inv, z_(iy,ix+1), z_(iy+1,ix+1), dz_dy_(iy,ix+1), dz_dy_(iy+1,ix+1));

  double dz_dx_0 = cubic_value(ty, dy, dy_inv, dz_dx_(iy,ix), dz_dx_(iy+1,ix), d2z_dx_dy_(iy,ix), d2z_dx_dy_(iy+1,ix));
  double dz_dx_1 = cubic_value(ty, dy, dy_inv, dz_dx_(iy,ix+1), dz_dx_(iy+1,ix+1), d2z_dx_dy_(iy,ix+1), d2z_dx_dy_(iy+1,ix+1));

  return cubic_value(tx, dx, dx_inv, z0, z1, dz_dx_0, dz_dx_1);
}
