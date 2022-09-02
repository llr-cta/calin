/*

   calin/math/b_spline.cpp -- Stephen Fegan -- 2018-10-23

   Functions for B-spline evaluation

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

#include <cmath>
#include <string>
#include <vector>
#include <stdexcept>

#include <math/b_spline.hpp>

namespace {

unsigned knot_index(double& x, double x0, double one_over_dx, unsigned n, unsigned k)
{
  if(x<x0)throw std::out_of_range("knot_index: x-value below start of spline knots: "
    + std::to_string(x) + " < " + std::to_string(x0));
  double t = (x-x0)*one_over_dx;
  double t_round = std::floor(t);
  unsigned ix = unsigned(t_round) + k;
  if(ix>=n)throw std::out_of_range("knot_index: x-value beyond end of spline knots: "
    + std::to_string(x) + " >= " + std::to_string(x0 + (n-k)/one_over_dx));
  x = t-t_round;
  return ix;
}

}

double calin::math::b_spline::
value0(double x, const double* c, unsigned nc, double x0, double one_over_dx)
{
  unsigned ix = knot_index(x, x0, one_over_dx, nc, 0);
  return interval_value0(x, c[ix]);
}

double calin::math::b_spline::
value1(double x, const double* c, unsigned nc, double x0, double one_over_dx)
{
  unsigned ix = knot_index(x, x0, one_over_dx, nc, 1);
  return interval_value1(x, c[ix], c[ix-1]);
}

double calin::math::b_spline::
value2(double x, const double* c, unsigned nc, double x0, double one_over_dx)
{
  unsigned ix = knot_index(x, x0, one_over_dx, nc, 2);
  return interval_value2(x, c[ix], c[ix-1], c[ix-2]);
}

double calin::math::b_spline::
value3(double x, const double* c, unsigned nc, double x0, double one_over_dx)
{
  unsigned ix = knot_index(x, x0, one_over_dx, nc, 3);
  return interval_value3(x, c[ix], c[ix-1], c[ix-2], c[ix-3]);
}

double calin::math::b_spline::
value0(double x, const std::vector<double>& c, double x0, double one_over_dx)
{
  return value0(x, c.data(), c.size(), x0, one_over_dx);
}

double calin::math::b_spline::
value1(double x, const std::vector<double>& c, double x0, double one_over_dx)
{
  return value1(x, c.data(), c.size(), x0, one_over_dx);
}

double calin::math::b_spline::
value2(double x, const std::vector<double>& c, double x0, double one_over_dx)
{
  return value2(x, c.data(), c.size(), x0, one_over_dx);
}

double calin::math::b_spline::
value3(double x, const std::vector<double>& c, double x0, double one_over_dx)
{
  return value3(x, c.data(), c.size(), x0, one_over_dx);
}

double calin::math::b_spline::
integral0(double x, const double* c, unsigned nc, double x0, double one_over_dx)
{
  unsigned ix = knot_index(x, x0, one_over_dx, nc, 0);
  return interval_integral0(x, c[ix]);
}

double calin::math::b_spline::
integral1(double x, const double* c, unsigned nc, double x0, double one_over_dx)
{
  unsigned ix = knot_index(x, x0, one_over_dx, nc, 1);
  return interval_integral1(x, c[ix], c[ix-1]);
}

double calin::math::b_spline::
integral2(double x, const double* c, unsigned nc, double x0, double one_over_dx)
{
  unsigned ix = knot_index(x, x0, one_over_dx, nc, 2);
  return interval_integral2(x, c[ix], c[ix-1], c[ix-2]);
}

double calin::math::b_spline::
integral3(double x, const double* c, unsigned nc, double x0, double one_over_dx)
{
  unsigned ix = knot_index(x, x0, one_over_dx, nc, 3);
  return interval_integral3(x, c[ix], c[ix-1], c[ix-2], c[ix-3]);
}

double calin::math::b_spline::
integral0(double x, const std::vector<double>& c, double x0, double one_over_dx)
{
  return integral0(x, c.data(), c.size(), x0, one_over_dx);
}

double calin::math::b_spline::
integral1(double x, const std::vector<double>& c, double x0, double one_over_dx)
{
  return integral1(x, c.data(), c.size(), x0, one_over_dx);
}

double calin::math::b_spline::
integral2(double x, const std::vector<double>& c, double x0, double one_over_dx)
{
  return integral2(x, c.data(), c.size(), x0, one_over_dx);
}

double calin::math::b_spline::
integral3(double x, const std::vector<double>& c, double x0, double one_over_dx)
{
  return integral3(x, c.data(), c.size(), x0, one_over_dx);
}


double calin::math::b_spline::
derivative0(double x, const double* c, unsigned nc, double x0, double one_over_dx)
{
  unsigned ix = knot_index(x, x0, one_over_dx, nc, 0);
  return interval_derivative0(x, c[ix]);
}

double calin::math::b_spline::
derivative1(double x, const double* c, unsigned nc, double x0, double one_over_dx)
{
  unsigned ix = knot_index(x, x0, one_over_dx, nc, 1);
  return interval_derivative1(x, c[ix], c[ix-1]);
}

double calin::math::b_spline::
derivative2(double x, const double* c, unsigned nc, double x0, double one_over_dx)
{
  unsigned ix = knot_index(x, x0, one_over_dx, nc, 2);
  return interval_derivative2(x, c[ix], c[ix-1], c[ix-2]);
}

double calin::math::b_spline::
derivative3(double x, const double* c, unsigned nc, double x0, double one_over_dx)
{
  unsigned ix = knot_index(x, x0, one_over_dx, nc, 3);
  return interval_derivative3(x, c[ix], c[ix-1], c[ix-2], c[ix-3]);
}

double calin::math::b_spline::
derivative0(double x, const std::vector<double>& c, double x0, double one_over_dx)
{
  return derivative0(x, c.data(), c.size(), x0, one_over_dx);
}

double calin::math::b_spline::
derivative1(double x, const std::vector<double>& c, double x0, double one_over_dx)
{
  return derivative1(x, c.data(), c.size(), x0, one_over_dx);
}

double calin::math::b_spline::
derivative2(double x, const std::vector<double>& c, double x0, double one_over_dx)
{
  return derivative2(x, c.data(), c.size(), x0, one_over_dx);
}

double calin::math::b_spline::
derivative3(double x, const std::vector<double>& c, double x0, double one_over_dx)
{
  return derivative3(x, c.data(), c.size(), x0, one_over_dx);
}
