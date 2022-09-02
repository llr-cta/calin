/*

   calin/math/special.cpp -- Stephen Fegan -- 2015-08-06

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

#include <limits>

#include <math/special.hpp>

using namespace calin::math::special;

void calin::math::special::gaussian(double* vec, unsigned n, double mean, double sigma)
{
#if INSTRSET >= 7
  calin::math::special::gaussian_vcl<calin::util::vcl::VCLDoubleReal<calin::util::vcl::VCL256Architecture> >(vec, n, mean, sigma);
#else
  calin::math::special::gaussian_scalar(vec, n, mean, sigma);
#endif
}

void calin::math::special::two_gaussian(double* vec, unsigned n, double mean, double sigma, double split)
{
#if INSTRSET >= 7
  calin::math::special::two_gaussian_vcl<calin::util::vcl::VCLDoubleReal<calin::util::vcl::VCL256Architecture> >(vec, n, mean, sigma, split);
#else
  calin::math::special::two_gaussian_scalar(vec, n, mean, sigma, split);
#endif
}

Eigen::VectorXd calin::math::special::gaussian(unsigned n, double mean, double sigma, bool vcl)
{
  Eigen::VectorXd vec(n);
  if(vcl) {
    gaussian(vec.data(), n, mean, sigma);
  } else {
    calin::math::special::gaussian_scalar(vec.data(), n, mean, sigma);
  }
  return vec;
}

Eigen::VectorXd calin::math::special::two_gaussian(unsigned n, double mean, double sigma, double scale, bool vcl)
{
  Eigen::VectorXd vec(n);
  if(vcl) {
    two_gaussian(vec.data(), n, mean, sigma, scale);
  } else {
    calin::math::special::two_gaussian_scalar(vec.data(), n, mean, sigma, scale);
  }
  return vec;
}

unsigned calin::math::special::
solve_cubic_equation(double& x0, double&x1, double& x2, double a, double b, double c)
{
  // Roots of cubic from Numerical Recipies

  double Q1 = (SQR(a) - 3*b)/9;
  double R1 = (a*(2*SQR(a) - 9*b) + 27*c)/54;

  double Q3 = CUBE(Q1);
  double R2 = SQR(R1);

  const double a_3 = a/3;

  if(R2 < Q3) {
    double theta = std::acos(R1/std::sqrt(Q3));
    const double sqrt_q1 = std::sqrt(Q1);
    x0 = -2*sqrt_q1*std::cos((theta-2*M_PI)/3) - a_3;
    x1 = -2*sqrt_q1*std::cos(theta/3)          - a_3;
    x2 = -2*sqrt_q1*std::cos((theta+2*M_PI)/3) - a_3;
    return 3;
  } else {
    double A = -((R1<0)?-1.0:1.0) * std::cbrt(std::abs(R1) + std::sqrt(R2-Q3));
    double B = (A==0)?0:Q1/A;
    x0 = (A+B) - a_3;
    if(std::fabs(A-B) < (std::fabs(A)+std::fabs(B)) * std::numeric_limits<double>::epsilon()) {
      x1 = -0.5*(A+B) - a_3;
      x2 = std::numeric_limits<double>::quiet_NaN();
      return 2;
    } else {
      x1 = std::numeric_limits<double>::quiet_NaN();
      x2 = std::numeric_limits<double>::quiet_NaN();
      return 1;
    }
  }
}
