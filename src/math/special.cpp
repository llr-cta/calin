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
