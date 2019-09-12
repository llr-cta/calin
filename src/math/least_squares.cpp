/*

   calin/math/least_squares.cpp -- Stephen Fegan -- 2019-03-07

   Implemantion of the clasic polynomial fit by least-squares

   Copyright 2019, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <exception>

#include <math/least_squares.hpp>
#include <util/log.hpp>

using namespace calin::math::least_squares;
using namespace calin::util::log;

Eigen::VectorXd
calin::math::least_squares::polyfit(const Eigen::VectorXd& x, const Eigen::VectorXd& y, unsigned norder)
{
  if(x.size() != y.size()) {
    throw std::runtime_error("polyfit: x and y arrays must have same size");
  }
  if(x.size() < norder+1) {
    throw std::runtime_error("polyfit: x and y arrays must have at least norder+1 points");
  }

  // See https://eigen.tuxfamily.org/dox/group__LeastSquares.html
  Eigen::MatrixXd X(x.size(), norder+1);
  for(unsigned ix = 0;ix<x.size();ix++) {
    double xi =  x[ix];
    double xn = 1;
    X(ix, 0) = xn;
    for(unsigned iorder=0;iorder<norder;++iorder) {
      xn *= xi;
      X(ix, iorder+1) = xn;
    }
  }
  return X.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(y);
}
