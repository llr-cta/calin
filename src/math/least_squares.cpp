/*

   calin/math/least_squares.cpp -- Stephen Fegan -- 2019-03-07

   Implemantion of the clasic polynomial fit by least-squares

   Copyright 2019, Stephen Fegan <sfegan@llr.in2p3.fr>
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

  Eigen::VectorXd xxs(2*norder+1);
  Eigen::VectorXd yxs(norder+1);

  xxs.setZero();
  yxs.setZero();

  for(unsigned ixy = 0;ixy<x.size();ixy++) {
    double xx = 1;
    xxs(0) += xx;
    for(unsigned iorder=0;iorder<2*norder;++iorder) {
      xx *= x[ixy];
      xxs(iorder+1) += xx;
    }

    double yx = y[ixy];
    yxs(0) += yx;
    for(unsigned iorder=0;iorder<norder;++iorder) {
      yx *= x[ixy];
      yxs(iorder+1) += yx;
    }
  }

  Eigen::MatrixXd XXs(norder+1,norder+1);
  for(unsigned iorder=0;iorder<norder+1;++iorder) {
    for(unsigned jorder=0;jorder<norder+1;++jorder) {
      XXs(jorder,iorder) = xxs(iorder+jorder);
    }
  }

  //LOG(INFO) << XXs << '\n' << yxs;

  return XXs.llt().solve(yxs);
}
