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

void I64LinearRegressionAccumulator::accumulate(int64_t x, int64_t y)
{
  if(not zero_set_) {
    x0_ = x;
    y0_ = y;
    W_ = 1;
    zero_set_ = true;
  } else {
    x -= x0_;
    y -= y0_;
    W_ += 1;
    X_ += x;
    Y_ += y;
    XX_ += x*x;
    XY_ += x*y;
    YY_ += y*y;
  }
}

void I64LinearRegressionAccumulator::integrate_into(I64LinearRegressionAccumulator& o)
{
  if(not o.zero_set_) {
    o = *this;
  } else {
    int64_t dx0 = x0_ - o.x0_;
    int64_t dy0 = y0_ - o.y0_;
    o.W_  += W_;
    o.X_  += X_ + W_*dx0;
    o.Y_  += Y_ + W_*dy0;
    o.XX_ += XX_ + 2*X_*dx0 + W_*dx0*dx0;
    o.XY_ += XY_ + X_*dy0 + dx0*Y_ + W_*dx0*dy0;
    o.YY_ += YY_ + 2*Y_*dy0 + W_*dy0*dy0;
  }
}

void I64LinearRegressionAccumulator::fit_parameters(double& a, double& b) const
{
  double D2 = 0;
  fit_parameters_and_d2(a, b, D2);
}

namespace {
  void mul256(__uint128_t a, __uint128_t b, __uint128_t& mh, __uint128_t& ml) {
    __uint128_t a_h = a>>64;
    __uint128_t a_l = a & ((__uint128_t(1)<<64)-1);
    __uint128_t b_h = b>>64;
    __uint128_t b_l = b & ((__uint128_t(1)<<64)-1);

    ml = b_l * a_l;
    mh = 0;

    __uint128_t m = a_l * b_h;
    mh += m>>64;
    m <<= 64;
    __uint128_t s = m + ml;
    if(s<m or s<ml)mh += 1;
    ml = s;

    m = a_h * b_l;
    mh += m>>64;
    m <<= 64;
    s = m + ml;
    if(s<m or s<ml)mh += 1;
    ml = s;

    mh += a_h * b_h;
  }
}

void I64LinearRegressionAccumulator::fit_parameters_and_d2(double& a, double& b, double& D2) const
{
  calin::util::log::LOG(calin::util::log::INFO)
    << "XX: " << double(XX_) << " X: " << double(X_) << " W: " << double(W_) << ' '
    << " XY: " << double(XY_) << " Y: " << double(Y_) << " YY: " << double(YY_);

  __int128_t W = W_;
  __int128_t X = X_;
  __int128_t Y = Y_;
  __int128_t XX = XX_;
  __int128_t XY = XY_;
  __int128_t YY = YY_;

  __int128_t sxy = XY*W - X*Y;
  __int128_t sxx = XX*W - X*X;
  __int128_t syy = YY*W - Y*Y;

  calin::util::log::LOG(calin::util::log::INFO)
    << "Sxx: " << double(sxx) << " Sxy: " << double(sxy) << " Syy: " << double(syy);

  a = double(sxy)/double(sxx);

  calin::util::log::LOG(calin::util::log::INFO)
    << "a_num: " << double(sxy) << " a_den: " << double(sxx) << " a: " << a;

  __int128_t b_num = Y_*XX_ - X_*XY_;

  b = double(b_num)/double(sxx);

  calin::util::log::LOG(calin::util::log::INFO)
    << "b_num: " << double(b_num) << " b: " << b;

  // This is where the fun starts !!

  __uint128_t sxy_p = (sxy<0) ? (-sxy) : sxy;

  __uint128_t syy_sxx_h;
  __uint128_t syy_sxx_l;
  __uint128_t sxy_sxy_h;
  __uint128_t sxy_sxy_l;

  mul256(sxx, syy, syy_sxx_h, syy_sxx_l);
  mul256(sxy_p, sxy_p, sxy_sxy_h, sxy_sxy_l);

  calin::util::log::LOG(calin::util::log::INFO)
    << "syy_sxx_h: " << double(syy_sxx_h) << " syy_sxx_l: " << double(syy_sxx_l)
    << " sxy_sxy_h: " << double(sxy_sxy_h) << " sxy_sxy_l: " << double(sxy_sxy_l);

  __uint128_t num_h = syy_sxx_h - sxy_sxy_h;
  __uint128_t num_l;
  if(syy_sxx_l >= sxy_sxy_l) {
    num_l = syy_sxx_l - sxy_sxy_l;
  } else {
    num_l = ~(sxy_sxy_l - syy_sxx_l);
    num_h -= 1;
  }

  calin::util::log::LOG(calin::util::log::INFO)
    << "num_h: " << double(num_h) << " num_l: " << double(num_l);

  D2 = (ldexp(double(num_h),128) + double(num_l))/double(sxx * W);
}


void KahanLinearRegressionAccumulator::accumulate(double x, double y)
{
  if(not zero_set) {
    x0 = x;
    y0 = y;
    W.accumulate(1.0);
    zero_set = true;
  } else {
    x -= x0;
    y -= y0;
    W.accumulate(1.0);
    X.accumulate(x);
    Y.accumulate(y);
    XX.accumulate(x*x);
    XY.accumulate(x*y);
    YY.accumulate(y*y);
  }
}

void KahanLinearRegressionAccumulator::integrate_into(KahanLinearRegressionAccumulator& o)
{
  if(not o.zero_set) {
    o.x0 = x0;
    o.y0 = y0;
    o.W  = W;
    o.X  = X;
    o.Y  = Y;
    o.XX = XX;
    o.XY = XY;
    o.YY = YY;
    o.zero_set = zero_set;
  } else {
    double dx0 = x0 - o.x0;
    double dy0 = y0 - o.y0;
    o.W.accumulate_from(W);
    o.X.accumulate_from(X);
    o.X.accumulate(W.total() * dx0);
    o.Y.accumulate_from(Y);
    o.Y.accumulate(W.total() * dy0);
    o.XX.accumulate_from(XX);
    o.XX.accumulate(2 * X.total() * dx0);
    o.XX.accumulate(W.total() * dx0*dx0);
    o.XY.accumulate_from(XY);
    o.XY.accumulate(X.total() * dy0);
    o.XY.accumulate(dx0 * Y.total());
    o.XY.accumulate(W.total() * dx0*dy0);
    o.YY.accumulate_from(YY);
    o.YY.accumulate(2 * Y.total() * dy0);
    o.YY.accumulate(W.total() * dy0*dy0);
  }
}

void KahanLinearRegressionAccumulator::fit_parameters_and_d2(double& a, double& b, double& D2)
{
  calin::util::log::LOG(calin::util::log::INFO)
    << "XX: " << XX.total() << " X: " << X.total() << " W: " << W.total() << ' '
    << " XY: " << XY.total() << " Y: " << Y.total() << " YY: " << YY.total();

  calin::math::accumulator::KahanAccumulator a_num { XY };
  a_num.accumulate(-X.total()*Y.total()/W.total());

  calin::math::accumulator::KahanAccumulator a_den { XX };
  a_den.accumulate(-X.total()*X.total()/W.total());

  double a_dir = a_num.total()/a_den.total();
  double b_dir = (Y.total() - a_dir*X.total())/W.total();

  calin::util::log::LOG(calin::util::log::INFO)
    << "Direct A : " << a_dir;

  calin::util::log::LOG(calin::util::log::INFO)
    << "Direct B : " << b_dir;

  calin::math::accumulator::KahanAccumulator d2 { YY };
  d2.accumulate_from_with_scaling(XX, a_dir*a_dir);
  d2.accumulate_from_with_scaling(W, b_dir*b_dir);
  d2.accumulate_from_with_scaling(XY, -2*a_dir);
  d2.accumulate_from_with_scaling(X, 2*a_dir*b_dir);
  d2.accumulate_from_with_scaling(Y, -2*b_dir);

  calin::util::log::LOG(calin::util::log::INFO)
    << "Direct D2 : " << d2.total();

  Eigen::Matrix2d m;
  m << XX.total(), X.total(), X.total(), W.total();
  Eigen::Vector2d v;
  v << XY.total(), Y.total();
  Eigen::Vector2d ab = m.fullPivHouseholderQr().solve(v);
  a = ab(0);
  b = ab(1);

  D2 = YY.total() + a*a*XX.total() + b*b*W.total() - 2*a*XY.total() - 2*b*Y.total() + 2*a*b*X.total();
  // b = b - a*x0 + y0;
}

void KahanLinearRegressionAccumulator::fit_parameters(double& a, double& b)
{
  double D2;
  fit_parameters_and_d2(a, b, D2);
}


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
