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
#include <math/special.hpp>

using namespace calin::math::least_squares;
using namespace calin::util::log;

using namespace calin::util::log;
using calin::math::special::SQR;

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
    __int128_t x128 = x;
    __int128_t y128 = y;
    W_ += 1;
    X_ += x128;
    Y_ += y128;
    XX_ += x128*x128;
    XY_ += x128*y128;
    YY_ += y128*y128;
  }
}

void I64LinearRegressionAccumulator::integrate_into(I64LinearRegressionAccumulator& o) const
{
  if(not o.zero_set_) {
    o = *this;
  } else {
    __int128_t dx0 = x0_ - o.x0_;
    __int128_t dy0 = y0_ - o.y0_;
    o.W_  += W_;
    o.X_  += X_                    + W_*dx0;
    o.Y_  += Y_                    + W_*dy0;
    o.XX_ += XX_ + 2*X_*dx0        + W_*dx0*dx0;
    o.XY_ += XY_ + X_*dy0 + dx0*Y_ + W_*dx0*dy0;
    o.YY_ += YY_ + 2*Y_*dy0        + W_*dy0*dy0;
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

void I64LinearRegressionAccumulator::rebalance()
{
  shift_origin(x0_ + X_/W_, y0_ + Y_/W_);
}

void I64LinearRegressionAccumulator::shift_origin(int64_t x0, int64_t y0)
{
  if(not zero_set_) {
    x0_ = x0;
    y0_ = y0;
    zero_set_ = true;
  } else {
    __int128_t dx0 = x0_ - x0;
    __int128_t dy0 = y0_ - y0;
    x0_  = x0;
    y0_  = y0;
    // Watch out for the order we do these shifts.. finish with X_ and Y_
    XX_ += 2*X_*dx0 + W_*dx0*dx0;
    XY_ += X_*dy0 + dx0*Y_ + W_*dx0*dy0;
    YY_ += 2*Y_*dy0 + W_*dy0*dy0;
    X_  += W_*dx0;
    Y_  += W_*dy0;
  }
}

double I64LinearRegressionAccumulator::mean_x() const
{
  return double(X_)/double(W_);
}

double I64LinearRegressionAccumulator::mean_y() const
{
  return double(Y_)/double(W_);
}

void I64LinearRegressionAccumulator::moments(
  double& entries, double& mean_x, double& mean_y,
  double& sigma_xx, double& sigma_yy, double& sigma_xy) const
{
  __int128_t sxy = XY_*W_ - X_*Y_;
  __int128_t sxx = XX_*W_ - X_*X_;
  __int128_t syy = YY_*W_ - Y_*Y_;

  entries = double(W_);
  mean_x = double(X_)/entries;
  mean_y = double(Y_)/entries;
  sigma_xx = double(sxx)/SQR(entries);
  sigma_yy = double(syy)/SQR(entries);
  sigma_xy = double(sxy)/SQR(entries);
}

void I64LinearRegressionAccumulator::fit_parameters_and_d2(double& a, double& b, double& D2) const
{
#ifdef DEBUG_I64LinearRegressionAccumulator
  calin::util::log::LOG(calin::util::log::INFO)
    << "XX: " << double(XX_) << " X: " << double(X_) << " W: " << double(W_) << ' '
    << " XY: " << double(XY_) << " Y: " << double(Y_) << " YY: " << double(YY_);
#endif

  // Possible limitation on accumulation here - overflow of X*X and XX*W etc
  __int128_t sxy = XY_*W_ - X_*Y_;
  __int128_t sxx = XX_*W_ - X_*X_;
  __int128_t syy = YY_*W_ - Y_*Y_;

#ifdef DEBUG_I64LinearRegressionAccumulator
  calin::util::log::LOG(calin::util::log::INFO)
    << "Sxx: " << double(sxx) << " Sxy: " << double(sxy) << " Syy: " << double(syy);
#endif

  a = double(sxy)/double(sxx);

#ifdef DEBUG_I64LinearRegressionAccumulator
  calin::util::log::LOG(calin::util::log::INFO)
    << "a_num: " << double(sxy) << " a_den: " << double(sxx) << " a: " << a;
#endif

  __int128_t b_num = Y_*XX_ - X_*XY_;

  b = double(b_num)/double(sxx);

#ifdef DEBUG_I64LinearRegressionAccumulator
  calin::util::log::LOG(calin::util::log::INFO)
    << "b_num: " << double(b_num) << " b: " << b;
#endif

  // This is where the fun starts - two 128bit multiplications into 256bits and
  // then subtraction of terms to maintain full precision in SYY*SXX-SXY*SXY

  __uint128_t sxy_p = (sxy<0) ? (-sxy) : sxy;

  __uint128_t syy_sxx_h;
  __uint128_t syy_sxx_l;
  __uint128_t sxy_sxy_h;
  __uint128_t sxy_sxy_l;

  mul256(sxx, syy, syy_sxx_h, syy_sxx_l);
  mul256(sxy_p, sxy_p, sxy_sxy_h, sxy_sxy_l);

#ifdef DEBUG_I64LinearRegressionAccumulator
  calin::util::log::LOG(calin::util::log::INFO)
    << "syy_sxx_h: " << double(syy_sxx_h) << " syy_sxx_l: " << double(syy_sxx_l)
    << " sxy_sxy_h: " << double(sxy_sxy_h) << " sxy_sxy_l: " << double(sxy_sxy_l);
#endif

  __uint128_t num_h = syy_sxx_h - sxy_sxy_h;
  __uint128_t num_l;
  if(syy_sxx_l >= sxy_sxy_l) {
    num_l = syy_sxx_l - sxy_sxy_l;
  } else {
    num_l = ~(sxy_sxy_l - syy_sxx_l - 1); // Two's complement subtraction
    num_h -= 1;                           // with carry from high-bits
  }

#ifdef DEBUG_I64LinearRegressionAccumulator
  calin::util::log::LOG(calin::util::log::INFO)
    << "num_h: " << double(num_h) << " num_l: " << double(num_l);
#endif

  D2 = (ldexp(double(num_h),128) + double(num_l))/double(sxx * W_);
}

void I64LinearRegressionAccumulator::dump_to_log() const
{
  LOG(INFO)
    << "x0: " << x0_ << '\n'
    << "y0: " << y0_ << '\n'
    << "W: " << int64_t(W_ >> 64) << ' ' << uint64_t(W_ & ((__int128_t(1)<<64)-1)) << '\n'
    << "X: " << int64_t(X_ >> 64) << ' ' << uint64_t(X_ & ((__int128_t(1)<<64)-1)) << '\n'
    << "Y: " << int64_t(Y_ >> 64) << ' ' << uint64_t(Y_ & ((__int128_t(1)<<64)-1)) << '\n'
    << "XX: " << int64_t(XX_ >> 64) << ' ' << uint64_t(XX_ & ((__int128_t(1)<<64)-1)) << '\n'
    << "XY: " << int64_t(XY_ >> 64) << ' ' << uint64_t(XY_ & ((__int128_t(1)<<64)-1)) << '\n'
    << "YY: " << int64_t(YY_ >> 64) << ' ' << uint64_t(YY_ & ((__int128_t(1)<<64)-1)) << '\n';
}

void I64LinearRegressionAccumulatorIgnoringFirstDatum::
accumulate(int64_t x, int64_t y)
{
  if(has_x0_) {
    if(x >= x0_) {
      accumulator_.accumulate(x, y);
    } else {
      accumulator_.accumulate(x0_, y0_);
      x0_ = x;
      y0_ = y;
    }
  } else {
    x0_ = x;
    y0_ = y;
    has_x0_ = true;
  }
}

void I64LinearRegressionAccumulatorIgnoringFirstDatum::
integrate_into(I64LinearRegressionAccumulatorIgnoringFirstDatum& o) const
{
  accumulator_.integrate_into(o.accumulator_);

  if(has_x0_) {
    if(o.has_x0_) {
      if(x0_ >= o.x0_) {
        o.accumulator_.accumulate(x0_, y0_);
      } else {
        o.accumulator_.accumulate(o.x0_, o.y0_);
        o.x0_ = x0_;
        o.y0_ = y0_;
      }
    } else {
      o.x0_ = x0_;
      o.y0_ = y0_;
      o.has_x0_ = true;
    }
  }
}

void I64LinearRegressionAccumulatorIgnoringFirstDatum::integrate_first_event()
{
  if(has_x0_) {
    accumulator_.accumulate(x0_, y0_);
    has_x0_ = false;
  }
}

void I64LinearRegressionAccumulatorIgnoringFirstDatum::rebalance()
{
  accumulator_.rebalance();
}

void KahanLinearRegressionAccumulator::accumulate(double x, double y)
{
  if(not zero_set_) {
    x0_ = x;
    y0_ = y;
    W_.accumulate(1.0);
    zero_set_ = true;
  } else {
    x -= x0_;
    y -= y0_;
    W_.accumulate(1.0);
    X_.accumulate(x);
    Y_.accumulate(y);
    XX_.accumulate(x*x);
    XY_.accumulate(x*y);
    YY_.accumulate(y*y);
  }
}

void KahanLinearRegressionAccumulator::integrate_into(KahanLinearRegressionAccumulator& o)
{
  if(not o.zero_set_) {
    o.x0_ = x0_;
    o.y0_ = y0_;
    o.W_  = W_;
    o.X_  = X_;
    o.Y_  = Y_;
    o.XX_ = XX_;
    o.XY_ = XY_;
    o.YY_ = YY_;
    o.zero_set_ = zero_set_;
  } else {
    double dx0 = x0_ - o.x0_;
    double dy0 = y0_ - o.y0_;
    o.W_.accumulate_from(W_);
    o.X_.accumulate_from(X_);
    o.X_.accumulate(W_.total() * dx0);
    o.Y_.accumulate_from(Y_);
    o.Y_.accumulate(W_.total() * dy0);
    o.XX_.accumulate_from(XX_);
    o.XX_.accumulate(2 * X_.total() * dx0);
    o.XX_.accumulate(W_.total() * dx0*dx0);
    o.XY_.accumulate_from(XY_);
    o.XY_.accumulate(X_.total() * dy0);
    o.XY_.accumulate(dx0 * Y_.total());
    o.XY_.accumulate(W_.total() * dx0*dy0);
    o.YY_.accumulate_from(YY_);
    o.YY_.accumulate(2 * Y_.total() * dy0);
    o.YY_.accumulate(W_.total() * dy0*dy0);
  }
}

void KahanLinearRegressionAccumulator::fit_parameters_and_d2(double& a, double& b, double& D2)
{
  calin::util::log::LOG(calin::util::log::INFO)
    << "XX: " << XX_.total() << " X: " << X_.total() << " W: " << W_.total() << ' '
    << " XY: " << XY_.total() << " Y: " << Y_.total() << " YY: " << YY_.total();

  calin::math::accumulator::KahanAccumulator a_num { XY_ };
  a_num.accumulate(-X_.total()*Y_.total()/W_.total());

  calin::math::accumulator::KahanAccumulator a_den { XX_ };
  a_den.accumulate(-X_.total()*X_.total()/W_.total());

  double a_dir = a_num.total()/a_den.total();
  double b_dir = (Y_.total() - a_dir*X_.total())/W_.total();

  calin::util::log::LOG(calin::util::log::INFO)
    << "Direct A : " << a_dir;

  calin::util::log::LOG(calin::util::log::INFO)
    << "Direct B : " << b_dir;

  calin::math::accumulator::KahanAccumulator d2 { YY_ };
  d2.accumulate_from_with_scaling(XX_, a_dir*a_dir);
  d2.accumulate_from_with_scaling(W_, b_dir*b_dir);
  d2.accumulate_from_with_scaling(XY_, -2*a_dir);
  d2.accumulate_from_with_scaling(X_, 2*a_dir*b_dir);
  d2.accumulate_from_with_scaling(Y_, -2*b_dir);

  calin::util::log::LOG(calin::util::log::INFO)
    << "Direct D2 : " << d2.total();

  Eigen::Matrix2d m;
  m << XX_.total(), X_.total(), X_.total(), W_.total();
  Eigen::Vector2d v;
  v << XY_.total(), Y_.total();
  Eigen::Vector2d ab = m.fullPivHouseholderQr().solve(v);
  a = ab(0);
  b = ab(1);

  D2 = YY_.total() + a*a*XX_.total() + b*b*W_.total() - 2*a*XY_.total() - 2*b*Y_.total() + 2*a*b*X_.total();
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

double calin::math::least_squares::polyval(const double* p, unsigned np, double x)
{
  double y = p[--np];
  while(np) {
    y = y*x + p[--np];
  }
  return y;
}

void calin::math::least_squares::polyval_and_derivative(double& y, double& dy_dx,
  const double* p, unsigned np, double x)
{
  y = p[--np];
  dy_dx = 0.0;
  while(np) {
    dy_dx = dy_dx*x + double(np)*p[np];
    y = y*x + p[--np];
  }
}
