/*

   calin/math/moments_calc.hpp -- Stephen Fegan -- 2017-01-24

   Calculate 1D and 2D moments

   Copyright 2017, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <math/special.hpp>
#include <math/accumulator.hpp>

namespace calin { namespace math { namespace moments_calc {

class SecondMomentsCalc1D
{
public:
  SecondMomentsCalc1D() { /* nothing to see here */ }
  void accumulate(double x, double w = 1.0) {
    acc_w.accumulate(w);
    w *= x;
    acc_wx.accumulate(w);
    w *= x;
    acc_wxx.accumulate(w);
  }
  double mean() { return sum_wx()/sum_w(); }
  double var() { return sum_wxx()/sum_w() - calin::math::special::SQR(mean()); }
  double sum_w() { return acc_w.total(); }
  double sum_wx() { return acc_wx.total(); }
  double sum_wxx() { return acc_wxx.total(); }
  void reset() { acc_w.reset(); acc_wx.reset(); acc_wxx.reset(); }
protected:
  calin::math::accumulator::RecommendedAccumulator acc_w {};
  calin::math::accumulator::RecommendedAccumulator acc_wx {};
  calin::math::accumulator::RecommendedAccumulator acc_wxx {};
};

class SecondMomentsCalc2D
{
public:
  SecondMomentsCalc2D() { /* nothing to see here */ }
  void accumulate(double x, double y, double w = 1.0) {
    acc_w.accumulate(w);
    double wx = w * x;
    acc_wx.accumulate(wx);
    acc_wxy.accumulate(wx*y);
    wx *= x;
    acc_wxx.accumulate(wx);
    double wy = w * y;
    acc_wy.accumulate(wy);
    wy *= y;
    acc_wyy.accumulate(wy);
  }
  double mean_x() { return sum_wx()/sum_w(); }
  double mean_y() { return sum_wy()/sum_w(); }
  double var_x() { return sum_wxx()/sum_w() - calin::math::special::SQR(mean_x()); }
  double var_y() { return sum_wyy()/sum_w() - calin::math::special::SQR(mean_y()); }
  double covar() { return sum_wxy()/sum_w() - mean_x()*mean_y(); }
  double sum_w() { return acc_w.total(); }
  double sum_wx() { return acc_wx.total(); }
  double sum_wy() { return acc_wy.total(); }
  double sum_wxx() { return acc_wxx.total(); }
  double sum_wxy() { return acc_wxy.total(); }
  double sum_wyy() { return acc_wyy.total(); }
  void reset() { acc_w.reset(); acc_wx.reset(); acc_wy.reset();
    acc_wxx.reset(); acc_wxy.reset(); acc_wyy.reset();}
protected:
  calin::math::accumulator::RecommendedAccumulator acc_w {};
  calin::math::accumulator::RecommendedAccumulator acc_wx {};
  calin::math::accumulator::RecommendedAccumulator acc_wy {};
  calin::math::accumulator::RecommendedAccumulator acc_wxx {};
  calin::math::accumulator::RecommendedAccumulator acc_wxy {};
  calin::math::accumulator::RecommendedAccumulator acc_wyy {};
};

class ThirdMomentsCalc2D
{
public:
  ThirdMomentsCalc2D() { /* nothing to see here */ }
  void accumulate(double x, double y, double w = 1.0) {
    acc_w.accumulate(w);
    double wx = w * x;
    acc_wx.accumulate(wx);
    acc_wxy.accumulate(wx*y);
    wx *= x;
    acc_wxx.accumulate(wx);
    acc_wxxy.accumulate(wx*y);
    acc_wxxx.accumulate(wx*x);
    double wy = w * y;
    acc_wy.accumulate(wy);
    wy *= y;
    acc_wyy.accumulate(wy);
    acc_wxyy.accumulate(wy*x);
    acc_wyyy.accumulate(wy*y);
  }
  double mean_x() { return sum_wx()/sum_w(); }
  double mean_y() { return sum_wy()/sum_w(); }
  double var_x() { return sum_wxx()/sum_w() - calin::math::special::SQR(mean_x()); }
  double var_y() { return sum_wyy()/sum_w() - calin::math::special::SQR(mean_y()); }
  double covar() { return sum_wxy()/sum_w() - mean_x()*mean_y(); }
  double skew_x() {
    return (sum_wxxx()/sum_w() - 3.0*var_x()*mean_x()
              - calin::math::special::CUBE(mean_x()))/
      calin::math::special::CUBE(std::sqrt(var_x())); }
  double skew_y() {
    return (sum_wyyy()/sum_w() - 3.0*var_y()*mean_y()
              - calin::math::special::CUBE(mean_y()))/
      calin::math::special::CUBE(std::sqrt(var_y())); }
  double sum_w() { return acc_w.total(); }
  double sum_wx() { return acc_wx.total(); }
  double sum_wy() { return acc_wy.total(); }
  double sum_wxx() { return acc_wxx.total(); }
  double sum_wxy() { return acc_wxy.total(); }
  double sum_wyy() { return acc_wyy.total(); }
  double sum_wxxx() { return acc_wxxx.total(); }
  double sum_wxxy() { return acc_wxxy.total(); }
  double sum_wxyy() { return acc_wxyy.total(); }
  double sum_wyyy() { return acc_wyyy.total(); }
  void reset() { acc_w.reset(); acc_wx.reset(); acc_wy.reset();
    acc_wxx.reset(); acc_wxy.reset(); acc_wyy.reset();
    acc_wxxx.reset(); acc_wxxy.reset(); acc_wxyy.reset(); acc_wyyy.reset(); }
protected:
  calin::math::accumulator::RecommendedAccumulator acc_w {};
  calin::math::accumulator::RecommendedAccumulator acc_wx {};
  calin::math::accumulator::RecommendedAccumulator acc_wy {};
  calin::math::accumulator::RecommendedAccumulator acc_wxx {};
  calin::math::accumulator::RecommendedAccumulator acc_wxy {};
  calin::math::accumulator::RecommendedAccumulator acc_wyy {};
  calin::math::accumulator::RecommendedAccumulator acc_wxxx {};
  calin::math::accumulator::RecommendedAccumulator acc_wxxy {};
  calin::math::accumulator::RecommendedAccumulator acc_wxyy {};
  calin::math::accumulator::RecommendedAccumulator acc_wyyy {};
};

} } }
