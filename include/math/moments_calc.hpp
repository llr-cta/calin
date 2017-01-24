/*

   calin/math/moments_calc.hpp -- Stephen Fegan -- 2017-01-24

   Calculate 1D and 2D moments

   Copyright 2017, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, Ecole polytechnique, CNRS/IN2P3, Universite Paris-Saclay

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

#include <math/accumulator.hpp>

namespace calin { namespace math { namespace moments_calc {

class MomentsCalc1D
{
public:
  MomentsCalc1D();
  double mean();
  double sum_w() { return acc_w.total(); }
  double sum_wx() { return acc_wx.total(); }
protected:
  calin::math::accumulator::RecommendedAccumulator acc_w;
  calin::math::accumulator::RecommendedAccumulator acc_wx;
};

} } }
