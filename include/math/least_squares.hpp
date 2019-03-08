/*

   calin/math/polyfit.hpp -- Stephen Fegan -- 2019-03-07

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

#pragma once

#include <Eigen/Dense>

namespace calin { namespace math { namespace least_squares {

Eigen::VectorXd polyfit(const Eigen::VectorXd& x, const Eigen::VectorXd& y, unsigned order);

} } } // namespace calin::math::least_squares
