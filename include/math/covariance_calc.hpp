/*

   calin/math/covariance_calc.hpp -- Stephen Fegan -- 2016-04-24

   Utility functions for covariance calculation

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
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

namespace calin { namespace math { namespace covariance_calc {

double cov_i64_gen(int64_t sij, int64_t nij,
  int64_t si, int64_t ni, int64_t sj, int64_t nj);

} } } // namespace calin::math::covariance_calc
