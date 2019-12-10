/*

   calin/util/algorithm.hpp -- Stephen Fegan -- 2019-12-10

   Various useful algorithms

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

#pragma once

#include <iostream>
#include <vector>
#include <numeric>      // std::iota
#include <algorithm>    // std::sort

namespace calin { namespace util { namespace algorithm {

template <typename Iterator> std::vector<size_t> argsort(Iterator begin, Iterator end) {
  // initialize original index locations
  std::vector<size_t> index(end-begin);
  std::iota(index.begin(), index.end(), 0);

  // sort indexes based on comparing values in v
  sort(index.begin(), index.end(),
       [&begin](size_t i, size_t j) { return *(begin+i) < *(begin+j); });

  return index;
}

} } } // namespace calin::util::algorithm
