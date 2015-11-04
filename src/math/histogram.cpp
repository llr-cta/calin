/* 

   calin/math/histogram.cpp -- Stephen Fegan -- 2015-02-16

   Simple histogramming classes

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include "math/accumulator.hpp"
#include "math/histogram.hpp"

namespace calin { namespace math { namespace histogram {

template class BasicHistogram1D<accumulator::SimpleAccumulator>;

} } } // namespace calin::math::histogram
