/*

   calin/math/histogram.cpp -- Stephen Fegan -- 2015-02-16

   Simple histogramming classes

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <math/accumulator.hpp>
#include <math/histogram.hpp>

namespace calin { namespace math { namespace histogram {

template class BasicHistogram1D<accumulator::SimpleAccumulator>;

} } } // namespace calin::math::histogram

calin::ix::math::histogram::Histogram1DData*
calin::math::histogram::rebin(
  const calin::ix::math::histogram::Histogram1DData& from, unsigned rebinning_factor)
{
  rebinning_factor = std::max(rebinning_factor, 1U);
  auto* to = new calin::ix::math::histogram::Histogram1DData(from);
  to->set_dxval(from.dxval() * rebinning_factor);
  to->clear_bins();
  for(int ibin=0,obin=-1;ibin<from.bins_size();++ibin) {
    if(ibin%rebinning_factor == 0) {
      ++obin;
    }
    to->increment_bins(obin, from.bins(ibin));
  }
  return to;
}
