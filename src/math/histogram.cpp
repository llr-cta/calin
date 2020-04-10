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

#include <algorithm>

#include <math/accumulator.hpp>
#include <math/histogram.hpp>
#include <util/log.hpp>

using namespace calin::util::log;

namespace calin { namespace math { namespace histogram {

template class BasicHistogram1D<accumulator::SimpleAccumulator>;

} } } // namespace calin::math::histogram

calin::ix::math::histogram::Histogram1DData*
calin::math::histogram::rebin(
  const calin::ix::math::histogram::Histogram1DData& original_hist, unsigned rebinning_factor,
  calin::ix::math::histogram::Histogram1DData* rebinned_hist)
{
  rebinning_factor = std::max(rebinning_factor, 1U);
  if(rebinned_hist == nullptr)
    rebinned_hist = new calin::ix::math::histogram::Histogram1DData;
  else
    rebinned_hist->Clear();
  rebinned_hist->CopyFrom(original_hist);
  rebinned_hist->set_dxval(original_hist.dxval() * rebinning_factor);
  rebinned_hist->clear_bins();
  for(int ibin=0,obin=-1;ibin<original_hist.bins_size();++ibin) {
    if(ibin%rebinning_factor == 0) {
      ++obin;
    }
    rebinned_hist->increment_bins(obin, original_hist.bins(ibin));
  }
  return rebinned_hist;
}

calin::ix::math::histogram::Histogram1DData*
calin::math::histogram::sparsify(
  const calin::ix::math::histogram::Histogram1DData& original_hist,
  calin::ix::math::histogram::Histogram1DData* sparsified_hist)
{
  constexpr int sparse_penelty = 1;
  int ihi = original_hist.bins_size();
  int ilo = 0;
  while(ihi>ilo and original_hist.bins(ihi-1) == 0)--ihi;
  while(ilo<ihi and original_hist.bins(ilo) == 0)++ilo;
  int size = ihi-ilo;;

  int best_size = size;
  int best_ihi = ihi;
  int best_ilo = ilo;

  int test_ihi = ihi;
  int test_size_hi = size;
  while(test_ihi>ilo and original_hist.bins(test_ihi-1) != 0)--test_ihi, test_size_hi += sparse_penelty;
  while(test_ihi>ilo and original_hist.bins(test_ihi-1) == 0)--test_ihi, test_size_hi -= 1;

  int test_ilo = ilo;
  int test_size_lo = size;
  while(test_ilo<ihi and original_hist.bins(test_ilo) != 0)++test_ilo, test_size_lo += sparse_penelty;
  while(test_ilo<ihi and original_hist.bins(test_ilo) == 0)++test_ilo, test_size_lo -= 1;

  while(ihi > ilo) {
#if 0
    LOG(INFO) << ilo << "," << ihi << "," << size << " -- "
      << ilo << "," << test_ihi << "," << test_size_hi << " -- "
      << test_ilo << "," << ihi << "," << test_size_lo << " -- "
      << best_ilo << "," << best_ihi << "," << best_size;
#endif

    if(test_size_hi <= test_size_lo) {
      ihi = test_ihi;
      test_size_lo += test_size_hi - size;
      size = test_size_hi;
      while(test_ihi>ilo and original_hist.bins(test_ihi-1) != 0)--test_ihi, test_size_hi += sparse_penelty;
      while(test_ihi>ilo and original_hist.bins(test_ihi-1) == 0)--test_ihi, test_size_hi -= 1;
    } else {
      ilo = test_ilo;
      test_size_hi += test_size_lo - size;
      size = test_size_lo;
      while(test_ilo<ihi and original_hist.bins(test_ilo) != 0)++test_ilo, test_size_lo += sparse_penelty;
      while(test_ilo<ihi and original_hist.bins(test_ilo) == 0)++test_ilo, test_size_lo -= 1;
    }

    if(size < best_size) {
      best_ihi = ihi;
      best_ilo = ilo;
      best_size = size;
    }
  }

  if(sparsified_hist == nullptr) {
    sparsified_hist = new calin::ix::math::histogram::Histogram1DData;
  } else {
    sparsified_hist->Clear();
  }
  sparsified_hist->CopyFrom(original_hist);

  sparsified_hist->clear_bins();
  sparsified_hist->mutable_sparse_bins()->clear();
  for(auto isparse : original_hist.sparse_bins()) {
    (*sparsified_hist->mutable_sparse_bins())[isparse.first - best_ilo] = isparse.second;
  }
  for(int ibin = 0; ibin<original_hist.bins_size(); ibin++) {
    if(ibin < best_ilo or ibin >= best_ihi) {
      if(original_hist.bins(ibin) != 0)
        (*sparsified_hist->mutable_sparse_bins())[ibin - best_ilo] += original_hist.bins(ibin);
    } else {
      sparsified_hist->add_bins(original_hist.bins(ibin));
    }
  }
  sparsified_hist->set_xval0(original_hist.xval0() + original_hist.dxval()*best_ilo);
  return sparsified_hist;
}

calin::ix::math::histogram::Histogram1DData*
calin::math::histogram::densify(
  const calin::ix::math::histogram::Histogram1DData& original_hist,
  calin::ix::math::histogram::Histogram1DData* densified_hist)
{
  if(densified_hist == nullptr) {
    densified_hist = new calin::ix::math::histogram::Histogram1DData;
  } else {
    densified_hist->Clear();
  }
  densified_hist->CopyFrom(original_hist);
  densified_hist->mutable_sparse_bins()->clear();

  if(not original_hist.sparse_bins().empty()) {
    int first_bin = original_hist.sparse_bins().begin()->first;
    int last_bin = original_hist.sparse_bins().begin()->first + 1;
    for(auto isparse : original_hist.sparse_bins()) {
      first_bin = std::min(first_bin, isparse.first);
      last_bin = std::max(last_bin, isparse.first + 1);
    }

    if(original_hist.bins_size()) {
      first_bin = std::min(0, first_bin);
      last_bin = std::max(original_hist.bins_size(), last_bin);
    }

    densified_hist->mutable_bins()->Resize(last_bin-first_bin,0);

    if(original_hist.bins_size() and first_bin!=0) {
      std::fill(densified_hist->mutable_bins()->mutable_data(),
        densified_hist->mutable_bins()->mutable_data()-first_bin, 0);
      std::copy(
        original_hist.bins().begin(), original_hist.bins().end(),
        densified_hist->mutable_bins()->mutable_data()-first_bin);
    }

    densified_hist->set_xval0(original_hist.xval0() + first_bin*original_hist.dxval());

    for(auto isparse : original_hist.sparse_bins()) {
      (*densified_hist->mutable_bins())[isparse.first - first_bin] += isparse.second;
    }
  }

  return densified_hist;
}
