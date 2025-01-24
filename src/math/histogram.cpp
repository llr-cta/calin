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

calin::math::histogram::Histogram1D* calin::math::histogram::new_histogram_if_enabled(
  const calin::ix::math::histogram::AccumulatedAndSerializedHistogram1DConfig& config)
{
  if(config.enable() and config.dxval()!=0)
    return new calin::math::histogram::Histogram1D(config);
  return nullptr;
}

calin::ix::math::histogram::Histogram1DData* calin::math::histogram::compactify(
  const calin::ix::math::histogram::Histogram1DData& original_hist,
  int max_dense_bins_in_output, int max_output_rebinning,
  calin::ix::math::histogram::Histogram1DData* compactified_hist)
{
  if(max_dense_bins_in_output<=0 or max_output_rebinning==1) {
    // No rebinning needed or allowed, so just sparsify into output hist
    return calin::math::histogram::sparsify(original_hist, compactified_hist);
  } else {
    std::unique_ptr<calin::ix::math::histogram::Histogram1DData> hist {
      calin::math::histogram::sparsify(original_hist) };
    int rebin = 1;
    if(hist->bins_size()>max_dense_bins_in_output) {
      rebin = (hist->bins_size()+max_dense_bins_in_output-1)/max_dense_bins_in_output;
      if(max_output_rebinning != 0) {
        rebin = std::min(rebin, max_output_rebinning);
      }
    }
    if(rebin != 1) {
      return calin::math::histogram::rebin(*hist, rebin, compactified_hist);
    } else if(compactified_hist) {
      compactified_hist->Clear();
      compactified_hist->CopyFrom(*hist);
      return compactified_hist;
    } else {
      return hist.release();
    }
  }
}

calin::ix::math::histogram::Histogram1DData*
calin::math::histogram::rebin(
  const calin::ix::math::histogram::Histogram1DData& original_hist, int rebinning_factor,
  calin::ix::math::histogram::Histogram1DData* rebinned_hist)
{
  rebinning_factor = std::max(rebinning_factor, 1);
  if(rebinned_hist == nullptr) {
    rebinned_hist = new calin::ix::math::histogram::Histogram1DData;
  } else {
    rebinned_hist->Clear();
  }
  rebinned_hist->CopyFrom(original_hist);
  rebinned_hist->set_dxval(original_hist.dxval() * rebinning_factor);
  rebinned_hist->clear_bins();
  for(int ibin=0,obin=-1;ibin<original_hist.bins_size();++ibin) {
    if(ibin%rebinning_factor == 0) {
      ++obin;
    }
    rebinned_hist->increment_bins(obin, original_hist.bins(ibin));
  }
  rebinned_hist->mutable_sparse_bins()->clear();
  for(auto isparse : original_hist.sparse_bins()) {
    int ibin = (isparse.first>=0) ? (isparse.first/rebinning_factor) :
      ((isparse.first - rebinning_factor + 1)/rebinning_factor);
    if(ibin >= 0 and ibin < rebinned_hist->bins_size()) {
      rebinned_hist->increment_bins(ibin, isparse.second);
    } else {
      (*rebinned_hist->mutable_sparse_bins())[ibin] += isparse.second;
    }
  }
  return rebinned_hist;
}

calin::ix::math::histogram::Histogram1DData*
calin::math::histogram::sparsify(
  const calin::ix::math::histogram::Histogram1DData& original_hist,
  calin::ix::math::histogram::Histogram1DData* sparsified_hist)
{
  constexpr int datum_cost = 8;
  constexpr int sparse_penelty = 4;
  int ihi = original_hist.bins_size();
  int ilo = 0;
  while(ihi>ilo and original_hist.bins(ihi-1) == 0)--ihi;
  while(ilo<ihi and original_hist.bins(ilo) == 0)++ilo;
  int size = (ihi-ilo)*datum_cost;

  int best_size = size;
  int best_ihi = ihi;
  int best_ilo = ilo;

  while(ihi > ilo) {
    int test_ihi = ihi;
    int test_size_hi = size;
    while(test_ihi>ilo and original_hist.bins(test_ihi-1) != 0)--test_ihi, test_size_hi += sparse_penelty;
    while(test_ihi>ilo and original_hist.bins(test_ihi-1) == 0)--test_ihi, test_size_hi -= datum_cost;

    int test_ilo = ilo;
    int test_size_lo = size;
    while(test_ilo<ihi and original_hist.bins(test_ilo) != 0)++test_ilo, test_size_lo += sparse_penelty;
    while(test_ilo<ihi and original_hist.bins(test_ilo) == 0)++test_ilo, test_size_lo -= datum_cost;

#if 0
    LOG(INFO) << ilo << "," << ihi << "," << size << " -- "
      << ilo << "," << test_ihi << "," << test_size_hi << " -- "
      << test_ilo << "," << ihi << "," << test_size_lo << " -- "
      << best_ilo << "," << best_ihi << "," << best_size;
#endif

    if(test_size_hi <= test_size_lo) {
      ihi = test_ihi;
      size = test_size_hi;
    } else {
      ilo = test_ilo;
      size = test_size_lo;
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

calin::ix::math::histogram::Histogram1DData*
calin::math::histogram::trim_fraction(
  const calin::ix::math::histogram::Histogram1DData& original_hist_data,
  double trim_fraction_left, double trim_fraction_right, bool winsonize, bool integer_weights,
  calin::ix::math::histogram::Histogram1DData* trimmed_hist_data)
{
  calin::math::histogram::Histogram1D hist(original_hist_data);
  hist.trim_fraction(trim_fraction_left, trim_fraction_right, winsonize, integer_weights);
  return hist.dump_as_proto(trimmed_hist_data);
}
