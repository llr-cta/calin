/*

   calin/math/histogram_mif.cpp -- Stephen Fegan -- 2017-03-14

   Message integration function for Histogram1DData

   Copyright 2017, Stephen Fegan <sfegan@llr.in2p3.fr>
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
#include <list>

#include "math/histogram.pb.h"

void calin::ix::math::histogram::Histogram1DData::IntegrateFrom(
  const calin::ix::math::histogram::Histogram1DData& from)
{
  calin::ix::math::histogram::Histogram1DData* to = this;
  if(from.dxval() == 0) {
    // Empty histogram - do nbothing !
    return;
  }
  if(to->dxval() and (to->dxval() != from.dxval() or
      to->xval_align() != from.xval_align())) {
    std::string warning = "merge_histogram1d_data: incompatible bin configurations: ";
    if(not from.name().empty()) {
      warning += from.name() + std::string(": ");
    }
    if(to->dxval() != from.dxval()) {
      warning += std::to_string(to->dxval()) + "!=" + std::to_string(from.dxval());
      if(to->xval_align() != from.xval_align()) {
        warning += " and " + std::to_string(to->xval_align()) + "!=" + std::to_string(from.xval_align());
      }
    } else {
      warning += std::to_string(to->xval_align()) + "!=" + std::to_string(from.xval_align());
    }

    throw std::invalid_argument(warning);
  }
  to->set_dxval(from.dxval());
  to->set_xval_align(from.xval_align());

  if(from.bins_size()) // If "from" has dense data then we must copy it into "to"
  {
    if(to->bins_size()) // If "to" also has dense data we must do an overlapping merge
    {
      double xval0 = std::min(to->xval0(), from.xval0());
      int to_offset = int((to->xval0() - xval0)/to->dxval());
      int from_offset = int((from.xval0() - xval0)/to->dxval());
      int to_n = to->bins_size();
      int N = std::max(to_n + to_offset, from.bins_size() + from_offset);

      to->mutable_bins()->Resize(N,0);
      if(to_offset) // Must move data in "to" array to later point
      {
        std::copy_backward(to->bins().data(), to->bins().data()+to_n,
          to->mutable_bins()->mutable_data()+to_offset+to_n);
        std::fill(to->mutable_bins()->mutable_data(),
          to->mutable_bins()->mutable_data()+to_offset, 0);
        to->set_xval0(xval0);
      }
      std::transform(from.bins().data(), from.bins().data()+from.bins_size(),
        to->bins().data()+from_offset,
        to->mutable_bins()->mutable_data()+from_offset,
        [](double a, double b) { return a+b; });

      std::list<std::pair<int, double> > new_sparse;
      for(auto isparse : to->sparse_bins()) {
        int ibin = isparse.first + to_offset;
        if(ibin<0 or ibin>=N) {
          new_sparse.emplace_back(ibin, isparse.second) ;
        } else {
          (*to->mutable_bins())[ibin] += isparse.second;
        }
      }
      to->mutable_sparse_bins()->clear();
      for(auto isparse : new_sparse) {
        (*to->mutable_sparse_bins())[isparse.first] = isparse.second;
      }
      for(auto isparse : from.sparse_bins()) {
        int ibin = isparse.first + from_offset;
        if(ibin<0 or ibin>=N) {
          (*to->mutable_sparse_bins())[ibin] += isparse.second;
        } else {
          (*to->mutable_bins())[ibin] += isparse.second;
        }
      }

      to->set_xval_min(std::min(to->xval_min(), from.xval_min()));
      to->set_xval_max(std::min(to->xval_max(), from.xval_max()));
    }
    else // "to" has no dense data in bins, just copy "from"
    {
      double xval0 = from.xval0();
      int to_offset = int((to->xval0() - xval0)/to->dxval());
      int from_offset = 0;
      int N = from.bins_size();
      to->mutable_bins()->Resize(N, 0);

      std::copy(from.bins().data(), from.bins().data()+N,
        to->mutable_bins()->mutable_data());
      to->set_xval0(from.xval0());

      std::list<std::pair<int, double> > new_sparse;
      for(auto isparse : to->sparse_bins()) {
        int ibin = isparse.first + to_offset;
        if(ibin<0 or ibin>=N) {
          new_sparse.emplace_back(ibin, isparse.second) ;
        } else {
          (*to->mutable_bins())[ibin] += isparse.second;
        }
      }
      to->mutable_sparse_bins()->clear();
      for(auto isparse : new_sparse) {
        (*to->mutable_sparse_bins())[isparse.first] = isparse.second;
      }
      for(auto isparse : from.sparse_bins()) {
        int ibin = isparse.first + from_offset;
        if(ibin<0 or ibin>=N) {
          (*to->mutable_sparse_bins())[ibin] += isparse.second;
        } else {
          (*to->mutable_bins())[ibin] += isparse.second;
        }
      }

      to->set_xval_min(from.xval_min());
      to->set_xval_max(from.xval_max());
    }
  } else { // "from" has no dense data
    double xval0 = to->xval0();
    int from_offset = int((from.xval0() - xval0)/to->dxval());
    int N = to->bins_size();
    for(auto isparse : from.sparse_bins()) {
      int ibin = isparse.first + from_offset;
      if(ibin<0 or ibin>=N) {
        (*to->mutable_sparse_bins())[ibin] += isparse.second;
      } else {
        (*to->mutable_bins())[ibin] += isparse.second;
      }
    }
  }

  to->set_limited(from.limited());
  to->set_xval_limit_lo(from.xval_limit_lo());
  to->set_xval_limit_hi(from.xval_limit_hi());
  to->set_overflow_lo(to->overflow_lo() + from.overflow_lo());
  to->set_overflow_hi(to->overflow_hi() + from.overflow_hi());
  to->set_sum_w(to->sum_w() + from.sum_w());
  to->set_sum_wx(to->sum_wx() + from.sum_wx());
  to->set_sum_wxx(to->sum_wxx() + from.sum_wxx());
  to->set_name(from.name());
  to->set_xval_units(from.xval_units());
  to->set_weight_units(from.weight_units());
}
