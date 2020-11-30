/*

   calin/math/nspace.hpp -- Stephen Fegan -- 2020-11-26

   (Re)Implemantion of the U Utah and UCLA N-Space algorithm

   Copyright 2020, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <vector>
#include <map>
#include <cstdlib>

#include <Eigen/Dense>

namespace calin { namespace math { namespace nspace {

class Axis
{
public:
  double xlo;
  double xhi;
  unsigned n;
};

class SparseNSpace
{
public:
  SparseNSpace(const Eigen::VectorXd& xlo, const Eigen::VectorXd& xhi,
    const Eigen::VectorXi& n);
  SparseNSpace(const std::vector<Axis> axes);

  std::vector<Axis> axes() const;

  unsigned naxes() const { return xlo_.size(); }
  Axis axis(unsigned iaxis) const;

  uint64_t num_occupied_cells() {
    return bins_.size();
  }

  uint64_t size() const {
    return N_;
  }

  int64_t index(const Eigen::VectorXd& x) const {
    if(x.size() != xlo_.size()) {
      throw std::runtime_error("SparseNSpace: dimensional mismatch");
    }
    int64_t indx = 0;
    for(int i=0; i<xlo_.size(); i++) {
      int ii = (x(i)-xlo_(i))*dx_inv_(i);
      if(ii<0 or ii>=n_(i))return -1;
      indx = indx*n_(i) + ii;
    }
    return indx;
  }

  Eigen::VectorXd x_corner(int64_t indx) const {
    if(indx < 0 or indx >= N_) {
      throw std::runtime_error("SparseNSpace: index out of range");
    }
    int n = xlo_.size();
    Eigen::VectorXd x(n);
    for(unsigned i=n; i>0;) {
      --i;
      auto qr = std::div(indx, int64_t(n_[i]));
      x[i] = xlo_[i] + dx_[i] * qr.rem;
      indx = qr.quot;
    }
    return x;
  }

  Eigen::VectorXd x_center(int64_t indx) const {
    if(indx < 0 or indx >= N_) {
      throw std::runtime_error("SparseNSpace: index out of range");
    }
    int n = xlo_.size();
    Eigen::VectorXd x(n);
    for(unsigned i=n; i>0;) {
      --i;
      auto qr = std::div(indx, int64_t(n_[i]));
      x[i] = xlo_[i] + dx_[i] * (qr.rem + 0.5);
      indx = qr.quot;
    }
    return x;
  }

  void x_center(Eigen::VectorXd& x, int64_t indx) const {
    if(indx < 0 or indx >= N_) {
      throw std::runtime_error("SparseNSpace: index out of range");
    }
    if(x.size() != xlo_.size()) {
      throw std::runtime_error("SparseNSpace: dimensional mismatch");
    }
    int n = xlo_.size();
    for(unsigned i=n; i>0;) {
      --i;
      auto qr = std::div(indx, int64_t(n_[i]));
      x[i] = xlo_[i] + dx_[i] * (qr.rem + 0.5);
      indx = qr.quot;
    }
  }

  void accumulate(const Eigen::VectorXd& x, double w = 1.0) {
    bins_[index(x)] += w;
  }

  double overflow_weight() const {
    auto ifind = bins_.find(-1);
    if(ifind == bins_.end())return 0;
    return (*ifind).second;
  }

  double weight(const Eigen::VectorXd& x) const {
    int64_t indx = index(x);
    if(indx<0)return 0;
    auto ifind = bins_.find(indx);
    if(ifind == bins_.end())return 0;
    return (*ifind).second;
  }

  typename std::map<int64_t, double>::const_iterator cells_begin() const {
    return bins_.begin();
  }

  typename std::map<int64_t, double>::const_iterator cells_end() const {
    return bins_.end();
  }

  std::vector<int64_t> index_keys() const {
    std::vector<int64_t> v;
    v.reserve(bins_.size());
    for(auto i : bins_) {
      v.emplace_back(i.first);
    }
    return v;
  }

  double weight_by_index(int64_t indx) const {
    if(indx<0 or indx>=N_)return 0;
    auto ifind = bins_.find(indx);
    if(ifind == bins_.end())return 0;
    return (*ifind).second;
  }

  SparseNSpace project_along_axis(unsigned iaxis, unsigned axis_cell_lo, unsigned axis_cell_hi);
  SparseNSpace project_along_axis(unsigned iaxis);

  Eigen::VectorXd as_vector() const;
  Eigen::MatrixXd as_matrix() const;

  double total_weight() const;

  Eigen::VectorXd mean_and_total_weight(double& w0) const;
  Eigen::VectorXd mean() const;

  Eigen::MatrixXd covar_mean_and_total_weight(Eigen::VectorXd& w1, double& w0) const;
  Eigen::MatrixXd covar() const;

private:
  Eigen::VectorXd xlo_;
  Eigen::VectorXd xhi_;
  Eigen::VectorXd dx_;
  Eigen::VectorXd dx_inv_;
  Eigen::VectorXi n_;
  int64_t N_;

  std::map<int64_t, double> bins_;
};

} } } // namespace calin::math::nspace
