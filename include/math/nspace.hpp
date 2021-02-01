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
#include <list>
#include <cstdlib>

#include <Eigen/Dense>
#include <calin_global_definitions.hpp>

namespace calin { namespace math { namespace nspace {

struct Axis
{
public:
  double xlo;
  double xhi;
  unsigned n;
};

class TreeSparseNSpace
{
public:
  TreeSparseNSpace(const Eigen::VectorXd& xlo, const Eigen::VectorXd& xhi,
    const Eigen::VectorXi& n);
  TreeSparseNSpace(const std::vector<Axis>& axes);

  void clear() { bins_.clear(); }
  void injest(const TreeSparseNSpace& o);

  std::vector<Axis> axes() const;

  unsigned naxes() const { return xlo_.size(); }
  Axis axis(unsigned iaxis) const;

  Eigen::VectorXd axis_bin_centers(unsigned iaxis) const;
  Eigen::VectorXd axis_bin_edges(unsigned iaxis) const;

  uint64_t num_occupied_cells() {
    return bins_.size();
  }

  uint64_t size() const {
    return N_;
  }

  int64_t index(const Eigen::VectorXd& x) const {
    if(x.size() != xlo_.size()) {
      throw std::runtime_error("TreeSparseNSpace: dimensional mismatch");
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
      throw std::runtime_error("TreeSparseNSpace: index out of range");
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
      throw std::runtime_error("TreeSparseNSpace: index out of range");
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
      throw std::runtime_error("TreeSparseNSpace: index out of range");
    }
    if(x.size() != xlo_.size()) {
      throw std::runtime_error("TreeSparseNSpace: dimensional mismatch");
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

  TreeSparseNSpace project_along_axis(unsigned iaxis, unsigned axis_cell_lo, unsigned axis_cell_hi);
  TreeSparseNSpace project_along_axis(unsigned iaxis);

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

class BlockSparseNSpace
{
public:
  BlockSparseNSpace(const Eigen::VectorXd& xlo, const Eigen::VectorXd& xhi,
    const Eigen::VectorXi& n);
  BlockSparseNSpace(const std::vector<Axis>& axes);

  // void clear() { bins_.clear(); }
  // void injest(const BlockSparseNSpace& o);

  // std::vector<Axis> axes() const;
  //
  // unsigned naxes() const { return xlo_.size(); }
  // Axis axis(unsigned iaxis) const;

  // Eigen::VectorXd axis_bin_centers(unsigned iaxis) const;
  // Eigen::VectorXd axis_bin_edges(unsigned iaxis) const;

  // uint64_t num_occupied_cells() {
  //   return array_occupancy_ * block_size_;
  // }
  //
  // uint64_t size() const {
  //   return N_;
  // }
  //
  void accumulate(const Eigen::VectorXd& x, double w = 1.0);

private:
  void index(const Eigen::VectorXd& x, int64_t& array_index, int64_t& block_index) const;
  double& cell_ref(int64_t array_index, int64_t block_index);
  double cell_val(int64_t array_index, int64_t block_index) const;

  Eigen::VectorXd xlo_;
  Eigen::VectorXd xhi_;
  Eigen::VectorXd dx_;
  Eigen::VectorXd dx_inv_;
  Eigen::VectorXi n_;
  int64_t N_;

  unsigned block_shift_;
  unsigned block_mask_;
  unsigned block_size_;

  double* alloc_next_;
  double* alloc_end_;
  unsigned alloc_size_;
  std::list<double*> alloc_free_list_;
  std::list<double*> alloc_all_list_;


  unsigned array_occupancy_;
  unsigned array_size_;
  double overflow_ = 0;
  std::vector<double*> array_;
  unsigned next_index_;
};

CALIN_TYPEALIAS(SparseNSpace, TreeSparseNSpace);

} } } // namespace calin::math::nspace
