/*

   calin/math/nspace.cpp -- Stephen Fegan -- 2020-11-27

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

#include <exception>

#include <math/nspace.hpp>

using namespace calin::math::nspace;

struct Axis
{
  double xlo;
  double xhi;
  unsigned n;
};

SparseNSpace::
SparseNSpace(const Eigen::VectorXd& xlo, const Eigen::VectorXd& xhi,
    const Eigen::VectorXi& n):
  xlo_(xlo.size()), xhi_(xlo.size()), dx_(xlo.size()), dx_inv_(xlo.size()),
  n_(xlo.size())
{
  if(std::min({xlo.size(), xhi.size(), n.size()}) !=
      std::max({xlo.size(), xhi.size(), n.size()})) {
    throw std::runtime_error("SparseNSpace: xhi, xlo and n must all have same size");
  }

  N_ = 1;
  for(unsigned i=0; i<xlo.size(); ++i) {
    const double dx = (xhi[i] - xlo[i])/n[i];
    xlo_[i] = xlo[i];
    xhi_[i] = xhi[i];
    n_[i] = n[i];
    N_ *= n[i];
    dx_[i] = dx;
    dx_inv_[i] = 1.0/dx;
  }
}

SparseNSpace::SparseNSpace(const std::vector<Axis> axes):
  xlo_(axes.size()), xhi_(axes.size()), dx_(axes.size()), dx_inv_(axes.size()),
  n_(axes.size())
{
  N_ = 1;
  for(unsigned i=0; i<axes.size(); ++i) {
    const Axis& ax ( axes[i] );
    const double dx = (ax.xhi - ax.xlo)/ax.n;
    xlo_[i] = ax.xlo;
    xhi_[i] = ax.xhi;
    n_[i] = ax.n;
    N_ *= ax.n;
    dx_[i] = dx;
    dx_inv_[i] = 1.0/dx;
  }
}

std::vector<calin::math::nspace::Axis> SparseNSpace::axes() const
{
  std::vector<calin::math::nspace::Axis> a;
  for(unsigned i=0; i<xlo_.size(); i++) {
    a.push_back({xlo_[i], xhi_[i], static_cast<unsigned int>(n_[i])});
  }
  return a;
}

double SparseNSpace::total_weight() const {
  double w0 = 0;
  for(auto i : bins_) {
    w0 += i.second;
  }
  return w0;
}

Eigen::VectorXd SparseNSpace::mean_and_total_weight(double& w0) const {
  Eigen::VectorXd w1(xlo_.size());
  Eigen::VectorXd x(xlo_.size());
  w0 = 0;
  w1.setZero();
  for(auto i : bins_) {
    x_center(x, i.first);
    w0 += i.second;
    w1 += x * i.second;
  }
  return w1/w0;
}

Eigen::VectorXd SparseNSpace::mean() const {
  double w0;
  return mean_and_total_weight(w0);
}

Eigen::MatrixXd SparseNSpace::covar_mean_and_total_weight(Eigen::VectorXd& w1, double& w0) const {
  Eigen::MatrixXd w2(xlo_.size(), xlo_.size());
  w1.resize(xlo_.size());
  Eigen::VectorXd x(xlo_.size());
  w0 = 0;
  w1.setZero();
  w2.setZero();
  for(auto i : bins_) {
    x_center(x, i.first);
    w0 += i.second;
    w1 += x * i.second;
    w2 += x * x.transpose() * i.second;
  }
  w1 /= w0;
  return w2/w0 - w1*w1.transpose();
}

Eigen::MatrixXd SparseNSpace::covar() const {
  double w0;
  Eigen::VectorXd w1(xlo_.size());
  return covar_mean_and_total_weight(w1, w0);
}
