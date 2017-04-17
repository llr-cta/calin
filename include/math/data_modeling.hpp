/*

   calin/math/data_modeling.hpp -- Stephen Fegan -- 2017-04-07

   Data modeling functions for various data types

   Copyright 2017, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <stdexcept>

#include "math/histogram.hpp"
#include "math/function.hpp"
#include "math/pdf_1d.hpp"
#include "math/m_estimate.hpp"
#include "math/optimizer.hpp"

namespace calin { namespace math { namespace data_modeling {

class IID1DDataLikelihoodFunction: public calin::math::function::MultiAxisFunction
{
public:
  IID1DDataLikelihoodFunction(calin::math::pdf_1d::Parameterizable1DPDF* pdf,
      const std::vector<double>& x, bool adopt_pdf = false):
    calin::math::function::MultiAxisFunction(),
    pdf_(pdf), adopt_pdf_(adopt_pdf), npar_(pdf_->num_parameters()),
    x_(x), w_(x.size(), 1.0) {
      // nothing to see here
  }

  IID1DDataLikelihoodFunction(calin::math::pdf_1d::Parameterizable1DPDF* pdf,
      const Eigen::VectorXd& x, bool adopt_pdf = false):
    calin::math::function::MultiAxisFunction(),
    pdf_(pdf), adopt_pdf_(adopt_pdf), npar_(pdf_->num_parameters()),
    x_(calin::eigen_to_stdvec(x)), w_(x.size(), 1.0) {
      // nothing to see here
  }

  IID1DDataLikelihoodFunction(calin::math::pdf_1d::Parameterizable1DPDF* pdf,
      const std::vector<double>& x, const std::vector<double>& w,
      bool adopt_pdf = false):
    calin::math::function::MultiAxisFunction(),
    pdf_(pdf), adopt_pdf_(adopt_pdf), npar_(pdf_->num_parameters()),
    x_(x), w_(w) {
      w_.resize(x_.size(), 1.0);
  }

  IID1DDataLikelihoodFunction(calin::math::pdf_1d::Parameterizable1DPDF* pdf,
      const Eigen::VectorXd& x, const Eigen::VectorXd& w,
      bool adopt_pdf = false):
    calin::math::function::MultiAxisFunction(),
    pdf_(pdf), adopt_pdf_(adopt_pdf), npar_(pdf_->num_parameters()),
    x_(calin::eigen_to_stdvec(x)), w_(calin::eigen_to_stdvec(w)) {
      w_.resize(x_.size(), 1.0);
  }

  template<typename Histogram1D> IID1DDataLikelihoodFunction(
      calin::math::pdf_1d::Parameterizable1DPDF* pdf,
      const Histogram1D& hist, bool adopt_pdf = false):
    calin::math::function::MultiAxisFunction(),
    pdf_(pdf), adopt_pdf_(adopt_pdf), npar_(pdf_->num_parameters()),
    x_(calin::eigen_to_stdvec(hist.all_xval_center())),
    w_(calin::eigen_to_stdvec(hist.all_weight())) {
      // nothing to see here
  }

  virtual ~IID1DDataLikelihoodFunction();

  unsigned num_domain_axes() override;
  std::vector<calin::math::function::DomainAxis> domain_axes() override;
  double value(ConstVecRef x) override;
  bool can_calculate_gradient() override;
  double value_and_gradient(ConstVecRef x, VecRef gradient) override;
  bool can_calculate_hessian() override;
  double value_gradient_and_hessian(ConstVecRef x, VecRef gradient,
                                    MatRef hessian) override;
  double error_up() override { return 0.5; }

protected:
  calin::math::pdf_1d::Parameterizable1DPDF* pdf_ = nullptr;
  bool adopt_pdf_ = false;
  unsigned npar_ = 0;
  std::vector<double> x_;
  std::vector<double> w_;
};

class IID1DDataMEstimateLikelihoodFunction: public calin::math::function::MultiAxisFunction
{
public:
  IID1DDataMEstimateLikelihoodFunction(calin::math::pdf_1d::Parameterizable1DPDF* pdf,
      calin::math::m_estimate::LikelihoodRhoFunction* rho,
      const std::vector<double>& x, bool adopt_pdf = false, bool adopt_rho = false):
    calin::math::function::MultiAxisFunction(),
    pdf_(pdf), adopt_pdf_(adopt_pdf), rho_(rho), adopt_rho_(adopt_rho),
    npar_(pdf_->num_parameters()),
    x_(x), w_(x.size(), 1.0) {
      // nothing to see here
  }

  IID1DDataMEstimateLikelihoodFunction(calin::math::pdf_1d::Parameterizable1DPDF* pdf,
      calin::math::m_estimate::LikelihoodRhoFunction* rho,
      const Eigen::VectorXd& x, bool adopt_pdf = false, bool adopt_rho = false):
    calin::math::function::MultiAxisFunction(),
    pdf_(pdf), adopt_pdf_(adopt_pdf), rho_(rho), adopt_rho_(adopt_rho),
    npar_(pdf_->num_parameters()),
    x_(calin::eigen_to_stdvec(x)), w_(x.size(), 1.0) {
      // nothing to see here
  }

  IID1DDataMEstimateLikelihoodFunction(calin::math::pdf_1d::Parameterizable1DPDF* pdf,
      calin::math::m_estimate::LikelihoodRhoFunction* rho,
      const std::vector<double>& x, const std::vector<double>& w,
      bool adopt_pdf = false, bool adopt_rho = false):
    calin::math::function::MultiAxisFunction(),
    pdf_(pdf), adopt_pdf_(adopt_pdf), rho_(rho), adopt_rho_(adopt_rho),
    npar_(pdf_->num_parameters()),
    x_(x), w_(w) {
      w_.resize(x_.size(), 1.0);
  }

  IID1DDataMEstimateLikelihoodFunction(calin::math::pdf_1d::Parameterizable1DPDF* pdf,
      calin::math::m_estimate::LikelihoodRhoFunction* rho,
      const Eigen::VectorXd& x, const Eigen::VectorXd& w,
      bool adopt_pdf = false, bool adopt_rho = false):
    calin::math::function::MultiAxisFunction(),
    pdf_(pdf), adopt_pdf_(adopt_pdf), rho_(rho), adopt_rho_(adopt_rho),
    npar_(pdf_->num_parameters()),
    x_(calin::eigen_to_stdvec(x)), w_(calin::eigen_to_stdvec(w)) {
      w_.resize(x_.size(), 1.0);
  }

  template<typename Histogram1D> IID1DDataMEstimateLikelihoodFunction(
      calin::math::pdf_1d::Parameterizable1DPDF* pdf,
      calin::math::m_estimate::LikelihoodRhoFunction* rho,
      const Histogram1D& hist, bool adopt_pdf = false, bool adopt_rho = false):
    calin::math::function::MultiAxisFunction(),
    pdf_(pdf), adopt_pdf_(adopt_pdf), rho_(rho), adopt_rho_(adopt_rho),
    npar_(pdf_->num_parameters()),
    x_(calin::eigen_to_stdvec(hist.all_xval_center())),
    w_(calin::eigen_to_stdvec(hist.all_weight())) {
      // nothing to see here
  }

  virtual ~IID1DDataMEstimateLikelihoodFunction();

  unsigned num_domain_axes() override;
  std::vector<calin::math::function::DomainAxis> domain_axes() override;
  double value(ConstVecRef x) override;
  bool can_calculate_gradient() override;
  double value_and_gradient(ConstVecRef x, VecRef gradient) override;
  bool can_calculate_hessian() override;
  double value_gradient_and_hessian(ConstVecRef x, VecRef gradient,
                                    MatRef hessian) override;
  double error_up() override { return 0.5; }

protected:
  calin::math::pdf_1d::Parameterizable1DPDF* pdf_ = nullptr;
  bool adopt_pdf_ = false;
  calin::math::m_estimate::LikelihoodRhoFunction* rho_;
  bool adopt_rho_ = false;
  unsigned npar_ = 0;
  std::vector<double> x_;
  std::vector<double> w_;
};


class IID1DDataChi2Function: public calin::math::function::MultiAxisFunction
{
public:
  IID1DDataChi2Function(calin::math::pdf_1d::Parameterizable1DPDF* pdf,
      const std::vector<double>& x, const std::vector<double>& w,
      bool adopt_pdf = false):
    calin::math::function::MultiAxisFunction(),
    pdf_(pdf), adopt_pdf_(adopt_pdf), npar_(pdf_->num_parameters()),
    x_(x), w_(w) {
      w_.resize(std::min(x_.size(), w_.size()));
      x_.resize(w_.size());
  }

  IID1DDataChi2Function(calin::math::pdf_1d::Parameterizable1DPDF* pdf,
      const Eigen::VectorXd& x, const Eigen::VectorXd& w,
      bool adopt_pdf = false):
    calin::math::function::MultiAxisFunction(),
    pdf_(pdf), adopt_pdf_(adopt_pdf), npar_(pdf_->num_parameters()),
    x_(calin::eigen_to_stdvec(x)), w_(calin::eigen_to_stdvec(w)) {
      w_.resize(std::min(x_.size(), w_.size()));
      x_.resize(w_.size());
  }

  template<typename Histogram1D> IID1DDataChi2Function(
      calin::math::pdf_1d::Parameterizable1DPDF* pdf,
      const Histogram1D& hist, bool adopt_pdf = false):
    calin::math::function::MultiAxisFunction(),
    pdf_(pdf), adopt_pdf_(adopt_pdf), npar_(pdf_->num_parameters()),
    x_(calin::eigen_to_stdvec(hist.all_xval_center())),
    w_(calin::eigen_to_stdvec(hist.all_weight())) {
      // nothing to see here
  }

  virtual ~IID1DDataChi2Function();

  unsigned num_domain_axes() override;
  std::vector<calin::math::function::DomainAxis> domain_axes() override;
  double value(ConstVecRef x) override;
  bool can_calculate_gradient() override;
  double value_and_gradient(ConstVecRef x, VecRef gradient) override;
  bool can_calculate_hessian() override;
  double value_gradient_and_hessian(ConstVecRef x, VecRef gradient,
                                    MatRef hessian) override;
  double error_up() override { return 1.0; }

protected:
  calin::math::pdf_1d::Parameterizable1DPDF* pdf_ = nullptr;
  bool adopt_pdf_ = false;
  unsigned npar_ = 0;
  std::vector<double> x_;
  std::vector<double> w_;
  double norm_ = 0;
};


} } } // namespace::math::data_modeling
