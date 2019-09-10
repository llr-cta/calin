/*
   calin/calib/pmt_ses_models.cpp -- Stephen Fegan -- 2017-04-24

   PMT single-electron spectrum models

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

#include <limits>

#include <Eigen/Dense>

#include <calib/pmt_ses_models.hpp>
#include <util/log.hpp>

using namespace calin::util::log;
using namespace calin::calib::pmt_ses_models;

static constexpr double inf = std::numeric_limits<double>::infinity();

TwoGaussianSES::TwoGaussianSES(double dx):
  pattern::delegation::Delegator<ParameterizableSingleAxisFunction>(
    new calin::math::pdf_1d::TwoComponent1DPDF(
      new calin::math::pdf_1d::LimitedGaussianPDF(0, inf, dx), "gauss_lo",
      new calin::math::pdf_1d::LimitedGaussianPDF(0, inf, dx), "gauss_hi",
      true, true), true),
  calin::math::function::BasicReducedSpaceParameterizable<calin::math::function::ParameterizableSingleAxisFunction>(delegate_)
{
  remove_parameter_from_subspace(1, 0.0); // Mean of gauss-lo
}

TwoGaussianSES::~TwoGaussianSES()
{
  // nothing to see here
}

double TwoGaussianSES::gain()
{
  double mx;
  double mxx;
  moments(mx, mxx);
  return mx;
}

double TwoGaussianSES::resolution()
{
  double mx;
  double mxx;
  moments(mx, mxx);
  return std::sqrt(mxx/(mx*mx)-1);
}

void TwoGaussianSES::moments(double& mx, double& mxx)
{
  const double sqrt2_pi = std::sqrt(2.0/M_PI);

  Eigen::VectorXd p = parameter_values();
  double Pl = p(0);
  double sl = p(1);
  double mh = p(2);
  double sh = p(3);

  mx = sqrt2_pi*Pl*sl + (1-Pl)*mh + sqrt2_pi*(1-Pl)/(1+std::erf(mh/(M_SQRT2*sh)))*sh*std::exp(-mh*mh/(2*sh*sh));
  mxx = Pl*(sl*sl - sqrt2_pi*sl*mh) + (1-Pl)*sh*sh + mx*mh;
}

ReparameterizedTwoGaussianSES::ReparameterizedTwoGaussianSES(double dx):
  calin::math::pdf_1d::Parameterizable1DPDF(),
  ses_(dx)
{
  update_cached_values();
}

ReparameterizedTwoGaussianSES::~ReparameterizedTwoGaussianSES()
{
  // nothing to see here
}

unsigned ReparameterizedTwoGaussianSES::num_parameters()
{
  return 4;
}

std::vector<calin::math::function::ParameterAxis> ReparameterizedTwoGaussianSES::parameters()
{
  constexpr double inf = std::numeric_limits<double>::infinity();
  constexpr double tiny_val = std::numeric_limits<double>::min();
  return { { "low_Q_fraction", "1", false, -inf, false, inf },
           { "low_Q_resolution ", "PE", true, tiny_val, false, inf },
           { "gain", "charge_units", true, tiny_val, false, inf },
           { "resolution ", "PE", true, tiny_val, false, inf }
          };
}

Eigen::VectorXd ReparameterizedTwoGaussianSES::parameter_values()
{
  Eigen::VectorXd pval(4);
  pval << prob_lo_, beta_lo_, gain_, beta_;
  return pval;
}

void ReparameterizedTwoGaussianSES::set_parameter_values(ConstVecRef values)
{
  verify_set_parameter_values(values, "ReparameterizedTwoGaussianSES");
  calin::math::function::assign_parameters(values, prob_lo_, beta_lo_, gain_, beta_);
  update_cached_values();
}

bool ReparameterizedTwoGaussianSES::can_calculate_gradient()
{
  return true;
}

bool ReparameterizedTwoGaussianSES::can_calculate_hessian()
{
  return true;
}

bool ReparameterizedTwoGaussianSES::can_calculate_parameter_gradient()
{
  return true;
}

bool ReparameterizedTwoGaussianSES::can_calculate_parameter_hessian()
{
  return false;
}

calin::math::function::DomainAxis ReparameterizedTwoGaussianSES::domain_axis()
{
  constexpr double inf = std::numeric_limits<double>::infinity();
  return { "signal", "signal-units", true, 0, false, inf };
}

double ReparameterizedTwoGaussianSES::value_1d(double x)
{
  return ses_.value_1d(x);
}

double ReparameterizedTwoGaussianSES::
value_and_gradient_1d(double x,  double& dfdx)
{
  return ses_.value_and_gradient_1d(x, dfdx);
}

double ReparameterizedTwoGaussianSES::
value_gradient_and_hessian_1d(double x, double& dfdx, double& d2fdx2)
{
  return ses_.value_gradient_and_hessian_1d(x, dfdx, d2fdx2);
}

double ReparameterizedTwoGaussianSES::
value_and_parameter_gradient_1d(double x,  VecRef gradient)
{
  double value = ses_.value_and_parameter_gradient_1d(x, gradient);
  gradient = jacobian_lu_.solve(gradient);
  return value;
}

double ReparameterizedTwoGaussianSES::
value_parameter_gradient_and_hessian_1d(double x, VecRef gradient, MatRef hessian)
{
  // Transforming the Hessian would require a 3-index tensor in addition to the Jacobian
  throw std::runtime_error("ReparameterizedTwoGaussianSES: cannot calculate Hessian");
}

void ReparameterizedTwoGaussianSES::update_cached_values()
{
  const double sqrt2_pi = std::sqrt(2.0/M_PI);

  const double Pl = prob_lo_;
  const double sl = beta_lo_ * gain_;
  const double mx = gain_;
  const double mxx = (beta_*beta_ + 1)*gain_*gain_;

  double mh = (mx - sqrt2_pi*Pl*sl)/(1-Pl);
  double sh = (mxx - Pl*(sl*sl - sqrt2_pi*sl*mh) - mx*mh)/(1-Pl);
  if(sh < 0)throw std::runtime_error(
    "ReparameterizedTwoGaussianSES: variance of high-Q component negative");

  sh = std::sqrt(sh);

  for(unsigned i = 0;; i++) {
    double mh_est = (mx - sqrt2_pi*Pl*sl - sqrt2_pi*(1-Pl)/(1+std::erf(mh/(M_SQRT2*sh)))*sh*exp(-mh*mh/(2.0*sh*sh)))/(1-Pl);
    double sh_est = (mxx - Pl*(sl*sl - sqrt2_pi*sl*mh_est) - mx*mh_est)/(1-Pl);

    if(sh_est < 0)throw std::runtime_error(
      "ReparameterizedTwoGaussianSES: variance of high-Q component negative");

    sh_est = std::sqrt(sh_est);

    if(std::abs(mh_est - mh)<1e-16*mx  and std::abs(sh_est - sh )<1e-16*mx) {
      mh = mh_est;
      sh = sh_est;
      break;
    }

    if(i>200)throw std::runtime_error(
        "ReparameterizedTwoGaussianSES: too many iterations");
    mh = mh_est;
    sh = sh_est;
  }

  Eigen::VectorXd p(4);
  p << Pl, sl, mh, sh;
  ses_.set_parameter_values(p);

  const double erfc_inv = 1/(1 + std::erf(mh/(M_SQRT2 * sh)));
  const double exp2 = std::exp(-mh*mh/(2*sh*sh));

  const double dg_dPl = sqrt2_pi*sl - mh - sqrt2_pi*sh*erfc_inv*exp2;
  const double dg_dsl = sqrt2_pi*Pl;
  const double dg_dmh =
    (1-Pl)
    - sqrt2_pi*sqrt2_pi*(1-Pl)*erfc_inv*erfc_inv*exp2*exp2
    - sqrt2_pi*(1-Pl)*erfc_inv*exp2*mh/sh;
  const double dg_dsh = sqrt2_pi*(1-Pl)*erfc_inv*exp2 *
    ( 1 + sqrt2_pi*erfc_inv*exp2*mh/sh + mh*mh/(sh*sh));

  const double dmxx_dPl = dg_dPl*mh + sl*sl - sqrt2_pi*sl*mh - sh*sh;
  const double dmxx_dsl = dg_dsl*mh + Pl*(2*sl - sqrt2_pi*mh);
  const double dmxx_dmh = dg_dmh*mh - sqrt2_pi*sl*Pl + mx;
  const double dmxx_dsh = dg_dsh*mh + 2*(1-Pl)*sh;

  const double db_dPl = 1/(2*beta_)*(dmxx_dPl - 2*mxx/mx*dg_dPl)/(mx*mx);
  const double db_dsl = 1/(2*beta_)*(dmxx_dsl - 2*mxx/mx*dg_dsl)/(mx*mx);
  const double db_dmh = 1/(2*beta_)*(dmxx_dmh - 2*mxx/mx*dg_dmh)/(mx*mx);
  const double db_dsh = 1/(2*beta_)*(dmxx_dsh - 2*mxx/mx*dg_dsh)/(mx*mx);

  jacobian_inv_ <<
    1, 0, 0, 0,
    -sl*dg_dPl/(mx*mx), 1/mx - sl*dg_dsl/(mx*mx), -sl*dg_dmh/(mx*mx), -sl*dg_dsh/(mx*mx),
    dg_dPl, dg_dsl, dg_dmh, dg_dsh,
    db_dPl, db_dsl, db_dmh, db_dsh;

  // jacobian_ = jacobian_inv_.inverse().transpose();
  jacobian_lu_ = Eigen::FullPivLU<Eigen::MatrixXd>(jacobian_inv_.transpose());
  // jacobian_lu_transpose_ = Eigen::FullPivLU<Eigen::MatrixXd>(jacobian_inv_);

  double mx_est;
  double mxx_est;
  ses_.moments(mx_est, mxx_est);

  LOG(INFO) << Pl << ' ' << sl << ' ' << mh << ' ' << sh << ' ' << mx << ' ' << mx_est << ' ' << mxx << ' ' << mxx_est;
}

// ********************************* OBSOLETE *********************************
// ********************************* OBSOLETE *********************************
// ********************************* OBSOLETE *********************************
// ********************************* OBSOLETE *********************************
// ********************************* OBSOLETE *********************************
// ********************************* OBSOLETE *********************************
// ********************************* OBSOLETE *********************************
// ********************************* OBSOLETE *********************************

TwoGaussianSESConstrained::TwoGaussianSESConstrained(double dx):
  pattern::delegation::Delegator<ParameterizableSingleAxisFunction>(
    new calin::math::pdf_1d::TwoComponent1DConstraintPDF(
      new calin::math::pdf_1d::TwoComponent1DPDF(
      new calin::math::pdf_1d::LimitedGaussianPDF(0, inf, dx), "gauss_lo",
      new calin::math::pdf_1d::LimitedGaussianPDF(0, inf, dx), "gauss_hi",
      true, true)), true),
  calin::math::function::BasicReducedSpaceParameterizable<calin::math::function::ParameterizableSingleAxisFunction>(delegate_)
{
  remove_parameter_from_subspace(1, 0.0); // Mean of gauss-lo
}

TwoGaussianSESConstrained::~TwoGaussianSESConstrained()
{
  // nothing to see here
}

TwoGaussianSESConstrained_Fast::TwoGaussianSESConstrained_Fast(double dx):
  pattern::delegation::Delegator<ParameterizableSingleAxisFunction>(
    new calin::math::pdf_1d::TwoComponent1DConstraintPDF(
      new calin::math::pdf_1d::TwoComponent1DPDF(
      new calin::math::pdf_1d::LimitedGaussianPDF(0, inf, dx), "gauss_lo",
      new calin::math::pdf_1d::LimitedGaussianPDF(0, inf, dx), "gauss_hi",
      true, true), true), true),
  calin::math::function::BasicReducedSpaceParameterizable<calin::math::function::ParameterizableSingleAxisFunction>(delegate_)
{
  remove_parameter_from_subspace(1, 0.0); // Mean of gauss-lo
}

TwoGaussianSESConstrained_Fast::~TwoGaussianSESConstrained_Fast()
{
  // nothing to see here
}
