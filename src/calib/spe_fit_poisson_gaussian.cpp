/*
   calin/calib/spe_fit_poisson_gaussian.cpp -- Stephen Fegan -- 2015-03-01

   Functions to do fit to multi-electron spectrum in the "single PE"
   domain. Poisson-Gaussian models.

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, Ecole Polytechnique, CNRS/IN2P3

   This file is part of "calin"

   "calin" is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 or
   later, as published by the Free Software Foundation.

   "calin" is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.
*/

#include <iostream>
#include <cmath>
#include <limits>
#include <cassert>
#include <memory>

#include <fftw3.h>

#include <math/fftw_util.hpp>
#include <math/accumulator.hpp>
#include <math/special.hpp>
#include <calib/spe_fit.hpp>
#include <calib/pmt_model_pg.hpp>
#include <io/log.hpp>

using namespace calin::math;
using namespace calin::calib::spe_fit;
using namespace calin::calib::pmt_model_pg;
using namespace calin::io::log;
using calin::math::special::SQR;
using calin::math::function::assign_parameters;

namespace {

constexpr double c_gauss_norm = 0.5*M_2_SQRTPI*M_SQRT1_2;

} // anonymous namespace

// ============================================================================
//
// PoissonGaussMES -Fast PG model for single PE regime where maximum
// number of PEs is specified
//
// ============================================================================

PoissonGaussianMES::PoissonGaussianMES(unsigned nmax, bool force_calc_hessian):
    nmax_(nmax), C1_(nmax), C2_(nmax), C3_(nmax),
    dC1_dF_(nmax), dC1_ds_(nmax), dC1_dg_(nmax), dC1_db_(nmax),
    dC2_ds_(nmax), dC2_dg_(nmax), dC2_db_(nmax), dC3_dg_(nmax),
    force_calc_hessian_(force_calc_hessian)
{
  calc_cached_vars();
}

PoissonGaussianMES::~PoissonGaussianMES()
{
  // nothing to see here
}

unsigned PoissonGaussianMES::num_parameters()
{
  return 5;
}

std::vector<function::ParameterAxis> PoissonGaussianMES::parameters()
{
  constexpr double tiny_val = std::numeric_limits<double>::min();
  constexpr double inf = std::numeric_limits<double>::infinity();
  //constexpr double phuge_val = std::numeric_limits<double>::max();
  //constexpr double nhuge_val = std::numeric_limits<double>::lowest();
  return { { "light_intensity", "PE", true, 0, false, inf },
    { "ped_zero", "DC", false, -inf, false, inf },
    { "ped_width", "DC", true, tiny_val, false, inf },
    { "gain", "DC/PE", false, -inf, false, inf }, // can be negative
    { "ses_width", "PE", true, tiny_val, false, inf }
  };
}

Eigen::VectorXd PoissonGaussianMES::parameter_values()
{
  Eigen::VectorXd vec(5);
  vec << intensity_pe_, ped_zero_dc_, ped_rms_dc_,
      gain_dc_pe_, ses_rms_pe_;
  return vec;
}

void PoissonGaussianMES::set_parameter_values(ConstVecRef values)
{
  verify_set_parameter_values(values, "PoissonGaussianMES");
  assign_parameters(values, intensity_pe_, ped_zero_dc_, ped_rms_dc_,
                    gain_dc_pe_, ses_rms_pe_);
  calc_cached_vars();
}

bool PoissonGaussianMES::can_calculate_parameter_gradient()
{
  return true;
}

bool PoissonGaussianMES::can_calculate_parameter_hessian()
{
  return true;
}

double PoissonGaussianMES::pdf_mes(double x)
{
  const double x0 = ped_zero_dc_;
  const double xc = x-x0;
  accumulator apdf;
  for(unsigned n=nmin_; n<nmax_; n++)
  {
    const double sqr_arg = xc-C3_[n];
    const double log_pn = C1_[n] - C2_[n]*SQR(sqr_arg);
    apdf.accumulate(std::exp(log_pn));
  }
  return c_gauss_norm * apdf.total();;
}

double PoissonGaussianMES::pdf_gradient_mes(double x, VecRef gradient)
{
  gradient.resize(5);

  const double x0 = ped_zero_dc_;
  const double xc = x-x0;

  accumulator apdf;
  accumulator agradient[5];

  for(unsigned n=nmin_; n<nmax_; n++)
  {
    const double xn = xc-C3_[n];
    const double xn2 = SQR(xn);
    const double log_pn = C1_[n] - C2_[n]*xn2;
    const double pdf_n = std::exp(log_pn);

    apdf.accumulate(pdf_n);

    const double dlog0 = dC1_dF_[n];
    const double dlog1 = C2_[n]*2.0*xn;
    const double dlog2 = dC1_ds_[n] - dC2_ds_[n]*xn2;
    const double dlog3 = dC1_dg_[n] - dC2_dg_[n]*xn2 + dC3_dg_[n]*xn;
    const double dlog4 = dC1_db_[n] - dC2_db_[n]*xn2;

    agradient[0].accumulate(pdf_n*dlog0);
    agradient[1].accumulate(pdf_n*dlog1);
    agradient[2].accumulate(pdf_n*dlog2);
    agradient[3].accumulate(pdf_n*dlog3);
    agradient[4].accumulate(pdf_n*dlog4);
  }

  gradient(0) = c_gauss_norm * agradient[0].total();
  gradient(1) = c_gauss_norm * agradient[1].total();
  gradient(2) = c_gauss_norm * agradient[2].total();
  gradient(3) = c_gauss_norm * agradient[3].total();
  gradient(4) = c_gauss_norm * agradient[4].total();
  return c_gauss_norm * apdf.total();
}

double PoissonGaussianMES::
pdf_gradient_hessian_mes(double x, VecRef gradient, MatRef hessian)
{
  gradient.resize(5);
  hessian.resize(5,5);

  if(!hessian_elements_good_)calc_cached_vars(true);

  const double x0 = ped_zero_dc_;
  const double xc = x-x0;

  accumulator apdf;
  accumulator agradient[5];
  accumulator ahessian[15];

  for(unsigned n=nmin_; n<nmax_; n++)
  {
    const double xn = xc-C3_[n];
    const double xn2 = SQR(xn);
    const double log_pn = C1_[n] - C2_[n]*xn2;
    const double pdf_n = std::exp(log_pn);

    apdf.accumulate(pdf_n);

    const double dln0 = dC1_dF_[n];
    const double dln1 = C2_[n]*2.0*xn;
    const double dln2 = dC1_ds_[n] - dC2_ds_[n]*xn2;
    const double dln3 = dC1_dg_[n] - dC2_dg_[n]*xn2 + dC3_dg_[n]*xn;
    const double dln4 = dC1_db_[n] - dC2_db_[n]*xn2;

    agradient[0].accumulate(pdf_n*dln0);
    agradient[1].accumulate(pdf_n*dln1);
    agradient[2].accumulate(pdf_n*dln2);
    agradient[3].accumulate(pdf_n*dln3);
    agradient[4].accumulate(pdf_n*dln4);

    const double d2ln00 = d2C1_dF2_[n];

    const double d2ln11 = C2_[n]*-2.0;
    const double d2ln12 = dC2_ds_[n]*2.0*xn;
    const double d2ln13 = dC2_dg_[n]*2.0*xn - dC3_dg_[n];
    const double d2ln14 = dC2_db_[n]*2.0*xn;

    const double d2ln22 = d2C1_ds2_[n] - d2C2_ds2_[n]*xn2;
    const double d2ln23 = d2C1_dsdg_[n] - d2C2_dsdg_[n]*xn2 + d2C3_dsdg_[n]*xn;
    const double d2ln24 = d2C1_dsdb_[n] - d2C2_dsdb_[n]*xn2;

    // caution how extra terms are included in d2C1_dg2_ and dC3_dg2_
    const double d2ln33 = d2C1_dg2_[n] - d2C2_dg2_[n]*xn2 + d2C3_dg2_[n]*xn;
    const double d2ln34 = d2C1_dgdb_[n] - d2C2_dgdb_[n]*xn2 + d2C3_dgdb_[n]*xn;

    const double d2ln44 = d2C1_db2_[n] - d2C2_db2_[n]*xn2;

    ahessian[0].accumulate(pdf_n*(dln0*dln0 + d2ln00));
    ahessian[1].accumulate(pdf_n*dln0*dln1);
    ahessian[2].accumulate(pdf_n*dln0*dln2);
    ahessian[3].accumulate(pdf_n*dln0*dln3);
    ahessian[4].accumulate(pdf_n*dln0*dln4);

    ahessian[5].accumulate(pdf_n*(dln1*dln1 + d2ln11));
    ahessian[6].accumulate(pdf_n*(dln1*dln2 + d2ln12));
    ahessian[7].accumulate(pdf_n*(dln1*dln3 + d2ln13));
    ahessian[8].accumulate(pdf_n*(dln1*dln4 + d2ln14));

    ahessian[9].accumulate(pdf_n*(dln2*dln2 + d2ln22));
    ahessian[10].accumulate(pdf_n*(dln2*dln3 + d2ln23));
    ahessian[11].accumulate(pdf_n*(dln2*dln4 + d2ln24));

    ahessian[12].accumulate(pdf_n*(dln3*dln3 + d2ln33));
    ahessian[13].accumulate(pdf_n*(dln3*dln4 + d2ln34));

    ahessian[14].accumulate(pdf_n*(dln4*dln4 + d2ln44));
  }

  gradient(0) = c_gauss_norm * agradient[0].total();
  gradient(1) = c_gauss_norm * agradient[1].total();
  gradient(2) = c_gauss_norm * agradient[2].total();
  gradient(3) = c_gauss_norm * agradient[3].total();
  gradient(4) = c_gauss_norm * agradient[4].total();

  hessian(0,0) = c_gauss_norm * ahessian[0].total();
  hessian(0,1) = c_gauss_norm * ahessian[1].total();
  hessian(0,2) = c_gauss_norm * ahessian[2].total();
  hessian(0,3) = c_gauss_norm * ahessian[3].total();
  hessian(0,4) = c_gauss_norm * ahessian[4].total();

  hessian(1,0) = hessian(0,1);
  hessian(1,1) = c_gauss_norm * ahessian[5].total();
  hessian(1,2) = c_gauss_norm * ahessian[6].total();
  hessian(1,3) = c_gauss_norm * ahessian[7].total();
  hessian(1,4) = c_gauss_norm * ahessian[8].total();

  hessian(2,0) = hessian(0,2);
  hessian(2,1) = hessian(1,2);
  hessian(2,2) = c_gauss_norm * ahessian[9].total();
  hessian(2,3) = c_gauss_norm * ahessian[10].total();
  hessian(2,4) = c_gauss_norm * ahessian[11].total();

  hessian(3,0) = hessian(0,3);
  hessian(3,1) = hessian(1,3);
  hessian(3,2) = hessian(2,3);
  hessian(3,3) = c_gauss_norm * ahessian[12].total();
  hessian(3,4) = c_gauss_norm * ahessian[13].total();

  hessian(4,0) = hessian(0,4);
  hessian(4,1) = hessian(1,4);
  hessian(4,2) = hessian(2,4);
  hessian(4,3) = hessian(3,4);
  hessian(4,4) = c_gauss_norm * ahessian[14].total();

  return c_gauss_norm * apdf.total();
}

double PoissonGaussianMES::pdf_ped(double x)
{
  const double xc = x - ped_zero_dc_;
  const double log_pn = -(log_s_ + C2_[0]*SQR(xc));
  return c_gauss_norm * std::exp(log_pn);
}

double PoissonGaussianMES::
pdf_gradient_ped(double x, VecRef gradient)
{
  gradient.resize(5);

  const double xc = x - ped_zero_dc_;
  const double xc2 = SQR(xc);
  const double log_pn = -(log_s_ + C2_[0]*xc2);
  const double pdf = c_gauss_norm * std::exp(log_pn);

  const double dln1 = 2.0*C2_[0]*xc;
  const double dln2 = -(1.0/ped_rms_dc_ + dC2_ds_[0]*xc2);

  gradient[0] = 0;
  gradient[1] = pdf * dln1;
  gradient[2] = pdf * dln2;
  gradient[3] = 0;
  gradient[4] = 0;

  return pdf;
}

double PoissonGaussianMES::
pdf_gradient_hessian_ped(double x, VecRef gradient, MatRef hessian)
{
  gradient.resize(5);
  hessian.resize(5,5);

  if(!hessian_elements_good_)calc_cached_vars(true);

  const double xc = x - ped_zero_dc_;
  const double xc2 = SQR(xc);
  const double log_pn = -(log_s_ + C2_[0]*xc2);
  const double pdf = c_gauss_norm * std::exp(log_pn);

  const double dln1 = 2.0*C2_[0]*xc;
  const double dln2 = -(1.0/ped_rms_dc_ + dC2_ds_[0]*xc2);

  gradient(0) = 0;
  gradient(1) = pdf * dln1;
  gradient(2) = pdf * dln2;
  gradient(3) = 0;
  gradient(4) = 0;

  hessian(0,0) = hessian(0,1) = hessian(0,2) = hessian(0,3) = hessian(0,4) = 0;

  hessian(1,0) = 0;
  hessian(1,1) = pdf * (dln1*dln1 - 2.0*C2_[0]);
  hessian(1,2) = pdf * (dln1*dln2 + 2.0*dC2_ds_[0]*xc);
  hessian(1,3) = hessian(1,4) = 0;

  hessian(2,0) = 0;
  hessian(2,1) = hessian(1,2);
  hessian(2,2) = pdf * (dln2*dln2 + 1.0/SQR(ped_rms_dc_) - d2C2_ds2_[0]*xc2);
  hessian(2,3) = hessian(2,4) = 0;

  hessian(3,0) = hessian(3,1) = hessian(3,2) = hessian(3,3) = hessian(3,4) = 0;
  hessian(4,0) = hessian(4,1) = hessian(4,2) = hessian(4,3) = hessian(4,4) = 0;

  return pdf;
}

void PoissonGaussianMES::calc_cached_vars(bool calc_hessian)
{
  const double F = intensity_pe_;
  const double g = gain_dc_pe_;
  const double s = ped_rms_dc_;
  const double b = ses_rms_pe_;

  s2_ = SQR(s);
  b2_ = SQR(b);
  g2_ = SQR(g);
  log_s_ = std::log(s);

  const double var_n = b2_*g2_;
  const double var_s = s2_;

  calc_hessian |= force_calc_hessian_;
  if(calc_hessian && d2C1_dF2_.empty())
  {
    d2C1_dF2_.resize(nmax_);
    d2C1_ds2_.resize(nmax_);
    d2C1_dg2_.resize(nmax_);
    d2C1_db2_.resize(nmax_);
    d2C1_dsdg_.resize(nmax_);
    d2C1_dsdb_.resize(nmax_);
    d2C1_dgdb_.resize(nmax_);
    d2C2_ds2_.resize(nmax_);
    d2C2_dg2_.resize(nmax_);
    d2C2_db2_.resize(nmax_);
    d2C2_dsdg_.resize(nmax_);
    d2C2_dsdb_.resize(nmax_);
    d2C2_dgdb_.resize(nmax_);
    d2C3_dg2_.resize(nmax_);
    d2C3_dsdg_.resize(nmax_);
    d2C3_dgdb_.resize(nmax_);
  }

  for(unsigned n=0; n<nmax_; n++) // always loop from n=0 here
    {
      const double dbl_n = double(n);
      const double var = dbl_n*var_n + var_s;

#if 0
      accumulator acc;
      acc.accumulate(dbl_n*std::log(F));
      acc.accumulate(-lgamma(dbl_n+1.0));
      acc.accumulate(-F);
      acc.accumulate(-0.5*std::log(var));
      C1_[n] = acc.total();
#else
      C1_[n] = dbl_n*std::log(F) - F - lgamma(dbl_n+1.0) - 0.5*std::log(var);
#endif
      C2_[n] = 0.5/var;
      C3_[n] = dbl_n*g;

      dC1_dF_[n] = dbl_n/F - 1.0;
      dC1_ds_[n] = -s/var;
      dC1_dg_[n] = -dbl_n*b2_*g/var;
      dC1_db_[n] = -dbl_n*g2_*b/var;

      dC2_ds_[n] = dC1_ds_[n]/var;
      dC2_dg_[n] = dC1_dg_[n]/var;
      dC2_db_[n] = dC1_db_[n]/var;

      dC3_dg_[n] = C2_[n]*2.0*dbl_n;

      if(!calc_hessian)
      {
        hessian_elements_good_ = false;
        continue;
      }

      d2C1_dF2_[n]  = -dbl_n/SQR(F);
      d2C1_ds2_[n]  = -(1.0/var + 2.0*s*dC2_ds_[n]);
      d2C1_dg2_[n]  = -dbl_n*b2_*(1.0/var + 2.0*g*dC2_dg_[n]);
      d2C1_db2_[n]  = -dbl_n*g2_*(1.0/var + 2.0*b*dC2_db_[n]);
      d2C1_dsdg_[n] = -2.0*s*dC2_dg_[n];
      d2C1_dsdb_[n] = -2.0*s*dC2_db_[n];
      d2C1_dgdb_[n] = -2.0*dbl_n*b*g*(1.0/var + b*dC2_db_[n]);

      d2C2_ds2_[n]  = d2C1_ds2_[n]/var + 2.0*dC1_ds_[n]*dC2_ds_[n];
      d2C2_dg2_[n]  = d2C1_dg2_[n]/var + 2.0*dC1_dg_[n]*dC2_dg_[n];
      d2C2_db2_[n]  = d2C1_db2_[n]/var + 2.0*dC1_db_[n]*dC2_db_[n];
      d2C2_dsdg_[n] = d2C1_dsdg_[n]/var + 2.0*dC1_ds_[n]*dC2_dg_[n];
      d2C2_dsdb_[n] = d2C1_dsdb_[n]/var + 2.0*dC1_ds_[n]*dC2_db_[n];
      d2C2_dgdb_[n] = d2C1_dgdb_[n]/var + 2.0*dC1_dg_[n]*dC2_db_[n];

      d2C3_dg2_[n]  = dC2_dg_[n]*2.0*dbl_n;
      d2C3_dsdg_[n] = dC2_ds_[n]*2.0*dbl_n;
      d2C3_dgdb_[n] = dC2_db_[n]*2.0*dbl_n;

      // Add extra terms needed for "d2ln33"
      d2C3_dg2_[n] *= 2.0;
      d2C1_dg2_[n] -= dC3_dg_[n]*dbl_n;

      hessian_elements_good_ = true;
    }
}

// ============================================================================
//
// PoissonGaussianMES_HighAccuracy - PG model for low and high
// intensity regime using legacy code. May be more flexable in certain
// circumstances but also likely more brittle. Definitely slower. NOT
// RECOMMENDED!!
//
// ============================================================================

PoissonGaussianMES_HighAccuracy::
PoissonGaussianMES_HighAccuracy(double tol):
    tol_(tol)
{
  // nothing to see here
}

PoissonGaussianMES_HighAccuracy::~PoissonGaussianMES_HighAccuracy()
{
  // nothing to see here
}

unsigned PoissonGaussianMES_HighAccuracy::num_parameters()
{
  return 5;
}

std::vector<function::ParameterAxis>
PoissonGaussianMES_HighAccuracy::parameters()
{
  return { { "light_intensity", "PE", true, 0, false, 0 },
    { "ped_zero", "DC", false, 0, false, 0 },
    { "ped_width", "DC", true, 0, false, 0 },
    { "gain", "DC/PE", false, 0, false, 0 }, // can be negative
    { "ses_width", "PE", true, 0, false, 0 }
  };
}

Eigen::VectorXd PoissonGaussianMES_HighAccuracy::parameter_values()
{
  Eigen::VectorXd vec(5);
  vec << intensity_pe_, ped_zero_dc_, ped_rms_dc_,
      gain_dc_pe_, ses_rms_pe_;
  return vec;
}

void PoissonGaussianMES_HighAccuracy::
set_parameter_values(ConstVecRef values)
{
  verify_set_parameter_values(values, "PoissonGaussianMES_HighAccuracy");
  assign_parameters(values, intensity_pe_, ped_zero_dc_, ped_rms_dc_,
                    gain_dc_pe_, ses_rms_pe_);
}

bool PoissonGaussianMES_HighAccuracy::can_calculate_parameter_gradient()
{
  return true;
}

bool PoissonGaussianMES_HighAccuracy::can_calculate_parameter_hessian()
{
  return false;
}

double PoissonGaussianMES_HighAccuracy::
pdf_gradient_hessian_ped(double x, VecRef gradient, MatRef hessian)
{
  assert(0);
  return 0;
}

double PoissonGaussianMES_HighAccuracy::
pdf_gradient_hessian_mes(double x, VecRef gradient, MatRef hessian)
{
  assert(0);
  return 0;
}

double PoissonGaussianMES_HighAccuracy::pdf_mes(double x)
{
  return std::exp(PMTModelPG::logL_multi_pe(x, intensity_pe_, 1.0,
                                            SQR(ses_rms_pe_), gain_dc_pe_,
                                            SQR(ped_rms_dc_), ped_zero_dc_,
                                            tol_));
}

double PoissonGaussianMES_HighAccuracy::pdf_ped(double x)
{
  return std::exp(PMTModelPG::logL_ped(x, SQR(ped_rms_dc_), ped_zero_dc_));
}

double PoissonGaussianMES_HighAccuracy::
pdf_gradient_mes(double x, VecRef gradient)
{
  gradient.resize(5);

  PMTModelPGDerivs log_gradient;
  double p =
      std::exp(PMTModelPG::logL_multi_pe_derivs(log_gradient, x,
           intensity_pe_, 1.0, SQR(ses_rms_pe_), gain_dc_pe_,
           SQR(ped_rms_dc_), ped_zero_dc_, tol_));

  gradient(0) = p * log_gradient.dfdF;
  gradient(1) = p * log_gradient.dfdx0;
  gradient(2) = p * 2.0 * ped_rms_dc_ * log_gradient.dfds2;
  gradient(3) = p * log_gradient.dfdg;
  gradient(4) = p * 2.0 * ses_rms_pe_ * log_gradient.dfdb2;

  return p;
}

double PoissonGaussianMES_HighAccuracy::
pdf_gradient_ped(double x, VecRef gradient)
{
  PMTModelPGDerivs log_gradient;
  double p =
      std::exp(PMTModelPG::logL_ped_derivs(log_gradient, x,
                                           SQR(ped_rms_dc_), ped_zero_dc_));

  gradient[0] = 0; // p * log_gradient.dfdF;
  gradient[1] = p * log_gradient.dfdx0;
  gradient[2] = p * 2.0 * ped_rms_dc_ * log_gradient.dfds2;
  gradient[3] = 0; // p * log_gradient.dfdg;
  gradient[4] = 0; // p * 2.0 * ses_rms_pe_ * log_gradient.dfdb2;

  return p;
}
