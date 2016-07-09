/*
   calin/calib/spe_fit.cpp -- Stephen Fegan -- 2015-03-01

   Functions to do fit to multi-electron spectrum in the "single PE"
   domain.

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#include <iostream>
#include <cmath>
#include <limits>
#include <cassert>
#include <memory>

#include <fftw3.h>

#include <math/fftw_util.hpp>
#include <math/accumulator.hpp>
#include <calib/spe_fit.hpp>
#include <calib/pmt_model_pg.hpp>
#include <io/log.hpp>

using namespace calin::math;
using namespace calin::math::fftw_util;
using namespace calin::calib::spe_fit;
using namespace calin::calib::pmt_model_pg;
using namespace calin::io::log;

using calin::math::function::assign_parameters;

using uptr_fftw_plan = std::unique_ptr<fftw_plan_s,void(*)(fftw_plan_s*)>;
using uptr_fftw_data = std::unique_ptr<double,void(*)(void*)>;

namespace {

static double SQR(double x)
{
  return x*x;
}

constexpr double c_gauss_norm = 0.5*M_2_SQRTPI*M_SQRT1_2;

} // anonymous namespace

#if 0
SingleElectronSpectrum::~SingleElectronSpectrum()
{
  // nothing to see here
}
#endif

// ============================================================================
//
// MultiElectronSpectrum base class
//
// ============================================================================

MultiElectronSpectrum::~MultiElectronSpectrum()
{
  // nothing to see here
}


double MultiElectronSpectrum::
pdf_gradient_hessian_ped(double x, VecRef gradient, MatRef hessian)
{
  assert(0);
  return 0;
}

double MultiElectronSpectrum::
pdf_gradient_hessian_mes(double x, VecRef gradient, MatRef hessian)
{
  assert(0);
  return 0;
}

bool MultiElectronSpectrum::can_calculate_parameter_hessian()
{
  return false;
}

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
  //constexpr double phuge_val = std::numeric_limits<double>::max();
  //constexpr double nhuge_val = std::numeric_limits<double>::lowest();
  return { { "light_intensity", "PE", true, 0, false, 0 },
    { "ped_zero", "DC", false, 0, false, 0 },
    { "ped_width", "DC", true, tiny_val, false, 0 },
    { "gain", "DC/PE", false, 0, false, 0 }, // can be negative
    { "ses_width", "PE", true, tiny_val, false, 0 }
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
  assign_parameters(values, intensity_pe_, ped_zero_dc_, ped_rms_dc_,
                    gain_dc_pe_, ses_rms_pe_);
}

bool PoissonGaussianMES_HighAccuracy::can_calculate_parameter_gradient()
{
  return true;
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

// ============================================================================
//
// GeneralPoissonMES - Poisson model using generic SES and pedestal
// distribution. Uses FFTs to do convolution
//
// ============================================================================

GeneralPoissonMES::
GeneralPoissonMES(double x0, double dx, unsigned npoint,
                  SingleElectronSpectrum* ses,
                  math::pdf_1d::Parameterizable1DPDF* ped,
                  unsigned nmax, bool adopt_ses, bool adopt_ped):
    MultiElectronSpectrum(),
    ses_pdf_(ses), ped_pdf_(ped),
    adopt_ses_pdf_(adopt_ses), adopt_ped_pdf_(adopt_ped),
    nmax_(nmax), x0_(x0), dx_(dx), nsample_(npoint), nes_fft_(nmax)
{
  ped_spec_ = fftw_alloc_real(nsample_);
  ped_fft_ = fftw_alloc_real(nsample_);
  for(unsigned ines=0; ines<nmax_;ines++)
    nes_fft_[ines] = fftw_alloc_real(nsample_);
  mes_spec_ = fftw_alloc_real(nsample_);

  ped_plan_fwd_ =
      fftw_plan_r2r_1d(nsample_, ped_fft_, ped_fft_, FFTW_R2HC, 0);
  ses_plan_fwd_ =
      fftw_plan_r2r_1d(nsample_, nes_fft_[0], nes_fft_[0], FFTW_R2HC, 0);
  mes_plan_rev_ =
      fftw_plan_r2r_1d(nsample_, mes_spec_, mes_spec_, FFTW_HC2R, 0);

  if(can_calculate_parameter_gradient())
  {
    unsigned ped_npar = ped_pdf_->num_parameters();

    ped_grad_.resize(ped_npar);
    ped_grad_fft_.resize(ped_npar);
    ped_grad_plan_fwd_.resize(ped_npar);
    for(unsigned ipar=0;ipar<ped_npar;ipar++)
    {
      ped_grad_[ipar] = fftw_alloc_real(nsample_);
      ped_grad_fft_[ipar] = fftw_alloc_real(nsample_);
      ped_grad_plan_fwd_[ipar] =
          fftw_plan_r2r_1d(nsample_, ped_grad_[ipar], ped_grad_fft_[ipar],
                           FFTW_R2HC, 0);
    }

    unsigned ses_npar = ses_pdf_->num_parameters();
    ses_grad_fft_.resize(ses_npar);
    ses_grad_plan_fwd_.resize(ses_npar);
    for(unsigned ipar=0;ipar<ses_npar;ipar++)
    {
      ses_grad_fft_[ipar] = fftw_alloc_real(nsample_);
      ses_grad_plan_fwd_[ipar] =
          fftw_plan_r2r_1d(nsample_, ses_grad_fft_[ipar], ses_grad_fft_[ipar],
                           FFTW_R2HC, 0);
    }

    unsigned mes_npar = 1+ses_npar+ped_npar;
    mes_grad_.resize(mes_npar);
    mes_grad_plan_rev_.resize(mes_npar);
    for(unsigned ipar=0;ipar<mes_npar;ipar++)
    {
      mes_grad_[ipar] = fftw_alloc_real(nsample_);
      mes_grad_plan_rev_[ipar] =
          fftw_plan_r2r_1d(nsample_, mes_grad_[ipar], mes_grad_[ipar],
                           FFTW_HC2R, 0);
    }
  }
}

GeneralPoissonMES::~GeneralPoissonMES()
{
  for(auto x : mes_grad_plan_rev_)fftw_destroy_plan(x);
  for(auto x : ses_grad_plan_fwd_)fftw_destroy_plan(x);
  for(auto x : ped_grad_plan_fwd_)fftw_destroy_plan(x);
  for(auto x : mes_grad_)fftw_free(x);
  for(auto x : ped_grad_)fftw_free(x);
  for(auto x : ses_grad_fft_)fftw_free(x);
  for(auto x : ped_grad_fft_)fftw_free(x);

  fftw_destroy_plan(mes_plan_rev_);
  fftw_destroy_plan(ses_plan_fwd_);
  fftw_destroy_plan(ped_plan_fwd_);

  fftw_free(ped_fft_);
  fftw_free(ped_spec_);
  for(auto x : nes_fft_)fftw_free(x);
  fftw_free(mes_spec_);

  if(adopt_ses_pdf_)delete ses_pdf_;
  if(adopt_ped_pdf_)delete ped_pdf_;
}

unsigned GeneralPoissonMES::num_parameters()
{
  return 1+ses_pdf_->num_parameters()+ped_pdf_->num_parameters();
}

auto GeneralPoissonMES::parameters() ->
    std::vector<math::function::ParameterAxis>
{
  std::vector<math::function::ParameterAxis> pvec;
  pvec.push_back({ "light_intensity", "PE", true, 0, false, 0 });
  std::vector<math::function::ParameterAxis> pped { ped_pdf_->parameters() };
  for(auto& ip : pped)ip.name = std::string("ped.") + ip.name;
  pvec.insert(pvec.end(), pped.begin(), pped.end());
  std::vector<math::function::ParameterAxis> pses { ses_pdf_->parameters() };
  for(auto& ip : pses)ip.name = std::string("ses.") + ip.name;
  pvec.insert(pvec.end(), pses.begin(), pses.end());
  return pvec;
}

Eigen::VectorXd GeneralPoissonMES::parameter_values()
{
  Eigen::VectorXd param(num_parameters());
  param[0] = intensity_pe_;
  unsigned ip = 1;
  unsigned num_ped_params = ped_pdf_->num_parameters();
  param.segment(ip,num_ped_params) = ped_pdf_->parameter_values();
  ip += num_ped_params;
  param.segment(ip,ses_pdf_->num_parameters()) = ses_pdf_->parameter_values();
  return param;
}

void GeneralPoissonMES::set_parameter_values(ConstVecRef values)
{
  assign_parameters(values, intensity_pe_);
  unsigned ip = 1;
  unsigned num_ped_params = ped_pdf_->num_parameters();
  ped_pdf_->set_parameter_values(values.segment(ip,num_ped_params));
  ip += num_ped_params;
  ses_pdf_->set_parameter_values(values.segment(ip,ses_pdf_->num_parameters()));
  set_cache();
}

bool GeneralPoissonMES::can_calculate_parameter_gradient()
{
  //return false;
  return ped_pdf_->can_calculate_parameter_gradient() &&
      ses_pdf_->can_calculate_parameter_gradient();
}

bool GeneralPoissonMES::can_calculate_parameter_hessian()
{
  return false; // for the moment we are lazy
}

double GeneralPoissonMES::pdf_ped(double x)
{
  return std::max(ped_spec_[ibin(x)],0.0);
}

double GeneralPoissonMES::pdf_gradient_ped(double x, VecRef gradient)
{
  assert(can_calculate_parameter_gradient());
  int thebin = ibin(x);
  //std::cout << "HELLO : " << x << ' ' << thebin << '\n';
  unsigned npar = ped_pdf_->num_parameters();
  gradient.resize(num_parameters());
  gradient.setZero();
  for(unsigned ipar=0;ipar<npar;ipar++)
    gradient[1+ipar] = ped_grad_[ipar][thebin];
  //std::cout << gradient << ' ' << std::max(ped_spec_[thebin],0.0) << '\n';
  return std::max(ped_spec_[thebin],0.0);
}

double GeneralPoissonMES::pdf_gradient_hessian_ped(double x, VecRef gradient,
                                                   MatRef hessian)
{
  throw std::logic_error("GeneralPoissonMES: cannot calculate ped hessian");
}

double GeneralPoissonMES::pdf_mes(double x)
{
  //std::cout << x << ' ' << ibin(x) << ' ' << mes_spec_[ibin(x)] << '\n';
  return std::max(mes_spec_[ibin(x)],0.0);
}

double GeneralPoissonMES::pdf_gradient_mes(double x, VecRef gradient)
{
  assert(can_calculate_parameter_gradient());
  int thebin = ibin(x);
  unsigned npar = num_parameters();
  gradient.resize(npar);
  for(unsigned ipar=0;ipar<npar;ipar++)
    gradient[ipar] = mes_grad_[ipar][thebin];
  return std::max(mes_spec_[thebin],0.0);
}

double GeneralPoissonMES::pdf_gradient_hessian_mes(double x, VecRef gradient,
                                                   MatRef hessian)
{
  throw std::logic_error("GeneralPoissonMES: cannot calculate ses hessian");
}

double GeneralPoissonMES::ped_rms_dc()
{
  assert(false);
  return 0;
}

double GeneralPoissonMES::ped_zero_dc()
{
  assert(false);
  return 0;
}

double GeneralPoissonMES::ses_mean_dc()
{
  assert(false);
  return 0;
}

double GeneralPoissonMES::ses_rms_pe()
{
  assert(false);
  return 0;
}

std::vector<double> GeneralPoissonMES::multi_electron_spectrum() const
{
  std::vector<double> spec(nsample_);
  std::copy(mes_spec_, mes_spec_+nsample_, spec.begin());
  return spec;
}

std::vector<double> GeneralPoissonMES::pedestal_spectrum() const
{
  std::vector<double> spec(nsample_);
  std::copy(ped_spec_, ped_spec_+nsample_, spec.begin());
  return spec;
}

std::vector<double> GeneralPoissonMES::
multi_electron_spectrum_gradient(unsigned iparam) const
{
  if(iparam >= const_cast<GeneralPoissonMES*>(this)->num_parameters())
    throw std::out_of_range(
        "GeneralPoissonMES::multi_electron_spectrum_gradient: "
        "iparam out of range");

  return std::vector<double>(mes_grad_[iparam], mes_grad_[iparam]+nsample_);
}

std::vector<double> GeneralPoissonMES::
pedestal_spectrum_gradient(unsigned iparam) const
{
  if(iparam >= const_cast<GeneralPoissonMES*>(this)->num_parameters())
    throw std::out_of_range("GeneralPoissonMES::pedestal_spectrum_gradient: "
                            "iparam out of range");

  if(iparam==0 || iparam>ped_pdf_->num_parameters())
    return std::vector<double>(nsample_, 0.0);
  else
    return std::vector<double>(ped_grad_[iparam], ped_grad_[iparam]+nsample_);
}

std::vector<double> GeneralPoissonMES::
single_electron_spectrum_gradient(unsigned iparam) const
{
  if(iparam >= const_cast<GeneralPoissonMES*>(this)->num_parameters())
    throw std::out_of_range("GeneralPoissonMES::pedestal_spectrum_gradient: "
                            "iparam out of range");

  unsigned npedpar = ped_pdf_->num_parameters();
  if(iparam<=npedpar)return std::vector<double>(nsample_, 0.0);
  iparam -= 1+npedpar;

  uptr_fftw_data spec_buffer { fftw_alloc_real(nsample_), fftw_free };
  assert(spec_buffer);

  uptr_fftw_plan spec_plan = {
    fftw_plan_r2r_1d(nsample_, spec_buffer.get(), spec_buffer.get(),
                     FFTW_HC2R, 0), fftw_destroy_plan };
  assert(spec_plan);

  std::copy(ses_grad_fft_[iparam], ses_grad_fft_[iparam]+nsample_,
            spec_buffer.get());
  fftw_execute(spec_plan.get());
  std::vector<double> spec_gradient(nsample_);
  double norm { 1.0/double(nsample_) };
  std::transform(spec_buffer.get(), spec_buffer.get()+nsample_,
                 spec_gradient.begin(),
                 [norm](double x){return x*norm;});
  return spec_gradient;
}

std::vector<double> GeneralPoissonMES::n_electron_spectrum(unsigned n) const
{
  if(n==0 or n>nmax_)
    throw std::out_of_range("GeneralPoissonMES::n_electron_spectrum: "
                            "number of PEs out of range");

  uptr_fftw_data spec_buffer { fftw_alloc_real(nsample_), fftw_free };
  assert(spec_buffer);

  uptr_fftw_plan spec_plan = {
    fftw_plan_r2r_1d(nsample_, spec_buffer.get(), spec_buffer.get(),
                     FFTW_HC2R, 0), fftw_destroy_plan };
  assert(spec_plan);

  std::copy(nes_fft_[n-1], nes_fft_[n-1]+nsample_, spec_buffer.get());
  fftw_execute(spec_plan.get());
  std::vector<double> spec(nsample_);
  double norm { 1.0/double(nsample_) };
  std::transform(spec_buffer.get(), spec_buffer.get()+nsample_, spec.begin(),
                 [norm](double x){return x*norm;});
  return spec;
}

std::vector<double> GeneralPoissonMES::mes_n_electron_cpt(unsigned n) const
{
  if(n>nmax_)
    throw std::out_of_range("GeneralPoissonMES::mes_n_electron_cpt: "
                            "number of PEs out of range");

  std::vector<double> spec(nsample_);
  double log_nsample = std::log(double(nsample_));

  if(n==0)
  {
    double scale = std::exp(-intensity_pe_);
    std::transform(ped_spec_, ped_spec_+nsample_, spec.begin(),
                   [scale](const double& x){ return x*scale; });
    return spec;
  }

  uptr_fftw_data spec_buffer { fftw_alloc_real(nsample_), fftw_free };
  assert(spec_buffer);

  uptr_fftw_plan spec_plan = {
    fftw_plan_r2r_1d(nsample_, spec_buffer.get(), spec_buffer.get(),
                     FFTW_HC2R, 0), fftw_destroy_plan };
  assert(spec_plan);

  double log_intensity = std::log(intensity_pe_);
  double dbl_n { double(n) };
  double poisson_factor {
    std::exp(dbl_n*log_intensity - intensity_pe_ - lgamma(dbl_n+1.0) -
             log_nsample) };

  hcvec_scale_and_multiply(spec_buffer.get(), nes_fft_[n-1], ped_fft_,
                           nsample_, poisson_factor*dx_);

  fftw_execute(spec_plan.get());
  std::copy(spec_buffer.get(), spec_buffer.get()+nsample_, spec.begin());
  return spec;
}

Eigen::VectorXd GeneralPoissonMES::
extract_ped_gradient_values(ConstVecRef gradient)
{
  return gradient.segment(1, ped_pdf_->num_parameters());
}

Eigen::VectorXd GeneralPoissonMES::
extract_ses_gradient_values(ConstVecRef gradient)
{
  return gradient.segment(1+ped_pdf_->num_parameters(),
                          ses_pdf_->num_parameters());
}

Eigen::MatrixXd GeneralPoissonMES::
extract_ped_hessian_values(ConstMatRef hessian)
{
  unsigned ped_npar = ped_pdf_->num_parameters();
  return hessian.block(1, 1, ped_npar, ped_npar);
}

Eigen::MatrixXd GeneralPoissonMES::
extract_ses_hessian_values(ConstMatRef hessian)
{
  unsigned ped_npar = ped_pdf_->num_parameters();
  unsigned ses_npar = ses_pdf_->num_parameters();
  return hessian.block(1+ped_npar, 1+ped_npar, ses_npar, ses_npar);
}

int GeneralPoissonMES::ibin(double x) const
{
  int thebin = std::floor((x-x0_)/dx_);
  if(thebin<0 or thebin>=(int)nsample_)
  {
    std::ostringstream str;
    str << "GeneralPoissonMES::ibin: x=" << x
        << " out of range (bin=" << thebin << ", nsample=" << nsample_ << ")";
    throw std::out_of_range(str.str());
  }
  return thebin;
}

void GeneralPoissonMES::set_cache()
{
  // THIS FUNCTION IS TOO LONG

  bool calc_gradient = can_calculate_parameter_gradient();
  unsigned ses_npar = ses_pdf_->num_parameters();
  Eigen::VectorXd ses_gradient(ses_npar);
  ses_gradient.setZero();
  calin::math::accumulator::LikelihoodAccumulator ses_acc;
  for(unsigned isample = 0;isample<nsample_;isample++)
  {
    // The SES describes charge delivered by the PMT and is assumed be
    // positive. Hence the function is sampled from zero here and not
    // from x0. We could allow negative charges by using the domain
    // axis limits to specify the minimum bound. We would need to
    // adjust x0 to compensate for where the pedestal distribution
    // would end up.
    const double x = ses_x(isample); // (0.5+double(isample))*dx_;
    double val;
    if(calc_gradient)
    {
      val = ses_pdf_->value_and_parameter_gradient_1d(x, ses_gradient);
      if(!std::isfinite(val)){ val = 0; ses_gradient.setZero(); }
      for(unsigned ipar=0;ipar<ses_npar;ipar++)
        ses_grad_fft_[ipar][isample] = ses_gradient[ipar];
    }
    else
    {
      val = ses_pdf_->value_1d(x);
      if(!std::isfinite(val))val = 0;
    }
    nes_fft_[0][isample] = val;
    ses_acc.accumulate(val);
  }

  if(std::abs(ses_acc.total() * dx_ - 1.0) > 1.0/double(nsample_)
     and n_ses_norm_warning_ < n_max_ses_norm_warning_)
  {
    LOG(WARNING) << "SES normalization is significantly different from 1.0: "
                 << ses_acc.total() * dx_ << '\n'
                 << "SES parameter values : "
                 << ses_pdf_->parameter_values().transpose();
    n_ses_norm_warning_++;
    if(n_ses_norm_warning_ == n_max_ses_norm_warning_)
      LOG(INFO) << "Further SES normalization warnings will be suppressed!";
  }

  fftw_execute(ses_plan_fwd_);
  for(unsigned ines=1;ines<nmax_;ines++)
    hcvec_scale_and_multiply(nes_fft_[ines], nes_fft_[ines-1], nes_fft_[0],
                             nsample_, dx_);

  unsigned ped_npar = ped_pdf_->num_parameters();
  Eigen::VectorXd ped_gradient(ped_npar);
  ped_gradient.setZero();
  for(unsigned isample = 0;isample<nsample_;isample++)
  {
    const double x = ped_x(isample); // x0_ + (0.5+double(isample))*dx_;
    double val;
    if(calc_gradient)
    {
      val = ped_pdf_->value_and_parameter_gradient_1d(x, ped_gradient);
      if(!std::isfinite(val)){ val = 0; ped_gradient.setZero(); }
      for(unsigned ipar=0;ipar<ped_npar;ipar++)
        ped_grad_[ipar][isample] = ped_gradient[ipar];
    }
    else
    {
      val = ped_pdf_->value_1d(x);
      if(!std::isfinite(val))val = 0;
    }
    ped_spec_[isample] = ped_fft_[isample] = val;
  }
  fftw_execute(ped_plan_fwd_);

  unsigned mes_npar = 1+ped_npar+ses_npar;
  if(calc_gradient)
  {
    for(unsigned ipar=0;ipar<ses_npar;ipar++)
      fftw_execute(ses_grad_plan_fwd_[ipar]);
    for(unsigned ipar=0;ipar<ped_npar;ipar++)
      fftw_execute(ped_grad_plan_fwd_[ipar]);
    for(unsigned ipar=0;ipar<mes_npar;ipar++)
      std::fill(mes_grad_[ipar], mes_grad_[ipar]+nsample_, 0.0);
  }

  std::fill(mes_spec_, mes_spec_+nsample_, 0.0);

  double log_intensity = std::log(intensity_pe_);
  double log_nsample = std::log(double(nsample_));
  for(unsigned ines = 0;ines<nmax_;ines++)
  {
    double dbl_n { double(ines+1) };
    double poisson_factor {
      std::exp(dbl_n*log_intensity - intensity_pe_ - lgamma(dbl_n+1.0) -
               log_nsample) };
    hcvec_scale_and_add(mes_spec_, nes_fft_[ines], nsample_, poisson_factor);

    if(calc_gradient)
    {
      hcvec_scale_and_add(mes_grad_[0], nes_fft_[ines], nsample_,
                          poisson_factor*(dbl_n/intensity_pe_ - 1.0));
      if(ines>0)
        for(unsigned ipar=0;ipar<ses_npar;ipar++)
          hcvec_scale_and_add(mes_grad_[1+ped_npar+ipar], nes_fft_[ines-1],
                              nsample_, dbl_n*poisson_factor);
    }
  }

  if(calc_gradient)
  {
    hcvec_scale_and_multiply(mes_grad_[0], mes_grad_[0], ped_fft_, nsample_,
                        dx_);
    hcvec_scale_and_add(mes_grad_[0], ped_fft_, nsample_,
                        -std::exp(-intensity_pe_-log_nsample));
    for(unsigned ipar=0;ipar<ped_npar;ipar++)
    {
      hcvec_scale_and_multiply(mes_grad_[1+ipar], mes_spec_,
        ped_grad_fft_[ipar], nsample_, dx_);
      hcvec_scale_and_add(mes_grad_[1+ipar], ped_grad_fft_[ipar], nsample_,
        std::exp(-intensity_pe_-log_nsample));
    }
    for(unsigned ipar=0;ipar<ses_npar;ipar++)
    {
      hcvec_scale_and_multiply(mes_grad_[1+ped_npar+ipar],
        mes_grad_[1+ped_npar+ipar], ses_grad_fft_[ipar], nsample_, dx_);
      hcvec_scale_and_add(mes_grad_[1+ped_npar+ipar], ses_grad_fft_[ipar],
        nsample_, std::exp(log_intensity -intensity_pe_ - log_nsample));
      hcvec_scale_and_multiply(mes_grad_[1+ped_npar+ipar],
        mes_grad_[1+ped_npar+ipar], ped_fft_, nsample_, dx_);
    }
  }

  hcvec_scale_and_multiply(mes_spec_, mes_spec_, ped_fft_, nsample_, dx_);
  hcvec_scale_and_add(mes_spec_, ped_fft_, nsample_,
                      std::exp(-intensity_pe_-log_nsample));
  fftw_execute(mes_plan_rev_);

  if(calc_gradient)
    for(unsigned ipar=0;ipar<mes_npar;ipar++)
      fftw_execute(mes_grad_plan_rev_[ipar]);
}

// ============================================================================
//
// SPELikelihood -- calculate likelihood for optimizer
//
// ============================================================================

SPELikelihood::SPELikelihood(MultiElectronSpectrum& mes_model,
                             const math::histogram::SimpleHist& mes_data):
    function::MultiAxisFunction(), mes_model_(&mes_model),
    npar_(mes_model.num_parameters()),
    mes_data_(mes_data), has_ped_data_(false), ped_data_(1.0)
{
  // nothing to see here
}

SPELikelihood::SPELikelihood(MultiElectronSpectrum& mes_model,
                             const math::histogram::SimpleHist& mes_data,
                             const math::histogram::SimpleHist& ped_data):
    function::MultiAxisFunction(), mes_model_(&mes_model),
    npar_(mes_model.num_parameters()),
    mes_data_(mes_data), has_ped_data_(true), ped_data_(ped_data)
{
  // nothing to see here
}

SPELikelihood::~SPELikelihood()
{
  // nothing to see here
}

unsigned SPELikelihood::num_domain_axes()
{
  return mes_model_->num_parameters();
}

auto SPELikelihood::domain_axes() -> std::vector<math::function::DomainAxis>
{
  return mes_model_->parameters();
}

double SPELikelihood::value(ConstVecRef x)
{
  mes_model_->set_parameter_values(x);
  math::accumulator::LikelihoodAccumulator acc;
  for(auto& ibin : mes_data_)
    if(ibin.weight())
    {
      double pdf = mes_model_->pdf_mes(ibin.xval_center());
      if(pdf<=0)continue;
      acc.accumulate(std::log(pdf)*ibin.weight());
    }

  if(has_ped_data_)
    for(auto& ibin : ped_data_)
      if(ibin.weight())
      {
        double pdf = mes_model_->pdf_ped(ibin.xval_center());
        if(pdf<=0)continue;
        acc.accumulate(std::log(pdf)*ibin.weight());
      }
  return -acc.total();
}

bool SPELikelihood::can_calculate_gradient()
{
  return mes_model_->can_calculate_parameter_gradient();
}

bool SPELikelihood::can_calculate_hessian()
{
  return mes_model_->can_calculate_parameter_hessian();
}

double SPELikelihood::value_and_gradient(ConstVecRef x, VecRef gradient)
{
  mes_model_->set_parameter_values(x);
  gradient.resize(npar_);
  math::accumulator::LikelihoodAccumulator acc;
  std::vector<math::accumulator::LikelihoodAccumulator> gradient_acc(npar_);
  for(auto& ibin : mes_data_)
    if(ibin.weight())
    {
      double pdf = mes_model_->pdf_gradient_mes(ibin.xval_center(), gradient);
      if(pdf<=0)continue;
      acc.accumulate(std::log(pdf)*ibin.weight());
      for(unsigned ipar=0;ipar<npar_;ipar++)
        gradient_acc[ipar].accumulate(gradient(ipar)/pdf*ibin.weight());
    }

  if(has_ped_data_)
    for(auto& ibin : ped_data_)
      if(ibin.weight())
      {
        double pdf = mes_model_->pdf_gradient_ped(ibin.xval_center(), gradient);
        if(pdf<=0)continue;
        acc.accumulate(std::log(pdf)*ibin.weight());
        for(unsigned ipar=0;ipar<npar_;ipar++)
          gradient_acc[ipar].accumulate(gradient(ipar)/pdf*ibin.weight());
      }

  for(unsigned ipar=0;ipar<npar_;ipar++)
    gradient(ipar) = -gradient_acc[ipar].total();
  return -acc.total();
}

double SPELikelihood::
value_gradient_and_hessian(ConstVecRef x, VecRef gradient, MatRef hessian)
{
  mes_model_->set_parameter_values(x);
  gradient.resize(npar_);
  hessian.resize(npar_,npar_);
  math::accumulator::LikelihoodAccumulator acc;
  std::vector<math::accumulator::LikelihoodAccumulator> gradient_acc(npar_);
  std::vector<math::accumulator::LikelihoodAccumulator>
      hessian_acc(npar_*(npar_+1)/2);
  for(auto& ibin : mes_data_)
    if(ibin.weight())
    {
      double pdf =
          mes_model_->pdf_gradient_hessian_mes(ibin.xval_center(),
                                               gradient, hessian);
      if(pdf<=0)continue;
      acc.accumulate(std::log(pdf)*ibin.weight());
      for(unsigned ipar=0;ipar<npar_;ipar++)
        gradient_acc[ipar].accumulate(gradient(ipar)/pdf*ibin.weight());
      unsigned itri = 0;
      for(unsigned icol=0;icol<npar_;icol++)
        for(unsigned irow=icol;irow<npar_;irow++)
        {
          double summand = (hessian(icol,irow)
                            - gradient[icol]*gradient[irow]/pdf)/pdf;
          hessian_acc[itri++].accumulate(summand*ibin.weight());
        }
    }

  if(has_ped_data_)
    for(auto& ibin : ped_data_)
      if(ibin.weight())
      {
        double pdf =
            mes_model_->pdf_gradient_hessian_ped(ibin.xval_center(),
                                                 gradient, hessian);
        if(pdf<=0)continue;
        acc.accumulate(std::log(pdf)*ibin.weight());
        for(unsigned ipar=0;ipar<npar_;ipar++)
          gradient_acc[ipar].accumulate(gradient(ipar)/pdf*ibin.weight());
        unsigned iacc = 0;
        for(unsigned icol=0;icol<npar_;icol++)
          for(unsigned irow=icol;irow<npar_;irow++)
          {
            double summand = (hessian(icol,irow)
                            - gradient[icol]*gradient[irow]/pdf)/pdf;
            hessian_acc[iacc++].accumulate(summand*ibin.weight());
          }
      }

  for(unsigned ipar=0;ipar<npar_;ipar++)
    gradient(ipar) = -gradient_acc[ipar].total();

  unsigned iacc = 0;
  for(unsigned icol=0;icol<npar_;icol++)
  {
    for(unsigned irow=icol;irow<npar_;irow++)
      hessian(icol,irow) = -hessian_acc[iacc++].total();
    for(unsigned irow=0;irow<icol;irow++)
      hessian(icol,irow) = hessian(irow,icol);
  }

  return -acc.total();
}
