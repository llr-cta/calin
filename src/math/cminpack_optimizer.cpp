/*

   calin/math/cminpack_optimizer.cpp -- Stephen Fegan -- 2015-05-21

   Interface to ANL CMinpack optimizer suite

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
#include <iomanip>
#include <cstring>

#include <cminpack/cminpack.h>

#include <io/log.hpp>
#include <math/cminpack_optimizer.hpp>
#include <math/hessian.hpp>

using namespace calin::math::optimizer;
using namespace calin::io::log;

CMinpackOptimizer::
CMinpackOptimizer(function::MultiAxisFunction* fcn, bool adopt_fcn):
    Optimizer(fcn, adopt_fcn)
{
  // nothing to see here
}

CMinpackOptimizer::~CMinpackOptimizer()
{
  // nothing to see here
}

bool CMinpackOptimizer::is_local_optimizer()
{
  return true;
}

bool CMinpackOptimizer::requires_gradient()
{
  return true;
}

bool CMinpackOptimizer::requires_hessian()
{
  return true; // mode_==CMinpackMode::REQUIRE_HESSIAN;
}

bool CMinpackOptimizer::requires_box_constraints()
{
  return false;
}

bool CMinpackOptimizer::can_estimate_error()
{
  return true;
}

bool CMinpackOptimizer::can_use_gradient()
{
  return true;
}

bool CMinpackOptimizer::can_use_hessian()
{
  return true; // mode_!=CMinpackMode::DONT_USE_HESSIAN
}

bool CMinpackOptimizer::can_impose_box_constraints()
{
  return false;
}

OptimizationStatus CMinpackOptimizer::minimize(VecRef xopt, double& fopt)
{
  const double xtol { std::sqrt(std::numeric_limits<double>::epsilon()) };
  const int maxfev = max_iterations()>0?max_iterations():0x7FFFFFFF;

  opt_starting("CMinpack", requires_gradient(),
               requires_hessian(), {}, {}, {});

  unsigned naxes = fcn_->num_domain_axes();
  xopt = xbest_;
  Eigen::VectorXd gvec(naxes);
  qmat_.resize(naxes,naxes);
  rvec_.resize(naxes*(naxes+1)/2);
  qtfvec_.resize(naxes);
  dxvec_.resize(naxes);

  Eigen::VectorXd diag(naxes);
  int ngrad { 0 };
  int nhess { 0 };

  Eigen::VectorXd wa1(naxes);
  Eigen::VectorXd wa2(naxes);
  Eigen::VectorXd wa3(naxes);
  Eigen::VectorXd wa4(naxes);

  int istat = hybrj(cminpack_callback_hess,     // fcn
                    this,                       // void* p
                    naxes,                      // int n
                    xopt.data(),                // double* x
                    gvec.data(),                // double* fvec
                    qmat_.data(),               // double& fjac
                    naxes,                      // int ldfjac
                    xtol,                       // double xtol
                    maxfev,                     // int maxfev
                    diag.data(),                // double* diag,
                    1,                          // int mode
                    100.0,                      // double factor
                    1,                          // int nprint
                    &ngrad,                     // int* nfev
                    &nhess,                     // int& njev
                    rvec_.data(),               // double* r
                    rvec_.size(),               // int lr
                    qtfvec_.data(),             // double* qtf
                    wa1.data(),                 // double* wa1
                    wa2.data(),                 // double* wa2
                    wa3.data(),                 // double* wa3
                    wa4.data());                // double* wa4

  gvec_ = gvec;

  switch(istat)
  {
    case -3:
      opt_status_ = OptimizationStatus::STOPPED_AT_MAXTIME;
      opt_message_ = "Time limit reached";
      break;
    case -2:
    case -1:
      opt_status_ = OptimizationStatus::TOLERANCE_REACHED;
      opt_message_ = "Tolerance reached";
      break;
    case 0:
      opt_status_ = OptimizationStatus::OPTIMIZER_FAILURE;
      opt_message_ = "Invalid arguments";
      break;
    case 1:
      opt_status_ = OptimizationStatus::TOLERANCE_REACHED;
      opt_message_ = "Tolerance reached";
      break;
    case 2:
      opt_status_ = OptimizationStatus::STOPPED_AT_MAXCALLS;
      opt_message_ = "Calls limit reached";
      break;
    case 3:
      opt_status_ = OptimizationStatus::LIMITED_BY_PRECISION;
      opt_message_ = "Minimization limited by machine precision";
      break;
    case 4:
    case 5:
      opt_status_ = OptimizationStatus::OPTIMIZER_FAILURE;
      opt_message_ = "Optimizer not making progress";
      break;
    default:
      opt_status_ = OptimizationStatus::OPTIMIZER_FAILURE;
      opt_message_ = "This should not happen";
      assert(0);
  }

  double f_edm = edm(naxes);

  fopt = fbest_;
  xopt = xbest_;

  opt_finished(opt_status_, fbest_, xbest_, &f_edm);

  return opt_status_;
}

ErrorMatrixStatus CMinpackOptimizer::error_matrix_estimate(MatRef error_matrix)
{
  // Use the CMinpack QR factorization of the Hessian to calculate the
  // error matrix: Sigma = ( Q R )^-1. We solve for each column the
  // equation: Q^-1 = R Sigma, noting that Q is orthogonal: Q^-1 = Q^T
  const double scale = 2.0*fcn_->error_up();
  unsigned n = fcn_->num_domain_axes();
  error_matrix.resize(n,n);
  for(unsigned ivec = 0; ivec<n; ++ivec)
  {
    for(unsigned iir = 0; iir<n; ++iir)
    {
      unsigned ir = n-iir-1;
      const double* rr = rvec_.data() + n*(n+1)/2 - iir*(iir+1)/2 - n;
      double sum = qmat_(ivec,ir);
      for(unsigned ic = ir+1; ic<n; ++ic)sum -= rr[ic]*error_matrix(ic,ivec);
      error_matrix(ir,ivec) = sum/rr[ir];
    }
  }

  // The calculated "error_matrix" may not be symmetric as the final rank-1
  // updates in CMinpack do not assume a symmetric Hessian (rather a more
  // generic Jacbian of the N functions). So we symmetrize the matrix here.
  // We don't worry about enforcing positive-definiteness.
  Eigen::MatrixXd ems = (0.5*scale)*(error_matrix + error_matrix.transpose());
  error_matrix = ems;

  return ErrorMatrixStatus::GOOD;
}

ErrorMatrixStatus CMinpackOptimizer::
calc_error_matrix_and_eigenvectors(MatRef error_matrix,
                                   VecRef eigenvalues, MatRef eigenvectors)
{
  // Duplicate of the code from nlopt_optimizer. Move to base class?
  const unsigned npar { fcn_->num_domain_axes() };
  error_matrix.resize(npar,npar);
  Eigen::VectorXd error_hint;
  if(error_matrix_estimate(error_matrix) != ErrorMatrixStatus::UNAVAILABLE)
    error_hint = error_matrix.diagonal().array().sqrt();
  Eigen::MatrixXd hessian(npar,npar);
  hessian::calculate_hessian(*fcn_, xbest_, hessian, error_hint);
  return hessian::hessian_to_error_matrix(*fcn_, hessian, error_matrix,
                                          eigenvalues, eigenvectors);
}

int CMinpackOptimizer::
cminpack_callback_grad(void* self, int n, const double* x,
                       double* grad, int iflag)
{
  CMinpackOptimizer* that = static_cast<CMinpackOptimizer*>(self);
  return that->eval_func(n,x,grad,nullptr,iflag);
}

int CMinpackOptimizer::
cminpack_callback_hess(void *self, int n, const double* x,
                       double* grad, double* hess, int ldfjac, int iflag)
{
  CMinpackOptimizer* that = static_cast<CMinpackOptimizer*>(self);
  (void)ldfjac;
  return that->eval_func(n,x,grad,hess,iflag);
}

int CMinpackOptimizer::
eval_func(unsigned n, const double* x, double* grad, double* hess, int iflag)
{
  xvec_ = Eigen::Map<const Eigen::VectorXd>(x,n);
  if(iflag == 2)
  {
    fval_ = fcn_->value_gradient_and_hessian(xvec_,gvec_,hmat_);
    this->opt_progress(fval_, xvec_, &gvec_, &hmat_);
    Eigen::Map<Eigen::MatrixXd>(hess,n,n) = hmat_;
    return iflag;
  }
  else if(iflag == 1)
  {
    fval_ = fcn_->value_and_gradient(xvec_,gvec_);
    this->opt_progress(fval_, xvec_, &gvec_, nullptr);
    Eigen::Map<Eigen::VectorXd>(grad,n) = gvec_;
    return iflag;
  }
  else if(iflag == 0)
  {
    if(abs_tolerance()>0 or rel_tolerance()>0)
    {
      double f_edm = edm(n);
      if(std::abs(f_edm)<abs_tolerance())return -1;
      if(std::abs(fval_*f_edm)<rel_tolerance())return -2;
    }
    if(max_walltime()>0)
    {
      TimeStamp ts = TimeStamp::now();
      double tss = ts.seconds_since(opt_start_time_);
      if(tss > max_walltime())return -3;
    }
    return iflag;
  }

  assert(0);
  return iflag;
}

double CMinpackOptimizer::edm(unsigned n)
{
  // Use the running QR factorization of the Hessian to calculate
  // EDM = Df * H^-1 * Df . The part Q^-1 Df is already calculated by
  // minpack in "qtfvec_" so it suffices to solve for R * qtvec_
  // and then finish with a scalar multiplication by Df.
  for(unsigned iir = 0; iir<n; ++iir)
  {
    unsigned ir = n-iir-1;
    const double* rr = rvec_.data() + n*(n+1)/2 - iir*(iir+1)/2 - n;
    double sum = qtfvec_(ir);
    for(unsigned ic = ir+1; ic<n; ++ic)sum -= rr[ic]*dxvec_(ic);
        dxvec_(ir) = sum/rr[ir];
  }
  return dxvec_.transpose()*gvec_;
}
