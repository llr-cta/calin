/*

   calin/math/optimizer.cpp -- Stephen Fegan -- 2015-03-12

   Base class providing interface to function optimizers (minimizers)

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
#include <iomanip>
#include <algorithm>
#include <cmath>

#include <util/log.hpp>
#include <math/optimizer.hpp>
#include <math/nlopt_optimizer.hpp>

using namespace calin::math::optimizer;
using namespace calin::util::log;

constexpr double Optimizer::inf;
constexpr double Optimizer::pos_inf;
constexpr double Optimizer::neg_inf;

Optimizer::~Optimizer()
{
  if(adopt_fcn_)delete fcn_;
}

std::vector<double> Optimizer::initial_values() const
{
  // Set initial values from caller's values with fallback to function defaults
  auto axes = fcn_->domain_axes();
  std::vector<double> x0;
  if(x0_.size() == axes.size())x0 = x0_;
  else for(const auto& ipar : axes)x0.push_back(ipar.initial_value);
  return x0;
}

std::vector<double> Optimizer::initial_stepsize() const
{
  // Set the step size from caller's values with fallback to function defaults
  auto axes = fcn_->domain_axes();
  std::vector<double> xscale;
  if(xscale_.size() == axes.size())xscale = xscale_;
  else for(const auto& ipar : axes)xscale.push_back(ipar.scale);
  for(auto& x : xscale)x *= stepsize_scale_;
  return xscale;
}

std::vector<double> Optimizer::limits_lo() const
{
  // Set lower limits from caller's values with fallback to function defaults
  auto axes = fcn_->domain_axes();
  std::vector<double> xlim(axes.size());
  std::transform(axes.begin(), axes.end(), xlim.begin(),
                 [](const decltype(axes)::value_type& ipar)
                 { return ipar.has_lo_bound?ipar.lo_bound:neg_inf;});
  if(xlim_lo_.size() <= axes.size())
    std::transform(xlim_lo_.begin(), xlim_lo_.end(), xlim.begin(),
                   xlim.begin(), [](double ipar1, double ipar2) {
                     return std::max(ipar1, ipar2); });
  return xlim;
}

std::vector<double> Optimizer::limits_hi() const
{
  // Set upper limits from caller's values with fallback to function defaults
  auto axes = fcn_->domain_axes();
  std::vector<double> xlim(axes.size());
  std::transform(axes.begin(), axes.end(), xlim.begin(),
                 [](const decltype(axes)::value_type& ipar)
                 { return ipar.has_hi_bound?ipar.hi_bound:pos_inf;});
  if(xlim_hi_.size() <= axes.size())
    std::transform(xlim_hi_.begin(), xlim_hi_.end(), xlim.begin(),
                   xlim.begin(), [](double ipar1, double ipar2) {
                     return std::min(ipar1, ipar2); });
  return xlim;
}

Optimizer* Optimizer::create_optimizer_for_function(
    function::MultiAxisFunction* fcn, bool adopt_fcn)
{
  if(fcn->can_calculate_gradient())
    return new NLOptOptimizer("LD_LBFGS", fcn, adopt_fcn);
  else
    return new NLOptOptimizer("LN_NELDERMEAD", fcn, adopt_fcn);
  return nullptr;
}

Optimizer* Optimizer::create_optimizer_by_name(const std::string& name,
    function::MultiAxisFunction* fcn, bool adopt_fcn)
{
  if(NLOptOptimizer::is_valid_algorithm(name))
    return new NLOptOptimizer(name, fcn, adopt_fcn);
  return nullptr;
}

void Optimizer::opt_starting(const std::string& opt_name,
                             bool requires_gradient, bool requires_hessian,
                             const std::vector<double>& lower_limit,
                             const std::vector<double>& upper_limit,
                             const std::vector<double>& step_size)
{
  // Set up common (base-class) variables used to keep track of optimization
  // progress and status

  opt_status_           = OptimizationStatus::OPTIMIZER_FAILURE;
  opt_message_          = "Optimizer did not run";
  opt_start_time_       = calin::util::timestamp::Timestamp::now();
  iterations_           = 0;
  fbest_                = inf;
  xbest_                = std_to_eigenvec(initial_values());
  progress_update_iter_ = 0;
  progress_update_time_ = 0;

  unsigned naxis = fcn_->num_domain_axes();
  unsigned wnaxis = std::max(1U, num_digits(naxis));

  // Print status message if required by verbosity
  if(verbose_ != OptimizerVerbosityLevel::SILENT and
     verbose_ != OptimizerVerbosityLevel::ALL_FCN_EVALS_ONLY)
  {
    auto L = LOG(INFO);
    L << "Optimization using " << opt_name;
    if(requires_gradient || requires_hessian)
    {
      L << " (requires: ";
      if(requires_gradient)
      {
        L << "gradient";
        if(requires_hessian)L << " and hessian)";
        else L << ')';
      }
      else
        L << "hessian)";
    }
    L << '\n';

    L << "- stopping criteria:";
    bool has_crit = false;
    if(abs_tolerance() > 0) {
      L << " dF < " << abs_tolerance(); has_crit = true; }
    if(rel_tolerance() > 0) { if(has_crit) L << " or";
      L << " dF/F < " << rel_tolerance(); has_crit = true; }
    if(max_iterations() > 0) { if(has_crit) L << " or";
      L << " N_eval > " << max_iterations(); }
    if(max_walltime() > 0) { if(has_crit) L << " or";
      L << " T_wall > " << max_walltime() << " sec"; }
    L << '\n';

    L << "- function: N_dim = " << naxis;
    if(fcn_->can_calculate_gradient() || fcn_->can_calculate_hessian())
    {
      L << ", provides: ";
      if(fcn_->can_calculate_gradient())
      {
        L << "gradient";
        if(fcn_->can_calculate_hessian())L << " and hessian";
      }
      else
        L << "hessian";
    }
    L << '\n';
  }

  if(this->can_impose_box_constraints())
  {
    bool xinit_inside_xlim_lo = true;
    // The STL function "equal" is badly named!!
    if((unsigned)xbest_.size() == lower_limit.size())
      xinit_inside_xlim_lo =
          std::equal(xbest_.data(), xbest_.data()+xbest_.size(),
                     lower_limit.begin(), std::greater_equal<double>());

    bool xinit_inside_xlim_hi = true;
    // The STL function "equal" is badly named!!
    if((unsigned)xbest_.size() == upper_limit.size())
      xinit_inside_xlim_hi =
          std::equal(xbest_.data(), xbest_.data()+xbest_.size(),
                     upper_limit.begin(), std::less_equal<double>());

    if(!xinit_inside_xlim_lo or !xinit_inside_xlim_hi)
    {
      std::string message;
      message += "Initial values outside ";
      if(!xinit_inside_xlim_lo) {
        message += "lower";
        if(!xinit_inside_xlim_hi)message += " and upper";
      } else {
        message += "upper";
      }
      message += " limits";
      LOG(WARNING) << message;
    }
  }

  if(verbose_ == OptimizerVerbosityLevel::SUMMARY_ONLY or
     verbose_ == OptimizerVerbosityLevel::SUMMARY_AND_PROGRESS)return;

  auto L = LOG(INFO);

  unsigned wname = 0;
  auto axes = fcn_->domain_axes();
  for(const auto& a : axes)wname = std::max(wname, (unsigned)a.name.size());
  wname = std::min(wname, 40U);

  L << "List of function dimensions: \n"
    << "- " << std::left
    << std::setw(wnaxis) << "#" << ' '
    << std::setw(wname) << "Name" << ' '
    << std::setw(12) << "Initial val" << ' ';
  if(this->can_impose_box_constraints())
    L << std::setw(8) << "Lo bound" << ' '
      << std::setw(8) << "Hi bound" << ' ';
  L << std::setw(8) << "Stepsize" << '\n';

  for(unsigned iaxis = 0; iaxis<axes.size(); iaxis++)
  {
    L << "- " << std::setw(wnaxis) << iaxis << ' '
      << std::setw(wname) << axes[iaxis].name << ' '
      << std::setw(12) << xbest_[iaxis] << ' ';
    if(this->can_impose_box_constraints())
    {
      if(iaxis < lower_limit.size())
        L << std::setw(8) << lower_limit[iaxis] << ' ';
      else
        L << std::setw(8) << '-' << ' ';
      if(iaxis < upper_limit.size())
        L << std::setw(8) << upper_limit[iaxis] << ' ';
      else
        L << std::setw(8) << '-' << ' ';
    }
    if(iaxis < step_size.size())
      L << std::setw(8) << step_size[iaxis] << '\n';
    else
      L << std::setw(8) << '-' << '\n';
  }
}

void Optimizer::opt_progress(double fval, const Eigen::VectorXd& x,
                             const Eigen::VectorXd* gradient,
                             const Eigen::MatrixXd* hessian)
{
  double flast = fbest_;
  if(iterations_ == 0 or fval < fbest_)fbest_ = fval, xbest_ = x;
  iterations_++;

  if(iterations_ == 1)
  {
    pfval_ = 0;
    if(abs_tolerance() > 0.0 and abs_tolerance() < 1.0)
      pfval_ = 1+unsigned(std::ceil(-std::log10(abs_tolerance())));
    if(rel_tolerance() > 0.0 and std::abs(fval)*rel_tolerance() < 1.0)
      pfval_ = std::max(pfval_, 1+unsigned(std::ceil(-std::log10(std::abs(fval)*
                                                     rel_tolerance()))));
    if(pfval_ == 0 and 0.001*fcn_->error_up() < 1.0)
      pfval_ = std::ceil(-std::log10(.001*fcn_->error_up()));

    std::ostringstream s;
    s << std::fixed << std::setprecision(pfval_) << fval;
    wfval_ = s.str().size();
  }

  if(verbose_ == OptimizerVerbosityLevel::SILENT or
     (verbose_ == OptimizerVerbosityLevel::SUMMARY_ONLY))
    return;

  calin::util::timestamp::Timestamp ts = calin::util::timestamp::Timestamp::now();
  double tss = ts.seconds_since(opt_start_time_);
  if(verbose_ == OptimizerVerbosityLevel::SUMMARY_AND_PROGRESS or
     verbose_ == OptimizerVerbosityLevel::ELEVATED)
  {
    if(fbest_ < flast
        or iterations_-progress_update_iter_ > 100
        or tss-progress_update_time_ > 15.0)
      progress_update_iter_ = iterations_, progress_update_time_ = tss;
    else
      return;
  }

  auto L = LOG(VERBOSE);

  double edm = fbest_-flast;

#if 0
  if(err_mat_est != nullptr and gradient != nullptr)
  {
    edm = ((gradient->transpose())*(*err_mat_est)*(*gradient));
    edm /= 2.0*fcn_->error_up();
  }
#endif

  L << std::left << std::setw(std::max(3U, num_digits(max_iterations())))
    << iterations_ << ' ';
  if(hessian)L << "2 ";
  else if(gradient)L << "1 ";
  else L << "0 ";
  L << std::fixed
    << std::setw(wfval_) << std::setprecision(pfval_) << fval << ' '
    << std::setw(wfval_) << std::setprecision(pfval_) << fbest_ << ' '
    << std::right
    << std::setw(wfval_) << std::setprecision(pfval_) << edm << ' '
    << std::scientific << std::setw(10) << std::setprecision(3)
    << edm/fbest_ << ' '
    << std::fixed << std::setw(9) << std::setprecision(3) << tss << '\n';
}

#if 0
enum class OptimizerVerbosityLevel { SILENT, SUMMARY_ONLY, ALL_FCN_EVALS_ONLY,
    SUMMARY_AND_PROGRESS, SUMMARY_AND_FCN_EVALS, ELEVATED, MAX };
#endif

bool Optimizer::is_near_lo_limit(unsigned iaxis, double value) const
{
  if(iaxis >= xlim_lo_.size())return false;
  if(iaxis < xlim_hi_.size())
    return value <= xlim_lo_[iaxis] + (xlim_hi_[iaxis]-xlim_lo_[iaxis])*0.001;
  else
    return value <= xlim_lo_[iaxis] * 1.001 + std::numeric_limits<double>::epsilon() * 1000.0;
}

bool Optimizer::is_near_hi_limit(unsigned iaxis, double value) const
{
  if(iaxis >= xlim_hi_.size())return false;
  if(iaxis < xlim_lo_.size())
    return value >= xlim_hi_[iaxis] - (xlim_hi_[iaxis]-xlim_lo_[iaxis])*0.001;
  else
    return value >= xlim_hi_[iaxis] / 1.001 - std::numeric_limits<double>::epsilon() * 1000.0;
}

void Optimizer::opt_finished(OptimizationStatus status, double fopt,
                             const Eigen::VectorXd& xopt, const double* edm)
{
  if(verbose_ == OptimizerVerbosityLevel::SILENT or
     (verbose_ == OptimizerVerbosityLevel::ALL_FCN_EVALS_ONLY))return;

  Level level = Level::FAILURE;
  switch(status)
  {
    case OptimizationStatus::TOLERANCE_REACHED:
    case OptimizationStatus::STOPPED_AT_MAXCALLS:
    case OptimizationStatus::STOPPED_AT_MAXTIME:
      level = Level::SUCCESS;
      break;
    case OptimizationStatus::LIMITED_BY_PRECISION:
      level = Level::WARNING;
      break;
    case OptimizationStatus::OPTIMIZER_FAILURE:
      level = Level::FAILURE;
      break;
  }

  calin::util::timestamp::Timestamp ts = calin::util::timestamp::Timestamp::now();
  double tss = ts.seconds_since(opt_start_time_);

  LOG(level)
      << "Optimization finished with status: " << opt_message_ << '\n'
      << "- " << iterations_ << " function evaluations in "
      << std::fixed << std::setprecision(3) << tss << " seconds\n"
      << "- Best function value: "
      << std::setw(wfval_) << std::setprecision(pfval_) << fopt
      << disable_logging_if(edm == nullptr)
      << " (EDM: " << std::scientific << (edm==nullptr?0.0:*edm) << ")"
      << enable_logging() << '\n';

  if(verbose_ != OptimizerVerbosityLevel::ELEVATED and
     verbose_ != OptimizerVerbosityLevel::MAX)return;

  bool has_err_mat = this->can_estimate_error();
  Eigen::MatrixXd err_mat;
  if(has_err_mat)
    has_err_mat =
        this->error_matrix_estimate(err_mat) != ErrorMatrixStatus::UNAVAILABLE;

  auto axes = fcn_->domain_axes();
  unsigned wnaxis = std::max(1U, num_digits(axes.size()));
  unsigned wname = 0;
  for(const auto& a : axes)wname = std::max(wname, (unsigned)a.name.size());
  wname = std::min(wname, 40U);

  auto L = LOG(INFO);
  L << "Final position values";
  if(has_err_mat) L << " and APPROXIMATE errors";
  L << ":\n";
  for(unsigned iaxis = 0; iaxis<axes.size(); iaxis++)
  {
    L << "- " << std::left << std::setw(wnaxis) << iaxis << ' '
      << std::setw(wname) << axes[iaxis].name << ' '
      << xopt(iaxis);
    if(has_err_mat)
      L << " +/- " << std::sqrt(err_mat(iaxis,iaxis));
    if(is_near_lo_limit(iaxis, xopt(iaxis))) {
      L << " (near lower";
      if(is_near_hi_limit(iaxis, xopt(iaxis)))
        L << " and upper";
      L << " limit)";
    } else if(is_near_hi_limit(iaxis, xopt(iaxis))) {
      L << " (near upper limit)";
    }
    L << '\n';
  }

  if(verbose_ != OptimizerVerbosityLevel::MAX)return;
  if(has_err_mat)
    L << "\nAPPROXIMATE error matrix:\n" << err_mat << '\n';
}

// =============================================================================
//
// Error matrix estimator
//
// =============================================================================

ErrorMatrixEstimator::~ErrorMatrixEstimator()
{
  // nothing to see here
}

IdentityErrorMatrixEstimator::~IdentityErrorMatrixEstimator()
{
  // nothing to see here
}

void IdentityErrorMatrixEstimator::reset(unsigned npar)
{
  npar_ = npar;
}

void IdentityErrorMatrixEstimator::
incorporate_func_value(ConstVecRef x, double f_val)
{
  // nothing to see here
}

void IdentityErrorMatrixEstimator::
incorporate_func_gradient(ConstVecRef x, double f_val,
                          ConstVecRef gradient)
{
  // nothing to see here
}

void IdentityErrorMatrixEstimator::
incorporate_func_hessian(ConstVecRef x, double f_val,
                         ConstVecRef gradient, ConstMatRef hessian)
{
  // nothing to see here
}

ErrorMatrixStatus IdentityErrorMatrixEstimator::
error_matrix(MatRef error_matrix)
{
  error_matrix.resize(npar_,npar_);
  error_matrix.setIdentity();
  return ErrorMatrixStatus::UNAVAILABLE;
}

BFGSErrorMatrixEstimator::~BFGSErrorMatrixEstimator()
{
  // nothing to see here
}

void BFGSErrorMatrixEstimator::reset(unsigned npar)
{
  npar_ = npar;
  last_good_ = false;
  Bk_.resize(npar_,npar_);
  Bk_.setIdentity();
  xk_.resize(npar);
  gk_.resize(npar);
}

void BFGSErrorMatrixEstimator::
incorporate_func_value(ConstVecRef x, double f_val)
{
  // BGFS requires derivative so ignore last point for next update
  last_good_ = false;
}

void BFGSErrorMatrixEstimator::
incorporate_func_gradient(ConstVecRef x, double f_val,
                          ConstVecRef gradient)
{
  // Form error matrix estimate using BFGS update method, see
  // http://en.wikipedia.org/wiki/Broyden-Fletcher-Goldfarb-Shanno_algorithm

  if(!std::isfinite(f_val) or
     !std::all_of(gradient.data(), gradient.data()+npar_,
                  [](const double& x){return std::isfinite(x);}))
  {
    last_good_ = false;
    return;
  }

  if(last_good_)
  {
    Eigen::VectorXd sk = x;
    sk -= xk_;

    Eigen::VectorXd yk = gradient;
    yk -= gk_;

    const double skyk = sk.dot(yk);
    if(skyk == 0)goto skip_rank1_update; // don't like goto? i don't care :-)
    sk /= skyk;

#ifdef BFGS_COMPUTE_WITH_LOOPS // todo: test if faster when npar->inf
    Eigen::VectorXd bkyk(n);
    for(unsigned i=0;i<n;i++)
    {
      bkyk(i) = Bk_(i,i)*yk(i);
      for(unsigned j=i+1;j<n;j++)bkyk(i) += Bk_(i,j)*yk(j);
      for(unsigned j=0;j<i;j++)bkyk(i) += Bk_(j,i)*yk(j);
    }
    const double C1 = skyk + yk.dot(bkyk);
    for(unsigned i=0;i<n;i++)
      for(unsigned j=i;j<n;j++)
        Bk_(i,j) += C1*sk(i)*sk(j) - bkyk(i)*sk(j) - sk(i)*bkyk(j);
#else
    const Eigen::VectorXd bkyk = Bk_*yk;
    Bk_.noalias() += // What will Eigen make of this?
        (skyk + yk.dot(bkyk))*sk*sk.transpose()
        - bkyk*sk.transpose()
        - sk*bkyk.transpose();
#endif
    //std::cout << std::scientific << std::setprecision(8) << Bk_ << "\n\n";
  }
skip_rank1_update:

  last_good_ = true;
  xk_ = x;
  gk_ = gradient;
}

void BFGSErrorMatrixEstimator::
incorporate_func_hessian(ConstVecRef x, double f_val,
                         ConstVecRef gradient, ConstMatRef hessian)
{
  if(!std::isfinite(f_val) or
     !std::all_of(gradient.data(), gradient.data()+npar_,
                  [](const double& x){return std::isfinite(x);}) or
     !std::all_of(hessian.data(), hessian.data()+npar_*npar_,
                  [](const double& x){return std::isfinite(x);}))
  {
    last_good_ = false;
    return;
  }

  Bk_ = hessian.inverse();
  last_good_ = true;
  xk_ = x;
  gk_ = gradient;
}

ErrorMatrixStatus BFGSErrorMatrixEstimator::error_matrix(MatRef error_matrix)
{
  if(error_up_ != 0.5)
    error_matrix = (2.0*error_up_)*Bk_;
  else
    error_matrix = Bk_;
#ifdef BFGS_COMPUTE_WITH_LOOPS
  for(unsigned i=0;i<n;i++)
    for(unsigned j=i+1;j<n;j++)
      error_matrix(j,i) = error_matrix(i,j);
#endif
  return ErrorMatrixStatus::GOOD;
}
