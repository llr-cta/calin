/*

   calin/math/optimizer.hpp -- Stephen Fegan -- 2015-03-12

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

#pragma once

#include <algorithm>
#include <vector>
#include <Eigen/Core>

#include "util/timestamp.hpp"
#include "io/log.hpp"
#include "function.hpp"
#include "hessian.hpp"

namespace calin { namespace math { namespace optimizer {

enum class OptimizerVerbosityLevel { SILENT, SUMMARY_ONLY, ALL_FCN_EVALS_ONLY,
    SUMMARY_AND_PROGRESS, SUMMARY_AND_ALL_FCN_EVALS, ELEVATED, MAX };

enum class OptimizationStatus { TOLERANCE_REACHED, STOPPED_AT_MAXCALLS,
    STOPPED_AT_MAXTIME, LIMITED_BY_PRECISION, OPTIMIZER_FAILURE };

CALIN_TYPEALIAS(ErrorMatrixStatus, hessian::ErrorMatrixStatus);

class Optimizer
{
 public:
  constexpr static double inf = std::numeric_limits<double>::infinity();
  constexpr static double pos_inf = inf;
  constexpr static double neg_inf = -inf;

  CALIN_TYPEALIAS(VerbosityLevel, OptimizerVerbosityLevel);

  Optimizer(function::MultiAxisFunction* fcn, bool adopt_fcn = false):
      fcn_(fcn), adopt_fcn_(adopt_fcn) { /* nothing to see here */ }
  virtual ~Optimizer();

  virtual bool is_local_optimizer() = 0;
  virtual bool requires_gradient() = 0;
  virtual bool requires_hessian() = 0;
  virtual bool requires_box_constraints() = 0;

  virtual bool can_estimate_error() = 0;
  virtual bool can_use_gradient() = 0;
  virtual bool can_use_hessian() = 0;
  virtual bool can_impose_box_constraints() = 0;

  virtual OptimizationStatus minimize(VecRef xopt, double& fopt) = 0;

  unsigned num_iterations() const { return iterations_; }
  double best_evaluation(VecRef xbest) const { xbest=xbest_; return fbest_; }
  OptimizationStatus last_minimize_status() const { return opt_status_; }
  const std::string& last_minimize_message() const { return opt_message_; }

  virtual ErrorMatrixStatus error_matrix_estimate(MatRef error_matrix) = 0;
  virtual ErrorMatrixStatus
  calc_error_matrix_and_eigenvectors(MatRef error_matrix,
                                 VecRef eigenvalues, MatRef eigenvectors) = 0;

  ErrorMatrixStatus calc_error_matrix(MatRef error_matrix)
  {
    Eigen::VectorXd eigenvalues;
    Eigen::MatrixXd eigenvectors;
    return
        calc_error_matrix_and_eigenvectors(error_matrix, eigenvalues,
                                           eigenvectors);
  }

  void set_verbosity_level(VerbosityLevel verbose = VerbosityLevel::SILENT) {
    verbose_=verbose; }
  VerbosityLevel verbosity_level() const { return verbose_; }

  void set_abs_tolerance(double tol) { abs_tol_ = std::max(0.0,tol); }
  void set_rel_tolerance(double tol) { rel_tol_ = std::max(0.0,tol); }
  void set_max_iterations(unsigned max_num) { max_iterations_ = max_num; }
  void set_max_walltime(double max_wtime) { max_walltime_ = max_wtime; }

  double abs_tolerance() const { return abs_tol_; }
  double rel_tolerance() const { return (abs_tol_<=0.0 and rel_tol_<=0.0 and
                                         max_iterations_==0 and
                                         max_walltime_<=0.0)?0.001:rel_tol_; }
  unsigned max_iterations() const { return max_iterations_; }
  double max_walltime() const { return max_walltime_; }

  double step_size_scale_factor() const { return stepsize_scale_; }
  void set_step_size_scale_factor(double sss = 1.0) { stepsize_scale_ = sss; }

  void set_initial_values(const std::vector<double>& x0 = {}) { x0_=x0; }
  void set_initial_values(ConstVecRef x0) {
    x0_.assign(x0.data(),x0.data()+x0.size()); }
  void set_scale(const std::vector<double>& xscale = {}) { xscale_=xscale; }
  void set_limit_lo(unsigned ipar, double lim) {
    if(ipar>=xlim_lo_.size())xlim_lo_.resize(ipar+1, neg_inf);
    xlim_lo_[ipar]=lim; }
  void set_limit_hi(unsigned ipar, double lim) {
    if(ipar>=xlim_hi_.size())xlim_hi_.resize(ipar+1, pos_inf);
    xlim_hi_[ipar]=lim; }
  void set_limits_lo(const std::vector<double>& xlim = {}) { xlim_lo_ = xlim; }
  void set_limits_hi(const std::vector<double>& xlim = {}) { xlim_hi_ = xlim; }

  void clear_initial_values() { x0_.clear(); }
  void clear_scale() { xscale_.clear(); }
  void clear_limit_lo(unsigned ipar) {
    if(ipar<=xlim_lo_.size())xlim_lo_[ipar] = neg_inf; }
  void clear_limit_hi(unsigned ipar) {
    if(ipar<=xlim_hi_.size())xlim_hi_[ipar] = pos_inf; }
  void clear_limits_lo() { xlim_lo_.clear(); }
  void clear_limits_hi() { xlim_hi_.clear(); }

  std::vector<double> initial_values() const;
  std::vector<double> initial_stepsize() const;
  std::vector<double> limits_lo() const;
  std::vector<double> limits_hi() const;

  static Optimizer* create_optimizer_for_function(
                    function::MultiAxisFunction* fcn, bool adopt_fcn = false);
  static Optimizer* create_optimizer_by_name(const std::string& name,
                    function::MultiAxisFunction* fcn, bool adopt_fcn = false);
 protected:
  bool is_near_lo_limit(unsigned iaxis, double value) const;
  bool is_near_hi_limit(unsigned iaxis, double value) const;

  void opt_starting(const std::string& opt_name,
                    bool requires_gradient, bool requires_hessian,
                    const std::vector<double>& lower_limit,
                    const std::vector<double>& upper_limit,
                    const std::vector<double>& step_size);
  void opt_progress(double fval, const Eigen::VectorXd& x,
                    const Eigen::VectorXd* gradient = nullptr,
                    const Eigen::MatrixXd* hessian = nullptr);
  void opt_finished(OptimizationStatus status, double fopt,
                    const Eigen::VectorXd& xopt, const double* edm = nullptr);

  function::MultiAxisFunction* fcn_  { nullptr };
  bool adopt_fcn_ { false };
  VerbosityLevel verbose_ { VerbosityLevel::SILENT };
  std::vector<double> x0_ { };
  std::vector<double> xscale_ { };
  std::vector<double> xlim_lo_ { };
  std::vector<double> xlim_hi_ { };
  double stepsize_scale_ { 1.0 };
  double abs_tol_ { 0.0 };
  double rel_tol_ { 0.0 };
  unsigned max_iterations_ { 0 };
  double max_walltime_ { 0.0 };

  OptimizationStatus opt_status_ { OptimizationStatus::OPTIMIZER_FAILURE };
  std::string opt_message_ { "Live long and prosper" };
  calin::util::timestamp::Timestamp opt_start_time_;
  unsigned iterations_ { 0 };
  double fbest_ { inf };
  Eigen::VectorXd xbest_;
  unsigned wfval_ { 0 };
  unsigned pfval_ { 0 };
  unsigned progress_update_iter_ { 0 };
  double progress_update_time_ { 0.0 };
};

class ErrorMatrixEstimator
{
 public:
  ErrorMatrixEstimator(double error_up):
      npar_(0), error_up_(error_up) { /* nothing to see here */ }
  virtual ~ErrorMatrixEstimator();
  virtual void reset(unsigned npar) = 0;
  virtual void incorporate_func_value(ConstVecRef x, double f_val) = 0;
  virtual void incorporate_func_gradient(ConstVecRef x, double f_val,
                                         ConstVecRef gradient) = 0;
  virtual void incorporate_func_hessian(ConstVecRef x, double f_val,
                                        ConstVecRef gradient,
                                        ConstMatRef hessian) = 0;
  virtual ErrorMatrixStatus error_matrix(MatRef error_matrix) = 0;
 protected:
  unsigned npar_ { 0 };
  double error_up_ { 0.5 };
};

class IdentityErrorMatrixEstimator: public ErrorMatrixEstimator
{
 public:
  using ErrorMatrixEstimator::ErrorMatrixEstimator;
  ~IdentityErrorMatrixEstimator();
  void reset(unsigned npar) override;
  void incorporate_func_value(ConstVecRef x, double f_val) override;
  void incorporate_func_gradient(ConstVecRef x, double f_val,
                                 ConstVecRef gradient) override;
  void incorporate_func_hessian(ConstVecRef x, double f_val,
                                ConstVecRef gradient,
                                ConstMatRef hessian) override;
  ErrorMatrixStatus error_matrix(MatRef error_matrix) override;
};

class BFGSErrorMatrixEstimator: public ErrorMatrixEstimator
{
 public:
  using ErrorMatrixEstimator::ErrorMatrixEstimator;
  ~BFGSErrorMatrixEstimator();
  void reset(unsigned npar) override;
  void incorporate_func_value(ConstVecRef x, double f_val) override;
  void incorporate_func_gradient(ConstVecRef x, double f_val,
                                 ConstVecRef gradient) override;
  void incorporate_func_hessian(ConstVecRef x, double f_val,
                                ConstVecRef gradient,
                                ConstMatRef hessian) override;
  ErrorMatrixStatus error_matrix(MatRef error_matrix) override;
 protected:
  bool last_good_ { false };
  Eigen::MatrixXd Bk_;  // last estimate of error matrix
  Eigen::VectorXd xk_;  // position of last evaluation
  Eigen::VectorXd gk_;  // gradient at last evaluation
};

} } } // namespace calin::math::optimizer
