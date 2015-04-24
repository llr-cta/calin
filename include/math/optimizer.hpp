/* 

   calin/math/optimizer.hpp -- Stephen Fegan -- 2015-03-12

   Base class providing interface to function optimizers (minimizers)

*/

#pragma once

#include <algorithm>
#include <vector>
#include <Eigen/Core>

#include "function.hpp"
#include "hessian.hpp"

namespace calin { namespace math { namespace optimizer {

enum class OptimizerVerbosityLevel { SILENT, FCN_EVALS_ONLY, ELEVATED, MAX };

CALIN_TYPEALIAS(ErrorMatrixStatus, hessian::ErrorMatrixStatus);

class Optimizer
{
 public:
  constexpr static double inf = std::numeric_limits<double>::infinity();
  constexpr static double pos_inf = inf;
  constexpr static double neg_inf = -inf;

  CALIN_TYPEALIAS(VerbosityLevel, OptimizerVerbosityLevel);

  Optimizer(function::MultiAxisFunction* fcn, bool adopt_fcn = false):
      fcn_(fcn), my_fcn_(adopt_fcn) { /* nothing to see here */ }
  virtual ~Optimizer();

  virtual bool requires_gradient() = 0;
  virtual bool requires_hessian() = 0;
  virtual bool requires_box_constraints() = 0;

  virtual bool can_estimate_error() = 0;
  virtual bool can_use_gradient() = 0;
  virtual bool can_use_hessian() = 0;
  virtual bool can_impose_box_constraints() = 0;
  
  virtual bool minimize(VecRef xopt, double& fopt) = 0;
  virtual ErrorMatrixStatus error_matrix_estimate(MatRef err_mat) = 0;
  virtual ErrorMatrixStatus calc_error_matrix_and_eigenvectors(MatRef err_mat,
                                 VecRef eigenvalues, MatRef eigenvectors) = 0;

  ErrorMatrixStatus calc_error_matrix(MatRef err_mat)
  {
    Eigen::VectorXd eigenvalues;
    Eigen::MatrixXd eigenvectors;
    return
        calc_error_matrix_and_eigenvectors(err_mat, eigenvalues, eigenvectors);
  }
  
  void set_verbosity_level(VerbosityLevel verbose = VerbosityLevel::SILENT) {
    verbose_=verbose; }
  VerbosityLevel verbosity_level() const { return verbose_; }

  void set_abs_tolerance(double tol) { abs_tol_ = std::max(0.0,tol); }
  void set_rel_tolerance(double tol) { rel_tol_ = std::max(0.0,tol); }
  void set_max_iterations(unsigned max_num) { max_iterations_ = max_num; }

  double abs_tolerance() const { return abs_tol_; }
  double rel_tolerance() const { return (abs_tol_==0.0 and rel_tol_==0.0 and
                                         max_iterations_==0)?0.001:rel_tol_; }
  unsigned max_iterations() const { return max_iterations_; }
  
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
  
 protected:
  bool my_fcn_ { false };
  function::MultiAxisFunction* fcn_  { nullptr };
  VerbosityLevel verbose_ { VerbosityLevel::SILENT };
  std::vector<double> x0_ { };
  std::vector<double> xscale_ { };
  std::vector<double> xlim_lo_ { };
  std::vector<double> xlim_hi_ { };
  double stepsize_scale_ { 1.0 };
  double abs_tol_ { 0.0 };
  double rel_tol_ { 0.0 };
  unsigned max_iterations_ { 0 };
};

class ErrorMatrixEstimator
{
 public:
  ErrorMatrixEstimator(bool error_up):
      npar_(0), error_up_(error_up) { /* nothing to see here */ }
  virtual ~ErrorMatrixEstimator();
  virtual void reset(unsigned npar) = 0;
  virtual void incorporate_func_value(ConstVecRef x, double f_val) = 0;
  virtual void incorporate_func_gradient(ConstVecRef x, double f_val,
                                         ConstVecRef gradient) = 0;
  virtual ErrorMatrixStatus error_matrix(MatRef err_mat) = 0;
 protected:
  unsigned npar_;
  bool error_up_;
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
  ErrorMatrixStatus error_matrix(MatRef err_mat) override;
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
  ErrorMatrixStatus error_matrix(MatRef err_mat) override;
 protected:
  bool last_good_ { false };
  Eigen::MatrixXd Bk_;  // last estimate of error matrix
  Eigen::VectorXd xk_;  // position of last evaluation
  Eigen::VectorXd gk_;  // gradient at last evaluation
};

} } } // namespace calin::math::optimizer
