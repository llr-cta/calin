/* 

   calin/math/optimizer.hpp -- Stephen Fegan -- 2015-03-12

   Base class providing interface to function optimizers (minimizers)

*/

#pragma once

#include <algorithm>
#include <vector>
#include <Eigen/Core>

#include "function.hpp"

namespace calin { namespace math { namespace optimizer {

enum class OptimizerVerbosityLevel { SILENT, FCN_EVALS_ONLY, ELEVATED, MAX };

class Optimizer
{
 public:
  using VerbosityLevel = OptimizerVerbosityLevel;
  Optimizer(MultiAxisFunction* fcn, bool adopt_fcn = false):
      fcn_(fcn), my_fcn_(adopt_fcn) { /* nothing to see here */ }
  virtual ~Optimizer();

  virtual bool requires_gradient() = 0;
  virtual bool requires_hessian() = 0;
  virtual bool requires_box_constraints() = 0;

  virtual bool can_estimate_error() = 0;
  virtual bool can_use_gradient() = 0;
  virtual bool can_use_hessian() = 0;
  virtual bool can_impose_box_constraints() = 0;
  
  virtual bool minimize(std::vector<double>& xopt, double& fopt) = 0;
  virtual bool error_matrix_estimate(Eigen::MatrixXd& err_mat) = 0;
  virtual bool calc_error_matrix(Eigen::MatrixXd& err_mat) = 0;

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
  void set_scale(const std::vector<double>& xscale = {}) { xscale_=xscale; }
  void set_limits_lo(const std::vector<double>& xlim = {}) { xlim_lo_ = xlim; }
  void set_limits_hi(const std::vector<double>& xlim = {}) { xlim_hi_ = xlim; }

  void clear_initial_values() { x0_.clear(); }
  void clear_scale() { xscale_.clear(); }
  void clear_limits_lo() { xlim_lo_.clear(); }
  void clear_limits_hi() { xlim_hi_.clear(); }

  std::vector<double> initial_values() const;
  std::vector<double> initial_stepsize() const;
  std::vector<double> limits_lo() const;
  std::vector<double> limits_hi() const;
  
 protected:
  bool my_fcn_ { false };
  MultiAxisFunction* fcn_  { nullptr };
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

} // namespace optimizer

using Optimizer = optimizer::Optimizer;

} } // namespace calin::math
