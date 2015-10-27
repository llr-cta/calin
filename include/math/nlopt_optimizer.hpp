/* 

   calin/math/nlopt_optimizer.hpp -- Stephen Fegan -- 2015-03-16

   Interface to MIT NLOpt optimizer suite

*/

#pragma once

#include <string>
#include <memory>
#include <Eigen/Dense>

#include "nlopt/nlopt_algorithm.h"

#include "function.hpp"
#include "optimizer.hpp"

namespace calin { namespace math { namespace optimizer {

class NLOptOptimizer: public Optimizer
{
 public:
  CALIN_TYPEALIAS(algorithm_type, nlopt_algorithm);
  
  NLOptOptimizer(algorithm_type algorithm,
                 function::MultiAxisFunction* fcn, bool adopt_fcn = false);
  NLOptOptimizer(const std::string& algorithm_name,
                 function::MultiAxisFunction* fcn, bool adopt_fcn = false);
  ~NLOptOptimizer();

  static algorithm_type algorithm_by_name(const std::string& algorithm_name);
  static bool is_valid_algorithm(const std::string& algorithm_name);

  static bool is_local_optimizer(algorithm_type algorithm);
  static bool requires_gradient(algorithm_type algorithm);
  static bool requires_hessian(algorithm_type algorithm);
  static bool requires_box_constraints(algorithm_type algorithm);
  static bool can_estimate_error(algorithm_type algorithm);
  static bool can_use_gradient(algorithm_type algorithm);
  static bool can_use_hessian(algorithm_type algorithm);
  static bool can_impose_box_constraints(algorithm_type algorithm);
  
  bool is_local_optimizer() override;
  bool requires_gradient() override;
  bool requires_hessian() override;
  bool requires_box_constraints() override;
  bool can_estimate_error() override;
  bool can_use_gradient() override;
  bool can_use_hessian() override;
  bool can_impose_box_constraints() override;

  OptimizationStatus minimize(VecRef xopt, double& fopt) override;
  
  ErrorMatrixStatus error_matrix_estimate(MatRef error_matrix) override;
  ErrorMatrixStatus calc_error_matrix_and_eigenvectors(MatRef error_matrix,
                     VecRef eigenvalues, MatRef eigenvectors) override;

 protected:
  static double nlopt_callback(unsigned n, const double* x, double* grad,
                               void* self);
  double eval_func(unsigned n, const double* x, double* grad);
  
  algorithm_type algorithm_;
  Eigen::VectorXd gvec_;
  std::unique_ptr<ErrorMatrixEstimator> err_est_;
};

} } } // namespace calin::math::optimizer
