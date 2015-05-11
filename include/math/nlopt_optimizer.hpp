/* 

   calin/math/nlopt_optimizer.hpp -- Stephen Fegan -- 2015-03-16

   Interface to MIT NLOpt optimizer suite

*/

#pragma once

#include <string>
#include <Eigen/Dense>

#include "function.hpp"
#include "optimizer.hpp"

namespace calin { namespace math { namespace optimizer {

class NLOptOptimizer: public Optimizer
{
 public:
  CALIN_TYPEALIAS(algorithm_type, std::string);

  NLOptOptimizer(const std::string& algorithm_name,
                 function::MultiAxisFunction* fcn, bool adopt_fcn = false);
  ~NLOptOptimizer();

  static bool is_valid_algorithm(const std::string& algorithm_name);
  static bool requires_gradient(const std::string& algorithm_name);
  static bool requires_hessian(const std::string& algorithm_name);
  static bool requires_box_constraints(const std::string& algorithm_name);
  static bool can_estimate_error(const std::string& algorithm_name);
  static bool can_use_gradient(const std::string& algorithm_name);
  static bool can_use_hessian(const std::string& algorithm_name);
  static bool can_impose_box_constraints(const std::string& algorithm_name);
  
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
  
  std::string algorithm_name_;

  std::unique_ptr<ErrorMatrixEstimator> err_est_;
};

} } } // namespace calin::math::optimizer
