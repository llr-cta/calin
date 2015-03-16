/* 

   calin/math/nlopt_optimizer.hpp -- Stephen Fegan -- 2015-03-16

   Interface to MIT NLOpt optimizer suite

*/

#pragma once

#include "function.hpp"
#include "optimizer.hpp"
#include "nlopt/nlopt.hpp"

namespace calin { namespace math { namespace optimizer {

class NLOptOptimizer: public Optimizer
{
 public:
  using algorithm_type = nlopt::algorithm;

  NLOptOptimizer(algorithm_type algorithm,
                 MultiAxisFunction* fcn, bool adopt_fcn = false);
  ~NLOptOptimizer();

  static bool requires_gradient(algorithm_type a);
  static bool requires_hessian(algorithm_type a);
  static bool can_estimate_error(algorithm_type a);
  static bool can_use_gradient(algorithm_type a);
  static bool can_use_hessian(algorithm_type a);
  
  bool requires_gradient() override;
  bool requires_hessian() override;
  bool can_estimate_error() override;
  bool can_use_gradient() override;
  bool can_use_hessian() override;

  bool minimize(std::vector<double>& xopt) override;
  bool error_matrix_estimate(Eigen::MatrixXd& err_mat) override;
  bool calc_error_matrix(Eigen::MatrixXd& err_mat) override;

 protected:
  algorithm_type algorithm_;
};

} } } // namespace calin::math::optimizer
