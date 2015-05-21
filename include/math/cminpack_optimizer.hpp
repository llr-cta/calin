/* 

   calin/math/cminpack_optimizer.hpp -- Stephen Fegan -- 2015-05-21

   Interface to ANL CMinpack optimizer suite

*/

#pragma once

#include <string>
#include <Eigen/Dense>

#include "function.hpp"
#include "optimizer.hpp"

namespace calin { namespace math { namespace optimizer {

class enum CMinpackMode { AUTO, REQUIRE_HESSIAN, DONT_USE_HESSIAN };

class CMinpackOptimizer: public Optimizer
{
 public:
  CMinpackOptimizer(CMinpackMode mode,
                    function::MultiAxisFunction* fcn, bool adopt_fcn = false);
  ~CMinpackOptimizer();

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

  int cminpack_callback_grad(void* self, int n, double* x, double* grad,
                             int iflag);

  int cminpack_callback_hess(void *self, int n, double* x, double* grad,
                             double* hess, int ldfjac, int iflag);
  
  int eval_func(unsigned n, const double* x, bool hessian);
  
  CMinpackMode mode_;
  Eigen::VectorXd gvec_;
  Eigen::MatrixXd hmat_;
  std::unique_ptr<ErrorMatrixEstimator> err_est_;
};

} } } // namespace calin::math::optimizer
