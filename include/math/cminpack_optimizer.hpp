/* 

   calin/math/cminpack_optimizer.hpp -- Stephen Fegan -- 2015-05-21

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

#pragma once

#include <string>
#include <memory>
#include <Eigen/Dense>

#include "function.hpp"
#include "optimizer.hpp"

namespace calin { namespace math { namespace optimizer {

class CMinpackOptimizer: public Optimizer
{
 public:
  CMinpackOptimizer(function::MultiAxisFunction* fcn, bool adopt_fcn = false);
  ~CMinpackOptimizer();

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

  static int cminpack_callback_grad(void* self, int n, const double* x,
                                    double* grad, int iflag);

  static int cminpack_callback_hess(void *self, int n, const double* x,
                                    double* grad, double* hess, int ldfjac,
                                    int iflag);
  
  int eval_func(unsigned n, const double* x, double* grad, double* hess,
                int iflag);
  double edm(unsigned n);
  
  double fval_;
  Eigen::VectorXd xvec_;
  Eigen::VectorXd gvec_;
  Eigen::MatrixXd hmat_;
  Eigen::MatrixXd qmat_;
  Eigen::VectorXd rvec_;
  Eigen::VectorXd qtfvec_;
  Eigen::VectorXd dxvec_;
  unsigned nbest_ { 0 };
  std::unique_ptr<ErrorMatrixEstimator> err_est_;
};

} } } // namespace calin::math::optimizer
