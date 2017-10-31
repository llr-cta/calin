/* 

   calin/math/minuit75_optimizer.hpp -- Stephen Fegan -- 2015-03-12

   Object oriented interface to the 1975 FORTRAN version of MINUIT
   from CERN. See http://www.adsabs.harvard.edu/abs/1975CoPhC..10..343J

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

#include "function.hpp"
#include "optimizer.hpp"

namespace calin { namespace math { namespace optimizer {

class Minuit75Optimizer: public Optimizer
{
 public:
  Minuit75Optimizer(function::MultiAxisFunction* fcn, bool adopt_fcn = false);
  ~Minuit75Optimizer();

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
  int do_command(const std::string& command);
  static long minuit_callback(long* npar, double* grad, double* fcnval,
                              const double* x, long* iflag, void* futil);
  void eval_function(long npar, double* grad, double& fcnval,
                     const double* x, long& iflag);
  void* fcb_;
};

} } } // namespace calin::math::optimizer
