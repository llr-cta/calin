/* 

   calin/math/hessian.hpp -- Stephen Fegan -- 2015-03-23

   Calculate the Hessian and/or error matrix of a function by finite
   differences

   Copyright 2015, Stephen Fegan <sfegan@gmail.com>

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
#include <vector>
#include <limits>

#include <Eigen/Core>

#include "function.hpp"

namespace calin { namespace math { namespace hessian {

enum class ErrorMatrixStatus { UNAVAILABLE, FORCED_POS_DEF, GOOD };

// Calculate error matrix by inverting the given Hessian, testing for
// positive-definiteness, and "correcting" negative eigenvalues if
// necessary. This version return the (corrected) eigenvalues and the
// eigenvectors.
ErrorMatrixStatus hessian_to_error_matrix(function::MultiAxisFunction& fcn,
                                    ConstMatRef hessian,
                                    MatRef error_matrix, VecRef eigenvalues,
                                    MatRef eigenvectors);

// Calculate error matrix by inverting the given Hessian, testing for
// positive-definiteness, and "correcting" negative eigenvalues if
// necessary.
ErrorMatrixStatus hessian_to_error_matrix(function::MultiAxisFunction& fcn,
                                          ConstMatRef hessian,
                                          MatRef error_matrix);

// Calculate error matrix by calling calculating the Hessiana and
// then inverting it.
ErrorMatrixStatus calculate_error_matrix(function::MultiAxisFunction& fcn,
                                         ConstVecRef x, MatRef error_matrix);

// Calculate step size as multiple of machine epslion
Eigen::VectorXd step_size_epsmult(ConstVecRef x, double eps_mult=100.0);

// Calculate step size to increase function value by fraction of its
// "error up", within a given tolerance. The vector "error hint"
// allows the user to pass in a hint on what size the errors are, if
// they are known
Eigen::VectorXd step_size_err_up(function::MultiAxisFunction& fcn,
                                 ConstVecRef x,
                                 ConstVecRef error_hint = Eigen::VectorXd(),
                                 double err_up_frac = 0.001, double tol = 0.03);

// Convenience function which chooses method depending on what Function
// supports
void calculate_hessian(function::MultiAxisFunction& fcn,
                       ConstVecRef x, MatRef hessian,
                       ConstVecRef error_hint = Eigen::VectorXd());

// Calculate Hessian by numerical differentiation of gradient using
// the two-point difference formula - O(2N). Step sizes in each of the
// directions should be given.
void calculate_hessian_1st_order_dx(function::MultiAxisFunction& fcn,
                                    ConstVecRef x,
                                    ConstVecRef dx, MatRef hessian);

// Calculate Hessian by numerical differentiation of gradient using
// the two-point difference formula - O(2N). The step sizes is given
// by multiple of machine epsilon.
void calculate_hessian_1st_order_eps(function::MultiAxisFunction& fcn,
                                     ConstVecRef x,
                                     MatRef hessian, double eps_mult=100.0);

// Calculate Hessian by numerical differentiation of gradient using
// the two-point difference formula - O(2N). The step size is chosen
// to increase function by some fraction of "error up" from its
// minimum position. The tol parameter sets the tolerance how much the
// function increases. 
void calculate_hessian_1st_order_err_up(function::MultiAxisFunction& fcn,
                                     ConstVecRef x, MatRef hessian,
                                     ConstVecRef error_hint = Eigen::VectorXd(),
                                     double err_up_frac = 0.001,
                                     double tol = 0.001);

// Calculate Hessian by numerical differentiation of the function
// using value using the 3/4-point difference formula - O(4N*N). Step
// sizes in each of the directions should be given.
void calculate_hessian_2nd_order_dx(function::MultiAxisFunction& fcn,
                                    ConstVecRef x,
                                    ConstVecRef dx, MatRef hessian);

// Calculate Hessian by numerical differentiation of the function
// using value using the 3/4-point difference formula - O(4N*N). The
// step sizes is given by multiple of machine epsilon.
void calculate_hessian_2nd_order_eps(function::MultiAxisFunction& fcn,
                                     ConstVecRef x,
                                     MatRef hessian, double eps_mult=100.0);

// Calculate Hessian by numerical differentiationof the function using
// value using the 3/4-point difference formula - O(4N*N). The step
// size is chosen to increase function by some fraction of "error up"
// from its minimum position. The tol parameter sets the tolerance how
// much the function increases.
void calculate_hessian_2nd_order_err_up(function::MultiAxisFunction& fcn,
                                     ConstVecRef x, MatRef hessian,
                                     ConstVecRef error_hint = Eigen::VectorXd(),
                                     double err_up_frac = 0.001,
                                     double tol = 0.001);

} } } // namespace calin::math::hessian
