/* 

   calin/math/hessian.hpp -- Stephen Fegan -- 2015-03-23

   Calculate the Hessian and/or error matrix of a function by finite
   differences

*/

#pragma once

#include <string>
#include <vector>
#include <limits>

#include <Eigen/Core>

#include "function.hpp"

namespace calin { namespace math { namespace hessian {

using function::MatRef;
using function::VecRef;
using function::ConstVecRef;
using function::ConstMatRef;

enum class ErrorStatus { FORCED_POS_DEF, GOOD };

// Calculate error matrix by inverting the given Hessian, testing for
// positive-definiteness, and "correcting" negative eigenvalues if
// necessary. This version return the (corrected) eigenvalues and the
// eigenvectors.
ErrorStatus hessian_to_error_matrix(MultiAxisFunction& fcn, ConstMatRef hessian,
                                    MatRef error_matrix, VecRef eigenvalues,
                                    MatRef eigenvectors);

// Calculate error matrix by inverting the given Hessian, testing for
// positive-definiteness, and "correcting" negative eigenvalues if
// necessary.
ErrorStatus hessian_to_error_matrix(MultiAxisFunction& fcn, ConstMatRef hessian,
                                    MatRef error_matrix);

// Calculate error matrix by calling calculating the Hessiana and
// then inverting it.
ErrorStatus calculate_error_matrix(MultiAxisFunction& fcn, ConstVecRef x,
                                   MatRef error_matrix);

// Calculate step size as multiple of machine epslion
Eigen::VectorXd step_size_epsmult(ConstVecRef x, double eps_mult=100.0);

// Calculate step size to increase function value by fraction of its
// "error up", within a given tolerance. The vector "error hint"
// allows the user to pass in a hint on what size the errors are, if
// they are known
Eigen::VectorXd step_size_err_up(MultiAxisFunction& fcn, ConstVecRef x,
                                 ConstVecRef error_hint = Eigen::VectorXd(),
                                 double err_up_frac = 0.01, double tol = 0.005);

// Convenience function which chooses method depending on what Function
// supports
void calculate_hessian(MultiAxisFunction& fcn, ConstVecRef x, MatRef hessian);

// Calculate Hessian by numerical differentiation of gradient using
// the two-point difference formula - O(2N). Step sizes in each of the
// directions should be given.
void calculate_hessian_gradient_dx(MultiAxisFunction& fcn, ConstVecRef x,
                                    ConstVecRef dx, MatRef hessian);

// Calculate Hessian by numerical differentiation of gradient using
// the two-point difference formula - O(2N). The step sizes is given
// by multiple of machine epsilon.
void calculate_hessian_gradient_eps(MultiAxisFunction& fcn, ConstVecRef x,
                                    MatRef hessian, double eps_mult=100.0);

// Calculate Hessian by numerical differentiation of gradient using
// the two-point difference formula - O(2N). The step size is chosen
// to increase function by some fraction of "error up" from its
// minimum position. The tol parameter sets the tolerance how much the
// function increases. 
void calculate_hessian_gradient_err_up(MultiAxisFunction& fcn, ConstVecRef x,
                                       MatRef hessian,
                                     ConstVecRef error_hint = Eigen::VectorXd(),
                                       double err_up_frac = 0.01,
                                       double tol = 0.01);

} // namespave hessian

} } // namespace calin::math
