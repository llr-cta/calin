/* 

   calin/math/hessian.cpp -- Stephen Fegan -- 2015-02-24

*/

#include <iostream>
#include <cassert>
#include <Eigen/Eigenvalues>

#include "math/hessian.hpp"

using namespace calin::math;

hessian::ErrorStatus
hessian::calculate_error_matrix(MultiAxisFunction& fcn, ConstMatRef hessian,
                                MatRef error_matrix, VecRef eigenvalues,
                                MatRef eigenvectors)
{
  static constexpr auto eps = std::numeric_limits<double>::epsilon();
  const unsigned npar = fcn.num_domain_axes();
  error_matrix.resize(npar,npar);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>
      es(hessian,Eigen::ComputeEigenvectors);
  eigenvalues = es.eigenvalues();
  eigenvectors = es.eigenvectors();
  hessian::ErrorStatus status = ErrorStatus::GOOD;
  for(unsigned ipar=0;ipar<npar;ipar++)
  {
    if(eigenvalues(ipar)<=0)
    {
      status = ErrorStatus::FORCED_POS_DEF;
      eigenvalues(ipar) = eps;
    }
    eigenvalues(ipar) = 1.0/eigenvalues(ipar);
  }
  error_matrix = eigenvectors*eigenvalues.asDiagonal()*eigenvectors.inverse();
  return status;
}

hessian::ErrorStatus
hessian::calculate_error_matrix(MultiAxisFunction& fcn, ConstMatRef hessian,
                                MatRef error_matrix)
{
  const unsigned npar = fcn.num_domain_axes();
  Eigen::VectorXd unused_eigenvalues(npar);
  Eigen::MatrixXd unused_eigenvectors(npar,npar);
  return calculate_error_matrix(fcn, hessian, error_matrix, unused_eigenvalues,
                                unused_eigenvectors);
}

hessian::ErrorStatus
hessian::calculate_error_matrix(MultiAxisFunction& fcn, ConstVecRef x)
{

}

void hessian::
calculate_hessian(MultiAxisFunction& fcn, ConstVecRef x, MatRef hessian)
{

}

void hessian::
calculate_hessian_gradient_eps(MultiAxisFunction& fcn, ConstVecRef x,
                               MatRef hessian, double eps_mult)
{
  static constexpr auto eps = std::numeric_limits<double>::epsilon();
  const unsigned npar = fcn.num_domain_axes();
  Eigen::VectorXd xx { x };
  Eigen::VectorXd grad_p(npar);
  Eigen::VectorXd grad_n(npar);
  hessian.resize(npar,npar);
  for(unsigned ipar=0;ipar<npar;ipar++)
  {
    double h = (x(ipar)*eps_mult)*eps;
    double xp = x(ipar)+h;
    double xn = x(ipar)-h;
    double h2 = xp-xn;
    xx(ipar) = xp;
    fcn.value_and_gradient(xx, grad_p);
    xx(ipar) = xn;
    fcn.value_and_gradient(xx, grad_n);
    xx(ipar) = x(ipar); // reset for next loop
    for(unsigned jpar=0;jpar<npar;jpar++)
      hessian(ipar,jpar) = (grad_p(jpar) - grad_n(jpar))/h2;
  }

  // Symmetrize by averaging off-diagonal elements
  for(unsigned ipar=0;ipar<npar;ipar++)
    for(unsigned jpar=0;jpar<ipar;jpar++)
      hessian(ipar,jpar) = hessian(jpar,ipar) =
                           0.5*(hessian(ipar,jpar) + hessian(jpar,ipar));
}

void hessian::
calculate_hessian_gradient_err_up(MultiAxisFunction& fcn, ConstVecRef x,
                                  MatRef hessian, ConstVecRef error_hint,
                                  double err_up_frac, double tol)
{
  
}
