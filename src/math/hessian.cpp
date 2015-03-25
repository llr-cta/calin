/* 

   calin/math/hessian.cpp -- Stephen Fegan -- 2015-02-24

*/

#include <iostream>
#include <cassert>
#include <Eigen/Eigenvalues>

#include "math/hessian.hpp"
#include "math/brent.hpp"

using namespace calin::math;

namespace {
inline double SQR(double x) { return x*x; }
}

hessian::ErrorStatus
hessian::hessian_to_error_matrix(MultiAxisFunction& fcn, ConstMatRef hessian,
                                 MatRef error_matrix, VecRef eigenvalues,
                                 MatRef eigenvectors)
{
  static constexpr auto eps = std::numeric_limits<double>::epsilon();
  const unsigned npar = fcn.num_domain_axes();
  const double scale = 2.0*fcn.error_up();
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
    eigenvalues(ipar) = scale/eigenvalues(ipar);
  }
  error_matrix = eigenvectors*eigenvalues.asDiagonal()*eigenvectors.inverse();
  return status;
}

hessian::ErrorStatus
hessian::hessian_to_error_matrix(MultiAxisFunction& fcn, ConstMatRef hessian,
                                 MatRef error_matrix)
{
  const unsigned npar = fcn.num_domain_axes();
  Eigen::VectorXd unused_eigenvalues(npar);
  Eigen::MatrixXd unused_eigenvectors(npar,npar);
  return hessian_to_error_matrix(fcn, hessian, error_matrix, unused_eigenvalues,
                                 unused_eigenvectors);
}

hessian::ErrorStatus
hessian::calculate_error_matrix(MultiAxisFunction& fcn, ConstVecRef x,
                                MatRef error_matrix)
{
  const unsigned npar = fcn.num_domain_axes();
  Eigen::MatrixXd hessian(npar,npar);
  calculate_hessian(fcn, x, hessian);
  return hessian_to_error_matrix(fcn, hessian, error_matrix);
}

Eigen::VectorXd hessian::
step_size_epsmult(ConstVecRef x, double eps_mult)
{
  constexpr double eps { std::numeric_limits<double>::epsilon() };
  unsigned npar = x.innerSize();
  Eigen::VectorXd dx(npar);
  for(unsigned ipar=0; ipar<npar; ipar++)
    dx(ipar) = (std::abs(x(ipar))*eps_mult)*eps;
  return dx;
}

Eigen::VectorXd hessian::
step_size_err_up(MultiAxisFunction& fcn, ConstVecRef x,
                 ConstVecRef error_hint, double err_up_frac, double tol)
{
  const double scale { 2.0*fcn.error_up() };
  double f0 { fcn.value(x) };
  double fup { f0+fcn.error_up()*err_up_frac };
  unsigned npar = x.innerSize();
  Eigen::VectorXd dx(npar);
  auto axes = fcn.domain_axes();
  if(error_hint.innerSize() == npar)
  {
    for(unsigned ipar=0;ipar<npar;ipar++)
      if(isfinite(error_hint(ipar)) and error_hint(ipar)>0)
        dx(ipar) = error_hint(ipar)*std::sqrt(scale*err_up_frac);
      else
        dx(ipar) = axes[ipar].scale;
  }
  else
  {
    for(unsigned ipar=0;ipar<npar;ipar++)
      dx(ipar) = axes[ipar].scale;
  }

  for(unsigned ipar=0;ipar<npar;ipar++)
  {
    Eigen::VectorXd xx { x };
    auto f_of_x = [ipar,&xx,&fcn,fup](double x){
      xx(ipar)=x; return fcn.value(xx)-fup; };
    double xlo = x(ipar);
    double flo = -fcn.error_up()*err_up_frac;
    double xhi = xlo + dx(ipar);
    double fhi = f_of_x(xhi);
    while(fhi<0){
      xlo=xhi; flo=fhi; dx(ipar)*=2.0; xhi=xlo+dx(ipar); fhi = f_of_x(xhi); };
    double xtol = std::abs((xhi-xlo)/(fhi-flo))*tol*fcn.error_up();
    double xroot = brent_zero(xlo,xhi,f_of_x,flo,fhi,xtol);
    dx(ipar) = xroot-x(ipar);
  }
  return dx;
}

void hessian::
calculate_hessian(MultiAxisFunction& fcn, ConstVecRef x, MatRef hessian)
{

}

void hessian::
calculate_hessian_gradient_dx(MultiAxisFunction& fcn, ConstVecRef x,
                              ConstVecRef dx, MatRef hessian)
{
  const unsigned npar = fcn.num_domain_axes();
  Eigen::VectorXd xx { x };
  Eigen::VectorXd grad_p(npar);
  Eigen::VectorXd grad_n(npar);
  hessian.resize(npar,npar);
  hessian.resize(npar,npar);
  for(unsigned ipar=0;ipar<npar;ipar++)
  {
    double h { dx(ipar) };
    double xp { x(ipar)+h };
    double xn { x(ipar)-h };
    double h2 = { xp-xn };
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
calculate_hessian_gradient_eps(MultiAxisFunction& fcn, ConstVecRef x,
                               MatRef hessian, double eps_mult)
{
  Eigen::VectorXd dx { step_size_epsmult(x, eps_mult) };
  calculate_hessian_gradient_dx(fcn, x, dx, hessian);
}

void hessian::
calculate_hessian_gradient_err_up(MultiAxisFunction& fcn, ConstVecRef x,
                                  MatRef hessian, ConstVecRef error_hint,
                                  double err_up_frac, double tol)
{
  
}
