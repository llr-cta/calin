/*

   calin/math/hessian.cpp -- Stephen Fegan -- 2015-02-24

   Calculate the Hessian and/or error matrix of a function by finite
   differences

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, CNRS, Ecole Polytechnique, Institut Polytechnique de Paris

   This file is part of "calin"

   "calin" is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 or
   later, as published by the Free Software Foundation.

   "calin" is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

*/

#include <iostream>
#include <cassert>
#include <cmath>
#include <Eigen/Eigenvalues>

#include <math/hessian.hpp>
#include <math/brent.hpp>

using namespace calin::math;

#if 0 // UNUSED
namespace {
inline double SQR(double x) { return x*x; }
}
#endif

hessian::ErrorMatrixStatus
hessian::hessian_to_error_matrix(function::MultiAxisFunction& fcn,
                                 ConstMatRef hessian,
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
  hessian::ErrorMatrixStatus status = ErrorMatrixStatus::GOOD;
  for(unsigned ipar=0;ipar<npar;ipar++)
  {
    if(eigenvalues(ipar)<=0)
    {
      status = ErrorMatrixStatus::FORCED_POS_DEF;
      eigenvalues(ipar) = eps;
    }
    eigenvalues(ipar) = scale/eigenvalues(ipar);
  }
  error_matrix = eigenvectors*eigenvalues.asDiagonal()*eigenvectors.inverse();
  return status;
}

hessian::ErrorMatrixStatus
hessian::hessian_to_error_matrix(function::MultiAxisFunction& fcn,
                                 ConstMatRef hessian,
                                 MatRef error_matrix)
{
  const unsigned npar = fcn.num_domain_axes();
  Eigen::VectorXd unused_eigenvalues(npar);
  Eigen::MatrixXd unused_eigenvectors(npar,npar);
  return hessian_to_error_matrix(fcn, hessian, error_matrix,
                                 unused_eigenvalues, unused_eigenvectors);
}

hessian::ErrorMatrixStatus
hessian::calculate_error_matrix(function::MultiAxisFunction& fcn,
                                ConstVecRef x,
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
  // Remember: at minimum dF~dx^2 so use sqrt(epsilon) as basis for step size!
  const double eps { std::sqrt(std::numeric_limits<double>::epsilon()) };
  unsigned npar = x.size();
  Eigen::VectorXd dx(npar);
  for(unsigned ipar=0; ipar<npar; ipar++)
    dx(ipar) = (std::abs(x(ipar))*eps_mult)*eps;
  return dx;
}

Eigen::VectorXd hessian::
step_size_err_up(function::MultiAxisFunction& fcn, ConstVecRef x,
                 ConstVecRef error_hint, double err_up_frac, double tol)
{
  constexpr static double inf = std::numeric_limits<double>::infinity();
  const double scale { 2.0*fcn.error_up() };
  double f0 { fcn.value(x) };
  double fup { f0+fcn.error_up()*err_up_frac };
  unsigned npar = x.size();
  Eigen::VectorXd dx(npar);
  auto axes = fcn.domain_axes();
  if(error_hint.size() == npar)
  {
    for(unsigned ipar=0;ipar<npar;ipar++)
      if(std::isfinite(error_hint(ipar)) and error_hint(ipar)>0)
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
    auto f_of_x = [ipar,&xx,&fcn,fup](double x){ xx(ipar)=x;
#if 1
      std::cout << "step_size_err_up C: " << ipar << ' ' << x << ' ' << fcn.value(xx)-fup << '\n';
#endif
      return fcn.value(xx)-fup; };
    double xlo = x(ipar);
    double flo = -fcn.error_up()*err_up_frac;
    double xhi = xlo + dx(ipar);
    double xlim = inf;
    if(fcn.domain_axes()[ipar].has_hi_bound)
      xlim = fcn.domain_axes()[ipar].hi_bound;
    if(xhi > xlim)xhi=xlim;
    double fhi = f_of_x(xhi);
    int counter = 32;
    while(fhi<0 and xhi<xlim and counter--) {
      xlo=xhi;
      flo=fhi;
      dx(ipar)*=2.0;
#if 1
      std::cout << "step_size_err_up A: " << ipar << ' ' << dx.transpose() << '\n';
#endif
      xhi=std::min(xlo+dx(ipar), xlim);
      fhi = f_of_x(xhi);
    };
    if(fhi>=0)
    {
      double xtol = std::abs((xhi-xlo)/(fhi-flo))*tol*fcn.error_up();
      double xroot = brent_zero(xlo,xhi,f_of_x,flo,fhi,xtol);
#if 1
      std::cout << "step_size_err_up B: " << xlo << ' ' << xhi << ' ' << flo << ' ' << fhi
        << ' ' << xtol << ' ' << xroot << ' ' << x(ipar) << '\n';
#endif
      if(xroot == xlo)dx(ipar) = xhi-x(ipar);
      else dx(ipar) = xroot-x(ipar);
    }
    else
    {
      //dx(ipar) = xhi-x(ipar);
    }
  }
  //std::cout << dx << '\n';
  return dx;
}

void hessian::
calculate_hessian(function::MultiAxisFunction& fcn,
                  ConstVecRef x, MatRef hessian,
                  ConstVecRef error_hint)
{
  if(fcn.can_calculate_hessian())
  {
    Eigen::VectorXd gradient(fcn.num_domain_axes());
    fcn.value_gradient_and_hessian(x, gradient, hessian);
  }
  else if(fcn.can_calculate_gradient())
  {
    calculate_hessian_1st_order_err_up(fcn, x, hessian, error_hint);
  }
  else
  {
    calculate_hessian_2nd_order_err_up(fcn, x, hessian, error_hint);
  }
}

void hessian::
calculate_hessian_1st_order_dx(function::MultiAxisFunction& fcn,
                               ConstVecRef x,
                               ConstVecRef dx, MatRef hessian)
{
  const unsigned npar = fcn.num_domain_axes();
  Eigen::VectorXd xx { x };
  Eigen::VectorXd grad_p(npar);
  Eigen::VectorXd grad_n(npar);
  hessian.resize(npar,npar);
  for(unsigned ipar=0;ipar<npar;ipar++)
  {
    double h { dx(ipar) };
    if(fcn.domain_axes()[ipar].has_hi_bound)
      h = std::min(h, fcn.domain_axes()[ipar].hi_bound - x(ipar));
    if(fcn.domain_axes()[ipar].has_lo_bound)
      h = std::min(h, x(ipar) - fcn.domain_axes()[ipar].lo_bound);
    if(h<=0)throw std::runtime_error("calculate_hessian_1st_order_dx: "
      "axis is at limit: " + fcn.domain_axes()[ipar].name);
    double xp { x(ipar)+h };
    double xn { x(ipar)-h };
    double h2 = { xp-xn };
    xx(ipar) = xp;
    fcn.value_and_gradient(xx, grad_p);
    xx(ipar) = xn;
    fcn.value_and_gradient(xx, grad_n);
#if 0
    std::cout << ipar << ' ' << xp << ' ' << xn << ' '
              << "[ " << grad_p.transpose() << " ] [ "
              << grad_n.transpose() << "]\n";
#endif
    xx(ipar) = x(ipar); // reset for next loop
    for(unsigned jpar=0;jpar<npar;jpar++)
      hessian(ipar,jpar) = (grad_p(jpar) - grad_n(jpar))/h2;
  }

  // Symmetrize by averaging off-diagonal elements
  for(unsigned ipar=0;ipar<npar;ipar++)
    for(unsigned jpar=0;jpar<ipar;jpar++)
      hessian(ipar,jpar) = hessian(jpar,ipar) =
                           0.5*(hessian(ipar,jpar) + hessian(jpar,ipar));

#if 0
  std::cout << hessian << '\n';
#endif
}

void hessian::
calculate_hessian_1st_order_eps(function::MultiAxisFunction& fcn,
                                ConstVecRef x,
                                MatRef hessian, double eps_mult)
{
  Eigen::VectorXd dx { step_size_epsmult(x, eps_mult) };
  calculate_hessian_1st_order_dx(fcn, x, dx, hessian);
}

void hessian::
calculate_hessian_1st_order_err_up(function::MultiAxisFunction& fcn,
                                   ConstVecRef x,
                                   MatRef hessian, ConstVecRef error_hint,
                                   double err_up_frac, double tol)
{
  Eigen::VectorXd dx {
    step_size_err_up(fcn, x, error_hint, err_up_frac, tol) };
  calculate_hessian_1st_order_dx(fcn, x, dx, hessian);
}

void hessian::
calculate_hessian_2nd_order_dx(function::MultiAxisFunction& fcn,
                               ConstVecRef x,
                               ConstVecRef dx, MatRef hessian)
{
  const unsigned npar = fcn.num_domain_axes();
  Eigen::VectorXd xx { x };
  hessian.resize(npar,npar);
  double v0 = fcn.value(xx);
  for(unsigned ipar=0;ipar<npar;ipar++)
  {
    double h_i { dx(ipar) };
    if(fcn.domain_axes()[ipar].has_hi_bound)
      h_i = std::min(h_i, fcn.domain_axes()[ipar].hi_bound - x(ipar));
    if(fcn.domain_axes()[ipar].has_lo_bound)
      h_i = std::min(h_i, x(ipar) - fcn.domain_axes()[ipar].lo_bound);
    if(h_i<=0)throw std::runtime_error("calculate_hessian_2nd_order_dx: "
      "axis is at limit: " + fcn.domain_axes()[ipar].name);
    double xp_i { x(ipar)+h_i };
    double xn_i { x(ipar)-h_i };
    double h2_i = { xp_i-xn_i };
    xx(ipar) = xp_i;
    double vp_i = fcn.value(xx);
    xx(ipar) = xn_i;
    double vn_i = fcn.value(xx);

    // Abromowitz & Stegun 25.3.23
    hessian(ipar,ipar) = 4.0*(vn_i - 2.0*v0 + vp_i)/(h2_i*h2_i);

    for(unsigned jpar=0;jpar<ipar;jpar++)
    {
      double h_j { dx(jpar) };
      if(fcn.domain_axes()[jpar].has_hi_bound)
        h_j = std::min(h_j, fcn.domain_axes()[jpar].hi_bound - x(jpar));
      if(fcn.domain_axes()[ipar].has_lo_bound)
        h_j = std::min(h_j, x(jpar) - fcn.domain_axes()[jpar].lo_bound);
      if(h_j<=0)throw std::runtime_error("calculate_hessian_2nd_order_dx: "
        "axis is at limit: " + fcn.domain_axes()[jpar].name);
      double xp_j { x(jpar)+h_j };
      double xn_j { x(jpar)-h_j };
      double h2_j = { xp_j-xn_j };
      xx(ipar) = xp_i; xx(jpar) = xp_j; double vpp = fcn.value(xx);
      xx(ipar) = xp_i; xx(jpar) = xn_j; double vpn = fcn.value(xx);
      xx(ipar) = xn_i; xx(jpar) = xp_j; double vnp = fcn.value(xx);
      xx(ipar) = xn_i; xx(jpar) = xn_j; double vnn = fcn.value(xx);

      // Abromowitz & Stegun 25.3.26
      hessian(ipar,jpar) = hessian(jpar,ipar)
                         = (vpp - vpn - vnp + vnn)/(h2_i*h2_j);

      xx(jpar) = x(jpar); // reset for next loop
    }

    xx(ipar) = x(ipar); // reset for next loop
  }
}

void hessian::
calculate_hessian_2nd_order_eps(function::MultiAxisFunction& fcn,
                                ConstVecRef x,
                                MatRef hessian, double eps_mult)
{
  Eigen::VectorXd dx { step_size_epsmult(x, eps_mult) };
  calculate_hessian_2nd_order_dx(fcn, x, dx, hessian);
}

void hessian::
calculate_hessian_2nd_order_err_up(function::MultiAxisFunction& fcn,
                                   ConstVecRef x,
                                   MatRef hessian, ConstVecRef error_hint,
                                   double err_up_frac, double tol)
{
  Eigen::VectorXd dx {
    step_size_err_up(fcn, x, error_hint, err_up_frac, tol) };
  calculate_hessian_2nd_order_dx(fcn, x, dx, hessian);
}
