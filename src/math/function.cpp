/* 

   calin/math/function.cpp -- Stephen Fegan -- 2015-02-24

*/

#include <iostream>

#include "math/function.hpp"

using namespace calin::math::function;

Parameterizable::~Parameterizable()
{
  // nothing to see here
}

MultiAxisFunction::~MultiAxisFunction()
{
  // nothing to see here
}

SingleAxisFunction::~SingleAxisFunction()
{
  // nothing to see here
}

std::vector<DomainAxis> SingleAxisFunction::domain_axes()
{
  return { domain_axis() };
}

double SingleAxisFunction::value(ConstVecRef x)
{
  return value(x(0));
}

double SingleAxisFunction::value_and_gradient(ConstVecRef x, VecRef gradient) 
{
  gradient.resize(1);
  return value_and_deriv(x(0),gradient(0));
}

double SingleAxisFunction::
value_gradient_and_hessian(ConstVecRef x, VecRef gradient, MatRef hessian) 
{
  gradient.resize(1);
  hessian.resize(1,1);
  return value_gradient_and_hessian(x(0),gradient(0),hessian(0,0));
}

bool calin::math::function::
gradient_check(MultiAxisFunction& fcn, ConstVecRef x, VecRef err,
               double eps_factor = 10.0)
{

}

bool calin::math::function::
gradient_check(MultiAxisFunction& fcn, ConstVecRef x, ConstVecRef dx,
               VecRef err)
{
  constexpr double eps = std::numeric_limits<double>::epsilon();
  unsigned fcn_num_axes = fcn.num_domain_axes();
  assert(fcn_num_axes == x.innerSize());
  assert(fcn_num_axes == dx.innerSize());
  double f0 = fcn.value(x);
  Eigen::VectorXd gradient(fcn_num_axes);
  double f0g = fcn.value_and_gradient(x, gradient);
  assert(fcn_num_axes == gradient.innerSize());
  if(std::abs(f0 - f0g)/std::abs(f0) > eps)return false;
  err.resize(fcn_num_axes);
  for(unsigned iaxis = 0;iaxis<fcn_num_axes;iaxis++)
  {
    const double h = dx(iaxis);
    Eigen::VectorXd xh = x;
    xh(iaxis) = x(iaxis) + h;
    const double fp1 = fcn.value(xh);
    xh(iaxis) = x(iaxis) + 2.0*h;
    const double fp2 = fcn.value(xh);
    xh(iaxis) = x(iaxis) - h;
    const double fm1 = fcn.value(xh);
    xh(iaxis) = x(iaxis) - 2.0*h;
    const double fm2 = fcn.value(xh);
    volatile double h2 = x(iaxis) + h;
    h2 -= (x(iaxis) - h);
    //double dfdx = (fm2 - 8.0*fm1 + 8.0*fp1 - fp2)/(6.0*h2);
    double dfdx = (-fm1 + fp1)/h2;
    double h2d3fdx3 = (-fm2 + 2.0*fm1 - 2.0*fp1 + fp2)/(6.0*h2);
    if(dfdx == gradient(iaxis))err(iaxis)=0;
    else err(iaxis) =
             std::abs(std::abs(dfdx - gradient(iaxis))/std::abs(h2d3fdx3)-1);
    std::cout << iaxis << ' ' << dfdx << ' ' << gradient(iaxis) << ' '
              << err(iaxis) << ' '
              << std::abs(dfdx - gradient(iaxis)) << ' '
              << h2d3fdx3 << ' '
              << std::abs(std::abs(dfdx - gradient(iaxis))-std::abs(h2d3fdx3))
              << '\n';
  }
  return true;
}

