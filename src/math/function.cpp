/* 

   calin/math/function.cpp -- Stephen Fegan -- 2015-02-24

*/

#include <iostream>
#include <iomanip>
#include <stdexcept>

#include "math/function.hpp"

using namespace calin::math::function;

namespace {

inline static double SQR(double x) { return x*x; }
constexpr double c_gauss_norm = 0.5*M_2_SQRTPI*M_SQRT1_2;

} // anonymous namespace

Parameterizable::~Parameterizable()
{
  // nothing to see here
}

MultiAxisFunction::~MultiAxisFunction()
{
  // nothing to see here
}

// *****************************************************************************
//
// SingleAxisFunction
//
// *****************************************************************************

SingleAxisFunction::~SingleAxisFunction()
{
  // nothing to see here
}

unsigned SingleAxisFunction::num_domain_axes()
{
  return 1;
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
  return value_and_gradient(x(0),gradient(0));
}

double SingleAxisFunction::
value_gradient_and_hessian(ConstVecRef x, VecRef gradient, MatRef hessian) 
{
  gradient.resize(1);
  hessian.resize(1,1);
  return value_gradient_and_hessian(x(0),gradient(0),hessian(0,0));
}

// *****************************************************************************
//
// ParameterizableSingleAxisFunction
//
// *****************************************************************************

ParameterizableMultiAxisFunction::~ParameterizableMultiAxisFunction()
{
  // nothing to see here
}

ParameterizableSingleAxisFunction::~ParameterizableSingleAxisFunction()
{
  // nothing to see here
}

double ParameterizableSingleAxisFunction::
value_and_parameter_gradient(ConstVecRef x, VecRef gradient)
{
  return value_and_parameter_gradient(x[0], gradient);
}

double ParameterizableSingleAxisFunction::
value_parameter_gradient_and_hessian(ConstVecRef x, VecRef gradient,
                                     MatRef hessian)
{
  return value_parameter_gradient_and_hessian(x[0], gradient, hessian);
}

// *****************************************************************************
//
// PMAFReverser
//
// *****************************************************************************

PMAFReverser::PMAFReverser(ParameterizableMultiAxisFunction* fcn_deligate,
                           bool adopt_fcn_deligate, double error_up):
    ParameterizableMultiAxisFunction(), fcn_deligate_(fcn_deligate),
    adopt_fcn_deligate_(adopt_fcn_deligate), error_up_(error_up),
    x_(fcn_deligate->num_domain_axes())
{
  // nothing to see here
}

PMAFReverser::~PMAFReverser()
{
  if(adopt_fcn_deligate_)delete fcn_deligate_;
}

unsigned PMAFReverser::num_parameters()
{
  return fcn_deligate_->num_domain_axes();
}

auto PMAFReverser::parameters() -> std::vector<function::ParameterAxis>
{
  return fcn_deligate_->domain_axes();
}

Eigen::VectorXd PMAFReverser::parameter_values()
{
  return x_;
}

void PMAFReverser::set_parameter_values(ConstVecRef values)
{
  x_ = values;
}

bool PMAFReverser::can_calculate_parameter_gradient()
{
  return fcn_deligate_->can_calculate_gradient();
}

bool PMAFReverser::can_calculate_parameter_hessian()
{
  return fcn_deligate_->can_calculate_hessian();
}

unsigned PMAFReverser::num_domain_axes()
{
  return fcn_deligate_->num_parameters();
}

auto PMAFReverser::domain_axes() -> std::vector<function::DomainAxis>
{
  return fcn_deligate_->parameters();
}

double PMAFReverser::value(ConstVecRef x)
{
  fcn_deligate_->set_parameter_values(x);
  return fcn_deligate_->value(x_);
}

bool PMAFReverser::can_calculate_gradient()
{
  return fcn_deligate_->can_calculate_parameter_gradient();
}

double PMAFReverser::value_and_gradient(ConstVecRef x, VecRef gradient)
{
  fcn_deligate_->set_parameter_values(x);
  return fcn_deligate_->value_and_parameter_gradient(x_, gradient);
}

bool PMAFReverser::can_calculate_hessian()
{
  return fcn_deligate_->can_calculate_parameter_hessian();
}

double PMAFReverser::value_gradient_and_hessian(ConstVecRef x, VecRef gradient,
                                                MatRef hessian)
{
  fcn_deligate_->set_parameter_values(x);
  return fcn_deligate_->value_parameter_gradient_and_hessian(x_, gradient,
                                                             hessian);
}

double PMAFReverser::error_up()
{
  return error_up_;
}

double PMAFReverser::value_and_parameter_gradient(ConstVecRef x,
                                                  VecRef gradient)
{
  fcn_deligate_->set_parameter_values(x);
  return fcn_deligate_->value_and_gradient(x_, gradient);
}

double PMAFReverser::value_parameter_gradient_and_hessian(ConstVecRef x,
                                           VecRef gradient, MatRef hessian)
{
  fcn_deligate_->set_parameter_values(x);
  return fcn_deligate_->value_gradient_and_hessian(x_, gradient, hessian);
}

// *****************************************************************************
//
// gradient_check
//
// *****************************************************************************

bool calin::math::function::
gradient_check(MultiAxisFunction& fcn, ConstVecRef x, VecRef good,
               double eps_factor)
{
  constexpr double eps = std::numeric_limits<double>::epsilon();
  constexpr double tiny_val = std::numeric_limits<double>::min();
  unsigned fcn_num_axes = fcn.num_domain_axes();
  Eigen::VectorXd dx(fcn_num_axes);
  for(unsigned idx=0;idx<fcn_num_axes;idx++)
    dx[idx] = std::max((std::abs(x[idx])*eps_factor)*eps, eps_factor*tiny_val);
  return gradient_check(fcn, x, dx, good);
}

bool calin::math::function::
gradient_check(MultiAxisFunction& fcn, ConstVecRef x, ConstVecRef dx,
               VecRef good)
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
  good.resize(fcn_num_axes);
  for(unsigned iaxis = 0;iaxis<fcn_num_axes;iaxis++)
  {
    // The philosophy used here is to calculate the gradient with the two-point
    // differnce formula, and compare that value to the expected numerical
    // error, which is: 1/6 h^2 d3f/dx3. We calculate the 3rd derivative with
    // the five point formula at the given point and at +/-h, and use the
    // maximum value found as the expected error (maxerr). If the difference
    // between the numerical and analytic gradient is less than maxerr we return
    // zero. If the difference is larger then we return 0.5*log10(diff/maxerr).
    // Users can test this to see how bad the computation is. Values up to 0.5
    // might be OK in practice.
    
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
    const double dfdx = (-fm1 + fp1)/h2;
    const double theerr = std::abs(dfdx - gradient(iaxis));

    // Coefficients from Abromowitz & Stegun Table 25.2
    const double h2d3fdx3 = (-fm2 + 2.0*fm1 - 2.0*fp1 + fp2)/h2;
    const double h2d3fdx3a = (-3.0*fm2 + 10.0*fm1 - 12.0*f0 + 6.0*fp1 - fp2)/h2;
    const double h2d3fdx3b = (fm2 - 6.0*fm1 + 12.0*f0 - 10.0*fp1 + 3.0*fp2)/h2;
    const double maxerr = std::max({std::abs(h2d3fdx3),std::abs(h2d3fdx3a),
                          std::abs(h2d3fdx3b)})/6.0;

    if(theerr < eps)
    {
      good(iaxis)=0;
    }
    else if(maxerr<eps)
    {
      // What to do if the 3rd derivative is itself zero?
      // Compare with epsilon - return zero if difference is less than
      // epsilon - grow to 0.5 (our nominal threshold) by sqrt(eps).
      if(theerr < eps)good(iaxis)=0.0;
      else good(iaxis)=std::log10(theerr)/std::log10(eps)-1;
    }
    else
    {
      // This is the main branch - compare difference to 3rd derivative,
      // return zero if smaller and grow to 0.5 (nominal threshold) by 10
      // times the expected error
      if(theerr < maxerr)good(iaxis)=0.0;
      else good(iaxis)=0.5*std::log10(theerr/(maxerr+eps)*(1.0+h+eps));
    }

    if(good(iaxis) > 0.5)
      std::cout << iaxis << ' '
                << "diff: " << dfdx << ' '
                << "grad: " << gradient(iaxis) << ' '
                << "err: " << std::abs(dfdx - gradient(iaxis)) << ' '
                << "maxerr: " << maxerr << ' '
                << "h2d3fdx3: " << h2d3fdx3 << ' '
                << "h2d3fdx3a: " << h2d3fdx3a << ' '
                << "h2d3fdx3b: " << h2d3fdx3b << ' '
                << "f(-2)-f(0): " << fm2-f0 <<  ' '
                << "f(-1)-f(0): " << fm1-f0 <<  ' '
                << "f(+1)-f(0): " << fp1-f0 <<  ' '
                << "f(+2)-f(0): " << fp1-f0 <<  ' '
                << "ret: " << good(iaxis) << ' '
                << '\n';
  }
  return true;
}

bool calin::math::function::
hessian_check(MultiAxisFunction& fcn, ConstVecRef x, ConstVecRef dx,
              MatRef good)
{
  constexpr double eps = std::numeric_limits<double>::epsilon();
  unsigned fcn_num_axes = fcn.num_domain_axes();
  assert(fcn_num_axes == x.innerSize());
  assert(fcn_num_axes == dx.innerSize());

  Eigen::VectorXd g0(fcn_num_axes);
  double f0 = fcn.value_and_gradient(x, g0);
  assert(fcn_num_axes == g0.innerSize());
  
  Eigen::VectorXd g0h(fcn_num_axes);
  Eigen::MatrixXd hessian(fcn_num_axes,fcn_num_axes);
  double f0h = fcn.value_gradient_and_hessian(x, g0h, hessian);
  assert(fcn_num_axes == g0h.innerSize());
  assert(fcn_num_axes == hessian.innerSize());
  assert(fcn_num_axes == hessian.outerSize());

  if(hessian != hessian.transpose())
  {
    std::cerr << "hessian_check: hessian matrix not symmetric\n";
    return false;
  }

  if(std::abs(f0 - f0h)/std::abs(f0) > eps)
  {
    std::cerr << "hessian_check: function value differs\n";
    return false;
  }
  
  for(unsigned iaxis = 0;iaxis<fcn_num_axes;iaxis++)
    if(std::abs(g0(iaxis) - g0h(iaxis))/std::abs(g0(iaxis)) > eps)
    {
      std::cerr << "hessian_check: gradient value differs (axis=" << iaxis
                << ")\n";
      return false;
    }

  good.resize(fcn_num_axes, fcn_num_axes);
  for(unsigned iaxis = 0;iaxis<fcn_num_axes;iaxis++)
  {
    // The philosophy used here is to calculate the hessian with the two-point
    // differnce formula, and compare that value to the expected numerical
    // error, which is: 1/6 h^2 d3f/dx3. We calculate the 3rd derivative with
    // the five point formula at the given point and at +/-h, and use the
    // maximum value found as the expected error (maxerr). If the difference
    // between the numerical and analytic hessian is less than maxerr we return
    // zero. If the difference is larger then we return 0.5*log10(diff/maxerr).
    // Users can test this to see how bad the computation is. Values up to 0.5
    // might be OK in practice.
    
    const double h = dx(iaxis);
    Eigen::VectorXd xh = x;
    xh(iaxis) = x(iaxis) + h;
    Eigen::VectorXd gp1(fcn_num_axes); fcn.value_and_gradient(xh,gp1);
    xh(iaxis) = x(iaxis) + 2.0*h;
    Eigen::VectorXd gp2(fcn_num_axes); fcn.value_and_gradient(xh,gp2);
    xh(iaxis) = x(iaxis) - h;
    Eigen::VectorXd gm1(fcn_num_axes); fcn.value_and_gradient(xh,gm1);
    xh(iaxis) = x(iaxis) - 2.0*h;
    Eigen::VectorXd gm2(fcn_num_axes); fcn.value_and_gradient(xh,gm2);
    volatile double h2v = x(iaxis) + h;
    h2v -= (x(iaxis) - h);
    const double h2 = h2v;
    
    Eigen::VectorXd dgdx = (-gm1 + gp1)/h2;
    Eigen::VectorXd theerr(fcn_num_axes);

    // Coefficients from Abromowitz & Stegun Table 25.2
    Eigen::VectorXd h2d3gdx3 = (-gm2 + 2.0*gm1 - 2.0*gp1 + gp2)/h2;
    Eigen::VectorXd h2d3gdx3a =
        (-3.0*gm2 + 10.0*gm1 - 12.0*g0 + 6.0*gp1 - gp2)/h2;
    Eigen::VectorXd h2d3gdx3b =
        (gm2 - 6.0*gm1 + 12.0*g0 - 10.0*gp1 + 3.0*gp2)/h2;

    for(unsigned jaxis = 0; jaxis<fcn_num_axes; jaxis++)
    {
      const double theerr = std::abs(dgdx(jaxis) - hessian(iaxis,jaxis));
      const double maxerr =
          std::max({std::abs(h2d3gdx3(jaxis)),std::abs(h2d3gdx3a(jaxis)),
                  std::abs(h2d3gdx3b(jaxis))})/6.0;
      
      if(theerr < eps)
      {
        good(iaxis,jaxis)=0;
      }
      else if(maxerr<eps)
      {
        // What to do if the 3rd derivative is itself zero?
        // Compare with epsilon - return zero if difference is less than
        // epsilon - grow to 0.5 (our nominal threshold) by sqrt(eps).
        if(theerr < eps)good(iaxis, jaxis)=0.0;
        else good(iaxis, jaxis)=std::log10(theerr)/std::log10(eps)-1;
      }
      else
      {
        // This is the main branch - compare difference to 3rd derivative,
        // return zero if smaller and grow to 0.5 (nominal threshold) by 10
        // times the expected error
        if(theerr < maxerr)good(iaxis, jaxis)=0.0;
        else good(iaxis, jaxis)=0.5*std::log10(theerr/(maxerr+eps)*(1.0+h+eps));
      }

      if(good(iaxis, jaxis) > 0.5)
        std::cout << iaxis << ' ' << jaxis << ' '
                  << "diff: " << dgdx(jaxis) << ' '
                  << "hess: " << hessian(iaxis, jaxis) << ' '
                  << "err: " << std::abs(dgdx(jaxis) - hessian(iaxis, jaxis)) << ' '
                  << "maxerr: " << maxerr << ' '
                  << "h2d3gdx3: " << h2d3gdx3(jaxis) << ' '
                  << "h2d3gdx3a: " << h2d3gdx3a(jaxis) << ' '
                  << "h2d3gdx3b: " << h2d3gdx3b(jaxis) << ' '
                  << "g(-2)-g(0): " << gm2(jaxis)-g0(jaxis) <<  ' '
                  << "g(-1)-g(0): " << gm1(jaxis)-g0(jaxis) <<  ' '
                  << "g(+1)-g(0): " << gp1(jaxis)-g0(jaxis) <<  ' '
                  << "g(+2)-g(0): " << gp1(jaxis)-g0(jaxis) <<  ' '
                  << "ret: " << good(iaxis, jaxis) << ' '
                  << '\n';
    }
  }
  return true;
}

