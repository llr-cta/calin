/* 

   calin/math/function.cpp -- Stephen Fegan -- 2015-02-24

*/

#include <iostream>

#include "math/function.hpp"

using namespace calin::math::function;

namespace {

inline double SQR(double x) { return x*x; }
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
// gradient_check
//
// *****************************************************************************

bool calin::math::function::
gradient_check(MultiAxisFunction& fcn, ConstVecRef x, VecRef err,
               double eps_factor = 10.0)
{
  constexpr double eps = std::numeric_limits<double>::epsilon();
  unsigned fcn_num_axes = fcn.num_domain_axes();
  Eigen::VectorXd dx(fcn_num_axes);
  for(unsigned idx=0;idx<fcn_num_axes;idx++)
    dx[idx] = (x[idx]*eps_factor)*eps;
  return gradient_check(fcn, x, dx, err);
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
      err(iaxis)=0;
    }
    else if(maxerr<eps)
    {
      // What to do if the 3rd derivative is itself zero?
      // Compare with epsilon - return zero if difference is less than
      // epsilon - grow to 0.5 (our nominal threshold) by sqrt(eps).
      if(theerr < eps)err(iaxis)=0.0;
      else err(iaxis)=std::log10(theerr)/std::log10(eps)-1;
    }
    else
    {
      // This is the main branch - compare difference to 3rd derivative,
      // return zero if smaller and grow to 0.5 (nominal threshold) by 10
      // times the expected error
      if(theerr < maxerr)err(iaxis)=0.0;
      else err(iaxis)=0.5*std::log10(theerr/(maxerr+eps)*(1.0+h+eps));
    }

    if(err(iaxis) > 0.5)
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
                << "ret: " << err(iaxis) << ' '
                << '\n';
  }
  return true;
}

// *****************************************************************************
//
// Miscellaneous functions
//
// *****************************************************************************

GaussianPDF::~GaussianPDF()
{
  // nothing to see here
}

unsigned GaussianPDF::num_parameters()
{
  return 2;
}

auto GaussianPDF::parameters() -> std::vector<math::ParameterAxis>
{
  constexpr double inf = std::numeric_limits<double>::infinity();
  constexpr double tiny_val = std::numeric_limits<double>::min();
  return { { "mean", "x-value units", false, -inf, false, inf },
    { "rms", "x-value units", true, tiny_val, false, inf } };
}

Eigen::VectorXd GaussianPDF::parameter_values()
{
  Eigen::VectorXd pval(2);
  pval << x0_, s_;
  return pval;
}

void GaussianPDF::set_parameter_values(ConstVecRef values)
{
  assign_parameters(values, x0_, s_);
}

bool GaussianPDF::can_calculate_gradient()
{
  return true;
}

bool GaussianPDF::can_calculate_hessian()
{
  return true;
}

bool GaussianPDF::can_calculate_parameter_gradient()
{
  return true;
}
  
bool GaussianPDF::can_calculate_parameter_hessian()
{
  return true;
}

DomainAxis GaussianPDF::domain_axis()
{
  constexpr double inf = std::numeric_limits<double>::infinity();
  return { "x-value", "x-value units", false, -inf, false, inf };
}

double GaussianPDF::value(double x) 
{
  return c_gauss_norm/s_*std::exp(-SQR(x-x0_)/(2.0*SQR(s_)));
}
    
double GaussianPDF::value_and_gradient(double x,  double& dfdx)
{
  const double xc = x-x0_;
  const double xs = xc/s_;
  const double xs2 = SQR(xs);
  double val = c_gauss_norm/s_*exp(-0.5*xs2);
  dfdx = -val*xs/s_;
  return val;
}

double GaussianPDF::
value_gradient_and_hessian(double x, double& dfdx, double& d2fdx2)
{
  const double xc = x-x0_;
  const double xs = xc/s_;
  const double xs2 = SQR(xs);
  double val = c_gauss_norm/s_*exp(-0.5*xs2);
  dfdx = -val*xs/s_;
  d2fdx2 = -(dfdx*xc + val)/SQR(s_);
  return val;
}

double GaussianPDF::
value_and_parameter_gradient(double x,  VecRef gradient)
{
  const double xc = x-x0_;
  const double xs = xc/s_;
  const double xs2 = SQR(xs);
  double val = c_gauss_norm/s_*exp(-0.5*xs2);
  gradient.resize(2);
  gradient[0] = val*xs/s_; // df/dx0
  gradient[1] = val*(xs2 - 1.0)/s_;
  return val;
}

double GaussianPDF::
value_parameter_gradient_and_hessian(double x, VecRef gradient, MatRef hessian)
{
  const double xc = x-x0_;
  const double xs = xc/s_;
  const double xs2 = SQR(xs);
  double val = c_gauss_norm/s_*exp(-0.5*xs2);
  gradient.resize(2);
  gradient[0] = val*xs/s_; // df/dx0
  gradient[1] = val*(xs2 - 1.0)/s_;
  hessian(0,0) = gradient[1]/s_; // val*(xs2 - 1.0)/SQR(s_);
  hessian(0,1) = val*(xs2 - 3.0)*xs/SQR(s_);
  hessian(1,0) = hessian(0,1);
  hessian(1,1) = val*(SQR(xs2) - 5.0*xs2 + 3.0)/SQR(s_);
  return val;
}
