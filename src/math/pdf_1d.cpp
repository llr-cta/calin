/* 

   calin/math/pdf_1d.cpp -- Stephen Fegan -- 2015-04-02

*/

#include <iostream>
#include <iomanip>
#include <stdexcept>

#include "math/pdf_1d.hpp"

using namespace calin::math;
using namespace calin::math::pdf_1d;

namespace {

inline static double SQR(double x) { return x*x; }
constexpr double c_gauss_norm = 0.5*M_2_SQRTPI*M_SQRT1_2;

} // anonymous namespace

using function::assign_parameters;

Parameterizable1DPDF::~Parameterizable1DPDF()
{
  // nothing to see here
}

// *****************************************************************************
//
// GaussianPDF
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
  const double xc = x-x0_;
  const double xs = xc/s_;
  const double xs2 = SQR(xs);
  double val = c_gauss_norm/s_*exp(-0.5*xs2);
  //std::cout << std::setprecision(17) << x << ' ' << val << '\n';
  return val;
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
  hessian(1,1) = val*(SQR(xs2) - 5.0*xs2 + 2.0)/SQR(s_);
  return val;
}

double GaussianPDF::error_up()
{
  return error_up_;
}

bool GaussianPDF::can_calculate_mean_and_variance()
{
  return true;
}

void GaussianPDF::get_mean_and_variance(double& mean, double& var)
{
  mean = x0_;
  var = SQR(s_);
}

// *****************************************************************************
//
// LimitedGaussianPDF
//
// *****************************************************************************

LimitedGaussianPDF::~LimitedGaussianPDF()
{
  // nothing to see here
}

DomainAxis LimitedGaussianPDF::domain_axis()
{
  DomainAxis a = GaussianPDF::domain_axis();
  a.has_lo_bound = xlo_>-inf;
  a.lo_bound = xlo_;
  a.has_hi_bound = xhi_<inf;
  a.hi_bound = xhi_;
  return a;
}

void LimitedGaussianPDF::set_parameter_values(ConstVecRef values)
{
  GaussianPDF::set_parameter_values(values);
  set_cache();
}

double LimitedGaussianPDF::value(double x)
{
  if(x<xlo_ or x>=xhi_)return 0;
  return norm_*GaussianPDF::value(x);
}
    
double LimitedGaussianPDF::value_and_gradient(double x,  double& dfdx)
{
  if(x<xlo_ or x>=xhi_)
  {
    dfdx = 0;
    return 0;
  }
  double val = norm_*GaussianPDF::value_and_gradient(x, dfdx);
  dfdx *= norm_;
  return val;
}

double LimitedGaussianPDF::value_gradient_and_hessian(double x, double& dfdx,
                                                      double& d2fdx2)
{
  if(x<xlo_ or x>=xhi_)
  {
    d2fdx2 = 0;
    dfdx = 0;
    return 0;
  }
  double val = GaussianPDF::value_gradient_and_hessian(x, dfdx, d2fdx2);
  dfdx *= norm_;
  d2fdx2 *= norm_;
  return norm_*val;
}

double LimitedGaussianPDF::
value_and_parameter_gradient(double x,  VecRef gradient)
{
  if(x<xlo_ or x>=xhi_)
  {
    gradient[0] = gradient[1] = 0;
    return 0;
  }  
  double val = GaussianPDF::value_and_parameter_gradient(x, gradient);
  gradient = norm_*gradient + val*norm_gradient_;
  return norm_*val;
}

double LimitedGaussianPDF::
value_parameter_gradient_and_hessian(double x, VecRef gradient, MatRef hessian)
{
  if(x<xlo_ or x>=xhi_)
  {
    gradient[0] = gradient[1] = 0;
    hessian(0,0) = hessian(0,1) = hessian(1,0) = hessian(1,1) = 0;
    return 0;
  }  
  double val =
      GaussianPDF::value_parameter_gradient_and_hessian(x, gradient, hessian);
  hessian = norm_*hessian + val*norm_hessian_;
  hessian(0,0) += 2.0*norm_gradient_[0]*gradient[0];
  hessian(0,1) += norm_gradient_[0]*gradient[1] + norm_gradient_[1]*gradient[0];
  hessian(1,0) = hessian(0,1);
  hessian(1,1) += 2.0*norm_gradient_[1]*gradient[1];
  gradient = norm_*gradient + val*norm_gradient_;
  return norm_*val;
}

bool LimitedGaussianPDF::can_calculate_mean_and_variance()
{
  return false;
}

void LimitedGaussianPDF::get_mean_and_variance(double& mean, double& var)
{
  assert(0);
}

void LimitedGaussianPDF::set_cache()
{
  double ehi         = 1.0;
  double dehi_dx0    = 0.0;
  double dehi_ds     = 0.0;
  double d2ehi_dx02  = 0.0;
  double d2ehi_ds2   = 0.0;
  double d2ehi_dx0ds = 0.0;
  
  if(xhi_!=inf)
  {
    double xc = xhi_-x0_;
    double xs = xc/s_;
    double xs2 = SQR(xs);

    ehi         = 0.5*(1.0+std::erf(M_SQRT1_2*xs));

    dehi_dx0    = -c_gauss_norm/s_*exp(-0.5*xs2);
    dehi_ds     = xs*dehi_dx0;

    d2ehi_dx02  = dehi_ds/s_;
    d2ehi_dx0ds = dehi_dx0*(xs2-1)/s_;
    d2ehi_ds2   = xs*d2ehi_dx0ds - d2ehi_dx02;
  }
    
  double elo         = 0.0;
  double delo_dx0    = 0.0;
  double delo_ds     = 0.0;
  double d2elo_dx02  = 0.0;
  double d2elo_ds2   = 0.0;
  double d2elo_dx0ds = 0.0;
  
  if(xlo_!=-inf)
  {
    double xc = xlo_-x0_;
    double xs = xc/s_;
    double xs2 = SQR(xs);

    elo         = 0.5*(1.0+std::erf(M_SQRT1_2*xs));

    delo_dx0    = -c_gauss_norm/s_*exp(-0.5*xs2);
    delo_ds     = xs*delo_dx0;

    d2elo_dx02  = delo_ds/s_;
    d2elo_dx0ds = delo_dx0*(xs2-1)/s_;
    d2elo_ds2   = xs*d2elo_dx0ds - d2elo_dx02;
  }

  norm_              = 1.0/(ehi - elo);
  const double norm2 = SQR(norm_);
  norm_gradient_(0)  = -norm2*(dehi_dx0 - delo_dx0);
  norm_gradient_(1)  = -norm2*(dehi_ds - delo_ds);
  norm_hessian_(0,0) = -norm2*(d2ehi_dx02 - d2elo_dx02)
                       + 2.0*SQR(norm_gradient_(0))/norm_;
  norm_hessian_(0,1) = -norm2*(d2ehi_dx0ds - d2elo_dx0ds)
                       + 2.0*norm_gradient_(0)*norm_gradient_(1)/norm_;
  norm_hessian_(1,0) = norm_hessian_(0,1);
  norm_hessian_(1,1) = -norm2*(d2ehi_ds2 - d2elo_ds2)
                       + 2.0*SQR(norm_gradient_(1))/norm_;
}

// *****************************************************************************
//
// LimitedExponentialPDF
//
// *****************************************************************************

LimitedExponentialPDF::~LimitedExponentialPDF()
{
  // nothing to see here
}

unsigned LimitedExponentialPDF::num_parameters()
{
  return 1;
}

std::vector<ParameterAxis> LimitedExponentialPDF::parameters()
{
  return { { "scale", "x-value units", false, -inf, false, inf } };
}

Eigen::VectorXd LimitedExponentialPDF::parameter_values()
{
  Eigen::VectorXd p(1);
  p << a_;
  return p;
}

void LimitedExponentialPDF::set_parameter_values(ConstVecRef values)
{
  assign_parameters(values, a_);
  set_cache();
}

DomainAxis LimitedExponentialPDF::domain_axis()
{
  return { "x-value", "x-value units", xlo_>-inf, xlo_, xhi_<inf, xhi_ };
}

bool LimitedExponentialPDF::can_calculate_gradient()
{
  return true;
}

bool LimitedExponentialPDF::can_calculate_hessian()
{
  return true;
}

bool LimitedExponentialPDF::can_calculate_parameter_gradient()
{
  return true;
}

bool LimitedExponentialPDF::can_calculate_parameter_hessian()
{
  return true;
}

double LimitedExponentialPDF::value(double x)
{
  if(x<xlo_ or x>=xhi_)return 0;
  const double xs = x/a_;
  return norm_ * std::exp(-xs);
}

double LimitedExponentialPDF::value_and_gradient(double x,  double& dfdx)
{
  if(x<xlo_ or x>=xhi_)
  {
    dfdx = 0;
    return 0;
  }
  const double xs = x/a_;
  const double val = norm_*std::exp(-xs);
  dfdx = -val/a_;
  return val;
}

double LimitedExponentialPDF::
value_gradient_and_hessian(double x, double& dfdx, double& d2fdx2)
{
  if(x<xlo_ or x>=xhi_)
  {
    dfdx = 0;
    d2fdx2 = 0;
    return 0;
  }
  const double xs = x/a_;
  const double val = norm_*std::exp(-xs);
  dfdx = -val/a_;
  d2fdx2 = val/SQR(a_);
  return val;
}

double LimitedExponentialPDF::
value_and_parameter_gradient(double x,  VecRef gradient)
{
  if(x<xlo_ or x>=xhi_)
  {
    gradient(0) = 0;
    return 0;
  }
  const double xs = x/a_;
  const double val = std::exp(-xs);
  gradient(0) = val*(norm_gradient_ + norm_*xs/a_);
  return norm_*val;
}

double LimitedExponentialPDF::
value_parameter_gradient_and_hessian(double x, VecRef gradient, MatRef hessian)
{
  if(x<xlo_ or x>=xhi_)
  {
    gradient(0) = 0;
    hessian(0,0) = 0;
    return 0;
  }
  const double xs = x/a_;
  const double val = std::exp(-xs);
  gradient(0) = val*(norm_gradient_ + norm_*xs/a_);  
  hessian(0,0) = val*(norm_hessian_ + 2.0*norm_gradient_*xs/a_
                      + norm_*xs*(xs - 2.0)/SQR(a_));
  return norm_*val;
}

double LimitedExponentialPDF::error_up()
{
  return error_up_;
}

bool LimitedExponentialPDF::can_calculate_mean_and_variance()
{
  return false;
}

void LimitedExponentialPDF::get_mean_and_variance(double& mean, double& var)
{
  assert(0);
}

void LimitedExponentialPDF::set_cache()
{
  double ehi         = 0.0;
  double dehi_da     = 0.0;
  double d2ehi_da2   = 0.0;
  
  if(xhi_!=inf)
  {
    double xs = xhi_/a_;
    double exs = exp(-xs);
    
    ehi         = exs*a_;
    dehi_da     = exs*(1.0 + xs);
    d2ehi_da2   = exs*SQR(xs)/a_;
  }
  else if(a_<=0)
  {
    throw std::out_of_range("LimitedExponentialPDF: scale must be strictly "
                            "positive when xhi=inf");
  }
    
  double elo         = 0.0;
  double delo_da     = 0.0;
  double d2elo_da2   = 0.0;
  
  if(xlo_!=-inf)
  {
    double xs = xlo_/a_;
    double exs = exp(-xs);
    
    elo         = exs*a_;
    delo_da     = exs*(1.0 + xs);
    d2elo_da2   = exs*SQR(xs)/a_;
  }
  else if(a_>=0)
  {
    throw std::out_of_range("LimitedExponentialPDF: scale must be strictly "
                            "negative when xlo=-inf");
  }

  norm_              = 1.0/(ehi - elo);
  const double norm2 = SQR(norm_);
  norm_gradient_     = -norm2*(dehi_da - delo_da);
  norm_hessian_      = -norm2*(d2ehi_da2 - d2elo_da2)
                       + 2.0*SQR(norm_gradient_)/norm_;
}
