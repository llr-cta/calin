/* 

   calin/math/log_quadratic_spline_pdf_1d.cpp -- Stephen Fegan -- 2015-07-30

   PDF base on quadratic spline in log space.

*/

#include <cmath>
#include <iostream>

#include "math/log_quadratic_spline_pdf_1d.hpp"

using namespace calin::math::pdf_1d;

LogQuadraticSpline1DPDF::
LogQuadraticSpline1DPDF(ConstVecRef xknot, double xlo, double xhi, bool normalize,
                        const std::string& yunits, const std::string& xunits,
                        double error_up):
    Parameterizable1DPDF(), yunits_(yunits), xunits_(xunits), error_up_(error_up),
    xlo_(xlo), xhi_(xhi), normalize_(normalize),
    nknot_(xknot.size()),
    xknot_(xknot), yknot_(Eigen::VectorXd::Zero(nknot_)),
    dx_(std::max(1U,nknot_)-1), dy_(std::max(1U,nknot_)-1),
    a_(std::max(1U,nknot_)-1), b_(std::max(1U,nknot_)-1),
    a_gradient_(nknot_+1, std::max(1U,nknot_)-1),
    b_gradient_(nknot_+1, std::max(1U,nknot_)-1),
    norm_gradient_(Eigen::VectorXd::Zero(nknot_+1))
{
  if(nknot_<2)
  {
    std::ostringstream stream;
    stream << "LogQuadraticSpline1DPDF - minimum of 2 knots required, " << nknot_
           << " supplied.";
    throw(std::runtime_error(stream.str()));
  }  
  
  std::sort(xknot_.data(), xknot_.data()+nknot_);
  dx_ = xknot_.tail(nknot_-1) - xknot_.head(nknot_-1);

  if(xknot_(0) < xlo_)
  {
    std::ostringstream stream;
    stream << "LogQuadraticSpline1DPDF - first knot is outside limits, " << xknot_(0)
           << " < " << xlo_;
    throw(std::out_of_range(stream.str()));
  }  

  if(xknot_(nknot_-1) > xhi_)
  {
    std::ostringstream stream;
    stream << "LogQuadraticSpline1DPDF - final knot is outside limits, " << xknot_(nknot_-1)
           << " > " << xhi_;
    throw(std::out_of_range(stream.str()));
  }  
  
  a_gradient_.setZero();
  b_gradient_.setZero();
  
  set_cache();
}
      
LogQuadraticSpline1DPDF::~LogQuadraticSpline1DPDF()
{
  // nothing to see here
}
  
unsigned LogQuadraticSpline1DPDF::num_parameters()
{
  return 1+nknot_;
}

std::vector<calin::math::function::ParameterAxis> LogQuadraticSpline1DPDF::parameters()
{
  std::vector<calin::math::function::ParameterAxis> axes;
  axes.push_back({ std::string("dy_dx(")+std::to_string(xknot_(0))+std::string(")"),
          yunits_+std::string("/")+xunits_, false, -inf, false, inf, 1, 0 });
  for(unsigned iknot=0; iknot<nknot_; iknot++)
    axes.push_back({ std::string("y(")+std::to_string(xknot_(iknot))+std::string(")"),
            yunits_, false, -inf, false, inf, 1, 0 });
  return axes;
}

Eigen::VectorXd LogQuadraticSpline1DPDF::parameter_values()
{
  Eigen::VectorXd p(nknot_+1);
  p[0] = slope0_;
  p.tail(nknot_) = yknot_;
  return p;
}

void LogQuadraticSpline1DPDF::set_parameter_values(ConstVecRef values)
{
  if(values.size() != num_parameters())
  {
    std::ostringstream stream;
    stream << "LogQuadraticSpline1DPDF - parameter vector has " << values.size()
           << " values, " << num_parameters() << " required.";
    throw(std::runtime_error(stream.str()));
  }  
  slope0_ = values(0);
  yknot_ = values.tail(nknot_);
  set_cache();
}

calin::math::function::DomainAxis LogQuadraticSpline1DPDF::domain_axis()
{
  constexpr double inf = std::numeric_limits<double>::infinity();
  return { "x-value", xunits_, xlo_>-inf, xlo_, xhi_<inf, xhi_ };
}

bool LogQuadraticSpline1DPDF::can_calculate_gradient()
{
  return true;
}

bool LogQuadraticSpline1DPDF::can_calculate_hessian()
{
  return false;
}

bool LogQuadraticSpline1DPDF::can_calculate_parameter_gradient()
{
  return true;
}

bool LogQuadraticSpline1DPDF::can_calculate_parameter_hessian()
{
  return false;
}

double LogQuadraticSpline1DPDF::value_1d(double x)
{
  if(x<xlo_ or x>=xhi_)return 0;
  unsigned isegment = find_segment(x);
  double xx = x-xknot_(isegment);
#if 0
  std::cout << ">> " << x << ' ' << isegment << ' ' << xx << ' ' << a_(isegment) << ' '
            << b_(isegment) << ' ' << yknot_(isegment) << "<< ";
#endif
  return norm_*std::exp((a_(isegment)*xx+b_(isegment))*xx+yknot_(isegment));
}

double LogQuadraticSpline1DPDF::value_and_gradient_1d(double x,  double& dfdx)
{
  if(x<xlo_ or x>=xhi_)
  {
    dfdx = 0;
    return 0;
  }
  unsigned isegment = find_segment(x);
  double xx = x-xknot_(isegment);
  double val = (a_(isegment)*xx+b_(isegment))*xx+yknot_(isegment);
  dfdx = (2.0*a_(isegment)*xx+b_(isegment));
  val = norm_*std::exp(val);
  dfdx *= val;
  return val;
}

double LogQuadraticSpline1DPDF::value_gradient_and_hessian_1d(double x, double& dfdx,
                                                              double& d2fdx2)
{
  throw std::runtime_error("LogQuadraticSpline1DPDF::value_gradient_and_hessian_1d not implemented");
}

double LogQuadraticSpline1DPDF::value_and_parameter_gradient_1d(double x, VecRef gradient)
{
  gradient.resize(num_parameters());
  if(x<xlo_ or x>=xhi_)
  {
    gradient.setZero();
    return 0;
  }
  unsigned isegment = find_segment(x);
  double xx = x-xknot_(isegment);
  double val = (a_(isegment)*xx+b_(isegment))*xx+yknot_(isegment);
  gradient = (a_gradient_.col(isegment)*xx + b_gradient_.col(isegment))*xx;
  gradient(isegment+1) += 1.0;
  val = norm_ * std::exp(val);
#if 0
  std::cout << "AAA: [" 
            << a_gradient_.col(isegment).transpose() << " ] [ "
            << b_gradient_.col(isegment).transpose() << " ] [ "
            << gradient.transpose() << " ] [ "
            << norm_gradient_.transpose() << " ] \n";
#endif
  gradient += norm_gradient_;
  gradient *= val;
  return val;
}

double LogQuadraticSpline1DPDF::
value_parameter_gradient_and_hessian_1d(double x, VecRef gradient, MatRef hessian)
{
  throw std::runtime_error("LogQuadraticSpline1DPDF::value_parameter_gradient_and_hessian_1d not implemented");  
}

double LogQuadraticSpline1DPDF::error_up()
{
  return error_up_;
}

bool LogQuadraticSpline1DPDF::can_calculate_mean_and_variance()
{
  return false;
}

void LogQuadraticSpline1DPDF::mean_and_variance(double& mean, double& var)
{
  throw std::runtime_error("LogQuadraticSpline1DPDF::mean_and_variance not implemented");
}

unsigned LogQuadraticSpline1DPDF::find_segment(double x) const
{
  unsigned ix = std::lower_bound(xknot_.data(), xknot_.data()+nknot_-1, x)-xknot_.data();
  return std::max(1U,ix)-1;
}

#include <gsl/gsl_sf_dawson.h>

namespace {

inline double SQR(double x) { return x*x; }
inline double CUB(double x) { return x*x*x; }

double dawson(double x)
{
  return gsl_sf_dawson(x);
}

} // anonymous namespace

void LogQuadraticSpline1DPDF::set_cache()
{
  dy_ = yknot_.tail(nknot_-1)-yknot_.head(nknot_-1);
  b_(0) = slope0_;
  a_(0) = (dy_(0) - dx_(0) * b_(0))/SQR(dx_(0));

  b_gradient_(0,0) = 1.0;
  a_gradient_(0,0) = -1.0/dx_(0);
  a_gradient_(1,0) = -1.0/SQR(dx_(0));
  a_gradient_(2,0) = 1.0/SQR(dx_(0));
  
  for(unsigned isegment=1; isegment<nknot_-1; isegment++)
  {
    double dx2 = SQR(dx_(0));
    b_(isegment) = 2.0*dx_(isegment-1)*a_(isegment-1) + b_(isegment-1);
    a_(isegment) = (dy_(isegment) - dx_(isegment)*b_(isegment))/dx2;

    b_gradient_.col(isegment) =
        2.0*dx_(isegment-1)*a_gradient_.col(isegment-1) + b_gradient_.col(isegment-1);
    a_gradient_.col(isegment) =
        -b_gradient_.col(isegment)/dx_(isegment);
    a_gradient_(isegment+2, isegment) += 1.0/dx2;
    a_gradient_(isegment+1, isegment) -= 1.0/dx2;
  }

  if(normalize_)
  {
    norm_ = 0;
    norm_gradient_.setZero();
    for(unsigned isegment = 0; isegment<nknot_-1; isegment++)
    {
      double xl = xknot_[isegment];
      double xr = xknot_[isegment+1];
      if(isegment == 0)xl=std::min(xl,xlo_);
      if(isegment == nknot_-2)xr=std::max(xr,xhi_);
      xl -= xknot_[isegment];
      xr -= xknot_[isegment];
      const double a = a_[isegment];
      const double b = b_[isegment];
      const double c = yknot_[isegment];
#if 0
      std::cout << a << ' ' << b << ' ' << c << '\n';
#endif
      double I;
      double dI_da;
      double dI_db;
      
      if(a == 0)
      {
        if(b == 0)
        {
          const double F_exp = std::exp(c);
          I = F_exp*(xr-xl);
          dI_da = F_exp * (CUB(xr) - CUB(xl))/3;
          dI_db = F_exp * (SQR(xr) - SQR(xl))/2;
        }
        else
        {
          const double F_exp_r = std::exp(b*xr+c);
          const double F_exp_l = std::exp(b*xl+c);
          I = (F_exp_r - F_exp_l)/b;
          dI_da = ((b*xr*(b*xr-2)+2)*F_exp_r
                   - (b*xl*(b*xl-2)+2)*F_exp_l)/CUB(b);
          dI_db = ((b*xr-1)*F_exp_r - (b*xl-1)*F_exp_l)/SQR(b);
        }
      }
      else // a!=0
      {
        const double F_exp = std::exp(-SQR(b)/(4*a)+c);
        const double F_exp_r = std::exp(a*SQR(xr+b/(2*a)));
        const double F_exp_l = std::exp(a*SQR(xl+b/(2*a)));

        if(a < 0) // negative curvature - integral is erf
        {
          const double s = std::sqrt(-a);
          I = 1/(M_2_SQRTPI*s)*F_exp*
              (std::erf(s*(xr+b/(2*a)))-std::erf(s*(xl+b/(2*a))));
        }
        else // positive curvature - integral is e^(x**2) dawson(x) 
        {
          const double s = std::sqrt(a);
          I = 1/s*F_exp*(F_exp_r*dawson(s*(xr+b/(2*a)))
                         - F_exp_l*dawson(s*(xl+b/(2*a))));
        }

        dI_da = I*(SQR(b)/(4*SQR(a)) - 1/(2*a))
                + 1/(2*a)*F_exp*(F_exp_r*(xr-b/(2*a)) - F_exp_l*(xl-b/(2*a)));
        dI_db = (-I*b + F_exp*(F_exp_r - F_exp_l))/(2*a);
      }

      norm_ += I;
      norm_gradient_ += dI_da * a_gradient_.col(isegment);
      norm_gradient_ += dI_db * b_gradient_.col(isegment);
      norm_gradient_(1+isegment) += I;
    }
    
    norm_ = 1/norm_;
    norm_gradient_ *= -norm_;
  }
}
