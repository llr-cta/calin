/*

   calin/math/log_quadratic_spline_pdf_1d.cpp -- Stephen Fegan -- 2015-07-30

   PDF based on quadratic spline in log space.

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
   Laboratoire Leprince-Ringuet, CNRS/IN2P3, Ecole Polytechnique, Institut Polytechnique de Paris

   This file is part of "calin"

   "calin" is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 or
   later, as published by the Free Software Foundation.

   "calin" is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

*/

#include <cmath>
#include <iostream>

#include <math/special.hpp>
#include <math/log_quadratic_spline_pdf_1d.hpp>

using namespace calin::math::special;
using namespace calin::math::pdf_1d;

LogQuadraticSpline1DPDF::
LogQuadraticSpline1DPDF(ConstVecRef xknot, double xlo, double xhi,
                        double bin_dx,
                        ParamZeroType p0_type, ParamZeroLocation p0_loc,
                        bool normalize,
                        const std::string& yunits, const std::string& xunits):
    Parameterizable1DPDF(), yunits_(yunits), xunits_(xunits),
    xlo_(xlo), xhi_(xhi), bin_dx_(bin_dx),
    p0_type_(p0_type), p0_loc_(p0_loc),
    nknot_(xknot.size()), xknot_(xknot), yknot_(Eigen::VectorXd::Zero(nknot_)),
    dx_(std::max(1U,nknot_)-1), dy_(std::max(1U,nknot_)-1),
    a_(std::max(1U,nknot_)-1), b_(std::max(1U,nknot_)-1),
    a_gradient_(nknot_+1, std::max(1U,nknot_)-1),
    b_gradient_(nknot_+1, std::max(1U,nknot_)-1),
    normalize_(normalize), norm_gradient_(Eigen::VectorXd::Zero(nknot_+1))
{
  if(nknot_<2)
  {
    std::ostringstream stream;
    stream << "LogQuadraticSpline1DPDF - minimum of 2 knots required, "
           << nknot_ << " supplied.";
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

std::vector<calin::math::function::ParameterAxis>
LogQuadraticSpline1DPDF::parameters()
{
  std::vector<calin::math::function::ParameterAxis> axes;
  std::string p0_name;
  if(p0_type_ == ParamZeroType::SLOPE)p0_name = "dy/dx(";
  else p0_name = "d^2y/2dx^2(";
  if(p0_loc_ == ParamZeroLocation::RIGHT)
    p0_name += std::to_string(xknot_(xknot_.size()-1));
  else p0_name += std::to_string(xknot_(0));
  p0_name += std::string(")");
  std::string p0_denom_sq;
  if(p0_type_ != ParamZeroType::SLOPE)p0_denom_sq = "**2";
  axes.push_back({ p0_name, yunits_+std::string("_")+xunits_+p0_denom_sq,
          false, -inf, false, inf, 1, 0 });
  for(unsigned iknot=0; iknot<nknot_; iknot++)
    axes.push_back({ std::string("y(")+std::to_string(xknot_(iknot))+std::string(")"),
            yunits_, false, -inf, false, inf, 1, 0 });
  return axes;
}

Eigen::VectorXd LogQuadraticSpline1DPDF::parameter_values()
{
  Eigen::VectorXd p(nknot_+1);
  p[0] = param0_;
  p.tail(nknot_) = yknot_;
  return p;
}

void LogQuadraticSpline1DPDF::set_parameter_values(ConstVecRef values)
{
  verify_set_parameter_values(values, "LogQuadraticSpline1DPDF");
  param0_ = values(0);
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
  return bin_dx_ == 0;
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
  const double xl = std::max(x-0.5*bin_dx_, xlo_);
  const double xr = std::min(x+0.5*bin_dx_, xhi_);

  if(xr<xlo_ or xl>xhi_ or norm_<=0.0 or !std::isfinite(norm_))return 0;
  unsigned isegment = find_segment(x);

  if(bin_dx_ == 0)
  {
    double xx = x-xknot_(isegment);
#if 0
    std::cout << ">> " << x << ' ' << isegment << ' ' << xx << ' '
              << a_(isegment) << ' ' << b_(isegment) << ' '
              << yknot_(isegment) << "<< ";
#endif
    return norm_*std::exp((a_(isegment)*xx+b_(isegment))*xx+yknot_(isegment));
  }
  else
  {
    double I = 0;
    double dI_da = 0;
    double dI_db = 0;
    if(isegment>0 and xl<xknot_(isegment))
    {
      // left edge of interval crosses segment boundary
      integral(a_(isegment-1), b_(isegment-1), yknot_(isegment-1),
               xl-xknot_(isegment-1), dx_(isegment-1), I, dI_da, dI_db);
      double I2 = 0;
      integral(a_(isegment), b_(isegment), yknot_(isegment),
               0, xr-xknot_(isegment), I2, dI_da, dI_db);
      I += I2;
    }
    else if(isegment<nknot_-2 and xr>xknot_(isegment+1))
    {
      // right edge of interval crosses segment boundary
      integral(a_(isegment), b_(isegment), yknot_(isegment),
               xl-xknot_(isegment), dx_(isegment), I, dI_da, dI_db);
      double I2 = 0;
      integral(a_(isegment+1), b_(isegment+1), yknot_(isegment+1),
               0, xr-xknot_(isegment+1), I2, dI_da, dI_db);
      I += I2;
    }
    else
    {
      integral(a_(isegment), b_(isegment), yknot_(isegment),
               xl-xknot_(isegment), xr-xknot_(isegment), I, dI_da, dI_db);
    }
    return norm_ * I / bin_dx_;
  }
}

double LogQuadraticSpline1DPDF::value_and_gradient_1d(double x,  double& dfdx)
{
  if(x<xlo_ or x>=xhi_ or norm_<=0.0 or !std::isfinite(norm_))
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

double LogQuadraticSpline1DPDF::
value_and_parameter_gradient_1d(double x, VecRef gradient)
{
  gradient.resize(num_parameters());

  const double xl = std::max(x-0.5*bin_dx_, xlo_);
  const double xr = std::min(x+0.5*bin_dx_, xhi_);

  if(xr<xlo_ or xl>xhi_ or norm_<=0.0 or !std::isfinite(norm_))
  {
    gradient.setZero();
    return 0;
  }

  unsigned isegment = find_segment(x);

  if(bin_dx_ == 0)
  {
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
  else
  {
    double I = 0;
    double dI_da = 0;
    double dI_db = 0;
    if(isegment>0 and xl<xknot_(isegment))
    {
      // left edge of interval crosses segment boundary
      integral(a_(isegment-1), b_(isegment-1), yknot_(isegment-1),
               xl-xknot_(isegment-1), dx_(isegment-1), I, dI_da, dI_db);
      gradient = dI_da * a_gradient_.col(isegment-1);
      gradient += dI_db * b_gradient_.col(isegment-1);
      gradient(isegment) += I;

      double I2 = 0;
      integral(a_(isegment), b_(isegment), yknot_(isegment),
               0, xr-xknot_(isegment), I2, dI_da, dI_db);
      I += I2;
      gradient += dI_da * a_gradient_.col(isegment);
      gradient += dI_db * b_gradient_.col(isegment);
      gradient(1+isegment) += I2;
    }
    else if(isegment<nknot_-2 and xr>xknot_(isegment+1))
    {
      // right edge of interval crosses segment boundary
      integral(a_(isegment), b_(isegment), yknot_(isegment),
               xl-xknot_(isegment), dx_(isegment), I, dI_da, dI_db);
      gradient = dI_da * a_gradient_.col(isegment);
      gradient += dI_db * b_gradient_.col(isegment);
      gradient(1+isegment) += I;

      double I2 = 0;
      integral(a_(isegment+1), b_(isegment+1), yknot_(isegment+1),
               0, xr-xknot_(isegment+1), I2, dI_da, dI_db);
      I += I2;
      gradient += dI_da * a_gradient_.col(isegment+1);
      gradient += dI_db * b_gradient_.col(isegment+1);
      gradient(2+isegment) += I2;
    }
    else
    {
      integral(a_(isegment), b_(isegment), yknot_(isegment),
               xl-xknot_(isegment), xr-xknot_(isegment), I, dI_da, dI_db);
      gradient = dI_da * a_gradient_.col(isegment);
      gradient += dI_db * b_gradient_.col(isegment);
      gradient(1+isegment) += I;
    }


    gradient *= norm_/bin_dx_;
    I *= norm_/bin_dx_;
    gradient += norm_gradient_*I;
    return I;
  }
}

double LogQuadraticSpline1DPDF::
value_parameter_gradient_and_hessian_1d(double x, VecRef gradient, MatRef hessian)
{
  throw std::runtime_error("LogQuadraticSpline1DPDF::value_parameter_gradient_and_hessian_1d not implemented");
}

Eigen::VectorXd LogQuadraticSpline1DPDF::a_gradient(unsigned isegment) const
{
  if(isegment > nknot_-2)
  {
    std::ostringstream stream;
    stream << "LogQuadraticSpline1DPDF::a_gradient: isegment out of range: "
           << isegment << " > " << nknot_-2;
    throw(std::out_of_range(stream.str()));
  }
  return a_gradient_.col(isegment);
}

Eigen::VectorXd LogQuadraticSpline1DPDF::b_gradient(unsigned isegment) const
{
  if(isegment > nknot_-2)
  {
    std::ostringstream stream;
    stream << "LogQuadraticSpline1DPDF::b_gradient: isegment out of range: "
           << isegment << " > " << nknot_-2;
    throw(std::out_of_range(stream.str()));
  }
  return b_gradient_.col(isegment);
}


unsigned LogQuadraticSpline1DPDF::find_segment(double x) const
{
  unsigned ix = std::lower_bound(xknot_.data(), xknot_.data()+nknot_-1, x)-xknot_.data();
  return std::max(1U,ix)-1;
}

void LogQuadraticSpline1DPDF::set_cache()
{
  dy_ = yknot_.tail(nknot_-1)-yknot_.head(nknot_-1);

  if(p0_loc_ == ParamZeroLocation::RIGHT)
    set_spline_coeffs_right_to_left();
  else
    set_spline_coeffs_left_to_right();

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

      integral(a, b, c, xl, xr, I, dI_da, dI_db);

      if(!std::isfinite(I) or I==0)
      {
        std::cout << "Segment has zero norm: "
                  << isegment << ' ' << xl << ' ' << xr << ' '
                  <<  a << ' ' << b << ' ' << c << ' '
                  << I << ' ' << dI_da << ' ' << dI_db <<'\n';
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

void LogQuadraticSpline1DPDF::set_spline_coeffs_left_to_right()
{
  if(p0_type_ == ParamZeroType::SLOPE)
  {
    b_(0) = param0_;
    a_(0) = (dy_(0) - dx_(0) * b_(0))/SQR(dx_(0));

    b_gradient_(0,0) = 1.0;
    a_gradient_(0,0) = -1.0/dx_(0);
    a_gradient_(1,0) = -1.0/SQR(dx_(0));
    a_gradient_(2,0) = 1.0/SQR(dx_(0));
  }
  else
  {
    a_(0) = param0_;
    b_(0) = dy_(0)/dx_(0) - a_(0) * dx_(0);

    a_gradient_(0,0) = 1.0;
    b_gradient_(0,0) = -dx_(0);
    b_gradient_(1,0) = -1.0/dx_(0);
    b_gradient_(2,0) = 1.0/dx_(0);
  }

  for(unsigned isegment=1; isegment<nknot_-1; isegment++)
  {
    double dx2 = SQR(dx_(isegment));
    b_(isegment) = 2.0*dx_(isegment-1)*a_(isegment-1) + b_(isegment-1);
    a_(isegment) = (dy_(isegment) - dx_(isegment)*b_(isegment))/dx2;

    b_gradient_.col(isegment) =
        2.0*dx_(isegment-1)*a_gradient_.col(isegment-1)
        + b_gradient_.col(isegment-1);
    a_gradient_.col(isegment) =
        -b_gradient_.col(isegment)/dx_(isegment);
    a_gradient_(isegment+2, isegment) += 1.0/dx2;
    a_gradient_(isegment+1, isegment) -= 1.0/dx2;
  }
}

void LogQuadraticSpline1DPDF::set_spline_coeffs_right_to_left()
{
  if(p0_type_ == ParamZeroType::SLOPE)
  {
    unsigned iseg = nknot_-2;
    b_(iseg) = 2*dy_(iseg)/dx_(iseg) - param0_;
    a_(iseg) = (dy_(iseg) - dx_(iseg) * b_(iseg))/SQR(dx_(iseg));

    b_gradient_(0,iseg) = -1.0;
    b_gradient_(iseg+1,iseg) -= 2.0/dx_(iseg);
    b_gradient_(iseg+2,iseg) += 2.0/dx_(iseg);
    a_gradient_.col(iseg) = -b_gradient_.col(iseg)/dx_(iseg);
    a_gradient_(iseg+1,iseg) -= 1.0/SQR(dx_(iseg));
    a_gradient_(iseg+2,iseg) += 1.0/SQR(dx_(iseg));
  }
  else
  {
    unsigned iseg = nknot_-2;
    a_(iseg) = param0_;
    b_(iseg) = dy_(iseg)/dx_(iseg) - a_(iseg) * dx_(iseg);

    a_gradient_(0,iseg) = 1.0;
    b_gradient_(0,iseg) = -dx_(iseg);
    b_gradient_(iseg+1,iseg) = -1.0/dx_(iseg);
    b_gradient_(iseg+2,iseg) = 1.0/dx_(iseg);
  }

  for(unsigned iloop=1; iloop<nknot_-1; iloop++)
  {
    unsigned iseg = nknot_-2-iloop;
    double dx2 = SQR(dx_(iseg));
    b_(iseg) = 2*dy_(iseg)/dx_(iseg) - b_(iseg+1);
    a_(iseg) = (dy_(iseg) - dx_(iseg) * b_(iseg))/dx2;

    b_gradient_.col(iseg) = -b_gradient_.col(iseg+1);
    b_gradient_(iseg+1,iseg) -= 2.0/dx_(iseg);
    b_gradient_(iseg+2,iseg) += 2.0/dx_(iseg);
    a_gradient_.col(iseg) = -b_gradient_.col(iseg)/dx_(iseg);
    a_gradient_(iseg+1,iseg) -= 1.0/dx2;
    a_gradient_(iseg+2,iseg) += 1.0/dx2;
  }
}

void LogQuadraticSpline1DPDF::
integral(double a, double b, double c, double xl, double xr,
         double& I, double& dI_da, double& dI_db)
{
  if(a == 0)
  {
    if(b == 0)
    {
      const double F_exp = std::exp(c);
      I = F_exp*(xr-xl);
      dI_da = F_exp * (CUBE(xr) - CUBE(xl))/3;
      dI_db = F_exp * (SQR(xr) - SQR(xl))/2;
    }
    else
    {
      const double F_exp_r = std::exp(b*xr+c);
      const double F_exp_l = std::exp(b*xl+c);
      I = (F_exp_r - F_exp_l)/b;
      dI_da = ((b*xr*(b*xr-2)+2)*F_exp_r
               - (b*xl*(b*xl-2)+2)*F_exp_l)/CUBE(b);
      dI_db = ((b*xr-1)*F_exp_r - (b*xl-1)*F_exp_l)/SQR(b);
    }
  }
  else // a!=0
  {
    const double arg_exp = -SQR(b)/(4*a)+c;
    const double arg_exp_r = a*SQR(xr+b/(2*a));
    const double arg_exp_l = a*SQR(xl+b/(2*a));

    const double FF_exp_r = std::exp(arg_exp + arg_exp_r);
    const double FF_exp_l = std::exp(arg_exp + arg_exp_l);

    if(a < 0) // negative curvature - integral is erf
    {
      const double s = std::sqrt(-a);
      const double xxr = s*(xr+b/(2*a));
      const double xxl = s*(xl+b/(2*a));
      if(xxl>20)
        I = 1/(M_2_SQRTPI*s)*(std::exp(arg_exp + lerfc(xxl))
                              - std::exp(arg_exp + lerfc(xxr)));
      else if(xxr<-20)
        I = 1/(M_2_SQRTPI*s)*(std::exp(arg_exp + lerfc(-xxr))
                              - std::exp(arg_exp + lerfc(-xxl)));
      else if(xxl>5)
        I = 1/(M_2_SQRTPI*s)*std::exp(arg_exp)*(std::erfc(xxl)-std::erfc(xxr));
      else if(xxr<-5)
        I = 1/(M_2_SQRTPI*s)*std::exp(arg_exp)*(std::erfc(-xxr)-std::erfc(-xxl));
      else
        I = 1/(M_2_SQRTPI*s)*std::exp(arg_exp)*(std::erf(xxr)-std::erf(xxl));
      //std::cout << "DDD: " << xxl << ' ' << xxr << ' ' << erf(xxr)-erf(xxl) << ' ' << erfc(xxl)-erfc(xxr) << ' ' << erfc(-xxr)-erfc(-xxl) << ' ' << F_exp << ' ' << I << '\n';

    }
    else // positive curvature - integral is e^(x**2) dawson(x)
    {
      const double s = std::sqrt(a);
      I = 1/s*(FF_exp_r*dawson(s*(xr+b/(2*a))) -
               FF_exp_l*dawson(s*(xl+b/(2*a))));
    }

    dI_da = I*(SQR(b)/(4*SQR(a)) - 1/(2*a))
            + 1/(2*a)*(FF_exp_r*(xr-b/(2*a)) - FF_exp_l*(xl-b/(2*a)));
    dI_db = (-I*b + (FF_exp_r - FF_exp_l))/(2*a);
  }

#if 0
  std::cout << "       " << xl << ' ' << xr << ' '
            <<  a << ' ' <<  b << ' ' << c << ' '
            << I << ' ' << dI_da << ' ' << dI_db << '\n';
#endif
}
