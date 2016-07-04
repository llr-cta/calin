/*

   calin/math/pdf_1d_two_component.cpp -- Stephen Fegan -- 2015-04-02

   Base classes for one-dimensional PDF functions. PDF functions are
   based on ParameterizableSingleAxisFunction, and are assumed to
   integrate out to 1.0. They also (optionally) provide analytic
   moments.

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, Ecole polytechnique, CNRS/IN2P3, Universite Paris-Saclay

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
#include <iomanip>
#include <stdexcept>

#include <math/pdf_1d.hpp>

using namespace calin::math;
using namespace calin::math::pdf_1d;
using function::assign_parameters;

// ============================================================================
//
// TwoComponentPDF
//
// ============================================================================

TwoComponentPDF::
TwoComponentPDF(Parameterizable1DPDF* pdf1, const std::string& cpt1_name,
                Parameterizable1DPDF* pdf2, const std::string& cpt2_name,
                bool adopt_pdf1, bool adopt_pdf2, double error_up):
    Parameterizable1DPDF(),
    pdf1_(pdf1), adopt_pdf1_(adopt_pdf1), cpt1_name_(cpt1_name),
    pdf2_(pdf2), adopt_pdf2_(adopt_pdf2), cpt2_name_(cpt2_name),
    error_up_(error_up)
{
  // nothing to see here
}

TwoComponentPDF::~TwoComponentPDF()
{
  // nothing to see here
}

unsigned TwoComponentPDF::num_parameters()
{
  return 1 + pdf1_->num_parameters() + pdf2_->num_parameters();
}

std::vector<function::ParameterAxis> TwoComponentPDF::parameters()
{
  std::vector<function::ParameterAxis> pvec {
    { cpt1_name_+std::string("_probability"), "1", true, 0, true, 1 } };
  std::vector<function::ParameterAxis> pdf1_pvec = pdf1_->parameters();
  for(auto& p : pdf1_pvec)p.name = cpt1_name_ + std::string(".") + p.name;
  pvec.insert(pvec.end(), pdf1_pvec.begin(), pdf1_pvec.end());
  std::vector<function::ParameterAxis> pdf2_pvec = pdf2_->parameters();
  for(auto& p : pdf2_pvec)p.name = cpt2_name_ + std::string(".") + p.name;
  pvec.insert(pvec.end(), pdf2_pvec.begin(), pdf2_pvec.end());
  return pvec;
}

Eigen::VectorXd TwoComponentPDF::parameter_values()
{
  Eigen::VectorXd param(num_parameters());
  param[0] = prob_cpt1_;
  unsigned num_cpt1_params = pdf1_->num_parameters();
  param.segment(1,num_cpt1_params)
      = pdf1_->parameter_values();
  param.segment(num_cpt1_params+1,pdf2_->num_parameters())
      = pdf2_->parameter_values();
  return param;
}

void TwoComponentPDF::set_parameter_values(ConstVecRef values)
{
  assign_parameters(values, prob_cpt1_);
  unsigned num_cpt1_params = pdf1_->num_parameters();
  pdf1_->set_parameter_values(values.segment(1,num_cpt1_params));
  pdf2_->set_parameter_values(values.segment(num_cpt1_params+1,
                                             pdf2_->num_parameters()));
}

function::DomainAxis TwoComponentPDF::domain_axis()
{
  return pdf1_->domain_axis();
}

bool TwoComponentPDF::can_calculate_gradient()
{
  return pdf1_->can_calculate_gradient() and pdf2_->can_calculate_gradient();
}

bool TwoComponentPDF::can_calculate_hessian()
{
  return pdf1_->can_calculate_hessian() and pdf2_->can_calculate_hessian();
}

bool TwoComponentPDF::can_calculate_parameter_gradient()
{
  return pdf1_->can_calculate_parameter_gradient() and
      pdf2_->can_calculate_parameter_gradient();
}

bool TwoComponentPDF::can_calculate_parameter_hessian()
{
  return pdf1_->can_calculate_parameter_hessian() and
      pdf2_->can_calculate_parameter_hessian();
}

double TwoComponentPDF::value_1d(double x)
{
  //std::cout << pdf1_->value_1d(x) << ' ' << pdf2_->value_1d(x) << '\n';
  return prob_cpt1_*pdf1_->value_1d(x) + (1.0-prob_cpt1_)*pdf2_->value_1d(x);
}

double TwoComponentPDF::value_and_gradient_1d(double x,  double& dfdx)
{
  double dfdx1;
  double dfdx2;
  double val = prob_cpt1_*pdf1_->value_and_gradient_1d(x,dfdx1) +
               (1.0-prob_cpt1_)*pdf2_->value_and_gradient_1d(x,dfdx2);
  dfdx = prob_cpt1_*dfdx1 + (1.0-prob_cpt1_)*dfdx2;
  return val;
}

double TwoComponentPDF::value_gradient_and_hessian_1d(double x, double& dfdx,
                                           double& d2fdx2)
{
  double dfdx1;
  double dfdx2;
  double d2fdx21;
  double d2fdx22;
  double val = prob_cpt1_*pdf1_->value_gradient_and_hessian_1d(x,dfdx1,d2fdx21)+
        (1.0-prob_cpt1_)*pdf2_->value_gradient_and_hessian_1d(x,dfdx2,d2fdx22);
  dfdx = prob_cpt1_*dfdx1 + (1.0-prob_cpt1_)*dfdx2;
  d2fdx2 = prob_cpt1_*d2fdx21 + (1.0-prob_cpt1_)*d2fdx22;
  return val;
}

double
TwoComponentPDF::value_and_parameter_gradient_1d(double x,  VecRef gradient)
{
  const double omp = 1.0-prob_cpt1_;
  const unsigned npar1 = pdf1_->num_parameters();
  const unsigned npar2 = pdf2_->num_parameters();
  gradient.resize(1+npar1+npar2);

#ifndef CALIN_USE_EIGEN_REF
  Eigen::VectorXd grad1(npar1);
  Eigen::VectorXd grad2(npar2);
#else
  auto grad1 = gradient.segment(1,npar1);
  auto grad2 = gradient.segment(1+npar1,npar2);
#endif

  double val1 = pdf1_->value_and_parameter_gradient_1d(x, grad1);
  double val2 = pdf2_->value_and_parameter_gradient_1d(x, grad2);

  grad1 *= prob_cpt1_;
  grad2 *= omp;
  gradient[0] = val1 - val2;
#ifndef CALIN_USE_EIGEN_REF
  gradient.segment(1,npar1) = grad1;
  gradient.segment(1+npar1,npar2) = grad2;
#endif
  return prob_cpt1_*val1 + omp*val2;
}

double TwoComponentPDF::
value_parameter_gradient_and_hessian_1d(double x, VecRef gradient,
                                        MatRef hessian)
{
  const double omp = 1.0-prob_cpt1_;
  const unsigned npar1 = pdf1_->num_parameters();
  const unsigned npar2 = pdf2_->num_parameters();
  gradient.resize(1+npar1+npar2);
  hessian.resize(1+npar1+npar2,1+npar1+npar2);
  hessian.setZero();

#ifndef CALIN_USE_EIGEN_REF
  Eigen::VectorXd grad1(npar1);
  Eigen::VectorXd grad2(npar2);
  Eigen::MatrixXd hess1(npar1, npar1);
  Eigen::MatrixXd hess2(npar2, npar2);
#else
  auto grad1 = gradient.segment(1,npar1);
  auto grad2 = gradient.segment(1+npar1,npar2);
  auto hess1 = hessian.block(1,1,npar1,npar1);
  auto hess2 = hessian.block(1+npar1,1+npar1,npar2,npar2);
#endif

  double val1 = pdf1_->value_parameter_gradient_and_hessian_1d(x, grad1, hess1);
  double val2 = pdf2_->value_parameter_gradient_and_hessian_1d(x, grad2, hess2);

  hess1 *= prob_cpt1_;
  hess2 *= omp;
#ifndef CALIN_USE_EIGEN_REF
  hessian.block(1,1,npar1,npar1) = hess1;
  hessian.block(1+npar1,1+npar1,npar2,npar2) = hess2;
#endif
  hessian.block(1,0,npar1,1) = grad1;
  hessian.block(1+npar1,0,npar2,1) = -grad2;
  hessian.block(0,1,1,npar1) = grad1.transpose();
  hessian.block(0,1+npar1,1,npar2) = -grad2.transpose();

  grad1 *= prob_cpt1_;
  grad2 *= omp;
  gradient[0] = val1 - val2;
#ifndef CALIN_USE_EIGEN_REF
  gradient.segment(1,npar1) = grad1;
  gradient.segment(1+npar1,npar2) = grad2;
#endif

  return prob_cpt1_*val1 + omp*val2;
}

double TwoComponentPDF::error_up()
{
  return error_up_;
}

bool TwoComponentPDF::can_calculate_mean_and_variance()
{
  return false;
}

void TwoComponentPDF::mean_and_variance(double& mean, double& var)
{
  assert(0);
}
