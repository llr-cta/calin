/*

   calin/math/pdf_1d.hpp -- Stephen Fegan -- 2015-04-02

   Base classes for one-dimensional PDF functions. PDF functions are
   based on ParameterizableSingleAxisFunction, and are assumed to
   integrate out to 1.0. They also (optionally) provide analytic
   moments.

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

#pragma once

#include <string>
#include <vector>
#include <limits>

#include <Eigen/Core>

#include "calin_global_definitions.hpp"
#include "function.hpp"

namespace calin { namespace math { namespace pdf_1d {

#if 0
class Parameterizable1DPDF: public function::ParameterizableSingleAxisFunction
{
public:
  virtual ~Parameterizable1DPDF();

  // Reiterate functions from ParameterizableSingleAxisFunction

  unsigned num_parameters() override = 0;
  std::vector<function::ParameterAxis> parameters() override = 0;
  Eigen::VectorXd parameter_values() override = 0;
  void set_parameter_values(ConstVecRef values) override = 0;

  function::DomainAxis domain_axis() override = 0;

  bool can_calculate_gradient() override = 0;
  bool can_calculate_hessian() override = 0;
  bool can_calculate_parameter_gradient() override = 0;
  bool can_calculate_parameter_hessian() override = 0;

  double value_1d(double x) override = 0;
  double value_and_gradient_1d(double x,  double& dfdx) override = 0;
  double value_gradient_and_hessian_1d(double x, double& dfdx,
                            double& d2fdx2) override = 0;
  double value_and_parameter_gradient_1d(double x,
                                         VecRef gradient) override = 0;
  double value_parameter_gradient_and_hessian_1d(double x, VecRef gradient,
                                              MatRef hessian) override = 0;
};
#else
CALIN_TYPEALIAS(Parameterizable1DPDF, function::ParameterizableSingleAxisFunction);
#endif

// *****************************************************************************
//
// Miscellaneous PDFs
//
// *****************************************************************************

class GaussianPDF: public Parameterizable1DPDF
{
 public:
  GaussianPDF(const std::string& xunits = "x-value units"):
      Parameterizable1DPDF(), xunits_(xunits)
  { /* nothing to see here */ }

  virtual ~GaussianPDF();

  unsigned num_parameters() override;
  std::vector<function::ParameterAxis> parameters() override;
  Eigen::VectorXd parameter_values() override;
  void set_parameter_values(ConstVecRef values) override;
  function::DomainAxis domain_axis() override;

  bool can_calculate_gradient() override;
  bool can_calculate_hessian() override;
  bool can_calculate_parameter_gradient() override;
  bool can_calculate_parameter_hessian() override;

  double value_1d(double x) override;
  double value_and_gradient_1d(double x,  double& dfdx) override;
  double value_gradient_and_hessian_1d(double x, double& dfdx,
                            double& d2fdx2) override;
  double value_and_parameter_gradient_1d(double x,  VecRef gradient) override;
  double value_parameter_gradient_and_hessian_1d(double x, VecRef gradient,
                                                 MatRef hessian) override;

 protected:
  std::string xunits_;
  double x0_ = 0;
  double s_ = 1;
};

class BinnedGaussianPDF: public Parameterizable1DPDF
{
 public:
  BinnedGaussianPDF(double dx = 0, const std::string& xunits = "x-value units"):
      Parameterizable1DPDF(),
      dx_2_(0.5*dx), xunits_(xunits)
  { /* nothing to see here */ }

  virtual ~BinnedGaussianPDF();

  unsigned num_parameters() override;
  std::vector<function::ParameterAxis> parameters() override;
  Eigen::VectorXd parameter_values() override;
  void set_parameter_values(ConstVecRef values) override;
  function::DomainAxis domain_axis() override;

  bool can_calculate_gradient() override;
  bool can_calculate_hessian() override;
  bool can_calculate_parameter_gradient() override;
  bool can_calculate_parameter_hessian() override;

  double value_1d(double x) override;
  double value_and_gradient_1d(double x,  double& dfdx) override;
  double value_gradient_and_hessian_1d(double x, double& dfdx,
                            double& d2fdx2) override;
  double value_and_parameter_gradient_1d(double x,  VecRef gradient) override;
  double value_parameter_gradient_and_hessian_1d(double x, VecRef gradient,
                                                 MatRef hessian) override;

protected:
  double dx_2_ = 0;
  std::string xunits_;
  double x0_ = 0;
  double s_ = 1;
};

class LimitedGaussianPDF: public BinnedGaussianPDF
{
 public:
  constexpr static double inf = std::numeric_limits<double>::infinity();

  LimitedGaussianPDF(double xlo, double xhi, double dx = 0,
    const std::string& xunits = "x-value units"):
      BinnedGaussianPDF(dx, xunits), xlo_(xlo), xhi_(xhi),
      norm_gradient_(2), norm_hessian_(2,2) { set_cache(); }

  virtual ~LimitedGaussianPDF();

  function::DomainAxis domain_axis() override;

  void set_parameter_values(ConstVecRef values) override;

  double value_1d(double x) override;
  double value_and_gradient_1d(double x,  double& dfdx) override;
  double value_gradient_and_hessian_1d(double x, double& dfdx,
                                    double& d2fdx2) override;
  double value_and_parameter_gradient_1d(double x,  VecRef gradient) override;
  double value_parameter_gradient_and_hessian_1d(double x, VecRef gradient,
                                                 MatRef hessian) override;

protected:
  void set_cache();
  double xlo_;
  double xhi_;
  double norm_ = 0;
  Eigen::VectorXd norm_gradient_;
  Eigen::MatrixXd norm_hessian_;
};

class LimitedExponentialPDF: public Parameterizable1DPDF
{
 public:
  constexpr static double inf = std::numeric_limits<double>::infinity();

  LimitedExponentialPDF(double xlo=0.0, double xhi=inf, double dx = 0,
			const std::string& xunits = "x-value units"):
    Parameterizable1DPDF(), xunits_(xunits),
    xlo_(xlo), xhi_(xhi), dx_(dx) { set_cache(); }

  virtual ~LimitedExponentialPDF();

  void limit_scale(double scale_lo, double scale_hi) {
    limit_a_lo_ = scale_lo;
    limit_a_hi_ = scale_hi;
  }

  unsigned num_parameters() override;
  std::vector<function::ParameterAxis> parameters() override;
  Eigen::VectorXd parameter_values() override;
  void set_parameter_values(ConstVecRef values) override;
  function::DomainAxis domain_axis() override;

  bool can_calculate_gradient() override;
  bool can_calculate_hessian() override;
  bool can_calculate_parameter_gradient() override;
  bool can_calculate_parameter_hessian() override;

  double value_1d(double x) override;
  double value_and_gradient_1d(double x,  double& dfdx) override;
  double value_gradient_and_hessian_1d(double x, double& dfdx,
                            double& d2fdx2) override;
  double value_and_parameter_gradient_1d(double x,  VecRef gradient) override;
  double value_parameter_gradient_and_hessian_1d(double x, VecRef gradient,
                                                 MatRef hessian) override;

 protected:
  void set_cache();
  std::string xunits_;
  double limit_a_lo_    = -inf;
  double limit_a_hi_    = inf;
  double a_             = 1.0;
  double xlo_;
  double xhi_;
  double dx_;
  double norm_          = 1.0;
  double norm_gradient_ = 0.0;
  double norm_hessian_  = 0.0;
};

class TwoComponent1DPDF: public Parameterizable1DPDF
{
 public:
  TwoComponent1DPDF(Parameterizable1DPDF* pdf1, const std::string& cpt1_name,
                  Parameterizable1DPDF* pdf2, const std::string& cpt2_name,
                  bool adopt_pdf1 = false, bool adopt_pdf2 = false);
  virtual ~TwoComponent1DPDF();

  unsigned num_parameters() override;
  std::vector<function::ParameterAxis> parameters() override;
  Eigen::VectorXd parameter_values() override;
  void set_parameter_values(ConstVecRef values) override;

  function::DomainAxis domain_axis() override;

  bool can_calculate_gradient() override;
  bool can_calculate_hessian() override;
  bool can_calculate_parameter_gradient() override;
  bool can_calculate_parameter_hessian() override;

  double value_1d(double x) override;
  double value_and_gradient_1d(double x,  double& dfdx) override;
  double value_gradient_and_hessian_1d(double x, double& dfdx,
                            double& d2fdx2) override;
  double value_and_parameter_gradient_1d(double x,  VecRef gradient) override;
  double value_parameter_gradient_and_hessian_1d(double x, VecRef gradient,
                                                 MatRef hessian) override;

 protected:
  double prob_cpt1_;
  Parameterizable1DPDF* pdf1_;
  bool adopt_pdf1_;
  std::string cpt1_name_;
  Parameterizable1DPDF* pdf2_;
  bool adopt_pdf2_;
  std::string cpt2_name_;
};

class TwoComponent1DConstraintPDF: public Parameterizable1DPDF
{
 public:
  TwoComponent1DConstraintPDF(Parameterizable1DPDF* pdfTest, bool fast_mode = false);
  virtual ~TwoComponent1DConstraintPDF();
  unsigned num_parameters() override;
  std::vector<function::ParameterAxis> parameters() override;
  Eigen::VectorXd parameter_values() override;
  void set_parameter_values(ConstVecRef values) override;

  function::DomainAxis domain_axis() override;

  bool can_calculate_gradient() override;
  bool can_calculate_hessian() override;
  bool can_calculate_parameter_gradient() override;
  bool can_calculate_parameter_hessian() override;

  double value_1d(double x) override;
  double value_and_gradient_1d(double x,  double& dfdx) override;
  double value_gradient_and_hessian_1d(double x, double& dfdx,
                            double& d2fdx2) override;
  double value_and_parameter_gradient_1d(double x,  VecRef gradient) override;
  double value_parameter_gradient_and_hessian_1d(double x, VecRef gradient,
                                                 MatRef hessian) override;
  Eigen::MatrixXd Jacobian(double x, double mu1,double mu2,double pp,double res,double n);


 protected:
  
  double prob_cpt1_;
  double pp_;
  double mu1_;
  double mu2_;
  double res_;
  double n_;
  Parameterizable1DPDF* pdfTest_;
  bool fast_mode_;
  
  
};

} } } // namespace calin::math::pdf_1d
