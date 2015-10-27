/* 

   calin/math/pdf_1d.hpp -- Stephen Fegan -- 2015-04-02

   Base classes for one-dimensional PDF functions. PDF functions are
   based on ParameterizableSingleAxisFunction, and are assumed to
   integrate out to 1.0. They also (optionally) provide analytic
   moments.

*/

#pragma once

#include <string>
#include <vector>
#include <limits>

#include <Eigen/Core>

#include "package_wide_definitions.hpp"
#include "function.hpp"

namespace calin { namespace math { namespace pdf_1d {

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

  double error_up() override = 0;

  // Moments

  virtual bool can_calculate_mean_and_variance() = 0;
  virtual void mean_and_variance(double& mean, double& var) = 0;
};

// *****************************************************************************
//
// Miscellaneous PDFs
//
// *****************************************************************************

class GaussianPDF: public Parameterizable1DPDF
{
 public:
  GaussianPDF(const std::string& xunits = "x-value units", double error_up = 0.5):
      Parameterizable1DPDF(), xunits_(xunits), error_up_(error_up)
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

  double error_up() override;

  bool can_calculate_mean_and_variance() override;
  void mean_and_variance(double& mean, double& var) override;

 protected:
  std::string xunits_;
  double error_up_ = 0.5;
  double x0_ = 0;
  double s_ = 1;
};

class LimitedGaussianPDF: public GaussianPDF
{
 public:
#ifndef SWIG
  constexpr static double inf = std::numeric_limits<double>::infinity();
#endif
  
  LimitedGaussianPDF(double xlo, double xhi, const std::string& xunits = "x-value units",
		     double error_up = 0.5):
      GaussianPDF(xunits, error_up), xlo_(xlo), xhi_(xhi),
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

  bool can_calculate_mean_and_variance() override;
  void mean_and_variance(double& mean, double& var) override;
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
#ifndef SWIG
  constexpr static double inf = std::numeric_limits<double>::infinity();
#endif
  
  LimitedExponentialPDF(double xlo, double xhi, double dx = 0,
			const std::string& xunits = "x-value units",
                        double error_up = 0.5):
    Parameterizable1DPDF(), xunits_(xunits), error_up_(error_up),
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

  double error_up() override;

  bool can_calculate_mean_and_variance() override;
  void mean_and_variance(double& mean, double& var) override;  
 protected:
  void set_cache();
  std::string xunits_;
  double limit_a_lo_    = -inf;
  double limit_a_hi_    = inf;
  double error_up_      = 0.5;
  double a_             = 1.0;
  double xlo_;
  double xhi_;
  double dx_;
  double norm_          = 1.0;
  double norm_gradient_ = 0.0;
  double norm_hessian_  = 0.0;
};

class TwoComponentPDF: public Parameterizable1DPDF
{
 public:
  TwoComponentPDF(Parameterizable1DPDF* pdf1, const std::string& cpt1_name,
                  Parameterizable1DPDF* pdf2, const std::string& cpt2_name,
                  bool adopt_pdf1 = false, bool adopt_pdf2 = false,
                  double error_up = 0.5);
  virtual ~TwoComponentPDF();
  
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

  double error_up() override;

  // Moments

  bool can_calculate_mean_and_variance() override;
  void mean_and_variance(double& mean, double& var) override;

 protected:
  double prob_cpt1_;
  Parameterizable1DPDF* pdf1_;
  Parameterizable1DPDF* pdf2_;
  bool adopt_pdf1_;
  bool adopt_pdf2_;
  std::string cpt1_name_;
  std::string cpt2_name_;
  double error_up_;
};

} } } // namespace calin::math::pdf_1d
