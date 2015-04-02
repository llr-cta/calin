/* 

   calin/math/function.hpp -- Stephen Fegan -- 2015-02-24

   Base classes for functions and general parameterizable objects that
   can be used with optimizers, root finders, the MCMC algorithm etc.

*/

#pragma once

#include <string>
#include <vector>
#include <limits>

#include <Eigen/Core>

namespace calin { namespace math { namespace function {

#if 1
using VecRef = Eigen::Ref<Eigen::VectorXd>;
using ConstVecRef = const Eigen::Ref<const Eigen::VectorXd>&;
using MatRef = Eigen::Ref<Eigen::MatrixXd>;
using ConstMatRef = const Eigen::Ref<const Eigen::MatrixXd>&;
#else
using VecRef = Eigen::VectorXd&;
using MatRef = Eigen::MatrixXd&;
using ConstVecRef = const Eigen::Ref<const Eigen::VectorXd>&;
using ConstMatRef = const Eigen::Ref<const Eigen::MatrixXd>&;
#endif

struct ParameterAxis
{
  ParameterAxis(const std::string& _name = { },
                const std::string& _units = { },
                bool _has_lo_bound = false,
                double _lo_bound = -std::numeric_limits<double>::infinity(),
                bool _has_hi_bound  = false,
                double _hi_bound = std::numeric_limits<double>::infinity(),
                double _scale = 1.0, double _initial_value = 0):
      name(_name), units(_units),
      has_lo_bound(_has_lo_bound), lo_bound(_lo_bound),
      has_hi_bound(_has_hi_bound), hi_bound(_hi_bound),
      scale(_scale), initial_value(_initial_value)
  { /* nothing to see here */ } 

  std::string name;
  std::string units;
  bool has_lo_bound;
  double lo_bound;
  bool has_hi_bound;
  double hi_bound;
  double scale;
  double initial_value;
};

using DomainAxis = ParameterAxis;

// Variadic template.. what fun!

template<typename T>
void assign_parameters(const double* values, T& x)
{
  x = *values;
}

template<typename T>
void assign_parameters(ConstVecRef values, T& x)
{
  x = values[0];
}

template<typename T, typename... Params>
void assign_parameters(const double* values, T& x, Params & ... params)
{
  x = *values;
  assign_parameters(++values, params...);
}

template<typename T, typename... Params>
void assign_parameters(ConstVecRef values, T& x, Params & ... params)
{
  const double* values_data = values.data();
  x = *values_data;
  assign_parameters(++values_data, params...);
}

class Parameterizable
{
 public:
  virtual ~Parameterizable();
  virtual unsigned num_parameters() = 0;
  virtual std::vector<ParameterAxis> parameters() = 0;
  virtual Eigen::VectorXd parameter_values() = 0;
  virtual void set_parameter_values(ConstVecRef values) = 0;
  virtual bool can_calculate_parameter_gradient() = 0;
  virtual bool can_calculate_parameter_hessian() = 0;
};

class MultiAxisFunction
{
 public:
  virtual ~MultiAxisFunction();
  virtual unsigned num_domain_axes() = 0;
  virtual std::vector<DomainAxis> domain_axes() = 0;
  virtual double value(ConstVecRef x) = 0;
  virtual bool can_calculate_gradient() = 0;
  virtual double value_and_gradient(ConstVecRef x, VecRef gradient) = 0;
  virtual bool can_calculate_hessian() = 0;
  virtual double value_gradient_and_hessian(ConstVecRef x, VecRef gradient,
                                            MatRef hessian) = 0;
  virtual double error_up() = 0;  
};

class SingleAxisFunction: virtual public MultiAxisFunction
{
 public:
  virtual ~SingleAxisFunction();
  virtual DomainAxis domain_axis() = 0;
  virtual double value(double x) = 0;
  virtual double value_and_gradient(double x,  double& dfdx) = 0;
  virtual double value_gradient_and_hessian(double x, double& dfdx,
                                          double& d2fdx2) = 0;

  // Members from MultiAxisFunction that we override
  unsigned num_domain_axes() override;
  std::vector<DomainAxis> domain_axes() override;
  double value(ConstVecRef x) override;
  double value_and_gradient(ConstVecRef x, VecRef gradient) override;
  double value_gradient_and_hessian(ConstVecRef x, VecRef gradient,
                                    MatRef hessian) override;
};

class ParameterizableMultiAxisFunction: public Parameterizable,
                                        virtual public MultiAxisFunction
{
 public:
  virtual ~ParameterizableMultiAxisFunction();
  virtual double value_and_parameter_gradient(ConstVecRef x,
                                         VecRef gradient) = 0;
  virtual double value_parameter_gradient_and_hessian(ConstVecRef x,
                                         VecRef gradient, MatRef hessian) = 0;
};

class ParameterizableSingleAxisFunction:
      public ParameterizableMultiAxisFunction,
      public SingleAxisFunction
{
 public:
  virtual ~ParameterizableSingleAxisFunction();
  virtual double value_and_parameter_gradient(double x,  VecRef gradient) = 0;
  virtual double value_parameter_gradient_and_hessian(double x,
                                         VecRef gradient, MatRef hessian) = 0;

  // Members from ParameterizableMultiAxisFunction that we override
  double value_and_parameter_gradient(ConstVecRef x,
                                      VecRef gradient) override;
  double value_parameter_gradient_and_hessian(ConstVecRef x,
                                      VecRef gradient, MatRef hessian) override;
};

class PMAFReverser: public ParameterizableMultiAxisFunction
{
 public:
  PMAFReverser(ParameterizableMultiAxisFunction* fcn_deligate,
               bool adopt_fcn_deligate = false, double error_up = 0.5);
  virtual ~PMAFReverser();

  // Parameterizable interface
  unsigned num_parameters() override;
  std::vector<ParameterAxis> parameters() override;
  Eigen::VectorXd parameter_values() override;
  void set_parameter_values(ConstVecRef values) override;
  bool can_calculate_parameter_gradient() override;
  bool can_calculate_parameter_hessian() override;

  // MultiAxisFunction interface
  unsigned num_domain_axes() override;
  std::vector<DomainAxis> domain_axes() override;
  double value(ConstVecRef x) override;
  bool can_calculate_gradient() override;
  double value_and_gradient(ConstVecRef x, VecRef gradient) override;
  bool can_calculate_hessian() override;
  double value_gradient_and_hessian(ConstVecRef x, VecRef gradient,
                                    MatRef hessian) override;
  double error_up() override;

  // ParameterizableMultiAxisFunction interface
  double value_and_parameter_gradient(ConstVecRef x,
                                      VecRef gradient) override;
  double value_parameter_gradient_and_hessian(ConstVecRef x,
                                      VecRef gradient, MatRef hessian) override;

 protected:
  ParameterizableMultiAxisFunction* fcn_deligate_;
  bool adopt_fcn_deligate_ = false;
  double error_up_ = 0.5;
  Eigen::VectorXd x_;
};

template<typename ParamType>
class MultiParameterSet: virtual public Parameterizable
{
 public:
  virtual ~MultiParameterSet() { /* nothing to see here */ }
  virtual std::vector<ParameterAxis> parameters() override;
  virtual Eigen::VectorXd parameter_values() override;
  virtual void set_parameter_values(ConstVecRef values) override;
 private:
  unsigned find_parameter_set(unsigned iparam);
  std::vector<std::pair<unsigned, ParamType*> > parameter_sets_;
};

bool gradient_check(MultiAxisFunction& fcn, ConstVecRef x, VecRef good,
                    double eps_factor = 10.0);
bool gradient_check(MultiAxisFunction& fcn, ConstVecRef x, ConstVecRef dx,
                    VecRef good);
bool hessian_check(MultiAxisFunction& fcn, ConstVecRef x, ConstVecRef dx,
                   MatRef good);

template<typename ParameterizableType> using ValGetter =
    std::function<double(ParameterizableType*)>;
template<typename ParameterizableType> using GradGetter =
    std::function<double(ParameterizableType*, VecRef)>;
template<typename ParameterizableType> using HessGetter =
    std::function<double(ParameterizableType*, VecRef, MatRef)>;
    
template<typename ParameterizableType>
class ParameterizableToMultiAxisFunctionAdapter: public MultiAxisFunction
{
 public:
  using PTValGetter = ValGetter<ParameterizableType>;
  using PTGradGetter = GradGetter<ParameterizableType>;
  using PTHessGetter = HessGetter<ParameterizableType>;

  ParameterizableToMultiAxisFunctionAdapter(ParameterizableType* par,
                                            PTValGetter val_getter,
                                            PTGradGetter grad_getter,
                                            PTHessGetter hess_getter,
                                            double error_up = 0.5):
      MultiAxisFunction(), par_(par), val_getter_(val_getter),
      grad_getter_(grad_getter), error_up_(error_up) { }
  ~ParameterizableToMultiAxisFunctionAdapter() { }
  unsigned num_domain_axes() override { return par_->num_parameters(); }
  std::vector<DomainAxis> domain_axes() override { return par_->parameters(); }
  double value(ConstVecRef x) {
    par_->set_parameter_values(x);
    return val_getter_(par_); }
  bool can_calculate_gradient() { 
    return par_->can_calculate_parameter_gradient(); }
  double value_and_gradient(ConstVecRef x, VecRef gradient) {
    par_->set_parameter_values(x);
    return grad_getter_(par_, gradient); }
  bool can_calculate_hessian() {
    return par_->can_calculate_parameter_hessian(); }
  double value_gradient_and_hessian(ConstVecRef x, VecRef gradient,
                                    MatRef hessian) {
    par_->set_parameter_values(x);
    return hess_getter_(par_, gradient, hessian); }
  double error_up() { return error_up_; }
 protected:
  ParameterizableType* par_;
  PTValGetter val_getter_;
  PTGradGetter grad_getter_;
  PTHessGetter hess_getter_;
  double error_up_;
};

template<typename ParameterizableType>
bool gradient_check_par(ParameterizableType& par_fcn,
                        ConstVecRef p, ConstVecRef dp, VecRef good,
                        ValGetter<ParameterizableType> val_getter,
                        GradGetter<ParameterizableType> grad_getter,
                        HessGetter<ParameterizableType> hess_getter,
                        double error_up = 0.5)
{
  function::ParameterizableToMultiAxisFunctionAdapter<ParameterizableType>
      fcn(&par_fcn, val_getter, grad_getter, hess_getter, error_up);
  return function::gradient_check(fcn, p, dp, good);
}

// *****************************************************************************
//
// Miscellaneous functions
//
// *****************************************************************************

class GaussianPDF: public ParameterizableSingleAxisFunction
{
 public:
  GaussianPDF(double error_up = 0.5):
      ParameterizableSingleAxisFunction(), error_up_(error_up)
  { /* nothing to see here */ }
  
  virtual ~GaussianPDF();

  unsigned num_parameters() override;
  std::vector<ParameterAxis> parameters() override;
  Eigen::VectorXd parameter_values() override;
  void set_parameter_values(ConstVecRef values) override;
  DomainAxis domain_axis() override;

  bool can_calculate_gradient() override;
  bool can_calculate_hessian() override;
  bool can_calculate_parameter_gradient() override;
  bool can_calculate_parameter_hessian() override;

  double value(double x) override;
  double value_and_gradient(double x,  double& dfdx) override;
  double value_gradient_and_hessian(double x, double& dfdx,
                                    double& d2fdx2) override;
  double value_and_parameter_gradient(double x,  VecRef gradient) override;
  double value_parameter_gradient_and_hessian(double x, VecRef gradient,
                                              MatRef hessian) override;

  double error_up() override;
 protected:
  double error_up_ = 0.5;
  double x0_ = 0;
  double s_ = 1;
};

class LimitedGaussianPDF: public GaussianPDF
{
 public:
  constexpr static double inf = std::numeric_limits<double>::infinity();

  LimitedGaussianPDF(double xlo, double xhi, double error_up = 0.5):
      GaussianPDF(error_up), xlo_(xlo), xhi_(xhi),
      norm_gradient_(2), norm_hessian_(2,2) { set_cache(); }
  
  virtual ~LimitedGaussianPDF();

  void set_parameter_values(ConstVecRef values) override;

  double value(double x) override;
  double value_and_gradient(double x,  double& dfdx) override;
  double value_gradient_and_hessian(double x, double& dfdx,
                                    double& d2fdx2) override;
  double value_and_parameter_gradient(double x,  VecRef gradient) override;
  double value_parameter_gradient_and_hessian(double x, VecRef gradient,
                                              MatRef hessian) override;
 protected:
  void set_cache();
  double xlo_;
  double xhi_;
  double norm_ = 0;
  Eigen::VectorXd norm_gradient_;
  Eigen::MatrixXd norm_hessian_;
};

class LimitedExponentialPDF: public ParameterizableSingleAxisFunction
{
 public:
  constexpr static double inf = std::numeric_limits<double>::infinity();

  LimitedExponentialPDF(double xlo=0.0, double xhi=inf, double error_up = 0.5):
      error_up_(error_up), xlo_(xlo), xhi_(xhi) { set_cache(); }

  virtual ~LimitedExponentialPDF();

  unsigned num_parameters() override;
  std::vector<ParameterAxis> parameters() override;
  Eigen::VectorXd parameter_values() override;
  void set_parameter_values(ConstVecRef values) override;
  DomainAxis domain_axis() override;

  bool can_calculate_gradient() override;
  bool can_calculate_hessian() override;
  bool can_calculate_parameter_gradient() override;
  bool can_calculate_parameter_hessian() override;

  double value(double x) override;
  double value_and_gradient(double x,  double& dfdx) override;
  double value_gradient_and_hessian(double x, double& dfdx,
                                    double& d2fdx2) override;
  double value_and_parameter_gradient(double x,  VecRef gradient) override;
  double value_parameter_gradient_and_hessian(double x, VecRef gradient,
                                              MatRef hessian) override;

  double error_up() override;
 protected:
  void set_cache();
  double error_up_      = 0.5;
  double a_             = 1.0;
  double xlo_;
  double xhi_;
  double norm_          = 1.0;
  double norm_gradient_ = 0.0;
  double norm_hessian_  = 0.0;
};

  
} // namespace function

using function::ParameterAxis;
using function::DomainAxis;
using function::Parameterizable;
using function::MultiAxisFunction;
using function::SingleAxisFunction;
using function::ParameterizableMultiAxisFunction;
using function::ParameterizableSingleAxisFunction;

using function::GaussianPDF;
using function::LimitedGaussianPDF;
using function::LimitedExponentialPDF;

} } // namespace calin::math
