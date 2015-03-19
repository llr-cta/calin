/* 

   calin/math/optimizer.cpp -- Stephen Fegan -- 2015-03-12

*/

#include <iostream>
#include <iomanip>

#include "math/optimizer.hpp"

using namespace calin::math::optimizer;

Optimizer::~Optimizer()
{
  if(my_fcn_)delete fcn_;
}

std::vector<double> Optimizer::initial_values() const
{
  // Set initial values from caller's values with fallback to function defaults
  auto axes = fcn_->domain_axes();
  std::vector<double> x0;
  if(x0_.size() == axes.size())x0 = x0_;
  else for(const auto& ipar : axes)x0.push_back(ipar.initial_value);
  return x0;
}

std::vector<double> Optimizer::initial_stepsize() const
{
  // Set the step size from caller's values with fallback to function defaults
  auto axes = fcn_->domain_axes();
  std::vector<double> xscale;
  if(xscale_.size() == axes.size())xscale = xscale_;
  else for(const auto& ipar : axes)xscale.push_back(ipar.scale);
  for(auto& x : xscale)x *= stepsize_scale_;
  return xscale;
}

std::vector<double> Optimizer::limits_lo() const
{
  // Set lower limits from caller's values with fallback to function defaults
  constexpr auto inf = std::numeric_limits<double>::infinity();
  auto axes = fcn_->domain_axes();
  std::vector<double> xlim_lo;
  if(xlim_lo_.size() == axes.size())xlim_lo = xlim_lo_;
  else for(const auto& ipar : axes) {
      xlim_lo.push_back(ipar.has_lo_bound?ipar.lo_bound:-inf); }
  if(std::all_of(xlim_lo.begin(),xlim_lo.end(),[](double x){return x==-inf;}))
    return std::vector<double>{};
  return xlim_lo;
}

std::vector<double> Optimizer::limits_hi() const
{
  // Set upper limits from caller's values with fallback to function defaults
  constexpr auto inf = std::numeric_limits<double>::infinity();
  auto axes = fcn_->domain_axes();
  std::vector<double> xlim_hi;
  if(xlim_hi_.size() == axes.size())xlim_hi = xlim_hi_;
  else for(const auto& ipar : axes) {
      xlim_hi.push_back(ipar.has_hi_bound?ipar.hi_bound:inf); }
  if(std::all_of(xlim_hi.begin(),xlim_hi.end(),[](double x){return x==inf;}))
    return std::vector<double>{};
  return xlim_hi;
}

// =============================================================================
//
// Error matrix estimator
//
// =============================================================================

ErrorMatrixEstimator::~ErrorMatrixEstimator()
{
  // nothing to see here
}

IdentityErrorMatrixEstimator::~IdentityErrorMatrixEstimator()
{
  // nothing to see here
}

void IdentityErrorMatrixEstimator::reset(unsigned npar)
{
  npar_ = npar;
}

void IdentityErrorMatrixEstimator::
invalid_func_value(function::ConstVecRef x)
{
  // nothing to see here
}

void IdentityErrorMatrixEstimator::
incorporate_func_value(function::ConstVecRef x, double f_val)
{
  // nothing to see here
}

void IdentityErrorMatrixEstimator::
incorporate_func_gradient(function::ConstVecRef x, double f_val,
                          function::ConstVecRef gradient)
{
  // nothing to see here
}

auto IdentityErrorMatrixEstimator::
error_matrix(Eigen::MatrixXd& err_mat) -> Status
{
  err_mat.resize(npar_,npar_);
  err_mat.setIdentity();
  return Status::UNAVAILABLE;
}

BFGSErrorMatrixEstimator::~BFGSErrorMatrixEstimator()
{
  // nothing to see here
}

void BFGSErrorMatrixEstimator::reset(unsigned npar)
{
  npar_ = npar;
  last_good_ = false;
  Bk_.resize(npar_,npar_);
  Bk_.setIdentity();
  xk_.resize(npar);
  gk_.resize(npar);
}

void BFGSErrorMatrixEstimator::invalid_func_value(function::ConstVecRef x)
{
  // Ignore last point for next update
  last_good_ = false;
}

void BFGSErrorMatrixEstimator::
incorporate_func_value(function::ConstVecRef x, double f_val)
{
  // BGFS requires derivative so ignore last point for next update
  last_good_ = false;
}

void BFGSErrorMatrixEstimator::
incorporate_func_gradient(function::ConstVecRef x, double f_val,
                          function::ConstVecRef gradient)
{
  // Form error matrix estimate using BFGS update method, see
  // http://en.wikipedia.org/wiki/Broyden-Fletcher-Goldfarb-Shanno_algorithm

  if(last_good_)
  {
    Eigen::VectorXd sk = x;
    sk -= xk_;
      
    Eigen::VectorXd yk = gradient;
    yk -= gk_;

    const double skyk = sk.dot(yk);
    sk /= skyk;

    std::cout << "B: " << gradient.data() << '\n';
    
#ifdef BFGS_COMPUTE_WITH_LOOPS
    Eigen::VectorXd bkyk(n);
    for(unsigned i=0;i<n;i++)
    {
      bkyk(i) = Bk_(i,i)*yk(i);
      for(unsigned j=i+1;j<n;j++)bkyk(i) += Bk_(i,j)*yk(j);
      for(unsigned j=0;j<i;j++)bkyk(i) += Bk_(j,i)*yk(j);
    }
    const double C1 = skyk + yk.dot(bkyk);
    for(unsigned i=0;i<n;i++)
      for(unsigned j=i;j<n;j++)
        Bk_(i,j) += C1*sk(i)*sk(j) - bkyk(i)*sk(j) - sk(i)*bkyk(j);
#else
    const Eigen::VectorXd bkyk = Bk_*yk;
    Bk_.noalias() += // What will Eigen make of this?
        (skyk + yk.dot(bkyk))*sk*sk.transpose()
        - bkyk*sk.transpose()
        - sk*bkyk.transpose();
#endif
    std::cout << std::scientific << std::setprecision(8) << Bk_ << "\n\n";
  }
      
  last_good_ = true;
  xk_ = x;
  gk_ = gradient;
}

auto BFGSErrorMatrixEstimator::error_matrix(Eigen::MatrixXd& err_mat) -> Status
{
  err_mat = 2.0*error_up_*Bk_;
#ifdef BFGS_COMPUTE_WITH_LOOPS
    for(unsigned i=0;i<n;i++)
      for(unsigned j=i+1;j<n;j++)
        err_mat(j,i) = err_mat(i,j);
#endif
  return Status::GOOD;
}


