/* 

   calin/math/optimizer.cpp -- Stephen Fegan -- 2015-03-12

*/

#include <iostream>
#include <iomanip>

#include "math/optimizer.hpp"

using namespace calin::math::optimizer;

const double Optimizer::inf;
const double Optimizer::pos_inf;
const double Optimizer::neg_inf;

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
  auto axes = fcn_->domain_axes();
  std::vector<double> xlim(axes.size());
  std::transform(axes.begin(), axes.end(), xlim.begin(),
                 [](const decltype(axes)::value_type& ipar)
                 { return ipar.has_lo_bound?ipar.lo_bound:neg_inf;});
  if(xlim_lo_.size() <= axes.size())
    std::transform(xlim_lo_.begin(), xlim_lo_.end(), xlim.begin(),
                   xlim.begin(), [](double ipar1, double ipar2) {
                     return std::max(ipar1, ipar2); });
  return xlim;
}

std::vector<double> Optimizer::limits_hi() const
{
  // Set upper limits from caller's values with fallback to function defaults
  auto axes = fcn_->domain_axes();
  std::vector<double> xlim(axes.size());
  std::transform(axes.begin(), axes.end(), xlim.begin(),
                 [](const decltype(axes)::value_type& ipar)
                 { return ipar.has_hi_bound?ipar.hi_bound:pos_inf;});
  if(xlim_hi_.size() <= axes.size())
    std::transform(xlim_hi_.begin(), xlim_hi_.end(), xlim.begin(),
                   xlim.begin(), [](double ipar1, double ipar2) {
                     return std::min(ipar1, ipar2); });
  return xlim;
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
invalid_func_value(ConstVecRef x)
{
  // nothing to see here
}

void IdentityErrorMatrixEstimator::
incorporate_func_value(ConstVecRef x, double f_val)
{
  // nothing to see here
}

void IdentityErrorMatrixEstimator::
incorporate_func_gradient(ConstVecRef x, double f_val,
                          ConstVecRef gradient)
{
  // nothing to see here
}

ErrorMatrixStatus IdentityErrorMatrixEstimator::
error_matrix(MatRef err_mat)
{
  err_mat.resize(npar_,npar_);
  err_mat.setIdentity();
  return ErrorMatrixStatus::UNAVAILABLE;
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

void BFGSErrorMatrixEstimator::invalid_func_value(ConstVecRef x)
{
  // Ignore last point for next update
  last_good_ = false;
}

void BFGSErrorMatrixEstimator::
incorporate_func_value(ConstVecRef x, double f_val)
{
  // BGFS requires derivative so ignore last point for next update
  last_good_ = false;
}

void BFGSErrorMatrixEstimator::
incorporate_func_gradient(ConstVecRef x, double f_val,
                          ConstVecRef gradient)
{
  // Form error matrix estimate using BFGS update method, see
  // http://en.wikipedia.org/wiki/Broyden-Fletcher-Goldfarb-Shanno_algorithm

  if(!isfinite(f_val) or !std::all_of(gradient.data(), gradient.data()+npar_, 
                                      [](const double& x){return isfinite(x);}))
  {
    last_good_ = false;
    return;
  }
  
  if(last_good_)
  {
    Eigen::VectorXd sk = x;
    sk -= xk_;
      
    Eigen::VectorXd yk = gradient;
    yk -= gk_;

    const double skyk = sk.dot(yk);
    if(skyk == 0)goto skip_rank1_update; // don't like goto? i don't care :-)
    sk /= skyk;

#ifdef BFGS_COMPUTE_WITH_LOOPS // todo: test if faster when npar->inf
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
    //std::cout << std::scientific << std::setprecision(8) << Bk_ << "\n\n";
  }
skip_rank1_update:
  
  last_good_ = true;
  xk_ = x;
  gk_ = gradient;
}

ErrorMatrixStatus BFGSErrorMatrixEstimator::error_matrix(MatRef err_mat)
{
  err_mat = 2.0*error_up_*Bk_;
#ifdef BFGS_COMPUTE_WITH_LOOPS
    for(unsigned i=0;i<n;i++)
      for(unsigned j=i+1;j<n;j++)
        err_mat(j,i) = err_mat(i,j);
#endif
  return ErrorMatrixStatus::GOOD;
}


