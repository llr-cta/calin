/* 

   calin/math/optimizer.cpp -- Stephen Fegan -- 2015-03-12

*/

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
