/* 

   calin/math/function.cpp -- Stephen Fegan -- 2015-02-24

*/

#include "math/function.hpp"

// namespace calin { namespace math { namespace function {

using namespace calin::math::function;

Parameterizable::~Parameterizable()
{
  // nothing to see here
}

MultiAxisFunction::~MultiAxisFunction()
{
  // nothing to see here
}

SingleAxisFunction::~SingleAxisFunction()
{
  // nothing to see here
}

std::vector<DomainAxis> SingleAxisFunction::domain_axes()
{
  return { domain_axis() };
}

double SingleAxisFunction::value(ConstVecRef x)
{
  return value(x(0));
}

double SingleAxisFunction::value_and_gradient(ConstVecRef x, VecRef gradient) 
{
  gradient.resize(1);
  return value_and_deriv(x(0),gradient(0));
}

double SingleAxisFunction::
value_gradient_and_hessian(ConstVecRef x, VecRef gradient, MatRef hessian) 
{
  gradient.resize(1);
  hessian.resize(1,1);
  return value_gradient_and_hessian(x(0),gradient(0),hessian(0,0));
}

#if 0
DoubleAxisFunction::~DoubleAxisFunction()
{
  // nothing to see here
}
#endif

// } } } // namespace calin::math::function
