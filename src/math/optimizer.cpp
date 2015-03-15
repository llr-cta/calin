/* 

   calin/math/optimizer.cpp -- Stephen Fegan -- 2015-03-12

*/

#include "math/optimizer.hpp"

using namespace calin::math::optimizer;

Optimizer::~Optimizer()
{
  if(my_fcn_)delete fcn_;
}
