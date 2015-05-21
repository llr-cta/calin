/* 

   calin/math/cminpack_optimizer.cpp -- Stephen Fegan -- 2015-05-21

   Interface to ANL CMinpack optimizer suite

*/

#include <iostream>
#include <iomanip>
#include <cstring>

#include "cminpack/cminpack.h"

#include "io/log.hpp"
#include "math/cminpack_optimizer.hpp"
#include "math/hessian.hpp"

using namespace calin::math::optimizer;
using namespace calin::io::log;

class enum CMinpackMode { AUTO, REQUIRE_HESSIAN, DONT_USE_HESSIAN };

CMinpackOptimizer::
CMinpackOptimizer(CMinpackMode mode,
                  function::MultiAxisFunction* fcn, bool adopt_fcn)
{

}

CMinpackOptimizer::~CMinpackOptimizer()
{

}

bool CMinpackOptimizer::requires_gradient()
{
  return true;
}

bool CMinpackOptimizer::requires_hessian()
{
  return mode_==CMinpackMode::REQUIRE_HESSIAN;
}

bool CMinpackOptimizer::requires_box_constraints()
{
  return false;
}

bool CMinpackOptimizer::can_estimate_error()
{
  return true;
}

bool CMinpackOptimizer::can_use_gradient()
{
  return true;
}

bool CMinpackOptimizer::can_use_hessian()
{
  return mode_!=CMinpackMode::DONT_USE_HESSIAN
}

bool CMinpackOptimizer::can_impose_box_constraints()
{
  return false;
}
