/* 

   calin/math/nlopt_optimizer.cpp -- Stephen Fegan -- 2015-03-16

   Interface to MIT NLOpt optimizer suite

*/

#include "math/nlopt_optimizer.hpp"

using namespace calin::math::optimizer;

NLOptOptimizer::
NLOptOptimizer(algorithm_type algorithm, MultiAxisFunction* fcn, bool adopt_fcn)
    : Optimizer(fcn, adopt_fcn), algorithm_(algorithm)
{
  // nothing to see here
}

NLOptOptimizer::~NLOptOptimizer()
{
  // nothing to see here
}

bool requires_gradient(algorithm_type a)
{
  switch(a)
  {
    case GN_DIRECT:
    case GN_DIRECT_L:
    case GN_DIRECT_L_RAND:
    case GN_DIRECT_NOSCAL:
    case GN_DIRECT_L_NOSCAL:
    case GN_DIRECT_L_RAND_NOSCAL:
    case GN_ORIG_DIRECT:
    case GN_ORIG_DIRECT_L:
    case LN_PRAXIS:
    case GN_CRS2_LM:
    case GN_MLSL:
    case GN_MLSL_LDS:
    case LN_COBYLA:
    case LN_NEWUOA:
    case LN_NEWUOA_BOUND:
    case LN_NELDERMEAD:
    case LN_SBPLX:
    case LN_AUGLAG:
    case LN_AUGLAG_EQ:
    case LN_BOBYQA:
    case GN_ISRES:
    case GN_ESCH:
      return false;

    case GD_STOGO:
    case GD_STOGO_RAND:
    case LD_LBFGS_NOCEDAL:
    case LD_LBFGS:
    case LD_VAR1:
    case LD_VAR2:
    case LD_TNEWTON:
    case LD_TNEWTON_RESTART:
    case LD_TNEWTON_PRECOND:
    case LD_TNEWTON_PRECOND_RESTART:
    case GD_MLSL:
    case GD_MLSL_LDS:
    case LD_MMA:
    case LD_AUGLAG:
    case LD_AUGLAG_EQ:
    case AUGLAG:
    case AUGLAG_EQ:
    case G_MLSL:
    case G_MLSL_LDS:
    case LD_SLSQP:
    case LD_CCSAQ:
      return true;
  }
}

bool requires_hessian(algorithm_type a)
{
  return false;
}

bool can_estimate_error(algorithm_type a)
{
  return false;
}

bool can_use_gradient(algorithm_type a)
{
  return requires_gradient(a);
}

bool can_use_hessian(algorithm_type a)
{
  return false; // (a==LD_CCSAQ);
}
    
bool NLOptOptimizer::requires_gradient()
{
  return requires_gradient(algorithm_);
}

bool NLOptOptimizer::requires_hessian()
{
  return requires_hessian(algorithm_);
}

bool NLOptOptimizer::can_estimate_error()
{
  return can_estimate_error(algorithm_);
}

bool NLOptOptimizer::can_use_gradient()
{
  return can_use_gradient(algorithm_);
}

bool NLOptOptimizer::can_use_hessian()
{
  return can_use_hessian(algorithm_);
}


bool NLOptOptimizer::minimize(std::vector<double>& xopt)
{

}

bool NLOptOptimizer::error_matrix_estimate(Eigen::MatrixXd& err_mat)
{

}

bool NLOptOptimizer::calc_error_matrix(Eigen::MatrixXd& err_mat)
{

}

