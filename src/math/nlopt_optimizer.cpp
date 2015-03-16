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

bool NLOptOptimizer::requires_gradient(algorithm_type a)
{
  switch(a)
  {
    case nlopt::GN_DIRECT:
    case nlopt::GN_DIRECT_L:
    case nlopt::GN_DIRECT_L_RAND:
    case nlopt::GN_DIRECT_NOSCAL:
    case nlopt::GN_DIRECT_L_NOSCAL:
    case nlopt::GN_DIRECT_L_RAND_NOSCAL:
    case nlopt::GN_ORIG_DIRECT:
    case nlopt::GN_ORIG_DIRECT_L:
    case nlopt::LN_PRAXIS:
    case nlopt::GN_CRS2_LM:
    case nlopt::GN_MLSL:
    case nlopt::GN_MLSL_LDS:
    case nlopt::LN_COBYLA:
    case nlopt::LN_NEWUOA:
    case nlopt::LN_NEWUOA_BOUND:
    case nlopt::LN_NELDERMEAD:
    case nlopt::LN_SBPLX:
    case nlopt::LN_AUGLAG:
    case nlopt::LN_AUGLAG_EQ:
    case nlopt::LN_BOBYQA:
    case nlopt::GN_ISRES:
    case nlopt::GN_ESCH:
      return false;

    case nlopt::GD_STOGO:
    case nlopt::GD_STOGO_RAND:
    case nlopt::LD_LBFGS_NOCEDAL:
    case nlopt::LD_LBFGS:
    case nlopt::LD_VAR1:
    case nlopt::LD_VAR2:
    case nlopt::LD_TNEWTON:
    case nlopt::LD_TNEWTON_RESTART:
    case nlopt::LD_TNEWTON_PRECOND:
    case nlopt::LD_TNEWTON_PRECOND_RESTART:
    case nlopt::GD_MLSL:
    case nlopt::GD_MLSL_LDS:
    case nlopt::LD_MMA:
    case nlopt::LD_AUGLAG:
    case nlopt::LD_AUGLAG_EQ:
    case nlopt::AUGLAG:
    case nlopt::AUGLAG_EQ:
    case nlopt::G_MLSL:
    case nlopt::G_MLSL_LDS:
    case nlopt::LD_SLSQP:
    case nlopt::LD_CCSAQ:
      return true;

    case nlopt::NUM_ALGORITHMS:
    default:
      assert(0);
  }

  assert(0);
}

bool NLOptOptimizer::requires_hessian(algorithm_type a)
{
  return false;
}

bool NLOptOptimizer::can_estimate_error(algorithm_type a)
{
  return false;
}

bool NLOptOptimizer::can_use_gradient(algorithm_type a)
{
  return requires_gradient(a);
}

bool NLOptOptimizer::can_use_hessian(algorithm_type a)
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

