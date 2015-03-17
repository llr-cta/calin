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

bool NLOptOptimizer::requires_box_constraints(algorithm_type a)
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
    case nlopt::GN_CRS2_LM:
    case nlopt::GN_MLSL:
    case nlopt::GN_MLSL_LDS:
    case nlopt::GN_ISRES:
    case nlopt::GN_ESCH:
    case nlopt::GD_STOGO:
    case nlopt::GD_STOGO_RAND:
    case nlopt::GD_MLSL:
    case nlopt::GD_MLSL_LDS:
    case nlopt::G_MLSL:
    case nlopt::G_MLSL_LDS:
      return true;

    case nlopt::LN_PRAXIS:
    case nlopt::LN_COBYLA:
    case nlopt::LN_NEWUOA:
    case nlopt::LN_NEWUOA_BOUND:
    case nlopt::LN_NELDERMEAD:
    case nlopt::LN_SBPLX:
    case nlopt::LN_AUGLAG:
    case nlopt::LN_AUGLAG_EQ:
    case nlopt::LN_BOBYQA:
    case nlopt::LD_LBFGS_NOCEDAL:
    case nlopt::LD_LBFGS:
    case nlopt::LD_VAR1:
    case nlopt::LD_VAR2:
    case nlopt::LD_TNEWTON:
    case nlopt::LD_TNEWTON_RESTART:
    case nlopt::LD_TNEWTON_PRECOND:
    case nlopt::LD_TNEWTON_PRECOND_RESTART:
    case nlopt::LD_MMA:
    case nlopt::LD_AUGLAG:
    case nlopt::LD_AUGLAG_EQ:
    case nlopt::AUGLAG:
    case nlopt::AUGLAG_EQ:
    case nlopt::LD_SLSQP:
    case nlopt::LD_CCSAQ:
      return false;

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

bool NLOptOptimizer::can_impose_box_constraints(algorithm_type a)
{
  return true;
}

bool NLOptOptimizer::requires_gradient()
{
  return requires_gradient(algorithm_);
}

bool NLOptOptimizer::requires_hessian()
{
  return requires_hessian(algorithm_);
}

bool NLOptOptimizer::requires_box_constraints()
{
  return requires_box_constraints(algorithm_);
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

bool NLOptOptimizer::can_impose_box_constraints()
{
  return can_impose_box_constraints(algorithm_);
}

bool NLOptOptimizer::minimize(std::vector<double>& xopt, double& fopt)
{
  constexpr auto inf = std::numeric_limits<double>::infinity();
  auto axes = fcn_->domain_axes();
  nlopt::opt opt(algorithm_, axes.size());
  opt.set_min_objective(nlopt_callback, this);

  // Set lower limits from caller's values with fallback to function defaults
  std::vector<double> xlim_lo;
  if(xlim_lo_.size() == axes.size())xlim_lo = xlim_lo_;
  else for(const auto& ipar : axes) {
      xlim_lo.push_back(ipar.has_lo_bound?ipar.lo_bound:-inf); }
  if(std::any_of(xlim_lo.begin(),xlim_lo.end(),[](double x){return x!=-inf;}))
    opt.set_lower_bounds(xlim_lo);

  // Set upper limits from caller's values with fallback to function defaults
  std::vector<double> xlim_hi;
  if(xlim_hi_.size() == axes.size())xlim_hi = xlim_hi_;
  else for(const auto& ipar : axes) {
      xlim_hi.push_back(ipar.has_hi_bound?ipar.hi_bound:inf); }
  if(std::any_of(xlim_hi.begin(),xlim_hi.end(),[](double x){return x!=inf;}))
    opt.set_upper_bounds(xlim_hi);

  // Set the step size from caller's values with fallback to function defaults
  std::vector<double> xscale;
  if(xscale_.size() == axes.size())xscale = xscale_;
  else for(const auto& ipar : axes)xscale.push_back(ipar.scale);
  opt.set_initial_step(xscale);

  // Set initial values from caller's values with fallback to function defaults
  xopt.clear();
  if(x0_.size() == axes.size())xopt = x0_;
  else for(const auto& ipar : axes)xopt.push_back(ipar.initial_value);
  
  opt.set_ftol_abs(0.001);
  opt.optimize(xopt, fopt);
}

bool NLOptOptimizer::error_matrix_estimate(Eigen::MatrixXd& err_mat)
{

}

bool NLOptOptimizer::calc_error_matrix(Eigen::MatrixXd& err_mat)
{

}

double NLOptOptimizer::nlopt_callback(unsigned n, const double* x, double* grad,
                                      void* self)
{
  return static_cast<NLOptOptimizer*>(self)->eval_func(x,grad);
}

double NLOptOptimizer::eval_func(const double* x, double* grad)
{
  double fcn_eval { grad?fcn_->value_and_gradient(x,grad):fcn_->value(x) };
  return fcn_eval;
}
