/* 

   calin/math/nlopt_optimizer.cpp -- Stephen Fegan -- 2015-03-16

   Interface to MIT NLOpt optimizer suite

*/

#include <iostream>
#include <iomanip>

#include "math/nlopt_optimizer.hpp"

using namespace calin::math::optimizer;

NLOptOptimizer::
NLOptOptimizer(algorithm_type algorithm, MultiAxisFunction* fcn, bool adopt_fcn)
    : Optimizer(fcn, adopt_fcn), algorithm_(algorithm)
{
  if(can_use_gradient(algorithm))
    err_est_.reset(new BFGSErrorMatrixEstimator(fcn->error_up()));
  else
    err_est_.reset(new IdentityErrorMatrixEstimator(fcn->error_up()));
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

bool NLOptOptimizer::minimize(VecRef xopt, double& fopt)
{
  constexpr auto inf = std::numeric_limits<double>::infinity();
  auto axes = fcn_->domain_axes();
  nlopt::opt opt(algorithm_, axes.size());
  opt.set_min_objective(nlopt_callback, this);

  std::vector<double> xlim_lo { limits_lo() };
  if(!xlim_lo.empty())opt.set_lower_bounds(xlim_lo);
  std::vector<double> xlim_hi { limits_hi() };
  if(!xlim_hi.empty())opt.set_upper_bounds(xlim_hi);
  opt.set_initial_step(initial_stepsize());

  if(abs_tolerance()>0)opt.set_ftol_abs(abs_tolerance());
  if(rel_tolerance()>0)opt.set_ftol_rel(rel_tolerance());
  if(max_iterations()>0)opt.set_maxeval(max_iterations());

  err_est_->reset(fcn_->domain_axes().size());
  iter_ = 0;

  std::vector<double> x { initial_values() };
  opt.optimize(x, fopt);

  xopt.resize(x.size());
  xopt = Eigen::Map<Eigen::VectorXd>(x.data(), x.size());
}

bool NLOptOptimizer::error_matrix_estimate(MatRef err_mat)
{
  return err_est_->error_matrix(err_mat) != ErrorMatrixStatus::UNAVAILABLE;
}

bool NLOptOptimizer::calc_error_matrix(MatRef err_mat)
{

}

double NLOptOptimizer::nlopt_callback(unsigned n, const double* x, double* grad,
                                      void* self)
{
  return static_cast<NLOptOptimizer*>(self)->eval_func(n,x,grad);
}

double NLOptOptimizer::eval_func(unsigned n, const double* x, double* grad)
{
  Eigen::Map<const Eigen::VectorXd> xvec(x,n);
  Eigen::Map<Eigen::VectorXd> gvec(grad,n);
  double fcn_value { grad?fcn_->value_and_gradient(xvec,gvec):fcn_->value(xvec) };
  if(verbose_ != VerbosityLevel::SILENT)
  {
    std::cout << std::left << std::setw(4) << iter_+1 << ' '
              << std::fixed << std::setprecision(4) << fcn_value;
    for(unsigned ipar=0;ipar<n;ipar++)std::cout << ' ' << x[ipar];
    std::cout << '\n';
  }
  
  if(!isfinite(fcn_value))
  {
    fcn_value = std::numeric_limits<double>::infinity();
    err_est_->invalid_func_value(xvec);
  }
  else if(grad)
  {
    err_est_->incorporate_func_gradient(xvec, fcn_value, gvec);
  }
  else
  {
    err_est_->incorporate_func_value(xvec, fcn_value);
  }
  
  iter_++;
  return fcn_value;
}
