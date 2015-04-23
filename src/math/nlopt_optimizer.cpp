/* 

   calin/math/nlopt_optimizer.cpp -- Stephen Fegan -- 2015-03-16

   Interface to MIT NLOpt optimizer suite

*/

#include <iostream>
#include <iomanip>

#include "math/nlopt_optimizer.hpp"
#include "math/hessian.hpp"

using namespace calin::math::optimizer;

NLOptOptimizer::
NLOptOptimizer(algorithm_type algorithm,
               function::MultiAxisFunction* fcn, bool adopt_fcn)
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
  auto axes = fcn_->domain_axes();
  nlopt::opt opt(algorithm_, axes.size());
  opt.set_min_objective(nlopt_callback, this);

#if 1
  std::vector<double> xlim_lo { limits_lo() };
  if(!xlim_lo.empty() &&
     std::count(xlim_lo.begin(), xlim_lo.end(), neg_inf)!=xlim_lo.size())
    opt.set_lower_bounds(xlim_lo);

  std::vector<double> xlim_hi { limits_hi() };
  if(!xlim_hi.empty() &&
     std::count(xlim_hi.begin(), xlim_hi.end(), pos_inf)!=xlim_hi.size())
    opt.set_upper_bounds(xlim_hi);
#endif
  
  std::vector<double> stepsize { initial_stepsize() };
  opt.set_initial_step(stepsize);

  if(abs_tolerance()>0)opt.set_ftol_abs(abs_tolerance());
  if(rel_tolerance()>0)opt.set_ftol_rel(rel_tolerance());
  if(max_iterations()>0)opt.set_maxeval(max_iterations());

  err_est_->reset(fcn_->domain_axes().size());
  iter_ = 0;

  std::vector<double> x { initial_values() };

  for(unsigned ipar=0; ipar<axes.size(); ipar++)
  {
    std::cout << axes[ipar].name << ' '
              << ((ipar<x.size())?x[ipar]:0.0) << ' '
              << ((ipar<stepsize.size())?stepsize[ipar]:0.0) << ' '
              << (ipar<xlim_lo.size()?xlim_lo[ipar]:0.0) << ' '
              << (ipar<xlim_hi.size()?xlim_hi[ipar]:0.0) << '\n';
  }

  opt.optimize(x, fopt);

  xopt.resize(x.size());
  xopt = Eigen::Map<Eigen::VectorXd>(x.data(), x.size());
  xopt_ = xopt;
}

ErrorMatrixStatus NLOptOptimizer::error_matrix_estimate(MatRef err_mat)
{
  return err_est_->error_matrix(err_mat);
}

ErrorMatrixStatus NLOptOptimizer::
calc_error_matrix_and_eigenvectors(MatRef err_mat,
                                   VecRef eigenvalues, MatRef eigenvectors)
{
  const unsigned npar { fcn_->num_domain_axes() };
  err_mat.resize(npar,npar);
  Eigen::VectorXd error_hint;
  if(error_matrix_estimate(err_mat) != ErrorMatrixStatus::UNAVAILABLE)
    error_hint = err_mat.diagonal().array().sqrt();
  Eigen::MatrixXd hessian(npar,npar);
  hessian::calculate_hessian(*fcn_, xopt_, hessian, error_hint);
  return hessian::hessian_to_error_matrix(*fcn_, hessian, err_mat,
                                          eigenvalues, eigenvectors);
}

double NLOptOptimizer::nlopt_callback(unsigned n, const double* x, double* grad,
                                      void* self)
{
  return static_cast<NLOptOptimizer*>(self)->eval_func(n,x,grad);
}

double NLOptOptimizer::eval_func(unsigned n, const double* x, double* grad)
{
  Eigen::VectorXd xvec = Eigen::Map<const Eigen::VectorXd>(x,n);

  double fcn_value;
  try
  {
    if(grad)
    {
      Eigen::VectorXd gvec(n);
      fcn_value = fcn_->value_and_gradient(xvec,gvec);
      Eigen::Map<Eigen::VectorXd>(grad,n) = gvec;
      err_est_->incorporate_func_gradient(xvec, fcn_value, gvec);
    }
    else
    {
      fcn_value = fcn_->value(xvec);
      err_est_->incorporate_func_value(xvec, fcn_value);
    }
  }
  catch(std::exception& x)
  {
    std::cout
        << "NLOptOptimizer::eval_func: exception in user function with x = "
        << xvec.transpose() << '\n'
        << "NLOptOptimizer::eval_func: exception is: " << x.what() << '\n';
    throw;
  }

  if(!isfinite(fcn_value))fcn_value = std::numeric_limits<double>::infinity();

  if(verbose_ != VerbosityLevel::SILENT)
  {
    std::cout << std::left << std::setw(4) << iter_+1 << ' '
              << std::fixed << std::setprecision(8) << fcn_value;
    for(unsigned ipar=0;ipar<n;ipar++)std::cout << ' ' << x[ipar];
    if(grad)
      for(unsigned ipar=0;ipar<n;ipar++)std::cout << ' ' << grad[ipar];
    std::cout << '\n';
  }
    
  iter_++;
  return fcn_value;
}
