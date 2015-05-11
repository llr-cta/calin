/* 

   calin/math/nlopt_optimizer.cpp -- Stephen Fegan -- 2015-03-16

   Interface to MIT NLOpt optimizer suite

*/

#include <iostream>
#include <iomanip>
#include <cstring>

#include "nlopt/nlopt.h"

#include "math/nlopt_optimizer.hpp"
#include "math/hessian.hpp"

using namespace calin::math::optimizer;

namespace {

bool does_algorithm_match(const std::string& name, const std::string& algo_name)
{
  if(strcasecmp(name.c_str(), algo_name.c_str()) == 0)return true;
  if(strncasecmp(algo_name.c_str(), "NLOPT_", 6)==0 and
     strcasecmp(name.c_str(), algo_name.c_str()+6) == 0)return true;
  if(strncasecmp(algo_name.c_str(), "NLOPP::", 7)==0 and
     strcasecmp(name.c_str(), algo_name.c_str()+7) == 0)return true;
  return false;
}

#define MATCH_ALGORITHM(x) if(does_algorithm_match(name, #x))return x

nlopt_algorithm name_to_algorithm(const std::string& name)
{
  MATCH_ALGORITHM(NLOPT_GN_DIRECT);
  MATCH_ALGORITHM(NLOPT_GN_DIRECT_L);
  MATCH_ALGORITHM(NLOPT_GN_DIRECT_L_RAND);
  MATCH_ALGORITHM(NLOPT_GN_DIRECT_NOSCAL);
  MATCH_ALGORITHM(NLOPT_GN_DIRECT_L_NOSCAL);
  MATCH_ALGORITHM(NLOPT_GN_DIRECT_L_RAND_NOSCAL);

  MATCH_ALGORITHM(NLOPT_GN_ORIG_DIRECT);
  MATCH_ALGORITHM(NLOPT_GN_ORIG_DIRECT_L);

  MATCH_ALGORITHM(NLOPT_GD_STOGO);
  MATCH_ALGORITHM(NLOPT_GD_STOGO_RAND);

  MATCH_ALGORITHM(NLOPT_LD_LBFGS_NOCEDAL);

  MATCH_ALGORITHM(NLOPT_LD_LBFGS);

  MATCH_ALGORITHM(NLOPT_LN_PRAXIS);
  MATCH_ALGORITHM(NLOPT_LD_VAR1);
  MATCH_ALGORITHM(NLOPT_LD_VAR2);

  MATCH_ALGORITHM(NLOPT_LD_TNEWTON);
  MATCH_ALGORITHM(NLOPT_LD_TNEWTON_RESTART);
  MATCH_ALGORITHM(NLOPT_LD_TNEWTON_PRECOND);
  MATCH_ALGORITHM(NLOPT_LD_TNEWTON_PRECOND_RESTART);

  MATCH_ALGORITHM(NLOPT_GN_CRS2_LM);

  MATCH_ALGORITHM(NLOPT_GN_MLSL);
  MATCH_ALGORITHM(NLOPT_GD_MLSL);
  MATCH_ALGORITHM(NLOPT_GN_MLSL_LDS);
  MATCH_ALGORITHM(NLOPT_GD_MLSL_LDS);

  MATCH_ALGORITHM(NLOPT_LD_MMA);

  MATCH_ALGORITHM(NLOPT_LN_COBYLA);

  MATCH_ALGORITHM(NLOPT_LN_NEWUOA);
  MATCH_ALGORITHM(NLOPT_LN_NEWUOA_BOUND);
  MATCH_ALGORITHM(NLOPT_LN_NELDERMEAD);
  MATCH_ALGORITHM(NLOPT_LN_SBPLX);

  MATCH_ALGORITHM(NLOPT_LN_AUGLAG);
  MATCH_ALGORITHM(NLOPT_LD_AUGLAG);
  MATCH_ALGORITHM(NLOPT_LN_AUGLAG_EQ);
  MATCH_ALGORITHM(NLOPT_LD_AUGLAG_EQ);

  MATCH_ALGORITHM(NLOPT_LN_BOBYQA);

  MATCH_ALGORITHM(NLOPT_GN_ISRES);

  /* new variants that require local_optimizer to be set,
     not with older constants for backwards compatibility */
  MATCH_ALGORITHM(NLOPT_AUGLAG);
  MATCH_ALGORITHM(NLOPT_AUGLAG_EQ);
  MATCH_ALGORITHM(NLOPT_G_MLSL);
  MATCH_ALGORITHM(NLOPT_G_MLSL_LDS);

  MATCH_ALGORITHM(NLOPT_LD_SLSQP);
  MATCH_ALGORITHM(NLOPT_LD_CCSAQ);

  MATCH_ALGORITHM(NLOPT_GN_ESCH);

  return NLOPT_NUM_ALGORITHMS;
}

}

NLOptOptimizer::
NLOptOptimizer(const std::string& algorithm_name,
               function::MultiAxisFunction* fcn, bool adopt_fcn)
    : Optimizer(fcn, adopt_fcn), algorithm_name_(algorithm_name)
{
  if(can_use_gradient(algorithm_name))
    err_est_.reset(new BFGSErrorMatrixEstimator(fcn->error_up()));
  else
    err_est_.reset(new IdentityErrorMatrixEstimator(fcn->error_up()));
}

NLOptOptimizer::~NLOptOptimizer()
{
  // nothing to see here
}

bool NLOptOptimizer::is_valid_algorithm(const std::string& algorithm_name)
{
  return name_to_algorithm(algorithm_name) != NLOPT_NUM_ALGORITHMS;
}

bool NLOptOptimizer::requires_gradient(const std::string& algorithm_name)
{
  switch(name_to_algorithm(algorithm_name))
  {
    case NLOPT_GN_DIRECT:
    case NLOPT_GN_DIRECT_L:
    case NLOPT_GN_DIRECT_L_RAND:
    case NLOPT_GN_DIRECT_NOSCAL:
    case NLOPT_GN_DIRECT_L_NOSCAL:
    case NLOPT_GN_DIRECT_L_RAND_NOSCAL:
    case NLOPT_GN_ORIG_DIRECT:
    case NLOPT_GN_ORIG_DIRECT_L:
    case NLOPT_LN_PRAXIS:
    case NLOPT_GN_CRS2_LM:
    case NLOPT_GN_MLSL:
    case NLOPT_GN_MLSL_LDS:
    case NLOPT_LN_COBYLA:
    case NLOPT_LN_NEWUOA:
    case NLOPT_LN_NEWUOA_BOUND:
    case NLOPT_LN_NELDERMEAD:
    case NLOPT_LN_SBPLX:
    case NLOPT_LN_AUGLAG:
    case NLOPT_LN_AUGLAG_EQ:
    case NLOPT_LN_BOBYQA:
    case NLOPT_GN_ISRES:
    case NLOPT_GN_ESCH:
      return false;

    case NLOPT_GD_STOGO:
    case NLOPT_GD_STOGO_RAND:
    case NLOPT_LD_LBFGS_NOCEDAL:
    case NLOPT_LD_LBFGS:
    case NLOPT_LD_VAR1:
    case NLOPT_LD_VAR2:
    case NLOPT_LD_TNEWTON:
    case NLOPT_LD_TNEWTON_RESTART:
    case NLOPT_LD_TNEWTON_PRECOND:
    case NLOPT_LD_TNEWTON_PRECOND_RESTART:
    case NLOPT_GD_MLSL:
    case NLOPT_GD_MLSL_LDS:
    case NLOPT_LD_MMA:
    case NLOPT_LD_AUGLAG:
    case NLOPT_LD_AUGLAG_EQ:
    case NLOPT_AUGLAG:
    case NLOPT_AUGLAG_EQ:
    case NLOPT_G_MLSL:
    case NLOPT_G_MLSL_LDS:
    case NLOPT_LD_SLSQP:
    case NLOPT_LD_CCSAQ:
      return true;

    case NLOPT_NUM_ALGORITHMS:
    default:
      assert(0);
  }

  assert(0);
}

bool NLOptOptimizer::requires_box_constraints(const std::string& algorithm_name)
{
  switch(name_to_algorithm(algorithm_name))
  {
    case NLOPT_GN_DIRECT:
    case NLOPT_GN_DIRECT_L:
    case NLOPT_GN_DIRECT_L_RAND:
    case NLOPT_GN_DIRECT_NOSCAL:
    case NLOPT_GN_DIRECT_L_NOSCAL:
    case NLOPT_GN_DIRECT_L_RAND_NOSCAL:
    case NLOPT_GN_ORIG_DIRECT:
    case NLOPT_GN_ORIG_DIRECT_L:
    case NLOPT_GN_CRS2_LM:
    case NLOPT_GN_MLSL:
    case NLOPT_GN_MLSL_LDS:
    case NLOPT_GN_ISRES:
    case NLOPT_GN_ESCH:
    case NLOPT_GD_STOGO:
    case NLOPT_GD_STOGO_RAND:
    case NLOPT_GD_MLSL:
    case NLOPT_GD_MLSL_LDS:
    case NLOPT_G_MLSL:
    case NLOPT_G_MLSL_LDS:
      return true;

    case NLOPT_LN_PRAXIS:
    case NLOPT_LN_COBYLA:
    case NLOPT_LN_NEWUOA:
    case NLOPT_LN_NEWUOA_BOUND:
    case NLOPT_LN_NELDERMEAD:
    case NLOPT_LN_SBPLX:
    case NLOPT_LN_AUGLAG:
    case NLOPT_LN_AUGLAG_EQ:
    case NLOPT_LN_BOBYQA:
    case NLOPT_LD_LBFGS_NOCEDAL:
    case NLOPT_LD_LBFGS:
    case NLOPT_LD_VAR1:
    case NLOPT_LD_VAR2:
    case NLOPT_LD_TNEWTON:
    case NLOPT_LD_TNEWTON_RESTART:
    case NLOPT_LD_TNEWTON_PRECOND:
    case NLOPT_LD_TNEWTON_PRECOND_RESTART:
    case NLOPT_LD_MMA:
    case NLOPT_LD_AUGLAG:
    case NLOPT_LD_AUGLAG_EQ:
    case NLOPT_AUGLAG:
    case NLOPT_AUGLAG_EQ:
    case NLOPT_LD_SLSQP:
    case NLOPT_LD_CCSAQ:
      return false;

    case NLOPT_NUM_ALGORITHMS:
    default:
      assert(0);
  }

  assert(0);
}


bool NLOptOptimizer::requires_hessian(const std::string& algorithm_name)
{
  return false;
}

bool NLOptOptimizer::can_estimate_error(const std::string& algorithm_name)
{
  return false;
}

bool NLOptOptimizer::can_use_gradient(const std::string& algorithm_name)
{
  return requires_gradient(algorithm_name);
}

bool NLOptOptimizer::can_use_hessian(const std::string& algorithm_name)
{
  return false; // (a==LD_CCSAQ);
}

bool NLOptOptimizer::can_impose_box_constraints(const std::string& algorithm_name)
{
  return true;
}

bool NLOptOptimizer::requires_gradient()
{
  return requires_gradient(algorithm_name_);
}

bool NLOptOptimizer::requires_hessian()
{
  return requires_hessian(algorithm_name_);
}

bool NLOptOptimizer::requires_box_constraints()
{
  return requires_box_constraints(algorithm_name_);
}

bool NLOptOptimizer::can_estimate_error()
{
  return can_estimate_error(algorithm_name_);
}

bool NLOptOptimizer::can_use_gradient()
{
  return can_use_gradient(algorithm_name_);
}

bool NLOptOptimizer::can_use_hessian()
{
  return can_use_hessian(algorithm_name_);
}

bool NLOptOptimizer::can_impose_box_constraints()
{
  return can_impose_box_constraints(algorithm_name_);
}

OptimizationStatus NLOptOptimizer::minimize(VecRef xopt, double& fopt)
{
  auto axes = fcn_->domain_axes();
  unsigned naxes { fcn_->num_domain_axes() };

  std::unique_ptr<nlopt_opt_s, void(*)(nlopt_opt)>
      opt(nlopt_create(name_to_algorithm(algorithm_name_), naxes),
          nlopt_destroy);
  nlopt_set_min_objective(opt.get(), nlopt_callback, this);

#if 1
  std::vector<double> xlim_lo { limits_lo() };
  if(!xlim_lo.empty() &&
     std::count(xlim_lo.begin(), xlim_lo.end(), neg_inf)!=xlim_lo.size())
    nlopt_set_lower_bounds(opt.get(), &xlim_lo.front());

  std::vector<double> xlim_hi { limits_hi() };
  if(!xlim_hi.empty() &&
     std::count(xlim_hi.begin(), xlim_hi.end(), pos_inf)!=xlim_hi.size())
    nlopt_set_upper_bounds(opt.get(), &xlim_hi.front());
#endif
  
  std::vector<double> stepsize { initial_stepsize() };
  nlopt_set_initial_step(opt.get(), &stepsize.front());
  
  if(abs_tolerance()>0)nlopt_set_ftol_abs(opt.get(), abs_tolerance());
  if(rel_tolerance()>0)nlopt_set_ftol_rel(opt.get(), rel_tolerance());
  if(max_iterations()>0)nlopt_set_maxeval(opt.get(), max_iterations());
  if(max_walltime()>0)nlopt_set_maxtime(opt.get(), max_walltime());

  std::vector<double> x { initial_values() };

  err_est_->reset(naxes);
  iterations_ = 0;
  fbest_ = inf;
  xbest_.resize(x.size());
  xbest_ = Eigen::Map<Eigen::VectorXd>(x.data(), x.size());
  
  for(unsigned ipar=0; ipar<axes.size(); ipar++)
  {
    std::cout << axes[ipar].name << ' '
              << ((ipar<x.size())?x[ipar]:0.0) << ' '
              << ((ipar<stepsize.size())?stepsize[ipar]:0.0) << ' '
              << (ipar<xlim_lo.size()?xlim_lo[ipar]:0.0) << ' '
              << (ipar<xlim_hi.size()?xlim_hi[ipar]:0.0) << '\n';
  }

  nlopt_result nlopt_status = nlopt_optimize(opt.get(), &x.front(), &fopt);

  OptimizationStatus return_status;
  switch(nlopt_status)
  {
    case NLOPT_SUCCESS:
    case NLOPT_FTOL_REACHED:
      return_status = OptimizationStatus::TOLERANCE_REACHED;
      break;
    case NLOPT_MAXEVAL_REACHED:
      return_status = OptimizationStatus::STOPPED_AT_MAXCALLS;
      break;
    case NLOPT_MAXTIME_REACHED:
      return_status = OptimizationStatus::STOPPED_AT_MAXTIME;
      break;
    case NLOPT_ROUNDOFF_LIMITED:
      return_status = OptimizationStatus::LIMITED_BY_PRECISION;
      break;
    case NLOPT_FAILURE:
    case NLOPT_INVALID_ARGS:
    case NLOPT_OUT_OF_MEMORY:
      return_status = OptimizationStatus::OPTIMIZER_FAILURE;
      break;
    case NLOPT_STOPVAL_REACHED:
    case NLOPT_XTOL_REACHED:
    case NLOPT_FORCED_STOP:
      assert(0);
  }

  if(return_status != OptimizationStatus::OPTIMIZER_FAILURE)
  {
    fbest_ = fopt;
    xbest_ = Eigen::Map<Eigen::VectorXd>(x.data(), x.size());
  }
  
  fopt = fbest_;
  xopt = xbest_;

  return return_status;
}

ErrorMatrixStatus NLOptOptimizer::error_matrix_estimate(MatRef error_matrix)
{
  return err_est_->error_matrix(error_matrix);
}

ErrorMatrixStatus NLOptOptimizer::
calc_error_matrix_and_eigenvectors(MatRef error_matrix,
                                   VecRef eigenvalues, MatRef eigenvectors)
{
  const unsigned npar { fcn_->num_domain_axes() };
  error_matrix.resize(npar,npar);
  Eigen::VectorXd error_hint;
  if(error_matrix_estimate(error_matrix) != ErrorMatrixStatus::UNAVAILABLE)
    error_hint = error_matrix.diagonal().array().sqrt();
  Eigen::MatrixXd hessian(npar,npar);
  hessian::calculate_hessian(*fcn_, xbest_, hessian, error_hint);
  return hessian::hessian_to_error_matrix(*fcn_, hessian, error_matrix,
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

  if(!isfinite(fcn_value))fcn_value = inf;
  
  if(verbose_ != VerbosityLevel::SILENT)
  {
    std::cout << std::left << std::setw(4) << iterations_+1 << ' '
              << std::fixed << std::setprecision(8) << fcn_value;
    for(unsigned ipar=0;ipar<n;ipar++)std::cout << ' ' << x[ipar];
    if(grad)
      for(unsigned ipar=0;ipar<n;ipar++)std::cout << ' ' << grad[ipar];
    std::cout << '\n';
  }
    
  if(iterations_ == 0 or fcn_value < fbest_)fbest_ = fcn_value, xbest_ = xvec;
  iterations_++;
  
  return fcn_value;
}
