/*

   calin/math/nlopt_optimizer.cpp -- Stephen Fegan -- 2015-03-16

   Interface to MIT NLOpt optimizer suite

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, Ecole Polytechnique, CNRS/IN2P3

   This file is part of "calin"

   "calin" is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 or
   later, as published by the Free Software Foundation.

   "calin" is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

*/

#include <iostream>
#include <iomanip>
#include <cstring>
#include <cmath>

#include <nlopt/nlopt.h>

#include <io/log.hpp>
#include <math/nlopt_optimizer.hpp>
#include <math/hessian.hpp>

using namespace calin::math::optimizer;
using namespace calin::io::log;

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

static nlopt_algorithm name_to_algorithm(const std::string& name)
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

#undef MATCH_ALGORITHM

#if 0
#define MATCH_ALGORITHM(x) case x: return(std::string(#x))

static std::string ZZZ_algorithm_to_name(nlopt_algorithm algo)
{
  switch(algo)
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

    case NLOPT_NUM_ALGORITHMS:
    default:
      return std::string();
  };

  return std::string();
}

#undef MATCH_ALGORITHM
#endif // #if 0 around ZZZ_algorithm_to_name
}

NLOptOptimizer::
NLOptOptimizer(NLOptOptimizer::algorithm_type algorithm,
               function::MultiAxisFunction* fcn, bool adopt_fcn)
    : Optimizer(fcn, adopt_fcn), algorithm_(algorithm)
{
  if(is_local_optimizer(algorithm) and can_use_gradient(algorithm))
    err_est_.reset(new BFGSErrorMatrixEstimator(fcn->error_up()));
  else
    err_est_.reset(new IdentityErrorMatrixEstimator(fcn->error_up()));
}

NLOptOptimizer::
NLOptOptimizer(const std::string& algorithm_name,
               function::MultiAxisFunction* fcn, bool adopt_fcn)
    : NLOptOptimizer(algorithm_by_name(algorithm_name),
                     fcn, adopt_fcn)
{
  // nothing to see here
}

NLOptOptimizer::~NLOptOptimizer()
{
  // nothing to see here
}

NLOptOptimizer::algorithm_type NLOptOptimizer::
algorithm_by_name(const std::string& algorithm_name)
{
  return name_to_algorithm(algorithm_name);
}


bool NLOptOptimizer::is_valid_algorithm(const std::string& algorithm_name)
{
  return name_to_algorithm(algorithm_name) != NLOPT_NUM_ALGORITHMS;
}

bool NLOptOptimizer::is_local_optimizer(algorithm_type algorithm)
{
  switch(algorithm)
  {
    case NLOPT_GN_DIRECT:
    case NLOPT_GN_DIRECT_L:
    case NLOPT_GN_DIRECT_L_RAND:
    case NLOPT_GN_DIRECT_NOSCAL:
    case NLOPT_GN_DIRECT_L_NOSCAL:
    case NLOPT_GN_DIRECT_L_RAND_NOSCAL:
    case NLOPT_GN_ORIG_DIRECT:
    case NLOPT_GN_ORIG_DIRECT_L:
    case NLOPT_GD_STOGO:
    case NLOPT_GD_STOGO_RAND:
    case NLOPT_GN_CRS2_LM:
    case NLOPT_GN_MLSL:
    case NLOPT_GD_MLSL:
    case NLOPT_GN_MLSL_LDS:
    case NLOPT_GD_MLSL_LDS:
    case NLOPT_GN_ISRES:
    case NLOPT_G_MLSL:
    case NLOPT_G_MLSL_LDS:
    case NLOPT_GN_ESCH:
      return false;

    case NLOPT_LD_LBFGS_NOCEDAL:
    case NLOPT_LD_LBFGS:
    case NLOPT_LN_PRAXIS:
    case NLOPT_LD_VAR1:
    case NLOPT_LD_VAR2:
    case NLOPT_LD_TNEWTON:
    case NLOPT_LD_TNEWTON_RESTART:
    case NLOPT_LD_TNEWTON_PRECOND:
    case NLOPT_LD_TNEWTON_PRECOND_RESTART:
    case NLOPT_LD_MMA:
    case NLOPT_LN_COBYLA:
    case NLOPT_LN_NEWUOA:
    case NLOPT_LN_NEWUOA_BOUND:
    case NLOPT_LN_NELDERMEAD:
    case NLOPT_LN_SBPLX:
    case NLOPT_LN_AUGLAG:
    case NLOPT_LD_AUGLAG:
    case NLOPT_LN_AUGLAG_EQ:
    case NLOPT_LD_AUGLAG_EQ:
    case NLOPT_LN_BOBYQA:
    case NLOPT_AUGLAG:
    case NLOPT_AUGLAG_EQ:
    case NLOPT_LD_SLSQP:
    case NLOPT_LD_CCSAQ:
      return true;

    case NLOPT_NUM_ALGORITHMS:
    default:
      assert(0);
  }

  assert(0);
  return false;
}

bool NLOptOptimizer::requires_gradient(algorithm_type algorithm)
{
  switch(algorithm)
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
  return false;
}

bool NLOptOptimizer::requires_box_constraints(algorithm_type algorithm)
{
  switch(algorithm)
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
  return false;
}


bool NLOptOptimizer::requires_hessian(algorithm_type algorithm)
{
  return false;
}

bool NLOptOptimizer::can_estimate_error(algorithm_type algorithm)
{
  return requires_gradient(algorithm);
}

bool NLOptOptimizer::can_use_gradient(algorithm_type algorithm)
{
  return requires_gradient(algorithm);
}

bool NLOptOptimizer::can_use_hessian(algorithm_type algorithm)
{
  return false; // (a==LD_CCSAQ);
}

bool NLOptOptimizer::can_impose_box_constraints(algorithm_type algorithm)
{
  return true;
}

bool NLOptOptimizer::is_local_optimizer()
{
  return is_local_optimizer(algorithm_);
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

OptimizationStatus NLOptOptimizer::minimize(VecRef xopt, double& fopt)
{
  nlopt_result nlopt_status;

  unsigned naxes { fcn_->num_domain_axes() };

  err_est_->reset(naxes);
  gvec_.resize(naxes);

  std::vector<double> xlim_lo { limits_lo() };
  std::vector<double> xlim_hi { limits_hi() };
  std::vector<double> stepsize { initial_stepsize() };

  opt_starting(nlopt_algorithm_name(algorithm_), requires_gradient(),
               requires_hessian(), xlim_lo, xlim_hi, stepsize);

  fopt = fbest_;
  xopt = xbest_;

  std::unique_ptr<nlopt_opt_s, void(*)(nlopt_opt)>
      opt(nlopt_create(algorithm_, naxes), nlopt_destroy);
  if(opt.get() == NULL)
  {
    opt_message_ = "Could not create NLOpt object";
    return opt_status_;
  }

  nlopt_status = nlopt_set_min_objective(opt.get(), nlopt_callback, this);
  if(nlopt_status != NLOPT_SUCCESS)
  {
    opt_message_ = "Could not set objective function";
    return opt_status_;
  }

  if(!xlim_lo.empty() &&
     std::count(xlim_lo.begin(), xlim_lo.end(), neg_inf) !=
     (unsigned)xlim_lo.size())
  {
    nlopt_status = nlopt_set_lower_bounds(opt.get(), &xlim_lo.front());
    if(nlopt_status != NLOPT_SUCCESS)
    {
      opt_message_ = "Could not set lower bounds";
      return opt_status_;
    }
  }

  if(!xlim_hi.empty() &&
     std::count(xlim_hi.begin(), xlim_hi.end(), pos_inf) !=
     (unsigned)xlim_hi.size())
  {
    nlopt_status = nlopt_set_upper_bounds(opt.get(), &xlim_hi.front());
    if(nlopt_status != NLOPT_SUCCESS)
    {
      opt_message_ = "Could not set upper bounds";
      return opt_status_;
    }
  }

  nlopt_status = nlopt_set_initial_step(opt.get(), &stepsize.front());
  if(nlopt_status != NLOPT_SUCCESS)
  {
    opt_message_ = "Could not set initial stepsize";
    return opt_status_;
  }

  if((abs_tolerance()>0 and
      nlopt_set_ftol_abs(opt.get(), abs_tolerance()) != NLOPT_SUCCESS) or
     (rel_tolerance()>0 and
      nlopt_set_ftol_rel(opt.get(), rel_tolerance()) != NLOPT_SUCCESS) or
     (max_iterations()>0 and
      nlopt_set_maxeval(opt.get(), max_iterations()) != NLOPT_SUCCESS) or
     (max_walltime()>0 and
      nlopt_set_maxtime(opt.get(), max_walltime()) != NLOPT_SUCCESS))
  {
    opt_message_ = "Could not set stopping criteria";
    return opt_status_;
  }

  nlopt_status = nlopt_optimize(opt.get(), xopt.data(), &fopt);

  switch(nlopt_status)
  {
    case NLOPT_SUCCESS:
    case NLOPT_FTOL_REACHED:
      opt_status_ = OptimizationStatus::TOLERANCE_REACHED;
      opt_message_ = "Tolerance reached";
      break;
    case NLOPT_MAXEVAL_REACHED:
      opt_status_ = OptimizationStatus::STOPPED_AT_MAXCALLS;
      opt_message_ = "Calls limit reached";
      break;
    case NLOPT_MAXTIME_REACHED:
      opt_status_ = OptimizationStatus::STOPPED_AT_MAXTIME;
      opt_message_ = "Time limit reached";
      break;
    case NLOPT_ROUNDOFF_LIMITED:
      opt_status_ = OptimizationStatus::LIMITED_BY_PRECISION;
      opt_message_ = "Minimization limited by machine precision";
      break;
    case NLOPT_FAILURE:
      opt_status_ = OptimizationStatus::OPTIMIZER_FAILURE;
      opt_message_ = "Optimizer failed for unknown reason";
      break;
    case NLOPT_INVALID_ARGS:
      opt_status_ = OptimizationStatus::OPTIMIZER_FAILURE;
      opt_message_ = "Invalid arguments";
      break;
    case NLOPT_OUT_OF_MEMORY:
      opt_status_ = OptimizationStatus::OPTIMIZER_FAILURE;
      opt_message_ = "Out of memory";
      break;
    case NLOPT_STOPVAL_REACHED:
    case NLOPT_XTOL_REACHED:
    case NLOPT_FORCED_STOP:
      opt_status_ = OptimizationStatus::OPTIMIZER_FAILURE;
      opt_message_ = "This should not happen";
      break;
  }

  if(opt_status_ == OptimizationStatus::OPTIMIZER_FAILURE)
  {
    fopt = fbest_;
    xopt = xbest_;
  }
  else
  {
    fbest_ = fopt;
    xbest_ = xopt;
  }

  opt_finished(opt_status_, fbest_, xbest_);

  return opt_status_;
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
  if(xscale_.size() == fcn_->num_domain_axes())
    error_hint = calin::std_to_eigenvec(xscale_);
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
      fcn_value = fcn_->value_and_gradient(xvec,gvec_);
      Eigen::Map<Eigen::VectorXd>(grad,n) = gvec_;
      err_est_->incorporate_func_gradient(xvec, fcn_value, gvec_);
    }
    else
    {
      fcn_value = fcn_->value(xvec);
      err_est_->incorporate_func_value(xvec, fcn_value);
    }
  }
  catch(std::exception& x)
  {
    LOG(ERROR)
        << "NLOptOptimizer::eval_func: exception in user function with x = "
        << xvec.transpose() << '\n'
        << "NLOptOptimizer::eval_func: exception is: " << x.what() << '\n';
    throw;
  }

  if(!std::isfinite(fcn_value))fcn_value = inf;

#if 0
  if(dynamic_cast<BFGSErrorMatrixEstimator*>(err_est_.get()))
  {
    Eigen::MatrixXd err_mat_est;
    err_est_->error_matrix(err_mat_est);
    opt_progress(fcn_value, xvec, grad?&gvec_:nullptr, nullptr, &err_mat_est);
  }
  else;
#endif

  opt_progress(fcn_value, xvec, grad?&gvec_:nullptr, nullptr);

  return fcn_value;
}
