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

CMinpackOptimizer::
CMinpackOptimizer(function::MultiAxisFunction* fcn, bool adopt_fcn):
    Optimizer(fcn, adopt_fcn)
{
  // nothing to see here
}

CMinpackOptimizer::~CMinpackOptimizer()
{
  // nothing to see here
}

bool CMinpackOptimizer::is_local_optimizer()
{
  return true;
}

bool CMinpackOptimizer::requires_gradient()
{
  return true;
}

bool CMinpackOptimizer::requires_hessian()
{
  return true; // mode_==CMinpackMode::REQUIRE_HESSIAN;
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
  return true; // mode_!=CMinpackMode::DONT_USE_HESSIAN
}

bool CMinpackOptimizer::can_impose_box_constraints()
{
  return false;
}

OptimizationStatus CMinpackOptimizer::minimize(VecRef xopt, double& fopt)
{
  const double xtol { std::sqrt(std::numeric_limits<double>::epsilon()) };
  const int maxfev = max_iterations()>0?max_iterations():0x7FFFFFFF;

  opt_starting("CMinpack", requires_gradient(),
               requires_hessian(), {}, {}, {});

  unsigned naxes = fcn_->num_domain_axes();
  Eigen::VectorXd xvec = xbest_;
  Eigen::VectorXd gvec(naxes);
  Eigen::MatrixXd hmat(naxes, naxes);

  Eigen::VectorXd diag(naxes);
  int ngrad { 0 };
  int nhess { 0 };
  Eigen::VectorXd rvec(naxes*(naxes+1)/2);
  Eigen::VectorXd qtfvec(naxes);
  Eigen::VectorXd wa1(naxes);
  Eigen::VectorXd wa2(naxes);
  Eigen::VectorXd wa3(naxes);
  Eigen::VectorXd wa4(naxes);
  
  int istat = hybrj(cminpack_callback_hess,     // fcn
                    this,                       // void* p
                    naxes,                      // int n
                    xvec.data(),                // double* x
                    gvec.data(),                // double* fvec
                    hmat.data(),                // double& fjac
                    naxes,                      // int ldfjac
                    xtol,                       // double xtol
                    maxfev,                     // int maxfev
                    diag.data(),                // double* diag,
                    1,                          // int mode
                    100.0,                      // double factor
                    0,                          // int nprint
                    &ngrad,                     // int* nfev
                    &nhess,                     // int& njev
                    rvec.data(),                // double* r
                    rvec.size(),                // int lr
                    qtfvec.data(),              // double* qtf
                    wa1.data(),                 // double* wa1
                    wa2.data(),                 // double* wa2
                    wa3.data(),                 // double* wa3
                    wa4.data());                // double* wa4

}
  
ErrorMatrixStatus CMinpackOptimizer::error_matrix_estimate(MatRef error_matrix)
{

}

ErrorMatrixStatus CMinpackOptimizer::
calc_error_matrix_and_eigenvectors(MatRef error_matrix,
                                   VecRef eigenvalues, MatRef eigenvectors)
{
  
}

int CMinpackOptimizer::
cminpack_callback_grad(void* self, int n, const double* x,
                       double* grad, int iflag)
{
  CMinpackOptimizer* that = static_cast<CMinpackOptimizer*>(self);
  return that->eval_func(n,x,grad,nullptr,iflag);
}

int CMinpackOptimizer::
cminpack_callback_hess(void *self, int n, const double* x,
                       double* grad, double* hess, int ldfjac, int iflag)
{
  CMinpackOptimizer* that = static_cast<CMinpackOptimizer*>(self);
  (void)ldfjac;
  return that->eval_func(n,x,grad,hess,iflag);
}
  
int CMinpackOptimizer::
eval_func(unsigned n, const double* x, double* grad, double* hess, int iflag)
{
  xvec_ = Eigen::Map<const Eigen::VectorXd>(x,n);
  double last_fbest = fbest_;
  if(iflag == 2)
  {
    fval_ = fcn_->value_gradient_and_hessian(xvec_,gvec_,hmat_);
    this->opt_progress(fval_, xvec_, &gvec_, &hmat_);
    Eigen::Map<Eigen::MatrixXd>(hess,n,n) = hmat_;    
  }
  else // (iflag == 1)
  {
    fval_ = fcn_->value_and_gradient(xvec_,gvec_);
    this->opt_progress(fval_, xvec_, &gvec_, nullptr);
    Eigen::Map<Eigen::VectorXd>(grad,n) = gvec_;
  }
}

