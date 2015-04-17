/* 

   calin/math/optimizer.cpp -- Stephen Fegan -- 2015-03-12

*/

#include <sstream>
#include <cstdio>
#include <cmath>
#include <limits>
#include <iostream>
#include <iomanip>

#include "math/minuit75_optimizer.hpp"

#include "f2c/f2c.h"

extern "C" {

  int s_stop(const char *s, ftnlen n);
  double d_lg10(doublereal *x);
  double pow_ri(real *ap, integer *bp);
  integer i_len(const char *s, ftnlen n);
  integer e_rsfi(void);
  integer f_open(olist *a);
  integer f_rew(alist *a);
  integer f_inqu(inlist *a);
  integer s_rsfe(cilist *a);
  integer e_rsfe(void);
  integer s_rsfi(icilist *a);
  void s_cat(char *lp, const char *rpp[], ftnint rnp[], ftnint *np, ftnlen ll);
  double d_sign(doublereal *a, doublereal *b);
  integer s_cmp(const char *a0, const char *b0, ftnlen la, ftnlen lb);
  integer e_wsfi(void);
  integer s_wsfi(icilist *a);
  integer e_wsfe(void);
  integer s_wsfe(cilist *a);
  integer do_fio(ftnint *number, const char *ptr, ftnlen len);
  void s_copy(register char *a, register const char *b, ftnlen la, ftnlen lb);
  integer i_indx(const char *a, const char *b, ftnlen la, ftnlen lb);
  integer e_wsfe(void);
    
};

#define MNE 400
#define MNI 200

namespace {

#include "minuit75_routines_f2c.cpp"

static mn7fcb* cast_fcb(void* x)
{
  return static_cast<mn7fcb*>(x);
}

};

#undef F2C_INCLUDE
#undef qbit_clear
#undef qbit_set
#undef TRUE_
#undef FALSE_
#undef Extern
#undef VOID
#undef abs
#undef dabs
#undef min
#undef max
#undef dmin
#undef dmax
#undef bit_test
#undef bit_clear
#undef bit_set
#undef F2C_proc_par_types

using namespace calin::math::optimizer;

Minuit75Optimizer::
Minuit75Optimizer(function::MultiAxisFunction* fcn, bool adopt_fcn):
    Optimizer(fcn, adopt_fcn), fcb_(new mn7fcb)
{
  integer i5=5, i6=6, i7=7;
  mninit_(&i5, &i6, &i7, cast_fcb(fcb_));
}

Minuit75Optimizer::~Minuit75Optimizer()
{
  delete cast_fcb(fcb_);
}

bool Minuit75Optimizer::requires_gradient()
{
  return false;
}

bool Minuit75Optimizer::requires_hessian()
{
  return false;
}

bool Minuit75Optimizer::requires_box_constraints()
{
  return false;
}

bool Minuit75Optimizer::can_estimate_error()
{
  return true;
}

bool Minuit75Optimizer::can_use_gradient()
{
  return true;
}

bool Minuit75Optimizer::can_use_hessian()
{
  return false;
}

bool Minuit75Optimizer::can_impose_box_constraints()
{
  return true;
}

bool Minuit75Optimizer::minimize(VecRef xopt, double& fopt)
{
  constexpr auto inf = std::numeric_limits<double>::infinity();
  mn7fcb* fcb { cast_fcb(fcb_) };
  integer error_flag;
  std::ostringstream pline;

  integer i5=5, i6=6, i7=7;
  mninit_(&i5, &i6, &i7, cast_fcb(fcb_));
  mintio_(&i5, &i6, &i7, fcb);

  integer verbose = -1;
  if(verbose_ == VerbosityLevel::ELEVATED)verbose=0;
  else if(verbose_ == VerbosityLevel::MAX)verbose=1;
  pline << "SET PRINT " << verbose;
  do_command(pline.str()); // Set verbosity of Minuit

  do_command("SET NOWARN");

  pline.str(std::string());
  pline << "SET ERR " << fcn_->error_up();
  do_command(pline.str());

  if(fcn_->can_calculate_gradient())
    do_command("SET GRAD 1");  // Use gradient calculated by fcn

  auto params = fcn_->domain_axes();
  for(unsigned iparam=0; iparam<params.size(); iparam++)
    {
      integer iparam_fortran { iparam+1 };
      double x0 = params[iparam].initial_value;
      if(x0_.size()==params.size()) x0 = x0_[iparam];
      double xscale = params[iparam].scale;
      if(xscale_.size()==params.size()) xscale = xscale_[iparam];
      double xlo = 0;
      // params[iparam].has_lo_bound?params[iparam].lo_bound:-inf;
      double xhi = 0;
      //params[iparam].has_hi_bound?params[iparam].hi_bound:inf;
      if(!params[iparam].has_lo_bound && !params[iparam].has_hi_bound)
        xlo = xhi = 0.0;
      mnparm_(&iparam_fortran, params[iparam].name.c_str(),
              &x0, &xscale, &xlo, &xhi, &error_flag,
              params[iparam].name.size(), fcb);
    }

  do_command("MIN 0 1.0");

  xopt.resize(params.size());
  for(unsigned iparam=0;iparam<params.size();iparam++)
    {
      char buffer[11];
      buffer[10]='\0';
      integer fortran_iparam = iparam+1;
      integer intvar;
      double xerr;
      double xlo;
      double xhi;
      mnpout_(&fortran_iparam, buffer, &xopt[iparam], &xerr, &xlo, &xhi,
              &intvar, sizeof(buffer)/sizeof(*buffer)-1, fcb);
    }
  
  return true;
}

ErrorMatrixStatus Minuit75Optimizer::error_matrix_estimate(MatRef err_mat)
{

}

ErrorMatrixStatus Minuit75Optimizer::calc_error_matrix(MatRef err_mat)
{

}

int Minuit75Optimizer::do_command(const std::string& command)
{
  // Pass a command string to Minuit
  integer errorFlag = 0;  
  mncomd_(minuit_callback, command.c_str(), &errorFlag,
          this, command.length(), cast_fcb(fcb_));
  return errorFlag;
}

integer Minuit75Optimizer::
minuit_callback(integer* npar, double* grad, double* fcnval,
                const double* x, long* iflag, void* fcn_data)
{
  // This is the function that Minuit minimizes.
  Minuit75Optimizer* m75 = static_cast<Minuit75Optimizer*>(fcn_data);
  m75->eval_function(*npar,grad,*fcnval,x,*iflag);
  return 0;
}

void Minuit75Optimizer::
eval_function(long npar, double* grad, double& fcnval,
              const double* x, long& iflag)
{
  Eigen::VectorXd xvec = Eigen::Map<const Eigen::VectorXd>(x,npar);
  
  if(iflag==2)
  {
    Eigen::VectorXd gvec(npar);
    fcnval = fcn_->value_and_gradient(xvec, gvec);
    Eigen::Map<Eigen::VectorXd>(grad,npar) = gvec;
  }
  else
  {
    fcnval = fcn_->value(xvec);
  }

  std::cout << fcnval;
  for(unsigned ipar=0;ipar<npar;ipar++)std::cout << ' ' << xvec(ipar);
  std::cout << '\n';
  
}
