/* 

   calin/math/minuit75_optimizer.hpp -- Stephen Fegan -- 2015-03-12

   Object oriented interface to the 1975 FORTRAN version of MINUIT
   from CERN. See http://www.adsabs.harvard.edu/abs/1975CoPhC..10..343J

*/

#pragma once

#include "function.hpp"
#include "optimizer.hpp"

namespace calin { namespace math { namespace optimizer {

class Minuit75Optimizer: public Optimizer
{
 public:
  Minuit75Optimizer(MultiAxisFunction* fcn, bool adopt_fcn = false);
  ~Minuit75Optimizer();

  bool find_min(std::vector<double>& xopt) override;
  bool get_optimizer_error_matrix_estimate(Eigen::MatrixXd& err_mat) override;
  bool calc_error_matrix(Eigen::MatrixXd& err_mat) override;

 protected:
  int do_command(const std::string& command);
  static long minuit_callback(long* npar, double* grad, double* fcnval,
                              const double* x, long* iflag, void* futil);
  void eval_function(long npar, double* grad, double& fcnval,
                     const double* x, long& iflag);
  void* fcb_;
};

} } } // namespace calin::math::optimizer
