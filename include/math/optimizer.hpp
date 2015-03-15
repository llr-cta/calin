/* 

   calin/math/optimizer.hpp -- Stephen Fegan -- 2015-03-12

   Base class providing interface to function optimizers (minimizers)

*/

#pragma once

#include <vector>
#include <Eigen/Core>

#include "function.hpp"

namespace calin { namespace math { namespace optimizer {

enum class OptimizerVerbosityLevel { SILENT, FCN_EVALS_ONLY, ELEVATED, MAX };

class Optimizer
{
 public:
  using VerbosityLevel = OptimizerVerbosityLevel;
  Optimizer(MultiAxisFunction* fcn, bool adopt_fcn = false):
      fcn_(fcn), my_fcn_(adopt_fcn) { /* nothing to see here */ }
  virtual ~Optimizer();

  virtual bool find_min(std::vector<double>& xopt) = 0;
  virtual bool get_optimizer_error_matrix_estimate(Eigen::MatrixXd& err_mat) = 0;
  virtual bool calc_error_matrix(Eigen::MatrixXd& err_mat) = 0;

  void set_verbosity_level(VerbosityLevel verbose = VerbosityLevel::SILENT) { verbose_=verbose; }
  void set_initial_values(const std::vector<double>& x0 = {}) { x0_=x0; }
  void set_scale(const std::vector<double>& xscale = {}) { xscale_=xscale; }

 protected:
  bool my_fcn_ { false };
  MultiAxisFunction* fcn_  { nullptr };
  VerbosityLevel verbose_ { VerbosityLevel::SILENT };
  std::vector<double> x0_ { };
  std::vector<double> xscale_ { };
};

} // namespace optimizer

using Optimizer = optimizer::Optimizer;

} } // namespace calin::math
