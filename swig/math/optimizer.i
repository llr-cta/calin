//-*-mode:swig;-*-

/* 

   calin/math/optimizer.i -- Stephen Fegan -- 2015-04-23

*/

%module (package="calin.math") optimizer

%{
#include "math/optimizer.hpp"
#include "math/nlopt_optimizer.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "package_wide_definitions.i"

%import "math/function.i"
%import "math/hessian.i"

//%include "nlopt/nlopt.hpp"

%include "typemaps.i"

%apply Eigen::VectorXd &OUTPUT { Eigen::VectorXd& xbest };
%apply Eigen::VectorXd &OUTPUT { Eigen::VectorXd& xopt };
%apply double &OUTPUT { double& fopt };
//%apply Eigen::MatrixXd &OUTPUT { Eigen::MatrixXd& error_matrix };
//%apply Eigen::MatrixXd &OUTPUT { Eigen::MatrixXd& eigenvectors };
//%apply Eigen::VectorXd &OUTPUT { Eigen::VectorXd& eigenvalues };

%include "math/optimizer.hpp"
%include "math/nlopt_optimizer.hpp"
