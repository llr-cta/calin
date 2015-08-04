//-*-mode:swig;-*-

/* 

   calin/math/function.i -- Stephen Fegan -- 2015-04-15

*/

%module (package="calin.math") function

%{
#include "math/function.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "package_wide_definitions.i"

%include "typemaps.i"
%apply Eigen::VectorXd &OUTPUT { Eigen::VectorXd& gradient };
%apply Eigen::MatrixXd &OUTPUT { Eigen::MatrixXd& hessian };
%apply double &OUTPUT { double& dfdx, double& d2fdx2 }

%include "math/function.hpp"

%template (VectorParameterAxis) std::vector<calin::math::function::ParameterAxis>;
