//-*-mode:swig;-*-

/* 

   calin/math/hessian.i -- Stephen Fegan -- 2015-04-27

*/

%module (package="calin.math") hessian

%{
#include "math/hessian.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "package_wide_definitions.i"

%import "math/function.i"

%include "typemaps.i"

%apply Eigen::MatrixXd &OUTPUT { Eigen::MatrixXd& error_matrix };
%apply Eigen::MatrixXd &OUTPUT { Eigen::MatrixXd& eigenvectors };
%apply Eigen::VectorXd &OUTPUT { Eigen::VectorXd& eigenvalues };
%apply Eigen::MatrixXd &OUTPUT { Eigen::MatrixXd& hessian };

%rename(hessian_to_error_matrix_no_eigenvectors)
hessian_to_error_matrix(function::MultiAxisFunction& fcn,
                        ConstMatRef hessian,
                        MatRef error_matrix);

%include "math/hessian.hpp"
