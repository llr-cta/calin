//-*-mode:swig;-*-

/*

   calin/math/hessian.i -- Stephen Fegan -- 2015-04-27

   SWIG interface file for calin.math.hessian

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, Ecole polytechnique, CNRS/IN2P3, Universite Paris-Saclay

   This file is part of "calin"

   "calin" is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 or
   later, as published by the Free Software Foundation.

   "calin" is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

*/

%module (package="calin.math") hessian

%{
#include "math/hessian.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

//%include "numpy.i"
%include "calin_typemaps.i"
%import "calin_global_definitions.i"

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
