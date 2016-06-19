//-*-mode:swig;-*-

/*

   calin/math/optimizer.i -- Stephen Fegan -- 2015-04-23

   SWIG interface file for calin.math.optimizer

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

%module (package="calin.math") optimizer

%{
#include "math/optimizer.hpp"
#include "math/nlopt_optimizer.hpp"
#include "math/cminpack_optimizer.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "calin_typemaps.i"
%import "calin_global_definitions.i"

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

%newobject calin::math::optimizer::Optimizer::create_optimizer_for_function;
%newobject calin::math::optimizer::Optimizer::create_optimizer_by_name;

%include "math/optimizer.hpp"
%include "math/nlopt_optimizer.hpp"
%include "math/cminpack_optimizer.hpp"
