//-*-mode:swig;-*-

/*

   calin/math/function.i -- Stephen Fegan -- 2015-04-15

   SWIG interface file for calin.math.function

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

%module (package="calin.math") function

%{
#include "math/function.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "calin_typemaps.i"
%import "calin_global_definitions.i"

%include "typemaps.i"

%import "pattern/delegation.hpp"

%apply Eigen::VectorXd &OUTPUT { Eigen::VectorXd& gradient };
%apply Eigen::MatrixXd &OUTPUT { Eigen::MatrixXd& hessian };
%apply double &OUTPUT { double& dfdx, double& d2fdx2 }

%include "math/function.hpp"

%template (VectorParameterAxis) std::vector<calin::math::function::ParameterAxis>;

%template (DelegatorSingleAxisFunction)
  calin::pattern::delegation::Delegator<
    calin::math::function::SingleAxisFunction>;
%template (DelegatorParameterizableMultiAxisFunction)
  calin::pattern::delegation::Delegator<
    calin::math::function::ParameterizableSingleAxisFunction>;

%template (BasicParameterizableDelegator_ParameterizableSingleAxisFunction)
  calin::math::function::BasicParameterizableDelegator<
    calin::math::function::ParameterizableSingleAxisFunction>;
%template (BasicSingleAxisFunctionDelegator_ParameterizableSingleAxisFunction)
  calin::math::function::BasicSingleAxisFunctionDelegator<
    calin::math::function::ParameterizableSingleAxisFunction>;

%template (SingleToMultiAxisFunctionAdapter)
  calin::math::function::BasicSingleToMultiAxisFunctionAdapter<
    calin::math::function::SingleAxisFunction>;
%template (SingleToMultiAxisFunctionAdapter_ParameterizableSingleAxisFunction)
  calin::math::function::BasicSingleToMultiAxisFunctionAdapter<
    calin::math::function::ParameterizableSingleAxisFunction>;
%template (SingleToParameterizableMultiAxisFunctionAdapter)
  calin::math::function::BasicSingleToParameterizableMultiAxisFunctionAdapter<
    calin::math::function::ParameterizableSingleAxisFunction>;

%template (ReducedSpaceParameterizable_ParameterizableSingleAxisFunction)
  calin::math::function::BasicReducedSpaceParameterizable<
    calin::math::function::ParameterizableSingleAxisFunction>;
%template (ReducedSpaceParameterizableSingleAxisFunction)
  calin::math::function::BasicReducedSpaceParameterizableSingleAxisFunction<
    calin::math::function::ParameterizableSingleAxisFunction>;
