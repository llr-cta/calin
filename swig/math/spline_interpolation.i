//-*-mode:swig;-*-

/*

   calin/math/spline_interpolation.i -- Stephen Fegan -- 2018-12-18

   SWIG interface file for calin.math.spline_interpolation

   Copyright 2018, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, Ecole Polytechnique, CNRS/IN2P3

   This file is part of "calin"

   "calin" is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 or
   later, as published by the Free Software Foundation.

   "calin" is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

*/

%module (package="calin.math") spline_interpolation

%{
#include "math/spline_interpolation.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "calin_typemaps.i"
%import "calin_global_definitions.i"

//%apply double &OUTPUT { double& x, double& y };
//%apply Eigen::VectorXd &OUTPUT { Eigen::VectorXd& xv, Eigen::VectorXd& yv };

%include "math/spline_interpolation.hpp"
