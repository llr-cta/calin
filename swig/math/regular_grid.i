//-*-mode:swig;-*-

/*

   calin/math/regular_grid.i -- Stephen Fegan -- 2017-08-28

   SWIG interface file for calin.math.regular_grid

   Copyright 2017, Stephen Fegan <sfegan@llr.in2p3.fr>
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

%module (package="calin.math") regular_grid

%{
#include "math/regular_grid.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "calin_typemaps.i"
%import "calin_global_definitions.i"

%apply double &OUTPUT { double& x, double& y };
%apply Eigen::VectorXd &OUTPUT { Eigen::VectorXd& xv, Eigen::VectorXd& yv };

%include "math/regular_grid.hpp"
