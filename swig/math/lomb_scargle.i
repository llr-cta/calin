//-*-mode:swig;-*-

/*

   calin/math/lomb_scargle.i -- Stephen Fegan -- 2017-12-05

   SWIG interface file for calin.math.lomb_scargle

   Copyright 2017, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, CNRS, Ecole Polytechnique, Institut Polytechnique de Paris

   This file is part of "calin"

   "calin" is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 or
   later, as published by the Free Software Foundation.

   "calin" is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

*/

%module (package="calin.math") lomb_scargle
%feature(autodoc,2);

%{
#include "math/lomb_scargle.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "calin_typemaps.i"
%import "calin_global_definitions.i"

//%apply double &OUTPUT { double& x, double& y };
//%apply Eigen::VectorXd &OUTPUT { Eigen::VectorXd& xv, Eigen::VectorXd& yv };

%include "math/lomb_scargle.hpp"
