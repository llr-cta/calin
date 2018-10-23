/*

   calin/math/b_spline.i -- Stephen Fegan -- 2018-10-22

   SWIG interface file for calin.math.b_spline

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

%module (package="calin.math") b_spline

%{
#include "math/b_spline.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "calin_typemaps.i"
%import "calin_global_definitions.i"

%apply double &OUTPUT { double& y0 };
%apply double &OUTPUT { double& y1 };
%apply double &OUTPUT { double& y2 };
%apply double &OUTPUT { double& y3 };

%include "math/b_spline.hpp"

%template(basis0) calin::math::b_spline::basis0<double>;
%template(basis1) calin::math::b_spline::basis1<double>;
%template(basis2) calin::math::b_spline::basis2<double>;
%template(basis3) calin::math::b_spline::basis3<double>;

%template(eval0) calin::math::b_spline::eval0<double>;
%template(eval1) calin::math::b_spline::eval1<double>;
%template(eval2) calin::math::b_spline::eval2<double>;
%template(eval3) calin::math::b_spline::eval3<double>;

%template(integrate0) calin::math::b_spline::integrate0<double>;
%template(integrate1) calin::math::b_spline::integrate1<double>;
%template(integrate2) calin::math::b_spline::integrate2<double>;
%template(integrate3) calin::math::b_spline::integrate3<double>;
