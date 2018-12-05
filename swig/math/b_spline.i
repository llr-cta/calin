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

%template(interval_basis0) calin::math::b_spline::interval_basis0<double>;
%template(interval_basis1) calin::math::b_spline::interval_basis1<double>;
%template(interval_basis2) calin::math::b_spline::interval_basis2<double>;
%template(interval_basis3) calin::math::b_spline::interval_basis3<double>;

%template(interval_value0) calin::math::b_spline::interval_value0<double>;
%template(interval_value1) calin::math::b_spline::interval_value1<double>;
%template(interval_value2) calin::math::b_spline::interval_value2<double>;
%template(interval_value3) calin::math::b_spline::interval_value3<double>;

%template(interval_integral0) calin::math::b_spline::interval_integral0<double>;
%template(interval_integral1) calin::math::b_spline::interval_integral1<double>;
%template(interval_integral2) calin::math::b_spline::interval_integral2<double>;
%template(interval_integral3) calin::math::b_spline::interval_integral3<double>;

%template(interval_derivative0) calin::math::b_spline::interval_derivative0<double>;
%template(interval_derivative1) calin::math::b_spline::interval_derivative1<double>;
%template(interval_derivative2) calin::math::b_spline::interval_derivative2<double>;
%template(interval_derivative3) calin::math::b_spline::interval_derivative3<double>;
