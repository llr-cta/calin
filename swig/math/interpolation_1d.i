//-*-mode:swig;-*-

/*

   calin/math/interpolation_1d.i -- Stephen Fegan -- 2016-10-17

   SWIG interface file for calin.math.interpolation_1d

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
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

%module (package="calin.math") interpolation_1d

%{
#include "math/interpolation_1d.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "calin_typemaps.i"
%import "calin_global_definitions.i"

%include "math/interpolation_1d.hpp"

%template(LinearInterpolatorDouble) calin::math::interpolation_1d::LinearInterpolator<double>;
%template(InterpLinear1DDouble) calin::math::interpolation_1d::Interpolation1D<double, calin::math::interpolation_1d::LinearInterpolator<double> >;
