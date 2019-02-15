//-*-mode:swig;-*-

/*

   calin/math/pdf_1d.i -- Stephen Fegan -- 2015-04-16

   SWIG interface file for calin.math.pdf_1d

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
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

%module (package="calin.math") pdf_1d
%feature(autodoc,2);

%{
#include "math/pdf_1d.hpp"
#include "math/log_quadratic_spline_pdf_1d.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "calin_typemaps.i"
%import "calin_global_definitions.i"

%include "typemaps.i"

%import "math/function.i"

%include "math/pdf_1d.hpp"
%include "math/log_quadratic_spline_pdf_1d.hpp"
