//-*-mode:swig;-*-

/*

   calin/calib/spe_fit.i -- Stephen Fegan -- 2015-04-24

   SWIG interface file for calin.calin.spe_fit

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
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

%module (package="calin.calib") spe_fit
%feature(autodoc,2);

%{
#include "calib/spe_fit.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "calin_typemaps.i"
%import "calin_global_definitions.i"

%import "math/function.i"
%import "math/histogram.i"

%import "calib/spe_fit.pb.i"
%include "calib/spe_fit.hpp"
