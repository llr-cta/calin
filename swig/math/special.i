//-*-mode:swig;-*-

/* 

   calin/math/special.i -- Stephen Fegan -- 2015-08-06

   SWIG interface file for calin.math.special

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

%module (package="calin.math") special

%{
#include "math/special.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "numpy.i"
%import "package_wide_definitions.i"
%include "calin_typemaps.i"

%include "math/special.hpp"
