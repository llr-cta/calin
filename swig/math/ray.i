//-*-mode:swig;-*-

/*

   calin/math/ray.i -- Stephen Fegan -- 2016-10-24

   SWIG interface file for calin.math.ray

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

%module (package="calin.math") ray

%{
#include "math/ray.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "calin_typemaps.i"
%import "calin_global_definitions.i"

//%import "math/ray.pb.i"

//%newobject calin::math::rng::RNG::as_proto() const;
//%newobject calin::math::rng::RNGCore::as_proto() const;

%include "math/ray.hpp"
