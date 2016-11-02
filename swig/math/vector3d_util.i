//-*-mode:swig;-*-

/*

   calin/math/ray.i -- Stephen Fegan -- 2016-10-24

   SWIG interface file for calin.math.ray

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

%module (package="calin.math") vector3d_util

%{
#include "math/vector3d_util.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "calin_typemaps.i"
%import "calin_global_definitions.i"

%import "common_types.pb.i"

//%newobject calin::math::rng::RNG::as_proto() const;
//%newobject calin::math::rng::RNGCore::as_proto() const;

%include "math/vector3d_util.hpp"
