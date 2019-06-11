//-*-mode:swig;-*-

/*

   calin/math/rng.i -- Stephen Fegan -- 2015-04-15

   SWIG interface file for calin.math.rng

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

%module (package="calin.math") rng
%feature(autodoc,2);

%{
#include "math/rng.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "calin_typemaps.i"
%import "calin_global_definitions.i"

%import "math/rng.pb.i"

%newobject calin::math::rng::RNG::as_proto() const;
%newobject calin::math::rng::RNGCore::as_proto() const;
%include "math/rng.hpp"

%template(NR3_EmulateSIMD_RNGCore_4) calin::math::rng::NR3_EmulateSIMD_RNGCore<4>;
