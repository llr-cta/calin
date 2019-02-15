//-*-mode:swig;-*-

/*

   calin/simulation/detector_efficiency.i -- Stephen Fegan -- 2016-10-17

   SWIG interface file for calin.simulation.detector_efficiency

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

%module (package="calin.simulation") detector_efficiency
%feature(autodoc,2);

%{
#include "simulation/detector_efficiency.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "calin_typemaps.i"
%import "calin_global_definitions.i"

%import "math/interpolation_1d.i"

%include "simulation/detector_efficiency.hpp"
