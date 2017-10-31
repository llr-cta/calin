/*

   calin/iact_data/instrument_layout.i -- Stephen Fegan -- 2017-09-15

   SWIG interface file for calin instrument (camera) layout

   Copyright 2017, Stephen Fegan <sfegan@llr.in2p3.fr>
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

%module (package="calin.iact_data") instrument_layout

%{
#include "iact_data/instrument_layout.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "calin_typemaps.i"
%import "calin_global_definitions.i"

%newobject make_grid_from_instrument_layout;

%import "iact_data/instrument_layout.pb.i"
%include "iact_data/instrument_layout.hpp"
