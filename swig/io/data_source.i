/*

   calin/iact_data/data_source.i -- Stephen Fegan -- 2019-12-17

   SWIG interface file for calin.io.data_source. Mostly just for FragmentList
   interface

   Copyright 2019, Stephen Fegan <sfegan@llr.in2p3.fr>
   Laboratoire Leprince-Ringuet, CNRS/IN2P3, Ecole Polytechnique, Institut Polytechnique de Paris

   This file is part of "calin"

   "calin" is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 or
   later, as published by the Free Software Foundation.

   "calin" is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

*/

%module (package="calin.io") data_source
%feature(autodoc,2);

%{
#include <calin_global_config.hpp>
#include <io/data_source.hpp>
#include <io/chained_data_source.hpp>
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "calin_typemaps.i"
%include "calin_global_config.hpp"
%import "calin_global_definitions.i"

%include "io/data_source.hpp"
%include "io/chained_data_source.hpp"
