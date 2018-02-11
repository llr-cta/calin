/*

   calin/provenance/system_info.i -- Stephen Fegan -- 2016-03-21

   SWIG interface file for provenance system build and host info

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

%module (package="calin.provenance") system_info

%{
#include "provenance/system_info.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "calin_typemaps.i"
%import "calin_global_definitions.i"

%import "util/log.pb.i"
%import "provenance/system_info.pb.i"

%newobject *::copy_the_build_info();
%newobject *::copy_the_host_info();
%include "provenance/system_info.hpp"
