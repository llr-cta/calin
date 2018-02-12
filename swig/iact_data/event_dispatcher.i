/*

   calin/iact_data/event_dispatcher.i -- Stephen Fegan -- 2016-02-10

   SWIG interface file for calin event dispatcher and visitor

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

%module (package="calin.iact_data", threads="1") event_dispatcher
%nothread;

%{
#include "iact_data/event_dispatcher.hpp"
#include "iact_data/parallel_event_dispatcher.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "calin_typemaps.i"
%import "calin_global_definitions.i"
%include "calin_global_config.hpp"

%import "iact_data/event_visitor.i"

%thread; // Release Pyhton GIL for all functions here (since some use threads)
%include "iact_data/event_dispatcher.hpp"
%include "iact_data/parallel_event_dispatcher.hpp"
%nothread;
