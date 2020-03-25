/*

   calin/diagnostics/simple_charge_stats.i -- Stephen Fegan -- 2020-03-25

   SWIG interface file for simple charge stats diagnostics

   Copyright 2020, Stephen Fegan <sfegan@llr.in2p3.fr>
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

%module (package="calin.diagnostics") simple_charge_stats
%feature(autodoc,2);

%{
#include "iact_data/event_visitor.hpp"
#include "diagnostics/simple_charge_stats.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "calin_typemaps.i"
%import "calin_global_definitions.i"

%import "iact_data/event_visitor.i"
%import "diagnostics/simple_charge_stats.pb.i"

%newobject calin::diagnostics::simple_charge_stats::SimpleChargeStatsVisitor::simple_charge_stats();

%include "diagnostics/simple_charge_stats.hpp"
