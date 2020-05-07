/*

   calin/diagnostics/simple_charge_hists.i -- Stephen Fegan -- 2020-04-18

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

%module (package="calin.diagnostics") simple_charge_hists
%feature(autodoc,2);

%{
#include "iact_data/event_visitor.hpp"
#include "diagnostics/simple_charge_hists.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "calin_typemaps.i"
%import "calin_global_definitions.i"

%import "iact_data/event_visitor.i"
%import "diagnostics/simple_charge_hists.pb.i"

%newobject calin::diagnostics::simple_charge_hists::SimpleChargeHistsParallelEventVisitor::simple_charge_hists();

%include "diagnostics/simple_charge_hists.hpp"
