/*

   calin/diagnostics/clock_regression.i -- Stephen Fegan -- 2020-05-29

   SWIG interface file for clock regression tests

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

%module (package="calin.diagnostics") clock_regression
%feature(autodoc,2);

%{
#include "iact_data/event_visitor.hpp"
#include "diagnostics/clock_regression.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "calin_typemaps.i"
%import "calin_global_definitions.i"

%import "iact_data/event_visitor.i"
%import "diagnostics/clock_regression.pb.i"

%newobject calin::diagnostics::clock_regression::ClockRegressionParallelEventVisitor::clock_regression() const;

%include "diagnostics/clock_regression.hpp"
