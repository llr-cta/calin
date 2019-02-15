/*

   calin/diagnostics/waveform.i -- Stephen Fegan -- 2016-03-23

   SWIG interface file for calin waveform diagnostics

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

%module (package="calin.diagnostics") functional
%feature(autodoc,2);

%{
#include "diagnostics/functional.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "calin_typemaps.i"
%import "calin_global_definitions.i"

%import "iact_data/event_dispatcher.i"
%import "iact_data/functional_event_visitor.i"
%import "diagnostics/value_capture.i"
%import "diagnostics/functional.pb.i"

%apply double &OUTPUT { double& common_variance_out };

%include "diagnostics/functional.hpp"

%template(FunctionalIntStatsVisitor)
calin::diagnostics::functional::FunctionalStatsVisitor<
  calin::iact_data::functional_event_visitor::
    DualGainInt32FunctionalTelescopeEventVisitor,
  calin::ix::diagnostics::functional::CameraIntFunctionalRawStats>;

%template(FunctionalDoubleStatsVisitor)
calin::diagnostics::functional::FunctionalStatsVisitor<
  calin::iact_data::functional_event_visitor::
    DualGainDoubleFunctionalTelescopeEventVisitor,
  calin::ix::diagnostics::functional::CameraDoubleFunctionalRawStats>;

%define ONEGAIN_TEMPLATE_WRAP(function)
%template(function)
  calin::diagnostics::functional:: ## function ##<
    calin::ix::diagnostics::functional::OneGainIntFunctionalRawStats>;
%template(function)
  calin::diagnostics::functional:: ## function ##<
    calin::ix::diagnostics::functional::OneGainDoubleFunctionalRawStats>;
%enddef

%define DUALGAIN_TEMPLATE_WRAP(function)
%template(function)
  calin::diagnostics::functional:: ## function ##<
    calin::ix::diagnostics::functional::CameraIntFunctionalRawStats>;
%template(function)
  calin::diagnostics::functional:: ## function ##<
    calin::ix::diagnostics::functional::CameraDoubleFunctionalRawStats>;
%enddef

ONEGAIN_TEMPLATE_WRAP(channel_mean)
ONEGAIN_TEMPLATE_WRAP(channel_var)
ONEGAIN_TEMPLATE_WRAP(channel_cov)
ONEGAIN_TEMPLATE_WRAP(channel_cov_frac)
ONEGAIN_TEMPLATE_WRAP(mean_of_mean_over_channels)
ONEGAIN_TEMPLATE_WRAP(var_of_mean_over_channels)
ONEGAIN_TEMPLATE_WRAP(decompose_channel_independent_and_common_var)

DUALGAIN_TEMPLATE_WRAP(channel_high_to_low_gain_cov)
DUALGAIN_TEMPLATE_WRAP(channel_high_to_low_gain_cov_frac)
