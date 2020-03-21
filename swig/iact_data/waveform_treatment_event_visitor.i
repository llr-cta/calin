/*

   calin/iact_data/waveform_treatment_event_visitor.i -- Stephen Fegan -- 2018-01-12

   SWIG interface file for waveform treatment event visitor

   Copyright 2018, Stephen Fegan <sfegan@llr.in2p3.fr>
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

%module (package="calin.iact_data") waveform_treatment_event_visitor
%feature(autodoc,2);

%{
#include "iact_data/event_visitor.hpp"
#include "iact_data/waveform_treatment_event_visitor.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "calin_typemaps.i"
%import "calin_global_definitions.i"

%import "event_visitor.i"
%import "iact_data/waveform_treatment_event_visitor.pb.i"
%include "iact_data/waveform_treatment_event_visitor.hpp"

%template (VCL_OptimalWindowSumWaveformTreatmentEventVisitor256)
  calin::iact_data::waveform_treatment_event_visitor::VCL_OptimalWindowSumWaveformTreatmentEventVisitor<calin::util::vcl::VCL256Architecture>;
%template (VCL_OptimalWindowSumWaveformTreatmentEventVisitor128)
  calin::iact_data::waveform_treatment_event_visitor::VCL_OptimalWindowSumWaveformTreatmentEventVisitor<calin::util::vcl::VCL128Architecture>;

%template (VCLSingleGainDualWindowWaveformTreatmentEventVisitor256)
  calin::iact_data::waveform_treatment_event_visitor::VCL_SingleGainDualWindowWaveformTreatmentEventVisitor<calin::util::vcl::VCL256Architecture>;
%template (VCLSingleGainDualWindowWaveformTreatmentEventVisitor128)
  calin::iact_data::waveform_treatment_event_visitor::VCL_SingleGainDualWindowWaveformTreatmentEventVisitor<calin::util::vcl::VCL128Architecture>;
