/*

   calin/diagnostics/value_capture.i -- Stephen Fegan -- 2016-05-14

   SWIG interface file for calin value capture diagnostics

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, Ecole polytechnique, CNRS/IN2P3, Universite Paris-Saclay

   This file is part of "calin"

   "calin" is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 or
   later, as published by the Free Software Foundation.

   "calin" is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

*/

%module (package="calin.diagnostics") value_capture

%{
#include "diagnostics/value_capture.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "calin_typemaps.i"
%import "calin_global_definitions.i"

%import "iact_data/event_dispatcher.i"
%import "diagnostics/value_capture.pb.i"
%include "diagnostics/value_capture.hpp"

%template (IntValueSupplierVisitor)
  calin::diagnostics::value_capture::ValueSupplierVisitor<int>;
%template (IntSequentialValueCaptureVisitor)
  calin::diagnostics::value_capture::SequentialValueCaptureVisitor<
    calin::diagnostics::value_capture::ValueSupplierVisitor<int>,
    calin::ix::diagnostics::value_capture::CapturedInt32Values>;

%template (LongValueSupplierVisitor)
  calin::diagnostics::value_capture::ValueSupplierVisitor<long>;
%template (LongSequentialValueCaptureVisitor)
  calin::diagnostics::value_capture::SequentialValueCaptureVisitor<
    calin::diagnostics::value_capture::ValueSupplierVisitor<long>,
    calin::ix::diagnostics::value_capture::CapturedInt64Values>;

%template (LongLongValueSupplierVisitor)
  calin::diagnostics::value_capture::ValueSupplierVisitor<long long>;
%template (LongLongSequentialValueCaptureVisitor)
  calin::diagnostics::value_capture::SequentialValueCaptureVisitor<
    calin::diagnostics::value_capture::ValueSupplierVisitor<long long>,
    calin::ix::diagnostics::value_capture::CapturedInt64Values>;

%template (DoubleValueSupplierVisitor)
  calin::diagnostics::value_capture::ValueSupplierVisitor<double>;
%template (DoubleSequentialValueCaptureVisitor)
  calin::diagnostics::value_capture::SequentialValueCaptureVisitor<
    calin::diagnostics::value_capture::ValueSupplierVisitor<double>,
    calin::ix::diagnostics::value_capture::CapturedDoubleValues>;
