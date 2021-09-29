/*

   calin/diagnostics/waveform.i -- Stephen Fegan -- 2016-03-23

   SWIG interface file for calin waveform diagnostics

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
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

%module (package="calin.diagnostics") waveform
%feature(autodoc,2);

%{
#include "util/log.hpp"
#include "diagnostics/waveform.hpp"
#include "diagnostics/waveform_psd_vcl.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "calin_typemaps.i"
%import "calin_global_definitions.i"

%import "util/log.i"
%import "iact_data/event_visitor.i"
%import "diagnostics/waveform.pb.i"

%apply Eigen::MatrixXd &OUTPUT { Eigen::MatrixXd& mean_waveform_out };
%apply Eigen::VectorXi &OUTPUT { Eigen::VectorXi& count_out };
%apply Eigen::VectorXd &OUTPUT { Eigen::VectorXd& h };
%newobject calin::diagnostics::waveform::WaveformSumParallelEventVisitor::mean_waveforms() const;
%newobject calin::diagnostics::waveform::WaveformCodeHistParallelEventVisitor::waveform_code_hist() const;
%newobject calin::diagnostics::waveform::WaveformCodeHistParallelEventVisitor::compact_waveform_code_hist() const;
%newobject calin::diagnostics::waveform::WaveformCodeHistParallelEventVisitor::compact_waveform_code_hist() const;
%newobject calin::diagnostics::waveform::WaveformPSDParallelVisitor::psd() const;

%include "diagnostics/waveform.hpp"
%include "diagnostics/waveform_psd_vcl.hpp"

%template(VCL128_WaveformPSDParallelVisitor)
  calin::diagnostics::waveform::VCL_WaveformPSDParallelVisitor<calin::util::vcl::VCL128Architecture>;
%template(VCL256_WaveformPSDParallelVisitor)
  calin::diagnostics::waveform::VCL_WaveformPSDParallelVisitor<calin::util::vcl::VCL256Architecture>;
