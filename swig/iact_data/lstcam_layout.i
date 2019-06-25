/*

   calin/iact_data/lstcam_layout.i -- Stephen Fegan -- 2018-10-16

   SWIG interface file for calin LSTCam camera layout

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

%module (package="calin.iact_data") lstcam_layout
%feature(autodoc,2);

%{
#include "iact_data/lstcam_layout.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "calin_typemaps.i"
%import "calin_global_definitions.i"

%newobject lstcam_layout(
  calin::ix::iact_data::instrument_layout::CameraLayout* layout = nullptr);

%import "iact_data/instrument_layout.pb.i"
%include "iact_data/lstcam_layout.hpp"
