/*

   calin/iact_data/nectarcam_layout.i -- Stephen Fegan -- 2016-06-06

   SWIG interface file for calin NectarCam camera layout

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
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

%module (package="calin.iact_data") nectarcam_layout
%feature(autodoc,2);

%{
#include "iact_data/nectarcam_layout.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "calin_typemaps.i"
%import "calin_global_definitions.i"

%newobject nectarcam_19module_layout(
  calin::ix::iact_data::instrument_layout::CameraLayout* layout = nullptr);
%newobject nectarcam_layout(
  calin::ix::iact_data::instrument_layout::CameraLayout* layout = nullptr);

%import "iact_data/instrument_layout.pb.i"
%include "iact_data/nectarcam_layout.hpp"
