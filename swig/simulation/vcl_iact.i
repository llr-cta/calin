/*

   calin/simulation/vcl_iact.i -- Stephen Fegan -- 2019-02-20

   SWIG interface file for calin.simulation.vcl_iact

   Copyright 2019, Stephen Fegan <sfegan@llr.in2p3.fr>
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

%module (package="calin.simulation") vcl_iact
%feature(autodoc,2);

%{
#include "simulation/tracker.hpp"
#include "simulation/misc_trackers.hpp"
#include "simulation/vcl_iact.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "calin_typemaps.i"
%import "calin_global_definitions.i"

%import "simulation/tracker.i"

%include "simulation/vcl_iact.hpp"

%template (VCLIACTTrackVisitor128)
  calin::simulation::vcl_iact::VCLIACTTrackVisitor<calin::util::vcl::VCL128Architecture>;
%template (VCLIACTTrackVisitor256)
  calin::simulation::vcl_iact::VCLIACTTrackVisitor<calin::util::vcl::VCL256Architecture>;
%template (VCLIACTTrackVisitor512)
  calin::simulation::vcl_iact::VCLIACTTrackVisitor<calin::util::vcl::VCL512Architecture>;
