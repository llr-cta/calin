//-*-mode:swig;-*-

/*

   calin/simulation/waveform_processor.i -- Stephen Fegan -- 2022-09-19

   SWIG interface file for calin.simulation.waveform_processor

   Copyright 2022, Stephen Fegan <sfegan@llr.in2p3.fr>
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

%module (package="calin.simulation") waveform_processor
%feature(autodoc,2);

%{
#include "simulation/ray_processor.hpp"
#include "simulation/vcl_ray_processor.hpp"
#include "simulation/vcl_ray_propagator.hpp"
#include "simulation/pe_processor.hpp"
#include "simulation/waveform_processor.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "calin_typemaps.i"
%import "calin_global_definitions.i"

%import "simulation/ray_processor.i"

%newobject new_trigger_memory_buffers;
%newobject new_trigger_patch_waveform_processor;
%include "simulation/waveform_processor.hpp"
