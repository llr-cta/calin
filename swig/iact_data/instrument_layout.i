/*

   calin/iact_data/instrument_layout.i -- Stephen Fegan -- 2017-09-15

   SWIG interface file for calin instrument (camera) layout

   Copyright 2017, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, CNRS, Ecole Polytechnique, Institut Polytechnique de Paris

   This file is part of "calin"

   "calin" is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 or
   later, as published by the Free Software Foundation.

   "calin" is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

*/

%module (package="calin.iact_data") instrument_layout
%feature(autodoc,2);

%{
#include "iact_data/instrument_layout.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "calin_typemaps.i"
%import "calin_global_definitions.i"

%newobject make_grid_from_instrument_layout;
%newobject reduce_camera_channels;
%newobject reduce_camera_modules;
%newobject channel_outline;

%apply const std::vector<unsigned int> & { const std::vector<unsigned int>& channel_id };
%apply const std::vector<unsigned int> & { const std::vector<unsigned int>& module_id };

%apply Eigen::VectorXi &OUTPUT { Eigen::VectorXi& map };

%import "iact_data/instrument_layout.pb.i"
%include "iact_data/instrument_layout.hpp"
