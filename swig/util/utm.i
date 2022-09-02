/*

   calin/util/utm.i -- Stephen Fegan -- 2022-07-31

   SWIG interface file for calin.util.utm

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

%module (package="calin.util") utm
%feature(autodoc,2);

%{
#include "util/datum.hpp"
#include "util/utm.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "calin_typemaps.i"
%import "calin_global_definitions.i"

%include "util/datum.hpp"

%apply double& OUTPUT { double& lat_rad, double& lon_rad, double& N, double& E,
  double& grid_convergence_rad, double& scale };
%apply int& INOUT { calin::util::utm::GridZone& zone, calin::util::utm::Hemisphere& hemi };

%include "util/utm.hpp"
